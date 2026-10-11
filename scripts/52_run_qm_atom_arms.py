#!/usr/bin/env python
"""Notebook 52's driver: Chemprop + CheMeleon with the QM data at BOTH levels.

THE ARM (`qm_both`), one object, used identically in stages 1 and 2:
  * graph input -- `canonical_smiles`, chemprop's own atom/bond features, V2 featurizer
  * molecule level -- notebook 44's 13 genuinely-new QM columns via `--descriptors-columns`
  * atom level -- `q_hirshfeld` and `mayer_va` per atom via `--atom-features-path`, which
    chemprop 2.3.1 hstacks onto the atom feature matrix BEFORE message passing
  * encoder -- CheMeleon, widened by two zero-initialised columns so it can accept them;
    see `scripts/52_atom_arrays.py` for why the published route is closed

Notebook 45's four AGGREGATE atom descriptors are deliberately absent: they are summaries of
these same per-atom values, so feeding both would be feeding one quantity twice.
No bond features: bond ordering is a second ordering problem and is out of scope here.

Modes, in order:
  --sanity   Plain CheMeleon on the RANDOM repeat_0/fold_0, reproducing notebook 05's stored
             `chemprop_chemeleoninit__repeat0_fold0` scores. The tenth notebook in a row to
             run this check. Aborts everything downstream on failure.
  --wiring   The gate that decides whether stage 1 is worth 4-5 hours. Two fits on one Butina
             fold, differing ONLY in the VALUES in the atom .npz -- real against all-zeros --
             with everything else (CLI, widened encoder, molecule descriptors, seed, folds)
             held identical. If the blind predictions do not differ, the features are not
             connected and the CV would train an expensive copy of `plain`.
  --confirm  The 25 frozen cluster-disjoint folds. Compared against notebook 37's STORED
             `plain` scores on the identical folds, never a refit.
  --final    One fit on all 4,905 training compounds, predicting the blind 750.

The harness is notebook 44's (itself notebook 37's), with the differences the arm forces:
a local foundation file, an atom-features npz per CSV, `--no-atom-feature-scaling`, and
`--patience` always passed explicitly.

Usage:
  python scripts/52_run_qm_atom_arms.py --sanity
  python scripts/52_run_qm_atom_arms.py --wiring
  python scripts/52_run_qm_atom_arms.py --confirm [--repeats 0,1]
  python scripts/52_run_qm_atom_arms.py --final
"""
from __future__ import annotations

import os

# THREAD ENVIRONMENT. Notebook 34 established that restricting BLAS/OMP to one thread changes
# floating-point reduction order and so chemprop's training trajectory, breaking
# bit-reproducibility against notebook 05's stored scores. OMP_NUM_THREADS is asserted unset
# because nothing sets it implicitly; KMP_DUPLICATE_LIB_OK is only reported, because
# `import xgboost` sets it and an assertion after imports is unsatisfiable (notebook 40's bug).
assert "OMP_NUM_THREADS" not in os.environ, (
    "OMP_NUM_THREADS is set. It changes chemprop's training trajectory and breaks "
    "bit-reproducibility against notebook 05's stored scores. Unset it and re-run.")
_KMP_AT_START = os.environ.get("KMP_DUPLICATE_LIB_OK")

import argparse  # noqa: E402
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from datetime import datetime, timezone  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.chemprop_screen import (  # noqa: E402
    build_predict_csv, build_training_csv, run_chemprop_predict, run_chemprop_train,
    run_subprocess_streamed, setup_logging, verify_predictions,
)
from src.cv_bootstrap import per_fold_bootstrap_seed  # noqa: E402
from src.features import assign_screen_split  # noqa: E402
from src.vendor.openadmet_eval.config import ACTIVITY_METRICS, REGRESSION_ENDPOINTS  # noqa: E402
from src.vendor.openadmet_eval.evaluate_predictions import (  # noqa: E402
    add_macro_endpoint, score_activity_predictions,
)

METRIC_NAMES = [name for name, _ in ACTIVITY_METRICS]
ISOFORMS = [e.split("_")[0] for e in REGRESSION_ENDPOINTS]
BASE_COLS = ["Molecule_Name", "inchikey", "canonical_smiles"]

OUT = REPO_ROOT / "outputs" / "52_chemprop_qm_atom_features"
PRED_DIR = OUT / "predictions"
SCORE_DIR = OUT / "scores"
OOF_DIR = OUT / "oof_long"
RUNS_DIR = OUT / "chemprop_runs"
LOG_DIR = REPO_ROOT / "logs" / "52_chemprop_qm_atom_features"

FOLDS_PATH = REPO_ROOT / "data" / "folds" / "cv_folds.csv"
BUTINA5X5_PATH = REPO_ROOT / "data" / "folds" / "cv_folds_butina_5x5.csv"
CURATED_PATH = REPO_ROOT / "data" / "processed" / "train_inhibition_curated.csv"
BLIND_PATH = REPO_ROOT / "data" / "processed" / "test_blinded_curated.csv"
QM_FULL_PATH = REPO_ROOT / "outputs" / "44_qm_molecular_block" / "features" / "QM_FULL.csv"
NB05_OUT = REPO_ROOT / "outputs" / "05_cv_comparison"
MANIFEST_PATH = NB05_OUT / "manifest.csv"

NB34_FOLD_SHA256 = "52f39b1a1cea9200b58dce0a061d4a21ce94d9d58ec77067cf6508b0f4b3e87e"
VAL_FRACTION = 0.15
EPOCHS, PATIENCE = 50, 5
CHEMELEON_ARCH_ARGS = ["--from-foundation", "CHEMELEON",
                       "--multi-hot-atom-featurizer-mode", "V2"]
SANITY_THRESHOLD = 0.03
N_FOLDS = 5
SEEDS_3 = [2684470948, 4091952314, 233227757]
FINAL_SEED = 2684470948

# QM_FULL's five algebraic derivatives, excluded to leave notebook 44's 13 genuinely-new columns.
DERIVED5 = ["dft_gap_ev", "gfn2_gap_ev", "chi_electronegativity", "omega_electrophilicity",
            "dft_minus_gfn2_gap"]

ARM = "qm_both"
ARM_ZEROATOMS = "qm_both_zeroatoms"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def atom_arrays_module():
    """`scripts/52_atom_arrays.py`, imported rather than duplicated."""
    spec = importlib.util.spec_from_file_location(
        "nb52_atom_arrays", REPO_ROOT / "scripts" / "52_atom_arrays.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def thirteen_qm_columns(qm: pd.DataFrame) -> list[str]:
    cols = [c for c in qm.columns if c not in ("Molecule_Name", "inchikey", "set", "split",
                                               *DERIVED5)]
    assert len(cols) == 13, (len(cols), cols)
    return cols


# ------------------------------------------------------------------ CSV + npz builders
def build_training_csv_multi(population, target_cols, descriptor_cols, out_path, logger):
    """Notebook 44's builder. `src/chemprop_screen.py`'s own hardcodes its column list and has
    no hook for extra descriptors; it is deliberately left unmodified."""
    pooled = population[population["screen_split"].isin(
        ["screen_inner_train", "screen_inner_val"])].copy()
    pooled["chemprop_split"] = pooled["screen_split"].map(
        {"screen_inner_train": "train", "screen_inner_val": "val"})
    assert pooled["chemprop_split"].notna().all()
    cols = ["Molecule_Name", "canonical_smiles", *target_cols, *descriptor_cols, "chemprop_split"]
    pooled[cols].to_csv(out_path, index=False)
    logger.info(f"wrote {out_path} ({len(pooled)} rows, {len(descriptor_cols)} descriptors, "
                f"splits {pooled['chemprop_split'].value_counts().to_dict()})")
    return pooled


def build_predict_csv_multi(frame, descriptor_cols, out_path, logger):
    cols = [c for c in BASE_COLS if c in frame.columns] + list(descriptor_cols)
    frame[cols].to_csv(out_path, index=False)
    logger.info(f"wrote {out_path} ({len(frame)} compounds)")
    return out_path


def run_chemprop_predict_atoms(model_dir, predict_csv, descriptor_cols, npz, output_csv, logger):
    argv = ["chemprop", "predict", "-i", str(predict_csv), "-s", "canonical_smiles",
            "--descriptors-columns", *descriptor_cols,
            "--model-paths", str(model_dir), "--accelerator", "cpu", "--devices", "1",
            "-o", str(output_csv)]
    if npz is not None:
        argv += ["--atom-features-path", str(npz), "--no-atom-feature-scaling"]
    run_subprocess_streamed(argv, logger)


# ------------------------------------------------------------------ shared state
class Ctx:
    def __init__(self):
        for d in [OUT, PRED_DIR, SCORE_DIR, OOF_DIR, RUNS_DIR, LOG_DIR]:
            d.mkdir(parents=True, exist_ok=True)
        self.curated = pd.read_csv(CURATED_PATH)
        assert len(self.curated) == 4905
        self.blind = pd.read_csv(BLIND_PATH)
        assert len(self.blind) == 750

        sha = sha256(BUTINA5X5_PATH)
        assert sha == NB34_FOLD_SHA256, (
            f"cv_folds_butina_5x5.csv sha256 {sha} != {NB34_FOLD_SHA256} -- the partition is not "
            f"the file notebooks 34/35/37 ran on, and every comparison would be invalid.")
        self.butina = pd.read_csv(BUTINA5X5_PATH)
        self.butina_cols = [c for c in self.butina.columns if c.startswith("butina_repeat_")]
        assert len(self.butina_cols) == 5
        self.random_folds = pd.read_csv(FOLDS_PATH)

        manifest = pd.read_csv(MANIFEST_PATH)
        mf = manifest[manifest["config"] == "chemprop_chemeleoninit"].sort_values(
            ["repeat", "fold"])
        assert len(mf) == 25
        self.fold_seed = {(int(r.repeat), int(r.fold)): int(r.seed) for r in mf.itertuples()}
        assert [self.fold_seed[(0, f)] for f in range(3)] == SEEDS_3, \
            "repeat_0's first three manifest seeds are not this project's established SEEDS"
        for nb in ["34_butina_5x5", "35_butina_repeat3", "37_butina_5x5_complete"]:
            p = REPO_ROOT / "outputs" / nb / "run_scores.csv"
            if not p.exists():
                continue
            rs = pd.read_csv(p)[["repeat", "fold", "seed"]].drop_duplicates()
            bad = [(int(r.repeat), int(r.fold), int(r.seed)) for r in rs.itertuples()
                   if self.fold_seed[(int(r.repeat), int(r.fold))] != int(r.seed)]
            assert not bad, f"seed mapping differs from notebook {nb}'s: {bad}"

        qm = pd.read_csv(QM_FULL_PATH)
        self.desc_cols = thirteen_qm_columns(qm)
        self.qm = qm[["Molecule_Name"] + self.desc_cols]
        assert len(self.qm) == 5655 and self.qm[self.desc_cols].notna().all().all()

        aa = atom_arrays_module()
        self.aa = aa
        self.src = aa.AtomFeatureSource()
        self.foundation = OUT / "chemeleon_widened.pt"
        if not self.foundation.exists():
            info = aa.build_widened_chemeleon(aa.K_ATOM_FEATURES, self.foundation)
            (OUT / "widened_encoder.json").write_text(json.dumps(info, indent=2) + "\n")
        self.arm_args = ["--from-foundation", str(self.foundation),
                         "--multi-hot-atom-featurizer-mode", "V2",
                         "--loss-function", "mse"]

        self.pop_base = self.curated[BASE_COLS + list(REGRESSION_ENDPOINTS)].merge(
            self.qm, on="Molecule_Name", validate="one_to_one")
        assert len(self.pop_base) == 4905
        self.blind_frame = self.blind[BASE_COLS].merge(
            self.qm, on="Molecule_Name", validate="one_to_one")
        assert len(self.blind_frame) == 750

        self.scores_path = OUT / "run_scores.csv"
        self.timings_path = OUT / "run_timings.csv"

    def population(self, partition, repeat_col, fold, seed):
        folds = self.butina if partition == "butina" else self.random_folds
        split_df = assign_screen_split(folds, repeat_col=repeat_col, test_fold=fold,
                                       val_fraction=VAL_FRACTION, seed=seed)
        pop = self.pop_base.merge(split_df[["inchikey", "screen_split"]], on="inchikey",
                                  how="left")
        assert len(pop) == len(self.pop_base), "split merge changed row count"
        pop.loc[pop["screen_split"].isna(), "screen_split"] = "screen_inner_train"
        assert pop["screen_split"].isin(
            ["screen_test", "screen_inner_train", "screen_inner_val"]).all()
        return pop


def _append_rows(path, new_rows, key_cols):
    """Notebook 37's upsert, verbatim."""
    new_df = pd.DataFrame(new_rows)
    if path.exists():
        old = pd.read_csv(path)
        keys = set(map(tuple, new_df[key_cols].values.tolist()))
        mask = ~old[key_cols].apply(lambda r: tuple(r) in keys, axis=1)
        combined = pd.concat([old[mask], new_df], ignore_index=True)
    else:
        combined = new_df
    combined = combined.sort_values(key_cols).reset_index(drop=True)
    combined.to_csv(path, index=False)
    return combined


def cached(pred_path, expected_names):
    """Notebook 37's content-based cache check. Existence is never enough."""
    if not pred_path.exists():
        return None
    try:
        df = pd.read_csv(pred_path)
    except Exception:
        return None
    if len(df) != len(expected_names) or set(df["Molecule_Name"]) != expected_names:
        print(f"  [cache MISMATCH] {pred_path.name}: content check failed -- will retrain")
        return None
    if df[list(REGRESSION_ENDPOINTS)].isna().to_numpy().any():
        print(f"  [cache MISMATCH] {pred_path.name}: NaN in a predicted column -- will retrain")
        return None
    return df


def score_and_append(ctx, arm, repeat, fold, seed, pred_df, tag):
    gt = ctx.curated[ctx.curated["inchikey"].isin(pred_df["inchikey"])].copy()
    with per_fold_bootstrap_seed(seed):
        scored = score_activity_predictions(pred_df, gt, list(REGRESSION_ENDPOINTS))
        scored = add_macro_endpoint(scored, list(REGRESSION_ENDPOINTS), ACTIVITY_METRICS)
    scored.to_csv(SCORE_DIR / f"{tag}.csv", index=False)
    pe = scored.groupby("Endpoint")[METRIC_NAMES].mean()
    rows = []
    for endpoint in pe.index:
        row = {"repeat": repeat, "fold": fold, "arm": arm, "seed": seed, "Endpoint": endpoint}
        row.update(pe.loc[endpoint].to_dict())
        rows.append(row)
    _append_rows(ctx.scores_path, rows, ["repeat", "fold", "arm", "seed", "Endpoint"])
    return pe


def save_oof_long(ctx, arm, repeat, fold, seed, pred_df, tag):
    frames = []
    for iso in ISOFORMS:
        ep = f"{iso}_pIC50_direct_inhibition"
        m = pred_df[["inchikey", ep]].rename(columns={ep: "y_pred"}).copy()
        m["isoform"] = iso
        frames.append(m[["inchikey", "isoform", "y_pred"]])
    out = pd.concat(frames, ignore_index=True)
    for col, val in [("seed", seed), ("fold", fold), ("repeat", repeat), ("arm", arm)]:
        out.insert(0, col, val)
    out.to_csv(OOF_DIR / f"{tag}.csv", index=False)


def record_timing(ctx, arm, repeat, fold, seed, wall, last_epoch=None, best_epoch=None):
    """Notebook 37's anti-clobber rule: a cache hit must not overwrite a completed run's timing.
    Notebook 35 lost its whole timing file exactly that way."""
    row = {"repeat": repeat, "fold": fold, "arm": arm, "seed": seed, "wall_min": wall,
           "last_epoch": last_epoch, "best_epoch": best_epoch,
           "min_per_epoch": (wall / last_epoch) if last_epoch else np.nan}
    key = ["repeat", "fold", "arm", "seed"]
    already = False
    if ctx.timings_path.exists():
        t = pd.read_csv(ctx.timings_path)
        already = bool(((t["repeat"] == repeat) & (t["fold"] == fold) & (t["arm"] == arm)
                        & (t["seed"] == seed)).any())
    if wall < 0.02 and already:
        print(f"  [cache hit] keeping the existing timing record for {arm}")
        return
    _append_rows(ctx.timings_path, [row], key)


def epochs_of(run_dir: Path):
    """(last_epoch, best_epoch). `best_epoch` comes from chemprop's own checkpoint FILENAME,
    which is where notebook 47 read the stopping count from."""
    m = run_dir / "model_0" / "trainer_logs" / "version_0" / "metrics.csv"
    last = None
    if m.exists():
        eps = [int(r.split(",")[0]) for r in m.read_text().splitlines()[1:]
               if r and r.split(",")[0].isdigit()]
        last = max(eps) if eps else None
    best = None
    ck = sorted((run_dir / "model_0" / "checkpoints").glob("best-epoch=*.ckpt"))
    if ck:
        best = int(ck[0].name.split("best-epoch=")[1].split("-")[0])
    return last, best


# ------------------------------------------------------------------ the arm
def fit_arm(ctx, tag, pop, zero_atoms: bool, seed: int, predict_frames: dict,
            epochs: int = EPOCHS, patience: int = PATIENCE, split_col_values=None):
    """One fit of the arm, plus one prediction per entry of `predict_frames`.

    `predict_frames` maps a label to a DataFrame of compounds to predict. Returns
    {label: DataFrame}. Everything is keyed on Molecule_Name, so the npz row order is built
    from the CSV that is actually written rather than assumed.
    """
    run_dir = RUNS_DIR / tag
    run_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(LOG_DIR / f"{tag}.log", f"nb52_{tag}")
    logger.info(f"=== {tag}: 13 molecule descriptors + "
                f"{ctx.aa.K_ATOM_FEATURES} atom features "
                f"{ctx.aa.ATOM_FEATURE_COLS}, zero_atoms={zero_atoms} ===")

    t0 = time.time()
    train_csv = run_dir / "train_input.csv"
    pooled = build_training_csv_multi(pop, list(REGRESSION_ENDPOINTS), ctx.desc_cols,
                                      train_csv, logger)

    # The atom-feature scaler is fitted on `screen_inner_train` atoms ONLY -- the population
    # chemprop's own scaler would have used, excluding the inner validation set and the
    # held-out fold. See scripts/52_atom_arrays.py for why the scaling is done here.
    fit_names = pooled.loc[pooled["screen_split"] == "screen_inner_train", "Molecule_Name"]
    mean, scale = ctx.src.fit_scaler(fit_names)
    logger.info(f"atom-feature scaler fitted on {len(fit_names)} screen_inner_train compounds: "
                f"loc={np.array2string(mean, precision=4)}, "
                f"scale={np.array2string(scale, precision=4)}")
    (run_dir / "atom_scaler.json").write_text(json.dumps(
        {"columns": ctx.aa.ATOM_FEATURE_COLS, "n_fit_compounds": int(len(fit_names)),
         "loc": [float(v) for v in mean], "scale": [float(v) for v in scale]}, indent=2) + "\n")

    train_order = pd.read_csv(train_csv, usecols=["Molecule_Name"])["Molecule_Name"].tolist()
    train_npz = run_dir / "atom_features_train.npz"
    info = (ctx.src.write_npz_zeroed(train_order, train_npz) if zero_atoms
            else ctx.src.write_npz(train_order, mean, scale, train_npz))
    logger.info(f"train npz: {json.dumps(info)}")

    argv_extra = ctx.arm_args + ["--descriptors-columns", *ctx.desc_cols,
                                 "--atom-features-path", str(train_npz),
                                 "--no-atom-feature-scaling"]
    run_chemprop_train(train_csv, list(REGRESSION_ENDPOINTS), run_dir, logger, argv_extra,
                       epochs, patience, seed)

    results = {}
    for label, frame in predict_frames.items():
        pcsv = run_dir / f"predict_input_{label}.csv"
        build_predict_csv_multi(frame, ctx.desc_cols, pcsv, logger)
        order = pd.read_csv(pcsv, usecols=["Molecule_Name"])["Molecule_Name"].tolist()
        pnpz = run_dir / f"atom_features_{label}.npz"
        pinfo = (ctx.src.write_npz_zeroed(order, pnpz) if zero_atoms
                 else ctx.src.write_npz(order, mean, scale, pnpz))
        logger.info(f"{label} npz: {json.dumps(pinfo)}")
        raw = run_dir / f"raw_predictions_{label}.csv"
        run_chemprop_predict_atoms(run_dir / "model_0", pcsv, ctx.desc_cols, pnpz, raw, logger)
        pred = pd.read_csv(raw)
        assert set(pred["Molecule_Name"]) == set(order), f"{label}: Molecule_Name set changed"
        assert not pred[list(REGRESSION_ENDPOINTS)].isna().to_numpy().any(), \
            f"{label}: NaN in a predicted column"
        results[label] = pred

    wall = (time.time() - t0) / 60.0
    last, best = epochs_of(run_dir)
    logger.info(f"{tag}: {wall:.2f} min (last_epoch={last}, best_epoch={best})")
    return results, wall, last, best


def run_confirm_unit(ctx, repeat, fold):
    seed = ctx.fold_seed[(repeat, fold)]
    tag = f"{ARM}__r{repeat}f{fold}"
    pop = ctx.population("butina", ctx.butina_cols[repeat], fold, seed)
    held = pop.loc[pop["screen_split"] == "screen_test"].copy()
    expected = set(held["Molecule_Name"])
    pred_path = PRED_DIR / f"{tag}.csv"
    hit = cached(pred_path, expected)
    if hit is not None:
        score_and_append(ctx, ARM, repeat, fold, seed, hit, tag)
        save_oof_long(ctx, ARM, repeat, fold, seed, hit, tag)
        record_timing(ctx, ARM, repeat, fold, seed, 0.0)
        print(f"  [cache hit] {tag}")
        return hit, 0.0
    res, wall, last, best = fit_arm(ctx, tag, pop, False, seed, {"heldout": held})
    pred = res["heldout"]
    pred.to_csv(pred_path, index=False)
    score_and_append(ctx, ARM, repeat, fold, seed, pred, tag)
    save_oof_long(ctx, ARM, repeat, fold, seed, pred, tag)
    record_timing(ctx, ARM, repeat, fold, seed, wall, last, best)
    print(f"  {tag}: {wall:.2f} min (last_epoch={last}, best_epoch={best})")
    return pred, wall


# ------------------------------------------------------------------ modes
def mode_sanity(ctx):
    from src.chemprop_screen import load_screen_population
    print("=" * 78)
    print("SANITY GATE -- reproduce notebook 05's stored chemprop_chemeleoninit repeat0_fold0")
    print("=" * 78)
    tag = "sanity_plain__random_r0f0"
    seed = ctx.fold_seed[(0, 0)]
    logger = setup_logging(LOG_DIR / f"{tag}.log", f"nb52_{tag}")
    pop = load_screen_population(FOLDS_PATH, CURATED_PATH, "repeat_0", 0, VAL_FRACTION, seed,
                                 logger)
    expected = set(pop.loc[pop["screen_split"] == "screen_test", "Molecule_Name"])
    pred_path = PRED_DIR / f"{tag}.csv"
    pred = cached(pred_path, expected)
    if pred is not None:
        print(f"[cache hit] {tag}")
    else:
        run_dir = RUNS_DIR / tag
        run_dir.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        train_csv = run_dir / "train_input.csv"
        build_training_csv(pop, list(REGRESSION_ENDPOINTS), train_csv, logger,
                           require_all_targets=False)
        predict_csv = build_predict_csv(pop, run_dir, logger)
        run_chemprop_train(train_csv, list(REGRESSION_ENDPOINTS), run_dir, logger,
                           CHEMELEON_ARCH_ARGS + ["--loss-function", "mse"], EPOCHS, PATIENCE,
                           seed)
        raw = run_dir / "raw_predictions.csv"
        run_chemprop_predict(run_dir / "model_0", predict_csv, raw, logger)
        verify_predictions(raw, expected, list(REGRESSION_ENDPOINTS), logger)
        pred = pd.read_csv(raw)
        pred.to_csv(pred_path, index=False)
        print(f"sanity run wall time: {(time.time() - t0) / 60:.2f} min")

    gt = ctx.curated[ctx.curated["inchikey"].isin(pred["inchikey"])].copy()
    with per_fold_bootstrap_seed(seed):
        scored = score_activity_predictions(pred, gt, list(REGRESSION_ENDPOINTS))
    mine = scored.groupby("Endpoint")["ST-RAE"].mean()
    theirs = pd.read_csv(
        NB05_OUT / "scores" / "chemprop_chemeleoninit__repeat0_fold0.csv"
    ).groupby("Endpoint")["ST-RAE"].mean()
    print(f"\n{'endpoint':42s} {'this run':>14s} {'nb05':>14s} {'abs diff':>10s}")
    rows, max_diff = [], 0.0
    for ep in sorted(set(mine.index) & set(theirs.index)):
        d = abs(float(mine[ep]) - float(theirs[ep]))
        max_diff = max(max_diff, d)
        print(f"{ep:42s} {mine[ep]:14.10f} {theirs[ep]:14.10f} {d:10.2e}")
        rows.append({"Endpoint": ep, "this_run": float(mine[ep]), "nb05": float(theirs[ep]),
                     "abs_diff": d})
    pd.DataFrame(rows).to_csv(OUT / "sanity_check.csv", index=False)
    print(f"\nmax abs ST-RAE difference: {max_diff:.10f}   threshold {SANITY_THRESHOLD}")
    (OUT / "sanity_verdict.json").write_text(json.dumps(
        {"max_abs_strae_diff": max_diff, "threshold": SANITY_THRESHOLD,
         "bit_identical": max_diff == 0.0, "pass": bool(max_diff <= SANITY_THRESHOLD),
         "at_utc": datetime.now(timezone.utc).isoformat()}, indent=2) + "\n")
    if max_diff > SANITY_THRESHOLD:
        print("*** SANITY GATE FAILED -- nothing downstream is interpretable. STOPPING. ***")
        return 1
    print("PASS" + (" (bit-identical)" if max_diff == 0.0 else ""))
    return 0


def mode_wiring(ctx):
    print("=" * 78)
    print("WIRING CHECK -- do the per-atom VALUES reach the model?")
    print("two fits on butina repeat_0/fold_0, identical but for the npz CONTENTS")
    print("=" * 78)
    repeat, fold = 0, 0
    seed = ctx.fold_seed[(repeat, fold)]
    pop = ctx.population("butina", ctx.butina_cols[repeat], fold, seed)
    held = pop.loc[pop["screen_split"] == "screen_test"].copy()
    frames = {"heldout": held, "blind": ctx.blind_frame}

    out = {}
    for arm, zero in [(ARM, False), (ARM_ZEROATOMS, True)]:
        tag = f"wiring_{arm}__r{repeat}f{fold}"
        cache = PRED_DIR / f"{tag}_blind.csv"
        hcache = PRED_DIR / f"{tag}_heldout.csv"
        if cache.exists() and hcache.exists():
            print(f"  [cache hit] {tag}")
            out[arm] = {"blind": pd.read_csv(cache), "heldout": pd.read_csv(hcache)}
            continue
        print(f"\n--- {tag} (zero_atoms={zero}) ---")
        res, wall, last, best = fit_arm(ctx, tag, pop, zero, seed, frames)
        res["blind"].to_csv(cache, index=False)
        res["heldout"].to_csv(hcache, index=False)
        record_timing(ctx, f"wiring_{arm}", repeat, fold, seed, wall, last, best)
        print(f"  {tag}: {wall:.2f} min (last_epoch={last}, best_epoch={best})")
        out[arm] = res

    rows = []
    for label in ["blind", "heldout"]:
        a = out[ARM][label].set_index("Molecule_Name").sort_index()
        b = out[ARM_ZEROATOMS][label].set_index("Molecule_Name").sort_index()
        assert list(a.index) == list(b.index)
        for ep in REGRESSION_ENDPOINTS:
            d = (a[ep].to_numpy() - b[ep].to_numpy())
            rows.append({"set": label, "Endpoint": ep, "n": len(d),
                         "max_abs_diff": float(np.abs(d).max()),
                         "mean_abs_diff": float(np.abs(d).mean()),
                         "pearson": float(np.corrcoef(a[ep], b[ep])[0, 1])})
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "wiring_check.csv", index=False)
    print()
    print(table.to_string(index=False))

    blind_max = float(table.loc[table["set"] == "blind", "max_abs_diff"].max())
    verdict = {
        "at_utc": datetime.now(timezone.utc).isoformat(),
        "arms": [ARM, ARM_ZEROATOMS],
        "difference_is_only_the_npz_contents": True,
        "blind_max_abs_diff": blind_max,
        "heldout_max_abs_diff": float(table.loc[table["set"] == "heldout",
                                                "max_abs_diff"].max()),
        "wired": bool(blind_max > 0.0),
    }
    (OUT / "wiring_verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")
    print(f"\nmax abs difference on the blind 750: {blind_max:.6f}")
    if blind_max == 0.0:
        print("*** NOT WIRED -- the per-atom values are not reaching the model. STOPPING. ***")
        return 1
    print("WIRED -- the per-atom values change the predictions.")
    return 0


def mode_confirm(ctx, repeats):
    print("=" * 78)
    print(f"CONFIRM -- {ARM} on the frozen cluster-disjoint partition, repeats {repeats}")
    print("compared against notebook 37's STORED `plain` scores, never a refit")
    print("=" * 78)
    t_start = time.time()
    done, walls = 0, []
    total = len(repeats) * N_FOLDS
    for repeat in repeats:
        for fold in range(N_FOLDS):
            print(f"\n[{done + 1}/{total}] repeat {repeat} fold {fold}")
            _, wall = run_confirm_unit(ctx, repeat, fold)
            done += 1
            if wall > 0:
                walls.append(wall)
            if len(walls) == 3:
                rate = float(np.mean(walls))
                print(f"\n  >>> MEASURED RATE after 3 real folds: {rate:.2f} min/fold "
                      f"-> {rate * 25 / 60:.2f} h for 25 folds "
                      f"(the brief projected 9-12 min/fold, 4-5 h)")
            if walls:
                print(f"  elapsed {(time.time() - t_start) / 60:.1f} min, "
                      f"mean {np.mean(walls):.2f} min/fold")
    print(f"\n{ARM}: {done} units, {(time.time() - t_start) / 60:.1f} min total")
    return 0


def mode_final(ctx, epochs):
    print("=" * 78)
    print(f"FINAL -- one fit on all 4,905 training compounds, {epochs} epochs, "
          f"early stopping disabled")
    print("=" * 78)
    tag = f"{ARM}__final"
    # 100% of the training compounds, zero held out -- notebook 47's own choice. chemprop falls
    # its tracking metric back to train_loss when no row carries `val` (verified against
    # cli/train.py's build_datasets), and `patience = epochs` is byte-for-byte what chemprop's
    # own default resolves to when --patience is omitted.
    pop = ctx.pop_base.copy()
    pop["screen_split"] = "screen_inner_train"
    frames = {"blind": ctx.blind_frame}
    pred_path = PRED_DIR / f"{tag}_blind.csv"
    if pred_path.exists() and len(pd.read_csv(pred_path)) == 750:
        print(f"[cache hit] {tag}")
        pred = pd.read_csv(pred_path)
        wall, last, best = 0.0, None, None
    else:
        res, wall, last, best = fit_arm(ctx, tag, pop, False, FINAL_SEED, frames,
                                        epochs=epochs, patience=epochs)
        pred = res["blind"]
        pred.to_csv(pred_path, index=False)
    record_timing(ctx, f"{ARM}__final", -1, -1, FINAL_SEED, wall, last, best)
    print(f"\nwrote {pred_path.relative_to(REPO_ROOT)} ({len(pred)} rows)")
    print(pred[list(REGRESSION_ENDPOINTS)].describe().round(4).to_string())
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sanity", action="store_true")
    ap.add_argument("--wiring", action="store_true")
    ap.add_argument("--confirm", action="store_true")
    ap.add_argument("--final", action="store_true")
    ap.add_argument("--repeats", default="0,1,2,3,4")
    ap.add_argument("--epochs", type=int, default=None, help="--final only")
    args = ap.parse_args()
    if not (args.sanity or args.wiring or args.confirm or args.final):
        ap.error("pick a mode")

    print(f"started {datetime.now(timezone.utc).isoformat(timespec='seconds')}")
    print(f"KMP_DUPLICATE_LIB_OK at start: {_KMP_AT_START!r} (reported, never asserted)")
    print(f"chemprop: "
          f"{subprocess.run(['which', 'chemprop'], capture_output=True, text=True).stdout.strip()}")
    ctx = Ctx()
    print(f"13 molecule-level QM columns: {ctx.desc_cols}")
    print(f"{ctx.aa.K_ATOM_FEATURES} per-atom QM columns: {ctx.aa.ATOM_FEATURE_COLS}")
    print(f"widened foundation: {ctx.foundation.name} sha256 {sha256(ctx.foundation)[:16]}")

    rc = 0
    if args.sanity:
        rc = mode_sanity(ctx)
    if rc == 0 and args.wiring:
        rc = mode_wiring(ctx)
    if rc == 0 and args.confirm:
        rc = mode_confirm(ctx, [int(x) for x in args.repeats.split(",")])
    if rc == 0 and args.final:
        if args.epochs is None:
            ap.error("--final needs --epochs (the median per-fold best epoch)")
        rc = mode_final(ctx, args.epochs)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
