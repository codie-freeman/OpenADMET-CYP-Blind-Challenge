#!/usr/bin/env python
"""Train the molecule-level QM arms. Notebook 44's driver.

Three modes, meant to be run in this order:

  --sanity    Retrain chemprop_chemeleoninit on the RANDOM repeat_0/fold_0 and reproduce
              outputs/05_cv_comparison's stored scores. Must be 0.0000 on all four isoforms.
              Would be the tenth consecutive piece of work to pass this check (18/20/22/28/29/
              31/34/35/37). Aborts everything downstream on failure.

  --tier1     The cheap screen: RANDOM partition (data/folds/cv_folds.csv), repeat_0 / fold 0,
              three seeds, 14 blocks x 4 model classes, then the 20-draw descriptive null band
              at the selected model class. Deliberately NOT on the Butina partition, so the
              frozen 25-fold design stays unseen until confirmation -- a departure from
              notebooks 18/20/22/28/31, where notebook 29 then found CYP2C9's screen verdict was
              a single-fold selection artefact. Writes model_selection.json (M*, chosen on the
              CONTROL block before any QM arm is read).

  --confirm   The frozen cluster-disjoint 25-fold design: 18 tabular arms + plain_qm.

The harness below is copied from notebook 37's own (cells 13-18), which is the established
pattern -- notebooks 35 and 37 each copied 34's. Differences, all deliberate:
  * `_in_design` is open on all 25 cells (37's returns `repeat in [3, 4]`).
  * No AID arm, so no winsorisation and no augmented population.
  * A `ridge` model class, and a `--descriptors-columns` chemprop arm generalised from notebook
    18's single-column local builders to many columns.
  * The OOF-wide assembly is NOT done here. Notebook 37's does it by starting from notebook 35's
    frozen 15-fold table and asserting every present arm is a column of it, which a new arm
    fires. Notebook 44 assembles from this run's own per-run oof_long/ files instead.

Usage:
  python scripts/44_run_qm_arms.py --sanity
  python scripts/44_run_qm_arms.py --tier1
  python scripts/44_run_qm_arms.py --confirm [--arms a,b,c] [--repeats 0,1]
"""
from __future__ import annotations

import os

# THREAD ENVIRONMENT. Notebook 34 established that restricting BLAS/OMP to one thread changes
# floating-point reduction order and so chemprop's training trajectory, breaking
# bit-reproducibility against notebook 05's stored scores: its sanity check failed at 0.0378
# against a 0.03 threshold purely because its setup had copied notebook 33's OMP_NUM_THREADS=1.
#
# OMP_NUM_THREADS is asserted unset, because nothing sets it implicitly.
# KMP_DUPLICATE_LIB_OK is only REPORTED, never asserted: notebook 40's Part 8 established that
# `import xgboost` SETS it, so an assertion placed after imports is unsatisfiable -- that is the
# exact defect that broke notebook 40 on a top-to-bottom re-run.
assert "OMP_NUM_THREADS" not in os.environ, (
    "OMP_NUM_THREADS is set. It changes chemprop's training trajectory and breaks "
    "bit-reproducibility against notebook 05's stored scores. Unset it and re-run.")
_KMP_AT_START = os.environ.get("KMP_DUPLICATE_LIB_OK")

import argparse  # noqa: E402
import hashlib  # noqa: E402
import inspect  # noqa: E402
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

import lightgbm as lgb  # noqa: E402
from lightgbm import LGBMRegressor  # noqa: E402
from sklearn.ensemble import RandomForestRegressor  # noqa: E402
from sklearn.linear_model import RidgeCV  # noqa: E402
from sklearn.pipeline import Pipeline  # noqa: E402
from sklearn.preprocessing import FunctionTransformer, StandardScaler  # noqa: E402
from xgboost import XGBRegressor  # noqa: E402

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
PIC50_COL = {iso: f"{iso}_pIC50_direct_inhibition" for iso in ISOFORMS}
BASE_COLS = ["Molecule_Name", "inchikey", "canonical_smiles"]

OUT = REPO_ROOT / "outputs" / "44_qm_molecular_block"
FEAT_DIR = OUT / "features"
PRED_DIR = OUT / "predictions"
SCORE_DIR = OUT / "scores"
OOF_DIR = OUT / "oof_long"
CHEMPROP_RUNS_DIR = OUT / "chemprop_runs"
T1 = OUT / "tier1"
LOG_DIR = REPO_ROOT / "logs" / "44_qm_molecular_block"

FOLDS_PATH = REPO_ROOT / "data" / "folds" / "cv_folds.csv"
BUTINA5X5_PATH = REPO_ROOT / "data" / "folds" / "cv_folds_butina_5x5.csv"
CURATED_PATH = REPO_ROOT / "data" / "processed" / "train_inhibition_curated.csv"
NB05_OUT = REPO_ROOT / "outputs" / "05_cv_comparison"
MANIFEST_PATH = NB05_OUT / "manifest.csv"

NB34_FOLD_SHA256 = "52f39b1a1cea9200b58dce0a061d4a21ce94d9d58ec77067cf6508b0f4b3e87e"
VAL_FRACTION = 0.15
EPOCHS, PATIENCE = 50, 5
CHEMELEON_ARCH_ARGS = ["--from-foundation", "CHEMELEON", "--multi-hot-atom-featurizer-mode", "V2"]
XGB_LGBM_EARLY_STOPPING_ROUNDS = 10
N_FOLDS = 5
SANITY_THRESHOLD = 0.03
RIDGE_ALPHAS = np.logspace(-2, 4, 13)
RIDGE_CLIP_PCT = (1.0, 99.0)
MODEL_CLASSES = ["rf", "lightgbm", "xgboost", "ridge"]
SEEDS_3 = [2684470948, 4091952314, 233227757]

TIER1_BLOCKS = ["RDKIT2D", "RDKIT2D_QM", "RDKIT2D_3D", "RDKIT2D_NEW3", "RDKIT2D_ELEC",
                "QM_FULL", "QM_ELEC", "QM_3D", "QM_NEW3", "QM_DFT3", "QM_GFN23",
                "RDKIT_MATCH18", "RDKIT_MATCH3_ELEC", "RDKIT_MATCH3_NEW3"]
# (block, model) on the frozen partition. `None` means "the selected M*".
CONFIRM_SPEC = [("RDKIT2D", None), ("RDKIT2D_QM", None), ("RDKIT2D_3D", None),
                ("RDKIT2D_NEW3", None), ("RDKIT2D_ELEC", None),
                ("QM_FULL", "rf"), ("QM_FULL", "lightgbm"), ("QM_FULL", "xgboost"),
                ("QM_FULL", "ridge"),
                ("RDKIT_MATCH18", None), ("RDKIT_MATCH3_ELEC", None),
                ("RDKIT_MATCH3_NEW3", None),
                ("QM_ELEC", None), ("QM_3D", None), ("QM_NEW3", None), ("QM_FULL_HMIN", None),
                ("QM_DFT3", None), ("QM_GFN23", None)]
CHEMPROP_DESC_BLOCK = "QM_FULL"       # plain_qm's extra molecule-level descriptors


def arm_name(block: str, model: str) -> str:
    return f"{block.lower()}__{model}"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ------------------------------------------------------------------ feature blocks
_FEATURE_CACHE: dict[str, tuple[pd.DataFrame, np.ndarray, list[str]]] = {}


def load_feature_matrix(block: str):
    """Reproduces scripts/run_5x5_cv_comparison.py's own load_feature_matrix (CSV branch)."""
    if block in _FEATURE_CACHE:
        return _FEATURE_CACHE[block]
    df = pd.read_csv(FEAT_DIR / f"{block}.csv")
    if "split" in df.columns:
        df = df[df["split"] == "train"].reset_index(drop=True)
    index_df = df[["Molecule_Name", "inchikey"]].reset_index(drop=True)
    drop_cols = [c for c in ["Molecule_Name", "inchikey", "split"] if c in df.columns]
    feat_cols = [c for c in df.columns if c not in drop_cols]
    X = df.drop(columns=drop_cols).to_numpy(dtype=np.float64)
    assert len(index_df) == 4905, f"{block}: expected 4905 training rows, got {len(index_df)}"
    assert np.isfinite(X).all(), f"{block}: non-finite feature value"
    _FEATURE_CACHE[block] = (index_df, X, feat_cols)
    return _FEATURE_CACHE[block]


def make_model(algo: str, seed: int, X_fit: np.ndarray):
    """One estimator. rf/xgboost/lightgbm are notebook 37's verbatim; ridge is new.

    Ridge needs the Ipc guard: that column reaches 5.06e+14 and becomes a ~75-sigma leverage
    point after standardisation, which made an early diagnostic ridge return R^2 as low as
    -6,400. Each feature is therefore clipped to its FIT-PORTION 1st/99th percentile before
    scaling -- fit-portion only, so no held-out value enters the bounds. Trees do not get the
    clip: they are invariant to monotone transforms, so it could only merge their tails for no
    benefit, and keeping them byte-identical to notebook 37's estimators matters more.
    """
    if algo == "rf":
        return RandomForestRegressor(random_state=seed)
    if algo == "xgboost":
        return XGBRegressor(early_stopping_rounds=XGB_LGBM_EARLY_STOPPING_ROUNDS,
                            random_state=seed)
    if algo == "lightgbm":
        return LGBMRegressor(random_state=seed, verbosity=-1)
    if algo == "ridge":
        lo = np.percentile(X_fit, RIDGE_CLIP_PCT[0], axis=0)
        hi = np.percentile(X_fit, RIDGE_CLIP_PCT[1], axis=0)
        return Pipeline([
            ("clip", FunctionTransformer(lambda A, lo=lo, hi=hi: np.clip(A, lo, hi))),
            ("scale", StandardScaler()),
            ("ridge", RidgeCV(alphas=RIDGE_ALPHAS)),
        ])
    raise ValueError(algo)


def run_tabular_one(block, algo, population, seed):
    """One tabular fit per endpoint. rf/ridge fit the pooled inner train+val; xgboost/lightgbm
    fit inner-train with inner-val as the early-stopping eval set -- notebook 37's own split of
    responsibilities, reproduced rather than reinvented."""
    feat_index, X, _ = load_feature_matrix(block)
    pop_pos = feat_index.reset_index().merge(population, on=["Molecule_Name", "inchikey"],
                                            how="left")
    assert pop_pos["screen_split"].notna().all(), f"{block}: feature index did not align 1:1"

    pool_mask = pop_pos["screen_split"].isin(["screen_inner_train", "screen_inner_val"])
    inner_train_mask = pop_pos["screen_split"] == "screen_inner_train"
    inner_val_mask = pop_pos["screen_split"] == "screen_inner_val"
    test_sel = pop_pos["screen_split"] == "screen_test"
    test_positions = pop_pos.loc[test_sel, "index"].to_numpy()

    pred_df = pop_pos.loc[test_sel, ["Molecule_Name", "inchikey"]].reset_index(drop=True)
    for endpoint in REGRESSION_ENDPOINTS:
        has_label = pop_pos[endpoint].notna()
        if algo in ("rf", "ridge"):
            idx = pop_pos.loc[pool_mask & has_label, "index"].to_numpy()
            y = pop_pos.loc[pool_mask & has_label, endpoint].to_numpy()
            model = make_model(algo, seed, X[idx])
            model.fit(X[idx], y)
        else:
            tr_idx = pop_pos.loc[inner_train_mask & has_label, "index"].to_numpy()
            va_idx = pop_pos.loc[inner_val_mask & has_label, "index"].to_numpy()
            y_tr = pop_pos.loc[inner_train_mask & has_label, endpoint].to_numpy()
            y_va = pop_pos.loc[inner_val_mask & has_label, endpoint].to_numpy()
            model = make_model(algo, seed, X[tr_idx])
            if algo == "xgboost":
                model.fit(X[tr_idx], y_tr, eval_set=[(X[va_idx], y_va)], verbose=False)
            else:
                model.fit(X[tr_idx], y_tr, eval_set=[(X[va_idx], y_va)],
                          callbacks=[lgb.early_stopping(
                              stopping_rounds=XGB_LGBM_EARLY_STOPPING_ROUNDS, verbose=False)])
        pred_df[endpoint] = np.asarray(model.predict(X[test_positions])).reshape(-1)
    return pred_df


# ------------------------------------------------------------------ chemprop with descriptors
def build_training_csv_multi(population, target_cols, descriptor_cols, out_path, logger):
    """Notebook 18's build_training_csv_with_charge, generalised from one descriptor column to
    many. src/chemprop_screen.py's own builder hardcodes its column list and has no hook, and is
    deliberately left unmodified -- notebook 18's own decision."""
    pooled = population[population["screen_split"].isin(
        ["screen_inner_train", "screen_inner_val"])].copy()
    split_map = {"screen_inner_train": "train", "screen_inner_val": "val"}
    pooled["chemprop_split"] = pooled["screen_split"].map(split_map)
    assert pooled["chemprop_split"].notna().all()
    cols = ["canonical_smiles", *target_cols, *descriptor_cols, "chemprop_split"]
    pooled[cols].to_csv(out_path, index=False)
    logger.info(f"wrote {out_path} ({len(pooled)} rows, {len(descriptor_cols)} extra descriptors, "
                f"chemprop_split {pooled['chemprop_split'].value_counts().to_dict()})")
    return pooled


def build_predict_csv_multi(population, descriptor_cols, run_dir, logger):
    run_dir.mkdir(parents=True, exist_ok=True)
    test_df = population.loc[population["screen_split"] == "screen_test",
                             BASE_COLS + list(descriptor_cols)].copy()
    out_path = run_dir / "predict_input.csv"
    test_df.to_csv(out_path, index=False)
    logger.info(f"wrote {out_path} ({len(test_df)} screen_test compounds)")
    return out_path


def run_chemprop_predict_multi(model_dir, predict_csv, descriptor_cols, output_csv, logger):
    argv = ["chemprop", "predict", "-i", str(predict_csv), "-s", "canonical_smiles",
            "--descriptors-columns", *descriptor_cols,
            "--model-paths", str(model_dir), "--accelerator", "cpu", "--devices", "1",
            "-o", str(output_csv)]
    run_subprocess_streamed(argv, logger)


# ------------------------------------------------------------------ shared state
class Ctx:
    """Everything the run functions need, assembled once."""

    def __init__(self, mode):
        for d in [OUT, PRED_DIR, SCORE_DIR, OOF_DIR, CHEMPROP_RUNS_DIR, T1, LOG_DIR]:
            d.mkdir(parents=True, exist_ok=True)
        self.mode = mode
        self.curated = pd.read_csv(CURATED_PATH)
        assert len(self.curated) == 4905
        self.curated_ik = set(self.curated["inchikey"].dropna())

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
        # Cross-check the mapping against what notebooks 34/35/37 actually used.
        for nb in ["34_butina_5x5", "35_butina_repeat3", "37_butina_5x5_complete"]:
            p = REPO_ROOT / "outputs" / nb / "run_scores.csv"
            if not p.exists():
                continue
            rs = pd.read_csv(p)[["repeat", "fold", "seed"]].drop_duplicates()
            bad = [(int(r.repeat), int(r.fold), int(r.seed)) for r in rs.itertuples()
                   if self.fold_seed[(int(r.repeat), int(r.fold))] != int(r.seed)]
            assert not bad, f"seed mapping differs from notebook {nb}'s: {bad}"

        self.pop_base = self.curated[BASE_COLS + list(REGRESSION_ENDPOINTS)].copy()
        self.scores_path = (T1 if mode == "tier1" else OUT) / "run_scores.csv"
        self.timings_path = (T1 if mode == "tier1" else OUT) / "run_timings.csv"
        self.pred_dir = (T1 / "predictions") if mode == "tier1" else PRED_DIR
        self.score_dir = (T1 / "scores") if mode == "tier1" else SCORE_DIR
        self.oof_dir = (T1 / "oof_long") if mode == "tier1" else OOF_DIR
        for d in [self.pred_dir, self.score_dir, self.oof_dir]:
            d.mkdir(parents=True, exist_ok=True)

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
    """Notebook 37's upsert, verbatim. Keyed so a new arm's rows land beside the old ones."""
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


def tag_of(arm, repeat, fold, seed, mode):
    return f"{arm}__r{repeat}f{fold}" + (f"s{seed}" if mode == "tier1" else "")


def cached(ctx, arm, repeat, fold, seed, expected_names):
    """Notebook 37's content-based cache check. Existence is never enough."""
    p = ctx.pred_dir / f"{tag_of(arm, repeat, fold, seed, ctx.mode)}.csv"
    if not p.exists():
        return None
    try:
        df = pd.read_csv(p)
    except Exception:
        return None
    if len(df) != len(expected_names) or set(df["Molecule_Name"]) != expected_names:
        print(f"  [cache MISMATCH] {p.name}: content check failed -- will retrain")
        return None
    if df[list(REGRESSION_ENDPOINTS)].isna().to_numpy().any():
        print(f"  [cache MISMATCH] {p.name}: NaN in a predicted column -- will retrain")
        return None
    return df


def score_and_append(ctx, arm, repeat, fold, seed, pred_df):
    gt = ctx.curated[ctx.curated["inchikey"].isin(pred_df["inchikey"])].copy()
    with per_fold_bootstrap_seed(seed):
        scored = score_activity_predictions(pred_df, gt, list(REGRESSION_ENDPOINTS))
        scored = add_macro_endpoint(scored, list(REGRESSION_ENDPOINTS), ACTIVITY_METRICS)
    scored.to_csv(ctx.score_dir / f"{tag_of(arm, repeat, fold, seed, ctx.mode)}.csv", index=False)
    pe = scored.groupby("Endpoint")[METRIC_NAMES].mean()
    rows = []
    for endpoint in pe.index:
        row = {"repeat": repeat, "fold": fold, "arm": arm, "seed": seed, "Endpoint": endpoint}
        row.update(pe.loc[endpoint].to_dict())
        rows.append(row)
    _append_rows(ctx.scores_path, rows, ["repeat", "fold", "arm", "seed", "Endpoint"])
    return pe


def save_oof_long(ctx, arm, repeat, fold, seed, pred_df):
    frames = []
    for iso in ISOFORMS:
        ep = PIC50_COL[iso]
        m = pred_df[["inchikey", ep]].rename(columns={ep: "y_pred"}).copy()
        m["isoform"] = iso
        frames.append(m[["inchikey", "isoform", "y_pred"]])
    out = pd.concat(frames, ignore_index=True)
    out.insert(0, "seed", seed)
    out.insert(0, "fold", fold)
    out.insert(0, "repeat", repeat)
    out.insert(0, "arm", arm)
    out.to_csv(ctx.oof_dir / f"{tag_of(arm, repeat, fold, seed, ctx.mode)}.csv", index=False)


def record_timing(ctx, arm, repeat, fold, seed, wall, last_epoch=None):
    """Notebook 37's anti-clobber rule: a cache hit must not overwrite a completed run's timing.
    Notebook 35 lost its whole timing file exactly that way."""
    row = {"repeat": repeat, "fold": fold, "arm": arm, "seed": seed, "wall_min": wall,
           "last_epoch": last_epoch,
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


def run_tabular_unit(ctx, block, algo, partition, repeat, fold, seed):
    arm = arm_name(block, algo)
    pop = ctx.population(partition, (ctx.butina_cols[repeat] if partition == "butina"
                                     else f"repeat_{repeat}"), fold, seed)
    expected = set(pop.loc[pop["screen_split"] == "screen_test", "Molecule_Name"])
    hit = cached(ctx, arm, repeat, fold, seed, expected)
    if hit is not None:
        score_and_append(ctx, arm, repeat, fold, seed, hit)
        save_oof_long(ctx, arm, repeat, fold, seed, hit)
        record_timing(ctx, arm, repeat, fold, seed, 0.0)
        print(f"  [cache hit] {tag_of(arm, repeat, fold, seed, ctx.mode)}")
        return hit
    t0 = time.time()
    pred_df = run_tabular_one(block, algo, pop, seed)
    wall = (time.time() - t0) / 60.0
    assert set(pred_df["Molecule_Name"]) == expected
    pred_df.to_csv(ctx.pred_dir / f"{tag_of(arm, repeat, fold, seed, ctx.mode)}.csv", index=False)
    score_and_append(ctx, arm, repeat, fold, seed, pred_df)
    save_oof_long(ctx, arm, repeat, fold, seed, pred_df)
    record_timing(ctx, arm, repeat, fold, seed, wall)
    print(f"  {tag_of(arm, repeat, fold, seed, ctx.mode)}: {wall:.2f} min")
    return pred_df


def run_plain_qm_unit(ctx, partition, repeat, fold, seed, arm="plain_qm"):
    """chemprop_chemeleoninit plus the QM block as extra molecule-level descriptors."""
    _, _, desc_cols = load_feature_matrix(CHEMPROP_DESC_BLOCK)
    qm_tbl = pd.read_csv(FEAT_DIR / f"{CHEMPROP_DESC_BLOCK}.csv")
    pop = ctx.population(partition, (ctx.butina_cols[repeat] if partition == "butina"
                                     else f"repeat_{repeat}"), fold, seed)
    n_before = len(pop)
    pop = pop.merge(qm_tbl[["inchikey"] + desc_cols], on="inchikey", how="left")
    assert len(pop) == n_before and pop[desc_cols].notna().all().all()
    expected = set(pop.loc[pop["screen_split"] == "screen_test", "Molecule_Name"])
    tag = tag_of(arm, repeat, fold, seed, ctx.mode)
    hit = cached(ctx, arm, repeat, fold, seed, expected)
    if hit is not None:
        score_and_append(ctx, arm, repeat, fold, seed, hit)
        save_oof_long(ctx, arm, repeat, fold, seed, hit)
        record_timing(ctx, arm, repeat, fold, seed, 0.0)
        print(f"  [cache hit] {tag}")
        return hit
    run_dir = CHEMPROP_RUNS_DIR / tag
    run_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(LOG_DIR / f"{tag}.log", f"nb44_{tag}")
    logger.info(f"=== {tag}: {len(desc_cols)} extra descriptors: {desc_cols} ===")
    t0 = time.time()
    train_csv = run_dir / "train_input.csv"
    build_training_csv_multi(pop, list(REGRESSION_ENDPOINTS), desc_cols, train_csv, logger)
    predict_csv = build_predict_csv_multi(pop, desc_cols, run_dir, logger)
    arch = CHEMELEON_ARCH_ARGS + ["--loss-function", "mse",
                                  "--descriptors-columns", *desc_cols]
    run_chemprop_train(train_csv, list(REGRESSION_ENDPOINTS), run_dir, logger, arch,
                       EPOCHS, PATIENCE, seed)
    raw = run_dir / "raw_predictions.csv"
    run_chemprop_predict_multi(run_dir / "model_0", predict_csv, desc_cols, raw, logger)
    verify_predictions(raw, expected, list(REGRESSION_ENDPOINTS), logger)
    wall = (time.time() - t0) / 60.0
    pred_df = pd.read_csv(raw)
    assert "inchikey" in pred_df.columns
    pred_df.to_csv(ctx.pred_dir / f"{tag}.csv", index=False)
    score_and_append(ctx, arm, repeat, fold, seed, pred_df)
    save_oof_long(ctx, arm, repeat, fold, seed, pred_df)
    ep = last_epoch(run_dir)
    record_timing(ctx, arm, repeat, fold, seed, wall, ep)
    print(f"  {tag}: {wall:.2f} min (last_epoch={ep})")
    return pred_df


def last_epoch(run_dir: Path):
    m = run_dir / "model_0" / "trainer_logs" / "version_0" / "metrics.csv"
    if not m.exists():
        return None
    eps = [int(r.split(",")[0]) for r in m.read_text().splitlines()[1:]
           if r and r.split(",")[0].isdigit()]
    return max(eps) if eps else None


# ------------------------------------------------------------------ modes
def mode_sanity(ctx):
    from src.chemprop_screen import load_screen_population
    print("=" * 78)
    print("SANITY GATE -- reproduce notebook 05's stored chemprop_chemeleoninit repeat0_fold0")
    print("=" * 78)
    tag = "sanity_plain__random_r0f0"
    seed = ctx.fold_seed[(0, 0)]
    logger = setup_logging(LOG_DIR / f"{tag}.log", f"nb44_{tag}")
    pop = load_screen_population(FOLDS_PATH, CURATED_PATH, "repeat_0", 0, VAL_FRACTION,
                                 seed, logger)
    expected = set(pop.loc[pop["screen_split"] == "screen_test", "Molecule_Name"])
    pred_path = PRED_DIR / f"{tag}.csv"
    if pred_path.exists() and set(pd.read_csv(pred_path)["Molecule_Name"]) == expected:
        print(f"[cache hit] {tag}")
        pred = pd.read_csv(pred_path)
    else:
        run_dir = CHEMPROP_RUNS_DIR / tag
        run_dir.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        train_csv = run_dir / "train_input.csv"
        build_training_csv(pop, list(REGRESSION_ENDPOINTS), train_csv, logger,
                           require_all_targets=False)
        predict_csv = build_predict_csv(pop, run_dir, logger)
        run_chemprop_train(train_csv, list(REGRESSION_ENDPOINTS), run_dir, logger,
                           CHEMELEON_ARCH_ARGS + ["--loss-function", "mse"],
                           EPOCHS, PATIENCE, seed)
        raw = run_dir / "raw_predictions.csv"
        run_chemprop_predict(run_dir / "model_0", predict_csv, raw, logger)
        verify_predictions(raw, expected, list(REGRESSION_ENDPOINTS), logger)
        pred = pd.read_csv(raw)
        pred.to_csv(pred_path, index=False)
        print(f"sanity run wall time: {(time.time()-t0)/60:.2f} min")

    gt = ctx.curated[ctx.curated["inchikey"].isin(pred["inchikey"])].copy()
    with per_fold_bootstrap_seed(seed):
        scored = score_activity_predictions(pred, gt, list(REGRESSION_ENDPOINTS))
    mine = scored.groupby("Endpoint")["ST-RAE"].mean()
    theirs = pd.read_csv(
        NB05_OUT / "scores" / "chemprop_chemeleoninit__repeat0_fold0.csv"
    ).groupby("Endpoint")["ST-RAE"].mean()
    print(f"\n{'endpoint':42s} {'this run':>12s} {'nb05':>12s} {'abs diff':>12s}")
    max_diff = 0.0
    rows = []
    for ep in sorted(set(mine.index) & set(theirs.index)):
        d = abs(float(mine[ep]) - float(theirs[ep]))
        max_diff = max(max_diff, d)
        print(f"{ep:42s} {mine[ep]:12.10f} {theirs[ep]:12.10f} {d:12.2e}")
        rows.append({"Endpoint": ep, "this_run": float(mine[ep]), "nb05": float(theirs[ep]),
                     "abs_diff": d})
    pd.DataFrame(rows).to_csv(OUT / "sanity_check.csv", index=False)
    print(f"\nmax abs ST-RAE difference: {max_diff:.10f}   threshold {SANITY_THRESHOLD}")
    (OUT / "sanity_verdict.json").write_text(json.dumps(
        {"max_abs_strae_diff": max_diff, "threshold": SANITY_THRESHOLD,
         "pass": bool(max_diff <= SANITY_THRESHOLD),
         "at_utc": datetime.now(timezone.utc).isoformat()}, indent=2))
    if max_diff > SANITY_THRESHOLD:
        print("*** SANITY GATE FAILED -- nothing downstream is interpretable. STOPPING. ***")
        return 1
    print("PASS -- the pipeline reproduces notebook 05's own stored numbers.")
    return 0


def mode_tier1(ctx):
    print("=" * 78)
    print("TIER-1 SCREEN -- RANDOM partition repeat_0/fold 0, 3 seeds")
    print("the frozen Butina 25-fold stays UNSEEN until --confirm")
    print("=" * 78)
    t_start = time.time()
    print(f"\nphase A: {len(TIER1_BLOCKS)} blocks x {len(MODEL_CLASSES)} models x "
          f"{len(SEEDS_3)} seeds = {len(TIER1_BLOCKS)*len(MODEL_CLASSES)*len(SEEDS_3)} runs")
    for block in TIER1_BLOCKS:
        for algo in MODEL_CLASSES:
            for seed in SEEDS_3:
                run_tabular_unit(ctx, block, algo, "random", 0, 0, seed)

    # M* is selected on the CONTROL block, before any QM arm is read. Pre-registered rule.
    rs = pd.read_csv(ctx.scores_path)
    ctrl = rs[(rs["arm"].str.startswith("rdkit2d__")) & (rs["Endpoint"] == "MA")]
    by_model = ctrl.assign(model=ctrl["arm"].str.split("__").str[-1]).groupby(
        "model")["ST-RAE"].mean().sort_values()
    m_star = str(by_model.index[0])
    print(f"\nmacro ST-RAE on the CONTROL block (rdkit2d), mean over {len(SEEDS_3)} seeds:")
    for mdl, v in by_model.items():
        print(f"  {mdl:10s} {v:.4f}" + ("   <-- M*" if mdl == m_star else ""))
    (OUT / "model_selection.json").write_text(json.dumps(
        {"m_star": m_star, "rule": ("lowest macro ST-RAE on the CONTROL block RDKIT2D at "
                                    "tier-1, chosen before any QM arm was read"),
         "control_macro_strae": {k: float(v) for k, v in by_model.items()},
         "seeds": SEEDS_3, "partition": "random repeat_0/fold 0",
         "at_utc": datetime.now(timezone.utc).isoformat()}, indent=2))
    print(f"\nM* = {m_star}  -> outputs/44_qm_molecular_block/model_selection.json")

    n_draws = len([p for p in FEAT_DIR.glob("RDKIT_RAND18_*.csv")])
    print(f"\nphase B: the descriptive null band, {n_draws} draws x {len(SEEDS_3)} seeds at M*")
    for p in sorted(FEAT_DIR.glob("RDKIT_RAND18_*.csv")):
        for seed in SEEDS_3:
            run_tabular_unit(ctx, p.stem, m_star, "random", 0, 0, seed)

    print(f"\nphase C: plain_qm at tier-1 is DELIBERATELY SKIPPED -- revision 1 item 3 ungated "
          f"it, so the 25-fold run supersedes a screen of it.")
    print(f"\ntier-1 complete in {(time.time()-t_start)/60:.1f} min")
    return 0


def mode_confirm(ctx, arms_filter, repeats):
    sel = OUT / "model_selection.json"
    assert sel.exists(), "run --tier1 first: M* must be selected on the control block"
    m_star = json.loads(sel.read_text())["m_star"]
    print("=" * 78)
    print(f"CONFIRM -- frozen cluster-disjoint 25-fold design, M* = {m_star}")
    print("=" * 78)
    spec = [(b, mdl or m_star) for b, mdl in CONFIRM_SPEC]
    # De-duplicate: QM_FULL is listed once per model class, one of which is M*.
    seen, uniq = set(), []
    for b, mdl in spec:
        if (b, mdl) not in seen:
            seen.add((b, mdl))
            uniq.append((b, mdl))
    if arms_filter:
        want = set(arms_filter)
        uniq = [(b, mdl) for b, mdl in uniq if arm_name(b, mdl) in want]
    print(f"{len(uniq)} tabular arms + plain_qm, repeats {repeats} x {N_FOLDS} folds")
    for b, mdl in uniq:
        print(f"  {arm_name(b, mdl)}")
    t_start = time.time()
    for repeat in repeats:
        for fold in range(N_FOLDS):
            seed = ctx.fold_seed[(repeat, fold)]
            print(f"\n--- butina repeat {repeat} fold {fold} (seed {seed}) ---")
            for b, mdl in uniq:
                run_tabular_unit(ctx, b, mdl, "butina", repeat, fold, seed)
            if not arms_filter or "plain_qm" in set(arms_filter):
                run_plain_qm_unit(ctx, "butina", repeat, fold, seed)
    print(f"\nconfirm complete in {(time.time()-t_start)/60:.1f} min")
    return 0


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--sanity", action="store_true")
    g.add_argument("--tier1", action="store_true")
    g.add_argument("--confirm", action="store_true")
    ap.add_argument("--arms", default=None, help="comma-separated arm names, --confirm only")
    ap.add_argument("--repeats", default="0,1,2,3,4", help="comma-separated, --confirm only")
    a = ap.parse_args()

    print(f"KMP_DUPLICATE_LIB_OK at start: {_KMP_AT_START!r}  (REPORTED, never asserted -- "
          f"`import xgboost` sets it; notebook 40 Part 8)")
    print(f"KMP_DUPLICATE_LIB_OK after imports: "
          f"{os.environ.get('KMP_DUPLICATE_LIB_OK')!r}")
    env_bin = str(Path(sys.executable).parent)
    parts = os.environ.get("PATH", "").split(os.pathsep)
    if env_bin not in parts:
        os.environ["PATH"] = os.pathsep.join([env_bin, *parts])
    print(f"chemprop: {subprocess.run(['which','chemprop'],capture_output=True,text=True).stdout.strip()}")
    assert (OUT / "prereg.json").exists(), (
        "prereg.json is absent -- run scripts/44_build_qm_feature_blocks.py first. The "
        "pre-registration must be frozen before any training run.")

    mode = "sanity" if a.sanity else ("tier1" if a.tier1 else "confirm")
    ctx = Ctx(mode)
    if a.sanity:
        return mode_sanity(ctx)
    if a.tier1:
        return mode_tier1(ctx)
    return mode_confirm(ctx, a.arms.split(",") if a.arms else None,
                        [int(x) for x in a.repeats.split(",")])


if __name__ == "__main__":
    sys.exit(main())
