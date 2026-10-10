"""Training driver for notebook 46's QM arm: tier 1 (choose M*) and tier 2 (the formal 25 folds).

SELF-CONTAINED BY CHOICE, not by accident. `scripts/44_run_qm_arms.py` does nearly all of this
already, and importing it would have been shorter -- but every path it writes to is a module-level
constant pointing inside `outputs/44_qm_molecular_block/`, which is protected. Monkeypatching
those globals risks missing one and writing into notebook 44's tree, which is precisely the
scope violation the run is checked for afterwards. The functions below are therefore copied from
it verbatim except where a comment says otherwise, following the same precedent by which
notebooks 35 and 37 each copied notebook 34's harness.

WHAT IS NEW HERE, and it is only two things:
  * `load_feature_matrix` tolerates NaN for the blocks declared in the pre-registration, instead
    of asserting finiteness unconditionally.
  * `run_tabular_one` fits a per-fold median fill for those columns on that fold's TRAINING
    portion only. A whole-population fill would put held-out information into every fit, which
    is why the feature CSV ships the column with NaN rather than pre-filled.
"""

from __future__ import annotations

import os

# Notebook 34 established that restricting BLAS/OMP to one thread changes chemprop's training
# trajectory; nothing here trains chemprop, but the assertion is kept so the environment this
# project measures everything under stays the same one.
assert os.environ.get("OMP_NUM_THREADS") is None, "OMP_NUM_THREADS must be unset"
_KMP_AT_START = os.environ.get("KMP_DUPLICATE_LIB_OK")

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from datetime import datetime, timezone  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from lightgbm import LGBMRegressor  # noqa: E402
from sklearn.ensemble import RandomForestRegressor  # noqa: E402
from sklearn.linear_model import RidgeCV  # noqa: E402
from sklearn.pipeline import Pipeline  # noqa: E402
from sklearn.preprocessing import FunctionTransformer, StandardScaler  # noqa: E402
from xgboost import XGBRegressor  # noqa: E402

from src.cv_bootstrap import per_fold_bootstrap_seed  # noqa: E402
from src.features import assign_screen_split  # noqa: E402
from src.vendor.openadmet_eval.config import ACTIVITY_METRICS, REGRESSION_ENDPOINTS  # noqa: E402
from src.vendor.openadmet_eval.evaluate_predictions import (  # noqa: E402
    add_macro_endpoint, score_activity_predictions,
)

METRIC_NAMES = [name for name, _ in ACTIVITY_METRICS]
ISOFORMS = [e.split("_")[0] for e in REGRESSION_ENDPOINTS]
PIC50_COL = {e.split("_")[0]: e for e in REGRESSION_ENDPOINTS}
BASE_COLS = ["Molecule_Name", "inchikey", "canonical_smiles"]

OUT = REPO_ROOT / "outputs" / "46_qm_arm"
FEAT_DIR = OUT / "features"
T1 = OUT / "tier1"
PRED_DIR = OUT / "predictions"
SCORE_DIR = OUT / "scores"
OOF_DIR = OUT / "oof_long"

BUTINA5X5_PATH = REPO_ROOT / "data" / "folds" / "cv_folds_butina_5x5.csv"
FOLDS_PATH = REPO_ROOT / "data" / "folds" / "cv_folds.csv"
CURATED_PATH = REPO_ROOT / "data" / "processed" / "train_inhibition_curated.csv"
MANIFEST_PATH = REPO_ROOT / "outputs" / "05_cv_comparison" / "manifest.csv"
NB34_FOLD_SHA256 = "52f39b1a1cea9200b58dce0a061d4a21ce94d9d58ec77067cf6508b0f4b3e87e"

VAL_FRACTION = 0.15
RIDGE_ALPHAS = np.logspace(-2, 4, 13)
RIDGE_CLIP_PCT = (1.0, 99.0)
MODEL_CLASSES = ["rf", "lightgbm", "xgboost", "ridge"]
SEEDS_3 = [2684470948, 4091952314, 233227757]
XGB_LGBM_EARLY_STOPPING_ROUNDS = 10

BLOCKS = ["FREE", "FREE_QM13", "FREE_ATOM4", "FREE_QM13_ATOM4", "FREE_ATOM4_HAL", "FREE_CONTROL"]
BASELINE_BLOCK = "FREE"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def arm_name(block: str, model: str) -> str:
    return f"{block.lower()}__{model}"


def prereg() -> dict:
    return json.loads((OUT / "prereg.json").read_text())


# ------------------------------------------------------------------ feature blocks
_FEATURE_CACHE: dict = {}


def load_feature_matrix(block: str):
    """Notebook 44's loader, with one change: blocks the pre-registration declares may carry NaN
    are allowed to, because their fill is fitted per fold by `run_tabular_one`."""
    if block in _FEATURE_CACHE:
        return _FEATURE_CACHE[block]
    may_nan = set(prereg()["nan_fill_policy"]["blocks"])
    df = pd.read_csv(FEAT_DIR / f"{block}.csv")
    index_df = df[["Molecule_Name", "inchikey"]].reset_index(drop=True)
    feat_cols = [c for c in df.columns if c not in ("Molecule_Name", "inchikey")]
    X = df[feat_cols].to_numpy(dtype=np.float64)
    assert len(index_df) == 4905, f"{block}: expected 4905 training rows, got {len(index_df)}"
    if block in may_nan:
        assert np.isfinite(X[~np.isnan(X)]).all(), f"{block}: non-finite value that is not NaN"
    else:
        assert np.isfinite(X).all(), f"{block}: non-finite feature value"
    _FEATURE_CACHE[block] = (index_df, X, feat_cols)
    return _FEATURE_CACHE[block]


def make_model(algo: str, seed: int, X_fit: np.ndarray):
    """Notebook 44's `make_model`, verbatim. rf/xgboost/lightgbm are notebook 37's estimators;
    ridge carries the Ipc guard -- that column reaches 5.06e+14 and becomes a ~75-sigma leverage
    point after standardisation, so each feature is clipped to its FIT-PORTION 1st/99th
    percentile first. Fit-portion only, so no held-out value enters the bounds."""
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
    """Notebook 44's `run_tabular_one`, with the per-fold NaN fill added.

    THE FILL. Computed once per fold from the pooled training portion (inner train + inner val),
    identical for every endpoint and model, and applied to the held-out portion too. Using the
    training portion only is what keeps it leak-free; using the pooled portion rather than a
    per-endpoint subset keeps one fold's feature matrix single-valued, which is what makes the
    six methods comparable on that fold.
    """
    feat_index, X, feat_cols = load_feature_matrix(block)
    X = X.copy()
    pop_pos = feat_index.reset_index().merge(population, on=["Molecule_Name", "inchikey"],
                                             how="left")
    assert pop_pos["screen_split"].notna().all(), f"{block}: feature index did not align 1:1"

    pool_mask = pop_pos["screen_split"].isin(["screen_inner_train", "screen_inner_val"])
    inner_train_mask = pop_pos["screen_split"] == "screen_inner_train"
    inner_val_mask = pop_pos["screen_split"] == "screen_inner_val"
    test_sel = pop_pos["screen_split"] == "screen_test"
    test_positions = pop_pos.loc[test_sel, "index"].to_numpy()

    nan_cols = np.where(np.isnan(X).any(axis=0))[0]
    if len(nan_cols):
        fit_pos = pop_pos.loc[pool_mask, "index"].to_numpy()
        for j in nan_cols:
            col = X[fit_pos, j]
            present = col[~np.isnan(col)]
            assert len(present) >= 30, f"{block}/{feat_cols[j]}: only {len(present)} present in fit"
            X[np.isnan(X[:, j]), j] = float(np.median(present))
    assert np.isfinite(X).all(), f"{block}: NaN survived the fill"

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
                import lightgbm as lgb
                model.fit(X[tr_idx], y_tr, eval_set=[(X[va_idx], y_va)],
                          callbacks=[lgb.early_stopping(XGB_LGBM_EARLY_STOPPING_ROUNDS,
                                                        verbose=False)])
        pred_df[endpoint] = np.asarray(model.predict(X[test_positions])).reshape(-1)
    return pred_df


# ------------------------------------------------------------------ shared state
class Ctx:
    """Notebook 44's Ctx, with its paths repointed at this notebook's own output tree."""

    def __init__(self, mode):
        for d in [OUT, PRED_DIR, SCORE_DIR, OOF_DIR, T1]:
            d.mkdir(parents=True, exist_ok=True)
        self.mode = mode
        self.curated = pd.read_csv(CURATED_PATH)
        assert len(self.curated) == 4905

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
    """Content-based, never existence -- notebook 37's rule."""
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
    for col, val in [("seed", seed), ("fold", fold), ("repeat", repeat), ("arm", arm)]:
        out.insert(0, col, val)
    out.to_csv(ctx.oof_dir / f"{tag_of(arm, repeat, fold, seed, ctx.mode)}.csv", index=False)


def record_timing(ctx, arm, repeat, fold, seed, wall):
    """Notebook 37's anti-clobber rule: a cache hit must not overwrite a completed run's timing.
    Notebook 35 lost its whole timing file exactly that way."""
    row = {"repeat": repeat, "fold": fold, "arm": arm, "seed": seed, "wall_min": wall}
    key = ["repeat", "fold", "arm", "seed"]
    already = False
    if ctx.timings_path.exists():
        t = pd.read_csv(ctx.timings_path)
        already = bool(((t["repeat"] == repeat) & (t["fold"] == fold) & (t["arm"] == arm)
                        & (t["seed"] == seed)).any())
    if wall < 0.02 and already:
        return
    _append_rows(ctx.timings_path, [row], key)


def run_unit(ctx, block, algo, partition, repeat, fold, seed):
    arm = arm_name(block, algo)
    repeat_col = ctx.butina_cols[repeat] if partition == "butina" else f"repeat_{repeat}"
    pop = ctx.population(partition, repeat_col, fold, seed)
    expected = set(pop.loc[pop["screen_split"] == "screen_test", "Molecule_Name"])
    hit = cached(ctx, arm, repeat, fold, seed, expected)
    t0 = time.time()
    if hit is not None:
        pred = hit
    else:
        pred = run_tabular_one(block, algo, pop, seed)
        pred.to_csv(ctx.pred_dir / f"{tag_of(arm, repeat, fold, seed, ctx.mode)}.csv", index=False)
    wall = (time.time() - t0) / 60
    score_and_append(ctx, arm, repeat, fold, seed, pred)
    save_oof_long(ctx, arm, repeat, fold, seed, pred)
    record_timing(ctx, arm, repeat, fold, seed, wall)
    return wall


# ------------------------------------------------------------------ tiers
def tier1(ctx):
    """6 blocks x 4 models x 3 seeds, random partition repeat_0 / fold 0. No p-values are
    produced here and none should be reported: this tier exists only to choose M*."""
    print(f"TIER 1 -- random partition repeat_0/fold 0, {len(BLOCKS)} blocks x "
          f"{len(MODEL_CLASSES)} models x {len(SEEDS_3)} seeds")
    n = 0
    for block in BLOCKS:
        for algo in MODEL_CLASSES:
            for seed in SEEDS_3:
                w = run_unit(ctx, block, algo, "random", 0, 0, seed)
                n += 1
                print(f"  [{n:3d}] {arm_name(block, algo):28s} seed {seed}  {w:5.2f} min")
    return n


def select_mstar(ctx):
    """M* is chosen on the BASELINE block alone, so the model cannot be selected for liking the
    QM blocks. Rule frozen in prereg.json before tier 1 ran."""
    rs = pd.read_csv(ctx.scores_path)
    iso = rs[~rs.Endpoint.str.contains("MA|macro", case=False, na=False)]
    base = iso[iso.arm.str.startswith(BASELINE_BLOCK.lower() + "__")]
    per = base.groupby(["arm", "seed"])["ST-RAE"].mean().groupby("arm").mean().sort_values()
    mstar = per.index[0].split("__")[1]
    print("\nM* selection, on the BASELINE block alone (mean macro ST-RAE over 3 seeds):")
    print(per.round(4).to_string())
    print(f"  M* = {mstar}")
    return mstar, per


def confirm(ctx, mstar, repeats):
    print(f"\nTIER 2 -- frozen Butina partition, {len(BLOCKS)} methods x 25 folds, M* = {mstar}")
    n = 0
    for block in BLOCKS:
        for repeat in repeats:
            for fold in range(5):
                seed = ctx.fold_seed[(repeat, fold)]
                w = run_unit(ctx, block, mstar, "butina", repeat, fold, seed)
                n += 1
                print(f"  [{n:3d}] {arm_name(block, mstar):28s} r{repeat}f{fold}  {w:5.2f} min")
    return n


def main() -> int:
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--tier1", action="store_true")
    g.add_argument("--confirm", action="store_true")
    ap.add_argument("--repeats", default="0,1,2,3,4")
    args = ap.parse_args()

    t0 = time.time()
    if args.tier1:
        ctx = Ctx("tier1")
        n = tier1(ctx)
        mstar, per = select_mstar(ctx)
        (T1 / "mstar.json").write_text(json.dumps({
            "selected_at_utc": datetime.now(timezone.utc).isoformat(),
            "mstar": mstar, "rule": prereg()["tier1"]["selection_rule"],
            "baseline_block": BASELINE_BLOCK,
            "mean_macro_strae_by_model": {k: float(v) for k, v in per.items()},
            "n_runs": n,
        }, indent=2))
        print(f"\nwrote {T1 / 'mstar.json'}")
    else:
        ctx = Ctx("confirm")
        mstar = json.loads((T1 / "mstar.json").read_text())["mstar"]
        n = confirm(ctx, mstar, [int(r) for r in args.repeats.split(",")])
    print(f"\n{n} units, {(time.time() - t0) / 60:.1f} min total")
    print(f"KMP_DUPLICATE_LIB_OK at start: {_KMP_AT_START!r} (reported, not asserted -- "
          f"`import xgboost` sets it; notebook 40 Part 8)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
