"""CYP2D6 sample-weighting and hyperparameter-tuning levers (Roadmap Step 2,
notebook 07).

BASELINE FOR THIS STEP: 06's CV-residual-excluded (criterion 1) results for
`chemeleon__rf` and `chemeleon__lightgbm` -- i.e. `outputs/06_outlier_check/
{predictions,scores}/residual__<config>__repeat*_fold*.csv` -- NOT 05's original
unexcluded models. Both levers below reuse 06's fixed training-pool exclusion
(`outputs/06_outlier_check/flagged_compounds/criterion1_residual_flagged.csv`, the 75
compounds flagged by CV-residual p95) as their own starting training pool; CI-width
exclusion (criterion 2) is out of scope here.

SCOPE -- `chemprop_chemeleoninit` is explicitly OUT for this entire step; tabular
configs only (`chemeleon__rf`, `chemeleon__lightgbm`). Reason: 06's Sections 5-8 found
that `chemprop_chemeleoninit`'s residual-excluded result is entangled between a
genuine data-quality effect and a leading, unconfirmed hypothesis that masking CYP2D6
labels implicitly reweights CYP2D6's share of the shared multitask encoder's loss.
Adding a third mechanism (explicit sample weighting) on top of an already-unresolved
two-mechanism baseline would make attribution substantially harder, so
`chemprop_chemeleoninit` is deferred rather than tested here -- its 06 result
(excluded, unweighted, untuned) is its best-known state going into step 3 unless
revisited. `chemeleon__rf` and `chemeleon__lightgbm` are single-task models with no
shared encoder -- unaffected by this mechanism, in scope as originally planned.

Two independent levers, tested SEPARATELY against the residual-excluded baseline (not
combined into one variant yet):

1. WEIGHTING (both configs): inverse-local-density sample weighting on pIC50, applied
   via `sample_weight` in `.fit()`. REFIT PER FOLD -- density is estimated separately
   for each of the 25 folds, from that fold's own post-exclusion training pool
   (`screen_inner_train` + `screen_inner_val`, CYP2D6-labeled, criterion-1-excluded
   compounds already dropped) -- see `compute_fold_cyp2d6_weights`. This deviates from
   the task's original "frozen once for all folds" framing: a single global fit over
   all 1,493 labeled compounds would use each fold's own held-out `screen_test`
   compounds' pIC50 values to shape that very fold's training weights -- a leakage
   risk, confirmed with the user (2026-09-03) and corrected before any retraining ran.
   Model hyperparameters otherwise match the 06 baseline (sklearn defaults,
   `random_state=seed`) -- weighting is the only change.

2. TUNING (both configs): grid search (not random/Bayesian -- precedented by Goossens
   et al.'s Polaris ADME challenge RF/XGBoost tabular baselines), scored on plain RAE
   (`src.vendor.openadmet_eval.custom_scoring_functions.rae` -- undecorated, no CI
   soft-thresholding; a deliberate simplification for this internal selection step
   only, see `grid_search_rf`/`grid_search_lgbm`) evaluated on each fold's own
   `screen_inner_val` split -- "one iteration of the inner loop of nested CV" (Ash et
   al., 5x5 CV paper, Section 3.1.4), not a second nested CV layered across all 25
   outer folds. The final before/after comparison against 06's baseline still uses the
   project's real bootstrapped ST-RAE, exactly as 06 did -- only which hyperparameters
   win the search uses this simpler criterion. Sample weighting is NOT applied during
   tuning -- the two levers are tested independently, per task spec, not combined.

RF GRID REDUCED (2026-09-03, confirmed with the user, supersedes the task's original
64-combo RF grid): `n_estimators` dropped from {200, 500, 1000, 1500} to {200, 500,
1000} and `max_depth` dropped from {None, 10, 20, 30} to {10, 20} -- 24 combos, down
from 64. `min_samples_leaf` {1, 2, 3, 5} is unchanged. Reason: with only 3,335
inner-training compounds and `min_samples_leaf` already searched down to 1,
unrestricted-depth trees at 1,500 estimators sit deep in overfitting territory for a
dataset this size -- and, separately, `n_estimators=1500`/`max_depth=None` on
2048-dim CheMeleon embeddings was taking over an hour per fold's 64-combo search
(~24-30h projected for all 25 folds), which is what actually surfaced the concern.
Both reasons point the same direction: the weighting lever's own result (a dramatic
distribution-manipulation costing far more than it recovered, see Section 3 of the
notebook) is independent evidence that extreme parameter choices aren't where this
task's useful signal lives, making this a proportionate, evidence-consistent
narrowing rather than an arbitrary compute-driven cut. `n_jobs=-1` (see
`grid_search_rf`/`retrain_weighted`/`retrain_tuned`) was already set on every RF fit
before this change -- confirmed present, not the fix here; it does not change results
(sklearn's RF is bit-identical regardless of `n_jobs` for a fixed `random_state`),
only wall-clock time. `chemeleon__lightgbm`'s grid (also 64 combos, ~75s/fold) is
untouched -- it was never the runtime bottleneck.

CORRECTED ASYMMETRY DISCLOSURE (2026-09-03, supersedes the task's original framing):
weighting is NOT frozen once for all folds -- it is refit per fold, on that fold's own
post-exclusion training pool only (see Lever 1 above), for the same fold-isolation
reason tuning's inner-split search already respects. What both levers share instead:
neither ever touches a fold's own `screen_test` compounds when fitting anything
(weights or hyperparameters) for that fold -- this is the invariant actually being
protected, and is disclosed here and again in the notebook.

Run:
    python scripts/cyp2d6_weighting_tuning.py --weighting-check
        # computes + plots the proposed weighting scheme only, no training. Stop
        # here for user confirmation before running --retrain (task spec).
    python scripts/cyp2d6_weighting_tuning.py --retrain
        # full run: both levers, both configs, all 25 folds
    python scripts/cyp2d6_weighting_tuning.py --retrain --configs chemeleon__rf \\
        --levers weighting
        # narrow scope, for reruns
"""

import os

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import argparse
import itertools
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
for _p in (REPO_ROOT, SCRIPTS_DIR):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import run_5x5_cv_comparison as cv5x5  # noqa: E402  (needs sys.path set up first)

from src.chemprop_screen import load_screen_population, setup_logging  # noqa: E402
from src.cv_bootstrap import per_fold_bootstrap_seed  # noqa: E402
from src.vendor.openadmet_eval.evaluate_predictions import score_activity_predictions  # noqa: E402

CURATED_PATH = REPO_ROOT / "data" / "processed" / "train_inhibition_curated.csv"
RESIDUAL_FLAGGED_PATH = (
    REPO_ROOT / "outputs" / "06_outlier_check" / "flagged_compounds" / "criterion1_residual_flagged.csv"
)
BASELINE_SCORE_DIR = REPO_ROOT / "outputs" / "06_outlier_check" / "scores"  # 06's residual-excluded scores (read-only "before")

OUT = REPO_ROOT / "outputs" / "07_weighting_tuning"
WEIGHTING_CHECK_DIR = OUT / "weighting_check"
PRED_DIR = OUT / "predictions"
SCORE_DIR = OUT / "scores"
TUNING_SEARCH_DIR = OUT / "tuning_search"
LOG_DIR = REPO_ROOT / "logs"

CYP2D6_COL = "CYP2D6_pIC50_direct_inhibition"

N_REPEATS = 5
N_FOLDS = 5

# chemprop_chemeleoninit deliberately excluded from this step -- see module docstring
# ("SCOPE") for the full reasoning (06 Sections 5-8: masking-vs-multitask-loss-
# reweighting entanglement, not tested here to avoid a third confounded mechanism).
CONFIGS = ["chemeleon__rf", "chemeleon__lightgbm"]
LEVERS = ["weighting", "tuning"]

WEIGHT_CLIP_MAX = 5.0  # proposed cap -- see compute_cyp2d6_sample_weights docstring; confirmed with user before use

# Reduced from the task's original 64-combo grid ({200,500,1000,1500} x {None,10,20,30}
# x {1,2,3,5}) on 2026-09-03 -- see module docstring "RF GRID REDUCED" for the full,
# confirmed reasoning (runtime + overfitting-territory evidence, not a silent cut).
RF_GRID = {
    "n_estimators": [200, 500, 1000],
    "max_depth": [10, 20],
    "min_samples_leaf": [1, 2, 3, 5],
}
LGBM_GRID = {
    "num_leaves": [15, 31, 63, 127],
    "learning_rate": [0.01, 0.03, 0.05, 0.1],
    "min_child_samples": [5, 10, 20, 30],
}
XGB_LGBM_EARLY_STOPPING_ROUNDS = cv5x5.XGB_LGBM_EARLY_STOPPING_ROUNDS  # matches baseline precedent (10)


# ---------------------------------------------------------------------------
# Lever 1: weighting scheme
# ---------------------------------------------------------------------------


def load_residual_excluded_inchikeys() -> set:
    flagged = pd.read_csv(RESIDUAL_FLAGGED_PATH)
    return set(flagged["inchikey"])


def get_manifest_seed(manifest: pd.DataFrame, config: str, repeat: int, fold: int) -> int:
    """Same seed 06/05 used for this (config, repeat, fold) -- shared across configs
    within a (repeat, fold) by `generate_5x5_cv_manifest.py`'s own design (one seed per
    (repeat, fold), reused everywhere that (repeat, fold) is touched), so any `config`
    argument here returns the same value for a given (repeat, fold).
    """
    row = manifest[(manifest["config"] == config) & (manifest["repeat"] == repeat) & (manifest["fold"] == fold)]
    if len(row) != 1:
        raise ValueError(f"expected exactly one manifest row for {config} repeat={repeat} fold={fold}, got {len(row)}")
    return int(row["seed"].iloc[0])


def compute_fold_cyp2d6_weights(
    population: pd.DataFrame, excluded_inchikeys: set, clip_max: float = WEIGHT_CLIP_MAX
) -> pd.DataFrame:
    """Inverse-local-density sample weights for ONE fold, fit only on that fold's own
    post-exclusion training pool (`screen_inner_train` + `screen_inner_val`,
    CYP2D6-labeled, criterion-1-flagged compounds already dropped) -- never on that
    fold's `screen_test` compounds. Refitting per fold (rather than once globally)
    avoids leaking a fold's own held-out labels into its own training weights --
    confirmed with the user (2026-09-03), see module docstring "CORRECTED ASYMMETRY
    DISCLOSURE".

    CONFIRMED SCHEME: Gaussian KDE (`scipy.stats.gaussian_kde`, Scott's-rule
    bandwidth -- scipy's own default, `bw = n**(-1/(d+4))` with d=1) fit on this pool's
    pIC50 values. Raw weight = 1 / density(pIC50); normalized so the mean weight
    within this pool is exactly 1.0 (keeps sum(weights) ~= pool size, so the weighted
    loss stays on the same scale as the unweighted loss and RF/LightGBM's default
    hyperparameters -- unchanged for this lever -- remain comparable to the unweighted
    baseline). Clipped at `clip_max` (confirmed: 5.0x the pool's own mean) to stop the
    handful of most-extreme high/low-potency compounds in this pool from
    single-handedly dominating the fit.
    """
    pool = population.loc[
        population["screen_split"].isin(["screen_inner_train", "screen_inner_val"])
        & population[CYP2D6_COL].notna()
        & ~population["inchikey"].isin(excluded_inchikeys),
        ["inchikey", "Molecule_Name", CYP2D6_COL],
    ].reset_index(drop=True)
    values = pool[CYP2D6_COL].to_numpy()

    kde = gaussian_kde(values, bw_method="scott")
    density = kde(values)
    raw_weight = 1.0 / density
    normalized_weight = raw_weight / raw_weight.mean()
    n_clipped = int((normalized_weight > clip_max).sum())
    clipped_weight = np.clip(normalized_weight, a_min=None, a_max=clip_max)

    pool["density"] = density
    pool["weight_unclipped"] = normalized_weight
    pool["weight"] = clipped_weight
    return pool


def run_weighting_check() -> pd.DataFrame:
    """Refit the confirmed weighting scheme per fold (all 25 repeat/fold pairs) and
    write its distribution (CSV + PNG) to `outputs/07_weighting_tuning/
    weighting_check/` for a second user sanity-check, per-fold this time -- NO
    retraining happens in this mode.
    """
    import matplotlib.pyplot as plt

    curated = pd.read_csv(CURATED_PATH)
    print(f"loaded {CURATED_PATH.name}: {curated.shape}")
    excluded_inchikeys = load_residual_excluded_inchikeys()
    print(f"loaded {RESIDUAL_FLAGGED_PATH.name}: {len(excluded_inchikeys)} criterion-1-flagged compounds (train-only excluded)")
    manifest = pd.read_csv(cv5x5.MANIFEST_PATH)
    logger = setup_logging(LOG_DIR / "07_weighting_check.log", "cyp2d6_weighting_check")

    rows = []
    for repeat in range(N_REPEATS):
        for fold in range(N_FOLDS):
            seed = get_manifest_seed(manifest, "chemeleon__rf", repeat, fold)  # config-agnostic: seed is shared per (repeat, fold)
            repeat_col = f"repeat_{repeat}"
            population = load_screen_population(
                cv5x5.FOLDS_PATH, cv5x5.CURATED_PATH, repeat_col, fold, cv5x5.VAL_FRACTION, seed, logger
            )
            fold_weights = compute_fold_cyp2d6_weights(population, excluded_inchikeys)
            fold_weights["repeat"] = repeat
            fold_weights["fold"] = fold
            rows.append(fold_weights)

    all_weights = pd.concat(rows, ignore_index=True)
    per_fold_summary = all_weights.groupby(["repeat", "fold"])["weight"].agg(["size", "mean", "median", "min", "max"])
    print(f"\nper-fold pool size: min={per_fold_summary['size'].min()}, max={per_fold_summary['size'].max()}")
    print("per-fold weight mean: " f"mean-of-means={per_fold_summary['mean'].mean():.4f}, std-of-means={per_fold_summary['mean'].std():.4f}")
    print("\npooled weight distribution across all 25 folds (n = sum of 25 per-fold pool sizes, compounds repeated across folds):")
    print(all_weights["weight"].describe().to_string())

    WEIGHTING_CHECK_DIR.mkdir(parents=True, exist_ok=True)
    weights_path = WEIGHTING_CHECK_DIR / "cyp2d6_sample_weights_per_fold.csv"
    all_weights.to_csv(weights_path, index=False)
    summary_path = WEIGHTING_CHECK_DIR / "cyp2d6_sample_weights_per_fold_summary.csv"
    per_fold_summary.reset_index().to_csv(summary_path, index=False)
    print(f"\nwrote {weights_path}\nwrote {summary_path}")

    fig, axes = plt.subplots(1, 2, figsize=(15, 5))

    fold_order = [(r, f) for r in range(N_REPEATS) for f in range(N_FOLDS)]
    fold_labels = [f"r{r}f{f}" for r, f in fold_order]
    box_data = [all_weights.loc[(all_weights["repeat"] == r) & (all_weights["fold"] == f), "weight"] for r, f in fold_order]
    axes[0].boxplot(box_data, tick_labels=fold_labels, showfliers=False, patch_artist=True,
                     boxprops={"facecolor": "lightgray"}, medianprops={"color": "black"})
    for i, vals in enumerate(box_data, start=1):
        jitter = np.random.default_rng(0).normal(0, 0.06, size=len(vals))
        axes[0].scatter(np.full(len(vals), i) + jitter, vals, color="black", s=3, alpha=0.3)
    axes[0].axhline(1.0, color="red", linestyle="--", linewidth=1, label="mean weight = 1.0 (within-fold)")
    axes[0].set_xticklabels(fold_labels, rotation=90, fontsize=6)
    axes[0].set_ylabel("sample weight (post-clip)")
    axes[0].set_title("per-fold weight distribution (25 folds)")
    axes[0].legend(fontsize=8)

    axes[1].hist(all_weights["weight"], bins=40, color="lightgray", edgecolor="black", linewidth=0.5)
    axes[1].axvline(1.0, color="red", linestyle="--", linewidth=1, label="mean weight = 1.0 (within-fold)")
    axes[1].axvline(WEIGHT_CLIP_MAX, color="black", linestyle=":", linewidth=1, label=f"clip_max = {WEIGHT_CLIP_MAX}")
    axes[1].set_xlabel("sample weight (post-clip)")
    axes[1].set_ylabel("count (pooled across 25 folds)")
    axes[1].set_title("pooled weight distribution, all 25 folds")
    axes[1].legend(fontsize=8)

    fig.suptitle("CYP2D6 sample-weighting sanity check -- confirmed method, refit per fold (2nd review)")
    plt.tight_layout()
    hist_path = WEIGHTING_CHECK_DIR / "weight_distribution_check_per_fold.png"
    fig.savefig(hist_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {hist_path}")

    return all_weights


# ---------------------------------------------------------------------------
# Shared population/feature setup (both levers)
# ---------------------------------------------------------------------------


def prepare_fold_arrays(config: str, repeat: int, fold: int, seed: int, excluded_inchikeys: set, logger):
    """Load this fold's population, join it to the chemeleon feature matrix by
    position, and split into the index arrays both levers need. Mirrors
    `cyp2d6_outlier_check.retrain_tabular_cyp2d6`'s own setup exactly (same
    train-only exclusion: a flagged compound is dropped from the training pool but
    still scored normally if it lands in this fold's held-out test set) -- factored
    out here since both levers need it identically, only what happens with
    `has_label` inside `.fit()` differs between them.
    """
    repeat_col = f"repeat_{repeat}"
    population = load_screen_population(
        cv5x5.FOLDS_PATH, cv5x5.CURATED_PATH, repeat_col, fold, cv5x5.VAL_FRACTION, seed, logger
    )
    test_ids = population.loc[
        population["screen_split"] == "screen_test", ["Molecule_Name", "inchikey"]
    ].reset_index(drop=True)

    feature_name, algo = config.split("__")
    feat_index, X = cv5x5.load_feature_matrix(feature_name)
    pop_pos = feat_index.reset_index().merge(population, on=["Molecule_Name", "inchikey"], how="left")
    if pop_pos["screen_split"].isna().any():
        raise ValueError(
            f"{config} repeat={repeat} fold={fold}: feature index did not align 1:1 "
            "against the screen population -- stopping."
        )

    has_label_all = pop_pos[CYP2D6_COL].notna()
    is_excluded = pop_pos["inchikey"].isin(excluded_inchikeys)
    has_label = has_label_all & ~is_excluded
    logger.info(
        f"{config} repeat={repeat} fold={fold}: CYP2D6-labeled rows in this fold's population="
        f"{int(has_label_all.sum())}, excluded from training pool={int((has_label_all & is_excluded).sum())}, "
        f"remaining for training={int(has_label.sum())}"
    )

    return {
        "population": population,
        "test_ids": test_ids,
        "algo": algo,
        "X": X,
        "pop_pos": pop_pos,
        "has_label": has_label,
        "pool_mask": pop_pos["screen_split"].isin(["screen_inner_train", "screen_inner_val"]),
        "inner_train_mask": pop_pos["screen_split"] == "screen_inner_train",
        "inner_val_mask": pop_pos["screen_split"] == "screen_inner_val",
        "test_positions": pop_pos.loc[pop_pos["screen_split"] == "screen_test", "index"].to_numpy(),
    }


# ---------------------------------------------------------------------------
# Lever 1: weighting -- retrain
# ---------------------------------------------------------------------------


def retrain_weighted(
    config: str, repeat: int, fold: int, seed: int, excluded_inchikeys: set, logger
) -> pd.DataFrame:
    """CYP2D6-only retrain with inverse-local-density `sample_weight` applied to the
    training objective -- hyperparameters otherwise match the 06 baseline (sklearn
    defaults, `random_state=seed`; `n_jobs=-1` added for RF only, a performance-only
    change -- sklearn's RF is bit-identical regardless of `n_jobs` given a fixed
    `random_state`, matching this repo's own documented `n_jobs`-independent-result
    convention). Sample weight is applied to the training fit only, NOT to LightGBM's
    early-stopping `eval_set` metric -- early stopping is a model-selection device,
    not part of the loss being reweighted; only the training objective changes for
    this lever, matching the task's own framing ("count more in the loss").
    """
    prep = prepare_fold_arrays(config, repeat, fold, seed, excluded_inchikeys, logger)
    population, test_ids, algo, X, pop_pos, has_label = (
        prep["population"], prep["test_ids"], prep["algo"], prep["X"], prep["pop_pos"], prep["has_label"]
    )
    pool_mask, inner_train_mask, inner_val_mask, test_positions = (
        prep["pool_mask"], prep["inner_train_mask"], prep["inner_val_mask"], prep["test_positions"]
    )

    fold_weights = compute_fold_cyp2d6_weights(population, excluded_inchikeys)
    weight_lookup = fold_weights.set_index("inchikey")["weight"]

    if algo == "rf":
        idx = pop_pos.loc[pool_mask & has_label, "index"].to_numpy()
        y = pop_pos.loc[pool_mask & has_label, CYP2D6_COL].to_numpy()
        sw = pop_pos.loc[pool_mask & has_label, "inchikey"].map(weight_lookup).to_numpy()
        if np.isnan(sw).any():
            raise ValueError(f"{config} repeat={repeat} fold={fold}: missing sample weight for a pooled training compound -- stopping.")
        model = cv5x5.RandomForestRegressor(random_state=seed, n_jobs=-1)
        model.fit(X[idx], y, sample_weight=sw)
    elif algo == "lightgbm":
        tr_idx = pop_pos.loc[inner_train_mask & has_label, "index"].to_numpy()
        va_idx = pop_pos.loc[inner_val_mask & has_label, "index"].to_numpy()
        y_tr = pop_pos.loc[inner_train_mask & has_label, CYP2D6_COL].to_numpy()
        y_va = pop_pos.loc[inner_val_mask & has_label, CYP2D6_COL].to_numpy()
        sw_tr = pop_pos.loc[inner_train_mask & has_label, "inchikey"].map(weight_lookup).to_numpy()
        if np.isnan(sw_tr).any():
            raise ValueError(f"{config} repeat={repeat} fold={fold}: missing sample weight for an inner_train compound -- stopping.")
        model = cv5x5.LGBMRegressor(random_state=seed, verbosity=-1)
        model.fit(
            X[tr_idx], y_tr, sample_weight=sw_tr, eval_set=[(X[va_idx], y_va)],
            callbacks=[cv5x5.lgb.early_stopping(stopping_rounds=XGB_LGBM_EARLY_STOPPING_ROUNDS, verbose=False)],
        )
    else:
        raise ValueError(f"unknown algorithm: {algo}")

    pred_df = test_ids.copy()
    pred_df[CYP2D6_COL] = np.asarray(model.predict(X[test_positions])).reshape(-1)
    return pred_df


# ---------------------------------------------------------------------------
# Lever 2: tuning -- grid search + retrain
# ---------------------------------------------------------------------------


def grid_search_rf(X: np.ndarray, tr_idx: np.ndarray, y_tr: np.ndarray, va_idx: np.ndarray, y_va: np.ndarray, seed: int):
    """Grid search over `RF_GRID` (64 combos), fit on `screen_inner_train` only,
    scored by plain RAE (`src.vendor.openadmet_eval.custom_scoring_functions.rae` --
    undecorated, no CI soft-thresholding) on `screen_inner_val` -- "one iteration of
    the inner loop of nested CV" (Ash et al., Section 3.1.4), never touching this
    fold's `screen_test` compounds. Plain RAE (not the CI-aware ST-RAE) is a
    deliberate simplification for this internal selection step only -- the final
    before/after comparison against 06's baseline still uses the project's real
    bootstrapped ST-RAE exactly as 06 did; only which hyperparameters win the search
    uses this simpler, CI-agnostic criterion.

    Returns (best_params, best_inner_val_RAE, search_results_df) -- the RF grid does
    NOT return a fitted model (unlike LightGBM): the final RF model is refit on the
    pooled inner_train+inner_val set with the winning params, matching how the 06
    baseline's own RF call fits (RF has no early-stopping concept to justify holding
    inner_val out of the final fit).
    """
    from src.vendor.openadmet_eval.custom_scoring_functions import rae

    rows = []
    best_params, best_score = None, np.inf
    for n_estimators, max_depth, min_samples_leaf in itertools.product(
        RF_GRID["n_estimators"], RF_GRID["max_depth"], RF_GRID["min_samples_leaf"]
    ):
        model = cv5x5.RandomForestRegressor(
            n_estimators=n_estimators, max_depth=max_depth, min_samples_leaf=min_samples_leaf,
            random_state=seed, n_jobs=-1,
        )
        model.fit(X[tr_idx], y_tr)
        pred_va = model.predict(X[va_idx])
        score = float(rae(y_va, pred_va))
        rows.append({
            "n_estimators": n_estimators, "max_depth": max_depth, "min_samples_leaf": min_samples_leaf,
            "inner_val_RAE": score,
        })
        if score < best_score:
            best_score = score
            best_params = {"n_estimators": n_estimators, "max_depth": max_depth, "min_samples_leaf": min_samples_leaf}
    return best_params, best_score, pd.DataFrame(rows)


def grid_search_lgbm(X: np.ndarray, tr_idx: np.ndarray, y_tr: np.ndarray, va_idx: np.ndarray, y_va: np.ndarray, seed: int):
    """Grid search over `LGBM_GRID` (64 combos), each combo fit on `screen_inner_train`
    with `screen_inner_val` as its early-stopping `eval_set` (unweighted default loss,
    same as the 06 baseline) -- scored by plain RAE on `screen_inner_val` predictions,
    same simplification/rationale as `grid_search_rf`.

    Returns (best_params, best_inner_val_RAE, best_fitted_model, search_results_df).
    Unlike RF, the winning combo's ALREADY-FITTED model is reused directly as the
    final model -- matching the 06 baseline's own LightGBM convention (fit on
    inner_train, early-stopped against inner_val, used as-is; no separate pooled
    refit), so no second fit is needed here.
    """
    from src.vendor.openadmet_eval.custom_scoring_functions import rae

    rows = []
    best_params, best_score, best_model = None, np.inf, None
    for num_leaves, learning_rate, min_child_samples in itertools.product(
        LGBM_GRID["num_leaves"], LGBM_GRID["learning_rate"], LGBM_GRID["min_child_samples"]
    ):
        model = cv5x5.LGBMRegressor(
            num_leaves=num_leaves, learning_rate=learning_rate, min_child_samples=min_child_samples,
            random_state=seed, verbosity=-1,
        )
        model.fit(
            X[tr_idx], y_tr, eval_set=[(X[va_idx], y_va)],
            callbacks=[cv5x5.lgb.early_stopping(stopping_rounds=XGB_LGBM_EARLY_STOPPING_ROUNDS, verbose=False)],
        )
        pred_va = model.predict(X[va_idx])
        score = float(rae(y_va, pred_va))
        rows.append({
            "num_leaves": num_leaves, "learning_rate": learning_rate, "min_child_samples": min_child_samples,
            "inner_val_RAE": score,
        })
        if score < best_score:
            best_score = score
            best_params = {"num_leaves": num_leaves, "learning_rate": learning_rate, "min_child_samples": min_child_samples}
            best_model = model
    return best_params, best_score, best_model, pd.DataFrame(rows)


def retrain_tuned(
    config: str, repeat: int, fold: int, seed: int, excluded_inchikeys: set, logger
) -> pd.DataFrame:
    """CYP2D6-only retrain with grid-searched hyperparameters (no sample weighting --
    the two levers are tested independently, per task spec). Writes this fold's full
    grid-search result table to `outputs/07_weighting_tuning/tuning_search/` (used by
    the notebook's search-landscape plot), in addition to returning the final
    predictions.
    """
    prep = prepare_fold_arrays(config, repeat, fold, seed, excluded_inchikeys, logger)
    algo, X, pop_pos, has_label = prep["algo"], prep["X"], prep["pop_pos"], prep["has_label"]
    pool_mask, inner_train_mask, inner_val_mask, test_positions, test_ids = (
        prep["pool_mask"], prep["inner_train_mask"], prep["inner_val_mask"], prep["test_positions"], prep["test_ids"]
    )

    tr_idx = pop_pos.loc[inner_train_mask & has_label, "index"].to_numpy()
    va_idx = pop_pos.loc[inner_val_mask & has_label, "index"].to_numpy()
    y_tr = pop_pos.loc[inner_train_mask & has_label, CYP2D6_COL].to_numpy()
    y_va = pop_pos.loc[inner_val_mask & has_label, CYP2D6_COL].to_numpy()

    if algo == "rf":
        best_params, best_score, search_df = grid_search_rf(X, tr_idx, y_tr, va_idx, y_va, seed)
        pool_idx = pop_pos.loc[pool_mask & has_label, "index"].to_numpy()
        y_pool = pop_pos.loc[pool_mask & has_label, CYP2D6_COL].to_numpy()
        final_model = cv5x5.RandomForestRegressor(random_state=seed, n_jobs=-1, **best_params)
        final_model.fit(X[pool_idx], y_pool)
    elif algo == "lightgbm":
        best_params, best_score, final_model, search_df = grid_search_lgbm(X, tr_idx, y_tr, va_idx, y_va, seed)
    else:
        raise ValueError(f"unknown algorithm: {algo}")

    search_df["config"] = config
    search_df["repeat"] = repeat
    search_df["fold"] = fold
    TUNING_SEARCH_DIR.mkdir(parents=True, exist_ok=True)
    search_path = TUNING_SEARCH_DIR / f"{config}__repeat{repeat}_fold{fold}.csv"
    search_df.to_csv(search_path, index=False)
    logger.info(
        f"{config} repeat={repeat} fold={fold}: grid search ({len(search_df)} combos) "
        f"best_params={best_params}, best inner_val RAE={best_score:.4f}, wrote {search_path.name}"
    )

    pred_df = test_ids.copy()
    pred_df[CYP2D6_COL] = np.asarray(final_model.predict(X[test_positions])).reshape(-1)
    return pred_df


# ---------------------------------------------------------------------------
# Resumable driver (mirrors cyp2d6_outlier_check.py's own pattern)
# ---------------------------------------------------------------------------


def result_paths(lever: str, config: str, repeat: int, fold: int) -> tuple[Path, Path]:
    tag = f"{lever}__{config}__repeat{repeat}_fold{fold}"
    return PRED_DIR / f"{tag}.csv", SCORE_DIR / f"{tag}.csv"


def is_retrain_done(lever: str, config: str, repeat: int, fold: int) -> bool:
    pred_path, score_path = result_paths(lever, config, repeat, fold)
    if not (pred_path.exists() and score_path.exists()):
        return False
    try:
        return len(pd.read_csv(pred_path)) > 0 and len(pd.read_csv(score_path)) > 0
    except (pd.errors.EmptyDataError, OSError):
        return False


def score_and_save_cyp2d6(
    lever: str, config: str, repeat: int, fold: int, seed: int,
    pred_df: pd.DataFrame, curated: pd.DataFrame, logger
) -> None:
    pred_path, score_path = result_paths(lever, config, repeat, fold)
    PRED_DIR.mkdir(parents=True, exist_ok=True)
    SCORE_DIR.mkdir(parents=True, exist_ok=True)

    ground_truth = curated[curated["inchikey"].isin(pred_df["inchikey"])].copy()
    with per_fold_bootstrap_seed(seed):
        scored = score_activity_predictions(pred_df, ground_truth, [CYP2D6_COL])
    scored["config"] = config
    scored["lever"] = lever
    scored["repeat"] = repeat
    scored["fold"] = fold
    scored["bootstrap_seed"] = seed

    pred_df.to_csv(pred_path, index=False)
    scored.to_csv(score_path, index=False)
    logger.info(f"wrote {pred_path.name} ({len(pred_df)} rows), {score_path.name} ({len(scored)} rows)")


def run_retrain(
    configs: list[str], levers: list[str], curated: pd.DataFrame, excluded_inchikeys: set, logger
) -> None:
    manifest = pd.read_csv(cv5x5.MANIFEST_PATH)
    n_done, n_run, n_failed = 0, 0, 0
    for config in configs:
        for lever in levers:
            for repeat in range(N_REPEATS):
                for fold in range(N_FOLDS):
                    if is_retrain_done(lever, config, repeat, fold):
                        n_done += 1
                        continue
                    seed = get_manifest_seed(manifest, config, repeat, fold)
                    logger.info(f"--- retrain lever={lever} config={config} repeat={repeat} fold={fold} seed={seed} ---")
                    t0 = time.time()
                    try:
                        if lever == "weighting":
                            pred_df = retrain_weighted(config, repeat, fold, seed, excluded_inchikeys, logger)
                        else:
                            pred_df = retrain_tuned(config, repeat, fold, seed, excluded_inchikeys, logger)
                        score_and_save_cyp2d6(lever, config, repeat, fold, seed, pred_df, curated, logger)
                        n_run += 1
                        logger.info(f"{lever}/{config} repeat={repeat} fold={fold} done in {time.time() - t0:.1f}s")
                    except Exception:
                        n_failed += 1
                        logger.error(
                            f"{lever}/{config} repeat={repeat} fold={fold} FAILED after {time.time() - t0:.1f}s "
                            f"-- no output written, will retry on next invocation:\n{traceback.format_exc()}"
                        )
    logger.info(f"retrain finished: {n_done} already done, {n_run} completed this run, {n_failed} failed")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--weighting-check", action="store_true",
        help="compute + plot the confirmed weighting scheme's per-fold distribution only, no training",
    )
    parser.add_argument("--retrain", action="store_true", help="run both levers' full retrain")
    parser.add_argument("--configs", nargs="+", default=CONFIGS, choices=CONFIGS)
    parser.add_argument("--levers", nargs="+", default=LEVERS, choices=LEVERS)
    args = parser.parse_args()

    if args.weighting_check:
        run_weighting_check()
        return

    if not args.retrain:
        parser.error("specify --weighting-check or --retrain")

    curated = pd.read_csv(CURATED_PATH)
    print(f"loaded {CURATED_PATH.name}: {curated.shape}")
    excluded_inchikeys = load_residual_excluded_inchikeys()
    print(f"loaded {RESIDUAL_FLAGGED_PATH.name}: {len(excluded_inchikeys)} criterion-1-flagged compounds (train-only excluded)")

    logger = setup_logging(LOG_DIR / "07_weighting_tuning_retrain.log", "cyp2d6_weighting_tuning_retrain")
    logger.info("=" * 70)
    logger.info(f"starting retrain: configs={args.configs} levers={args.levers}")
    run_retrain(args.configs, args.levers, curated, excluded_inchikeys, logger)


if __name__ == "__main__":
    main()
