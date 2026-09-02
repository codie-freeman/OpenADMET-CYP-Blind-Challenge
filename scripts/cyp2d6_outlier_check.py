"""CYP2D6 outlier check (Roadmap Step 1, notebook 06).

Two independent, model-choice-agnostic flagging criteria, computed from existing
05_cv_comparison artifacts -- no new training for the flagging step itself:

1. CV-residual flagging: pools out-of-fold residuals (y_pred - y_true) for every
   CYP2D6-labeled compound across chemprop_chemeleoninit's 25 (5 repeat x 5 fold) OOF
   prediction files only -- deliberately not all 12 configs, so outlier status is
   model-independent once frozen (CLAUDE.md / task spec). Threshold CONFIRMED by user:
   flag mean-absolute-residual > the 95th percentile of the pooled distribution.

2. CI-width flagging: flags the top 5% of CYP2D6 compounds by (conf_high - conf_low),
   read directly from train_inhibition_curated.csv. No model predictions involved.

Both flagged-compound lists are written to disk (outputs/06_outlier_check/
flagged_compounds/) before any retraining, so the excluded set is fixed and auditable.

RETRAINING (--retrain), gated behind the confirmed thresholds above. For each of
chemeleon__rf, chemeleon__lightgbm, ecfp4_narrow__lightgbm, chemprop_chemeleoninit, on
the same frozen 25-fold partition (data/folds/cv_folds.csv, untouched):

- Exclusion is TRAIN-ONLY: a flagged compound is dropped from whichever fold's
  training pool it would have been part of, but still scored normally whenever it
  falls in a fold's held-out test set. This keeps the "before" and "after" evaluation
  population identical (same 1,493 CYP2D6-labeled compounds, same 25 folds) so any
  metric change reflects the retrained model, not a smaller/easier test set --
  confirmed with the user (the alternative, dropping flagged compounds from eval too,
  would conflate the two).
- For chemprop_chemeleoninit's multitask model, exclusion MASKS only the CYP2D6 label
  for flagged compounds (set to NaN in the training CSV -- chemprop already handles
  missing multitask targets, since every other endpoint has missing labels for most
  compounds too); the compound keeps contributing to CYP1A2/2C9/3A4 training.
  Confirmed with the user, to isolate the intervention to CYP2D6 -- matching how the
  tabular configs are naturally isolated already (one model per endpoint).
- Each (config, repeat, fold)'s retrain reuses the exact bootstrap/model seed
  `scripts/generate_5x5_cv_manifest.py` originally assigned that row (from
  `outputs/05_cv_comparison/manifest.csv`), so training-data exclusion is the only
  thing that differs between the original ("before") and retrained ("after") run --
  not incidental RNG noise.
- Reuses `src.chemprop_screen`/`scripts.run_5x5_cv_comparison`'s already-validated
  population loading, feature loading, and Chemprop CLI orchestration rather than
  reimplementing them (CLAUDE.md: one shared path per concern). Only scores the
  CYP2D6 endpoint (`score_activity_predictions(..., endpoints=["CYP2D6_..."])`) --
  the other 3 isoforms/the MA pseudo-endpoint are out of scope for this check.

Resumable the same way as `run_5x5_cv_comparison.py`: a (criterion, config, repeat,
fold) row is "done" purely by its own prediction+score CSVs existing and non-empty.

Run:
    python scripts/cyp2d6_outlier_check.py                  # flagging report only
    python scripts/cyp2d6_outlier_check.py --retrain         # + exclude/retrain/score
    python scripts/cyp2d6_outlier_check.py --retrain \\
        --configs chemeleon__rf --criteria residual          # narrow scope, for reruns
"""

import os

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import argparse
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
for _p in (REPO_ROOT, SCRIPTS_DIR):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import run_5x5_cv_comparison as cv5x5  # noqa: E402  (needs sys.path set up first)

from src.chemprop_screen import (  # noqa: E402
    build_predict_csv,
    build_training_csv,
    load_screen_population,
    run_chemprop_predict,
    run_chemprop_train,
    setup_logging,
    verify_predictions,
)
from src.cv_bootstrap import per_fold_bootstrap_seed  # noqa: E402
from src.vendor.openadmet_eval.config import REGRESSION_ENDPOINTS  # noqa: E402
from src.vendor.openadmet_eval.evaluate_predictions import score_activity_predictions  # noqa: E402

CURATED_PATH = REPO_ROOT / "data" / "processed" / "train_inhibition_curated.csv"
SRC_PRED_DIR = REPO_ROOT / "outputs" / "05_cv_comparison" / "predictions"  # 05's OOF preds (read-only input)

OUT = REPO_ROOT / "outputs" / "06_outlier_check"
FLAGGED_DIR = OUT / "flagged_compounds"
PRED_DIR = OUT / "predictions"
SCORE_DIR = OUT / "scores"
CHEMPROP_RUNS_DIR = OUT / "chemprop_runs"
LOG_DIR = REPO_ROOT / "logs"

CYP2D6_COL = "CYP2D6_pIC50_direct_inhibition"
CI_WIDTH_TOP_PCT = 0.05  # fixed by task spec
RESIDUAL_PERCENTILE_CONFIRMED = 95  # confirmed by user after reviewing the candidate-threshold report
N_REPEATS = 5
N_FOLDS = 5
RESIDUAL_CONFIG = "chemprop_chemeleoninit"
RETRAIN_CONFIGS = ["chemeleon__rf", "chemeleon__lightgbm", "ecfp4_narrow__lightgbm", "chemprop_chemeleoninit"]  # tabular (fast) first, chemprop (slow) last
CRITERIA = ["residual", "ci_width"]


# ---------------------------------------------------------------------------
# 1. Flagging
# ---------------------------------------------------------------------------


def load_cyp2d6_oof_residuals(pred_dir: Path, curated: pd.DataFrame) -> pd.DataFrame:
    """Pool (y_pred - y_true) for CYP2D6 across all 25 chemprop_chemeleoninit OOF
    prediction files. Returns one row per (inchikey, repeat, fold, residual) --
    aggregation to one row per compound happens in the caller so the raw per-repeat
    residual count actually observed (<=5) stays inspectable rather than assumed.
    """
    truth = curated.loc[
        curated[CYP2D6_COL].notna(), ["inchikey", "Molecule_Name", CYP2D6_COL]
    ].rename(columns={CYP2D6_COL: "y_true"})

    rows = []
    n_files = 0
    for repeat in range(N_REPEATS):
        for fold in range(N_FOLDS):
            path = pred_dir / f"{RESIDUAL_CONFIG}__repeat{repeat}_fold{fold}.csv"
            if not path.exists():
                raise FileNotFoundError(f"expected OOF prediction file missing: {path}")
            n_files += 1
            pred = pd.read_csv(path)[["inchikey", CYP2D6_COL]].rename(
                columns={CYP2D6_COL: "y_pred"}
            )
            merged = pred.merge(truth, on="inchikey", how="inner")
            merged["repeat"] = repeat
            merged["fold"] = fold
            merged["residual"] = merged["y_pred"] - merged["y_true"]
            rows.append(merged)

    print(f"loaded {n_files} {RESIDUAL_CONFIG} OOF prediction files (expected {N_REPEATS * N_FOLDS})")
    all_residuals = pd.concat(rows, ignore_index=True)
    print(f"row count before per-compound pooling: {len(all_residuals)} "
          f"(CYP2D6-labeled compounds: {truth['inchikey'].nunique()})")
    return all_residuals


def pool_mean_abs_residual(oof_residuals: pd.DataFrame) -> pd.DataFrame:
    """One row per compound: mean absolute residual + how many OOF observations
    (repeats) it was pooled from.
    """
    pooled = (
        oof_residuals.groupby(["inchikey", "Molecule_Name"])
        .agg(mean_abs_residual=("residual", lambda s: s.abs().mean()), n_oof_obs=("residual", "size"))
        .reset_index()
    )
    print(f"row count after per-compound pooling: {len(pooled)}")
    print(f"n_oof_obs distribution: {pooled['n_oof_obs'].value_counts().sort_index().to_dict()}")
    return pooled


def flag_ci_width_top_pct(curated: pd.DataFrame, top_pct: float = CI_WIDTH_TOP_PCT) -> pd.DataFrame:
    """Flag the top `top_pct` of CYP2D6-labeled compounds by conf_high - conf_low."""
    labeled = curated.loc[curated[CYP2D6_COL].notna()].copy()
    labeled["ci_width"] = (
        labeled[f"{CYP2D6_COL}_conf_high"] - labeled[f"{CYP2D6_COL}_conf_low"]
    )
    if labeled["ci_width"].isna().any():
        raise ValueError("CYP2D6-labeled compound(s) missing conf_low/conf_high -- stopping.")

    threshold = labeled["ci_width"].quantile(1 - top_pct)
    flagged = labeled.loc[labeled["ci_width"] >= threshold].copy()
    print(f"CI-width flagging: n={len(labeled)} CYP2D6-labeled compounds, "
          f"{top_pct:.0%} quantile threshold={threshold:.4f}, flagged={len(flagged)} "
          f"({len(flagged) / len(labeled):.2%})")
    return flagged[["inchikey", "Molecule_Name", "ci_width"]].sort_values(
        "ci_width", ascending=False
    ).reset_index(drop=True)


def report_candidate_residual_thresholds(pooled: pd.DataFrame, ci_flagged: pd.DataFrame) -> None:
    """Print the mean-abs-residual distribution and, for a handful of candidate
    thresholds, how many compounds each would flag and how that set overlaps with the
    CI-width flagged set. Informational context alongside the now-confirmed p95 cut.
    """
    vals = pooled["mean_abs_residual"]
    percentiles = [50, 75, 90, 95, 97.5, 99]
    print("\nmean_abs_residual percentiles:")
    for p in percentiles:
        print(f"  p{p}: {np.percentile(vals, p):.4f}")
    print(f"  mean={vals.mean():.4f}  std={vals.std():.4f}  max={vals.max():.4f}")

    ci_flagged_keys = set(ci_flagged["inchikey"])
    print("\ncandidate thresholds:")
    candidates = sorted(set(np.percentile(vals, p) for p in [90, 95, 97.5, 99]))
    for t in candidates:
        flagged = pooled.loc[pooled["mean_abs_residual"] > t]
        flagged_keys = set(flagged["inchikey"])
        overlap = flagged_keys & ci_flagged_keys
        smaller = min(len(flagged_keys), len(ci_flagged_keys)) or 1
        print(
            f"  threshold={t:.4f} (> ): flagged={len(flagged)} "
            f"({len(flagged) / len(pooled):.2%}); overlap with CI-width flagged: "
            f"{len(overlap)} ({len(overlap) / smaller:.2%} of smaller set)"
        )


def flag_residual_confirmed(pooled: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    """CONFIRMED criterion 1: flag mean_abs_residual strictly above the p95 of the
    pooled per-compound distribution.
    """
    threshold = float(np.percentile(pooled["mean_abs_residual"], RESIDUAL_PERCENTILE_CONFIRMED))
    flagged = pooled.loc[pooled["mean_abs_residual"] > threshold].copy()
    return flagged.sort_values("mean_abs_residual", ascending=False).reset_index(drop=True), threshold


def run_flagging(curated: pd.DataFrame) -> dict[str, pd.DataFrame]:
    oof_residuals = load_cyp2d6_oof_residuals(SRC_PRED_DIR, curated)
    pooled = pool_mean_abs_residual(oof_residuals)
    ci_flagged = flag_ci_width_top_pct(curated)
    report_candidate_residual_thresholds(pooled, ci_flagged)

    residual_flagged, residual_threshold = flag_residual_confirmed(pooled)
    overlap = set(residual_flagged["inchikey"]) & set(ci_flagged["inchikey"])
    smaller = min(len(residual_flagged), len(ci_flagged)) or 1
    print(
        f"\nCONFIRMED criterion 1 (residual, p{RESIDUAL_PERCENTILE_CONFIRMED}): "
        f"threshold={residual_threshold:.4f}, flagged={len(residual_flagged)}"
    )
    print(f"CONFIRMED overlap (criterion 1 vs criterion 2): {len(overlap)} ({len(overlap) / smaller:.2%} of smaller set)")

    FLAGGED_DIR.mkdir(parents=True, exist_ok=True)
    residual_path = FLAGGED_DIR / "criterion1_residual_flagged.csv"
    ci_path = FLAGGED_DIR / "criterion2_ci_width_flagged.csv"
    residual_flagged.to_csv(residual_path, index=False)
    ci_flagged.to_csv(ci_path, index=False)
    print(f"wrote {residual_path}\nwrote {ci_path}")

    return {"residual": residual_flagged, "ci_width": ci_flagged}


# ---------------------------------------------------------------------------
# 2. Retraining (train-only exclusion; see module docstring)
# ---------------------------------------------------------------------------


def retrain_tabular_cyp2d6(
    config: str, repeat: int, fold: int, seed: int, excluded_inchikeys: set, logger
) -> pd.DataFrame:
    """CYP2D6-only tabular retrain, mirroring `run_5x5_cv_comparison.run_tabular_or_naive`
    but restricted to the CYP2D6 endpoint and with flagged compounds dropped from the
    training pool only -- the screen_test population (and therefore what gets scored)
    is untouched.
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

    pool_mask = pop_pos["screen_split"].isin(["screen_inner_train", "screen_inner_val"])
    inner_train_mask = pop_pos["screen_split"] == "screen_inner_train"
    inner_val_mask = pop_pos["screen_split"] == "screen_inner_val"
    test_positions = pop_pos.loc[pop_pos["screen_split"] == "screen_test", "index"].to_numpy()

    if algo == "rf":
        idx = pop_pos.loc[pool_mask & has_label, "index"].to_numpy()
        y = pop_pos.loc[pool_mask & has_label, CYP2D6_COL].to_numpy()
        model = cv5x5.RandomForestRegressor(random_state=seed)
        model.fit(X[idx], y)
    elif algo in ("xgboost", "lightgbm"):
        tr_idx = pop_pos.loc[inner_train_mask & has_label, "index"].to_numpy()
        va_idx = pop_pos.loc[inner_val_mask & has_label, "index"].to_numpy()
        y_tr = pop_pos.loc[inner_train_mask & has_label, CYP2D6_COL].to_numpy()
        y_va = pop_pos.loc[inner_val_mask & has_label, CYP2D6_COL].to_numpy()
        if algo == "xgboost":
            model = cv5x5.XGBRegressor(
                early_stopping_rounds=cv5x5.XGB_LGBM_EARLY_STOPPING_ROUNDS, random_state=seed
            )
            model.fit(X[tr_idx], y_tr, eval_set=[(X[va_idx], y_va)], verbose=False)
        else:
            model = cv5x5.LGBMRegressor(random_state=seed, verbosity=-1)
            model.fit(
                X[tr_idx], y_tr, eval_set=[(X[va_idx], y_va)],
                callbacks=[cv5x5.lgb.early_stopping(stopping_rounds=cv5x5.XGB_LGBM_EARLY_STOPPING_ROUNDS, verbose=False)],
            )
    else:
        raise ValueError(f"unknown algorithm: {algo}")

    pred_df = test_ids.copy()
    pred_df[CYP2D6_COL] = np.asarray(model.predict(X[test_positions])).reshape(-1)
    return pred_df


def retrain_chemprop_cyp2d6(
    config: str, repeat: int, fold: int, seed: int, excluded_inchikeys: set, run_dir: Path, logger
) -> pd.DataFrame:
    """CYP2D6-excluded multitask Chemprop retrain, mirroring
    `run_5x5_cv_comparison.run_chemprop` but with flagged compounds' CYP2D6 label
    masked to NaN before building the training CSV -- chemprop already drops missing
    multitask targets from that task's loss, so this excludes them from CYP2D6
    training only, leaving their CYP1A2/2C9/3A4 contribution untouched.
    """
    repeat_col = f"repeat_{repeat}"
    population = load_screen_population(
        cv5x5.FOLDS_PATH, cv5x5.CURATED_PATH, repeat_col, fold, cv5x5.VAL_FRACTION, seed, logger
    )
    n_before = int(population[CYP2D6_COL].notna().sum())
    mask = population["inchikey"].isin(excluded_inchikeys)
    population.loc[mask, CYP2D6_COL] = np.nan
    n_after = int(population[CYP2D6_COL].notna().sum())
    logger.info(
        f"{config} repeat={repeat} fold={fold}: CYP2D6-labeled rows before exclusion-mask="
        f"{n_before}, after={n_after} (masked {n_before - n_after}); other 3 endpoints' labels untouched"
    )

    run_dir.mkdir(parents=True, exist_ok=True)
    predict_csv = build_predict_csv(population, run_dir, logger)
    train_csv = run_dir / "train_input.csv"
    build_training_csv(population, REGRESSION_ENDPOINTS, train_csv, logger, require_all_targets=False)
    run_chemprop_train(
        train_csv, REGRESSION_ENDPOINTS, run_dir, logger,
        cv5x5.CHEMPROP_ARCHITECTURE_ARGS[config], cv5x5.CHEMPROP_EPOCHS, cv5x5.CHEMPROP_PATIENCE, seed,
    )
    raw_pred_csv = run_dir / "raw_predictions.csv"
    run_chemprop_predict(run_dir / "model_0", predict_csv, raw_pred_csv, logger)

    expected_names = set(population.loc[population["screen_split"] == "screen_test", "Molecule_Name"])
    verify_predictions(raw_pred_csv, expected_names, REGRESSION_ENDPOINTS, logger)
    return pd.read_csv(raw_pred_csv)


def result_paths(criterion: str, config: str, repeat: int, fold: int) -> tuple[Path, Path]:
    tag = f"{criterion}__{config}__repeat{repeat}_fold{fold}"
    return PRED_DIR / f"{tag}.csv", SCORE_DIR / f"{tag}.csv"


def is_retrain_done(criterion: str, config: str, repeat: int, fold: int) -> bool:
    pred_path, score_path = result_paths(criterion, config, repeat, fold)
    if not (pred_path.exists() and score_path.exists()):
        return False
    try:
        return len(pd.read_csv(pred_path)) > 0 and len(pd.read_csv(score_path)) > 0
    except (pd.errors.EmptyDataError, OSError):
        return False


def score_and_save_cyp2d6(
    criterion: str, config: str, repeat: int, fold: int, seed: int,
    pred_df: pd.DataFrame, curated: pd.DataFrame, logger
) -> None:
    pred_path, score_path = result_paths(criterion, config, repeat, fold)
    PRED_DIR.mkdir(parents=True, exist_ok=True)
    SCORE_DIR.mkdir(parents=True, exist_ok=True)

    ground_truth = curated[curated["inchikey"].isin(pred_df["inchikey"])].copy()
    with per_fold_bootstrap_seed(seed):
        scored = score_activity_predictions(pred_df, ground_truth, [CYP2D6_COL])
    scored["config"] = config
    scored["criterion"] = criterion
    scored["repeat"] = repeat
    scored["fold"] = fold
    scored["bootstrap_seed"] = seed

    pred_df.to_csv(pred_path, index=False)
    scored.to_csv(score_path, index=False)
    logger.info(f"wrote {pred_path.name} ({len(pred_df)} rows), {score_path.name} ({len(scored)} rows)")


def get_manifest_seed(manifest: pd.DataFrame, config: str, repeat: int, fold: int) -> int:
    row = manifest[(manifest["config"] == config) & (manifest["repeat"] == repeat) & (manifest["fold"] == fold)]
    if len(row) != 1:
        raise ValueError(f"expected exactly one manifest row for {config} repeat={repeat} fold={fold}, got {len(row)}")
    return int(row["seed"].iloc[0])


def run_retrain(
    configs: list[str], criteria: list[str], curated: pd.DataFrame,
    flagged_by_criterion: dict[str, set], logger
) -> None:
    manifest = pd.read_csv(cv5x5.MANIFEST_PATH)
    n_done, n_run, n_failed = 0, 0, 0
    for config in configs:
        for criterion in criteria:
            excluded_inchikeys = flagged_by_criterion[criterion]
            for repeat in range(N_REPEATS):
                for fold in range(N_FOLDS):
                    if is_retrain_done(criterion, config, repeat, fold):
                        n_done += 1
                        continue
                    seed = get_manifest_seed(manifest, config, repeat, fold)
                    logger.info(f"--- retrain criterion={criterion} config={config} repeat={repeat} fold={fold} seed={seed} ---")
                    t0 = time.time()
                    try:
                        if config == "chemprop_chemeleoninit":
                            run_dir = CHEMPROP_RUNS_DIR / f"{criterion}__{config}__repeat{repeat}_fold{fold}"
                            pred_df = retrain_chemprop_cyp2d6(config, repeat, fold, seed, excluded_inchikeys, run_dir, logger)
                        else:
                            pred_df = retrain_tabular_cyp2d6(config, repeat, fold, seed, excluded_inchikeys, logger)
                        score_and_save_cyp2d6(criterion, config, repeat, fold, seed, pred_df, curated, logger)
                        n_run += 1
                        logger.info(f"{criterion}/{config} repeat={repeat} fold={fold} done in {time.time() - t0:.1f}s")
                    except Exception:
                        n_failed += 1
                        logger.error(
                            f"{criterion}/{config} repeat={repeat} fold={fold} FAILED after {time.time() - t0:.1f}s "
                            f"-- no output written, will retry on next invocation:\n{traceback.format_exc()}"
                        )
    logger.info(f"retrain finished: {n_done} already done, {n_run} completed this run, {n_failed} failed")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--retrain", action="store_true", help="run the exclude+retrain step after flagging")
    parser.add_argument("--configs", nargs="+", default=RETRAIN_CONFIGS, choices=RETRAIN_CONFIGS)
    parser.add_argument("--criteria", nargs="+", default=CRITERIA, choices=CRITERIA)
    args = parser.parse_args()

    curated = pd.read_csv(CURATED_PATH)
    print(f"loaded {CURATED_PATH.name}: {curated.shape}")

    flagged = run_flagging(curated)

    if not args.retrain:
        return

    logger = setup_logging(LOG_DIR / "06_outlier_check_retrain.log", "cyp2d6_outlier_check_retrain")
    logger.info("=" * 70)
    logger.info(f"starting retrain: configs={args.configs} criteria={args.criteria}")
    flagged_by_criterion = {name: set(df["inchikey"]) for name, df in flagged.items()}
    run_retrain(args.configs, args.criteria, curated, flagged_by_criterion, logger)


if __name__ == "__main__":
    main()
