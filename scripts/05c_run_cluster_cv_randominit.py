"""Resumable driver for notebook 05c's cluster-aware CV refit of `chemprop_randominit`
-- the one expensive config in 05c's 6-config leakage-sensitivity check (the tabular
configs are fast enough to fit directly in the notebook). Directly mirrors
`scripts/05b_run_cluster_cv.py`, which did the same thing for `chemprop_chemeleoninit`
-- see that script's own docstring for the full reasoning behind each piece; only what
differs is called out below.

WHAT DIFFERS FROM `05b_run_cluster_cv.py`:
- `CONFIG = "chemprop_randominit"` with its own architecture args (`--message-hidden-dim
  300 --depth 3 --aggregation mean --batch-norm`, copied verbatim from
  `run_5x5_cv_comparison.py`'s `CHEMPROP_ARCHITECTURE_ARGS["chemprop_randominit"]`)
  instead of CheMeleon-init's `--from-foundation CHEMELEON ...`.
- `FOLDS_PATH` points at **05b's own** `outputs/05b_cluster_cv_comparison/
  cluster_cv_folds.csv` -- reused read-only, exactly as-is. Notebook 05c does not
  rebuild clustering or fold assignment; this script never writes into
  `outputs/05b_cluster_cv_comparison/`, only reads one file from it.
- Output goes under `outputs/05c_cluster_cv_leakage_sensitivity/` (`predictions/`,
  `scores/`, `chemprop_runs/`) -- a directory notebook 05c's own tabular-config fits
  also write into, using the same `{config}__repeat{r}_fold{f}.csv` naming, so the
  notebook's aggregation step can glob all 6 configs uniformly regardless of which one
  trained in-notebook vs. via this script.
- **Training stays full 4-target multitask** (`REGRESSION_ENDPOINTS`, `require_all_
  targets=False`) -- identical to notebook 05's own `chemprop_randominit` run and to
  05b's `chemprop_chemeleoninit` run. Scoring is restricted to `SCORE_ENDPOINTS`
  (CYP1A2 + CYP3A4 only, notebook 05c's own stated scope) in `score_and_save` below --
  this is a scoring-time restriction, not a training-time one. Dropping CYP2C9/CYP2D6
  from *training* would risk exactly the multitask-loss-reweighting effect notebook 06
  flagged (excluding one isoform's labels shifted `chemprop_chemeleoninit`'s
  predictions on other, untouched isoforms) -- training scope must match notebook 05's
  real run exactly, so scoring is the only place this notebook's narrower isoform
  scope is applied.
- No macro ("MA") endpoint is added to the scored output -- a macro average over just
  2 of the 4 real isoforms would be a different, non-comparable quantity from the "MA"
  used everywhere else in this project, and 05c never asks for one.
- Seed cross-check at startup is against `05`'s manifest.csv filtered to
  `config == "chemprop_randominit"` (not `chemprop_chemeleoninit`).

USAGE:
    cd /path/to/OpenADMET-CYP-Blind-Challenge
    caffeinate -i nohup /Users/codiefreeman/miniconda3-arm64/envs/cyp-admet-v2/bin/python \\
        scripts/05c_run_cluster_cv_randominit.py > logs/05c_cluster_cv_randominit_stdout.log 2>&1 &
    tail -f logs/05c_cluster_cv_randominit.log

Expected cost: notebook 05's own `chemprop_randominit` run averaged ~40s/fold (vs.
CheMeleon-init's ~420s/fold) -- 25 folds should take on the order of 15-20 minutes, not
the ~3 hours 05b's CheMeleon-init refit took.
"""

import os

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

_env_bin = str(Path(sys.executable).parent)
_path_parts = os.environ.get("PATH", "").split(os.pathsep)
if _env_bin not in _path_parts:
    os.environ["PATH"] = os.pathsep.join([_env_bin, *_path_parts])

from src.chemprop_screen import (
    build_predict_csv,
    build_training_csv,
    load_screen_population,
    run_chemprop_predict,
    run_chemprop_train,
    setup_logging,
    verify_predictions,
)
from src.cv_bootstrap import per_fold_bootstrap_seed
from src.vendor.openadmet_eval.config import REGRESSION_ENDPOINTS
from src.vendor.openadmet_eval.evaluate_predictions import score_activity_predictions

CONFIG = "chemprop_randominit"
CURATED_PATH = REPO_ROOT / "data" / "processed" / "train_inhibition_curated.csv"
FOLDS_PATH = REPO_ROOT / "outputs" / "05b_cluster_cv_comparison" / "cluster_cv_folds.csv"  # read-only, reused as-is
ORIGINAL_MANIFEST_PATH = REPO_ROOT / "outputs" / "05_cv_comparison" / "manifest.csv"  # read-only cross-check
OUT = REPO_ROOT / "outputs" / "05c_cluster_cv_leakage_sensitivity"
PRED_DIR = OUT / "predictions"
SCORE_DIR = OUT / "scores"
CHEMPROP_RUNS_DIR = OUT / "chemprop_runs"
LOG_DIR = REPO_ROOT / "logs"

SCORE_ENDPOINTS = ["CYP1A2_pIC50_direct_inhibition", "CYP3A4_pIC50_direct_inhibition"]  # 05c's own scope

VAL_FRACTION = 0.15  # matches run_5x5_cv_comparison.py / 05b verbatim
CV_SEED_BASE = 42  # matches generate_5x5_cv_manifest.py's CV_SEED_BASE verbatim
N_REPEATS = 5
N_FOLDS = 5

CHEMPROP_ARCHITECTURE_ARGS = [
    "--message-hidden-dim", "300", "--depth", "3", "--aggregation", "mean", "--batch-norm",
]
CHEMPROP_EPOCHS = 50
CHEMPROP_PATIENCE = 5


def fold_seeds() -> dict:
    """Same formula, same enumeration order as `scripts/generate_5x5_cv_manifest.py`."""
    seed_sequences = np.random.SeedSequence(CV_SEED_BASE).spawn(N_REPEATS * N_FOLDS)
    seeds = {}
    i = 0
    for repeat in range(N_REPEATS):
        for fold in range(N_FOLDS):
            seeds[(repeat, fold)] = int(seed_sequences[i].generate_state(1)[0])
            i += 1
    return seeds


def cross_check_seeds_against_notebook_05(seeds: dict, logger) -> None:
    if not ORIGINAL_MANIFEST_PATH.exists():
        logger.warning(f"{ORIGINAL_MANIFEST_PATH} not found -- skipping seed cross-check against notebook 05")
        return
    manifest = pd.read_csv(ORIGINAL_MANIFEST_PATH)
    ref = manifest[manifest["config"] == CONFIG].set_index(["repeat", "fold"])["seed"].to_dict()
    mismatches = [k for k, v in seeds.items() if ref.get(k) != v]
    if mismatches:
        raise ValueError(
            f"seed mismatch vs. notebook 05's manifest.csv for {CONFIG} at (repeat, fold) = {mismatches} -- "
            "training randomness would silently differ from notebook 05, stopping."
        )
    logger.info(f"cross-checked all {len(seeds)} seeds against {ORIGINAL_MANIFEST_PATH} ({CONFIG} rows) -- identical.")


def result_paths(repeat: int, fold: int) -> tuple[Path, Path]:
    tag = f"{CONFIG}__repeat{repeat}_fold{fold}"
    return PRED_DIR / f"{tag}.csv", SCORE_DIR / f"{tag}.csv"


def is_done(repeat: int, fold: int) -> bool:
    pred_path, score_path = result_paths(repeat, fold)
    if not (pred_path.exists() and score_path.exists()):
        return False
    try:
        return len(pd.read_csv(pred_path)) > 0 and len(pd.read_csv(score_path)) > 0
    except (pd.errors.EmptyDataError, OSError):
        return False


def run_chemprop(repeat: int, fold: int, seed: int, logger) -> pd.DataFrame:
    repeat_col = f"repeat_{repeat}"
    population = load_screen_population(
        FOLDS_PATH, CURATED_PATH, repeat_col, fold, VAL_FRACTION, seed, logger
    )
    run_dir = CHEMPROP_RUNS_DIR / f"{CONFIG}__repeat{repeat}_fold{fold}"
    run_dir.mkdir(parents=True, exist_ok=True)

    predict_csv = build_predict_csv(population, run_dir, logger)
    train_csv = run_dir / "train_input.csv"
    # Full 4-target multitask training, matching notebook 05's real chemprop_randominit run --
    # see this file's own docstring for why CYP2C9/CYP2D6 stay in training despite not being scored.
    build_training_csv(population, REGRESSION_ENDPOINTS, train_csv, logger, require_all_targets=False)
    run_chemprop_train(
        train_csv, REGRESSION_ENDPOINTS, run_dir, logger,
        CHEMPROP_ARCHITECTURE_ARGS, CHEMPROP_EPOCHS, CHEMPROP_PATIENCE, seed,
    )
    raw_pred_csv = run_dir / "raw_predictions.csv"
    run_chemprop_predict(run_dir / "model_0", predict_csv, raw_pred_csv, logger)

    expected_names = set(population.loc[population["screen_split"] == "screen_test", "Molecule_Name"])
    verify_predictions(raw_pred_csv, expected_names, REGRESSION_ENDPOINTS, logger)
    return pd.read_csv(raw_pred_csv)


def score_and_save(repeat: int, fold: int, seed: int, pred_df: pd.DataFrame, curated: pd.DataFrame, logger) -> None:
    pred_path, score_path = result_paths(repeat, fold)
    PRED_DIR.mkdir(parents=True, exist_ok=True)
    SCORE_DIR.mkdir(parents=True, exist_ok=True)

    ground_truth = curated[curated["inchikey"].isin(pred_df["inchikey"])].copy()
    with per_fold_bootstrap_seed(seed):
        scored = score_activity_predictions(pred_df, ground_truth, SCORE_ENDPOINTS)  # 05c's 2-isoform scope, no macro
    scored["config"] = CONFIG
    scored["repeat"] = repeat
    scored["fold"] = fold
    scored["bootstrap_seed"] = seed

    pred_df.to_csv(pred_path, index=False)
    scored.to_csv(score_path, index=False)
    logger.info(f"wrote {pred_path.name} ({len(pred_df)} rows), {score_path.name} ({len(scored)} rows)")


def main():
    if not FOLDS_PATH.exists():
        raise FileNotFoundError(
            f"{FOLDS_PATH} not found -- run notebook 05b (Sections 1-4) first; 05c reuses its fold file as-is."
        )

    log_path = LOG_DIR / "05c_cluster_cv_randominit.log"
    logger = setup_logging(log_path, "run_05c_cluster_cv_randominit")
    logger.info("=" * 70)
    logger.info("starting scripts/05c_run_cluster_cv_randominit.py")
    logger.info(f"python executable: {sys.executable}")
    logger.info(f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS')}, KMP_DUPLICATE_LIB_OK={os.environ.get('KMP_DUPLICATE_LIB_OK')}")
    logger.info(f"FOLDS_PATH (read-only, reused from 05b): {FOLDS_PATH}")
    logger.info(f"SCORE_ENDPOINTS (scoring-time restriction only, training stays full multitask): {SCORE_ENDPOINTS}")

    seeds = fold_seeds()
    cross_check_seeds_against_notebook_05(seeds, logger)

    curated = pd.read_csv(CURATED_PATH)
    logger.info(f"loaded {CURATED_PATH.name}: {curated.shape}")
    folds_df = pd.read_csv(FOLDS_PATH)
    logger.info(f"loaded {FOLDS_PATH.name}: {folds_df.shape}")

    n_done, n_run, n_failed = 0, 0, 0
    for repeat in range(N_REPEATS):
        for fold in range(N_FOLDS):
            seed = seeds[(repeat, fold)]
            if is_done(repeat, fold):
                n_done += 1
                continue

            logger.info(f"--- {CONFIG} repeat={repeat} fold={fold} seed={seed} ---")
            t0 = time.time()
            try:
                pred_df = run_chemprop(repeat, fold, seed, logger)
                score_and_save(repeat, fold, seed, pred_df, curated, logger)
                n_run += 1
                logger.info(f"{CONFIG} repeat={repeat} fold={fold} done in {time.time() - t0:.1f}s")
            except Exception:
                n_failed += 1
                logger.error(
                    f"{CONFIG} repeat={repeat} fold={fold} FAILED after {time.time() - t0:.1f}s "
                    f"-- no output written, will retry on next invocation:\n{traceback.format_exc()}"
                )

    logger.info(f"finished: {n_done} already done, {n_run} completed this run, {n_failed} failed")


if __name__ == "__main__":
    main()
