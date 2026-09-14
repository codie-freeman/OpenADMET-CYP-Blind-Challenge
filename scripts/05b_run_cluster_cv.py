"""Resumable driver for notebook 05b's cluster-aware CV refit -- `chemprop_chemeleoninit`
only, over the 25 (repeat, fold) fold-partitions in `outputs/05b_cluster_cv_comparison/
cluster_cv_folds.csv` (built by the notebook's Sections 1-4, NOT this script). Directly
mirrors `scripts/run_5x5_cv_comparison.py`'s `run_chemprop`/`score_and_save`/`is_done`
functions and its per-(repeat,fold) seeding scheme -- see that script's own docstring
for the reasoning behind each; not re-explained here except where this script differs.

WHAT DIFFERS FROM `run_5x5_cv_comparison.py`, AND WHY:
- Only one config (`chemprop_chemeleoninit`) over one fold source (`cluster_cv_folds.csv`
  instead of `data/folds/cv_folds.csv`) -- notebook 05b's whole point is that the fold
  assignment is the only variable under test, so everything else (feature representation,
  architecture args, epochs, patience, val_fraction) is copied verbatim from that script's
  own `CHEMPROP_ARCHITECTURE_ARGS["chemprop_chemeleoninit"]` / `CHEMPROP_EPOCHS` /
  `CHEMPROP_PATIENCE` / `VAL_FRACTION` constants, not re-derived.
- Per-(repeat, fold) seeds are regenerated here via the exact same formula as
  `scripts/generate_5x5_cv_manifest.py` (`np.random.SeedSequence(42).spawn(25)`, same
  enumeration order) rather than read from that script's manifest.csv, so this script has
  no file dependency on notebook 05's own run -- and cross-checked once at startup against
  that manifest.csv (read-only) to confirm the two really do match, since an unnoticed
  mismatch here would silently make training-randomness a second variable under test
  alongside the fold assignment, defeating the notebook's whole comparison.
- No manifest/family/shard machinery -- only 25 rows total (vs. 300), run sequentially,
  single process. Resumability is still file-existence-based (`is_done`), same contract.
- All output goes under `outputs/05b_cluster_cv_comparison/` (`predictions/`, `scores/`,
  `chemprop_runs/`) -- `outputs/05_cv_comparison/` is never read or written by this script
  except for the one read-only manifest.csv cross-check above.

ENVIRONMENT: KMP_DUPLICATE_LIB_OK=TRUE set before RDKit/PyTorch/Chemprop import, and this
conda env's own bin/ prepended to PATH before the `chemprop` subprocess call -- both
copied verbatim from `run_5x5_cv_comparison.py` (see that script's docstring for why).
OMP_NUM_THREADS is left unset, same deliberate deviation as that script, same reasoning.

USAGE:
    cd /path/to/OpenADMET-CYP-Blind-Challenge
    caffeinate -i nohup /Users/codiefreeman/miniconda3-arm64/envs/cyp-admet-v2/bin/python \\
        scripts/05b_run_cluster_cv.py > logs/05b_cluster_cv_stdout.log 2>&1 &
    tail -f logs/05b_cluster_cv.log
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
from src.vendor.openadmet_eval.config import ACTIVITY_METRICS, REGRESSION_ENDPOINTS
from src.vendor.openadmet_eval.evaluate_predictions import (
    add_macro_endpoint,
    score_activity_predictions,
)

CONFIG = "chemprop_chemeleoninit"
CURATED_PATH = REPO_ROOT / "data" / "processed" / "train_inhibition_curated.csv"
OUT = REPO_ROOT / "outputs" / "05b_cluster_cv_comparison"
FOLDS_PATH = OUT / "cluster_cv_folds.csv"
ORIGINAL_MANIFEST_PATH = REPO_ROOT / "outputs" / "05_cv_comparison" / "manifest.csv"  # read-only cross-check
PRED_DIR = OUT / "predictions"
SCORE_DIR = OUT / "scores"
CHEMPROP_RUNS_DIR = OUT / "chemprop_runs"
LOG_DIR = REPO_ROOT / "logs"

VAL_FRACTION = 0.15  # matches run_5x5_cv_comparison.py verbatim
CV_SEED_BASE = 42  # matches generate_5x5_cv_manifest.py's CV_SEED_BASE verbatim
N_REPEATS = 5
N_FOLDS = 5

CHEMPROP_ARCHITECTURE_ARGS = ["--from-foundation", "CHEMELEON", "--multi-hot-atom-featurizer-mode", "V2"]
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
    logger.info(
        f"cross-checked all {len(seeds)} seeds against {ORIGINAL_MANIFEST_PATH} ({CONFIG} rows) -- identical."
    )


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
        scored = score_activity_predictions(pred_df, ground_truth, REGRESSION_ENDPOINTS)
        scored = add_macro_endpoint(scored, REGRESSION_ENDPOINTS, ACTIVITY_METRICS)
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
            f"{FOLDS_PATH} not found -- run notebooks/05b_cluster_cv_comparison.ipynb Sections 1-4 first."
        )

    log_path = LOG_DIR / "05b_cluster_cv.log"
    logger = setup_logging(log_path, "run_05b_cluster_cv")
    logger.info("=" * 70)
    logger.info("starting scripts/05b_run_cluster_cv.py")
    logger.info(f"python executable: {sys.executable}")
    logger.info(f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS')}, KMP_DUPLICATE_LIB_OK={os.environ.get('KMP_DUPLICATE_LIB_OK')}")

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
