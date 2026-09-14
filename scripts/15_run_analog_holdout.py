"""Driver for notebook 15's analog-holdout comparison -- `chemprop_chemeleoninit` ONLY,
trained once on the single train/test split built by notebook 15's Part 1 cells
(`data/folds/analog_holdout_split.csv`, NOT repeated CV, NOT `data/folds/cv_folds.csv` or
05b's `cluster_cv_folds.csv` -- both of those are read-only elsewhere in this project and
untouched here).

WHAT THIS IS: a single train/predict/score run, structurally like one fold of
`scripts/05b_run_cluster_cv.py` (same `src/chemprop_screen.py` helpers, same architecture
args/epochs/patience) but with no repeat/fold loop -- this split is one seed set of
compounds held out once, not resampled. Matches this project's own precedent for one-off,
non-CV final runs (`scripts/train_final_submission_multitask.py`,
`scripts/10_final_retrain_predict.py`): SEED = 42 throughout (inner train/val carve-out
for early stopping, chemprop `--data-seed`/`--pytorch-seed`, and the bootstrap evaluator),
not one of the per-(repeat,fold) varying manifest seeds -- there is no fold axis here to
vary across.

VAL_FRACTION = 0.15 for the internal early-stopping carve-out, matching every other
Chemprop run in this project (`run_5x5_cv_comparison.py`, `05b_run_cluster_cv.py`,
`train_final_submission_multitask.py`) -- not re-derived.

Reuses `src/chemprop_screen.py`'s `build_predict_csv`/`build_training_csv`/
`run_chemprop_train`/`run_chemprop_predict`/`verify_predictions` unchanged (same
functions 05b's script calls) -- this script only supplies a differently-shaped
`population` DataFrame (`screen_split` values `screen_test`/`screen_inner_train`/
`screen_inner_val`, built from `analog_holdout_split.csv`'s `split` column plus a plain
`train_test_split` carve-out of the `train` rows, rather than from `assign_screen_split`
over `cv_folds.csv`).

All output goes under `outputs/15_analog_holdout_comparison/` (`predictions/`, `scores/`,
`chemprop_runs/`). Resumability is file-existence-based (`is_done`), same contract as
every other training driver in this project.

ENVIRONMENT: KMP_DUPLICATE_LIB_OK=TRUE set before RDKit/PyTorch/Chemprop import, and this
conda env's own bin/ prepended to PATH before the `chemprop` subprocess call -- copied
verbatim from `run_5x5_cv_comparison.py`/`05b_run_cluster_cv.py` (see those scripts'
docstrings for why).

USAGE:
    cd /path/to/OpenADMET-CYP-Blind-Challenge
    caffeinate -i nohup /Users/codiefreeman/miniconda3-arm64/envs/cyp-admet-v2/bin/python \\
        scripts/15_run_analog_holdout.py > logs/15_analog_holdout_stdout.log 2>&1 &
    tail -f logs/15_analog_holdout.log
"""

import os

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

_env_bin = str(Path(sys.executable).parent)
_path_parts = os.environ.get("PATH", "").split(os.pathsep)
if _env_bin not in _path_parts:
    os.environ["PATH"] = os.pathsep.join([_env_bin, *_path_parts])

from src.chemprop_screen import (
    build_predict_csv,
    build_training_csv,
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
SPLIT_PATH = REPO_ROOT / "data" / "folds" / "analog_holdout_split.csv"
OUT = REPO_ROOT / "outputs" / "15_analog_holdout_comparison"
PRED_DIR = OUT / "predictions"
SCORE_DIR = OUT / "scores"
CHEMPROP_RUNS_DIR = OUT / "chemprop_runs"
LOG_DIR = REPO_ROOT / "logs"

VAL_FRACTION = 0.15  # matches every other Chemprop run in this project verbatim
SEED = 42  # matches this project's one-off final-run precedent (04b, train_final_submission_multitask.py, 10_final_retrain_predict.py)

CHEMPROP_ARCHITECTURE_ARGS = ["--from-foundation", "CHEMELEON", "--multi-hot-atom-featurizer-mode", "V2"]
CHEMPROP_EPOCHS = 50
CHEMPROP_PATIENCE = 5


def result_paths() -> tuple[Path, Path]:
    tag = f"{CONFIG}__analog_holdout"
    return PRED_DIR / f"{tag}.csv", SCORE_DIR / f"{tag}.csv"


def is_done() -> bool:
    pred_path, score_path = result_paths()
    if not (pred_path.exists() and score_path.exists()):
        return False
    try:
        return len(pd.read_csv(pred_path)) > 0 and len(pd.read_csv(score_path)) > 0
    except (pd.errors.EmptyDataError, OSError):
        return False


def build_population(logger) -> pd.DataFrame:
    """Mirrors `load_screen_population`'s contract (same `screen_split` value set,
    same row-count logging), but builds it from `analog_holdout_split.csv`'s `split`
    column instead of `assign_screen_split` over `cv_folds.csv`.
    """
    split_df = pd.read_csv(SPLIT_PATH)
    curated = pd.read_csv(CURATED_PATH)
    logger.info(f"loaded {SPLIT_PATH.name}: {split_df.shape}")
    logger.info(f"loaded {CURATED_PATH.name}: {curated.shape}")
    logger.info(f"analog_holdout_split.csv split counts: {split_df['split'].value_counts().to_dict()}")

    train_pool = split_df[split_df["split"] == "train"].copy()
    test_pool = split_df[split_df["split"] == "test"].copy()

    inner_train_ik, inner_val_ik = train_test_split(
        train_pool["inchikey"].to_numpy(), test_size=VAL_FRACTION, random_state=SEED
    )
    inner_val_set = set(inner_val_ik)
    train_pool["screen_split"] = np.where(
        train_pool["inchikey"].isin(inner_val_set), "screen_inner_val", "screen_inner_train"
    )
    test_pool["screen_split"] = "screen_test"
    population = pd.concat(
        [train_pool[["Molecule_Name", "inchikey", "screen_split"]],
         test_pool[["Molecule_Name", "inchikey", "screen_split"]]],
        ignore_index=True,
    )
    logger.info(
        f"train_test_split(val_fraction={VAL_FRACTION}, seed={SEED}) on the {len(train_pool)} "
        f"analog-holdout train compounds -> screen_split counts: "
        f"{population['screen_split'].value_counts().to_dict()}"
    )

    merged = population.merge(
        curated[["inchikey", "canonical_smiles", *REGRESSION_ENDPOINTS]], on="inchikey", how="left"
    )
    if len(merged) != len(population) or merged["canonical_smiles"].isna().any():
        raise ValueError(
            "analog-holdout split compounds did not align 1:1 against "
            "train_inhibition_curated.csv by inchikey -- stopping before training anything."
        )
    logger.info(f"rows after SMILES/target join: {len(merged)} (missing: 0)")

    for endpoint in REGRESSION_ENDPOINTS:
        counts = merged.dropna(subset=[endpoint]).groupby("screen_split").size().to_dict()
        logger.info(f"  {endpoint}: labeled-compound counts per split = {counts}")

    return merged


def run_chemprop(population: pd.DataFrame, logger) -> pd.DataFrame:
    run_dir = CHEMPROP_RUNS_DIR / f"{CONFIG}__analog_holdout"
    run_dir.mkdir(parents=True, exist_ok=True)

    predict_csv = build_predict_csv(population, run_dir, logger)
    train_csv = run_dir / "train_input.csv"
    build_training_csv(population, REGRESSION_ENDPOINTS, train_csv, logger, require_all_targets=False)
    run_chemprop_train(
        train_csv, REGRESSION_ENDPOINTS, run_dir, logger,
        CHEMPROP_ARCHITECTURE_ARGS, CHEMPROP_EPOCHS, CHEMPROP_PATIENCE, SEED,
    )
    raw_pred_csv = run_dir / "raw_predictions.csv"
    run_chemprop_predict(run_dir / "model_0", predict_csv, raw_pred_csv, logger)

    expected_names = set(population.loc[population["screen_split"] == "screen_test", "Molecule_Name"])
    verify_predictions(raw_pred_csv, expected_names, REGRESSION_ENDPOINTS, logger)
    return pd.read_csv(raw_pred_csv)


def score_and_save(pred_df: pd.DataFrame, curated: pd.DataFrame, logger) -> None:
    pred_path, score_path = result_paths()
    PRED_DIR.mkdir(parents=True, exist_ok=True)
    SCORE_DIR.mkdir(parents=True, exist_ok=True)

    ground_truth = curated[curated["inchikey"].isin(pred_df["inchikey"])].copy()
    with per_fold_bootstrap_seed(SEED):
        scored = score_activity_predictions(pred_df, ground_truth, REGRESSION_ENDPOINTS)
        scored = add_macro_endpoint(scored, REGRESSION_ENDPOINTS, ACTIVITY_METRICS)
    scored["config"] = CONFIG
    scored["split"] = "analog_holdout"
    scored["bootstrap_seed"] = SEED

    pred_df.to_csv(pred_path, index=False)
    scored.to_csv(score_path, index=False)
    logger.info(f"wrote {pred_path.name} ({len(pred_df)} rows), {score_path.name} ({len(scored)} rows)")


def main():
    if not SPLIT_PATH.exists():
        raise FileNotFoundError(
            f"{SPLIT_PATH} not found -- run notebooks/15_analog_holdout_comparison.ipynb Part 1 first."
        )

    log_path = LOG_DIR / "15_analog_holdout.log"
    logger = setup_logging(log_path, "run_15_analog_holdout")
    logger.info("=" * 70)
    logger.info("starting scripts/15_run_analog_holdout.py")
    logger.info(f"python executable: {sys.executable}")
    logger.info(f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS')}, KMP_DUPLICATE_LIB_OK={os.environ.get('KMP_DUPLICATE_LIB_OK')}")

    if is_done():
        logger.info(f"{CONFIG}__analog_holdout already done -- nothing to do.")
        return

    curated = pd.read_csv(CURATED_PATH)
    logger.info(f"loaded {CURATED_PATH.name}: {curated.shape}")
    population = build_population(logger)

    logger.info(f"--- {CONFIG} analog_holdout seed={SEED} ---")
    t0 = time.time()
    try:
        pred_df = run_chemprop(population, logger)
        score_and_save(pred_df, curated, logger)
        logger.info(f"{CONFIG} analog_holdout done in {time.time() - t0:.1f}s")
    except Exception:
        logger.error(
            f"{CONFIG} analog_holdout FAILED after {time.time() - t0:.1f}s -- no output "
            f"written, will retry on next invocation:\n{traceback.format_exc()}"
        )
        raise


if __name__ == "__main__":
    main()
