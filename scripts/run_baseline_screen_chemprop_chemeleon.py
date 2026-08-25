"""Train the CheMeleon-initialized counterpart of notebook 04a's 5 Chemprop
configurations: 4 single-task models (one per CYP isoform, trained independently) and
1 multitask model (all 4 isoforms jointly) -- alongside, not replacing,
`run_baseline_screen_chemprop.py`'s random-init results. Both stay in the notebook,
clearly labeled, for direct comparison (score AND measured training time).

WHY THIS SCRIPT EXISTS: `run_baseline_screen_chemprop.py`'s own docstring already
explains why this screen's Chemprop configs use a fast random-init stand-in rather
than CheMeleon initialization (CheMeleon's checkpoint forces the encoder to its own
real architecture, 2048-wide/8.7M params, regardless of `--epochs`/`--patience` --
measured directly at ~2:47/epoch, infeasible for that screen's original "fast, seconds
per fit" framing). Random-init is a stand-in used to unblock the cheap screen quickly;
CheMeleon-init is the architecture the project's actual plan (the full 25-fold
comparison) uses. This script runs that real architecture, once, at its real measured
cost, so the two can be compared directly rather than left as an assumption.

TIGHTLY SCOPED -- ONE ADDITION, NOT A RE-DO: same fold, same split, same 5 model
shapes (4 single-task + 1 multitask), same `--epochs 50 --patience 5
--data-seed 42 --pytorch-seed 42` as every other Chemprop config in this screen -- no
extra seeds, no hyperparameter variations, no re-running anything that already exists.
The only two differences from `run_baseline_screen_chemprop.py` are the encoder's
initialization/architecture (see ARCHITECTURE below) and the output paths (see OUTPUTS
below, chosen so nothing collides with the random-init run's own files).

SHARED INFRASTRUCTURE: `src/chemprop_screen.py` holds everything identical between
this script and `run_baseline_screen_chemprop.py` -- the `assign_screen_split` call
(same arguments, so both scripts see a byte-identical screen split), training/
predict-CSV construction, the per-epoch metrics.csv-polling progress logger, the
subprocess stdin fix, and disk-reread verification. Only `ARCHITECTURE_ARGS` below
differs between the two scripts' `run_chemprop_train` calls.

ARCHITECTURE: `--from-foundation CHEMELEON --multi-hot-atom-featurizer-mode V2`.
CheMeleon's checkpoint (already cached at `~/.chemprop/chemeleon_mp.pt`, so this
triggers no network call) fixes the message-passing encoder's own architecture
(2048-wide, 8.7M params) -- `--message-hidden-dim`/`--depth`/`--aggregation`/
`--batch-norm` are inapplicable here and are not passed (chemprop's `--from-foundation`
code path auto-configures the encoder from the checkpoint's own hyperparameters,
verified directly in `chemprop/cli/train.py`). `--multi-hot-atom-featurizer-mode V2` is
CheMeleon's own hard requirement (chemprop raises otherwise). Single-task and
multitask models get the *same* CheMeleon starting weights for the encoder -- only the
FFN head differs (1 output vs. 4) -- so the single-task-vs-multitask comparison
isolates head structure, not initialization, matching the full 25-fold comparison's
own design.

EXPECT THIS TO TAKE MULTIPLE HOURS -- AN OVERNIGHT RUN IS THE EXPECTED, ACCEPTABLE WAY
TO RUN THIS, NOT A SIGN SOMETHING'S WRONG: today's aborted first attempt at
CheMeleon-init for this screen measured ~2:47/epoch for single-task CYP1A2 (969
training compounds) before being interrupted. At up to 50 epochs (patience=5 may stop
earlier) across 4 single-task models (roughly 900-1900 pooled-train compounds each)
plus 1 multitask model (3,924 compounds, so proportionally slower per epoch), total
wall time is expected to run into multiple hours on this CPU-only hardware. Launch it,
walk away, and check back -- do not increase epochs/patience "to be safe": consistency
with every other Chemprop config in this screen matters more than squeezing extra
performance out of one config, and a longer run only pushes the compute-budget
question further out rather than answering it.

    cd /path/to/OpenADMET-CYP-Blind-Challenge
    caffeinate -i nohup python scripts/run_baseline_screen_chemprop_chemeleon.py > \\
        logs/04a_chemprop_chemeleon_stdout.log 2>&1 &

Per-model AND total wall-clock time are logged explicitly (see OUTPUTS) -- this is real
data for the compute-budget conversation (does CheMeleon-init's measured cost buy
enough of a score improvement over random-init to justify it at the full 25-fold
scale?), worth having on record regardless of how the scores themselves turn out.
Model-name keys in the elapsed-time log lines (`singletask_CYP1A2`, ..., `multitask`)
are deliberately identical to `run_baseline_screen_chemprop.py`'s own, so the two
scripts' log files can be compared line-for-line by the notebook -- only the output
file/directory paths are prefixed to avoid collisions (see OUTPUTS).

STILL THE CHEAP SCREEN'S SINGLE FOLD, NOT THE 25-FOLD COMPARISON: same `repeat_0`
fold-0 held-out test set as every other config in this screen. The wall-clock time
measured here is a same-fold, same-hardware data point for the compute-budget
question -- not yet a full-25-fold projection.

OUTPUTS:
    logs/04a_chemprop_screen_chemeleon.log                                  -- progress log
    outputs/04a_baseline_screen/chemprop_runs/predict_input.csv             -- shared with
                                                                                the random-init
                                                                                script (identical
                                                                                content, rebuilt
                                                                                here independently)
    outputs/04a_baseline_screen/chemprop_runs/chemeleon_{model}/            -- chemprop's own
                                                                                train_input.csv,
                                                                                checkpoints,
                                                                                config.toml
                                                                                ("chemeleon_"
                                                                                prefix avoids
                                                                                colliding with the
                                                                                random-init run's
                                                                                own model dirs)
    outputs/04a_baseline_screen/predictions/chemprop_chemeleon_singletask_{ISOFORM}.csv  (x4)
    outputs/04a_baseline_screen/predictions/chemprop_chemeleon_multitask.csv
        -- each: Molecule_Name, inchikey, canonical_smiles, and one column per isoform
           the model covers (named identically to `REGRESSION_ENDPOINTS`), for all 981
           `screen_test` compounds.
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import shutil
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))  # so `from src... import ...` works regardless of cwd

from src.chemprop_screen import (
    ISOFORM_SHORT_NAMES,
    build_predict_csv,
    build_training_csv,
    load_screen_population,
    run_chemprop_predict,
    run_chemprop_train,
    setup_logging,
    verify_predictions,
)
from src.vendor.openadmet_eval.config import REGRESSION_ENDPOINTS

FOLDS = REPO_ROOT / "data" / "folds"
PROCESSED = REPO_ROOT / "data" / "processed"
LOGS = REPO_ROOT / "logs"
OUTPUTS = REPO_ROOT / "outputs" / "04a_baseline_screen"
CHEMPROP_RUNS = OUTPUTS / "chemprop_runs"
PREDICTIONS = OUTPUTS / "predictions"

SEED = 42
REPEAT_COL = "repeat_0"
TEST_FOLD = 0
VAL_FRACTION = 0.15
EPOCHS = 50
PATIENCE = 5
ARCHITECTURE_ARGS = ["--from-foundation", "CHEMELEON", "--multi-hot-atom-featurizer-mode", "V2"]

LOG_PATH = LOGS / "04a_chemprop_screen_chemeleon.log"


def main():
    logger = setup_logging(LOG_PATH, "run_baseline_screen_chemprop_chemeleon")
    script_start = time.time()
    logger.info("=" * 70)
    logger.info("starting run_baseline_screen_chemprop_chemeleon.py")
    logger.info(
        "THIS IS THE CHEAP, SINGLE-FOLD BASELINE SCREEN -- NOT the 25-fold (5x5) "
        "comparison. Same repeat_0/fold-0 split as run_baseline_screen_chemprop.py "
        "(random-init) -- see that script's own predictions for direct comparison."
    )
    logger.info(
        "EXPECT MULTIPLE HOURS: CheMeleon's checkpoint forces the real 2048-wide/"
        "8.7M-param encoder, measured at ~2:47/epoch for single-task CYP1A2 on a "
        "prior attempt -- an overnight run is expected, not a sign something's wrong."
    )

    import chemprop
    import lightning
    import torch

    logger.info(f"python: {sys.version.split()[0]}")
    logger.info(f"torch: {torch.__version__}")
    logger.info(f"chemprop: {chemprop.__version__}")
    logger.info(f"lightning: {lightning.__version__}")
    logger.info(f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS')}")
    logger.info(f"seed: {SEED} (--data-seed and --pytorch-seed)")
    logger.info(
        f"hyperparameters: epochs={EPOCHS}, patience={PATIENCE}, "
        "from_foundation=CHEMELEON, multi_hot_atom_featurizer_mode=V2 -- identical "
        "epochs/patience/seeds to run_baseline_screen_chemprop.py, architecture fixed "
        "by the CheMeleon checkpoint (see module docstring)"
    )

    PREDICTIONS.mkdir(parents=True, exist_ok=True)
    CHEMPROP_RUNS.mkdir(parents=True, exist_ok=True)

    # ---- 1. screen population + shared predict input -----------------------------
    population = load_screen_population(
        FOLDS / "cv_folds.csv", PROCESSED / "train_inhibition_curated.csv",
        REPEAT_COL, TEST_FOLD, VAL_FRACTION, SEED, logger,
    )
    predict_csv = build_predict_csv(population, CHEMPROP_RUNS, logger)
    expected_molecule_names = set(
        population.loc[population["screen_split"] == "screen_test", "Molecule_Name"]
    )

    # ---- 2. four single-task models -----------------------------------------------
    for endpoint, short_name in zip(REGRESSION_ENDPOINTS, ISOFORM_SHORT_NAMES):
        model_name = f"singletask_{short_name}"  # log key -- matches the random-init script's own
        model_start = time.time()
        logger.info("-" * 70)
        logger.info(f"model: {model_name} (target: {endpoint})")

        model_dir = CHEMPROP_RUNS / f"chemeleon_{model_name}"
        shutil.rmtree(model_dir, ignore_errors=True)  # fresh dir -- keeps trainer_logs at version_0
        model_dir.mkdir(parents=True, exist_ok=True)
        train_csv = model_dir / "train_input.csv"
        build_training_csv(population, [endpoint], train_csv, logger, require_all_targets=True)

        run_chemprop_train(train_csv, [endpoint], model_dir, logger, ARCHITECTURE_ARGS, EPOCHS, PATIENCE, SEED)

        pred_out = PREDICTIONS / f"chemprop_chemeleon_{model_name}.csv"
        run_chemprop_predict(model_dir, predict_csv, pred_out, logger)
        verify_predictions(pred_out, expected_molecule_names, [endpoint], logger)

        logger.info(f"model {model_name} elapsed: {time.time() - model_start:.1f}s")

    # ---- 3. one multitask model ------------------------------------------------
    model_name = "multitask"  # log key -- matches the random-init script's own
    model_start = time.time()
    logger.info("-" * 70)
    logger.info(f"model: {model_name} (targets: {REGRESSION_ENDPOINTS})")

    model_dir = CHEMPROP_RUNS / f"chemeleon_{model_name}"
    shutil.rmtree(model_dir, ignore_errors=True)  # fresh dir -- keeps trainer_logs at version_0
    model_dir.mkdir(parents=True, exist_ok=True)
    train_csv = model_dir / "train_input.csv"
    build_training_csv(population, REGRESSION_ENDPOINTS, train_csv, logger, require_all_targets=False)

    run_chemprop_train(train_csv, REGRESSION_ENDPOINTS, model_dir, logger, ARCHITECTURE_ARGS, EPOCHS, PATIENCE, SEED)

    pred_out = PREDICTIONS / f"chemprop_chemeleon_{model_name}.csv"
    run_chemprop_predict(model_dir, predict_csv, pred_out, logger)
    verify_predictions(pred_out, expected_molecule_names, REGRESSION_ENDPOINTS, logger)

    logger.info(f"model {model_name} elapsed: {time.time() - model_start:.1f}s")

    total_elapsed = time.time() - script_start
    logger.info("=" * 70)
    logger.info(f"total script wall time: {total_elapsed:.1f}s ({total_elapsed / 60:.1f} min)")
    logger.info("done.")


if __name__ == "__main__":
    main()
