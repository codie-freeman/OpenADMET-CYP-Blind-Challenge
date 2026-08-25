"""Train the 5 Chemprop configurations for notebook 04a's cheap, single-fold baseline
screen: 4 single-task models (one per CYP isoform, trained independently) and 1
multitask model (all 4 isoforms jointly).

RANDOM-INIT, NOT CHEMELEON-INITIALIZED -- A DELIBERATE, TEMPORARY SWAP FOR THIS SCREEN
ONLY: an earlier version of this script used `--from-foundation CHEMELEON`, matching
the full 25-fold comparison's eventual methodology (CheMeleon-initialized encoders for
both single-task and multitask, so that comparison isolates head structure from
initialization -- see the full comparison's own design, unchanged). That defeated the
cheap screen's whole purpose: CheMeleon's checkpoint forces the message-passing
encoder to its own real architecture (2048-wide, 8.7M params) regardless of
`--epochs`/`--patience`, and a live run confirmed this is infeasible for a fast
screen on this hardware -- single-task CYP1A2 took ~2:47/epoch (measured directly,
`trainer_logs/version_0/metrics.csv` from that run, since deleted per the cleanup
below) against 03b's benchmarked ~13s/epoch for the same-sized random-init network, a
~13x slowdown per epoch. This script now uses the same small, random-init architecture
already validated as fast in `train_log2fc_encoder.py` (`d_h=300, depth=3, aggregation=
mean, batch_norm=True`) instead -- no `--from-foundation` flag at all.

CONSEQUENCE, STATED PLAINLY: this screen's Chemprop numbers are therefore NOT directly
representative of the eventual CheMeleon-initialized model -- a random-init encoder
learns representations from scratch on this screen's few-thousand-compound training
side, which is a different (and likely weaker) starting point than CheMeleon's
pretrained one. Treat these Chemprop results as a check that the training/eval
plumbing works end-to-end, not as a preview of CheMeleon-initialized performance.
`scripts/run_baseline_screen_chemprop_chemeleon.py` now runs the CheMeleon-initialized
variant too, at its real (much slower) measured cost, so both are directly comparable
in the notebook -- see that script and notebook 04a's Section 4 for the comparison.

STANDALONE SCRIPT, NOT A NOTEBOOK CELL: real training runs (not a single forward pass),
so per the project's agreed convention this belongs in a script the person launches
themselves -- e.g.

    cd /path/to/OpenADMET-CYP-Blind-Challenge
    caffeinate -i nohup python scripts/run_baseline_screen_chemprop.py > /dev/null 2>&1 &

-- not inside an agentic session or a notebook cell. `notebooks/04a_baseline_screen.
ipynb` loads and reports on this script's output prediction files; it does not train
anything itself.

THIS IS THE CHEAP SCREEN, NOT THE 25-FOLD COMPARISON: one held-out fold (`repeat_0`,
fold 0) is used as the test set, with the remaining 4 folds pooled as the training side
-- a fast first read on relative performance, not the final reportable result. Fold 0
was the first available option, not chosen after looking at results.

SHARED INFRASTRUCTURE: `src/chemprop_screen.py` holds everything that's identical
between this script and `run_baseline_screen_chemprop_chemeleon.py` -- the
`assign_screen_split` call, training/predict-CSV construction, the per-epoch
metrics.csv-polling progress logger, the subprocess stdin fix, and disk-reread
verification. Only `run_chemprop_train`'s `architecture_args` differ between the two
scripts (see below); everything else is called identically, so both scripts see a
byte-identical screen split and neither can silently drift from the other's
infrastructure.

ARCHITECTURE: `--message-hidden-dim 300 --depth 3 --aggregation mean --batch-norm` --
identical to `train_log2fc_encoder.py`'s validated-fast configuration (that script's
own Python-API equivalents: `BondMessagePassing(d_h=300, depth=3)`, `MeanAggregation()`,
`MPNN(..., batch_norm=True)`). `d_h=300`/`depth=3` already match chemprop's own CLI
defaults; `aggregation=mean` and `--batch-norm` do not (CLI defaults are `norm` and
off) and are passed explicitly to fully match 03b's architecture, not just its two
most-cited numbers. No `--multi-hot-atom-featurizer-mode` override is needed here --
that flag's CheMeleon-specific hard requirement (`must be used with ... V2`) is
enforced only inside chemprop's own `--from-foundation` code path (confirmed directly
in `chemprop/cli/train.py`), which never executes without that flag; chemprop's plain
CLI default (`V2`) applies unchanged.

EARLY STOPPING: `--patience 5` (of a 50-epoch cap) -- confirmed with the user, since the
task only specified that patience must be set explicitly (it defaults to None, silently
disabling early stopping) without a number. `--epochs 50` matches chemprop's own CLI
default and this project's established convention (`train_log2fc_encoder.py`).
Multitask NaN handling needs no special code: chemprop masks missing per-task targets
in its loss function automatically (verified directly in `chemprop/nn/metrics.py`), so
a compound measured on only some isoforms simply contributes to those isoforms' loss.

TRAINING/PREDICTION SPLIT: each model's training CSV contains only
`screen_inner_train` + `screen_inner_val` compounds (mapped to chemprop's expected
`train`/`val` values via `--splits-column`) -- `screen_test` is deliberately never
included in a training CSV. A single shared prediction CSV (`predict_input.csv`,
`chemprop_runs/`) covers *every* `screen_test` compound's SMILES regardless of
per-isoform label availability, decoupling "has features" from "has this isoform's
label" -- exactly matching the tabular configs' convention in the notebook. Each
single-task model's training CSV additionally drops rows lacking that one isoform's
label (logged before/after, per CLAUDE.md) -- a masked-away row contributes nothing to a
single-task loss anyway.

SEEDS: `SEED = 42` (project convention) for both `--data-seed` (unused here beyond
dataloader shuffling, since splits are pre-assigned via `--splits-column`, not
chemprop's own splitter) and `--pytorch-seed` (weight init / shuffling) -- logged below.

OUTPUTS:
    logs/04a_chemprop_screen.log                                    -- progress log
    outputs/04a_baseline_screen/chemprop_runs/predict_input.csv     -- shared predict input
    outputs/04a_baseline_screen/chemprop_runs/{model}/              -- chemprop's own
                                                                        train_input.csv,
                                                                        checkpoints,
                                                                        config.toml
    outputs/04a_baseline_screen/predictions/chemprop_singletask_{ISOFORM}.csv  (x4)
    outputs/04a_baseline_screen/predictions/chemprop_multitask.csv
        -- each: Molecule_Name, inchikey, canonical_smiles, and one column per isoform
           the model covers (named identically to `REGRESSION_ENDPOINTS`, e.g.
           `CYP1A2_pIC50_direct_inhibition`), for all 981 `screen_test` compounds.
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
D_H = 300  # matches train_log2fc_encoder.py's validated-fast random-init architecture
DEPTH = 3  # (also chemprop's own CLI defaults for these two -- see module docstring)
ARCHITECTURE_ARGS = ["--message-hidden-dim", str(D_H), "--depth", str(DEPTH), "--aggregation", "mean", "--batch-norm"]

LOG_PATH = LOGS / "04a_chemprop_screen.log"


def main():
    logger = setup_logging(LOG_PATH, "run_baseline_screen_chemprop")
    script_start = time.time()
    logger.info("=" * 70)
    logger.info("starting run_baseline_screen_chemprop.py")
    logger.info(
        "THIS IS THE CHEAP, SINGLE-FOLD BASELINE SCREEN -- NOT the 25-fold (5x5) "
        "comparison. Fold 0 of repeat_0 was the first available option, not chosen "
        "after looking at results."
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
        f"hyperparameters: epochs={EPOCHS}, patience={PATIENCE}, d_h={D_H}, depth={DEPTH}, "
        "aggregation=mean, batch_norm=True (random-init, NOT CheMeleon-initialized -- see "
        "module docstring for why this screen deliberately differs from the full 25-fold "
        "comparison's methodology here)"
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
        model_name = f"singletask_{short_name}"
        model_start = time.time()
        logger.info("-" * 70)
        logger.info(f"model: {model_name} (target: {endpoint})")

        model_dir = CHEMPROP_RUNS / model_name
        shutil.rmtree(model_dir, ignore_errors=True)  # fresh dir -- keeps trainer_logs at version_0
        model_dir.mkdir(parents=True, exist_ok=True)
        train_csv = model_dir / "train_input.csv"
        build_training_csv(population, [endpoint], train_csv, logger, require_all_targets=True)

        run_chemprop_train(train_csv, [endpoint], model_dir, logger, ARCHITECTURE_ARGS, EPOCHS, PATIENCE, SEED)

        pred_out = PREDICTIONS / f"chemprop_{model_name}.csv"
        run_chemprop_predict(model_dir, predict_csv, pred_out, logger)
        verify_predictions(pred_out, expected_molecule_names, [endpoint], logger)

        logger.info(f"model {model_name} elapsed: {time.time() - model_start:.1f}s")

    # ---- 3. one multitask model ------------------------------------------------
    model_name = "multitask"
    model_start = time.time()
    logger.info("-" * 70)
    logger.info(f"model: {model_name} (targets: {REGRESSION_ENDPOINTS})")

    model_dir = CHEMPROP_RUNS / model_name
    shutil.rmtree(model_dir, ignore_errors=True)  # fresh dir -- keeps trainer_logs at version_0
    model_dir.mkdir(parents=True, exist_ok=True)
    train_csv = model_dir / "train_input.csv"
    build_training_csv(population, REGRESSION_ENDPOINTS, train_csv, logger, require_all_targets=False)

    run_chemprop_train(train_csv, REGRESSION_ENDPOINTS, model_dir, logger, ARCHITECTURE_ARGS, EPOCHS, PATIENCE, SEED)

    pred_out = PREDICTIONS / f"chemprop_{model_name}.csv"
    run_chemprop_predict(model_dir, predict_csv, pred_out, logger)
    verify_predictions(pred_out, expected_molecule_names, REGRESSION_ENDPOINTS, logger)

    logger.info(f"model {model_name} elapsed: {time.time() - model_start:.1f}s")

    total_elapsed = time.time() - script_start
    logger.info("=" * 70)
    logger.info(f"total script wall time: {total_elapsed:.1f}s ({total_elapsed / 60:.1f} min)")
    logger.info("done.")


if __name__ == "__main__":
    main()
