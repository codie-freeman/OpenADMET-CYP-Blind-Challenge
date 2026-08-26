"""Train the final activity-track submission model: the 04a baseline screen's winning
config (`chemprop_multitask_chemeleoninit` -- a Chemprop model whose own encoder is
CheMeleon-initialized, NOT the tabular `chemeleon_{rf,xgboost,lightgbm}` configs, which
are unrelated and mediocre), retrained on the FULL labeled training set (all 4,905
curated compounds) instead of the screen's held-out-fold subset, then used to predict
on the real blinded test set.

SAME RECIPE, MORE DATA -- NOTHING ELSE: the Chemprop CLI invocation below is byte-
identical to `run_baseline_screen_chemprop_chemeleon.py`'s own multitask model (reused
directly via `src.chemprop_screen.run_chemprop_train`, not re-typed here) --
`--from-foundation CHEMELEON --multi-hot-atom-featurizer-mode V2 --epochs 50
--patience 5 --data-seed 42 --pytorch-seed 42`. No task-weighting, no architecture
change, no hyperparameter tuning -- any such idea is a distinct future experiment, not
this script. The only two differences from the screen: (1) this trains on all 4,905
labeled compounds rather than the screen's ~3,924-compound pooled screen_inner_train/
screen_inner_val subset -- there is no held-out "test" fold at all during training,
only train vs. an early-stopping validation slice -- and (2) predictions are made on
the REAL blinded test set (`data/processed/test_blinded_curated.csv`, 750 compounds,
confirmed against `ACTIVITY_DATASET_SIZE` in both
`src/vendor/openadmet_eval/config.py` and the vendored
`src/vendor/validation/activity_validation.py`) rather than a held-out training fold.

INTERNAL VALIDATION SPLIT: `assign_final_submission_split` (src/features.py) carves
`VAL_FRACTION` of the 4,905 compounds into a `final_val` slice purely for Chemprop's
`--patience`-based early stopping -- same fraction as the screen's own
`assign_screen_split` (established precedent for this exact purpose), same global
`SEED`. Entirely independent of `data/folds/cv_folds.csv` and the 25-fold CV backbone:
this split is computed only from `train_inhibition_curated.csv`, written to its own
file (`outputs/final_submission/final_train_val_split.csv`) under a
`final_submission_split` column name that cannot be confused with anything
screen_*-named or fold-related, and `cv_folds.csv` itself is never opened here.

VALIDATION MODULE: uses the tutorial repo's own `validate_activity_submission`,
vendored unmodified at `src/vendor/validation/` (see PROVENANCE.md there) rather than
reimplemented -- imported the same way the tutorial's own notebook does
(`from validation.activity_validation import validate_activity_submission`). If
validation fails, this script raises and stops before writing anything past the
submission CSV -- no further step (formatting was already the last one) is attempted.

EXPECT THIS TO TAKE SEVERAL HOURS on this CPU-only hardware -- more training compounds
than the already-multi-hour screen's own multitask model (~101 min for 7 epochs on
3,924 pooled train+val compounds, per `logs/04a_chemprop_screen_chemeleon.log`). An
overnight run is the expected, acceptable way to run this, matching
`run_baseline_screen_chemprop_chemeleon.py`'s own precedent -- do not increase
epochs/patience "to be safe": consistency with the screened recipe matters more here
than squeezing out extra performance.

    cd /path/to/OpenADMET-CYP-Blind-Challenge
    caffeinate -i nohup python scripts/train_final_submission_multitask.py > \\
        logs/train_final_submission_multitask_stdout.log 2>&1 &

OUTPUTS:
    logs/train_final_submission_multitask.log                        -- progress log
    outputs/final_submission/final_train_val_split.csv                -- Molecule_Name,
                                                                           inchikey,
                                                                           final_submission_split
                                                                           (NOT cv_folds.csv)
    outputs/final_submission/chemprop_run/                            -- chemprop's own
                                                                           train_input.csv,
                                                                           checkpoints,
                                                                           config.toml
    outputs/final_submission/predict_input.csv                        -- real blinded
                                                                           test-set predict
                                                                           input (750
                                                                           compounds)
    outputs/final_submission/predictions/
        chemprop_multitask_chemeleoninit_raw_predictions.csv          -- chemprop's raw
                                                                           predict output
                                                                           (Molecule_Name,
                                                                           inchikey,
                                                                           canonical_smiles,
                                                                           + 4 isoform cols)
    outputs/submissions/activity_submission_v1.csv                    -- final, schema-
                                                                           formatted,
                                                                           validated
                                                                           submission (750
                                                                           rows: SMILES,
                                                                           Molecule_Name,
                                                                           + 4 isoform cols)
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import sys
import time
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))  # so `from src... import ...` works regardless of cwd
sys.path.insert(0, str(REPO_ROOT / "src" / "vendor"))  # so `from validation... import ...`
# matches the tutorial's own import line exactly (see src/vendor/validation/PROVENANCE.md)

from src.chemprop_screen import (
    run_chemprop_predict,
    run_chemprop_train,
    setup_logging,
    verify_predictions,
)
from src.features import assign_final_submission_split
from src.vendor.openadmet_eval.config import REGRESSION_ENDPOINTS
from validation.activity_validation import validate_activity_submission

PROCESSED = REPO_ROOT / "data" / "processed"
LOGS = REPO_ROOT / "logs"
FINAL_OUT = REPO_ROOT / "outputs" / "final_submission"
SUBMISSIONS_OUT = REPO_ROOT / "outputs" / "submissions"
CHEMPROP_RUN_DIR = FINAL_OUT / "chemprop_run"
PREDICTIONS_DIR = FINAL_OUT / "predictions"

SEED = 42
VAL_FRACTION = 0.15  # matches assign_screen_split's own precedent for this exact purpose
EPOCHS = 50
PATIENCE = 5
ARCHITECTURE_ARGS = ["--from-foundation", "CHEMELEON", "--multi-hot-atom-featurizer-mode", "V2"]
# ^ byte-identical to run_baseline_screen_chemprop_chemeleon.py's ARCHITECTURE_ARGS --
# confirmed by reading that script and src/chemprop_screen.py directly before writing
# this file, not reconstructed from memory.

LOG_PATH = LOGS / "train_final_submission_multitask.log"
SUBMISSION_PATH = SUBMISSIONS_OUT / "activity_submission_v1.csv"


def main():
    logger = setup_logging(LOG_PATH, "train_final_submission_multitask")
    script_start = time.time()
    logger.info("=" * 70)
    logger.info("starting train_final_submission_multitask.py")
    logger.info(
        "Reproducing the 04a screen's winning config (chemprop_multitask_chemeleoninit) "
        "on the FULL labeled training set (4,905 compounds) -- same architecture/"
        "hyperparameters as the screen, more data, real blinded test-set predictions."
    )
    logger.info(
        "EXPECT SEVERAL HOURS: more training compounds than the screen's own multitask "
        "model (~101 min for 7 epochs on 3,924 compounds) -- an overnight run is "
        "expected, not a sign something's wrong."
    )

    import chemprop
    import lightning
    import torch

    logger.info(f"python: {sys.version.split()[0]}")
    logger.info(f"torch: {torch.__version__}")
    logger.info(f"chemprop: {chemprop.__version__}")
    logger.info(f"lightning: {lightning.__version__}")
    logger.info(f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS')}")
    logger.info(f"seed: {SEED} (--data-seed, --pytorch-seed, and assign_final_submission_split)")
    logger.info(
        f"hyperparameters: epochs={EPOCHS}, patience={PATIENCE}, "
        "from_foundation=CHEMELEON, multi_hot_atom_featurizer_mode=V2 -- identical to "
        "run_baseline_screen_chemprop_chemeleon.py's multitask config"
    )

    FINAL_OUT.mkdir(parents=True, exist_ok=True)
    SUBMISSIONS_OUT.mkdir(parents=True, exist_ok=True)
    PREDICTIONS_DIR.mkdir(parents=True, exist_ok=True)
    CHEMPROP_RUN_DIR.mkdir(parents=True, exist_ok=True)

    # ---- 1. load full labeled training set (ALL 4,905 compounds) -------------------
    curated_path = PROCESSED / "train_inhibition_curated.csv"
    curated = pd.read_csv(curated_path)
    logger.info(f"loaded {curated_path.name}: {curated.shape}")
    logger.info(
        f"rows before/after filtering: {len(curated)}/{len(curated)} (no filtering -- "
        "multitask training keeps every compound; chemprop NaN-masks any missing "
        "per-task target on its own, matching the screen's own multitask config)"
    )

    # ---- 2. internal train/val split for early stopping (NOT cv_folds.csv) ---------
    split_df = assign_final_submission_split(curated, val_fraction=VAL_FRACTION, seed=SEED)
    logger.info(
        f"assign_final_submission_split(val_fraction={VAL_FRACTION}, seed={SEED}) "
        f"split counts: {split_df['final_submission_split'].value_counts().to_dict()}"
    )
    split_out_path = FINAL_OUT / "final_train_val_split.csv"
    split_df.to_csv(split_out_path, index=False)
    logger.info(
        f"wrote {split_out_path} ({len(split_df)} rows) -- independent of "
        "data/folds/cv_folds.csv, which this script never opens"
    )

    population = split_df.merge(
        curated[["Molecule_Name", "canonical_smiles", *REGRESSION_ENDPOINTS]],
        on="Molecule_Name", how="left",
    )
    if len(population) != len(split_df) or population["canonical_smiles"].isna().any():
        raise ValueError(
            "final-submission split compounds did not align 1:1 against "
            "train_inhibition_curated.csv by Molecule_Name -- stopping before training."
        )
    logger.info(f"rows after SMILES/target join: {len(population)} (missing: 0)")

    for endpoint in REGRESSION_ENDPOINTS:
        counts = population.dropna(subset=[endpoint]).groupby("final_submission_split").size().to_dict()
        logger.info(f"  {endpoint}: labeled-compound counts per split = {counts}")

    # ---- 3. build training CSV (multitask, all 4 isoforms, NaN-masked) -------------
    split_map = {"final_train": "train", "final_val": "val"}
    population["chemprop_split"] = population["final_submission_split"].map(split_map)
    if population["chemprop_split"].isna().any():
        raise ValueError("unexpected final_submission_split value(s) -- stopping.")

    train_csv = CHEMPROP_RUN_DIR / "train_input.csv"
    cols = ["canonical_smiles", *REGRESSION_ENDPOINTS, "chemprop_split"]
    population[cols].to_csv(train_csv, index=False)
    logger.info(
        f"wrote {train_csv} ({len(population)} rows, chemprop_split counts: "
        f"{population['chemprop_split'].value_counts().to_dict()})"
    )

    # ---- 4. train (exact same CLI invocation as the screen's multitask config) -----
    run_chemprop_train(
        train_csv, REGRESSION_ENDPOINTS, CHEMPROP_RUN_DIR, logger,
        ARCHITECTURE_ARGS, EPOCHS, PATIENCE, SEED,
    )

    # ---- 5. predict on the REAL blinded test set ------------------------------------
    test_path = PROCESSED / "test_blinded_curated.csv"
    test_df = pd.read_csv(test_path)
    logger.info(f"loaded {test_path.name}: {test_df.shape}")
    if len(test_df) != 750:
        raise ValueError(f"expected 750 blinded test compounds, got {len(test_df)} -- stopping.")

    predict_csv = FINAL_OUT / "predict_input.csv"
    test_df[["Molecule_Name", "inchikey", "canonical_smiles"]].to_csv(predict_csv, index=False)
    logger.info(f"wrote {predict_csv} ({len(test_df)} blinded test compounds)")

    raw_pred_path = PREDICTIONS_DIR / "chemprop_multitask_chemeleoninit_raw_predictions.csv"
    run_chemprop_predict(CHEMPROP_RUN_DIR, predict_csv, raw_pred_path, logger)
    expected_molecule_names = set(test_df["Molecule_Name"])
    verify_predictions(raw_pred_path, expected_molecule_names, REGRESSION_ENDPOINTS, logger)

    # ---- 6. format to the required submission schema --------------------------------
    raw_pred = pd.read_csv(raw_pred_path)
    submission = test_df[["SMILES", "Molecule_Name"]].merge(
        raw_pred[["Molecule_Name", *REGRESSION_ENDPOINTS]], on="Molecule_Name", how="left"
    )
    if len(submission) != len(test_df) or submission[REGRESSION_ENDPOINTS].isna().any().any():
        raise ValueError("submission formatting produced missing rows/values -- stopping.")
    submission = submission[["SMILES", "Molecule_Name", *REGRESSION_ENDPOINTS]]
    submission.to_csv(SUBMISSION_PATH, index=False)
    logger.info(f"wrote {SUBMISSION_PATH} ({len(submission)} rows)")

    reread = pd.read_csv(SUBMISSION_PATH)
    logger.info(f"rows after reading {SUBMISSION_PATH.name} back from disk: {len(reread)}")
    if len(reread) != 750:
        raise ValueError(f"reread {SUBMISSION_PATH.name} has {len(reread)} rows, expected 750 -- stopping.")

    # ---- 7. validate via the tutorial's own validate_activity_submission ------------
    is_valid, validation_errors = validate_activity_submission(
        SUBMISSION_PATH, expected_ids=expected_molecule_names
    )
    if is_valid:
        logger.info("VALIDATION: PASS -- activity submission file is valid.")
    else:
        logger.error("VALIDATION: FAIL -- activity submission file is invalid:")
        for msg in validation_errors:
            logger.error(f"  - {msg}")
        raise ValueError(
            "submission failed validate_activity_submission -- stopping before anything else."
        )

    total_elapsed = time.time() - script_start
    logger.info("=" * 70)
    logger.info(f"total script wall time: {total_elapsed:.1f}s ({total_elapsed / 60:.1f} min)")
    logger.info("done.")


if __name__ == "__main__":
    main()
