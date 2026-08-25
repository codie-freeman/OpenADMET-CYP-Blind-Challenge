"""Retrain the log2fc Chemprop encoder (`train_log2fc_encoder.py`, notebook 03b) with
`screen_test` compounds excluded from pretraining -- a leakage check for notebook 04a's
baseline screen, not a replacement for 03b's own frozen artifacts.

WHY THIS SCRIPT EXISTS: 04a's screen found `log2fc_pretrained_{rf,xgboost,lightgbm}`
scoring best ST-RAE of all 12 tabular configs, on every isoform -- suspiciously good.
`train_log2fc_encoder.py` pretrained its encoder on the FULL 4,376-compound
primary-screen set, which includes compounds now sitting in 04a's `screen_test` split
(confirmed directly: 872 of the 4,376 log2fc compounds are `screen_test` compounds, by
`assign_screen_split` -- see `load_training_data` below). The downstream RF/XGBoost/
LightGBM fit in 04a is properly blind (fit only on `screen_inner_train`/`screen_inner_
val`), but the *feature itself* (the frozen embedding) may carry information leaked
from `screen_test` compounds through the encoder's own pretraining -- this script
tests that directly by retraining the identical architecture with those 872 compounds
held out of pretraining entirely, then re-running the same three tabular algorithms on
the resulting embeddings (`log2fc_pretrained_blind_{rf,xgboost,lightgbm}` in notebook
04a) for direct comparison against the original.

EVERYTHING BUT THE TRAINING POPULATION IS IDENTICAL TO `train_log2fc_encoder.py`: same
architecture (`d_h=300, depth=3, aggregation=mean, batch_norm=True`), same `SEED = 42`,
same `--epochs 50` (no validation split here either -- see that script's own EPOCHS
rationale, unchanged), same embedding-extraction step (frozen forward pass over the
full 5,655-compound index, `screen_test` compounds included -- excluding them from
*embedding extraction* would leave 04a with no features to predict on for its held-out
fold; only *pretraining supervision* excludes them, exactly mirroring how 04a's own
RF/XGBoost/LightGBM models are fit train-only but predict on the full test set).

TRAINING POPULATION -- computed directly via `src.features.assign_screen_split` (same
function, same default arguments as 04a's notebook and Chemprop script), not assumed:
of the 4,376 log2fc compounds, 3,504 are `screen_inner_train` (2,975) or `screen_inner_
val` (529) and are kept; the remaining 872 are `screen_test` and are dropped from this
run's training set entirely (logged before/after, per CLAUDE.md).

STANDALONE SCRIPT, NOT A NOTEBOOK CELL, NOT INSIDE AN AGENTIC SESSION: real training
(not a single forward pass), per this project's established convention -- e.g.

    cd /path/to/OpenADMET-CYP-Blind-Challenge
    caffeinate -i nohup python scripts/train_log2fc_encoder_blind.py > /dev/null 2>&1 &

`notebooks/04a_baseline_screen.ipynb` loads and reports on this script's output
artifacts; it does not run any training itself. At ~13s/epoch (03b's own benchmark) on
~3,504 compounds (80% of 03b's 4,376), expect roughly 03b's own ~11-minute runtime or
somewhat less.

DOES NOT TOUCH 03b's OWN FROZEN ARTIFACTS: `models/log2fc_encoder_mp.pt` and
`data/processed/log2fc_pretrained_embeddings*` (03b's originals) are left completely
untouched -- this script writes to `*_blind` -suffixed paths only (see OUTPUTS below),
so both the original (leakage-suspect) and blind (leakage-checked) versions exist
side by side for comparison, per notebook 04a's design.

OUTPUTS:
    logs/04a_train_log2fc_encoder_blind.log                      -- progress log
    models/log2fc_encoder_blind_mp.pt                             -- frozen encoder checkpoint
    data/processed/log2fc_pretrained_blind_embeddings.npy         -- (5655, 300) float32
    data/processed/log2fc_pretrained_blind_embeddings_index.csv   -- matching compound index
                                                                      (Molecule_Name, inchikey,
                                                                      split), same row order
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import logging
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))  # so `from src.features import ...` works regardless of cwd

import numpy as np
import pandas as pd
import torch
from lightning.pytorch.callbacks import Callback
from sklearn.preprocessing import StandardScaler

from src.features import assign_screen_split, frozen_encoder_embeddings

RAW = REPO_ROOT / "data" / "raw"
PROCESSED = REPO_ROOT / "data" / "processed"
FOLDS = REPO_ROOT / "data" / "folds"
MODELS = REPO_ROOT / "models"
LOGS = REPO_ROOT / "logs"

ISOFORMS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]

SEED = 42
REPEAT_COL = "repeat_0"  # same screen-split arguments as notebooks/04a_baseline_screen.ipynb
TEST_FOLD = 0            # and scripts/run_baseline_screen_chemprop.py, for an identical split
VAL_FRACTION = 0.15
N_EPOCHS = 50  # chemprop's own CLI default -- matches train_log2fc_encoder.py, see its docstring
BATCH_SIZE = 64  # chemprop's own default (`build_dataloader`'s default)
D_H = 300  # identical architecture to train_log2fc_encoder.py -- see that script's own rationale
DEPTH = 3
EMBED_BATCH_SIZE = 128  # matches chemeleon_embeddings' / train_log2fc_encoder.py's default

LOG_PATH = LOGS / "04a_train_log2fc_encoder_blind.log"


def setup_logging() -> logging.Logger:
    LOGS.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("train_log2fc_encoder_blind")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    file_handler = logging.FileHandler(LOG_PATH, mode="a")
    file_handler.setFormatter(fmt)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(fmt)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def load_training_data(logger: logging.Logger):
    """Load log2fc_targets.csv (notebook 03 Section 7), attach canonical_smiles (per
    notebook 01's InChIKey identity convention, identical to train_log2fc_encoder.py),
    then -- the one real difference from that script -- drop every compound whose
    `screen_split` (via `assign_screen_split`, same arguments as 04a's notebook/Chemprop
    script) is `screen_test`, so this encoder never sees 04a's held-out fold during
    pretraining. Row counts before/after this filter are logged explicitly (CLAUDE.md).
    """
    log2fc_path = PROCESSED / "log2fc_targets.csv"
    train_curated_path = PROCESSED / "train_inhibition_curated.csv"
    cv_folds_path = FOLDS / "cv_folds.csv"

    log2fc = pd.read_csv(log2fc_path)
    train_curated = pd.read_csv(train_curated_path)
    cv_folds = pd.read_csv(cv_folds_path)
    logger.info(f"loaded {log2fc_path.name}: {log2fc.shape}")
    logger.info(f"loaded {train_curated_path.name}: {train_curated.shape}")
    logger.info(f"loaded {cv_folds_path.name}: {cv_folds.shape}")

    merged = log2fc.merge(
        train_curated[["inchikey", "canonical_smiles"]], on="inchikey", how="left"
    )
    n_missing_smiles = int(merged["canonical_smiles"].isna().sum())
    logger.info(f"rows after SMILES join: {len(merged)} (missing SMILES: {n_missing_smiles})")
    if n_missing_smiles:
        raise ValueError(
            f"{n_missing_smiles} log2fc compounds have no canonical_smiles match in "
            "train_inhibition_curated.csv by InChIKey -- train_log2fc_encoder.py found "
            "0 such compounds, so this is a real change worth stopping for."
        )
    if len(merged) != len(log2fc):
        raise ValueError(
            f"row count changed on SMILES join ({len(log2fc)} -> {len(merged)}) -- "
            "the InChIKey join should be 1:many at most on the train_curated side, "
            "never row-multiplying; stopping."
        )

    target_cols = [f"{iso}_log2fc" for iso in ISOFORMS]
    n_nan = int(merged[target_cols].isna().sum().sum())
    logger.info(f"NaN values across the four log2fc target columns: {n_nan}")
    if n_nan:
        raise ValueError(
            "unexpected NaN in log2fc targets -- train_log2fc_encoder.py found zero "
            "nulls and a full 4-isoform panel per compound, so this should never "
            "trigger; stopping."
        )

    split_df = assign_screen_split(
        cv_folds, repeat_col=REPEAT_COL, test_fold=TEST_FOLD, val_fraction=VAL_FRACTION, seed=SEED
    )
    n_before_screen_filter = len(merged)
    merged = merged.merge(split_df[["inchikey", "screen_split"]], on="inchikey", how="left")
    n_missing_split = int(merged["screen_split"].isna().sum())
    logger.info(
        f"rows after screen_split join: {len(merged)} (missing screen_split: {n_missing_split})"
    )
    if n_missing_split or len(merged) != n_before_screen_filter:
        raise ValueError(
            f"screen_split join broken (missing={n_missing_split}, row count "
            f"{n_before_screen_filter} -> {len(merged)}) -- every log2fc compound is "
            "expected to be one of cv_folds.csv's 4,905 training compounds; stopping."
        )

    logger.info(f"screen_split counts before excluding screen_test: {merged['screen_split'].value_counts().to_dict()}")
    blind = merged[merged["screen_split"] != "screen_test"].copy()
    logger.info(
        f"rows before excluding screen_test: {n_before_screen_filter}, after: {len(blind)} "
        f"({n_before_screen_filter - len(blind)} screen_test compounds dropped from pretraining)"
    )

    return blind, target_cols


def build_compound_index_with_smiles(logger: logging.Logger) -> pd.DataFrame:
    """Identical to train_log2fc_encoder.py's own function -- embeddings are still
    extracted for the full 5,655-compound index (`screen_test` compounds included).
    Only pretraining *supervision* excludes `screen_test` (see `load_training_data`);
    the frozen feature-extraction step must still cover every compound 04a needs
    predictions for, exactly mirroring how 04a's own RF/XGBoost/LightGBM models are fit
    train-only but predict on the full held-out fold.
    """
    index_cols = ["Molecule_Name", "inchikey", "split"]
    compound_index = pd.read_csv(PROCESSED / "tabular_baseline_features.csv", usecols=index_cols)
    logger.info(f"loaded tabular_baseline_features.csv compound index: {compound_index.shape}")

    train_curated = pd.read_csv(PROCESSED / "train_inhibition_curated.csv")
    test_curated = pd.read_csv(PROCESSED / "test_blinded_curated.csv")
    smiles_lookup = pd.concat(
        [
            train_curated[["inchikey", "canonical_smiles"]],
            test_curated[["inchikey", "canonical_smiles"]],
        ],
        ignore_index=True,
    )

    n_before = len(compound_index)
    compound_index = compound_index.merge(smiles_lookup, on="inchikey", how="left")
    n_missing = int(compound_index["canonical_smiles"].isna().sum())
    logger.info(
        f"rows after SMILES join: {len(compound_index)} (before: {n_before}, "
        f"missing SMILES: {n_missing})"
    )
    if n_missing or len(compound_index) != n_before:
        raise ValueError(
            f"compound index alignment broken on SMILES join (missing={n_missing}, "
            f"row count {n_before} -> {len(compound_index)}) -- stopping."
        )
    return compound_index


class EpochLogger(Callback):
    """Identical to train_log2fc_encoder.py's own callback -- one log line per epoch."""

    def __init__(self, logger: logging.Logger, total_epochs: int):
        self.logger = logger
        self.total_epochs = total_epochs

    def on_train_epoch_end(self, trainer, pl_module):
        loss = trainer.callback_metrics.get("train_loss")
        loss_val = float(loss) if loss is not None else float("nan")
        self.logger.info(
            f"epoch {trainer.current_epoch + 1}/{self.total_epochs} done, train_loss={loss_val:.4f}"
        )


def main():
    logger = setup_logging()
    script_start = time.time()
    logger.info("=" * 70)
    logger.info("starting train_log2fc_encoder_blind.py")
    logger.info(
        "LEAKAGE CHECK FOR 04a: retrains the log2fc encoder with screen_test compounds "
        "excluded from pretraining -- 03b's own frozen artifacts are untouched."
    )

    import chemprop
    import lightning

    logger.info(f"python: {sys.version.split()[0]}")
    logger.info(f"torch: {torch.__version__}")
    logger.info(f"chemprop: {chemprop.__version__}")
    logger.info(f"lightning: {lightning.__version__}")
    logger.info(f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS')}")
    logger.info(f"seed: {SEED}")
    logger.info(
        f"hyperparameters: epochs={N_EPOCHS}, batch_size={BATCH_SIZE}, "
        f"d_h={D_H}, depth={DEPTH}, aggregation=mean, batch_norm=True "
        "(identical to train_log2fc_encoder.py -- only the training population differs)"
    )

    from chemprop import data as cp_data
    from chemprop import models as cp_models
    from chemprop import nn as cp_nn
    from lightning import pytorch as pl

    torch.manual_seed(SEED)
    pl.seed_everything(SEED, workers=True, verbose=False)

    # ---- 1. load training data (screen_test compounds excluded) -----------------
    log2fc, target_cols = load_training_data(logger)
    smiles_list = log2fc["canonical_smiles"].tolist()
    y = log2fc[target_cols].to_numpy(dtype=np.float64)

    scaler = StandardScaler()
    y_scaled = scaler.fit_transform(y)
    logger.info(
        "target scaler fit -- per-task mean/std used for training (predictions are "
        "un-scaled back to raw log2fc units by the model's own output_transform):"
    )
    for iso, mean, std in zip(ISOFORMS, scaler.mean_, scaler.scale_):
        logger.info(f"  {iso}_log2fc: mean={mean:.4f}, std={std:.4f}")

    datapoints = [
        cp_data.MoleculeDatapoint.from_smi(smi, y=y_scaled[i]) for i, smi in enumerate(smiles_list)
    ]
    dataset = cp_data.MoleculeDataset(datapoints)
    train_loader = cp_data.build_dataloader(
        dataset, batch_size=BATCH_SIZE, shuffle=True, num_workers=0, drop_last=False
    )
    logger.info(f"training dataset: {len(dataset)} compounds, {len(target_cols)} tasks")

    # ---- 2. build model, random init, no CheMeleon weights anywhere here --------
    mp = cp_nn.BondMessagePassing(d_h=D_H, depth=DEPTH)
    agg = cp_nn.MeanAggregation()
    output_transform = cp_nn.UnscaleTransform.from_standard_scaler(scaler)
    ffn = cp_nn.RegressionFFN(
        n_tasks=len(target_cols), input_dim=mp.output_dim, output_transform=output_transform
    )
    mpnn = cp_models.MPNN(mp, agg, ffn, batch_norm=True, metrics=[cp_nn.metrics.RMSE()])

    trainer = pl.Trainer(
        max_epochs=N_EPOCHS,
        accelerator="cpu",
        devices=1,
        enable_progress_bar=False,
        enable_model_summary=False,
        enable_checkpointing=False,
        logger=False,
        callbacks=[EpochLogger(logger, N_EPOCHS)],
    )

    # ---- 3. train ----------------------------------------------------------------
    train_start = time.time()
    logger.info(f"training start: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(train_start))}")
    trainer.fit(mpnn, train_loader)
    train_end = time.time()
    train_elapsed = train_end - train_start
    logger.info(f"training end:   {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(train_end))}")
    logger.info(f"training elapsed: {train_elapsed:.1f}s ({train_elapsed / 60:.1f} min)")

    # ---- 4. checkpoint the trained encoder ---------------------------------------
    MODELS.mkdir(parents=True, exist_ok=True)
    trained_mp = mpnn.message_passing
    trained_agg = mpnn.agg

    hparams = dict(trained_mp.hparams)
    hparams.pop("cls", None)
    if hasattr(hparams.get("activation"), "value"):
        hparams["activation"] = hparams["activation"].value
    checkpoint = {"hyper_parameters": hparams, "state_dict": trained_mp.state_dict()}
    encoder_path = MODELS / "log2fc_encoder_blind_mp.pt"
    torch.save(checkpoint, encoder_path)
    logger.info(f"wrote encoder checkpoint: {encoder_path}")

    # ---- 5. freeze encoder, extract embeddings for the full 5,655-compound index -
    trained_mp.eval()
    trained_agg.eval()

    compound_index = build_compound_index_with_smiles(logger)
    embed_smiles = compound_index["canonical_smiles"].tolist()

    embed_start = time.time()
    logger.info(
        f"embedding extraction start: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(embed_start))} "
        f"({len(embed_smiles)} compounds)"
    )
    embeddings = frozen_encoder_embeddings(
        embed_smiles, trained_mp, trained_agg, batch_size=EMBED_BATCH_SIZE, device="cpu"
    )
    embed_end = time.time()
    embed_elapsed = embed_end - embed_start
    logger.info(f"embedding extraction end:   {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(embed_end))}")
    logger.info(f"embedding extraction elapsed: {embed_elapsed:.1f}s")
    logger.info(f"embeddings shape: {embeddings.shape}, dtype: {embeddings.dtype}")

    n_nan = int(np.isnan(embeddings).sum())
    logger.info(f"NaN values in embeddings: {n_nan}")
    if n_nan:
        raise ValueError("unexpected NaN in extracted embeddings -- stopping before save.")

    # ---- 6. save + reread verification -------------------------------------------
    embeddings_path = PROCESSED / "log2fc_pretrained_blind_embeddings.npy"
    index_path = PROCESSED / "log2fc_pretrained_blind_embeddings_index.csv"

    logger.info(f"rows before save: {len(embeddings)}")
    np.save(embeddings_path, embeddings)
    compound_index[["Molecule_Name", "inchikey", "split"]].to_csv(index_path, index=False)
    logger.info(f"wrote {embeddings_path}")
    logger.info(f"wrote {index_path}")

    reread_embeddings = np.load(embeddings_path)
    reread_index = pd.read_csv(index_path)
    logger.info(f"rows after reading back from disk: {reread_embeddings.shape[0]} / {len(reread_index)}")
    if reread_embeddings.shape != embeddings.shape or len(reread_index) != len(compound_index):
        raise ValueError("save/reread row count or shape mismatch -- stopping.")
    if not (reread_index["inchikey"].values == compound_index["inchikey"].values).all():
        raise ValueError("embedding index row order does not match compound_index after reread -- stopping.")
    logger.info(
        "alignment confirmed: log2fc_pretrained_blind_embeddings.npy row i corresponds to "
        "log2fc_pretrained_blind_embeddings_index.csv row i."
    )

    total_elapsed = time.time() - script_start
    logger.info(f"total script wall time: {total_elapsed:.1f}s ({total_elapsed / 60:.1f} min)")
    logger.info("done.")


if __name__ == "__main__":
    main()
