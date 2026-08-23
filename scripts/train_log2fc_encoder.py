"""Train a fresh, random-init Chemprop encoder on the log2fc primary-screen data,
freeze it, and extract embeddings for every compound in the tabular feature index.

STANDALONE SCRIPT, NOT A NOTEBOOK CELL: this is long-running compute (training, not a
single forward pass), so per the project's agreed convention, it belongs in a script
the person launches themselves -- e.g.

    cd /path/to/OpenADMET-CYP-Blind-Challenge
    caffeinate -i nohup python scripts/train_log2fc_encoder.py > /dev/null 2>&1 &

-- not inside an agentic session or a notebook cell. `notebooks/03b_log2fc_pretrained_
encoder.ipynb` is a thin companion notebook that loads and reports on this script's
output artifacts; it does not run any training itself.

RECIPE: Wichrowski et al. (Chem. Res. Toxicol. 2026) -- pretrain a Chemprop encoder on a
related, larger-coverage task, then freeze it and use it purely as a fixed feature
extractor for the real, smaller target task. This is NOT the Dinh Pham-style variant
(fine-tuning the pretrained heads further), and it is NOT a CheMeleon continued-
pretraining variant -- the encoder here is trained completely from scratch (random
initialization), with no CheMeleon weights involved anywhere in this script.

WHY NOT CHEMELEON CONTINUED PRETRAINING -- stated plainly, not just asserted: OpenADMET's
own team (Burns & West, "Throwing Everything AND the Kitchen Sink at CheMeleon", Aug
2026) systematically tested continued pretraining of CheMeleon's own encoder across many
datasets/featurizations and found that pretraining-objective improvements do not
reliably transfer to downstream fine-tuning performance, with several conditions
actively hurting relative to the original checkpoint. That result is specifically about
modifying CheMeleon's own foundational training run -- it does not bear on training a
separate, freshly-initialized Chemprop encoder on a real, task-relevant endpoint
(log2fc), which is a different mechanism entirely. The two should not be confused: this
script does the latter, not the former.

LEAKAGE-AVOIDANCE DECISION (already agreed, not re-litigated here): notebook 03 (Section
5) found that the single-concentration/log2fc population (4,376 compounds) overlaps
almost completely with the pIC50 training set (0 log2fc compounds have no pIC50 label
at all, by InChIKey). Per the earlier agreement, this script trains on the FULL
4,376-compound log2fc set as one single global run, accepting that overlap as a
documented, CheMeleon-precedented assumption (large pretrained encoders routinely train
on data that overlaps downstream tasks -- CheMeleon itself is an instance of exactly
this), rather than the much more expensive per-fold-retraining alternative. That
alternative hasn't been shown to help anything yet, so it doesn't get the expensive
treatment until 04 shows it's needed.

ARCHITECTURE SIZE -- benchmarked, then decided by the user, not picked silently: two
options were benchmarked directly on the real 4,376-compound log2fc data before this
script was written: chemprop's own default (d_h=300, depth=3) measured at ~13s/epoch
(~11 min for 50 epochs), vs. an architecture matching CheMeleon's own size (d_h=2048,
depth=6) measured at ~758s/epoch (~10.5 HOURS for 50 epochs) -- a ~58x difference, all
CPU/single-threaded per the OMP constraint below. The user chose chemprop's own default
(smaller, well-tested, appropriately sized for a 4,376-compound run with no
validation-based early stopping, and far cheaper on CPU) over CheMeleon-matching
embedding dimensionality. Embeddings from this script are therefore 300-dim, not
2048-dim like `chemeleon_embeddings.npy` -- a different shape is expected and fine; nothing
downstream requires dimensional parity between feature options.

EPOCHS: 50, matching chemprop's own CLI default (`--epochs 50`) -- not an invented
number. There's no validation split here (see leakage decision above: full dataset, one
global run), so there's no validation-loss-based early-stopping signal to tune epochs
against; a fixed, tool-default budget is used instead. At the benchmarked ~13s/epoch
this is a ~11 minute run, cheap enough that this default isn't a large commitment either
way.

THREADING (reused directly from notebook 03's `chemeleon_embeddings` section, not
rediscovered here): rdkit and torch each bundle their own OpenMP runtime, and running
either training or inference multi-threaded in the same process segfaults on this
machine. `OMP_NUM_THREADS=1` is set below, before rdkit/torch/chemprop are first
imported anywhere in this process, and training/inference both run CPU-only --
Apple Silicon MPS was already ruled out for this exact chemprop message-passing
aggregation op (`scatter_reduce` not implemented for MPS in this torch version).

SEEDS: `SEED = 42` (matching the rest of this project's convention) controls PyTorch's
global RNG (weight initialization) and the training dataloader's shuffling; logged
below before use. This is a single global training run, not a per-CV-fold loop, so
`CLAUDE.md`'s "seeds must vary per fold" rule (about per-fold model-fitting seeds) does
not apply here -- there is only one fold.

CHEMPROP VERSION: this script was written and its API calls verified against whichever
chemprop version is installed in the `cyp-admet` conda environment at the time it was
run -- logged below at runtime rather than hardcoded, since chemprop's Python API has
changed across versions.

OUTPUTS:
    logs/03b_train_log2fc_encoder.log      -- progress log (tail -f while this runs)
    models/log2fc_encoder_mp.pt            -- frozen encoder checkpoint
                                               (chemprop BondMessagePassing state_dict +
                                               hyper_parameters, same on-disk format
                                               convention as CheMeleon's own checkpoint)
    data/processed/log2fc_pretrained_embeddings.npy       -- (5655, 300) float32
    data/processed/log2fc_pretrained_embeddings_index.csv -- matching compound index
                                                              (Molecule_Name, inchikey,
                                                              split), same row order

This is one candidate feature option among several for notebook 04's tabular baselines
(alongside ECFP4+descriptors, RDKit 2D descriptors, and raw CheMeleon embeddings) -- no
claim about its usefulness is made anywhere in this script. That's 04's job to test.
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

from src.features import frozen_encoder_embeddings

RAW = REPO_ROOT / "data" / "raw"
PROCESSED = REPO_ROOT / "data" / "processed"
MODELS = REPO_ROOT / "models"
LOGS = REPO_ROOT / "logs"

ISOFORMS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]

SEED = 42
N_EPOCHS = 50  # chemprop's own CLI default -- see module docstring
BATCH_SIZE = 64  # chemprop's own default (`build_dataloader`'s default)
D_H = 300  # chemprop's own default message-hidden-dim -- user-chosen over CheMeleon's 2048, see docstring
DEPTH = 3  # chemprop's own default message-passing depth -- user-chosen over CheMeleon's 6, see docstring
EMBED_BATCH_SIZE = 128  # matches chemeleon_embeddings' default

LOG_PATH = LOGS / "03b_train_log2fc_encoder.log"


def setup_logging() -> logging.Logger:
    LOGS.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("train_log2fc_encoder")
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
    """Load log2fc_targets.csv (notebook 03 Section 7) and attach canonical_smiles from
    the curated training data, per notebook 01's InChIKey identity convention. No
    canonicalization/InChIKey computation happens here -- both files already carry
    pre-computed, pre-verified identity columns; this only joins them.
    """
    log2fc_path = PROCESSED / "log2fc_targets.csv"
    train_curated_path = PROCESSED / "train_inhibition_curated.csv"

    log2fc = pd.read_csv(log2fc_path)
    train_curated = pd.read_csv(train_curated_path)
    logger.info(f"loaded {log2fc_path.name}: {log2fc.shape}")
    logger.info(f"loaded {train_curated_path.name}: {train_curated.shape}")

    merged = log2fc.merge(
        train_curated[["inchikey", "canonical_smiles"]], on="inchikey", how="left"
    )
    n_missing_smiles = int(merged["canonical_smiles"].isna().sum())
    logger.info(f"rows after SMILES join: {len(merged)} (missing SMILES: {n_missing_smiles})")
    if n_missing_smiles:
        raise ValueError(
            f"{n_missing_smiles} log2fc compounds have no canonical_smiles match in "
            "train_inhibition_curated.csv by InChIKey -- notebook 03 Section 5c found "
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
            "unexpected NaN in log2fc targets -- notebook 03 Section 5a/7 found zero "
            "nulls and a full 4-isoform panel per compound, so this should never "
            "trigger; stopping."
        )

    return merged, target_cols


def build_compound_index_with_smiles(logger: logging.Logger) -> pd.DataFrame:
    """Load tabular_baseline_features.csv's compound index (Molecule_Name, inchikey,
    split) -- the authoritative reference for which 5,655 compounds and what row order
    every feature file in this project uses -- and attach canonical_smiles from the
    curated train/test files (not carried in tabular_baseline_features.csv itself, so
    it has to be rejoined here rather than assumed present).
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
    """Lightning callback that logs one line per completed epoch to the shared logger
    -- "epoch progress if easily available" per the task brief; available directly via
    Lightning's own `on_train_epoch_end` hook, so no separate progress-tracking
    machinery is added."""

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
    logger.info("starting train_log2fc_encoder.py")

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
        f"d_h={D_H}, depth={DEPTH}, aggregation=mean"
    )

    from chemprop import data as cp_data
    from chemprop import models as cp_models
    from chemprop import nn as cp_nn
    from lightning import pytorch as pl

    torch.manual_seed(SEED)
    pl.seed_everything(SEED, workers=True, verbose=False)

    # ---- 1. load training data --------------------------------------------------
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
    encoder_path = MODELS / "log2fc_encoder_mp.pt"
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
    embeddings_path = PROCESSED / "log2fc_pretrained_embeddings.npy"
    index_path = PROCESSED / "log2fc_pretrained_embeddings_index.csv"

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
    logger.info("alignment confirmed: log2fc_pretrained_embeddings.npy row i corresponds to "
                "log2fc_pretrained_embeddings_index.csv row i.")

    total_elapsed = time.time() - script_start
    logger.info(f"total script wall time: {total_elapsed:.1f}s ({total_elapsed / 60:.1f} min)")
    logger.info("done.")


if __name__ == "__main__":
    main()
