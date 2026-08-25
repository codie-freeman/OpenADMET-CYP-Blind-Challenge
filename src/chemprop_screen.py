"""Shared Chemprop CLI orchestration for notebook 04a's baseline screen -- used by both
`scripts/run_baseline_screen_chemprop.py` (random-init) and
`scripts/run_baseline_screen_chemprop_chemeleon.py` (CheMeleon-init), which differ only
in the message-passing encoder's initialization/architecture. Factored out here, rather
than duplicated across the two scripts, so the per-epoch progress-logging mechanism,
the subprocess stdin fix, and the training-population construction only ever need
fixing in one place if either script needs changing again.

Not placed in `src/features.py`: that module is scoped to fingerprints/similarity/
splits (per its own docstring, and CLAUDE.md's "fingerprints/similarity/splits go
through one shared module" rule) -- `assign_screen_split` (a split) lives there and is
called from here (in `load_screen_population`), but Chemprop CLI subprocess
orchestration is a different kind of shared concern and gets its own module.

See each function's own docstring for why it exists; this module-level docstring
doesn't repeat that reasoning. `run_chemprop_train` is the one function whose CLI
arguments genuinely differ between the two calling scripts (architecture flags) --
callers pass their own `architecture_args`; everything else here is identical between
random-init and CheMeleon-init.
"""

import logging
import subprocess
import sys
import threading
from pathlib import Path

import pandas as pd

from src.features import assign_screen_split
from src.vendor.openadmet_eval.config import REGRESSION_ENDPOINTS

ISOFORM_SHORT_NAMES = [endpoint.split("_")[0] for endpoint in REGRESSION_ENDPOINTS]


def setup_logging(log_path: Path, logger_name: str) -> logging.Logger:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    file_handler = logging.FileHandler(log_path, mode="a")
    file_handler.setFormatter(fmt)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(fmt)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def run_subprocess_streamed(argv: list, logger: logging.Logger) -> None:
    """Run `argv`, streaming its combined stdout/stderr into `logger` line-by-line as it
    runs (rather than capturing and logging only at the end) -- so `tail -f` on the
    caller's log file shows live chemprop progress, matching this project's established
    "tail -f while this runs" convention.

    `stdin=subprocess.DEVNULL` is explicit -- confirmed necessary on this project's own
    first CheMeleon-init attempt, which crashed a `chemprop predict` child process with
    `Fatal Python error: init_sys_streams ... OSError: [Errno 9] Bad file descriptor`,
    most consistent with an inherited, invalid stdin under this script's own
    `nohup ... &` launch convention. This is the one place that fix lives -- both
    calling scripts get it automatically.
    """
    logger.info(f"running: {' '.join(str(a) for a in argv)}")
    process = subprocess.Popen(
        argv, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, bufsize=1,
    )
    for line in process.stdout:
        logger.info(line.rstrip())
    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(
            f"command failed with exit code {return_code}: {' '.join(str(a) for a in argv)}"
        )


def _tail_metrics_csv(
    metrics_csv_path: Path, logger: logging.Logger, stop_event: threading.Event, poll_interval: float = 10.0
) -> None:
    """Poll `metrics_csv_path` (Lightning's own `CSVLogger` output, written
    incrementally during training) every `poll_interval` seconds and log one clean
    "epoch N done" line per newly-completed epoch. Chemprop's CLI renders training
    progress via a `rich` progress bar that, once its stdout is a pipe rather than a
    real terminal (exactly our situation -- `run_subprocess_streamed` captures it),
    only prints its *final* rendered state once, not one line per epoch -- confirmed
    directly on this project's own first CheMeleon-init attempt, whose log had a single
    progress line for the whole run and was otherwise silent for ~24 minutes. This is
    the smallest fix that recovers real progress without abandoning the CLI-subprocess
    pattern (rewriting as a direct Python-API script instead would throw away
    `--splits-column`, `--patience`, and every other CLI convenience used here).

    Runs in a background thread alongside `run_subprocess_streamed`'s blocking stdout
    read loop; stops when `stop_event` is set. Train-loss and val-loss land on separate
    CSV rows for the same epoch (Lightning's own logging convention) -- an epoch is
    only reported once both are present, so a row is never reported with a missing
    value.
    """
    seen_epochs = set()
    while not stop_event.is_set():
        if metrics_csv_path.exists():
            try:
                metrics = pd.read_csv(metrics_csv_path)
            except (pd.errors.EmptyDataError, pd.errors.ParserError, OSError):
                metrics = None  # file mid-write or not yet flushed -- try again next poll
            if metrics is not None and len(metrics) and "epoch" in metrics.columns:
                for epoch, group in metrics.groupby("epoch"):
                    if epoch in seen_epochs:
                        continue
                    train_loss = group["train_loss_epoch"].dropna() if "train_loss_epoch" in group else pd.Series(dtype=float)
                    val_loss = group["val_loss"].dropna() if "val_loss" in group else pd.Series(dtype=float)
                    if train_loss.empty or val_loss.empty:
                        continue  # this epoch's train/val rows haven't both landed yet
                    seen_epochs.add(epoch)
                    logger.info(
                        f"  epoch {int(epoch)} done: train_loss={train_loss.iloc[-1]:.4f}, "
                        f"val_loss={val_loss.iloc[-1]:.4f}"
                    )
        stop_event.wait(poll_interval)


def load_screen_population(
    cv_folds_path: Path,
    train_curated_path: Path,
    repeat_col: str,
    test_fold: int,
    val_fraction: float,
    seed: int,
    logger: logging.Logger,
) -> pd.DataFrame:
    """Load `cv_folds.csv` + `train_inhibition_curated.csv`, assign this screen's split
    via `assign_screen_split`, and log overall + per-isoform compound counts per split
    bucket (CLAUDE.md: row counts logged before/after every filtering step). Calling
    this with the same arguments from both the random-init and CheMeleon-init scripts
    (and from `notebooks/04a_baseline_screen.ipynb`'s own tabular configs) yields a
    byte-identical split -- `assign_screen_split` is a pure, deterministic function.
    """
    cv_folds = pd.read_csv(cv_folds_path)
    curated = pd.read_csv(train_curated_path)
    logger.info(f"loaded {cv_folds_path.name}: {cv_folds.shape}")
    logger.info(f"loaded {train_curated_path.name}: {curated.shape}")

    split_df = assign_screen_split(
        cv_folds, repeat_col=repeat_col, test_fold=test_fold, val_fraction=val_fraction, seed=seed
    )
    logger.info(
        f"assign_screen_split(repeat_col={repeat_col!r}, test_fold={test_fold}, "
        f"val_fraction={val_fraction}, seed={seed}) split counts: "
        f"{split_df['screen_split'].value_counts().to_dict()}"
    )

    merged = split_df.merge(
        curated[["inchikey", "canonical_smiles", *REGRESSION_ENDPOINTS]], on="inchikey", how="left"
    )
    if len(merged) != len(split_df) or merged["canonical_smiles"].isna().any():
        raise ValueError(
            "screen split compounds did not align 1:1 against train_inhibition_curated.csv "
            "by inchikey -- stopping before training anything."
        )
    logger.info(f"rows after SMILES/target join: {len(merged)} (missing: 0)")

    for endpoint in REGRESSION_ENDPOINTS:
        counts = merged.dropna(subset=[endpoint]).groupby("screen_split").size().to_dict()
        logger.info(f"  {endpoint}: labeled-compound counts per split = {counts}")

    return merged


def build_predict_csv(population: pd.DataFrame, chemprop_runs_dir: Path, logger: logging.Logger) -> Path:
    """Write the one shared prediction-input CSV covering every `screen_test` compound
    (Molecule_Name, inchikey, canonical_smiles) -- reused, unchanged, across every
    model in either script, regardless of which isoform(s) that model was trained on.
    """
    chemprop_runs_dir.mkdir(parents=True, exist_ok=True)
    test_df = population.loc[
        population["screen_split"] == "screen_test", ["Molecule_Name", "inchikey", "canonical_smiles"]
    ].copy()
    out_path = chemprop_runs_dir / "predict_input.csv"
    test_df.to_csv(out_path, index=False)
    logger.info(f"wrote {out_path} ({len(test_df)} screen_test compounds)")
    return out_path


def build_training_csv(
    population: pd.DataFrame,
    target_cols: list,
    out_path: Path,
    logger: logging.Logger,
    require_all_targets: bool,
) -> pd.DataFrame:
    """Build one model's training CSV: pooled `screen_inner_train` + `screen_inner_val`
    compounds only (`screen_test` is never included here), mapped to chemprop's
    expected `train`/`val` split values via a `chemprop_split` column.

    `require_all_targets=True` (single-task models) additionally drops rows lacking
    that one isoform's label -- logged before/after, per CLAUDE.md. `require_all_
    targets=False` (multitask) keeps every pooled-train compound; chemprop masks any
    missing per-task target automatically.
    """
    pooled = population[population["screen_split"].isin(["screen_inner_train", "screen_inner_val"])].copy()
    n_before = len(pooled)
    if require_all_targets:
        pooled = pooled.dropna(subset=target_cols)
    logger.info(
        f"{out_path.name}: pooled-train rows before label filter={n_before}, "
        f"after={len(pooled)} (require_all_targets={require_all_targets}, targets={target_cols})"
    )

    split_map = {"screen_inner_train": "train", "screen_inner_val": "val"}
    pooled["chemprop_split"] = pooled["screen_split"].map(split_map)
    if pooled["chemprop_split"].isna().any():
        raise ValueError(f"unexpected screen_split value(s) building {out_path.name} -- stopping.")

    cols = ["canonical_smiles", *target_cols, "chemprop_split"]
    pooled[cols].to_csv(out_path, index=False)
    logger.info(
        f"wrote {out_path} ({len(pooled)} rows, chemprop_split counts: "
        f"{pooled['chemprop_split'].value_counts().to_dict()})"
    )
    return pooled


def run_chemprop_train(
    train_csv: Path,
    target_cols: list,
    output_dir: Path,
    logger: logging.Logger,
    architecture_args: list,
    epochs: int,
    patience: int,
    seed: int,
) -> Path:
    """Run `chemprop train`. `architecture_args` is the one place the two calling
    scripts differ (random-init's `--message-hidden-dim`/`--depth`/`--aggregation`/
    `--batch-norm` vs. CheMeleon-init's `--from-foundation CHEMELEON`/
    `--multi-hot-atom-featurizer-mode V2`) -- everything else about how a Chemprop
    training run is launched, logged, and progress-tracked is identical and lives here.
    """
    argv = [
        "chemprop", "train",
        "-i", str(train_csv),
        "-s", "canonical_smiles",
        "--target-columns", *target_cols,
        "--splits-column", "chemprop_split",
        "-t", "regression",
        *architecture_args,
        "--epochs", str(epochs),
        "--patience", str(patience),
        "--data-seed", str(seed),
        "--pytorch-seed", str(seed),
        "--accelerator", "cpu",
        "--devices", "1",
        "-o", str(output_dir),
    ]
    metrics_csv_path = output_dir / "model_0" / "trainer_logs" / "version_0" / "metrics.csv"
    stop_event = threading.Event()
    tail_thread = threading.Thread(
        target=_tail_metrics_csv, args=(metrics_csv_path, logger, stop_event), daemon=True
    )
    tail_thread.start()
    try:
        run_subprocess_streamed(argv, logger)
    finally:
        stop_event.set()
        tail_thread.join(timeout=15.0)
    return output_dir


def run_chemprop_predict(model_dir: Path, predict_csv: Path, output_csv: Path, logger: logging.Logger) -> None:
    argv = [
        "chemprop", "predict",
        "-i", str(predict_csv),
        "-s", "canonical_smiles",
        "--model-paths", str(model_dir),
        "--accelerator", "cpu",
        "--devices", "1",
        "-o", str(output_csv),
    ]
    run_subprocess_streamed(argv, logger)


def verify_predictions(
    csv_path: Path, expected_molecule_names: set, target_cols: list, logger: logging.Logger
) -> None:
    """Reread `csv_path` from disk and confirm: every expected `screen_test` compound is
    present exactly once, and every target column the model covers is fully populated
    (a forward pass should never itself produce a missing value for a curated,
    already-parseable SMILES).
    """
    reread = pd.read_csv(csv_path)
    logger.info(f"rows after reading {csv_path.name} back from disk: {len(reread)}")
    if set(reread["Molecule_Name"]) != expected_molecule_names:
        raise ValueError(f"{csv_path.name}: Molecule_Name set does not match screen_test compounds -- stopping.")
    if len(reread) != len(expected_molecule_names):
        raise ValueError(f"{csv_path.name}: row count {len(reread)} != {len(expected_molecule_names)} -- stopping.")
    n_nan = int(reread[target_cols].isna().sum().sum())
    logger.info(f"NaN values in {csv_path.name} target column(s): {n_nan}")
    if n_nan:
        raise ValueError(f"{csv_path.name}: unexpected NaN in predicted target column(s) -- stopping.")
