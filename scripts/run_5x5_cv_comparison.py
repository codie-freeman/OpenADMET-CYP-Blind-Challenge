"""Resumable driver for the 5x5 (25-fold) CV model comparison. Reads the static plan
written by `scripts/generate_5x5_cv_manifest.py` (`outputs/05_cv_comparison/
manifest.csv`, 300 rows = 12 configs x 5 repeats x 5 folds) and, for each row not
already completed, trains/predicts/scores it and writes results to disk immediately.

DEVIATION FROM PROJECT CONVENTION -- OMP_NUM_THREADS LEFT UNSET (2026-08-30): every
other script in this repo sets OMP_NUM_THREADS=1 before RDKit/PyTorch/Chemprop import
(CLAUDE.md convention, originally adopted out of general caution about arm64-native
behavior on the then-new `cyp-admet-v2` env). This script deliberately does NOT set it,
leaving it to whatever the launching shell provides. Evidence: an isolated timing check
run the same day, on this real repo's actual data and CheMeleon-init multitask config
(fold 0/repeat 0, `--epochs 2`, scratchpad-only, not part of this driver), measured a
clean, reproducible steady-state-epoch speedup with OMP_NUM_THREADS unset vs. =1 (56.4s
-> 27.2s, ~2.07x, 206% CPU, no crash/warning beyond a pre-existing unrelated
TensorBoardLogger-missing warning also present at OMP=1 and in the original 04a
baseline log) -- on top of an even larger ~13.2x speedup already present at OMP=1 from
the arm64-native migration itself (never previously isolated/measured). CAVEAT: this is
one fold, two epochs, one trial each -- not yet validated across repeats, full epoch
counts, or sustained multi-hour load. If the full run below shows instability
(crashes, degraded throughput, thermal/memory pressure) that this small check couldn't
reveal, revert to `os.environ.setdefault("OMP_NUM_THREADS", "1")` and re-launch --
resumability means nothing already completed is lost.

THIS STEP PRODUCES RAW PER-FOLD BOOTSTRAP SCORES ONLY -- no statistical comparison
(Levene's/ANOVA/Tukey/compact letter display); that's a separate step once results
exist.

RESUMABILITY: "done" for a given (config, repeat, fold) is determined purely by
checking whether that row's own prediction + score CSVs already exist and are
non-empty on disk (`is_done`) -- there is no shared, mutated manifest/status file. This
means multiple driver processes can run concurrently against the same manifest.csv
with zero shared mutable state and therefore no locking/race concerns -- the intended
parallelization mechanism (see `--family`/`--shard` below and the run-plan writeup).
If interrupted (Ctrl-C, disconnect, crash) mid-row, that row's outputs are never
written (predictions are only saved after a full, successful fit -- see
`run_tabular_or_naive`/`run_chemprop`), so it is correctly retried, in full, on the
next invocation. NOTE: this is fold-completion-granularity resumability, not
mid-epoch -- an interrupted Chemprop fold restarts that fold's training from epoch 0,
it does not resume from a checkpoint. The CheMeleon-init config's original ~100 min/
fold estimate (from the old `cyp-admet` Rosetta-environment 04a screen) is superseded
by the 2026-08-30 timing check referenced above (~3.5 min/fold extrapolated) -- see
that note for the caveats on trusting this before the full run confirms it.

A single row's failure (raised exception) is caught, logged at ERROR level with the
full traceback, and the driver moves on to the next row rather than aborting the whole
run -- appropriate for a many-hour unattended job where one flaky fold should not cost
the other 299. A failed row simply has no output files, so it naturally shows up as
"not done" and gets retried on the next invocation (e.g. after investigating the log).

PARALLELIZATION: `--family` selects which rows of the manifest this process is
responsible for (`naive`, `tabular`, `chemprop_randominit`, `chemprop_chemeleoninit`,
or `all` for every row regardless of family); `--shard i/n` further restricts to every
nth row (by manifest row order) within that selection, so e.g.
`--family chemprop_chemeleoninit --shard 0/2` and `--shard 1/2` run two independent
halves of the CheMeleon-init sweep concurrently, if ever needed again. No in-process
multiprocessing/pooling -- each invocation is single-purpose; concurrency (when used)
comes from running multiple OS processes, each its own `nohup`, each its own log file.
Per-row dispatch (tabular/naive fitting vs. Chemprop subprocess) reads each row's own
`family` column from the manifest, not `args.family` -- correct regardless of whether
this process is scoped to one family or `all`.

ENVIRONMENT: KMP_DUPLICATE_LIB_OK=TRUE is set before RDKit/PyTorch/Chemprop are first
imported, matching every existing script (CLAUDE.md convention). OMP_NUM_THREADS is
deliberately left unset -- see the dated note at the top of this docstring. This conda
env's own bin/ directory (holding the
`chemprop` console-script entry point) is also explicitly prepended to `PATH` before
any subprocess call -- confirmed necessary on this machine: pyenv shims intercept bare
`python`/`python3` ahead of `conda activate`'s own PATH entry in at least one shell
context tested here, and `run_chemprop_train`/`run_chemprop_predict`
(`src/chemprop_screen.py`) invoke `chemprop` as a bare command via `subprocess.Popen`,
which resolves purely from the *launching* shell's inherited PATH, not from
`sys.executable`. Deriving the prepended directory from `sys.executable` itself makes
this correct regardless of how the launching shell was set up, as long as this script
is itself invoked via this env's own python binary (see the run-plan writeup for the
exact launch command).
"""

import os

# OMP_NUM_THREADS intentionally left unset here -- see the dated deviation note at the
# top of this file's docstring. KMP_DUPLICATE_LIB_OK is unaffected, set as always.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import argparse
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

import lightgbm as lgb
from lightgbm import LGBMRegressor
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor

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

FOLDS_PATH = REPO_ROOT / "data" / "folds" / "cv_folds.csv"
CURATED_PATH = REPO_ROOT / "data" / "processed" / "train_inhibition_curated.csv"
PROCESSED = REPO_ROOT / "data" / "processed"
OUT = REPO_ROOT / "outputs" / "05_cv_comparison"
PRED_DIR = OUT / "predictions"
SCORE_DIR = OUT / "scores"
CHEMPROP_RUNS_DIR = OUT / "chemprop_runs"
LOG_DIR = REPO_ROOT / "logs"
MANIFEST_PATH = OUT / "manifest.csv"

VAL_FRACTION = 0.15  # matches assign_screen_split's/04a's own precedent
XGB_LGBM_EARLY_STOPPING_ROUNDS = 10  # matches 04a precedent

FEATURE_FILES = {
    "ecfp4_narrow": "tabular_baseline_features.csv",
    "chemeleon": "chemeleon_embeddings.npy",
    "mordred_pca": "tabular_mordred_pca.csv",
}

CHEMPROP_ARCHITECTURE_ARGS = {
    "chemprop_randominit": [
        "--message-hidden-dim", "300", "--depth", "3", "--aggregation", "mean", "--batch-norm",
    ],
    "chemprop_chemeleoninit": ["--from-foundation", "CHEMELEON", "--multi-hot-atom-featurizer-mode", "V2"],
}
CHEMPROP_EPOCHS = 50
CHEMPROP_PATIENCE = 5


def result_paths(config: str, repeat: int, fold: int) -> tuple[Path, Path]:
    tag = f"{config}__repeat{repeat}_fold{fold}"
    return PRED_DIR / f"{tag}.csv", SCORE_DIR / f"{tag}.csv"


def is_done(config: str, repeat: int, fold: int) -> bool:
    pred_path, score_path = result_paths(config, repeat, fold)
    if not (pred_path.exists() and score_path.exists()):
        return False
    try:
        return len(pd.read_csv(pred_path)) > 0 and len(pd.read_csv(score_path)) > 0
    except (pd.errors.EmptyDataError, OSError):
        return False


def load_feature_matrix(name: str) -> tuple[pd.DataFrame, np.ndarray]:
    """Load one of the three CV feature sources, restricted to `split == "train"` (the
    4,905 `cv_folds.csv` compounds) and re-indexed 0..N-1 in that restricted order.

    `tabular_baseline_features.csv` and `chemeleon_embeddings.npy` also carry the 750
    real blinded-test-set rows (used elsewhere for the final submission) --
    `tabular_mordred_pca.csv` does not (built train-only). Confirmed directly (not
    assumed from the 04a screen's own "row order is identical across all feature
    sources" note, which predates this Mordred+PCA file and does not apply to it
    as-is): once each source is filtered to `split == "train"`, all three -- plus
    `cv_folds.csv` itself -- share byte-identical (Molecule_Name, inchikey) row order.
    That shared order is what lets a single position index, built once per fold from
    `cv_folds.csv`, index correctly into any of the three.
    """
    path = PROCESSED / FEATURE_FILES[name]
    if path.suffix == ".csv":
        df = pd.read_csv(path)
        if "split" in df.columns:
            df = df[df["split"] == "train"].reset_index(drop=True)
        index_df = df[["Molecule_Name", "inchikey"]].reset_index(drop=True)
        drop_cols = [c for c in ["Molecule_Name", "inchikey", "split"] if c in df.columns]
        X = df.drop(columns=drop_cols).to_numpy(dtype=np.float64)
    else:
        X_full = np.load(path)
        index_full = pd.read_csv(path.with_name(f"{path.stem}_index.csv"))
        train_mask = (index_full["split"] == "train").to_numpy()
        index_df = index_full.loc[train_mask, ["Molecule_Name", "inchikey"]].reset_index(drop=True)
        X = X_full[train_mask]
    return index_df, X


def run_tabular_or_naive(
    config: str, repeat: int, fold: int, seed: int, logger
) -> pd.DataFrame:
    repeat_col = f"repeat_{repeat}"
    population = load_screen_population(
        FOLDS_PATH, CURATED_PATH, repeat_col, fold, VAL_FRACTION, seed, logger
    )
    test_ids = population.loc[
        population["screen_split"] == "screen_test", ["Molecule_Name", "inchikey"]
    ].reset_index(drop=True)

    if config == "naive_mean":
        pred_df = test_ids.copy()
        for endpoint in REGRESSION_ENDPOINTS:
            pool = population[
                population["screen_split"].isin(["screen_inner_train", "screen_inner_val"])
                & population[endpoint].notna()
            ]
            model = DummyRegressor(strategy="mean")
            model.fit(np.zeros((len(pool), 1)), pool[endpoint].to_numpy())
            pred_df[endpoint] = model.predict(np.zeros((len(test_ids), 1)))
        return pred_df

    feature_name, algo = config.split("__")
    feat_index, X = load_feature_matrix(feature_name)
    pop_pos = feat_index.reset_index().merge(population, on=["Molecule_Name", "inchikey"], how="left")
    if pop_pos["screen_split"].isna().any():
        raise ValueError(
            f"{config} repeat={repeat} fold={fold}: feature index did not align 1:1 "
            "against the screen population -- stopping."
        )

    pool_mask = pop_pos["screen_split"].isin(["screen_inner_train", "screen_inner_val"])
    inner_train_mask = pop_pos["screen_split"] == "screen_inner_train"
    inner_val_mask = pop_pos["screen_split"] == "screen_inner_val"
    test_positions = pop_pos.loc[pop_pos["screen_split"] == "screen_test", "index"].to_numpy()

    pred_df = test_ids.copy()
    for endpoint in REGRESSION_ENDPOINTS:
        has_label = pop_pos[endpoint].notna()
        if algo == "rf":
            idx = pop_pos.loc[pool_mask & has_label, "index"].to_numpy()
            y = pop_pos.loc[pool_mask & has_label, endpoint].to_numpy()
            model = RandomForestRegressor(random_state=seed)
            model.fit(X[idx], y)
        elif algo in ("xgboost", "lightgbm"):
            tr_idx = pop_pos.loc[inner_train_mask & has_label, "index"].to_numpy()
            va_idx = pop_pos.loc[inner_val_mask & has_label, "index"].to_numpy()
            y_tr = pop_pos.loc[inner_train_mask & has_label, endpoint].to_numpy()
            y_va = pop_pos.loc[inner_val_mask & has_label, endpoint].to_numpy()
            if algo == "xgboost":
                model = XGBRegressor(
                    early_stopping_rounds=XGB_LGBM_EARLY_STOPPING_ROUNDS, random_state=seed
                )
                model.fit(X[tr_idx], y_tr, eval_set=[(X[va_idx], y_va)], verbose=False)
            else:
                model = LGBMRegressor(random_state=seed, verbosity=-1)
                model.fit(
                    X[tr_idx], y_tr, eval_set=[(X[va_idx], y_va)],
                    callbacks=[lgb.early_stopping(stopping_rounds=XGB_LGBM_EARLY_STOPPING_ROUNDS, verbose=False)],
                )
        else:
            raise ValueError(f"unknown algorithm: {algo}")
        pred_df[endpoint] = np.asarray(model.predict(X[test_positions])).reshape(-1)

    return pred_df


def run_chemprop(config: str, repeat: int, fold: int, seed: int, logger) -> pd.DataFrame:
    repeat_col = f"repeat_{repeat}"
    population = load_screen_population(
        FOLDS_PATH, CURATED_PATH, repeat_col, fold, VAL_FRACTION, seed, logger
    )
    run_dir = CHEMPROP_RUNS_DIR / f"{config}__repeat{repeat}_fold{fold}"
    run_dir.mkdir(parents=True, exist_ok=True)

    predict_csv = build_predict_csv(population, run_dir, logger)
    train_csv = run_dir / "train_input.csv"
    build_training_csv(
        population, REGRESSION_ENDPOINTS, train_csv, logger, require_all_targets=False
    )
    run_chemprop_train(
        train_csv, REGRESSION_ENDPOINTS, run_dir, logger,
        CHEMPROP_ARCHITECTURE_ARGS[config], CHEMPROP_EPOCHS, CHEMPROP_PATIENCE, seed,
    )
    raw_pred_csv = run_dir / "raw_predictions.csv"
    run_chemprop_predict(run_dir / "model_0", predict_csv, raw_pred_csv, logger)

    expected_names = set(
        population.loc[population["screen_split"] == "screen_test", "Molecule_Name"]
    )
    verify_predictions(raw_pred_csv, expected_names, REGRESSION_ENDPOINTS, logger)
    return pd.read_csv(raw_pred_csv)


def score_and_save(config: str, repeat: int, fold: int, seed: int, pred_df: pd.DataFrame, curated: pd.DataFrame, logger) -> None:
    pred_path, score_path = result_paths(config, repeat, fold)
    PRED_DIR.mkdir(parents=True, exist_ok=True)
    SCORE_DIR.mkdir(parents=True, exist_ok=True)

    ground_truth = curated[curated["inchikey"].isin(pred_df["inchikey"])].copy()
    with per_fold_bootstrap_seed(seed):
        scored = score_activity_predictions(pred_df, ground_truth, REGRESSION_ENDPOINTS)
        scored = add_macro_endpoint(scored, REGRESSION_ENDPOINTS, ACTIVITY_METRICS)
    scored["config"] = config
    scored["repeat"] = repeat
    scored["fold"] = fold
    scored["bootstrap_seed"] = seed

    pred_df.to_csv(pred_path, index=False)
    scored.to_csv(score_path, index=False)
    logger.info(f"wrote {pred_path.name} ({len(pred_df)} rows), {score_path.name} ({len(scored)} rows)")


def parse_shard(shard: str) -> tuple[int, int]:
    i_str, n_str = shard.split("/")
    i, n = int(i_str), int(n_str)
    if not (0 <= i < n):
        raise ValueError(f"--shard must satisfy 0 <= i < n, got {shard!r}")
    return i, n


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--family", required=True,
        choices=["naive", "tabular", "chemprop_randominit", "chemprop_chemeleoninit", "all"],
        help="Which manifest rows this process handles. 'all' = every row regardless of family.",
    )
    parser.add_argument(
        "--shard", default="0/1",
        help="i/n: run only every nth manifest row (within --family), 0-indexed. Default 0/1 (all rows).",
    )
    args = parser.parse_args()
    shard_i, shard_n = parse_shard(args.shard)

    log_path = LOG_DIR / f"05_cv_comparison_{args.family}_shard{shard_i}of{shard_n}.log"
    logger = setup_logging(log_path, f"run_5x5_cv_comparison_{args.family}_{shard_i}of{shard_n}")
    logger.info("=" * 70)
    logger.info(f"starting run_5x5_cv_comparison.py --family {args.family} --shard {args.shard}")
    logger.info(f"python executable: {sys.executable}")
    logger.info(f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS')}, KMP_DUPLICATE_LIB_OK={os.environ.get('KMP_DUPLICATE_LIB_OK')}")

    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(
            f"{MANIFEST_PATH} not found -- run scripts/generate_5x5_cv_manifest.py first."
        )
    manifest = pd.read_csv(MANIFEST_PATH)
    family_rows = manifest if args.family == "all" else manifest[manifest["family"] == args.family]
    family_rows = family_rows.reset_index(drop=True)
    my_rows = family_rows.iloc[shard_i::shard_n].reset_index(drop=True)
    logger.info(
        f"manifest: {len(manifest)} total rows, {len(family_rows)} in family "
        f"{args.family!r}, {len(my_rows)} assigned to shard {shard_i}/{shard_n}"
    )

    curated = pd.read_csv(CURATED_PATH)
    logger.info(f"loaded {CURATED_PATH.name}: {curated.shape}")

    n_done, n_run, n_failed = 0, 0, 0
    for row in my_rows.itertuples(index=False):
        config, repeat, fold, seed, row_family = (
            row.config, int(row.repeat), int(row.fold), int(row.seed), row.family
        )
        if is_done(config, repeat, fold):
            n_done += 1
            continue

        logger.info(f"--- {config} repeat={repeat} fold={fold} seed={seed} ---")
        t0 = time.time()
        try:
            if row_family in ("naive", "tabular"):
                pred_df = run_tabular_or_naive(config, repeat, fold, seed, logger)
            else:
                pred_df = run_chemprop(config, repeat, fold, seed, logger)
            score_and_save(config, repeat, fold, seed, pred_df, curated, logger)
            n_run += 1
            logger.info(f"{config} repeat={repeat} fold={fold} done in {time.time() - t0:.1f}s")
        except Exception:
            n_failed += 1
            logger.error(
                f"{config} repeat={repeat} fold={fold} FAILED after {time.time() - t0:.1f}s "
                f"-- no output written, will retry on next invocation:\n{traceback.format_exc()}"
            )

    logger.info(
        f"shard {shard_i}/{shard_n} of family {args.family!r} finished: "
        f"{n_done} already done, {n_run} completed this run, {n_failed} failed"
    )


if __name__ == "__main__":
    main()
