"""Notebook 12 -- TRAIN stage: full retrain + blind predict for the Caruana-weighted
ensemble (`outputs/11b_caruana_selection/weights.csv`).

Trains, on the FULL 4,905-compound labeled training set (`train_inhibition_curated
.csv`, no held-out fold), every config that has non-zero Caruana weight for at least
one isoform AND isn't already reusable from notebook 10's existing full-data blind
predictions, then predicts on the real 750-compound blind test set
(`data/processed/test_blinded_curated.csv`). Combination (Caruana-weighted average,
NOT simple averaging), recalibration, validation and submission-file formatting are
notebook 12's job, not this script's (CLAUDE.md: scripts train, notebooks report) --
this script only trains and writes one raw prediction file per config to
`outputs/12_caruana_retrain_predict/blind_predictions/`.

TRAINING-POOL LINEAGE: unlike notebook 10 (which mixed 05's plain/unexcluded configs
with 06/07's CYP2D6-residual-excluded/tuned variants, per isoform), notebook 11b's
Caruana selection pool used 05's plain, UNEXCLUDED configs uniformly for EVERY
isoform -- including CYP2D6 itself (11_caruana_prep.ipynb pulled straight from
`outputs/05_cv_comparison/predictions/`, never from `06_outlier_check` or
`07_weighting_tuning`). So every config trained here uses the FULL, unexcluded
4,905-compound pool for every isoform, and hyperparameters are library defaults
throughout (RandomForestRegressor/LGBMRegressor/XGBRegressor with no tuning) --
CYP2D6's 07-tuned `chemeleon__rf` variant is NOT part of this pool and is not
trained here. This is a known, deliberate scope limit (see notebook 12's own intro
markdown): CYP2D6's result from this submission is a WEAKER test than the tuned
variant could give; that comparison is reserved for a later, separate submission.

DEDUPLICATION / REUSE, AT (config, isoform) GRANULARITY, NOT JUST CONFIG NAME. Tabular
configs (RF/LightGBM/XGBoost) are single-task -- one fit per isoform -- so a config's
existing `outputs/10_final_retrain_predict/blind_predictions/{config}__unexcluded.csv`
may cover some of this notebook's needed isoforms and not others (notebook 10 only
ever trained a tabular config for the isoforms where it was a step-3 winning-combo
member). Reuse is checked and reported per (config, isoform) pair -- see
`outputs/12_caruana_retrain_predict/config_provenance.csv`, which is this script's
answer to "reused or newly trained, and why" for every (config, isoform) this run
needed. Confirmed directly against the real `weights.csv` and the real
`outputs/10_final_retrain_predict/blind_predictions/*__unexcluded.csv` files (not
assumed from memory) -- as of this run:
  - `chemprop_chemeleoninit`: fully reusable, all 4 isoforms (multitask, one existing
    file already covers everything this pool needs from it).
  - `chemeleon__rf`, `ecfp4_narrow__rf`: reusable for CYP1A2 only; fresh elsewhere.
  - `ecfp4_narrow__lightgbm`: reusable for CYP1A2 and CYP3A4; fresh for CYP2C9/CYP2D6.
  - Everything else (`chemeleon__lightgbm`, `chemeleon__xgboost`, `chemprop_randominit`,
    `ecfp4_narrow__xgboost`, `mordred_pca__lightgbm`, `mordred_pca__rf`,
    `mordred_pca__xgboost`) has no existing full-data blind prediction at all --
    trained fresh here for every isoform it's needed on.
Only files literally named `*__unexcluded.csv` are ever treated as reusable --
`*__cyp2d6_residual_excluded*.csv` files exist in notebook 10's output but their
training-pool lineage (excluded pool, and for `chemeleon__rf` also tuned
hyperparameters) does not match this pool's plain/unexcluded lineage, so they are
never reused here, per this project's "exact lineage" reuse rule.

MORDRED_PCA BLIND-TEST FEATURES: `data/processed/tabular_mordred_pca.csv` is
train-only by design (`scripts/generate_mordred_pca_features.py`'s own docstring:
"the blind-test transform is explicit future work, not part of this task"). No config
in notebook 10's pool ever used `mordred_pca`, so that transform was never done. This
script does it once here, applying the ALREADY-FITTED pipeline bundle
(`models/mordred_pca_pipeline.joblib` -- impute+scale+PCA fit on the 4,905 training
compounds only) to the 750 blind-test compounds' freshly-computed Mordred descriptors
via `.transform()`, never `.fit_transform()` -- no refitting, no leakage. Written to
`outputs/12_caruana_retrain_predict/tabular_mordred_pca_blind_test.csv`, a NEW file --
the frozen `data/processed/tabular_mordred_pca.csv` is never modified.

XGBOOST SUPPORT: notebook 10's script only had RF/LightGBM branches (none of 08's
ensembles used an XGBoost member). This script adds an XGBoost branch, directly
extending 10's own LightGBM pattern (same `final_split`-derived train/val split for
early stopping, same `XGB_LGBM_EARLY_STOPPING_ROUNDS=10`) with 05's own XGBoost
hyperparameters (`XGBRegressor(early_stopping_rounds=10, random_state=seed)`, no other
kwargs -- library defaults, matching every other tabular algo's "no tuning" status
here) -- not a new modelling decision, just the same established pattern applied to a
third algorithm 10 never needed.

SEED: 42 throughout, matching 10's own precedent and reasoning for this "final,
full-data, no-CV" step category (one training run per config now, not 25 CV folds --
no fold axis for CLAUDE.md's "seeds vary per fold" rule to apply to).

REUSED, NOT RECOMPUTED: `outputs/10_final_retrain_predict/final_train_val_split.csv`
(deterministic function of curated+seed+val_fraction, 10's own script already reused
it across its own configs/reruns) and `.../predict_input.csv` (the shared 750-compound
blind predict input) are read directly from notebook 10's output, not regenerated.

EXPECT chemprop_randominit TO BE THE LONG POLE. `chemprop_chemeleoninit` is fully
reused (zero new chemprop training needed for it). `chemprop_randominit` has never
been trained on the full data before -- only per-CV-fold in `05_cv_comparison`. This
script logs its wall-clock start/end time explicitly so this project has a first real,
verified full-data timing number for a chemprop config, rather than continuing to rely
on the unverified README "~5hrs" estimate (itself for `chemeleoninit`, a specifically
called-out-as-much-slower architecture -- see README's own "28-83x slower than
random-init" note; `chemprop_randominit` is not expected to take anywhere near 5hrs,
but this script measures rather than assumes).

    cd /path/to/OpenADMET-CYP-Blind-Challenge
    caffeinate -i nohup python scripts/12_caruana_retrain_predict.py > \\
        logs/12_caruana_retrain_predict_stdout.log 2>&1 &

RESUMABILITY: each config's own blind-prediction CSV (once it covers every isoform
this run needs from it) is the completion marker; already-done configs are skipped on
re-invocation, matching 05/06/07/10's established pattern.

OUTPUTS:
    logs/12_caruana_retrain_predict.log
    outputs/12_caruana_retrain_predict/tabular_mordred_pca_blind_test.csv -- new,
        blind-test-only Mordred PCA features (see MORDRED_PCA note above)
    outputs/12_caruana_retrain_predict/chemprop_runs/chemprop_randominit/ -- chemprop
        CLI's own train_input.csv, checkpoints, config.toml
    outputs/12_caruana_retrain_predict/blind_predictions/<config>.csv -- one file per
        config this run needed (reused configs are copied in too, so notebook 12 can
        load everything from one consistent directory), columns = Molecule_Name,
        inchikey, one column per isoform this run needed from that config
    outputs/12_caruana_retrain_predict/config_provenance.csv -- (config, isoform,
        weight, source ["reused_from_10" | "trained_fresh"], blind_pred_file) --
        the full, explicit reuse/train report for notebook 12 to read directly
"""

import os

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
# OMP_NUM_THREADS intentionally left unset -- matches 05/06/07/08/10's validated ~2x
# speedup finding.

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

_env_bin = str(Path(sys.executable).parent)
_path_parts = os.environ.get("PATH", "").split(os.pathsep)
if _env_bin not in _path_parts:
    os.environ["PATH"] = os.pathsep.join([_env_bin, *_path_parts])

import joblib
import lightgbm as lgb
from lightgbm import LGBMRegressor
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor

from src.chemprop_screen import (
    run_chemprop_predict,
    run_chemprop_train,
    setup_logging,
    verify_predictions,
)
from src.features import mordred_2d_descriptors
from src.vendor.openadmet_eval.config import REGRESSION_ENDPOINTS

PROCESSED = REPO_ROOT / "data" / "processed"
LOGS = REPO_ROOT / "logs"
MODELS = REPO_ROOT / "models"
OUT10 = REPO_ROOT / "outputs" / "10_final_retrain_predict"
OUT11B = REPO_ROOT / "outputs" / "11b_caruana_selection"
OUT = REPO_ROOT / "outputs" / "12_caruana_retrain_predict"
BLIND_PRED_DIR = OUT / "blind_predictions"
CHEMPROP_RUNS_DIR = OUT / "chemprop_runs"

CURATED_PATH = PROCESSED / "train_inhibition_curated.csv"
TEST_PATH = PROCESSED / "test_blinded_curated.csv"
WEIGHTS_PATH = OUT11B / "weights.csv"
MORDRED_PIPELINE_PATH = MODELS / "mordred_pca_pipeline.joblib"
MORDRED_TRAIN_PATH = PROCESSED / "tabular_mordred_pca.csv"
MORDRED_BLIND_TEST_PATH = OUT / "tabular_mordred_pca_blind_test.csv"

SEED = 42
VAL_FRACTION = 0.15
CHEMPROP_EPOCHS = 50
CHEMPROP_PATIENCE = 5
XGB_LGBM_EARLY_STOPPING_ROUNDS = 10

CHEMPROP_ARCHITECTURE_ARGS = {
    "chemprop_randominit": [
        "--message-hidden-dim", "300", "--depth", "3", "--aggregation", "mean", "--batch-norm",
    ],
    "chemprop_chemeleoninit": ["--from-foundation", "CHEMELEON", "--multi-hot-atom-featurizer-mode", "V2"],
}
CHEMPROP_CONFIG_NAMES = set(CHEMPROP_ARCHITECTURE_ARGS)

FEATURE_FILES = {
    "ecfp4_narrow": "tabular_baseline_features.csv",
    "chemeleon": "chemeleon_embeddings.npy",
}

ISOFORMS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]  # matches REGRESSION_ENDPOINTS order
ISOFORM_COL = dict(zip(ISOFORMS, REGRESSION_ENDPOINTS))


def is_done(path: Path, needed_cols: list) -> bool:
    if not path.exists():
        return False
    try:
        df = pd.read_csv(path)
    except (pd.errors.EmptyDataError, OSError):
        return False
    if len(df) != 750 or not set(needed_cols).issubset(df.columns):
        return False
    return not df[needed_cols].isna().any().any()


def load_weights_and_needed_pairs(logger) -> dict:
    """Read `weights.csv` directly (never hardcoded/assumed from memory) and return
    {config: [isoform, ...]} for every non-zero-weight (config, isoform) pair."""
    weights = pd.read_csv(WEIGHTS_PATH)
    logger.info(f"loaded {WEIGHTS_PATH}: {weights.shape}")
    expected_isoforms = set(ISOFORMS)
    if set(weights["isoform"].unique()) != expected_isoforms:
        raise ValueError(f"weights.csv isoforms {sorted(weights['isoform'].unique())} != expected {sorted(expected_isoforms)}")

    nonzero = weights[weights["weight"] > 0]
    needed: dict[str, list[str]] = {}
    for row in nonzero.itertuples(index=False):
        needed.setdefault(row.config, []).append(row.isoform)

    logger.info(f"non-zero-weight (config, isoform) pairs: {int(nonzero.shape[0])} total, {len(needed)} distinct configs")
    for config, isos in sorted(needed.items()):
        logger.info(f"  {config}: needed for {sorted(isos)}")
    return needed


def build_mordred_blind_test_features(test_df: pd.DataFrame, logger) -> None:
    """Apply the ALREADY-FITTED mordred_pca pipeline (train-fit only) to the blind
    test set's freshly-computed Mordred descriptors. `.transform()` only, never
    `.fit_transform()` -- see module docstring."""
    if MORDRED_BLIND_TEST_PATH.exists():
        existing = pd.read_csv(MORDRED_BLIND_TEST_PATH)
        if len(existing) == 750:
            logger.info(f"[skip, already done] {MORDRED_BLIND_TEST_PATH.name} ({existing.shape})")
            return

    bundle = joblib.load(MORDRED_PIPELINE_PATH)
    pipeline, dropped, descriptor_cols = bundle["pipeline"], bundle["dropped_descriptors"], bundle["descriptor_columns"]
    train_n_components = pd.read_csv(MORDRED_TRAIN_PATH, nrows=0).columns
    n_components = sum(1 for c in train_n_components if c.startswith("mordred_pca_"))
    logger.info(
        f"mordred_pca pipeline bundle loaded from {MORDRED_PIPELINE_PATH.name}: "
        f"{len(descriptor_cols)} descriptor columns (after {len(dropped)} always-fail dropped), "
        f"transforming to match training file's {n_components} saved components"
    )

    t0 = time.time()
    raw = mordred_2d_descriptors(test_df["canonical_smiles"].tolist())
    logger.info(f"mordred_2d_descriptors on {len(test_df)} blind-test compounds: {raw.shape} in {time.time() - t0:.1f}s")
    if not set(descriptor_cols).issubset(raw.columns):
        missing = set(descriptor_cols) - set(raw.columns)
        raise ValueError(f"blind-test Mordred descriptors missing {len(missing)} column(s) the fitted pipeline expects: {sorted(missing)[:10]}")
    features = raw[descriptor_cols]  # exact same column set/order the pipeline was fit on

    transformed = pipeline.transform(features.to_numpy())  # transform ONLY -- never refit
    pc_cols = [f"mordred_pca_{i:04d}" for i in range(n_components)]
    out_df = pd.DataFrame(transformed[:, :n_components], columns=pc_cols)
    out_df.insert(0, "split", "test")
    out_df.insert(0, "inchikey", test_df["inchikey"].values)
    out_df.insert(0, "Molecule_Name", test_df["Molecule_Name"].values)

    OUT.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(MORDRED_BLIND_TEST_PATH, index=False)
    logger.info(f"wrote {MORDRED_BLIND_TEST_PATH} ({out_df.shape}) -- transform only, pipeline not refit")


def load_feature_matrix(name: str) -> tuple[pd.DataFrame, np.ndarray]:
    """Load a full feature matrix (both the 4,905 train rows and the 750 blind-test
    rows), tagged with `feat_pos`. `ecfp4_narrow`/`chemeleon` already carry both splits
    (04b/05's own precedent); `mordred_pca` is assembled here from the frozen
    train-only file plus this script's freshly-built blind-test file."""
    if name == "mordred_pca":
        train_df = pd.read_csv(MORDRED_TRAIN_PATH)
        test_df = pd.read_csv(MORDRED_BLIND_TEST_PATH)
        df = pd.concat([train_df, test_df], ignore_index=True)
        index_df = df[["Molecule_Name", "inchikey", "split"]].reset_index(drop=True)
        X = df.drop(columns=["Molecule_Name", "inchikey", "split"]).to_numpy(dtype=np.float64)
    else:
        path = PROCESSED / FEATURE_FILES[name]
        if path.suffix == ".csv":
            df = pd.read_csv(path)
            index_df = df[["Molecule_Name", "inchikey", "split"]].reset_index(drop=True)
            X = df.drop(columns=["Molecule_Name", "inchikey", "split"]).to_numpy(dtype=np.float64)
        else:
            X = np.load(path)
            index_df = pd.read_csv(path.with_name(f"{path.stem}_index.csv"))[["Molecule_Name", "inchikey", "split"]]
    index_df = index_df.reset_index().rename(columns={"index": "feat_pos"})
    return index_df, X


def check_reuse(config: str, isoform: str) -> Path | None:
    """A config is reusable for one isoform only if notebook 10's *unexcluded*
    (plain, no CYP2D6 exclusion, no tuning) blind prediction file exists AND already
    contains that isoform's column, fully populated for all 750 rows. Excluded/tuned
    lineage files are never matched here -- see module docstring."""
    candidate = OUT10 / "blind_predictions" / f"{config}__unexcluded.csv"
    if not candidate.exists():
        return None
    df = pd.read_csv(candidate)
    col = ISOFORM_COL[isoform]
    if col not in df.columns or len(df) != 750 or df[col].isna().any():
        return None
    return candidate


def fit_predict_tabular(
    feature_name: str, algo: str, target_col: str, curated: pd.DataFrame, test_df: pd.DataFrame,
    final_split: pd.DataFrame, seed: int, logger,
) -> pd.DataFrame:
    feat_index, X = load_feature_matrix(feature_name)
    train_feat = feat_index[feat_index["split"] == "train"]
    test_feat = feat_index[feat_index["split"] == "test"]

    train_merged = train_feat.merge(curated[["Molecule_Name", "inchikey", target_col]], on=["Molecule_Name", "inchikey"], how="left")
    if len(train_merged) != len(train_feat):
        raise ValueError(f"{feature_name}: train feature index did not align 1:1 against curated -- stopping.")
    test_merged = test_feat.merge(test_df[["Molecule_Name", "inchikey"]], on=["Molecule_Name", "inchikey"], how="inner")
    if len(test_merged) != 750:
        raise ValueError(f"{feature_name}: expected 750 blind-test rows after alignment, got {len(test_merged)}.")

    usable = train_merged[target_col].notna()  # no exclusion -- this pool is fully unexcluded, see module docstring
    logger.info(f"{feature_name}__{algo} [{target_col}]: labeled rows={int(usable.sum())} (no exclusion applied -- unexcluded pool)")

    if algo == "rf":
        idx = train_merged.loc[usable, "feat_pos"].to_numpy()
        y = train_merged.loc[usable, target_col].to_numpy()
        model = RandomForestRegressor(random_state=seed)
        model.fit(X[idx], y)
    elif algo in ("lightgbm", "xgboost"):
        split_merged = train_merged.merge(final_split[["Molecule_Name", "final_submission_split"]], on="Molecule_Name", how="left")
        tr_mask = usable & (split_merged["final_submission_split"] == "final_train")
        va_mask = usable & (split_merged["final_submission_split"] == "final_val")
        logger.info(f"{feature_name}__{algo} [{target_col}]: train={int(tr_mask.sum())}, val={int(va_mask.sum())} (early-stopping split)")
        tr_idx = split_merged.loc[tr_mask, "feat_pos"].to_numpy()
        va_idx = split_merged.loc[va_mask, "feat_pos"].to_numpy()
        y_tr = split_merged.loc[tr_mask, target_col].to_numpy()
        y_va = split_merged.loc[va_mask, target_col].to_numpy()
        if algo == "xgboost":
            model = XGBRegressor(early_stopping_rounds=XGB_LGBM_EARLY_STOPPING_ROUNDS, random_state=seed)
            model.fit(X[tr_idx], y_tr, eval_set=[(X[va_idx], y_va)], verbose=False)
        else:
            model = LGBMRegressor(random_state=seed, verbosity=-1)
            model.fit(
                X[tr_idx], y_tr, eval_set=[(X[va_idx], y_va)],
                callbacks=[lgb.early_stopping(stopping_rounds=XGB_LGBM_EARLY_STOPPING_ROUNDS, verbose=False)],
            )
    else:
        raise ValueError(f"unknown algorithm: {algo}")

    test_pos = test_merged["feat_pos"].to_numpy()
    preds = np.asarray(model.predict(X[test_pos])).reshape(-1)
    out = test_merged[["Molecule_Name", "inchikey"]].copy()
    out[target_col] = preds
    return out


def run_tabular_configs(needed: dict, curated: pd.DataFrame, test_df: pd.DataFrame, final_split: pd.DataFrame, logger) -> list[dict]:
    provenance = []
    tabular_configs = {c: isos for c, isos in needed.items() if c not in CHEMPROP_CONFIG_NAMES}
    for config, isos in sorted(tabular_configs.items()):
        feature_name, algo = config.split("__")
        out_path = BLIND_PRED_DIR / f"{config}.csv"

        reused_cols, fresh_isos = {}, []
        for isoform in isos:
            reuse_path = check_reuse(config, isoform)
            if reuse_path is not None:
                reused_cols[isoform] = reuse_path
            else:
                fresh_isos.append(isoform)

        needed_cols = [ISOFORM_COL[i] for i in isos]
        if is_done(out_path, needed_cols):
            logger.info(f"[skip, already done] {config} ({sorted(isos)})")
            for isoform in isos:
                provenance.append({
                    "config": config, "isoform": isoform,
                    "source": "reused_from_10" if isoform in reused_cols else "trained_fresh",
                    "blind_pred_file": out_path.name,
                })
            continue

        logger.info(f"--- tabular config: {config} (reused={sorted(reused_cols)}, fresh={sorted(fresh_isos)}) ---")
        t0 = time.time()
        combined = test_df[["Molecule_Name", "inchikey"]].copy()
        for isoform, reuse_path in reused_cols.items():
            col = ISOFORM_COL[isoform]
            src_df = pd.read_csv(reuse_path)[["Molecule_Name", "inchikey", col]]
            combined = combined.merge(src_df, on=["Molecule_Name", "inchikey"], how="left")
            logger.info(f"  {config} [{isoform}]: reused from {reuse_path.relative_to(REPO_ROOT)}")
        for isoform in fresh_isos:
            col = ISOFORM_COL[isoform]
            pred_df = fit_predict_tabular(feature_name, algo, col, curated, test_df, final_split, SEED, logger)
            combined = combined.merge(pred_df, on=["Molecule_Name", "inchikey"], how="left")

        if len(combined) != 750 or combined[needed_cols].isna().any().any():
            raise ValueError(f"{config}: incomplete after combining reused+fresh columns -- stopping.")
        BLIND_PRED_DIR.mkdir(parents=True, exist_ok=True)
        combined.to_csv(out_path, index=False)
        logger.info(f"wrote {out_path} ({len(combined)} rows, cols={needed_cols}) in {time.time() - t0:.1f}s")

        for isoform in isos:
            provenance.append({
                "config": config, "isoform": isoform,
                "source": "reused_from_10" if isoform in reused_cols else "trained_fresh",
                "blind_pred_file": out_path.name,
            })
    return provenance


def run_chemprop_config(config: str, isos: list, curated: pd.DataFrame, test_df: pd.DataFrame, final_split: pd.DataFrame, predict_csv: Path, logger) -> list[dict]:
    out_path = BLIND_PRED_DIR / f"{config}.csv"
    needed_cols = [ISOFORM_COL[i] for i in isos]

    reuse_path = OUT10 / "blind_predictions" / f"{config}__unexcluded.csv"
    if reuse_path.exists():
        df = pd.read_csv(reuse_path)
        if len(df) == 750 and set(needed_cols).issubset(df.columns) and not df[needed_cols].isna().any().any():
            BLIND_PRED_DIR.mkdir(parents=True, exist_ok=True)
            df[["Molecule_Name", "inchikey", *needed_cols]].to_csv(out_path, index=False)
            logger.info(f"[fully reused from notebook 10] {config}: copied {reuse_path.relative_to(REPO_ROOT)} -> {out_path}")
            return [{"config": config, "isoform": i, "source": "reused_from_10", "blind_pred_file": out_path.name} for i in isos]

    if is_done(out_path, needed_cols):
        logger.info(f"[skip, already done] {config} ({sorted(isos)})")
        return [{"config": config, "isoform": i, "source": "trained_fresh", "blind_pred_file": out_path.name} for i in isos]

    logger.info(f"--- chemprop config: {config} (fresh, no reusable lineage found) -- isoforms: {sorted(isos)} ---")
    t0 = time.time()
    t0_str = time.strftime("%Y-%m-%d %H:%M:%S")
    logger.info(f"{config}: training start = {t0_str}")

    population = curated.copy()  # unexcluded pool -- no masking, matches this pool's lineage throughout
    population = population.merge(final_split[["Molecule_Name", "final_submission_split"]], on="Molecule_Name", how="left")
    if population["final_submission_split"].isna().any():
        raise ValueError(f"{config}: final_submission_split did not cover every curated row -- stopping.")
    split_map = {"final_train": "train", "final_val": "val"}
    population["chemprop_split"] = population["final_submission_split"].map(split_map)

    run_dir = CHEMPROP_RUNS_DIR / config
    run_dir.mkdir(parents=True, exist_ok=True)
    train_csv = run_dir / "train_input.csv"
    cols = ["canonical_smiles", *REGRESSION_ENDPOINTS, "chemprop_split"]
    population[cols].to_csv(train_csv, index=False)
    logger.info(f"wrote {train_csv} ({len(population)} rows, chemprop_split counts: {population['chemprop_split'].value_counts().to_dict()})")

    run_chemprop_train(
        train_csv, REGRESSION_ENDPOINTS, run_dir, logger,
        CHEMPROP_ARCHITECTURE_ARGS[config], CHEMPROP_EPOCHS, CHEMPROP_PATIENCE, SEED,
    )

    raw_pred_csv = run_dir / "raw_predictions.csv"
    run_chemprop_predict(run_dir / "model_0", predict_csv, raw_pred_csv, logger)
    expected_names = set(test_df["Molecule_Name"])
    verify_predictions(raw_pred_csv, expected_names, REGRESSION_ENDPOINTS, logger)

    BLIND_PRED_DIR.mkdir(parents=True, exist_ok=True)
    raw_pred = pd.read_csv(raw_pred_csv)
    raw_pred.to_csv(out_path, index=False)

    t1_str = time.strftime("%Y-%m-%d %H:%M:%S")
    elapsed = time.time() - t0
    logger.info(f"{config}: training end = {t1_str}, elapsed = {elapsed:.1f}s ({elapsed / 60:.1f} min, {elapsed / 3600:.2f} hr)")
    logger.info(f"wrote {out_path} ({len(raw_pred)} rows)")
    return [{"config": config, "isoform": i, "source": "trained_fresh", "blind_pred_file": out_path.name} for i in isos]


def main() -> None:
    logger = setup_logging(LOGS / "12_caruana_retrain_predict.log", "12_caruana_retrain_predict")
    script_start = time.time()
    logger.info("=" * 70)
    logger.info("starting scripts/12_caruana_retrain_predict.py")

    OUT.mkdir(parents=True, exist_ok=True)
    BLIND_PRED_DIR.mkdir(parents=True, exist_ok=True)
    CHEMPROP_RUNS_DIR.mkdir(parents=True, exist_ok=True)

    needed = load_weights_and_needed_pairs(logger)

    curated = pd.read_csv(CURATED_PATH)
    logger.info(f"loaded {CURATED_PATH.name}: {curated.shape}")
    test_df = pd.read_csv(TEST_PATH)
    logger.info(f"loaded {TEST_PATH.name}: {test_df.shape}")
    if len(test_df) != 750:
        raise ValueError(f"expected 750 blinded test compounds, got {len(test_df)} -- stopping.")

    split_path = OUT10 / "final_train_val_split.csv"
    final_split = pd.read_csv(split_path)
    logger.info(f"reusing existing {split_path} ({len(final_split)} rows) from notebook 10 -- deterministic function of (curated, seed={SEED}, val_fraction={VAL_FRACTION}), not regenerated")

    predict_csv = OUT10 / "predict_input.csv"
    logger.info(f"reusing existing {predict_csv} from notebook 10 -- deterministic function of test_blinded_curated.csv, not regenerated")

    if "mordred_pca" in {c.split("__")[0] for c in needed if c not in CHEMPROP_CONFIG_NAMES}:
        logger.info("=" * 70)
        logger.info("stage: mordred_pca blind-test feature transform (one-time, pipeline.transform() only)")
        build_mordred_blind_test_features(test_df, logger)

    logger.info("=" * 70)
    logger.info("stage: tabular configs (RF/LightGBM/XGBoost, fast)")
    provenance = run_tabular_configs(needed, curated, test_df, final_split, logger)

    logger.info("=" * 70)
    logger.info("stage: chemprop configs (slow -- wall time logged explicitly per config)")
    for config, isos in sorted(needed.items()):
        if config in CHEMPROP_CONFIG_NAMES:
            provenance.extend(run_chemprop_config(config, isos, curated, test_df, final_split, predict_csv, logger))

    provenance_df = pd.DataFrame(provenance).sort_values(["isoform", "config"]).reset_index(drop=True)
    provenance_path = OUT / "config_provenance.csv"
    provenance_df.to_csv(provenance_path, index=False)
    logger.info(f"wrote {provenance_path} ({provenance_df.shape})")
    logger.info(f"provenance summary: {provenance_df['source'].value_counts().to_dict()}")

    total_elapsed = time.time() - script_start
    logger.info("=" * 70)
    logger.info(f"total script wall time: {total_elapsed:.1f}s ({total_elapsed / 60:.1f} min)")
    logger.info("done.")


if __name__ == "__main__":
    main()
