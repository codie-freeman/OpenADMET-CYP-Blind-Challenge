"""Roadmap Step "final retrain and predict" (notebook 10) -- TRAIN stage.

Retrains, on the FULL 4,905-compound labeled training set (`train_inhibition_curated
.csv`, no held-out fold -- CV's job of model/ensemble selection is already done by
notebooks 05-08), every distinct model config referenced by any isoform's step-3
winning source (`outputs/08_ensemble_selection/step3_summary.csv`), then predicts on
the real 750-compound blind test set (`data/processed/test_blinded_curated.csv`).
Combination (simple-averaging per isoform), recalibration, and final-submission
formatting/validation are notebook 10's job, not this script's (CLAUDE.md: scripts
train, notebooks report) -- this script only trains and writes one raw prediction
file per deduplicated config to `outputs/10_final_retrain_predict/blind_predictions/`.

DEDUPLICATION KEY IS (config, training-pool lineage), NOT JUST CONFIG NAME. Reading
`step3_summary.csv`'s `winning_combo`/`single_best_label` columns directly (the
literal "step-3 winning sources", not the broader candidate pools 08 also tested)
gives 4 isoforms' worth of members:

    CYP1A2: chemprop_chemeleoninit + ecfp4_narrow__rf + chemeleon__rf + ecfp4_narrow__lightgbm
    CYP2C9: chemprop_chemeleoninit (single model, no ensemble)
    CYP2D6: chemeleon__rf (07 tuned) + chemprop_chemeleoninit
    CYP3A4: chemprop_chemeleoninit + ecfp4_narrow__lightgbm

Cross-referencing each occurrence against 08's own saved pool-source definitions
(`notebooks/08_ensemble_selection.ipynb`, cells defining CYP1A2_POOL/CYP2C9 single-
best load/CYP2D6_POOL/CYP3A4_POOL -- read directly, not from memory, confirming the
task's own instruction not to assume symmetry across isoforms) shows the lineage is
NOT symmetric for either name shared across isoforms:

    chemprop_chemeleoninit:
      - CYP1A2, CYP2C9, CYP3A4 all use 05's ORIGINAL run (`("05", "chemprop_
        chemeleoninit")` / `("05", ...)` via `load_existing_fold_scores("05", ...)`
        for CYP2C9) -- full labeled pool, no exclusion.
      - CYP2D6 uses 06's CYP2D6-residual-excluded run (`("06", "residual__
        chemprop_chemeleoninit")`) -- SAME multitask architecture, but CYP2D6's
        label is NaN-masked for the 75 compounds in `outputs/06_outlier_check/
        flagged_compounds/criterion1_residual_flagged.csv` before training (06's own
        confirmed exclusion design: train-only, label-masked not row-dropped, so the
        other 3 isoforms' training contribution from those 75 compounds is untouched).
      => two distinct chemprop trainings, not one.
    chemeleon__rf:
      - CYP1A2 uses 05's ORIGINAL run (`("05", "chemeleon__rf")`) -- sklearn RF
        defaults, no exclusion.
      - CYP2D6 uses 07's TUNED run (`("07", "tuning__chemeleon__rf")`) -- hyper-
        parameters n_estimators=1000, max_depth=10, min_samples_leaf=5 (confirmed
        against 07's own written outcome in `notebooks/07_weighting_tuning.ipynb`,
        cells 20/22/25 -- the winning grid point, not retyped from the task prompt
        alone), trained on 06's CYP2D6-residual-excluded pool (07's tuning was built
        directly on top of 06's excluded baseline -- confirmed via 07's own
        parameter-table cell: "fit scope: refit per fold, on that fold's own
        post-exclusion training pool").
      => two distinct chemeleon__rf trainings, not one.
    ecfp4_narrow__lightgbm:
      - CYP1A2 AND CYP3A4 both use 05's ORIGINAL run (`("05", "ecfp4_narrow__
        lightgbm")` in both CYP1A2_POOL and CYP3A4_POOL) -- same lineage, so this one
        genuinely IS trained once and its predictions reused for both isoforms.
    ecfp4_narrow__rf:
      - CYP1A2 only (`("05", "ecfp4_narrow__rf")`).

NOT RETRAINED: `chemeleon__lightgbm`. It was a candidate in CYP2D6's step-3 pool
(`CYP2D6_POOL` in notebook 08) but is NOT a member of CYP2D6's winning combo
(`chemeleon__rf+chemprop_chemeleoninit`, per `step3_summary.csv`'s own
`winning_combo` field) or any other isoform's winning source. The task prompt's own
illustrative config list names it, but `step3_summary.csv` -- the file the task
explicitly says to read the winning sources from -- does not. Following the file
over the prompt's illustrative list here: retraining a config with zero winning-combo
membership would be exactly the kind of "add a model variant while you're at it"
CLAUDE.md asks not to do. Flagged prominently (this docstring + printed at runtime)
rather than silently decided, per the task's own "confirm explicitly, don't assume"
instruction for pool-lineage calls generally.

TRAINING-POOL MECHANICS, per algorithm family (mirrors 05/06's own established
patterns, reused not reinvented):
  - chemprop (multitask): trained on ALL 4,905 compounds every time (chemprop
    NaN-masks any missing per-task target on its own -- `require_all_targets=False`
    precedent from `src.chemprop_screen.build_training_csv` / 06's
    `retrain_chemprop_cyp2d6`). The "excluded" chemprop run additionally NaN-masks
    CYP2D6's label (not the row) for the 75 flagged compounds before writing the
    training CSV -- other 3 isoforms' labels for those compounds are untouched,
    exactly mirroring 06's `retrain_chemprop_cyp2d6`.
  - RF (single-task, no early-stopping concept): pools ALL available labeled rows
    for that target (05's own pattern: `screen_inner_train`+`screen_inner_val`
    pooled together for `algo == "rf"`) -- here, that pooling is simply "every
    curated row with a non-null label for that target", since there is no held-out
    fold this time. The excluded RF (CYP2D6, tuned) additionally drops the 75
    flagged compounds from that pool (06's own `retrain_tabular_cyp2d6` pattern:
    `has_label & ~is_excluded`).
  - LightGBM (single-task, early-stopping via `eval_set`): needs an internal
    train/val split even though there's no CV fold to draw one from. Reuses
    `src.features.assign_final_submission_split` (already established in `04b`/
    `scripts/train_final_submission_multitask.py` for exactly this "final, full-
    data, no-CV" early-stopping-split role for chemprop) for LightGBM's `eval_set`
    too -- same purpose (an internal val slice for the one final run), same
    function, not a new invention. `early_stopping_rounds=10` matches 05's own
    `XGB_LGBM_EARLY_STOPPING_ROUNDS` precedent.

SEED: 42 throughout (chemprop `--data-seed`/`--pytorch-seed`, `assign_final_
submission_split`, every sklearn/LightGBM `random_state`) -- matches `train_final_
submission_multitask.py`'s own precedent for this same "final, full-data, no-CV"
step category (05/06/07's per-(config,repeat,fold) varying manifest seeds don't
apply here: there is exactly one training run per config now, not 25, so CLAUDE.md's
"seeds must vary per CV fold" rule has no fold axis to vary across -- the closest
real precedent, 04b, used one fixed seed for its one final run, and this step
follows that same precedent for the same reason).

OMP_NUM_THREADS: deliberately left unset, matching 05/06/07/08's validated ~2x
speedup finding (see `run_5x5_cv_comparison.py`'s own dated deviation note) rather
than reverting to 04b's older `OMP_NUM_THREADS=1`. Purely a wall-clock setting --
06's Section 6 determinism check found chemprop training bit-identical regardless,
so this has no effect on correctness/reproducibility.

EXPECT THIS TO TAKE SEVERAL HOURS for the two chemprop configs (comparable to or
longer than 04b's own multi-hour final run, on similarly-sized ~4,900-compound
training pools) -- an overnight/background run is expected, matching every prior
Chemprop-training script's own precedent in this repo.

    cd /path/to/OpenADMET-CYP-Blind-Challenge
    caffeinate -i nohup python scripts/10_final_retrain_predict.py > \\
        logs/10_final_retrain_predict_stdout.log 2>&1 &

RESUMABILITY: each config's own blind-prediction CSV is the completion marker (exists
+ non-empty); already-done configs are skipped on re-invocation, matching 05/06/07's
established pattern.

OUTPUTS:
    logs/10_final_retrain_predict.log
    outputs/10_final_retrain_predict/final_train_val_split.csv   -- Molecule_Name,
                                                                      inchikey,
                                                                      final_submission_split
                                                                      (assign_final_
                                                                      submission_split,
                                                                      independent of
                                                                      cv_folds.csv)
    outputs/10_final_retrain_predict/chemprop_runs/<config>/      -- chemprop CLI's own
                                                                      train_input.csv,
                                                                      checkpoints,
                                                                      config.toml
    outputs/10_final_retrain_predict/predict_input.csv             -- shared 750-compound
                                                                      blind predict input
    outputs/10_final_retrain_predict/blind_predictions/<config>.csv -- one file per
                                                                      deduplicated config
    outputs/10_final_retrain_predict/config_manifest.csv            -- isoform ->
                                                                      (config, blind
                                                                      prediction file,
                                                                      column) lookup, for
                                                                      notebook 10 to read
                                                                      directly rather than
                                                                      re-deriving this
                                                                      lineage itself
"""

import os

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
# OMP_NUM_THREADS intentionally left unset -- see module docstring.

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

import lightgbm as lgb
from lightgbm import LGBMRegressor
from sklearn.ensemble import RandomForestRegressor

from src.chemprop_screen import (
    run_chemprop_predict,
    run_chemprop_train,
    setup_logging,
    verify_predictions,
)
from src.features import assign_final_submission_split
from src.vendor.openadmet_eval.config import REGRESSION_ENDPOINTS

PROCESSED = REPO_ROOT / "data" / "processed"
LOGS = REPO_ROOT / "logs"
OUT = REPO_ROOT / "outputs" / "10_final_retrain_predict"
BLIND_PRED_DIR = OUT / "blind_predictions"
CHEMPROP_RUNS_DIR = OUT / "chemprop_runs"

CURATED_PATH = PROCESSED / "train_inhibition_curated.csv"
TEST_PATH = PROCESSED / "test_blinded_curated.csv"
FLAGGED_PATH = (
    REPO_ROOT / "outputs" / "06_outlier_check" / "flagged_compounds" / "criterion1_residual_flagged.csv"
)

SEED = 42
VAL_FRACTION = 0.15
CHEMPROP_EPOCHS = 50
CHEMPROP_PATIENCE = 5
CHEMPROP_ARCHITECTURE_ARGS = ["--from-foundation", "CHEMELEON", "--multi-hot-atom-featurizer-mode", "V2"]
XGB_LGBM_EARLY_STOPPING_ROUNDS = 10

FEATURE_FILES = {
    "ecfp4_narrow": "tabular_baseline_features.csv",
    "chemeleon": "chemeleon_embeddings.npy",
}

CYP1A2_COL, CYP2C9_COL, CYP2D6_COL, CYP3A4_COL = REGRESSION_ENDPOINTS

CHEMPROP_CONFIGS = [
    {"name": "chemprop_chemeleoninit__unexcluded", "mask_cyp2d6": False},
    {"name": "chemprop_chemeleoninit__cyp2d6_residual_excluded", "mask_cyp2d6": True},
]

TABULAR_CONFIGS = [
    {
        "name": "ecfp4_narrow__rf__unexcluded", "feature": "ecfp4_narrow", "algo": "rf",
        "targets": [CYP1A2_COL], "excluded": False, "hyperparams": {},
    },
    {
        "name": "chemeleon__rf__unexcluded", "feature": "chemeleon", "algo": "rf",
        "targets": [CYP1A2_COL], "excluded": False, "hyperparams": {},
    },
    {
        "name": "ecfp4_narrow__lightgbm__unexcluded", "feature": "ecfp4_narrow", "algo": "lightgbm",
        "targets": [CYP1A2_COL, CYP3A4_COL], "excluded": False, "hyperparams": {},
    },
    {
        "name": "chemeleon__rf__cyp2d6_residual_excluded_tuned", "feature": "chemeleon", "algo": "rf",
        "targets": [CYP2D6_COL], "excluded": True,
        "hyperparams": {"n_estimators": 1000, "max_depth": 10, "min_samples_leaf": 5},
    },
]

# isoform -> [(step3 combo-member label, blind-prediction filename)] -- derived from
# step3_summary.csv's winning_combo/single_best_label columns, cross-referenced
# against 08's own saved pool-source lineage (see module docstring). Written to disk
# as config_manifest.csv so notebook 10 reads this mapping directly rather than
# re-deriving/hardcoding it a second time.
ISOFORM_COMBO_MEMBERS = {
    "CYP1A2": [
        ("chemprop_chemeleoninit", "chemprop_chemeleoninit__unexcluded.csv", CYP1A2_COL),
        ("ecfp4_narrow__rf", "ecfp4_narrow__rf__unexcluded.csv", CYP1A2_COL),
        ("chemeleon__rf", "chemeleon__rf__unexcluded.csv", CYP1A2_COL),
        ("ecfp4_narrow__lightgbm", "ecfp4_narrow__lightgbm__unexcluded.csv", CYP1A2_COL),
    ],
    "CYP2C9": [
        ("chemprop_chemeleoninit", "chemprop_chemeleoninit__unexcluded.csv", CYP2C9_COL),
    ],
    "CYP2D6": [
        ("chemeleon__rf (07 tuned)", "chemeleon__rf__cyp2d6_residual_excluded_tuned.csv", CYP2D6_COL),
        ("chemprop_chemeleoninit", "chemprop_chemeleoninit__cyp2d6_residual_excluded.csv", CYP2D6_COL),
    ],
    "CYP3A4": [
        ("chemprop_chemeleoninit", "chemprop_chemeleoninit__unexcluded.csv", CYP3A4_COL),
        ("ecfp4_narrow__lightgbm", "ecfp4_narrow__lightgbm__unexcluded.csv", CYP3A4_COL),
    ],
}


def is_done(path: Path) -> bool:
    if not path.exists():
        return False
    try:
        return len(pd.read_csv(path)) > 0
    except (pd.errors.EmptyDataError, OSError):
        return False


def load_feature_matrix(name: str) -> tuple[pd.DataFrame, np.ndarray]:
    """Load a full feature matrix (both the 4,905 train rows and the 750 blind-test
    rows), tagged with `feat_pos` -- the row's fixed position in `X` -- so callers can
    filter/reorder the index freely (by `split`, by label availability, by exclusion)
    while still indexing correctly into `X`.
    """
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


def fit_predict_tabular(
    feature_name: str, algo: str, target_col: str, hyperparams: dict,
    curated: pd.DataFrame, test_df: pd.DataFrame, excluded_inchikeys: set,
    final_split: pd.DataFrame, seed: int, logger,
) -> pd.DataFrame:
    feat_index, X = load_feature_matrix(feature_name)
    train_feat = feat_index[feat_index["split"] == "train"]
    test_feat = feat_index[feat_index["split"] == "test"]

    train_merged = train_feat.merge(
        curated[["Molecule_Name", "inchikey", target_col]], on=["Molecule_Name", "inchikey"], how="left"
    )
    if len(train_merged) != len(train_feat):
        raise ValueError(f"{feature_name}: train feature index did not align 1:1 against curated -- stopping.")
    test_merged = test_feat.merge(test_df[["Molecule_Name", "inchikey"]], on=["Molecule_Name", "inchikey"], how="inner")
    if len(test_merged) != 750:
        raise ValueError(f"{feature_name}: expected 750 blind-test rows after alignment, got {len(test_merged)}.")

    has_label = train_merged[target_col].notna()
    is_excluded = train_merged["inchikey"].isin(excluded_inchikeys)
    usable = has_label & ~is_excluded
    logger.info(
        f"{feature_name}__{algo} [{target_col}]: labeled rows={int(has_label.sum())}, "
        f"excluded={int((has_label & is_excluded).sum())}, usable for training={int(usable.sum())}"
    )

    if algo == "rf":
        idx = train_merged.loc[usable, "feat_pos"].to_numpy()
        y = train_merged.loc[usable, target_col].to_numpy()
        model = RandomForestRegressor(random_state=seed, **hyperparams)
        model.fit(X[idx], y)
    elif algo == "lightgbm":
        split_merged = train_merged.merge(final_split[["Molecule_Name", "final_submission_split"]], on="Molecule_Name", how="left")
        tr_mask = usable & (split_merged["final_submission_split"] == "final_train")
        va_mask = usable & (split_merged["final_submission_split"] == "final_val")
        logger.info(f"{feature_name}__{algo} [{target_col}]: train={int(tr_mask.sum())}, val={int(va_mask.sum())} (early-stopping split)")
        tr_idx = split_merged.loc[tr_mask, "feat_pos"].to_numpy()
        va_idx = split_merged.loc[va_mask, "feat_pos"].to_numpy()
        y_tr = split_merged.loc[tr_mask, target_col].to_numpy()
        y_va = split_merged.loc[va_mask, target_col].to_numpy()
        model = LGBMRegressor(random_state=seed, verbosity=-1, **hyperparams)
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


def run_tabular_configs(curated: pd.DataFrame, test_df: pd.DataFrame, excluded_inchikeys: set, final_split: pd.DataFrame, logger) -> None:
    for cfg in TABULAR_CONFIGS:
        out_path = BLIND_PRED_DIR / f"{cfg['name']}.csv"
        if is_done(out_path):
            logger.info(f"[skip, already done] {cfg['name']}")
            continue
        logger.info(f"--- tabular config: {cfg['name']} ---")
        t0 = time.time()
        excluded = excluded_inchikeys if cfg["excluded"] else set()
        combined = None
        for target_col in cfg["targets"]:
            pred_df = fit_predict_tabular(
                cfg["feature"], cfg["algo"], target_col, cfg["hyperparams"],
                curated, test_df, excluded, final_split, SEED, logger,
            )
            combined = pred_df if combined is None else combined.merge(pred_df, on=["Molecule_Name", "inchikey"], how="inner")
        if len(combined) != 750:
            raise ValueError(f"{cfg['name']}: expected 750 rows, got {len(combined)} -- stopping.")
        BLIND_PRED_DIR.mkdir(parents=True, exist_ok=True)
        combined.to_csv(out_path, index=False)
        logger.info(f"wrote {out_path} ({len(combined)} rows) in {time.time() - t0:.1f}s")


def run_chemprop_config(cfg: dict, curated: pd.DataFrame, test_df: pd.DataFrame, excluded_inchikeys: set, final_split: pd.DataFrame, predict_csv: Path, logger) -> None:
    out_path = BLIND_PRED_DIR / f"{cfg['name']}.csv"
    if is_done(out_path):
        logger.info(f"[skip, already done] {cfg['name']}")
        return
    logger.info(f"--- chemprop config: {cfg['name']} ---")
    t0 = time.time()

    population = curated.copy()
    if cfg["mask_cyp2d6"]:
        n_before = int(population[CYP2D6_COL].notna().sum())
        mask = population["inchikey"].isin(excluded_inchikeys)
        population.loc[mask, CYP2D6_COL] = np.nan
        n_after = int(population[CYP2D6_COL].notna().sum())
        logger.info(
            f"{cfg['name']}: CYP2D6-labeled rows before exclusion-mask={n_before}, "
            f"after={n_after} (masked {n_before - n_after}); other 3 endpoints' labels untouched"
        )
    else:
        logger.info(f"{cfg['name']}: no exclusion mask applied (05's original, unexcluded lineage)")

    population = population.merge(final_split[["Molecule_Name", "final_submission_split"]], on="Molecule_Name", how="left")
    if population["final_submission_split"].isna().any():
        raise ValueError(f"{cfg['name']}: final_submission_split did not cover every curated row -- stopping.")
    split_map = {"final_train": "train", "final_val": "val"}
    population["chemprop_split"] = population["final_submission_split"].map(split_map)

    run_dir = CHEMPROP_RUNS_DIR / cfg["name"]
    run_dir.mkdir(parents=True, exist_ok=True)
    train_csv = run_dir / "train_input.csv"
    cols = ["canonical_smiles", *REGRESSION_ENDPOINTS, "chemprop_split"]
    population[cols].to_csv(train_csv, index=False)
    logger.info(
        f"wrote {train_csv} ({len(population)} rows, chemprop_split counts: "
        f"{population['chemprop_split'].value_counts().to_dict()})"
    )

    run_chemprop_train(
        train_csv, REGRESSION_ENDPOINTS, run_dir, logger,
        CHEMPROP_ARCHITECTURE_ARGS, CHEMPROP_EPOCHS, CHEMPROP_PATIENCE, SEED,
    )

    raw_pred_csv = run_dir / "raw_predictions.csv"
    run_chemprop_predict(run_dir / "model_0", predict_csv, raw_pred_csv, logger)
    expected_names = set(test_df["Molecule_Name"])
    verify_predictions(raw_pred_csv, expected_names, REGRESSION_ENDPOINTS, logger)

    BLIND_PRED_DIR.mkdir(parents=True, exist_ok=True)
    raw_pred = pd.read_csv(raw_pred_csv)
    raw_pred.to_csv(out_path, index=False)
    logger.info(f"wrote {out_path} ({len(raw_pred)} rows) in {time.time() - t0:.1f}s")


def write_config_manifest(logger) -> None:
    rows = []
    for isoform, members in ISOFORM_COMBO_MEMBERS.items():
        for combo_member_label, blind_pred_file, pred_column in members:
            rows.append({
                "isoform": isoform, "combo_member_label": combo_member_label,
                "blind_pred_file": blind_pred_file, "pred_column": pred_column,
            })
    manifest = pd.DataFrame(rows)
    out_path = OUT / "config_manifest.csv"
    manifest.to_csv(out_path, index=False)
    logger.info(f"wrote {out_path} ({len(manifest)} rows)")


def main() -> None:
    logger = setup_logging(LOGS / "10_final_retrain_predict.log", "10_final_retrain_predict")
    script_start = time.time()
    logger.info("=" * 70)
    logger.info("starting scripts/10_final_retrain_predict.py")
    logger.info(
        "chemeleon__lightgbm is NOT retrained -- confirmed not a step3_summary.csv "
        "winning-combo member for any isoform (candidate-pool member for CYP2D6 only, "
        "excluded from the winning combo there). See module docstring."
    )

    OUT.mkdir(parents=True, exist_ok=True)
    BLIND_PRED_DIR.mkdir(parents=True, exist_ok=True)
    CHEMPROP_RUNS_DIR.mkdir(parents=True, exist_ok=True)

    curated = pd.read_csv(CURATED_PATH)
    logger.info(f"loaded {CURATED_PATH.name}: {curated.shape}")
    test_df = pd.read_csv(TEST_PATH)
    logger.info(f"loaded {TEST_PATH.name}: {test_df.shape}")
    if len(test_df) != 750:
        raise ValueError(f"expected 750 blinded test compounds, got {len(test_df)} -- stopping.")

    flagged = pd.read_csv(FLAGGED_PATH)
    excluded_inchikeys = set(flagged["inchikey"])
    logger.info(f"loaded {FLAGGED_PATH.name}: {len(flagged)} CYP2D6-flagged compounds (06's confirmed residual criterion)")

    split_out_path = OUT / "final_train_val_split.csv"
    if split_out_path.exists():
        final_split = pd.read_csv(split_out_path)
        logger.info(f"reusing existing {split_out_path} ({len(final_split)} rows) -- deterministic function of (curated, seed, val_fraction), safe to reuse across configs/reruns")
    else:
        final_split = assign_final_submission_split(curated, val_fraction=VAL_FRACTION, seed=SEED)
        final_split.to_csv(split_out_path, index=False)
        logger.info(
            f"wrote {split_out_path} ({len(final_split)} rows) via "
            f"assign_final_submission_split(val_fraction={VAL_FRACTION}, seed={SEED}); "
            f"split counts: {final_split['final_submission_split'].value_counts().to_dict()}"
        )

    predict_csv = OUT / "predict_input.csv"
    test_df[["Molecule_Name", "inchikey", "canonical_smiles"]].to_csv(predict_csv, index=False)
    logger.info(f"wrote {predict_csv} ({len(test_df)} blind test compounds, shared across both chemprop configs)")

    logger.info("=" * 70)
    logger.info("stage: tabular configs (RF/LightGBM, fast)")
    run_tabular_configs(curated, test_df, excluded_inchikeys, final_split, logger)

    logger.info("=" * 70)
    logger.info("stage: chemprop configs (slow -- expect hours per config)")
    for cfg in CHEMPROP_CONFIGS:
        run_chemprop_config(cfg, curated, test_df, excluded_inchikeys, final_split, predict_csv, logger)

    write_config_manifest(logger)

    total_elapsed = time.time() - script_start
    logger.info("=" * 70)
    logger.info(f"total script wall time: {total_elapsed:.1f}s ({total_elapsed / 60:.1f} min)")
    logger.info("done.")


if __name__ == "__main__":
    main()
