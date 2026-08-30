"""Compute Mordred 2D descriptors + PCA-reduced features for the tabular arm of the
upcoming 5x5 CV, replacing full RDKit2D there -- per Mauricio's suggestion.

COMPATIBILITY (checked separately, before this script was written, per the explicit
request not to silently work around problems): `mordred-community` (PyPI package
`mordredcommunity==2.0.7`, already pinned in `environment.yml`; import name is still
`mordred` -- NOT the original, unmaintained `mordred` PyPI package) installs cleanly and
runs without any exceptions, Python warnings, or RDKit stderr output against this
project's pinned RDKit 2026.3.3, confirmed on an 8-compound sample of
`train_inhibition_curated.csv` before this script was written.

NaN HANDLING -- a user decision, not a silent pick: the full run produces 1613 2D
descriptors, of which 82 are 100%-NaN (zero information -- exotic elements/bonding
patterns absent from this dataset) and 245 partially NaN (9.18% of all cells overall).
Presented to the user with the full per-descriptor/per-compound breakdown; user chose:
drop the 82 always-fail descriptors, median-impute the remaining 245 partial-NaN ones,
medians computed from these training compounds only (this script never touches any test
data). Kept: 1531 of 1613 descriptors.

STANDARDIZATION: applied (mean 0, std 1 per descriptor, fit on train only) before PCA --
not explicitly requested in the prompt, but flagged here rather than silently baked in:
Mordred descriptors span wildly different natural scales (e.g. molecular weight ~100s vs.
some topological indices ~0.01), and PCA is not scale-invariant -- fitting PCA on raw
unscaled descriptors would make the explained-variance curve (and the components
themselves) dominated by whichever descriptors happen to have the largest raw magnitude,
not the largest genuine variance in shape space. This is standard practice for PCA, not a
tuned hyperparameter, but noted explicitly since it changes the numbers the user is being
shown.

SINGLE FIT, MAXIMALLY REUSABLE: imputer + scaler + PCA are fit ONCE, together, as one
sklearn Pipeline, on the training compounds only (`train_inhibition_curated.csv`, 4905
rows) -- never on the 750-compound blind test set, which this script does not open at
all (leakage avoidance; the blind-test transform is explicit future work, not part of
this task). PCA itself is fit with the maximum number of components (1531 = min(n_samples,
n_features)) rather than pre-truncated to a guessed count, so: (a) the full
explained-variance-vs-component-count curve can be shown for the user's own call, and
(b) the ONE saved pipeline object can produce any component-count slice later (including
for the blind test set) without ever being refit -- "the same fitted transform, not a
separately-fit one," per the user's own requirement. The saved training feature CSV uses
a 95%-cumulative-variance default slice, explicitly flagged as a suggestion (see log
output / explained_variance_curve.csv) -- not a locked-in decision.

OUTPUTS:
    logs/generate_mordred_pca_features.log                    -- progress log
    outputs/mordred_pca/descriptor_nan_report.csv              -- per-descriptor NaN
                                                                    count/rate (all 1613,
                                                                    pre-drop)
    outputs/mordred_pca/explained_variance_curve.csv           -- component, individual
                                                                    and cumulative
                                                                    explained variance
                                                                    ratio, all 1531
                                                                    components
    outputs/mordred_pca/explained_variance_curve.png           -- plot of the above
    models/mordred_pca_pipeline.joblib                          -- fitted
                                                                    {imputer, scaler, pca,
                                                                    dropped_descriptors,
                                                                    descriptor_columns}
                                                                    bundle -- reuse this
                                                                    (never refit) to
                                                                    transform the blind
                                                                    test set later
    data/processed/tabular_mordred_pca.csv                     -- Molecule_Name,
                                                                    inchikey, split="train",
                                                                    mordred_pca_0000..N
                                                                    (N = suggested 95%-
                                                                    variance default)
"""

import sys
import time
from pathlib import Path

import joblib
import matplotlib

matplotlib.use("Agg")  # headless -- this script runs non-interactively
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))  # so `from src... import ...` works regardless of cwd

from src.chemprop_screen import setup_logging
from src.features import mordred_2d_descriptors

PROCESSED = REPO_ROOT / "data" / "processed"
LOGS = REPO_ROOT / "logs"
MODELS = REPO_ROOT / "models"
RUN_OUT = REPO_ROOT / "outputs" / "mordred_pca"

LOG_PATH = LOGS / "generate_mordred_pca_features.log"
FEATURES_OUT_PATH = PROCESSED / "tabular_mordred_pca.csv"
PIPELINE_OUT_PATH = MODELS / "mordred_pca_pipeline.joblib"

VARIANCE_THRESHOLD_SUGGESTED = 0.95  # user's own example default -- flagged, not final


def main():
    RUN_OUT.mkdir(parents=True, exist_ok=True)
    MODELS.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(LOG_PATH, "generate_mordred_pca_features")
    script_start = time.time()
    logger.info("=" * 70)
    logger.info("starting generate_mordred_pca_features.py")

    # ---- 1. load training compounds ONLY -- never test_blinded_curated.csv ---------
    curated_path = PROCESSED / "train_inhibition_curated.csv"
    curated = pd.read_csv(curated_path)
    logger.info(f"loaded {curated_path.name}: {curated.shape}")
    logger.info(f"rows before/after (no filtering, this script only computes descriptors): "
                f"{len(curated)}/{len(curated)}")

    # ---- 2. compute full Mordred 2D descriptor set ----------------------------------
    t0 = time.time()
    raw = mordred_2d_descriptors(curated["canonical_smiles"].tolist())
    logger.info(f"mordred_2d_descriptors: {raw.shape} in {time.time() - t0:.1f}s "
                f"(nproc={__import__('os').cpu_count()})")
    assert len(raw) == len(curated)

    # ---- 3. NaN diagnostics (report BEFORE dropping/imputing anything) --------------
    nan_per_desc = raw.isna().sum(axis=0)
    nan_rate_per_desc = nan_per_desc / len(raw)
    nan_report = pd.DataFrame({
        "descriptor": nan_rate_per_desc.index,
        "nan_count": nan_per_desc.values,
        "nan_rate": nan_rate_per_desc.values,
    }).sort_values("nan_rate", ascending=False)
    nan_report_path = RUN_OUT / "descriptor_nan_report.csv"
    nan_report.to_csv(nan_report_path, index=False)

    always_fail = nan_rate_per_desc[nan_rate_per_desc == 1.0].index.tolist()
    clean = (nan_rate_per_desc == 0.0).sum()
    partial = raw.shape[1] - len(always_fail) - clean
    nan_per_compound_rate = raw.isna().sum(axis=1) / raw.shape[1]
    logger.info(f"descriptor NaN summary: {clean} always-clean, {len(always_fail)} always-fail "
                f"(100% NaN), {partial} partial-fail -- full breakdown at {nan_report_path}")
    logger.info(f"overall cell-level NaN rate: {raw.isna().values.mean() * 100:.2f}%")
    logger.info(f"per-compound NaN rate: mean={nan_per_compound_rate.mean() * 100:.2f}%, "
                f"max={nan_per_compound_rate.max() * 100:.2f}%, "
                f">10% NaN: {(nan_per_compound_rate > 0.10).sum()} compounds")

    # ---- 4. drop always-fail, keep partial-fail for imputation (user's own decision) --
    kept_cols = [c for c in raw.columns if c not in always_fail]
    features = raw[kept_cols].copy()
    logger.info(f"dropped {len(always_fail)} always-fail (100% NaN) descriptors -- "
                f"kept {len(kept_cols)} of {raw.shape[1]}")
    logger.info(f"rows before/after descriptor drop: {len(features)}/{len(features)} "
                "(row count unaffected -- only columns dropped)")

    # ---- 5. impute (median, train-only) -> standardize -> PCA (full components) -----
    n_components = min(features.shape[0], features.shape[1])
    pipeline = Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("pca", PCA(n_components=n_components, svd_solver="full", random_state=42)),
    ])
    t0 = time.time()
    transformed = pipeline.fit_transform(features.to_numpy())
    logger.info(f"fit_transform (median-impute -> standardize -> PCA, "
                f"n_components={n_components}): {time.time() - t0:.1f}s")
    logger.info("seed: 42 (PCA random_state -- unused by svd_solver='full', logged for the record)")

    pca = pipeline.named_steps["pca"]
    cum_var = np.cumsum(pca.explained_variance_ratio_)
    curve_df = pd.DataFrame({
        "n_components": np.arange(1, n_components + 1),
        "explained_variance_ratio": pca.explained_variance_ratio_,
        "cumulative_explained_variance": cum_var,
    })
    curve_path = RUN_OUT / "explained_variance_curve.csv"
    curve_df.to_csv(curve_path, index=False)

    thresholds = {}
    for t in (0.90, 0.95, 0.99):
        n_needed = int(np.searchsorted(cum_var, t) + 1)
        thresholds[t] = n_needed
        logger.info(f"components needed for {t * 100:.0f}% cumulative variance: {n_needed}")

    n_suggested = thresholds[VARIANCE_THRESHOLD_SUGGESTED]
    logger.info(f"SUGGESTED default (not a locked-in decision): {n_suggested} components "
                f"covering {VARIANCE_THRESHOLD_SUGGESTED * 100:.0f}% cumulative variance")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(curve_df["n_components"], curve_df["cumulative_explained_variance"])
    for t, n_needed in thresholds.items():
        ax.axhline(t, color="gray", linestyle="--", linewidth=0.8)
        ax.axvline(n_needed, color="gray", linestyle="--", linewidth=0.8)
        ax.annotate(f"{t * 100:.0f}% @ n={n_needed}", (n_needed, t),
                    textcoords="offset points", xytext=(5, -12), fontsize=8)
    ax.set_xlabel("number of PCA components")
    ax.set_ylabel("cumulative explained variance")
    ax.set_title("Mordred 2D descriptors (standardized) -- PCA explained variance, "
                  "fit on 4905 training compounds only")
    fig.tight_layout()
    plot_path = RUN_OUT / "explained_variance_curve.png"
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    logger.info(f"wrote {plot_path}")

    # ---- 6. save the fitted pipeline (single fit, reusable for the blind test set) --
    bundle = {
        "pipeline": pipeline,
        "dropped_descriptors": always_fail,
        "descriptor_columns": kept_cols,
        "n_components_fit": n_components,
        "suggested_n_components": n_suggested,
    }
    joblib.dump(bundle, PIPELINE_OUT_PATH)
    logger.info(f"wrote fitted pipeline bundle to {PIPELINE_OUT_PATH} "
                f"(impute+scale+pca, {n_components} components, all fit on these 4905 "
                "training compounds only)")

    # ---- 7. save the suggested-default PCA-reduced training feature matrix ----------
    pc_cols = [f"mordred_pca_{i:04d}" for i in range(n_suggested)]
    out_df = pd.DataFrame(transformed[:, :n_suggested], columns=pc_cols)
    out_df.insert(0, "split", "train")
    out_df.insert(0, "inchikey", curated["inchikey"].values)
    out_df.insert(0, "Molecule_Name", curated["Molecule_Name"].values)
    out_df.to_csv(FEATURES_OUT_PATH, index=False)
    logger.info(f"wrote {FEATURES_OUT_PATH} ({out_df.shape}) -- {n_suggested} PCA components "
                f"({VARIANCE_THRESHOLD_SUGGESTED * 100:.0f}% variance, suggested default)")
    logger.info(f"rows before/after: {len(curated)}/{len(out_df)}")

    logger.info(f"total runtime: {time.time() - script_start:.1f}s")
    logger.info("done.")


if __name__ == "__main__":
    main()
