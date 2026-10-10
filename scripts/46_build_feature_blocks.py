"""Feature blocks and the frozen pre-registration for notebook 46's QM arm.

Builds the six blocks the arm compares and writes `prereg.json` BEFORE anything is trained.
Nothing here reads a fold file or a prediction, so the criteria cannot be tuned to a result.

THE SIX BLOCKS. One variable across the family: the feature block. FREE is notebook 45's own
axis A -- rdkit2d's 217 plus ecfp4_narrow's 9 non-fingerprint columns -- because that is the
block the four descriptors were screened for redundancy against, so the increment this arm
measures is the increment over the thing the screen controlled for.

    FREE              226  the baseline
    FREE_QM13         239  + notebook 44's 13 genuinely-new molecule-level QM columns
    FREE_ATOM4        230  + notebook 45's 4 full-coverage strongly-resistant atom descriptors
    FREE_QM13_ATOM4   243  + both: the global-plus-local combination, never tested here
    FREE_ATOM4_HAL    236  + the 4, plus Cl/F values with presence flags
    FREE_CONTROL      227  + hirshfeld_range, notebook 45's pre-registered REJECTED control

NOT REUSED FROM NOTEBOOK 44, and the reason is worth recording. Notebook 44 ran the same frozen
partition, so its arms looked reusable at zero cost -- but its modelling baseline `RDKIT2D` is
rdkit2d's 217 columns alone (no ecfp4), and its `RDKIT2D_QM` adds all 18 QM columns rather than
the 13. Neither matches a block here, so reuse would have meant adopting a baseline that is not
the one notebook 45 screened against. All six are trained fresh instead; at roughly 0.2 min per
arm per 25 folds that costs minutes.

MISSINGNESS. Only FREE_ATOM4_HAL has any. Cl is defined on 669 of 4,905 training compounds and
F on 1,372, so each element-conditioned value ships with a presence flag and the value column
carries NaN where the element is absent. The fill is NOT applied here: it is fitted per fold on
that fold's training portion by the runner, so no held-out value can enter it. Writing a filled
column to this file would bake a whole-population statistic into every fold.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT = REPO_ROOT / "outputs" / "46_qm_arm"
FEAT_DIR = OUT / "features"

CURATED = REPO_ROOT / "data" / "processed" / "train_inhibition_curated.csv"
RDKIT2D = REPO_ROOT / "data" / "processed" / "tabular_baseline_features_rdkit2d.csv"
ECFP4 = REPO_ROOT / "data" / "processed" / "tabular_baseline_features.csv"
QM_FULL = REPO_ROOT / "outputs" / "44_qm_molecular_block" / "features" / "QM_FULL.csv"
AGGREGATES = REPO_ROOT / "outputs" / "qm_parsed" / "qm_atom_aggregates.csv"
NB45_VERDICT = REPO_ROOT / "outputs" / "45_atom_level_screen" / "verdict_table.csv"
NB37_OUT = REPO_ROOT / "outputs" / "37_butina_5x5_complete"
BUTINA = REPO_ROOT / "data" / "folds" / "cv_folds_butina_5x5.csv"
NB34_FOLD_SHA256 = "52f39b1a1cea9200b58dce0a061d4a21ce94d9d58ec77067cf6508b0f4b3e87e"

NINE = ["mol_wt", "hba", "hbd", "ring_count", "stereocenter_count", "formal_charge",
        "logp", "logd_proxy", "vsa_acc_proxy"]
# QM_FULL's 5 algebraic derivatives, excluded to leave the 13 genuinely-new columns.
DERIVED5 = ["dft_gap_ev", "gfn2_gap_ev", "chi_electronegativity", "omega_electrophilicity",
            "dft_minus_gfn2_gap"]
FOUR = ["hirshfeld_C_min", "hirshfeld_H_min", "bond_order_min", "mayer_va_min"]
HALOGEN = ["hirshfeld_Cl_min", "hirshfeld_Cl_max", "hirshfeld_F_min", "hirshfeld_F_max"]
CONTROL = "hirshfeld_range"
ISOFORMS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
METRICS = ["ST-RAE", "Spearman_R"]
BAR_POOL = ["34_butina_5x5", "35_butina_repeat3", "37_butina_5x5_complete"]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def nb37_bars(metric: str) -> dict:
    """Notebook 37's bar construction, reproduced exactly: 2 x the mean standard error of the
    25-fold paired difference against `plain`, pooled across its TEN arms. Never 2 x the SE of
    the comparison under test, which is the significance threshold restated.

    The pool MUST span notebooks 34+35+37. Notebook 37's own run_scores.csv holds only its
    repeats 3-4, so reading it alone gives 10 folds per comparison instead of 25 and inflates
    every SE by about sqrt(25/10) -- the bug notebook 45 caught by cross-checking against
    notebook 44's published Spearman figures.
    """
    rs = (pd.concat([pd.read_csv(REPO_ROOT / "outputs" / n / "run_scores.csv") for n in BAR_POOL],
                    ignore_index=True)
          .drop_duplicates(["repeat", "fold", "arm", "Endpoint"]))
    arms = [a for a in rs["arm"].unique() if a != "plain"]
    out = {}
    for iso in ISOFORMS:
        ses = []
        for a in arms:
            pa = rs[(rs.arm == a) & (rs.Endpoint.str.contains(iso))].set_index(["repeat", "fold"])[metric]
            pb = rs[(rs.arm == "plain") & (rs.Endpoint.str.contains(iso))].set_index(["repeat", "fold"])[metric]
            idx = pa.index.intersection(pb.index)
            if len(idx) >= 10:
                d = (pa.loc[idx] - pb.loc[idx]).to_numpy()
                ses.append(d.std(ddof=1) / np.sqrt(len(d)))
        out[iso] = round(2 * float(np.mean(ses)), 4)
    return out


def main() -> int:
    FEAT_DIR.mkdir(parents=True, exist_ok=True)
    assert sha256(BUTINA) == NB34_FOLD_SHA256, "frozen partition is not the file 34/35/37 ran on"

    curated = pd.read_csv(CURATED)
    assert len(curated) == 4905
    keys = curated[["Molecule_Name", "inchikey"]].copy()

    r2d = pd.read_csv(RDKIT2D)
    rc = [c for c in r2d.columns if c not in ("Molecule_Name", "inchikey", "split")]
    ecfp = pd.read_csv(ECFP4, usecols=["Molecule_Name"] + NINE)
    free = (keys.merge(r2d[["Molecule_Name"] + rc], on="Molecule_Name", validate="1:1")
                .merge(ecfp.rename(columns={c: "n9_" + c for c in NINE}),
                       on="Molecule_Name", validate="1:1"))
    FREE_COLS = rc + ["n9_" + c for c in NINE]
    assert len(FREE_COLS) == 226, len(FREE_COLS)

    qm = pd.read_csv(QM_FULL)
    THIRTEEN = [c for c in qm.columns
                if c not in ("Molecule_Name", "inchikey", "set", "split", *DERIVED5)]
    assert len(THIRTEEN) == 13, (len(THIRTEEN), THIRTEEN)

    agg = pd.read_csv(AGGREGATES)
    verdict = pd.read_csv(NB45_VERDICT)
    # read the four back from notebook 45's own table rather than trusting the literal above
    four_check = sorted(verdict[(verdict.tier == "strong") & (verdict.pct_blind > 99)].descriptor)
    assert four_check == sorted(FOUR), (four_check, FOUR)
    assert verdict.loc[verdict.descriptor == CONTROL, "resists_both"].iloc[0] == False, \
        "the control must be a descriptor notebook 45's screen REJECTED"

    base = free.merge(qm[["Molecule_Name"] + THIRTEEN], on="Molecule_Name", validate="1:1") \
               .merge(agg[["Molecule_Name"] + FOUR + HALOGEN + [CONTROL]],
                      on="Molecule_Name", validate="1:1")
    assert len(base) == 4905

    # presence flags, from the value columns' own missingness
    for el in ["Cl", "F"]:
        base[f"has_{el}"] = base[f"hirshfeld_{el}_min"].notna().astype(float)
    HAL_COLS = HALOGEN + ["has_Cl", "has_F"]

    blocks = {
        "FREE": FREE_COLS,
        "FREE_QM13": FREE_COLS + THIRTEEN,
        "FREE_ATOM4": FREE_COLS + FOUR,
        "FREE_QM13_ATOM4": FREE_COLS + THIRTEEN + FOUR,
        "FREE_ATOM4_HAL": FREE_COLS + FOUR + HAL_COLS,
        "FREE_CONTROL": FREE_COLS + [CONTROL],
    }
    MAY_HAVE_NAN = {"FREE_ATOM4_HAL"}

    manifest = {}
    for name, cols in blocks.items():
        df = base[["Molecule_Name", "inchikey"] + cols]
        X = df[cols].to_numpy(dtype=np.float64)
        nan_cols = [c for c in cols if df[c].isna().any()]
        if name in MAY_HAVE_NAN:
            assert set(nan_cols) <= set(HALOGEN), (name, nan_cols)
        else:
            assert not nan_cols, f"{name}: unexpected NaN in {nan_cols}"
            assert np.isfinite(X).all(), f"{name}: non-finite value"
        const = [c for c in cols if df[c].nunique(dropna=False) <= 1]
        path = FEAT_DIR / f"{name}.csv"
        df.to_csv(path, index=False)
        manifest[name] = {"n_cols": len(cols), "n_rows": len(df), "columns": cols,
                          "nan_columns": nan_cols, "n_constant": len(const),
                          "sha256": sha256(path)}
        print(f"  {name:18s} {len(cols):4d} cols  {len(df)} rows  "
              f"nan_cols={len(nan_cols)}  constant={len(const)}")

    prereg = {
        "written_at_utc": datetime.now(timezone.utc).isoformat(),
        "written_by": "scripts/46_build_feature_blocks.py",
        "status": "FROZEN BEFORE ANY TRAINING. No fold file and no prediction is read here.",
        "question": ("Do notebook 45's screened QM descriptors improve pIC50 prediction over the "
                     "free 2D block they were screened against?"),
        "partition": {"path": "data/folds/cv_folds_butina_5x5.csv", "sha256": NB34_FOLD_SHA256,
                      "n_folds": 25, "identical_across_every_method": True},
        "blocks": manifest,
        "nan_fill_policy": {
            "blocks": sorted(MAY_HAVE_NAN),
            "rule": ("median of the PRESENT subset of that fold's training portion "
                     "(screen_inner_train + screen_inner_val), applied to both the training and "
                     "held-out portions of that fold"),
            "why_not_baked_into_the_csv": ("a whole-population median would put held-out "
                                           "information into every fold's fit"),
            "presence_flags": ["has_Cl", "has_F"],
            "coverage_train": {"Cl": int(base.has_Cl.sum()), "F": int(base.has_F.sum())},
        },
        "statistical_route": {
            "name": "Ash et al. 2025 parametric default, main text Section 3.2.2",
            "steps": ["repeated-measures ANOVA per metric",
                      "Bonferroni the ANOVA alpha by the number of metrics",
                      "Tukey HSD post-hoc only where the ANOVA clears",
                      "compact letter display from the Tukey pairwise matrix"],
            "alpha": 0.05, "n_metrics": len(METRICS), "alpha_adjusted": 0.05 / len(METRICS),
            "n_methods": len(blocks),
            "expected_posthoc": "tukey",
            "why_not_BH": ("BH is Ash's exception above 10 methods and its justification lives in "
                           "SI Section A.2, which is not available in this project. The family is "
                           "held at 6 so the verifiable parametric route applies, and so Tukey's "
                           "confidence intervals -- which the materiality bars need -- exist."),
            "branch_rule": ("measured per (isoform, metric): parametric iff Levene p > 0.05 AND "
                            "max/min variance ratio <= 9. If a cell routes nonparametric, "
                            "Conover-Friedman + Holm gives no CI, and the paired-t CI is reported "
                            "beside it as an effect-size estimate, not as the test."),
            "no_correction_stacked_on_tukey": ("Tukey controls FWER within each (isoform, metric) "
                                               "family by construction"),
        },
        "co_primary_metrics": METRICS,
        "materiality_bars": {
            "ST-RAE": nb37_bars("ST-RAE"), "Spearman_R": nb37_bars("Spearman_R"),
            "construction": ("2 x the mean paired SE against `plain`, pooled across notebook 37's "
                             "ten arms, from the 34+35+37 run_scores pool"),
        },
        "tier1": {
            "purpose": "choose M*, no p-values presented as evidence",
            "partition": "data/folds/cv_folds.csv, repeat_0, fold 0 -- the RANDOM partition, so "
                         "the frozen 25-fold stays unseen until the formal tier",
            "grid": "6 blocks x 4 models x 3 seeds",
            "selection_rule": ("lowest mean macro ST-RAE on the FREE baseline block ALONE, so the "
                               "model cannot be selected for liking the QM blocks"),
        },
        "verdict_rule": {
            "MATERIAL": ("ANOVA clears alpha_adjusted AND Tukey significant vs FREE AND "
                         "|difference| >= that isoform's bar, on >=1 isoform, and the passing set "
                         "is not {CYP2D6} alone"),
            "DETECTED_NOT_MATERIAL": "significant but |difference| below the bar",
            "NOT_DIFFERENT": "ANOVA fails or Tukey does not clear",
            "expected_outcome": ("DETECTED, NOT MATERIAL. Every one of notebook 44's 10 "
                                 "BH-significant cells was under its bar, the largest at 53%."),
            "consequence_of_DETECTED_NOT_MATERIAL": ("recorded as real but immaterial; does NOT "
                                                     "justify a submission, and does NOT justify "
                                                     "the chemprop atom-feature arm"),
            "cyp2d6": "cannot satisfy the criteria alone -- notebook 21",
        },
        "conditioning_disclosure": ("the formal tier's p-values are conditional on tier 1's choice "
                                    "of M*. Mitigation: M* is chosen on the baseline block alone, "
                                    "on a different partition."),
        "not_a_submission_decision": ("four CV-to-board transfer failures are on record. No arm "
                                      "from this registration is submitted on CV evidence alone."),
    }
    (OUT / "prereg.json").write_text(json.dumps(prereg, indent=2))
    print(f"\nwrote {OUT / 'prereg.json'}")
    print(f"  ST-RAE bars   : {prereg['materiality_bars']['ST-RAE']}")
    print(f"  Spearman bars : {prereg['materiality_bars']['Spearman_R']}")
    print(f"  methods: {len(blocks)}  alpha_adjusted: {prereg['statistical_route']['alpha_adjusted']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
