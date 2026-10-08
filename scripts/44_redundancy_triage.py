#!/usr/bin/env python
"""Characterise the QM columns against size and against the free 2D blocks. No labels.

Reproduces every number in the plan's Part 2.2/2.3 and revision 2's R2.1 into two versioned
tables, so notebook 44 plots a file rather than recomputing 17 minutes of fits inline:

  outputs/44_qm_molecular_block/redundancy_triage.csv
  outputs/44_qm_molecular_block/dft_vs_gfn2_increment.csv

NOTHING HERE TOUCHES A pIC50 LABEL. Every quantity is a property of the feature matrices alone,
which is what lets the pre-registration stay blind to the outcome: notebook 36 had to admit its
own criteria were "NOT a true pre-registration" because it read the relevant diagnostics before
freezing them.

Three recoverability measures are reported per column, because they disagree and the
disagreements are informative:
  1. max |Spearman| against any single free descriptor (the plan's own 0.9 bar);
  2. out-of-sample R^2 from the WHOLE free block, rank-space RidgeCV -- linear, robust,
     scale-free (an early raw-space attempt returned R^2 as low as -6,400: `Ipc` reaches
     5.06e+14 and becomes a ~75-sigma leverage point after standardisation);
  3. the same with HistGradientBoostingRegressor -- nonlinear, i.e. the strongest cheap recovery.

The honest headline number is the MAX of 2 and 3: if any cheap model reproduces it, it is cheaply
reproducible. That max is upward-biased by construction, which is stated rather than hidden --
and reproducibility of a FEATURE is not predictive redundancy, which only the increment arms can
settle.

Usage:  python scripts/44_redundancy_triage.py
"""
from __future__ import annotations

import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, rankdata, spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import RidgeCV
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold, cross_val_predict

warnings.filterwarnings("ignore")
REPO_ROOT = Path(__file__).resolve().parent.parent
OUT = REPO_ROOT / "outputs" / "44_qm_molecular_block"
FEAT_DIR = OUT / "features"
QM_DESC = REPO_ROOT / "outputs" / "qm_descriptors" / "descriptors.csv"
QM_MOL = REPO_ROOT / "outputs" / "qm_parsed" / "qm_molecular_raw.csv"
RDKIT2D_PATH = REPO_ROOT / "data" / "processed" / "tabular_baseline_features_rdkit2d.csv"
ECFP4_PATH = REPO_ROOT / "data" / "processed" / "tabular_baseline_features.csv"
MORDRED_PATH = REPO_ROOT / "data" / "processed" / "tabular_mordred_pca.csv"

NINE = ["mol_wt", "hba", "hbd", "ring_count", "stereocenter_count", "formal_charge",
        "logp", "logd_proxy", "vsa_acc_proxy"]
KF = KFold(5, shuffle=True, random_state=42)
HGB = lambda: HistGradientBoostingRegressor(random_state=0, max_iter=200)  # noqa: E731
RIDGE = lambda: RidgeCV(alphas=np.logspace(-2, 4, 13))                     # noqa: E731


def rank_z(a):
    r = rankdata(a) / len(a)
    return (r - r.mean()) / r.std()


def main() -> int:
    d = pd.read_csv(QM_DESC)
    m = pd.read_csv(QM_MOL)
    r2d = pd.read_csv(RDKIT2D_PATH)
    ecfp = pd.read_csv(ECFP4_PATH, usecols=["Molecule_Name"] + NINE)

    Q = d.merge(m[["Molecule_Name", "homo1_ev", "homo2_ev", "lumo1_ev",
                   "gfn2_relaxation_strain_eh", "rot_const_a_cm1", "rot_const_b_cm1",
                   "rot_const_c_cm1"]], on="Molecule_Name", validate="1:1")
    # the shape block, straight from the frozen feature file rather than recomputed
    shp = pd.read_csv(FEAT_DIR / "QM_3D.csv")
    Q = Q.merge(shp[["Molecule_Name", "npr1", "npr2", "radius_of_gyration", "I1"]],
                on="Molecule_Name", validate="1:1")
    Q["chi_electronegativity"] = -(Q.dft_homo_ev + Q.dft_lumo_ev) / 2
    Q["omega_electrophilicity"] = Q.chi_electronegativity ** 2 / (Q.dft_lumo_ev - Q.dft_homo_ev)
    Q["dft_minus_gfn2_gap"] = Q.dft_gap_ev - Q.gfn2_gap_ev
    Q["energy_per_atom"] = Q.dft_energy_eh / Q.total_atoms_with_h
    Q["dispersion_per_atom"] = Q.dft_dispersion_eh / Q.total_atoms_with_h
    Q["hirshfeld_range"] = Q.hirshfeld_max - Q.hirshfeld_min

    TARGETS = [
        # size / composition proxies
        "heavy_atoms", "total_atoms_with_h", "n_basis", "gfn2_energy_eh", "dft_dispersion_eh",
        # reproducible from the free block
        "formal_charge", "dft_energy_eh", "energy_per_atom", "dispersion_per_atom",
        "hirshfeld_max", "hirshfeld_abs_mean", "hirshfeld_range", "hirshfeld_min",
        # genuinely new: electronic
        "dft_homo_ev", "dft_lumo_ev", "dft_gap_ev", "homo1_ev", "homo2_ev", "lumo1_ev",
        "gfn2_homo_ev", "gfn2_lumo_ev", "gfn2_gap_ev", "dft_dipole_debye",
        "chi_electronegativity", "omega_electrophilicity", "dft_minus_gfn2_gap",
        "gfn2_relaxation_strain_eh",
        # genuinely new: 3D shape
        "npr1", "npr2", "radius_of_gyration", "I1",
        # the reparameterisations, reported so the identity is visible in the table
        "rot_const_a_cm1", "rot_const_b_cm1", "rot_const_c_cm1",
    ]
    GROUP = {}
    for c in TARGETS:
        GROUP[c] = ("b_size_proxy" if c in {"heavy_atoms", "total_atoms_with_h", "n_basis",
                                            "gfn2_energy_eh", "dft_dispersion_eh"}
                    else "c_reproducible" if c in {"formal_charge", "dft_energy_eh",
                                                   "energy_per_atom", "hirshfeld_max",
                                                   "hirshfeld_abs_mean", "hirshfeld_range",
                                                   "hirshfeld_min"}
                    else "reparameterisation" if c.startswith("rot_const_")
                    else "a_genuine")

    rc = [c for c in r2d.columns if c not in ("Molecule_Name", "inchikey", "split")]
    base = Q[["Molecule_Name", "set"] + TARGETS].merge(
        r2d[["Molecule_Name"] + rc], on="Molecule_Name", validate="1:1").merge(
        ecfp.rename(columns={c: "n9_" + c for c in NINE}), on="Molecule_Name", validate="1:1")
    cheap = rc + ["n9_" + c for c in NINE]
    print(f"free block: rdkit2d {len(rc)} + ecfp4_narrow non-fp {len(NINE)} = {len(cheap)}"
          f"   n = {len(base)}")
    Xr = np.column_stack([rank_z(base[c].to_numpy(float)) for c in cheap])
    Xraw = base[cheap].to_numpy(float)

    mp = pd.read_csv(MORDRED_PATH)
    pcs = [c for c in mp.columns if c.startswith("mordred_pca_")]
    mbase = Q[["Molecule_Name"] + TARGETS].merge(mp[["Molecule_Name"] + pcs],
                                                on="Molecule_Name", validate="1:1")
    print(f"mordred_pca: {len(pcs)} PCs, n = {len(mbase)} (TRAIN ONLY -- that file has no "
          f"blind rows, which is why it cannot be the control)")
    Xm = np.column_stack([rank_z(mbase[c].to_numpy(float)) for c in pcs])

    rows, t0 = [], time.time()
    for t in TARGETS:
        y = base[t].to_numpy(float)
        yz = rank_z(y)
        ridge_r2 = r2_score(yz, cross_val_predict(RIDGE(), Xr, yz, cv=KF))
        hgb_r2 = r2_score(y, cross_val_predict(HGB(), Xraw, y, cv=KF))
        best, bc = 0.0, None
        for c in cheap:
            v = base[c].to_numpy(float)
            if np.nanstd(v) == 0:
                continue
            rr = abs(spearmanr(y, v).statistic)
            if rr > best:
                best, bc = rr, c
        ym = rank_z(mbase[t].to_numpy(float))
        mord_r2 = r2_score(ym, cross_val_predict(RIDGE(), Xm, ym, cv=KF))
        rows.append({
            "column": t, "group": GROUP[t],
            "rho_heavy_atoms": spearmanr(y, base["heavy_atoms"]).statistic,
            "rho_molwt": spearmanr(y, base["MolWt"]).statistic,
            "max_abs_size_corr": max(abs(spearmanr(y, base["heavy_atoms"]).statistic),
                                     abs(pearsonr(y, base["heavy_atoms"]).statistic),
                                     abs(spearmanr(y, base["MolWt"]).statistic),
                                     abs(pearsonr(y, base["MolWt"]).statistic)),
            "max_abs_rho_single_free": best, "best_single_free_descriptor": bc,
            "cv_R2_rank_ridge_free": ridge_r2, "cv_R2_histgb_free": hgb_r2,
            "cv_R2_max_free": max(ridge_r2, hgb_r2),
            "cv_R2_rank_ridge_mordred_pca": mord_r2,
        })
        print(f"  {t:28s} size={rows[-1]['max_abs_size_corr']:.3f}  "
              f"ridge={ridge_r2:+.4f}  histgb={hgb_r2:+.4f}  max={rows[-1]['cv_R2_max_free']:+.4f}")
    tri = pd.DataFrame(rows)
    tri.to_csv(OUT / "redundancy_triage.csv", index=False)
    print(f"\nwrote redundancy_triage.csv ({len(tri)} rows, {time.time()-t0:.0f}s)")

    # ----------------------------------------------------------- revision 2 R2.1
    print("\n=== DFT vs GFN2: is the semi-empirical column redundant with the free block? ===")
    print("Three numbers per quantity, ALL out-of-sample R^2, same model, same folds, so the")
    print("two kinds of estimate are never set against each other again.")
    inc = []
    for q, gq in [("dft_homo_ev", "gfn2_homo_ev"), ("dft_lumo_ev", "gfn2_lumo_ev"),
                  ("dft_gap_ev", "gfn2_gap_ev")]:
        y = base[q].to_numpy(float)
        g = base[gq].to_numpy(float).reshape(-1, 1)
        a = r2_score(y, cross_val_predict(HGB(), Xraw, y, cv=KF))
        b = r2_score(y, cross_val_predict(HGB(), g, y, cv=KF))
        c = r2_score(y, cross_val_predict(HGB(), np.hstack([Xraw, g]), y, cv=KF))
        tr = base["set"] == "train"
        inc.append({"quantity": q.replace("dft_", "").replace("_ev", ""),
                    "free_226_oos_r2": a, "gfn2_alone_oos_r2": b, "free_plus_gfn2_227_oos_r2": c,
                    "increment_c_minus_a": c - a,
                    "univariate_pearson2_in_sample_NOT_COMPARABLE": pearsonr(y, g.ravel())[0] ** 2,
                    "spearman_dft_vs_gfn2_all": spearmanr(y, g.ravel()).statistic,
                    "spearman_train": spearmanr(y[tr.to_numpy()], g.ravel()[tr.to_numpy()]).statistic,
                    "spearman_blind": spearmanr(y[~tr.to_numpy()], g.ravel()[~tr.to_numpy()]).statistic})
        print(f"  {inc[-1]['quantity']:5s} free={a:.4f}  gfn2_alone={b:.4f}  "
              f"free+gfn2={c:.4f}  increment={c-a:+.4f}   "
              f"(train rho {inc[-1]['spearman_train']:.4f} / blind {inc[-1]['spearman_blind']:.4f})")
    pd.DataFrame(inc).to_csv(OUT / "dft_vs_gfn2_increment.csv", index=False)
    print("\nwrote dft_vs_gfn2_increment.csv")
    print("\nDONE -- no label read, no model trained for prediction, nothing in data/ written.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
