#!/usr/bin/env python
"""Build the molecule-level QM feature blocks and freeze the pre-registration.

Notebook 44's companion builder. Produces one CSV per feature block under
`outputs/44_qm_molecular_block/features/`, in the same schema the 5x5 harness expects
(`Molecule_Name`, `inchikey`, `split`, then features), plus `prereg.json`.

NOTHING IS TRAINED HERE and no label is read. Every block is a function of the feature
matrices alone. `prereg.json` is written by this script specifically so it exists before
`44_run_qm_arms.py` is ever invoked.

Sources, all read-only:
  outputs/qm_descriptors/descriptors.csv      (5,655 x 35)
  outputs/qm_parsed/qm_molecular_raw.csv      (5,655 x 23)
  outputs/qm_parsed/qm_atoms.csv              (252,237 atom rows, optimised coordinates)
  data/processed/tabular_baseline_features_rdkit2d.csv   (217 free RDKit 2D descriptors)
  data/processed/train_inhibition_curated.csv / test_blinded_curated.csv  (row counts only)

`data/processed/` is never written: it is on the protected list of every recent notebook's
scope check, so writing there would make this work indistinguishable from a violation.

Usage:  python scripts/44_build_qm_feature_blocks.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.vendor.openadmet_eval.config import REGRESSION_ENDPOINTS  # noqa: E402

ISOFORMS = [e.split("_")[0] for e in REGRESSION_ENDPOINTS]

OUT = REPO_ROOT / "outputs" / "44_qm_molecular_block"
FEAT_DIR = OUT / "features"
QM_DESC = REPO_ROOT / "outputs" / "qm_descriptors" / "descriptors.csv"
QM_MOL = REPO_ROOT / "outputs" / "qm_parsed" / "qm_molecular_raw.csv"
QM_ATOMS = REPO_ROOT / "outputs" / "qm_parsed" / "qm_atoms.csv"
RDKIT2D_PATH = REPO_ROOT / "data" / "processed" / "tabular_baseline_features_rdkit2d.csv"
CURATED_PATH = REPO_ROOT / "data" / "processed" / "train_inhibition_curated.csv"
BLIND_PATH = REPO_ROOT / "data" / "processed" / "test_blinded_curated.csv"
NB37_OUT = REPO_ROOT / "outputs" / "37_butina_5x5_complete"

N_TOTAL, N_TRAIN, N_BLIND, N_ATOMS = 5655, 4905, 750, 252237

# Atomic masses for the mass-weighted inertia tensor. Only the nine elements notebook 41
# established are present; anything else raises rather than defaulting.
MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
        "P": 30.974, "S": 32.06, "Cl": 35.45, "Br": 79.904, "I": 126.904}

# ---------------------------------------------------------------- block definitions
# Group (a) electronic, 9 columns: not reproducible from the free 2D block at R^2 >= 0.85 by
# either a linear or a nonlinear cheap model (plan Part 2.3).
QM_ELEC_PRIMARY = ["dft_homo_ev", "dft_lumo_ev", "homo1_ev", "homo2_ev", "lumo1_ev",
                   "gfn2_homo_ev", "gfn2_lumo_ev", "dft_dipole_debye",
                   "gfn2_relaxation_strain_eh"]
# Group (a'), 5 columns: exact algebraic derivatives of the above. Zero extra information, a
# real reparameterisation for ridge and a mild one for trees. eta (hardness) is excluded --
# it is gap/2, a pure rescale.
QM_ELEC_DERIVED = ["dft_gap_ev", "gfn2_gap_ev", "chi_electronegativity",
                   "omega_electrophilicity", "dft_minus_gfn2_gap"]
QM_ELEC = QM_ELEC_PRIMARY + QM_ELEC_DERIVED
# The four independent shape axes. Ten candidate shape descriptors collapse to these: NPR1 ~
# normalised asphericity ~ relative shape anisotropy (|rho| 0.999-1.000), NPR2 ~ normalised
# acylindricity (0.998), Rg ~ I2 ~ I3 ~ max atom radius (0.90-0.98), I1 near-independent of Rg.
QM_3D = ["npr1", "npr2", "radius_of_gyration", "I1"]
QM_FULL = QM_ELEC + QM_3D
# The three quantities no cheap 2D model reproduces above R^2 0.32. Revision 1 item 1.
QM_NEW3 = ["dft_dipole_debye", "npr1", "npr2"]
# 3-vs-3 cost attribution: the same three quantities at each level of theory.
QM_DFT3 = ["dft_homo_ev", "dft_lumo_ev", "dft_gap_ev"]
QM_GFN23 = ["gfn2_homo_ev", "gfn2_lumo_ev", "gfn2_gap_ev"]
# The borderline column: R^2 0.885 against a 0.90 bar, assigned to (c) and carried as a
# sensitivity rather than resolved by fiat.
HMIN = ["hirshfeld_min"]

# Mechanism-matched cheap controls, chosen on chemistry and frozen here before any training.
RDKIT_MATCH18 = ["MaxPartialCharge", "MinPartialCharge", "BCUT2D_CHGHI", "BCUT2D_CHGLO",
                 "MaxEStateIndex", "MinEStateIndex", "TPSA", "FractionCSP3",
                 "NumAromaticRings", "MolMR", "LabuteASA", "NumRotatableBonds",
                 "BertzCT", "Chi1n", "Kappa2", "Kappa3", "HallKierAlpha", "PEOE_VSA6"]
# Each the top single free correlate of its DFT counterpart in the plan's own Part 2.3 table.
RDKIT_MATCH3_ELEC = ["fr_aniline", "FractionCSP3", "NumAromaticRings"]
# Measured top free correlates of QM_NEW3's own columns (revision 2 R2.2): dipole ->
# NumHeteroatoms 0.217; NPR1 -> BalabanJ 0.285; NPR2 -> BalabanJ 0.225, so the third slot
# takes PEOE_VSA8 (shape's next, distinct) rather than repeating BalabanJ.
RDKIT_MATCH3_NEW3 = ["NumHeteroatoms", "BalabanJ", "PEOE_VSA8"]

N_RAND_DRAWS = 20
RAND_BASE_SEED = 44_000

# Columns that must never appear in any block: size/composition proxies, columns reproducible
# from the free 2D block, and run provenance. The exclusion is part of the pre-registration,
# not a modelling convenience -- a tree splitting on n_basis has measured molecular size by DFT.
FORBIDDEN = {
    # (b) size / composition proxies
    "heavy_atoms", "total_atoms_with_h", "n_atoms", "n_basis", "n_basis_estimated",
    "gfn2_energy_eh", "dft_dispersion_eh", "dispersion_eh",
    # (c) reproducible from the free 2D block at R^2 >= 0.90
    "formal_charge", "dft_energy_eh", "total_energy_eh", "energy_per_atom",
    "hirshfeld_max", "hirshfeld_abs_mean", "hirshfeld_range",
    # (d) run provenance
    "status", "error", "cached", "gfn2_opt_converged", "dft_not_fully_converged",
    "dft_scf_trouble", "stage1_s", "jobA_s", "jobB_s", "gfn2_orca_wall_s",
    "dft_orca_wall_s", "total_s", "freed_bytes", "dft_scf_cycles", "scf_cycles",
    "gfn2_opt_cycles", "gfn2_final_grad_norm",
    # monotone reparameterisations of I1/I2/I3 -- zero extra information for a tree
    "rot_const_a_cm1", "rot_const_b_cm1", "rot_const_c_cm1",
    # unit duplicate
    "gfn2_relaxation_strain_kcal",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def shape_descriptors(atoms: pd.DataFrame) -> pd.DataFrame:
    """Mass-weighted inertia / gyration shape block from the GFN2-optimised coordinates.

    Mirror-invariant throughout, which is why this block is usable at all: notebook 41
    established 93.2% of centre-bearing blind compounds never state configuration, so the
    geometry's handedness is arbitrary-but-consistent and must not reach a descriptor.
    """
    unknown = sorted(set(atoms["element"]) - set(MASS))
    assert not unknown, f"element(s) with no tabulated mass: {unknown}"
    atoms = atoms.assign(m=atoms["element"].map(MASS))
    assert atoms["m"].notna().all()
    rows = []
    for name, g in atoms.groupby("Molecule_Name", sort=False):
        m = g["m"].to_numpy()
        X = g[["x", "y", "z"]].to_numpy()
        com = (m[:, None] * X).sum(0) / m.sum()
        P = X - com
        I = np.zeros((3, 3))
        for mi, p in zip(m, P):
            I += mi * (np.dot(p, p) * np.eye(3) - np.outer(p, p))
        ev = np.sort(np.linalg.eigvalsh(I))                      # I1 <= I2 <= I3
        rg = float(np.sqrt((m * (P ** 2).sum(1)).sum() / m.sum()))
        rows.append((name, rg, float(ev[0]), float(ev[1]), float(ev[2]),
                     float(ev[0] / ev[2]), float(ev[1] / ev[2])))
    return pd.DataFrame(rows, columns=["Molecule_Name", "radius_of_gyration",
                                       "I1", "I2", "I3", "npr1", "npr2"])


def main() -> int:
    FEAT_DIR.mkdir(parents=True, exist_ok=True)
    print("=" * 78)
    print("PART 0 -- provenance")
    print("=" * 78)
    prov = {}
    for p in [QM_DESC, QM_MOL, QM_ATOMS, RDKIT2D_PATH, CURATED_PATH, BLIND_PATH]:
        prov[str(p.relative_to(REPO_ROOT))] = sha256(p)
        print(f"  {p.relative_to(REPO_ROOT)}  sha256 {prov[str(p.relative_to(REPO_ROOT))][:16]}...")

    d = pd.read_csv(QM_DESC)
    m = pd.read_csv(QM_MOL)
    atoms = pd.read_csv(QM_ATOMS)
    r2d = pd.read_csv(RDKIT2D_PATH)
    assert len(d) == len(m) == N_TOTAL, (len(d), len(m))
    assert len(atoms) == N_ATOMS, len(atoms)
    assert len(r2d) == N_TOTAL, len(r2d)
    print(f"\n  descriptors {d.shape}  molecular_raw {m.shape}  atoms {atoms.shape}  "
          f"rdkit2d {r2d.shape}")

    # ---------------------------------------------------------------- Part 1: identities
    print("\n" + "=" * 78)
    print("PART 1 -- the identities the block definitions depend on")
    print("=" * 78)
    j = d.merge(m, on="Molecule_Name", suffixes=("_D", "_M"), validate="1:1")
    for a, b, tol in [("dft_homo_ev", "homo_ev", 0.0), ("dft_lumo_ev", "lumo_ev", 0.0),
                      ("dft_gap_ev", "gap_ev", 0.0), ("dft_dipole_debye", "dipole_debye", 1e-3)]:
        diff = float(np.nanmax(np.abs(j[a] - j[b])))
        assert diff <= tol, f"{a} vs {b}: {diff:.3e} exceeds {tol}"
        print(f"  {a:26s} == {b:18s} max abs diff {diff:.3e}")
    gap_err = float(np.nanmax(np.abs(j["dft_gap_ev"] - (j["dft_lumo_ev"] - j["dft_homo_ev"]))))
    assert gap_err < 1e-12
    print(f"  dft_gap_ev == lumo - homo exactly ({gap_err:.1e}) -- gap is DERIVED, not independent")

    shp = shape_descriptors(atoms)
    assert len(shp) == N_TOTAL, len(shp)
    print(f"\n  shape block computed for {len(shp)} compounds")
    # The rotational constants are const/I on every axis, so they carry nothing a tree can use
    # beyond the inertia eigenvalues. Asserted, then excluded via FORBIDDEN.
    chk = shp.merge(m[["Molecule_Name", "rot_const_a_cm1", "rot_const_b_cm1",
                       "rot_const_c_cm1"]], on="Molecule_Name", validate="1:1")
    for rc, Icol in [("rot_const_a_cm1", "I1"), ("rot_const_b_cm1", "I2"),
                     ("rot_const_c_cm1", "I3")]:
        rho = float(spearmanr(chk[rc], chk[Icol]).statistic)
        prod = chk[rc] * chk[Icol]
        rel = float(prod.std() / prod.mean())
        assert rho < -0.99999 and rel < 1e-3, (rc, Icol, rho, rel)
        print(f"  {rc:18s} == const/{Icol}  (Spearman {rho:+.6f}, product rel. sd {rel:.1e})"
              f" -> EXCLUDED as a reparameterisation")

    # ---------------------------------------------------------------- Part 2: assemble
    print("\n" + "=" * 78)
    print("PART 2 -- assemble the QM columns")
    print("=" * 78)
    qm = d[["Molecule_Name", "inchikey", "set", "dft_homo_ev", "dft_lumo_ev", "dft_gap_ev",
            "gfn2_homo_ev", "gfn2_lumo_ev", "gfn2_gap_ev", "dft_dipole_debye",
            "hirshfeld_min"]].copy()
    qm = qm.merge(m[["Molecule_Name", "homo1_ev", "homo2_ev", "lumo1_ev",
                     "gfn2_relaxation_strain_eh"]], on="Molecule_Name", validate="1:1")
    qm = qm.merge(shp[["Molecule_Name"] + QM_3D], on="Molecule_Name", validate="1:1")
    qm["chi_electronegativity"] = -(qm["dft_homo_ev"] + qm["dft_lumo_ev"]) / 2.0
    qm["omega_electrophilicity"] = (qm["chi_electronegativity"] ** 2
                                   / (qm["dft_lumo_ev"] - qm["dft_homo_ev"]))
    qm["dft_minus_gfn2_gap"] = qm["dft_gap_ev"] - qm["gfn2_gap_ev"]
    qm["split"] = qm["set"].map({"train": "train", "blind": "test"})
    assert qm["split"].notna().all(), "unexpected value in the QM table's `set` column"
    assert (qm["split"] == "train").sum() == N_TRAIN
    assert (qm["split"] == "test").sum() == N_BLIND
    print(f"  QM table {qm.shape}  train={N_TRAIN} test={N_BLIND}")

    r2d_cols = [c for c in r2d.columns if c not in ("Molecule_Name", "inchikey", "split")]
    assert len(r2d_cols) == 217, len(r2d_cols)
    # rdkit2d carries SIXTEEN zero-variance columns, not one -- found by this assertion rather
    # than assumed, and more than the revision note recorded. All sixteen are identically 0:
    # NumRadicalElectrons (notebook 41's zero-radicals finding), two empty VSA bins (SMR_VSA8,
    # SlogP_VSA9), and thirteen rare-substructure fragment counters absent from this library
    # (azide, azo, barbiturate, diazo, epoxide, isocyanate, isothiocyanate, nitroso, phosphoric
    # acid/ester, primary sulfonamide, thiocyanate, C-S). So the "217 free descriptors" block is
    # 201 informative columns plus 16 constants. Harmless in a fit, but it matters twice: no
    # hand-picked control may draw one, and the random 18-column draws must sample from the
    # informative pool only, or a draw landing several constants would be artificially weak and
    # would bias the null band downward.
    const_cols = sorted(c for c in r2d_cols if r2d[c].nunique(dropna=False) <= 1)
    nonzero_const = [c for c in const_cols if float(r2d[c].iloc[0]) != 0.0]
    assert not nonzero_const, f"a constant column is not zero, which needs explaining: {nonzero_const}"
    print(f"  rdkit2d: {len(r2d_cols)} descriptor columns, {len(const_cols)} zero-variance "
          f"(all identically 0), {len(r2d_cols) - len(const_cols)} informative")
    print(f"    constant: {const_cols}")

    # Zero variance is not the only way a column can misbehave in a linear arm. `Ipc` spans
    # 9.65 to 5.06e+14 -- a 5e14 dynamic range, four orders of magnitude wider than the next
    # widest column (BertzCT, 2.7e+03) -- and without a guard it becomes a ~63-sigma leverage
    # point after standardisation. The ridge arm clips every column to its fit-portion p1/p99
    # before scaling, which brings Ipc's worst row to ~7 sigma; this assertion exists so that a
    # future column with that profile is REPORTED rather than silently absorbed, and so no
    # hand-picked control draws one.
    rng = {c: float(r2d[c].max() - r2d[c].min()) for c in r2d_cols}
    EXTREME_RANGE = 1e6
    extreme = sorted((c for c, v in rng.items() if v > EXTREME_RANGE), key=lambda c: -rng[c])
    print(f"  rdkit2d columns with dynamic range > {EXTREME_RANGE:.0e}: "
          f"{[(c, f'{rng[c]:.3e}') for c in extreme]}")
    assert set(extreme) <= {"Ipc"}, (
        f"a column other than Ipc has extreme dynamic range: {extreme}. The ridge arm's "
        f"fit-portion p1/p99 clip handles Ipc, but a new such column must be inspected before "
        f"it enters a linear arm, not discovered afterwards.")
    for nm, cols in [("RDKIT_MATCH18", RDKIT_MATCH18),
                     ("RDKIT_MATCH3_ELEC", RDKIT_MATCH3_ELEC),
                     ("RDKIT_MATCH3_NEW3", RDKIT_MATCH3_NEW3)]:
        missing = [c for c in cols if c not in r2d_cols]
        assert not missing, f"{nm}: columns absent from rdkit2d: {missing}"
        bad = [c for c in cols if c in const_cols]
        assert not bad, f"{nm}: draws a constant column: {bad}"
        wild = [c for c in cols if c in extreme]
        assert not wild, f"{nm}: draws an extreme-dynamic-range column: {wild}"
        assert len(cols) == len(set(cols)), f"{nm}: duplicate column"
        print(f"  {nm:18s} {len(cols):3d} cols, all present, none constant")

    rng = np.random.default_rng(RAND_BASE_SEED)
    draw_pool = [c for c in r2d_cols if c not in const_cols]
    rand_draws = {}
    for i in range(N_RAND_DRAWS):
        seed = RAND_BASE_SEED + i
        rand_draws[f"RDKIT_RAND18_{i:02d}"] = {
            "seed": seed,
            "columns": sorted(np.random.default_rng(seed).choice(
                draw_pool, size=18, replace=False).tolist()),
        }
    print(f"  {N_RAND_DRAWS} random 18-column draws, seeds {RAND_BASE_SEED}"
          f"..{RAND_BASE_SEED + N_RAND_DRAWS - 1} (descriptive band only -- see prereg 4.4)")

    # ---------------------------------------------------------------- Part 3: write blocks
    print("\n" + "=" * 78)
    print("PART 3 -- write the feature blocks")
    print("=" * 78)
    keyed = qm[["Molecule_Name", "inchikey", "split"]].copy()
    r2d_feat = r2d[["Molecule_Name"] + r2d_cols]

    def qm_block(cols):
        return keyed.merge(qm[["Molecule_Name"] + cols], on="Molecule_Name", validate="1:1")

    def mixed_block(cheap_cols, qm_cols):
        out = keyed.merge(r2d_feat[["Molecule_Name"] + cheap_cols], on="Molecule_Name",
                          validate="1:1")
        if qm_cols:
            out = out.merge(qm[["Molecule_Name"] + qm_cols], on="Molecule_Name", validate="1:1")
        return out

    blocks = {
        "QM_ELEC": qm_block(QM_ELEC),
        "QM_3D": qm_block(QM_3D),
        "QM_FULL": qm_block(QM_FULL),
        "QM_NEW3": qm_block(QM_NEW3),
        "QM_DFT3": qm_block(QM_DFT3),
        "QM_GFN23": qm_block(QM_GFN23),
        "QM_FULL_HMIN": qm_block(QM_FULL + HMIN),
        "RDKIT2D": mixed_block(r2d_cols, []),
        "RDKIT2D_QM": mixed_block(r2d_cols, QM_FULL),
        "RDKIT2D_3D": mixed_block(r2d_cols, QM_3D),
        "RDKIT2D_NEW3": mixed_block(r2d_cols, QM_NEW3),
        "RDKIT2D_ELEC": mixed_block(r2d_cols, QM_ELEC),
        "RDKIT_MATCH18": mixed_block(RDKIT_MATCH18, []),
        "RDKIT_MATCH3_ELEC": mixed_block(RDKIT_MATCH3_ELEC, []),
        "RDKIT_MATCH3_NEW3": mixed_block(RDKIT_MATCH3_NEW3, []),
    }
    for nm, spec in rand_draws.items():
        blocks[nm] = mixed_block(spec["columns"], [])

    manifest = {}
    for nm, df in sorted(blocks.items()):
        feat_cols = [c for c in df.columns if c not in ("Molecule_Name", "inchikey", "split")]
        hit = sorted(set(feat_cols) & FORBIDDEN)
        assert not hit, f"{nm} contains forbidden column(s): {hit}"
        assert len(df) == N_TOTAL, f"{nm}: {len(df)} rows"
        assert df["Molecule_Name"].is_unique and df["inchikey"].is_unique
        nulls = int(df[feat_cols].isna().to_numpy().sum())
        assert nulls == 0, f"{nm}: {nulls} null feature values"
        assert np.isfinite(df[feat_cols].to_numpy(dtype=np.float64)).all(), f"{nm}: non-finite"
        path = FEAT_DIR / f"{nm}.csv"
        df.to_csv(path, index=False)
        manifest[nm] = {"n_cols": len(feat_cols), "columns": feat_cols,
                        "path": str(path.relative_to(REPO_ROOT)), "sha256": sha256(path)}
        tier1_only = nm.startswith("RDKIT_RAND18_")
        print(f"  {nm:22s} {len(feat_cols):4d} cols  {len(df)} rows  0 nulls"
              f"{'   [tier-1 only]' if tier1_only else ''}")

    for nm, want in [("QM_ELEC", 14), ("QM_3D", 4), ("QM_FULL", 18), ("QM_NEW3", 3),
                     ("QM_DFT3", 3), ("QM_GFN23", 3), ("QM_FULL_HMIN", 19),
                     ("RDKIT2D", 217), ("RDKIT2D_QM", 235), ("RDKIT2D_3D", 221),
                     ("RDKIT2D_NEW3", 220), ("RDKIT2D_ELEC", 231), ("RDKIT_MATCH18", 18),
                     ("RDKIT_MATCH3_ELEC", 3), ("RDKIT_MATCH3_NEW3", 3)]:
        got = manifest[nm]["n_cols"]
        assert got == want, f"{nm}: expected {want} columns, built {got}"
    print("\n  all block widths match the plan's own table")

    # ---------------------------------------------------------------- Part 4: prereg
    print("\n" + "=" * 78)
    print("PART 4 -- freeze prereg.json (BEFORE any training run)")
    print("=" * 78)
    pt = pd.read_csv(NB37_OUT / "paired_tests_vs_plain_25fold.csv")
    se = (pt[pt["metric"] == "ST-RAE"].groupby("isoform")["se"].mean().round(4).to_dict())
    bars = {iso: round(2 * se[iso], 3) for iso in ISOFORMS}
    blend = pd.read_csv(NB37_OUT / "blend_results.csv").set_index("isoform")
    nested = {iso: float(blend.loc[iso, "raw_nested"]) for iso in ISOFORMS}
    optimism = {iso: round(float(blend.loc[iso, "raw_nested"]
                                - blend.loc[iso, "raw_in_sample"]), 4) for iso in ISOFORMS}
    print(f"  measured paired SE (nb37, mean over 10 arms): {se}")
    print(f"  materiality bar = 2 x SE:                     {bars}")
    print(f"  11-arm nested blend ST-RAE:                   "
          f"{ {k: round(v, 4) for k, v in nested.items()} }")
    print(f"  measured blend optimism:                      {optimism}")

    prereg = {
        "written_at_utc": datetime.now(timezone.utc).isoformat(),
        "written_by": "scripts/44_build_qm_feature_blocks.py",
        "status": "FROZEN before any training run",
        "question": ("Does the molecule-level QM block carry information that improves pIC50 "
                     "prediction over what this repo already computes for free, and if so does "
                     "it belong as features in a single model or as a decorrelated ensemble "
                     "member?"),
        "source_sha256": prov,
        "blocks": manifest,
        "random_draws": rand_draws,
        "forbidden_columns": sorted(FORBIDDEN),
        "rdkit2d_zero_variance_columns": {
            "n_descriptor_columns": 217,
            "n_zero_variance": None,          # filled below
            "n_informative": None,            # filled below
            "columns": None,                  # filled below
            "note": ("All are identically 0: NumRadicalElectrons (notebook 41's zero-radicals "
                     "finding), two empty VSA bins, and thirteen rare-substructure fragment "
                     "counters absent from this 5,655-compound library. The free block is "
                     "therefore narrower than its column count suggests. No hand-picked control "
                     "draws one, and the 20 random draws sample the informative pool only -- "
                     "otherwise a draw landing several constants would be artificially weak and "
                     "would bias the descriptive null band downward."),
        },
        "model_selection_rule": (
            "The primary model type M* is selected on the CONTROL block RDKIT2D at tier-1, "
            "before any QM arm is examined. Recorded to model_selection.json by "
            "44_run_qm_arms.py --tier1."),
        "resolution": {"paired_se_nb37": se, "materiality_bar_2se": bars,
                       "basis": ("2 x that isoform's own mean paired ST-RAE SE measured across "
                                 "notebook 37's 10 arms x 25 folds")},
        "criterion_I": {
            "I.1_PRIMARY_headline": {
                "test": "rdkit2d_qm__M* vs rdkit2d__M*, paired over the 25 common folds",
                "question": "does the whole QM block add to what is already free?",
                "procedure": ("notebook 37's own paired_test: Levene/Brown-Forsythe homogeneity "
                              "check -> paired t or Wilcoxon signed-rank; BH across the 4 "
                              "isoforms, in this family only"),
                "pass_per_isoform": "mean_diff < 0 AND p_BH < 0.05 AND |mean_diff| >= bar",
                "overall_pass": ("at least 1 isoform passes AND the passing set is not {CYP2D6} "
                                 "alone (notebook 21: OOF ST-RAE is not a valid ranking signal "
                                 "for CYP2D6; notebook 38 reproduced the mismatch on CYP3A4)"),
            },
            "I.2_PRIMARY": {
                "test": "rdkit2d_3d__M* vs rdkit2d__M*", "bh_family": "its own, 4 isoforms",
                "question": "does the 3D shape block add to what is free?"},
            "I.3_PRIMARY": {
                "test": "rdkit2d_new3__M* vs rdkit2d__M*", "bh_family": "its own, 4 isoforms",
                "question": ("do the three columns no cheap model reproduces above R^2 0.32 "
                             "(dipole, NPR1, NPR2) add to what is free?")},
            "three_primaries_note": (
                "Three questions, so three primary tests, each with its own BH family of 4 -- "
                "consistent with the one-test-per-question rule. Three independent families at "
                "alpha=0.05 put the chance of at least one spurious PASS across the three "
                "questions at UP TO ~14% rather than 5% -- up to, because the three blocks "
                "overlap heavily (QM_FULL contains both others), share all 25 folds and share "
                "the rdkit2d comparator arm, so the tests are positively dependent and the true "
                "figure is below 14%. I.1 is the headline; a PASS on I.2 or I.3 with I.1 failing "
                "is a NARROWER claim and may not be reported as 'the QM block helps'."),
            "I.4_sensitivity_cannot_PASS": {
                "test": "qm_full__M* vs rdkit_match18__M*",
                "why_weaker": ("every tabular arm in notebook 37's design is significantly worse "
                               "than plain on CYP2C9/CYP3A4 by +0.062 to +0.297, so 'QM beats a "
                               "cheap block' can pass while the arm stays far worse than the "
                               "model that matters")},
            "I.5_co_reported_metric": ("Spearman, same folds, same test -- the one CV quantity "
                                       "whose board counterpart no affine correction can move. "
                                       "In principle only: notebook 30's aid arm had a better "
                                       "OOF CYP3A4 Spearman than plain and a worse board "
                                       "Spearman (0.7651 vs 0.8177)."),
            "I.6_model_class_robustness": ("rf / lightgbm / xgboost / ridge all reported; a "
                                           "primary PASS holding in <=1 of 4 is reported as "
                                           "model-class-specific and cannot change the verdict"),
        },
        "criterion_II_ensemble": {
            "II.1_PRIMARY": {
                "test": ("add the QM arm to notebook 37's 11-arm pool and re-run its own nested "
                         "greedy selection: fit members AND weights on 24 folds, score the "
                         "held-out one, 25 times"),
                "nested_11arm_strae": nested, "measured_optimism": optimism,
                "required_improvement": {iso: round(optimism[iso] + bars[iso], 4)
                                         for iso in ISOFORMS},
                "pass": ("bar cleared on >=2 of 4 isoforms, at least one of which is not "
                         "CYP2D6")},
            "II.2_membership_is_not_the_criterion": (
                "The third blend slot is unidentified at n=25: it swapped member on 2 of 4 "
                "isoforms between the 15- and 25-fold versions while the nested score moved "
                "<0.005. Non-zero weight is reported as weak evidence, never as the criterion. "
                "CYP2D6 is also assessed under notebook 39's shared-weight objective, since 39 "
                "showed its collapse to deadzone 1.0 is a property of the per-isoform objective."),
            "II.3_decorrelation_descriptive": {
                "quantity": ("residual (prediction - true) correlation, never raw-prediction "
                             "correlation -- CLAUDE.md's own rule, Dietterich 2000"),
                "plain_reseed_floor": {"CYP1A2": 0.9227, "CYP2C9": 0.8949,
                                       "CYP2D6": 0.9210, "CYP3A4": 0.8847},
                "floor_sd": {"CYP1A2": 0.0070, "CYP2C9": 0.0067,
                             "CYP2D6": 0.0065, "CYP3A4": 0.0082},
                "tree_arm_band": {"CYP1A2": [0.8125, 0.8683], "CYP2C9": [0.7375, 0.7994],
                                  "CYP2D6": [0.8502, 0.8896], "CYP3A4": [0.7093, 0.7751]},
                "nb24_pool_band_per_isoform_mean": {"CYP1A2": 0.8845, "CYP2C9": 0.8422,
                                                    "CYP2D6": 0.9078, "CYP3A4": 0.7918},
                "criterion": ("'inside or below' that isoform's tree-arm band -- 'below' is a "
                              "bar no existing arm meets, because the tree arms define the "
                              "band. Descriptive, not a gate. Margins under ~0.013 (2 SD) are "
                              "not separable from re-seed noise."),
            },
        },
        "criterion_III_negative_result_predeclared": (
            "If I.1 and I.2 both fail, the recorded conclusion is: the molecule-level QM block "
            "does not carry pIC50-relevant information beyond free 2D descriptors on this "
            "dataset at this level of theory -- with the cost attached (281.6 core-hours, 47.03 h "
            "wall). A publishable negative; notebook 18 is the precedent for reporting exactly "
            "this and not promoting. Part 2.3 makes it the more likely outcome for the "
            "electronic block, which is why it is written down now rather than after."),
        "null_band_is_descriptive_only": (
            "With B draws the smallest achievable permutation p is 1/(B+1). BH across 4 isoforms "
            "puts the rank-1 threshold at 0.0125, so a formal permutation test needs B >= 79. "
            f"B = {N_RAND_DRAWS} gives a floor of {1/(N_RAND_DRAWS+1):.4f} and CANNOT support a "
            "BH-corrected permutation test. The band shows only whether RDKIT_MATCH18 was a "
            "lucky or unlucky draw, and it runs at tier-1 on the random partition so it does "
            "not occupy 20 of the frozen partition's permanent arms."),
        "frozen_display_subset": {
            "arms": ["plain", "deadzone", "aid", "rdkit2d", "rdkit2d_qm", "rdkit2d_3d",
                     "rdkit2d_new3", "rdkit_match18", "plain_qm"],
            "why_frozen_here": (
                "Membership is frozen by name before any result exists. This is STRICTER than "
                "notebook 05 rather than inherited from it: 05's curated_subset() picks 3 of its "
                "6 slots by argmin over already-computed mean ST-RAE, recomputed per isoform. "
                "Plot ORDERING may be result-derived (notebooks 34/35/37 do that); membership "
                "may not. The full-pool table and full-arm figures are saved alongside."),
        },
        "bh_family_rule": (
            "The new arms form their OWN BH family. Notebook 37's p_BH values are read from its "
            "saved paired_tests_vs_plain_25fold.csv and NEVER recomputed: 37's verdicts were "
            "BH-corrected within a family of 10 arms, so pooling new arms into that family would "
            "retroactively change already-published verdicts."),
        "not_a_submission_decision": (
            "No arm from this design will be submitted on CV evidence. Four CV-to-board transfer "
            "failures are on record; the diagnostic that did predict board movement is "
            "prediction-spread compression of a submitted column, which this design cannot "
            "produce and does not claim to."),
        "decisions_taken": [
            {"decision": "hirshfeld_min assigned to group (c), not (a)",
             "reason": ("R^2 0.885 against the rule's 'either model >= 0.90' bar -- within "
                        "measurement noise of it. Carried as the QM_FULL_HMIN sensitivity arm "
                        "rather than resolved by fiat."),
             "alternative": "assign to (a) and include it in QM_ELEC"},
            {"decision": "group (a') algebraic derivatives included in QM_ELEC",
             "reason": ("zero extra information but a real reparameterisation for ridge and a "
                        "mild one for trees"),
             "alternative": "the 9 raw electronic columns only (the stated sensitivity)"},
            {"decision": f"B = {N_RAND_DRAWS} random draws",
             "reason": ("the brief asked for a band without specifying B; this value is the "
                        "plan's own choice, and the plan then caught its own permutation floor. "
                        "That is the plan correcting itself, not a correction of the brief."),
             "alternative": "B >= 79, which would support a formal permutation test (~5.6 h)"},
            {"decision": "feature files written under outputs/, not data/processed/",
             "reason": ("data/processed/ is on the protected list of every recent notebook's "
                        "scope check; writing there would make this work indistinguishable from "
                        "a violation"),
             "alternative": "add them to data/processed/ and the shared FEATURE_FILES registry"},
        ],
    }
    prereg["rdkit2d_zero_variance_columns"].update(
        {"n_zero_variance": len(const_cols), "n_informative": len(r2d_cols) - len(const_cols),
         "columns": const_cols})
    (OUT / "prereg.json").write_text(json.dumps(prereg, indent=2))
    print(f"\n  wrote {(OUT / 'prereg.json').relative_to(REPO_ROOT)}")
    print(f"  sha256 {sha256(OUT / 'prereg.json')}")
    print("\nDONE -- no model trained, no label read, no submission artefact built.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
