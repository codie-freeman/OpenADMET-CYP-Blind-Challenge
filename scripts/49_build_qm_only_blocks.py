"""QM-only feature blocks, and blind-extended copies of notebook 46's six.

WHY THIS EXISTS. Notebook 46 compared six blocks and every one of them is notebook 45's 226 free
2D descriptors with 13 to 17 QM columns bolted on. Notebook 49's diversity analysis found those
six arms correlate 0.995 to 0.997 with `free__ridge`, which carries no quantum content at all, so
what notebook 46 measured is the free descriptors' correlation profile and not the quantum data's.
A QM-ONLY arm has never been run anywhere in this project. It is the clean reading of whether the
quantum data is decorrelated from notebook 48's ensemble pool, which is the only question left
open about this block.

TWO NEW BLOCKS, no free descriptors in either:

    QM13         13  notebook 44's genuinely-new molecule-level QM columns, alone
    QM13_ATOM4   17  + notebook 45's four full-coverage strongly-resistant atom descriptors

SCOPE. Nothing here modifies `scripts/46_build_feature_blocks.py` or anything under
`outputs/46_qm_arm/`. Notebook 46's own six blocks are REBUILT here only as a faithfulness check:
the rebuild must reproduce its saved files byte for byte on the 4,905 training rows, which is what
licenses the blind-extended copies written alongside.

BLIND EXTENSION. All four upstream sources already cover all 5,655 compounds, so extending costs
nothing but changing which keys the merge starts from. The 4,905 restriction in notebook 46's
builder is one line (`keys = curated[...]`). Two outputs are written per block:

    features/{BLOCK}.csv                 4,905 training rows, the contract notebook 46's runner
                                         asserts, used for the 25-fold out-of-fold run
    features_blind_extended/{BLOCK}.csv  5,655 rows, used ONLY if the decorrelation gate passes

No model is fitted here and no blind prediction is produced. The gate in
`outputs/49_qm_only_arm/gate.json` has to pass first.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
NB46 = REPO_ROOT / "outputs" / "46_qm_arm"
OUT = REPO_ROOT / "outputs" / "49_qm_only_arm"
FEAT_DIR = OUT / "features"
FEAT_ALL_DIR = OUT / "features_blind_extended"

CURATED = REPO_ROOT / "data" / "processed" / "train_inhibition_curated.csv"
BLIND = REPO_ROOT / "data" / "processed" / "test_blinded_curated.csv"
RDKIT2D = REPO_ROOT / "data" / "processed" / "tabular_baseline_features_rdkit2d.csv"
ECFP4 = REPO_ROOT / "data" / "processed" / "tabular_baseline_features.csv"
QM_FULL = REPO_ROOT / "outputs" / "44_qm_molecular_block" / "features" / "QM_FULL.csv"
AGGREGATES = REPO_ROOT / "outputs" / "qm_parsed" / "qm_atom_aggregates.csv"
NB45_VERDICT = REPO_ROOT / "outputs" / "45_atom_level_screen" / "verdict_table.csv"
BUTINA = REPO_ROOT / "data" / "folds" / "cv_folds_butina_5x5.csv"
NB34_FOLD_SHA256 = "52f39b1a1cea9200b58dce0a061d4a21ce94d9d58ec77067cf6508b0f4b3e87e"

# all three literals below are notebook 46's own, reused verbatim
NINE = ["mol_wt", "hba", "hbd", "ring_count", "stereocenter_count", "formal_charge",
        "logp", "logd_proxy", "vsa_acc_proxy"]
DERIVED5 = ["dft_gap_ev", "gfn2_gap_ev", "chi_electronegativity", "omega_electrophilicity",
            "dft_minus_gfn2_gap"]
FOUR = ["hirshfeld_C_min", "hirshfeld_H_min", "bond_order_min", "mayer_va_min"]
HALOGEN = ["hirshfeld_Cl_min", "hirshfeld_Cl_max", "hirshfeld_F_min", "hirshfeld_F_max"]
CONTROL = "hirshfeld_range"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(keys: pd.DataFrame, expect_rows: int):
    """Notebook 46's own merge chain, with `keys` as the only variable."""
    r2d = pd.read_csv(RDKIT2D)
    rc = [c for c in r2d.columns if c not in ("Molecule_Name", "inchikey", "split")]
    ecfp = pd.read_csv(ECFP4, usecols=["Molecule_Name"] + NINE)
    free = (keys.merge(r2d[["Molecule_Name"] + rc], on="Molecule_Name", validate="1:1")
                .merge(ecfp.rename(columns={c: "n9_" + c for c in NINE}),
                       on="Molecule_Name", validate="1:1"))
    free_cols = rc + ["n9_" + c for c in NINE]
    assert len(free_cols) == 226, len(free_cols)

    qm = pd.read_csv(QM_FULL)
    thirteen = [c for c in qm.columns
                if c not in ("Molecule_Name", "inchikey", "set", "split", *DERIVED5)]
    assert len(thirteen) == 13, (len(thirteen), thirteen)

    agg = pd.read_csv(AGGREGATES)
    verdict = pd.read_csv(NB45_VERDICT)
    four_check = sorted(verdict[(verdict.tier == "strong") & (verdict.pct_blind > 99)].descriptor)
    assert four_check == sorted(FOUR), (four_check, FOUR)
    assert verdict.loc[verdict.descriptor == CONTROL, "resists_both"].iloc[0] == False, \
        "the control must be a descriptor notebook 45's screen REJECTED"

    base = free.merge(qm[["Molecule_Name"] + thirteen], on="Molecule_Name", validate="1:1") \
               .merge(agg[["Molecule_Name"] + FOUR + HALOGEN + [CONTROL]],
                      on="Molecule_Name", validate="1:1")
    assert len(base) == expect_rows, (len(base), expect_rows)
    for el in ["Cl", "F"]:
        base[f"has_{el}"] = base[f"hirshfeld_{el}_min"].notna().astype(float)
    hal_cols = HALOGEN + ["has_Cl", "has_F"]

    blocks = {
        # notebook 46's six, rebuilt verbatim so the rebuild can be checked against its files
        "FREE": free_cols,
        "FREE_QM13": free_cols + thirteen,
        "FREE_ATOM4": free_cols + FOUR,
        "FREE_QM13_ATOM4": free_cols + thirteen + FOUR,
        "FREE_ATOM4_HAL": free_cols + FOUR + hal_cols,
        "FREE_CONTROL": free_cols + [CONTROL],
        # the two new ones: no free descriptors at all
        "QM13": thirteen,
        "QM13_ATOM4": thirteen + FOUR,
    }
    return base, blocks


def write_blocks(base, blocks, out_dir: Path, may_have_nan=("FREE_ATOM4_HAL",)):
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for name, cols in blocks.items():
        df = base[["Molecule_Name", "inchikey"] + cols]
        X = df[cols].to_numpy(dtype=np.float64)
        nan_cols = [c for c in cols if df[c].isna().any()]
        if name in may_have_nan:
            assert set(nan_cols) <= set(HALOGEN), (name, nan_cols)
        else:
            assert not nan_cols, f"{name}: unexpected NaN in {nan_cols}"
            assert np.isfinite(X).all(), f"{name}: non-finite value"
        path = out_dir / f"{name}.csv"
        df.to_csv(path, index=False)
        manifest[name] = {"n_cols": len(cols), "n_rows": len(df), "columns": cols,
                          "nan_columns": nan_cols,
                          "n_constant": int(sum(df[c].nunique(dropna=False) <= 1 for c in cols)),
                          "sha256": sha256(path)}
        print(f"  {name:18s} {len(cols):4d} cols  {len(df):5d} rows  "
              f"nan_cols={len(nan_cols)}  constant={manifest[name]['n_constant']}")
    return manifest


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    assert sha256(BUTINA) == NB34_FOLD_SHA256, "frozen partition is not the file 34/35/37 ran on"
    print(f"frozen Butina partition sha256 matches notebook 34's record")

    curated = pd.read_csv(CURATED)
    blind = pd.read_csv(BLIND)
    assert len(curated) == 4905 and len(blind) == 750
    keys_train = curated[["Molecule_Name", "inchikey"]].copy()
    keys_all = pd.concat([keys_train, blind[["Molecule_Name", "inchikey"]]], ignore_index=True)
    assert len(keys_all) == 5655 and keys_all.Molecule_Name.is_unique

    print(f"\n--- 4,905 training rows (notebook 46's own contract) ---")
    base_tr, blocks_tr = build(keys_train, 4905)
    man_tr = write_blocks(base_tr, blocks_tr, FEAT_DIR)

    # FAITHFULNESS CHECK. The rebuild of notebook 46's own six must reproduce its saved files
    # exactly, otherwise nothing downstream of this builder is comparable to notebook 46.
    print("\n--- faithfulness: rebuilt blocks against notebook 46's saved files ---")
    identical = {}
    for name in ["FREE", "FREE_QM13", "FREE_ATOM4", "FREE_QM13_ATOM4", "FREE_ATOM4_HAL",
                 "FREE_CONTROL"]:
        a, b = sha256(FEAT_DIR / f"{name}.csv"), sha256(NB46 / "features" / f"{name}.csv")
        identical[name] = (a == b)
        print(f"  {name:18s} sha256 {'IDENTICAL' if a == b else 'DIFFERS -- ' + a[:12]}")
        assert a == b, f"{name}: rebuild does not reproduce notebook 46's saved block"

    print(f"\n--- 5,655 rows, blind-extended (step 4 only, if the gate passes) ---")
    base_all, blocks_all = build(keys_all, 5655)
    man_all = write_blocks(base_all, blocks_all, FEAT_ALL_DIR)

    # the blind extension must not have changed a single training value
    print("\n--- blind extension changed nothing for the 4,905 training compounds ---")
    for name in blocks_tr:
        a = pd.read_csv(FEAT_DIR / f"{name}.csv")
        b = pd.read_csv(FEAT_ALL_DIR / f"{name}.csv")
        b = b[b.Molecule_Name.isin(set(a.Molecule_Name))].reset_index(drop=True)
        assert list(a.columns) == list(b.columns) and len(a) == len(b), name
        num = [c for c in a.columns if c not in ("Molecule_Name", "inchikey")]
        d = np.nanmax(np.abs(a[num].to_numpy(float) - b[num].to_numpy(float)))
        nan_same = bool((a[num].isna().to_numpy() == b[num].isna().to_numpy()).all())
        print(f"  {name:18s} max abs diff {d:.3e}  NaN pattern identical: {nan_same}")
        assert d == 0.0 and nan_same, name

    # blind coverage of the two new blocks, reported because notebook 43's bind was coverage
    print("\n--- blind coverage of the new QM-only blocks ---")
    bset = set(blind.Molecule_Name)
    for name in ["QM13", "QM13_ATOM4"]:
        d = pd.read_csv(FEAT_ALL_DIR / f"{name}.csv")
        db = d[d.Molecule_Name.isin(bset)]
        cols = [c for c in d.columns if c not in ("Molecule_Name", "inchikey")]
        print(f"  {name:12s} blind rows {len(db)}/750, complete on "
              f"{int(db[cols].notna().all(axis=1).sum())}/750")
        assert len(db) == 750 and db[cols].notna().all().all()

    json.dump({"written_at_utc": datetime.now(timezone.utc).isoformat(),
               "written_by": "scripts/49_build_qm_only_blocks.py",
               "partition_sha256": NB34_FOLD_SHA256,
               "nb46_rebuild_byte_identical": identical,
               "blocks_train_4905": man_tr, "blocks_all_5655": man_all},
              open(OUT / "feature_manifest.json", "w"), indent=1)
    print(f"\nwrote {OUT / 'feature_manifest.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
