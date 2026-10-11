#!/usr/bin/env python
"""Build and verify the per-compound ORCA-to-chemprop atom permutation. Notebook 52's stage 0.

Why this exists. The QM run (`scripts/run_qm_descriptors.py`) featurised every compound as
`Chem.AddHs(Chem.MolFromSmiles(row.SMILES))` -- the RAW `SMILES` column -- and wrote atoms in
RDKit's native index order, so an ORCA atom index is an index into that molecule. Chemprop, as
this project runs it, is given `-s canonical_smiles`, and `canonical_smiles` reproduces ORCA's
element sequence on only 4,212 of 5,655 compounds (notebook 43, reproduced here as a check).
Feeding a per-atom array in ORCA's order to a graph built from `canonical_smiles` would put
Hirshfeld charges on the wrong atoms on a quarter of the library, silently.

Two further facts the permutation has to absorb, both established here rather than assumed:

  * `Chem.AddHs` appends hydrogens AFTER the heavy atoms, so heavy-atom indices 0..n_heavy-1 are
    preserved between `MolFromSmiles(SMILES)` and `AddHs(MolFromSmiles(SMILES))`. Asserted per
    compound, so ORCA's first n_heavy rows are exactly the heavy atoms.
  * Chemprop does NOT add hydrogens (`add_h: False`, `keep_h: False` in its own resolved args),
    so its graph carries heavy atoms only. The permutation therefore maps chemprop's n_heavy
    atom indices into ORCA's first n_heavy rows and the hydrogen rows are dropped. That is a real
    loss of information and notebook 52 says so; it is not something this script can fix.

How the permutation is built. Both molecules are canonicalised through the same
`Chem.MolToSmiles` call, which records `_smilesAtomOutputOrder`: position k of the shared
canonical SMILES string is raw-mol atom `o_raw[k]` and canonical-mol atom `o_can[k]`. So
chemprop atom `o_can[k]` is ORCA atom `o_raw[k]`. This routes both orderings through one
canonicalisation rather than trying to match them against each other, which is why it is robust
to symmetry: RDKit's own tie-breaking decides, identically for both sides.

Verification, all of it a hard failure rather than a warning (the brief's requirement, and the
reason notebook 43's own assertions existed):

  (1) AddHs preserves heavy-atom indices.
  (2) The permutation is a bijection on 0..n_heavy-1.
  (3) Element at chemprop index j == element at ORCA index perm[j]. Notebook 43's check.
  (4) Every chemprop bond, mapped through the permutation, is a bond of the ORCA molecule.
      Notebook 43's bond check, in its subset form -- the element sequence alone cannot see a
      swap of two same-element atoms, which the bond check can.
  (5) A NEGATIVE CONTROL: a deliberately corrupted permutation must fail (3) or (4). Without
      this the assertions could be vacuous, which is the failure mode that makes an ordering bug
      invisible.

Usage:
  python scripts/52_build_atom_ordering.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

OUT = REPO_ROOT / "outputs" / "52_chemprop_qm_atom_features"
CURATED = REPO_ROOT / "data" / "processed" / "train_inhibition_curated.csv"
BLIND = REPO_ROOT / "data" / "processed" / "test_blinded_curated.csv"
ATOMS = REPO_ROOT / "outputs" / "qm_parsed" / "qm_atoms.csv"

ORDERING_CSV = OUT / "atom_ordering.csv"
REPORT_JSON = OUT / "ordering_verification.json"

# Notebook 43's own figure, cross-checked rather than quoted: the number of compounds whose
# `canonical_smiles` reproduces ORCA's element sequence outright.
NB43_ELEMENT_SEQUENCE_MATCHES = 4212


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def output_order(mol) -> list[int]:
    """RDKit's `_smilesAtomOutputOrder` for `mol`, as a list of atom indices."""
    Chem.MolToSmiles(mol)
    raw = mol.GetProp("_smilesAtomOutputOrder")
    return [int(t) for t in raw.strip()[1:-1].rstrip(",").split(",")]


def bonds_of(mol) -> set:
    return {(min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
             max(b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in mol.GetBonds()}


def permutation_for(raw_smiles: str, canonical_smiles: str) -> tuple[np.ndarray, dict]:
    """Return (perm, facts). `perm[j]` is the ORCA atom index of chemprop atom `j`.

    Raises on any failure of checks (1)-(4). Never returns a partially-verified permutation.
    """
    raw = Chem.MolFromSmiles(raw_smiles)
    can = Chem.MolFromSmiles(canonical_smiles)
    if raw is None or can is None:
        raise ValueError("RDKit could not parse a curated SMILES")
    n_heavy = raw.GetNumAtoms()
    if can.GetNumAtoms() != n_heavy:
        raise ValueError(f"heavy-atom count differs: raw {n_heavy} vs canonical {can.GetNumAtoms()}")

    raw_h = Chem.AddHs(raw)
    # (1) AddHs must preserve heavy-atom indices, so ORCA rows 0..n_heavy-1 are the heavy atoms.
    for i in range(n_heavy):
        if raw_h.GetAtomWithIdx(i).GetSymbol() != raw.GetAtomWithIdx(i).GetSymbol():
            raise ValueError(f"AddHs did not preserve heavy-atom index {i}")

    o_raw, o_can = output_order(raw), output_order(can)
    perm = np.empty(n_heavy, dtype=np.int64)
    perm[np.asarray(o_can, dtype=np.int64)] = np.asarray(o_raw, dtype=np.int64)

    # (2) bijection
    if sorted(perm.tolist()) != list(range(n_heavy)):
        raise ValueError("permutation is not a bijection on the heavy atoms")

    # (3) element sequence under the permutation -- notebook 43's check
    for j in range(n_heavy):
        if can.GetAtomWithIdx(j).GetSymbol() != raw_h.GetAtomWithIdx(int(perm[j])).GetSymbol():
            raise ValueError(f"element mismatch at chemprop index {j}")

    # (4) bond subset under the permutation -- notebook 43's check
    mapped = {(min(int(perm[i]), int(perm[j])), max(int(perm[i]), int(perm[j])))
              for i, j in bonds_of(can)}
    if not mapped <= bonds_of(raw_h):
        raise ValueError(f"{len(mapped - bonds_of(raw_h))} mapped bonds are not bonds of the "
                         f"ORCA molecule")

    raw_h_elements = [a.GetSymbol() for a in raw_h.GetAtoms()]
    can_h_elements = [a.GetSymbol() for a in Chem.AddHs(can).GetAtoms()]
    return perm, {
        "n_heavy": n_heavy,
        "n_with_h": raw_h.GetNumAtoms(),
        "n_hydrogens": raw_h.GetNumAtoms() - n_heavy,
        "is_identity": bool((perm == np.arange(n_heavy)).all()),
        "element_sequence_already_matched": raw_h_elements == can_h_elements,
    }


def negative_control(raw_smiles: str, canonical_smiles: str) -> bool:
    """Corrupt a verified permutation by swapping two entries and confirm the checks catch it.

    Returns True if the corrupted permutation is rejected. A compound whose heavy atoms are all
    one element AND whose swap happens to preserve the bond set would legitimately pass, so the
    caller tries several compounds and requires the control to bite on at least one.
    """
    raw = Chem.MolFromSmiles(raw_smiles)
    can = Chem.MolFromSmiles(canonical_smiles)
    raw_h = Chem.AddHs(raw)
    n = raw.GetNumAtoms()
    if n < 2:
        return False
    perm, _ = permutation_for(raw_smiles, canonical_smiles)
    bad = perm.copy()
    bad[0], bad[1] = bad[1], bad[0]
    elem_fails = any(can.GetAtomWithIdx(j).GetSymbol()
                     != raw_h.GetAtomWithIdx(int(bad[j])).GetSymbol() for j in range(n))
    mapped = {(min(int(bad[i]), int(bad[j])), max(int(bad[i]), int(bad[j])))
              for i, j in bonds_of(can)}
    return bool(elem_fails or not mapped <= bonds_of(raw_h))


def main() -> int:
    RDLogger.DisableLog("rdApp.*")
    OUT.mkdir(parents=True, exist_ok=True)

    curated = pd.read_csv(CURATED)
    blind = pd.read_csv(BLIND)
    assert len(curated) == 4905, len(curated)
    assert len(blind) == 750, len(blind)

    frames = [curated[["Molecule_Name", "SMILES", "canonical_smiles"]].assign(set="train"),
              blind[["Molecule_Name", "SMILES", "canonical_smiles"]].assign(set="blind")]
    both = pd.concat(frames, ignore_index=True)
    assert both["Molecule_Name"].is_unique
    print(f"compounds: {len(both)} ({len(curated)} train + {len(blind)} blind)")

    # ORCA's own per-atom table, for the one check that ties the permutation to the real file
    # rather than to a re-derivation of what the QM run should have written.
    atoms = pd.read_csv(ATOMS, usecols=["Molecule_Name", "atom_index", "element"])
    print(f"qm_atoms.csv: {len(atoms)} atom rows, sha256 {sha256(ATOMS)[:16]}")
    atom_elements = {name: g.sort_values("atom_index")["element"].tolist()
                     for name, g in atoms.groupby("Molecule_Name", sort=False)}
    assert len(atom_elements) == len(both), (len(atom_elements), len(both))

    rows, facts = [], []
    for r in both.itertuples():
        perm, f = permutation_for(r.SMILES, r.canonical_smiles)
        # (3b) the element check again, this time against the REAL qm_atoms.csv rows rather than
        # a re-parse of the SMILES -- so a silent change to the QM table would be caught here.
        orca_elements = atom_elements[r.Molecule_Name]
        assert len(orca_elements) == f["n_with_h"], (r.Molecule_Name, len(orca_elements))
        can = Chem.MolFromSmiles(r.canonical_smiles)
        for j in range(f["n_heavy"]):
            assert can.GetAtomWithIdx(j).GetSymbol() == orca_elements[int(perm[j])], \
                f"{r.Molecule_Name}: qm_atoms.csv element mismatch at chemprop index {j}"
        rows.append({"Molecule_Name": r.Molecule_Name, "set": r.set, **f,
                     "perm": " ".join(str(int(v)) for v in perm)})
        facts.append(f)

    table = pd.DataFrame(rows)
    table.to_csv(ORDERING_CSV, index=False)

    n = len(table)
    n_identity = int(table["is_identity"].sum())
    n_elseq = int(table["element_sequence_already_matched"].sum())
    print()
    print(f"AddHs preserves heavy-atom indices:            {n}/{n}")
    print(f"permutation is a bijection:                    {n}/{n}")
    print(f"element sequence matches under permutation:    {n}/{n}")
    print(f"  (checked against qm_atoms.csv itself:        {n}/{n})")
    print(f"chemprop bonds map into ORCA's bond set:       {n}/{n}")
    print(f"permutation is the IDENTITY on:                {n_identity}/{n} "
          f"({100 * n_identity / n:.1f}%) -- it does real work on {n - n_identity}")
    print(f"canonical_smiles already matched ORCA's element sequence on: {n_elseq}/{n}")
    assert n_elseq == NB43_ELEMENT_SEQUENCE_MATCHES, (n_elseq, NB43_ELEMENT_SEQUENCE_MATCHES)
    print(f"  cross-check against notebook 43's recorded 4,212: MATCH")

    # (5) negative control, over the first 50 non-identity compounds
    probes = table[~table["is_identity"]].head(50)["Molecule_Name"].tolist()
    lookup = both.set_index("Molecule_Name")
    caught = sum(negative_control(lookup.loc[p, "SMILES"], lookup.loc[p, "canonical_smiles"])
                 for p in probes)
    print(f"negative control (swap two entries): rejected on {caught}/{len(probes)} probes")
    assert caught >= 1, "the ordering assertions are vacuous -- a corrupted permutation passed"

    hydrogens = int(table["n_hydrogens"].sum())
    print()
    print(f"hydrogen rows DROPPED (chemprop has no explicit H): {hydrogens} of "
          f"{int(table['n_with_h'].sum())} ORCA atoms "
          f"({100 * hydrogens / table['n_with_h'].sum():.1f}%)")

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_compounds": n,
        "n_train": int((table["set"] == "train").sum()),
        "n_blind": int((table["set"] == "blind").sum()),
        "addhs_preserves_heavy_indices": n,
        "bijection": n,
        "element_sequence_under_permutation": n,
        "element_sequence_against_qm_atoms_csv": n,
        "bond_subset_under_permutation": n,
        "permutation_is_identity": n_identity,
        "permutation_does_work_on": n - n_identity,
        "canonical_smiles_element_sequence_matches": n_elseq,
        "nb43_recorded_element_sequence_matches": NB43_ELEMENT_SEQUENCE_MATCHES,
        "negative_control_probes": len(probes),
        "negative_control_rejected": int(caught),
        "orca_atom_rows_total": int(table["n_with_h"].sum()),
        "hydrogen_rows_dropped": hydrogens,
        "heavy_atom_rows_used": int(table["n_heavy"].sum()),
        "qm_atoms_sha256": sha256(ATOMS),
        "ordering_csv_sha256": sha256(ORDERING_CSV),
    }
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {ORDERING_CSV.relative_to(REPO_ROOT)} and "
          f"{REPORT_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
