"""Shared SMILES canonicalization, InChIKey generation, fingerprints, similarity,
SALI activity-cliff scoring, structural descriptors, and scaffolds.

Centralized here per CLAUDE.md's "fingerprints/similarity/splits go through one
shared module" rule -- these are used identically across notebooks (canonicalization/
InChIKey in curation; fingerprints/similarity/SALI/descriptors/scaffolds in chemical
space exploration and beyond) rather than redefined per notebook.
"""

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold


def canonicalize_smiles(smiles: str) -> str | None:
    """Parse `smiles` with RDKit and return its canonical form, or None if it fails to parse."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def smiles_to_inchikey(smiles: str) -> str | None:
    """Parse `smiles` with RDKit and return its InChIKey, or None if parsing fails."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    inchikey = Chem.MolToInchiKey(mol)
    return inchikey or None


def add_canonical_smiles_and_inchikey(
    df: pd.DataFrame, smiles_col: str = "SMILES"
) -> pd.DataFrame:
    """Return a copy of `df` with `canonical_smiles` and `inchikey` columns added.

    Both are derived directly from `smiles_col` via RDKit. Rows whose SMILES fail to
    parse get None in both new columns rather than being dropped -- callers are
    responsible for surfacing/handling parse failures.
    """
    out = df.copy()
    out["canonical_smiles"] = out[smiles_col].apply(canonicalize_smiles)
    out["inchikey"] = out[smiles_col].apply(smiles_to_inchikey)
    return out


def ecfp4_fingerprints(
    smiles_list: list[str], n_bits: int = 2048, include_chirality: bool = False
) -> list[DataStructs.ExplicitBitVect | None]:
    """Compute ECFP4 (Morgan, radius=2, diameter 4) fingerprints for a list of SMILES.

    `include_chirality` defaults to False -- see notebook 02 for why (a starting-point
    default, not yet tested either way on this data). Returns a list the same length and
    order as `smiles_list`; entries are None where the SMILES fails to parse -- callers
    are responsible for handling parse failures.
    """
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=n_bits, includeChirality=include_chirality
    )
    fps = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        fps.append(generator.GetFingerprint(mol) if mol is not None else None)
    return fps


def pairwise_tanimoto_matrix(fps: list) -> np.ndarray:
    """Full symmetric pairwise Tanimoto similarity matrix (NxN) for a list of fingerprints.

    Diagonal is 1.0 (self-similarity). O(N^2) -- intended for the single-set,
    few-thousand-compound scale used in notebook 02, not for arbitrarily large N.
    """
    n = len(fps)
    matrix = np.eye(n, dtype=float)
    for i in range(n - 1):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1 :])
        matrix[i, i + 1 :] = sims
        matrix[i + 1 :, i] = sims
    return matrix


def nearest_neighbor_similarity(
    query_fps: list, reference_fps: list, exclude_self_index: bool = False
) -> np.ndarray:
    """Max Tanimoto similarity from each query fingerprint to the reference set.

    If `exclude_self_index` is True, `query_fps` and `reference_fps` must be the same
    list in the same order, and each query's similarity to itself (index i vs i) is
    excluded from its own max -- used for leave-one-out nearest-neighbour within one set.
    """
    nn_sims = np.empty(len(query_fps))
    for i, query_fp in enumerate(query_fps):
        sims = np.array(DataStructs.BulkTanimotoSimilarity(query_fp, reference_fps))
        if exclude_self_index:
            sims[i] = -1.0
        nn_sims[i] = sims.max()
    return nn_sims


def compute_sali(
    fps: list, activities, similarity_threshold: float, activity_threshold: float
) -> pd.DataFrame:
    """SALI (Structure-Activity Landscape Index, Guha & Van Drie 2008) per compound pair.

    Restricted to pairs at or above `similarity_threshold` ("similar pairs" -- the
    denominator population for cliff prevalence). `activities` must be a NaN-free
    array-like aligned index-for-index with `fps` (callers pre-filter to a single
    isoform's non-null pIC50 subset -- SALI is not computed across isoforms).

    SALI_ij = |A_i - A_j| / (1 - Sim_ij); pairs with similarity == 1.0 (identical
    fingerprints) are excluded to avoid division by zero. Returns one row per similar
    pair: `i`, `j` (positional indices into fps/activities), `similarity`,
    `delta_activity`, `sali`, and `is_cliff` (delta_activity >= activity_threshold).
    """
    activities = np.asarray(activities, dtype=float)
    n = len(fps)
    rows = []
    for i in range(n - 1):
        sims = np.array(DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1 :]))
        js = np.arange(i + 1, n)
        mask = (sims >= similarity_threshold) & (sims < 1.0)
        if not mask.any():
            continue
        sim_pairs = sims[mask]
        j_pairs = js[mask]
        delta = np.abs(activities[i] - activities[j_pairs])
        sali = delta / (1 - sim_pairs)
        rows.append(
            pd.DataFrame(
                {
                    "i": i,
                    "j": j_pairs,
                    "similarity": sim_pairs,
                    "delta_activity": delta,
                    "sali": sali,
                }
            )
        )
    if not rows:
        return pd.DataFrame(
            columns=["i", "j", "similarity", "delta_activity", "sali", "is_cliff"]
        )
    out = pd.concat(rows, ignore_index=True)
    out["is_cliff"] = out["delta_activity"] >= activity_threshold
    return out


def isoform_structural_descriptors(smiles: str) -> dict:
    """RDKit-derived descriptors for the isoform-specific SAR checks in notebook 02,
    from the Kiani & Jabeen (2019) isoform-specific determinant lists.

    `logd_proxy` uses Crippen `MolLogP` -- RDKit has no native logD calculator (logD is
    pH/ionization-state-dependent; logP is not), and no dedicated logD tool is available
    in this environment, so logP is used as a stated approximation, not a true logD.

    `vsa_acc_proxy` sums `PEOE_VSA1` + `PEOE_VSA2` (surface area of the atoms with the
    most negative partial charges, i.e. O/N-type acceptor character) as an open-source
    stand-in for QikProp's `vsa_acc` (accessible surface area of H-bond-acceptor atoms),
    which RDKit does not compute under that name.

    Returns a dict of NaN for every field if the SMILES fails to parse.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "mol_wt": np.nan,
            "hba": np.nan,
            "hbd": np.nan,
            "ring_count": np.nan,
            "stereocenter_count": np.nan,
            "formal_charge": np.nan,
            "logp": np.nan,
            "logd_proxy": np.nan,
            "vsa_acc_proxy": np.nan,
        }
    logp = Crippen.MolLogP(mol)
    return {
        "mol_wt": Descriptors.MolWt(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "ring_count": rdMolDescriptors.CalcNumRings(mol),
        "stereocenter_count": rdMolDescriptors.CalcNumAtomStereoCenters(mol),
        "formal_charge": Chem.GetFormalCharge(mol),
        "logp": logp,
        "logd_proxy": logp,
        "vsa_acc_proxy": Descriptors.PEOE_VSA1(mol) + Descriptors.PEOE_VSA2(mol),
    }


def bemis_murcko_scaffold(smiles: str) -> str | None:
    """Return the Bemis-Murcko scaffold (rings + linkers, atom/bond types retained) for
    `smiles` as a canonical SMILES string, or None if the SMILES fails to parse."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(scaffold)
