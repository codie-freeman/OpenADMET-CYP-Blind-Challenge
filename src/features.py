"""Shared SMILES canonicalization, InChIKey generation, fingerprints, similarity,
SALI activity-cliff scoring, structural descriptors, and scaffolds.

Centralized here per CLAUDE.md's "fingerprints/similarity/splits go through one
shared module" rule -- these are used identically across notebooks (canonicalization/
InChIKey in curation; fingerprints/similarity/SALI/descriptors/scaffolds in chemical
space exploration and beyond) rather than redefined per notebook.

Three descriptor options are provided, for three different purposes: `isoform_structural_
descriptors` is a narrow, pharmacologically-motivated 9-descriptor set built for
notebook 02's isoform-specific SAR interpretability work; `rdkit_2d_descriptors` is the
full RDKit 2D descriptor set (200+ descriptors), added in notebook 03 (Section 2) as the
broader tabular-baseline feature option matching this challenge's official baseline
models (RDKit 2D + ECFP4); `mordred_2d_descriptors` is the full Mordred 2D descriptor
set (~1600 descriptors, via the maintained `mordred-community` fork), added to replace
RDKit2D in the tabular arm of the 5x5 CV comparison, per Mauricio's suggestion. None
replaces the others -- pick per the notebook's purpose.

`frozen_encoder_embeddings` is the shared forward-pass/batching loop for extracting
molecule embeddings from an already-loaded, already-frozen chemprop encoder + agg pair
-- used by both `chemeleon_embeddings` (CheMeleon's own checkpoint) and
`scripts/train_log2fc_encoder.py` (a freshly-trained-then-frozen encoder, notebook 03b),
so that loop isn't duplicated between the two.

`assign_screen_split`, added for notebook 04a's cheap single-fold baseline screen, is
the other half of the "splits go through one shared module" rule -- it's called
identically (same args) from both `notebooks/04a_baseline_screen.ipynb` (tabular
configs) and `scripts/run_baseline_screen_chemprop.py` (Chemprop configs), so both
processes see byte-identical train/inner-val/test assignments without needing to freeze
an intermediate split file or coordinate launch order.

`butina_clusters`, added for notebook 05b's cluster-aware CV split, is the clustering
half of the same rule -- Tanimoto-similarity-threshold clustering (via RDKit's Butina
algorithm), distinct from `bemis_murcko_scaffold`'s identity-based grouping.

`qm_atom_aggregates`, added for notebook 45's label-free redundancy screen, is the only
function here that does not start from SMILES: it aggregates `scripts/parse_orca_output.py`'s
per-atom and per-bond ORCA tables to one row per compound. It is centralized for the same
reason as the fingerprints above -- notebook 45 screens the block and a later modelling arm
would consume it, and the two must enumerate identical descriptors. Note it needs no atom-index
join and so is immune to the hazard notebook 43 documented, that `canonical_smiles` reproduces
ORCA's atom order on only 4,212 of 5,655 compounds; it groups on `Molecule_Name` and reads
`element` from the ORCA table itself.
"""

import os

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina
from sklearn.model_selection import train_test_split


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


# Reference descriptor-name ordering for `rdkit_2d_descriptors`, computed once from a
# trivial valid molecule at module load (confirmed identical across a range of
# molecules -- including disconnected/ionic species -- during development of that
# function, not assumed) and reused as the NaN-fallback key set for unparseable SMILES,
# so the failure case always returns the same shape as the success case.
_RDKIT_2D_DESCRIPTOR_NAMES = tuple(Descriptors.CalcMolDescriptors(Chem.MolFromSmiles("C")).keys())


def rdkit_2d_descriptors(smiles: str) -> dict:
    """All RDKit-computed 2D descriptors (via `Descriptors.CalcMolDescriptors`) for
    `smiles`, as a dict keyed by descriptor name.

    This is the full RDKit 2D descriptor set (200+ descriptors), used here as the
    broader companion to `isoform_structural_descriptors` -- see that function's
    docstring for the narrower, pharmacologically-motivated alternative. Matches the
    feature recipe used by this challenge's official baseline models (RDKit 2D
    descriptors + ECFP4), per OpenADMET's own description of XGB-baseline/LGBM-baseline
    (confirmed via the project's OpenADMET Discord thread).

    Returns a dict of NaN for every field if the SMILES fails to parse. Field set is
    determined by the installed RDKit version -- callers needing a stable, versioned
    column set should pin rdkit's version and log it, not assume this list is fixed
    across environments.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {name: np.nan for name in _RDKIT_2D_DESCRIPTOR_NAMES}
    return Descriptors.CalcMolDescriptors(mol)


def mordred_2d_descriptors(smiles_list: list[str], nproc: int | None = None) -> pd.DataFrame:
    """Full Mordred 2D descriptor set (~1600 descriptors) for a list of SMILES.

    Uses the actively-maintained `mordred-community` fork (PyPI package
    `mordredcommunity`, already pinned in `environment.yml`; import name is still
    `mordred`) -- NOT the original PyPI `mordred` package, which is unmaintained and
    known to break on recent RDKit. Compatibility with this project's pinned RDKit
    2026.3.3 was confirmed separately before this function was written: clean install,
    clean run on a small sample of `train_inhibition_curated.csv`, no exceptions, no
    Python warnings, no RDKit stderr output.

    Mordred computes some descriptors (e.g. high-topological-distance autocorrelation
    descriptors like `AATS8*`) that are structurally undefined for small/simple
    molecules, and returns its own `Missing`/`Error` wrapper objects for those cells
    instead of raising. This function coerces every such wrapper to NaN via
    `pandas.to_numeric(errors="coerce")` per column, then casts every column to
    float64 (some Mordred descriptors are boolean-valued, which would otherwise leave
    mixed dtypes across columns) -- it does NOT drop or impute anything; NaN cells are
    left as NaN in the returned DataFrame for the caller to inspect and decide how to
    handle. SMILES that fail to parse produce an all-NaN row
    (same descriptor columns) rather than being dropped, matching `rdkit_2d_descriptors`'s
    contract.

    Descriptor computation is farmed out to `nproc` worker processes internally by
    Mordred (default: `os.cpu_count()`) for speed -- row order is preserved exactly
    under multiprocessing, confirmed via direct comparison of `nproc=1` vs
    `nproc=os.cpu_count()` output on a 200-compound sample during development.

    Args:
        smiles_list: SMILES strings (parsed independently here -- need not already be
            canonical).
        nproc: Mordred worker process count. Defaults to `os.cpu_count()`.

    Returns:
        DataFrame with `len(smiles_list)` rows in the same order, one column per
        Mordred 2D descriptor (column names are Mordred's own descriptor names).
    """
    from mordred import Calculator
    from mordred import descriptors as mordred_descriptors

    calc = Calculator(mordred_descriptors, ignore_3D=True)
    descriptor_names = [str(d) for d in calc.descriptors]

    mols = [Chem.MolFromSmiles(s) for s in smiles_list]
    valid_idx = [i for i, m in enumerate(mols) if m is not None]
    valid_mols = [mols[i] for i in valid_idx]

    if nproc is None:
        nproc = os.cpu_count() or 1

    out = pd.DataFrame(np.nan, index=range(len(smiles_list)), columns=descriptor_names)
    if valid_mols:
        computed = calc.pandas(valid_mols, nproc=nproc, quiet=True)
        computed = computed[descriptor_names].apply(pd.to_numeric, errors="coerce").astype(np.float64)
        out.iloc[valid_idx] = computed.to_numpy()
    return out


def chemeleon_embeddings(
    smiles_list: list[str],
    checkpoint_path: str | None = None,
    batch_size: int = 128,
    device: str = "cpu",
) -> np.ndarray:
    """Mean-pooled molecule embeddings from the pretrained CheMeleon foundation-model
    encoder (chemprop `BondMessagePassing`, Zenodo record 15460715).

    A single frozen forward pass -- the checkpoint's weights are loaded and never
    updated here, so this is feature extraction, not training, matching the "fast
    forward pass" framing it's used under in notebook 03. `checkpoint_path` defaults to
    `~/.chemprop/chemeleon_mp.pt`, chemprop's own foundation-model cache location
    (downloaded automatically by chemprop/openadmet tooling on first use, or manually
    from https://zenodo.org/records/15460715/files/chemeleon_mp.pt); raises
    `FileNotFoundError` rather than downloading it here, so this function has no
    network dependency. `device` defaults to CPU rather than MPS for exact
    run-to-run determinism, since this output is meant to be frozen to disk once and
    never recomputed.

    torch/chemprop are imported lazily inside this function so notebooks that only
    need canonicalization/fingerprints/descriptors aren't forced to load them.

    Caller must set `OMP_NUM_THREADS=1` (via `os.environ`, before rdkit/torch/chemprop
    are first imported anywhere in the process -- e.g. at the very top of a notebook)
    to avoid a segfault from rdkit and torch each bundling their own conflicting
    OpenMP runtime on macOS; confirmed via direct reproduction during development of
    this function (crashes intermittently, sometimes silently, without it). This
    function does not set it itself since it must be set before those modules'
    top-level imports elsewhere in the process, which may already have happened by the
    time this function is called.

    Every SMILES in `smiles_list` must already be RDKit-parseable (callers should
    pre-filter with `canonicalize_smiles`) -- chemprop's own SMILES parser raises on a
    failure rather than returning None, so there is no per-row failure mode to handle
    here. Returns a `(len(smiles_list), 2048)` float32 array in the same row order as
    `smiles_list`; order is preserved by disabling dataloader shuffling and forcing
    `drop_last=False` (chemprop's own default can silently drop the final row when
    `len(dataset) % batch_size == 1`).
    """
    from pathlib import Path

    import torch
    from chemprop import data as cp_data
    from chemprop import nn as cp_nn

    if checkpoint_path is None:
        checkpoint_path = Path.home() / ".chemprop" / "chemeleon_mp.pt"
    checkpoint_path = Path(checkpoint_path)
    if not checkpoint_path.exists():
        raise FileNotFoundError(
            f"CheMeleon checkpoint not found at {checkpoint_path} -- download it from "
            "https://zenodo.org/records/15460715/files/chemeleon_mp.pt first."
        )

    ckpt = torch.load(checkpoint_path, weights_only=True, map_location="cpu")
    encoder = cp_nn.BondMessagePassing(**ckpt["hyper_parameters"])
    encoder.load_state_dict(ckpt["state_dict"])
    encoder = encoder.to(device).eval()
    agg = cp_nn.MeanAggregation().to(device)

    datapoints = [cp_data.MoleculeDatapoint.from_smi(smi) for smi in smiles_list]
    dataset = cp_data.MoleculeDataset(datapoints)
    loader = cp_data.build_dataloader(
        dataset, batch_size=batch_size, shuffle=False, num_workers=0, drop_last=False
    )

    embeddings = []
    with torch.no_grad():
        for batch in loader:
            bmg = batch.bmg
            bmg.to(device)
            h = encoder(bmg)
            embeddings.append(agg(h, bmg.batch).cpu().numpy())
    return np.concatenate(embeddings, axis=0).astype(np.float32)


def frozen_encoder_embeddings(
    smiles_list: list[str],
    encoder,
    agg,
    batch_size: int = 128,
    device: str = "cpu",
) -> np.ndarray:
    """Mean-pooled molecule embeddings from an already-loaded, already-frozen chemprop
    message-passing encoder (`encoder`) + aggregation module (`agg`).

    This is the shared batching/forward-pass loop behind `chemeleon_embeddings` (which
    loads CheMeleon's specific Zenodo checkpoint) and `scripts/train_log2fc_encoder.py`
    (which extracts embeddings from a freshly-trained-then-frozen encoder) -- added here
    as a standalone function rather than duplicating the loop in the training script,
    per the shared-module rule. `chemeleon_embeddings` is left as-is (not refactored to
    call this) since it's already a shipped, verified artifact-producing function --
    this is purely additive.

    Callers are responsible for putting `encoder`/`agg` in eval mode and on `device`
    before calling this -- this function does not change their mode or device
    placement, only runs the forward pass under `torch.no_grad()`.

    Same row-order and parse-failure contract as `chemeleon_embeddings`: dataloader
    shuffling is disabled and `drop_last=False` is forced, so no row is silently
    reordered or dropped; every SMILES must already be RDKit-parseable (callers should
    pre-filter with `canonicalize_smiles`). Returns a `(len(smiles_list), d_h)` float32
    array in the same row order as `smiles_list`, where `d_h` is `encoder`'s own output
    dimension.
    """
    import torch
    from chemprop import data as cp_data

    datapoints = [cp_data.MoleculeDatapoint.from_smi(smi) for smi in smiles_list]
    dataset = cp_data.MoleculeDataset(datapoints)
    loader = cp_data.build_dataloader(
        dataset, batch_size=batch_size, shuffle=False, num_workers=0, drop_last=False
    )

    embeddings = []
    with torch.no_grad():
        for batch in loader:
            bmg = batch.bmg
            bmg.to(device)
            h = encoder(bmg)
            embeddings.append(agg(h, bmg.batch).cpu().numpy())
    return np.concatenate(embeddings, axis=0).astype(np.float32)


def butina_clusters(fps: list, cutoff: float = 0.4) -> list[tuple[int, ...]]:
    """Cluster `fps` via RDKit's Butina algorithm (Butina, J. Chem. Inf. Comput. Sci.
    1999) at Tanimoto similarity `cutoff` -- pairs within `1 - cutoff` Tanimoto
    distance can land in the same cluster. Added for notebook 05b's cluster-aware CV
    split (finer-grained than the Bemis-Murcko scaffold split notebook 02 already ruled
    out) -- notebook 02's own SALI/nearest-neighbour analysis never needed clustering
    itself, only the raw pairwise/NN similarity functions above.

    `fps` must contain no `None` entries -- callers pre-filter (e.g. via
    `ecfp4_fingerprints`'s own parse-failure contract) before calling this. O(N^2)
    pairwise distance computation, same scale caveat as `pairwise_tanimoto_matrix`.

    Returns one tuple of positional indices into `fps` per cluster, in Butina's own
    output order (each tuple's first element is that cluster's centroid, the point
    with the most neighbours within `cutoff`); singleton clusters (no other point
    within `cutoff`) are single-element tuples. Every index in `range(len(fps))`
    appears in exactly one cluster.
    """
    n = len(fps)
    distances = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        distances.extend(1.0 - s for s in sims)
    return Butina.ClusterData(distances, n, 1.0 - cutoff, isDistData=True)


def bemis_murcko_scaffold(smiles: str) -> str | None:
    """Return the Bemis-Murcko scaffold (rings + linkers, atom/bond types retained) for
    `smiles` as a canonical SMILES string, or None if the SMILES fails to parse."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(scaffold)


def assign_screen_split(
    cv_folds: pd.DataFrame,
    repeat_col: str = "repeat_0",
    test_fold: int = 0,
    val_fraction: float = 0.15,
    seed: int = 42,
) -> pd.DataFrame:
    """Compound-level train/inner-val/test assignment for notebook 04a's cheap,
    single-fold baseline screen -- NOT the 25-fold (5x5) comparison's own splitting.

    `test_fold` of `repeat_col` (default: fold 0 of `repeat_0`) is held out as
    `"screen_test"`; the remaining compounds (the pooled training side for this screen)
    are further split `1 - val_fraction` / `val_fraction` into `"screen_inner_train"` /
    `"screen_inner_val"` via `sklearn.train_test_split(random_state=seed)`, for
    early-stopping-capable configurations (XGBoost, LightGBM, both Chemprop variants)
    to validate against without ever touching `screen_test`. Configurations with no
    early-stopping concept (RF, TabICLv2) are expected to fit on
    `screen_inner_train` + `screen_inner_val` pooled back together -- this function
    doesn't do that pooling itself, callers select the rows they need by
    `screen_split`.

    Deliberately named `screen_*` rather than plain `"train"`/`"val"`/"test"` -- these
    values must never be confused with the `split` column already used in the tabular
    feature files (official train/test, i.e. which compounds are in the real blinded
    leaderboard test set) or with that real test set itself. This screen's `"screen_
    test"` compounds are ordinary held-out *training* compounds; the official blinded
    test set is never touched anywhere in this function.

    A pure, deterministic function of `cv_folds` + these fixed arguments (no I/O, no
    hidden state) -- calling it with the same arguments from two different processes
    (as `notebooks/04a_baseline_screen.ipynb` and
    `scripts/run_baseline_screen_chemprop.py` both do) yields byte-identical output, so
    the split never needs to be frozen to disk or coordinated by launch order.

    Args:
        cv_folds: `data/folds/cv_folds.csv`, loaded as-is -- must contain
            `Molecule_Name`, `inchikey`, and `repeat_col`.
        repeat_col: Which of the 5 repeat columns to draw `test_fold` from.
        test_fold: Which fold (0-4) of `repeat_col` is this screen's held-out test set.
        val_fraction: Fraction of the non-test compounds set aside as
            `"screen_inner_val"`.
        seed: Passed to `train_test_split` for the inner train/val split. Logged by
            callers, not here (this function has no logging of its own).

    Returns:
        DataFrame with `Molecule_Name`, `inchikey`, `screen_split` (one of
        `"screen_test"`, `"screen_inner_train"`, `"screen_inner_val"`), one row per
        compound in `cv_folds`.
    """
    test_mask = cv_folds[repeat_col] == test_fold
    test_df = cv_folds.loc[test_mask, ["Molecule_Name", "inchikey"]].copy()
    test_df["screen_split"] = "screen_test"

    pool_df = cv_folds.loc[~test_mask, ["Molecule_Name", "inchikey"]].copy()
    inner_train_ik, inner_val_ik = train_test_split(
        pool_df["inchikey"].to_numpy(), test_size=val_fraction, random_state=seed
    )
    inner_val_set = set(inner_val_ik)
    pool_df["screen_split"] = np.where(
        pool_df["inchikey"].isin(inner_val_set), "screen_inner_val", "screen_inner_train"
    )

    return pd.concat([test_df, pool_df], ignore_index=True)


def assign_final_submission_split(
    curated: pd.DataFrame,
    val_fraction: float = 0.15,
    seed: int = 42,
) -> pd.DataFrame:
    """Compound-level train/val assignment for the final activity-track submission
    model (`scripts/train_final_submission_multitask.py`) -- an internal validation
    slice used purely for Chemprop's `--patience`-based early stopping when training on
    the FULL labeled set (all 4,905 compounds in `train_inhibition_curated.csv`), not a
    CV fold.

    Entirely independent of `data/folds/cv_folds.csv` and `assign_screen_split` above:
    this function only ever reads `curated`, never touches `cv_folds.csv`, and returns
    a `final_submission_split` column (values `"final_train"`/`"final_val"`) --
    deliberately distinct naming from `screen_split`'s `"screen_inner_train"`/
    `"screen_inner_val"` so the two can never be confused or accidentally merged.

    `val_fraction=0.15` matches this project's own precedent for the same purpose
    (`assign_screen_split`'s own `VAL_FRACTION`, used for the same "hold out a slice for
    Chemprop's --patience early stopping" role) -- not an independently invented cutoff.

    A pure, deterministic function of `curated` + these fixed arguments (no I/O) --
    the caller writes the returned split to disk itself.

    Args:
        curated: e.g. `data/processed/train_inhibition_curated.csv`, loaded as-is --
            must contain `Molecule_Name` and `inchikey`.
        val_fraction: Fraction of compounds set aside as `"final_val"`.
        seed: Passed to `train_test_split`. Logged by callers, not here.

    Returns:
        DataFrame with `Molecule_Name`, `inchikey`, `final_submission_split` (one of
        `"final_train"`, `"final_val"`), one row per compound in `curated`.
    """
    out = curated[["Molecule_Name", "inchikey"]].copy()
    train_ik, val_ik = train_test_split(
        out["inchikey"].to_numpy(), test_size=val_fraction, random_state=seed
    )
    val_set = set(val_ik)
    out["final_submission_split"] = np.where(
        out["inchikey"].isin(val_set), "final_val", "final_train"
    )
    return out


# Aggregation specification for `qm_atom_aggregates`, kept module-level so notebook 45's
# screen and any later modelling arm enumerate the same candidate block without either
# re-deriving it. Group labels match notebook 45's own G1-G8 table.
QM_AGG_CHARGE_SCHEMES = ("q_mulliken", "q_loewdin", "q_hirshfeld")
QM_AGG_CORE_ELEMENTS = ("C", "H", "N", "O")
QM_AGG_SECONDARY_ELEMENTS = ("Cl", "S", "F")
QM_AGG_BOND_PAIRS = (("C", "C"), ("C", "H"), ("C", "N"), ("C", "O"))
QM_AGG_HIRSHFELD_ANION_CUT = -0.25
QM_AGG_BOND_ORDER_CUT = 1.5


def qm_atom_aggregates(
    atoms: pd.DataFrame,
    bonds: pd.DataFrame,
    key: str = "Molecule_Name",
) -> pd.DataFrame:
    """Aggregate ORCA per-atom and per-bond output to one row per compound.

    Added for notebook 45's label-free redundancy screen on the atom-level QM block.
    `atoms` is `outputs/qm_parsed/qm_atoms.csv` and `bonds` is `qm_bonds.csv`, both
    produced by `scripts/parse_orca_output.py`. Lives here rather than in that notebook
    so a later modelling arm reuses one definition -- descriptor drift between a screen
    and the arm it motivates is exactly what this module exists to prevent. Distinct
    from `shape_descriptors()` in `scripts/44_build_qm_feature_blocks.py`, which
    aggregates the same table's *coordinates* for notebook 44's molecule-level 3D block
    and never touches charges; that function is notebook 44's and is left alone.

    Returns 54 descriptors in eight groups. Three deliberate design choices, each of
    which a caller may need to know about:

    1. **Intensive statistics only; no sums.** A sum over atoms is a size descriptor in
       a chemistry costume, and pIC50 is intensive (the QTAIM warning,
       10.1039/d4dd00057a). The two exceptions are named `n_atoms_*`/`n_bonds_*` and
       are *deliberately* extensive -- notebook 45 uses them as positive controls for
       its size-confound check, so that a clean pass means something.
    2. **Global mean charge is excluded by construction**, because it is exactly
       `formal_charge / n_atoms`: the schemes sum to formal charge to 8e-6 (Mulliken,
       Loewdin) and 1.1e-3 (Hirshfeld, integration noise), `formal_charge` is an exact
       duplicate of a free descriptor, and `n_atoms` is a size proxy. Per-element mean
       charge is kept -- it is not charge-constrained -- as is `mayer_va` mean.
    3. **Shell populations are computed over non-hydrogen atoms only.** ORCA reports no
       p or d population for hydrogen, so `mulliken_p`/`mulliken_d`/`loewdin_p`/
       `loewdin_d` are null on exactly the 113,596 H rows (verified both ways: every
       null is an H, every H is null). That is structural, not missing data, so the
       restriction is principled rather than imputation.

    Element-conditioned columns are **NaN for compounds not containing that element**,
    never sentinel-filled. Notebook 45 screens each on its present-subset for a stated
    reason: with a sentinel standing in for most of the population, a model predicting
    the descriptor mostly learns to predict the sentinel from "does this molecule
    contain sulfur", which 2D descriptors encode perfectly, so the descriptor looks
    recoverable while its informative subset was never tested. Callers building a
    modelling matrix must decide a sentinel-and-presence-flag policy themselves.

    Bond descriptors use `is_rdkit_bond == True` rows only. That flag separates
    something physical rather than drawing an arbitrary line: notebook 43 found the
    Mayer bond-order distribution cleanly bimodal with an empty gap (real bonds
    0.604-3.125, weak contacts -0.157-0.293, zero of 269,256 pairs between), and ORCA's
    own 0.1 print threshold sits inside the weak-contact population rather than in the
    gap. Pair-conditioned means cover the four pairs above 90% population coverage.

    Spread statistics use **ddof=0**. A molecule's atoms are a complete enumeration
    rather than a sample from a larger population, and ddof=1 would additionally return
    NaN for any single-atom group. Stated because it differs from the pandas default.

    No row is dropped and no compound is invented: the returned frame is indexed by the
    unique `key` values of `atoms`, so a compound present in `atoms` but absent from
    `bonds` keeps its atom-derived columns and carries NaN in the bond-derived ones.
    """
    missing_atom_cols = {
        "element", "mayer_na", "mayer_va", "mulliken_p", "mulliken_d",
        "loewdin_p", "loewdin_d", key, *QM_AGG_CHARGE_SCHEMES,
    } - set(atoms.columns)
    if missing_atom_cols:
        raise ValueError(f"`atoms` is missing columns: {sorted(missing_atom_cols)}")
    missing_bond_cols = {
        key, "element_i", "element_j", "mayer_bond_order", "is_rdkit_bond",
    } - set(bonds.columns)
    if missing_bond_cols:
        raise ValueError(f"`bonds` is missing columns: {sorted(missing_bond_cols)}")

    out = pd.DataFrame(index=pd.Index(atoms[key].unique(), name=key))
    grp = atoms.groupby(key, sort=False)

    # G1 (15): global charge extremes and mean absolute charge, per scheme. The SIGNED
    # mean is excluded -- see docstring (2) -- but mean|q| is not charge-constrained, so
    # that argument does not reach it, and it is the one column of this family with
    # measured predictive value inside the QM block (notebook 44's post-hoc probe).
    for scheme in QM_AGG_CHARGE_SCHEMES:
        short = scheme.removeprefix("q_")
        g = grp[scheme]
        lo, hi = g.min(), g.max()
        out[f"{short}_min"] = lo
        out[f"{short}_max"] = hi
        out[f"{short}_range"] = hi - lo
        out[f"{short}_std"] = g.std(ddof=0)
        out[f"{short}_abs_mean"] = atoms.assign(_a=atoms[scheme].abs()).groupby(key, sort=False)["_a"].mean()

    # G2 (4): Mayer atomic valence. Not charge-constrained, so the mean is meaningful.
    g = grp["mayer_va"]
    out["mayer_va_min"] = g.min()
    out["mayer_va_max"] = g.max()
    out["mayer_va_mean"] = g.mean()
    out["mayer_va_std"] = g.std(ddof=0)

    # G3 (2): Mayer atomic population. Suspected element-identity proxy; kept so
    # notebook 45's size-confound and redundancy checks can test exactly that.
    out["mayer_na_max"] = grp["mayer_na"].max()
    out["mayer_na_std"] = grp["mayer_na"].std(ddof=0)

    # G4 (12) + G5 (6): element-conditioned Hirshfeld charge. NaN where absent.
    for element in (*QM_AGG_CORE_ELEMENTS, *QM_AGG_SECONDARY_ELEMENTS):
        sub = atoms.loc[atoms["element"] == element].groupby(key, sort=False)["q_hirshfeld"]
        out[f"hirshfeld_{element}_min"] = sub.min()
        out[f"hirshfeld_{element}_max"] = sub.max()
        if element in QM_AGG_CORE_ELEMENTS:
            out[f"hirshfeld_{element}_mean"] = sub.mean()

    # G6 (5): shell populations, non-hydrogen atoms only -- see docstring (3).
    heavy = atoms.loc[atoms["element"] != "H"].groupby(key, sort=False)
    out["mulliken_p_mean"] = heavy["mulliken_p"].mean()
    out["mulliken_p_std"] = heavy["mulliken_p"].std(ddof=0)
    out["mulliken_d_mean"] = heavy["mulliken_d"].mean()
    out["loewdin_p_mean"] = heavy["loewdin_p"].mean()
    out["loewdin_d_mean"] = heavy["loewdin_d"].mean()

    # G7 (8): Mayer bond order over real bonds only -- see docstring.
    real = bonds.loc[bonds["is_rdkit_bond"].astype(bool)]
    bg = real.groupby(key, sort=False)["mayer_bond_order"]
    out["bond_order_min"] = bg.min()
    out["bond_order_max"] = bg.max()
    out["bond_order_mean"] = bg.mean()
    out["bond_order_std"] = bg.std(ddof=0)
    pair_key = pd.Series(
        [
            "".join(sorted((a, b)))
            for a, b in zip(real["element_i"].to_numpy(), real["element_j"].to_numpy())
        ],
        index=real.index,
    )
    for a, b in QM_AGG_BOND_PAIRS:
        want = "".join(sorted((a, b)))
        out[f"bond_order_{a}{b}_mean"] = (
            real.loc[pair_key == want].groupby(key, sort=False)["mayer_bond_order"].mean()
        )

    # G8 (2): DELIBERATELY EXTENSIVE. Notebook 45's size-confound positive controls --
    # these should fail its |rho| >= 0.80 check, and if they pass, the check is suspect.
    out[f"n_atoms_hirshfeld_lt_neg{abs(QM_AGG_HIRSHFELD_ANION_CUT):.2f}".replace(".", "p")] = (
        atoms.assign(_f=atoms["q_hirshfeld"] < QM_AGG_HIRSHFELD_ANION_CUT)
        .groupby(key, sort=False)["_f"]
        .sum()
    )
    out[f"n_bonds_order_gt_{QM_AGG_BOND_ORDER_CUT:.1f}".replace(".", "p")] = (
        real.assign(_f=real["mayer_bond_order"] > QM_AGG_BOND_ORDER_CUT)
        .groupby(key, sort=False)["_f"]
        .sum()
    )

    return out.reset_index()
