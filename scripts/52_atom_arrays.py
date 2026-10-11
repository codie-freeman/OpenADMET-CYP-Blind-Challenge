#!/usr/bin/env python
"""Per-atom QM feature arrays and the widened CheMeleon encoder. Notebook 52.

Two jobs, both of which notebook 52 needs before it can train anything.

================================================================================
1. THE .npz ARRAYS
================================================================================
`chemprop train --atom-features-path p.npz` expects `np.savez(p, *V_fs)` -- one 2D array per
datapoint, `n_atoms x n_features`, in the SAME ROW ORDER as the input CSV (confirmed against
`chemprop/cli/utils/parsing.py:load_input_feats_and_descs`, which reads them back as
`arr_0 ... arr_{n-1}`). So an npz is tied to one CSV, and the train and predict npz files are
built separately.

THREE columns (stage 1a; stage 0 ran the first two):

  q_hirshfeld            per-atom Hirshfeld partial charge
  mayer_va               per-atom Mayer valence
  min_attached_H_charge  per HEAVY atom, the minimum Hirshfeld charge among ITS OWN attached
                         hydrogens -- see below

The first two were chosen because both are defined on every atom, so no hydrogen-null policy is
needed: notebook 45's finding 5 recorded 113,596 nulls in the shell-population columns, exactly
the hydrogen count, which is why `mulliken_s/p/d` and `loewdin_s/p/d` are not used here.

THE THIRD COLUMN EXISTS BECAUSE CHEMPROP DROPS THE HYDROGENS. Chemprop adds no explicit
hydrogens (`add_h=False`, `keep_h=False` in its own resolved args), so 113,596 of 252,237 ORCA
atom rows have no vertex to attach to and are lost. One of the casualties is `hirshfeld_H_min`,
which notebook 45's label-free screen placed in its strong tier at full blind coverage
(axis-A R2 0.6752, i.e. poorly recoverable from the free 2D descriptors). This column recovers
that quantity INSIDE the heavy-atom graph by pushing each hydrogen's charge onto the heavy atom
it hangs off. Verified rather than asserted: the molecule-wise minimum of this column over the
heavy atoms that carry hydrogens reproduces `hirshfeld_H_min` exactly, because every hydrogen is
attached to exactly one heavy atom.

THE NO-HYDROGEN CONVENTION IS 0.0, and it is chosen on four pieces of evidence rather than
preference. 62,562 of 138,641 heavy atoms (45.1%) carry no attached hydrogen, so the choice is
not a corner case.

  1. It is CHEMPROP'S OWN convention. `MoleculeDatapoint.__post_init__` sets `NAN_TOKEN = 0` and
     fills `V_f`'s NaNs with it, so passing NaN would produce exactly 0.0 anyway.
  2. It CANNOT create a confound. The V2 atom featurizer already one-hot encodes `total_num_Hs`
     (a 6-wide block, verified directly), so the model can distinguish a heavy atom with no
     attached hydrogen from one with hydrogens WITHOUT this column. The fill can only distort a
     scale, not hide a distinction.
  3. It is physically the right filler rather than an arbitrary one. Attached-hydrogen charges
     run -0.0027 to 0.2033 with p1 at 0.0092, so 0.0 sits at the extreme low edge of the real
     range and reads as "no hydrogen-borne positive charge", which is literally true. Stated
     honestly: 7 of 113,596 hydrogens (0.006%) are at or below 0, so the fill overlaps the real
     range on a negligible number of atoms rather than on none.
  4. A FAR sentinel would destroy the column, measured not assumed. At 0.0 the between-group
     (has-H against not) share of the column's variance is 44.9%, leaving 55% genuine within-
     group variation. At -1.0 it is 99.8% and at -0.5 it is 99.3%: the column would become a
     has-H indicator in disguise, which is the failure mode notebook 45 documented for sentinel
     fills. Mean-imputation (0.0407) goes the other way, to 0.0%, erasing the distinction.

Values are reordered into chemprop's atom order by `scripts/52_build_atom_ordering.py`'s
permutation and restricted to heavy atoms.

SCALING IS DONE HERE, NOT BY CHEMPROP, and that is a correction rather than a preference.
chemprop 2.3.1 computes an atom-feature StandardScaler from the training split and logs it, but
on the `--from-foundation` path it never attaches the resulting `GraphTransform` to the loaded
message passing block (`cli/train.py` assigns `graph_transform=` only in the non-foundation
branch, lines 1529/1555). `MoleculeDataset.normalize_inputs` mutates the training dataset's
arrays in place, so TRAINING sees scaled features while the saved checkpoint carries
`nn.Identity` -- and `chemprop predict` therefore feeds the same features UNSCALED. Verified
directly: a foundation-path checkpoint's `message_passing.graph_transform` is `Identity` where an
otherwise identical non-foundation checkpoint's is a `GraphTransform` wrapping a `ScaleTransform`.
That is a silent train/predict mismatch on every blind prediction, so notebook 52 passes
`--no-atom-feature-scaling` and standardises here, from the fold's `screen_inner_train` atoms
only -- which is the population chemprop's own scaler would have used, and excludes both the
early-stopping validation set and the held-out fold.

================================================================================
2. THE WIDENED CheMeleon ENCODER
================================================================================
chemprop 2.3.1 REFUSES `--atom-features-path` outright when CheMeleon is named as the foundation
model: `cli/train.py:698-706` raises `CheMeleon does not support passing --atom-features-path`.
The reason is structural rather than arbitrary. CheMeleon's published checkpoint fixes `d_v = 72`
(the V2 atom featurizer's width) and its `W_i` is `(2048, 86) = (d_h, d_v + d_e)`, so extra atom
features widen the input that the pretrained projection can no longer accept.

The refusal is skipped for a LOCAL foundation file (`if not local_foundation:` guards it), and
the local path takes only `message_passing` and `agg` from the file. So this builder writes a
local foundation file holding CheMeleon's own weights with `k` EXTRA COLUMNS OF ZEROS inserted
into `W_i` and `W_o` at index `d_v`, which is where the featurizer appends the extra features
(`featurizers/molgraph/molecule.py:73`, `V = np.hstack((V, atom_features_extra))`) and where both
matrices' V block ends (`W_i` sees `[V || E]` per `mixins.py:9`; `W_o` sees `[V || M]` per
`base.py:finalize`).

Two consequences, both asserted in `build_widened_chemeleon` rather than argued for:

  * At initialisation the widened encoder is NUMERICALLY IDENTICAL to CheMeleon for any value of
    the extra features, because their weights are exactly zero. The pretrained encoder is
    preserved bit-for-bit; the only new parameters are the `k` columns, which training can move.
  * `agg` is `MeanAggregation`, matching what the remote CheMeleon path builds
    (`cli/train.py`: `agg = Factory.build(AggregationRegistry["mean"])`). The local path reads
    `agg` from the file, so this has to be set deliberately.

The file is named by `k` (`chemeleon_widened_k{k}.pt`) and its own `d_v` is asserted on load, so
a stale file from a different `k` can never be picked up silently.

`--multi-hot-atom-featurizer-mode V2` must still be passed by the caller: the local path also
skips CheMeleon's own check that it is set.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

OUT = REPO_ROOT / "outputs" / "52_chemprop_qm_atom_features"
ORDERING_CSV = OUT / "atom_ordering.csv"
ATOMS_CSV = REPO_ROOT / "outputs" / "qm_parsed" / "qm_atoms.csv"
CURATED = REPO_ROOT / "data" / "processed" / "train_inhibition_curated.csv"
BLIND = REPO_ROOT / "data" / "processed" / "test_blinded_curated.csv"
CHEMELEON_CKPT = Path.home() / ".chemprop" / "chemeleon_mp.pt"

# Stage 0's set, kept so that notebook's own record stays readable against what it ran.
ATOM_FEATURE_COLS_STAGE0 = ["q_hirshfeld", "mayer_va"]
# The current set. `min_attached_H_charge` is DERIVED here, not a column of qm_atoms.csv.
ATOM_FEATURE_COLS = ["q_hirshfeld", "mayer_va", "min_attached_H_charge"]
DERIVED_COLS = {"min_attached_H_charge"}
NO_HYDROGEN_FILL = 0.0
K_ATOM_FEATURES = len(ATOM_FEATURE_COLS)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def foundation_path(k: int) -> Path:
    return OUT / f"chemeleon_widened_k{k}.pt"


# ------------------------------------------------------------------ the .npz arrays
class AtomFeatureSource:
    """Per-compound heavy-atom QM feature matrices, already in chemprop's atom order."""

    def __init__(self, cols: list[str] | None = None) -> None:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")

        self.cols = list(ATOM_FEATURE_COLS if cols is None else cols)
        self.k = len(self.cols)
        if not ORDERING_CSV.exists():
            raise FileNotFoundError(
                f"{ORDERING_CSV} is missing -- run scripts/52_build_atom_ordering.py first")
        order = pd.read_csv(ORDERING_CSV)
        self.n_heavy = dict(zip(order["Molecule_Name"], order["n_heavy"]))
        perms = {r.Molecule_Name: np.fromstring(r.perm, dtype=np.int64, sep=" ")
                 for r in order.itertuples()}

        raw_cols = [c for c in self.cols if c not in DERIVED_COLS]
        need_h = bool(DERIVED_COLS & set(self.cols))
        usecols = ["Molecule_Name", "atom_index"] + sorted(set(raw_cols) | ({"q_hirshfeld"}
                                                                            if need_h else set()))
        atoms = pd.read_csv(ATOMS_CSV, usecols=usecols)
        assert atoms[[c for c in usecols if c not in ("Molecule_Name", "atom_index")]] \
            .notna().all().all(), "a source per-atom column carries nulls"
        grouped = {n: g.sort_values("atom_index") for n, g in
                   atoms.groupby("Molecule_Name", sort=False)}

        smiles = pd.concat([pd.read_csv(CURATED, usecols=["Molecule_Name", "SMILES"]),
                            pd.read_csv(BLIND, usecols=["Molecule_Name", "SMILES"])],
                           ignore_index=True)
        assert smiles["Molecule_Name"].is_unique and len(smiles) == len(perms)
        smi = dict(zip(smiles["Molecule_Name"], smiles["SMILES"]))

        self.matrices: dict[str, np.ndarray] = {}
        self.no_h_mask: dict[str, np.ndarray] = {}
        for name, perm in perms.items():
            g = grouped[name]
            nh = self.n_heavy[name]
            assert len(perm) == nh
            built = {}
            for c in raw_cols:
                built[c] = g[c].to_numpy(dtype=np.float64)       # ORCA order, all atoms
            if need_h:
                q = g["q_hirshfeld"].to_numpy(dtype=np.float64)
                # The RAW SMILES column, exactly as the QM run used, so RDKit's atom indices ARE
                # ORCA's. scripts/52_build_atom_ordering.py asserts that correspondence already.
                mol = Chem.AddHs(Chem.MolFromSmiles(smi[name]))
                assert mol.GetNumAtoms() == len(q)
                vals = np.full(len(q), NO_HYDROGEN_FILL, dtype=np.float64)
                no_h = np.zeros(nh, dtype=bool)
                for i in range(nh):
                    hs = [nb.GetIdx() for nb in mol.GetAtomWithIdx(i).GetNeighbors()
                          if nb.GetAtomicNum() == 1]
                    if hs:
                        vals[i] = min(q[j] for j in hs)
                    else:
                        no_h[i] = True
                built["min_attached_H_charge"] = vals
                self.no_h_mask[name] = no_h[perm]
            full = np.column_stack([built[c] for c in self.cols])
            self.matrices[name] = full[perm]                      # chemprop order, heavy only
        assert len(self.matrices) == len(perms), (len(self.matrices), len(perms))

    def fit_scaler(self, names) -> tuple[np.ndarray, np.ndarray]:
        """Mean/scale over the atoms of `names` only. Chemprop's own `StandardScaler` convention:
        population standard deviation, and a zero-variance column maps to scale 1.0."""
        stacked = np.concatenate([self.matrices[n] for n in names], axis=0)
        mean = stacked.mean(axis=0)
        scale = stacked.std(axis=0)
        scale[scale == 0.0] = 1.0
        return mean, scale

    def write_npz(self, names_in_csv_order, mean, scale, path: Path) -> dict:
        """Write one npz for a CSV whose rows are `names_in_csv_order`, in that order."""
        arrays = [(self.matrices[n] - mean) / scale for n in names_in_csv_order]
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(path, *arrays)
        stacked = np.concatenate(arrays, axis=0)
        return {"path": str(path), "n_datapoints": len(arrays), "n_atoms": int(len(stacked)),
                "n_features": self.k, "columns": self.cols,
                "scaled_mean": [float(v) for v in stacked.mean(axis=0)],
                "scaled_std": [float(v) for v in stacked.std(axis=0)]}

    def write_npz_zeroed(self, names_in_csv_order, path: Path) -> dict:
        """An all-zeros npz of the right shape. The control arm's input: same CLI, same widened
        encoder, same code path, only the VALUES differ -- so any difference in predictions can
        only come from the values having reached the model."""
        arrays = [np.zeros_like(self.matrices[n]) for n in names_in_csv_order]
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(path, *arrays)
        return {"path": str(path), "n_datapoints": len(arrays), "n_features": self.k,
                "columns": self.cols, "zeroed": True}

    def recovers_hirshfeld_H_min(self, aggregates_csv: Path) -> dict:
        """The check that `min_attached_H_charge` is what it claims to be: the molecule-wise
        minimum over the heavy atoms that CARRY hydrogens must reproduce notebook 45's own
        `hirshfeld_H_min`, because every hydrogen hangs off exactly one heavy atom."""
        assert "min_attached_H_charge" in self.cols
        j = self.cols.index("min_attached_H_charge")
        agg = pd.read_csv(aggregates_csv, usecols=["Molecule_Name", "hirshfeld_H_min"])
        ref = dict(zip(agg["Molecule_Name"], agg["hirshfeld_H_min"]))
        worst, n, n_skipped = 0.0, 0, 0
        for name, m in self.matrices.items():
            has_h = ~self.no_h_mask[name]
            if not has_h.any() or pd.isna(ref.get(name, np.nan)):
                n_skipped += 1
                continue
            worst = max(worst, abs(float(m[has_h, j].min()) - float(ref[name])))
            n += 1
        return {"n_compounds_checked": n, "n_skipped_no_hydrogen_at_all": n_skipped,
                "max_abs_diff_vs_nb45_hirshfeld_H_min": worst}


# ------------------------------------------------------------------ the widened encoder
def build_widened_chemeleon(k: int, path: Path | None = None) -> dict:
    """Write a local foundation file: CheMeleon with `k` zero-initialised extra atom-feature
    columns. Asserts at-init numerical identity with CheMeleon before saving."""
    import torch
    from rdkit import Chem

    from chemprop.data import BatchMolGraph
    from chemprop.featurizers import MultiHotAtomFeaturizer, SimpleMoleculeMolGraphFeaturizer
    from chemprop.models import MPNN
    from chemprop.models.utils import save_model
    from chemprop.nn import BondMessagePassing, MeanAggregation, RegressionFFN

    path = foundation_path(k) if path is None else path
    if not CHEMELEON_CKPT.exists():
        raise FileNotFoundError(
            f"{CHEMELEON_CKPT} is missing. It is the file `--from-foundation CHEMELEON` caches; "
            f"run any CheMeleon training once to populate it.")
    ck = torch.load(CHEMELEON_CKPT, weights_only=True)
    hp = dict(ck["hyper_parameters"])
    sd = {key: v.clone() for key, v in ck["state_dict"].items()}
    d_v, d_e, d_h = hp["d_v"], hp["d_e"], hp["d_h"]
    assert (d_v, d_e) == (72, 14), (d_v, d_e)
    assert tuple(sd["W_i.weight"].shape) == (d_h, d_v + d_e)
    assert tuple(sd["W_o.weight"].shape) == (d_h, d_v + d_h)
    assert hp["d_vd"] is None and hp["graph_transform"] is None

    def widen(W, at, n):
        return torch.cat([W[:, :at], torch.zeros(W.shape[0], n, dtype=W.dtype), W[:, at:]], dim=1)

    sd["W_i.weight"] = widen(sd["W_i.weight"], d_v, k)
    sd["W_o.weight"] = widen(sd["W_o.weight"], d_v, k)
    keep = ("d_e", "d_h", "bias", "depth", "dropout", "activation", "undirected", "d_vd")
    hp_wide = {key: hp[key] for key in keep}
    hp_wide["d_v"] = d_v + k
    mp = BondMessagePassing(**hp_wide)
    mp.load_state_dict(sd)

    # At-init identity, on a real molecule, with RANDOM extra feature values -- so the check
    # cannot pass merely because the features happened to be zero.
    mp0 = BondMessagePassing(**hp)
    mp0.load_state_dict(ck["state_dict"])
    f0 = SimpleMoleculeMolGraphFeaturizer(atom_featurizer=MultiHotAtomFeaturizer.v2())
    fk = SimpleMoleculeMolGraphFeaturizer(atom_featurizer=MultiHotAtomFeaturizer.v2(),
                                          extra_atom_fdim=k)
    probe = Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1")
    extra = np.random.RandomState(20261011).randn(probe.GetNumAtoms(), k)
    mp.eval()
    mp0.eval()
    with torch.no_grad():
        h0 = mp0(BatchMolGraph([f0(probe)]))
        hk = mp(BatchMolGraph([fk(probe, atom_features_extra=extra)]))
    max_abs = float((hk - h0).abs().max())
    assert max_abs == 0.0, (
        f"the widened encoder is not identical to CheMeleon at init (max abs {max_abs:.3e}) -- "
        f"the zero columns were inserted in the wrong place")

    model = MPNN(mp, MeanAggregation(), RegressionFFN(input_dim=d_h, n_tasks=1))
    path.parent.mkdir(parents=True, exist_ok=True)
    save_model(path, model, ["placeholder"])
    return {"path": str(path), "k": k, "d_v_original": d_v, "d_v_widened": d_v + k,
            "W_i_shape": list(mp.W_i.weight.shape), "W_o_shape": list(mp.W_o.weight.shape),
            "max_abs_diff_vs_chemeleon_at_init": max_abs,
            "chemeleon_ckpt_sha256": sha256(CHEMELEON_CKPT), "foundation_sha256": sha256(path)}


def load_foundation_d_v(path: Path) -> int:
    """Read a local foundation file's own `d_v`, so a stale file from a different k is caught.

    Lightning stores the message passing block's constructor arguments as a plain `AttributeDict`
    under `hyper_parameters`, not as a live module, so `d_v` is read straight out of that mapping
    rather than off a `.hparams` attribute.
    """
    import torch
    d = torch.load(path, map_location="cpu", weights_only=False)
    return int(d["hyper_parameters"]["message_passing"]["d_v"])


if __name__ == "__main__":
    import json

    src = AtomFeatureSource()
    print(f"per-atom matrices: {len(src.matrices)} compounds, "
          f"{sum(m.shape[0] for m in src.matrices.values())} heavy atoms, "
          f"{src.k} features {src.cols}")
    print(json.dumps(src.recovers_hirshfeld_H_min(
        REPO_ROOT / "outputs" / "qm_parsed" / "qm_atom_aggregates.csv"), indent=2))
    print(json.dumps(build_widened_chemeleon(K_ATOM_FEATURES), indent=2))
