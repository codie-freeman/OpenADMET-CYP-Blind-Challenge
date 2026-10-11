#!/usr/bin/env python
"""Per-atom QM feature arrays and the widened CheMeleon encoder. Notebook 52's stage 0.

Two jobs, both of which notebook 52 needs before it can train anything.

================================================================================
1. THE .npz ARRAYS
================================================================================
`chemprop train --atom-features-path p.npz` expects `np.savez(p, *V_fs)` -- one 2D array per
datapoint, `n_atoms x n_features`, in the SAME ROW ORDER as the input CSV (confirmed against
`chemprop/cli/utils/parsing.py:load_input_feats_and_descs`, which reads them back as
`arr_0 ... arr_{n-1}`). So an npz is tied to one CSV, and the train and predict npz files are
built separately.

Two columns: `q_hirshfeld` and `mayer_va`. Both are defined on every atom, which is why they were
chosen -- notebook 45's finding 5 recorded 113,596 nulls in the shell-population columns, exactly
the hydrogen count, so `mulliken_s/p/d` and `loewdin_s/p/d` would need a hydrogen-null policy
this notebook deliberately does not take on.

Values are reordered into chemprop's atom order by `scripts/52_build_atom_ordering.py`'s
permutation and restricted to heavy atoms, because chemprop does not add hydrogens.

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
CHEMELEON_CKPT = Path.home() / ".chemprop" / "chemeleon_mp.pt"

# The two per-atom columns. Both defined on every atom; see the module docstring.
ATOM_FEATURE_COLS = ["q_hirshfeld", "mayer_va"]
K_ATOM_FEATURES = len(ATOM_FEATURE_COLS)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ------------------------------------------------------------------ the .npz arrays
class AtomFeatureSource:
    """Per-compound heavy-atom QM feature matrices, already in chemprop's atom order."""

    def __init__(self) -> None:
        if not ORDERING_CSV.exists():
            raise FileNotFoundError(
                f"{ORDERING_CSV} is missing -- run scripts/52_build_atom_ordering.py first")
        order = pd.read_csv(ORDERING_CSV)
        self.n_heavy = dict(zip(order["Molecule_Name"], order["n_heavy"]))
        perms = {r.Molecule_Name: np.fromstring(r.perm, dtype=np.int64, sep=" ")
                 for r in order.itertuples()}

        atoms = pd.read_csv(ATOMS_CSV, usecols=["Molecule_Name", "atom_index", *ATOM_FEATURE_COLS])
        assert atoms[ATOM_FEATURE_COLS].notna().all().all(), \
            "a selected per-atom column carries nulls -- both were chosen because neither does"

        self.matrices: dict[str, np.ndarray] = {}
        for name, g in atoms.groupby("Molecule_Name", sort=False):
            perm = perms[name]
            # ORCA order, then permuted into chemprop order and cut to the heavy atoms.
            full = g.sort_values("atom_index")[ATOM_FEATURE_COLS].to_numpy(dtype=np.float64)
            assert len(perm) == self.n_heavy[name]
            self.matrices[name] = full[perm]
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
                "n_features": K_ATOM_FEATURES,
                "scaled_mean": [float(v) for v in stacked.mean(axis=0)],
                "scaled_std": [float(v) for v in stacked.std(axis=0)]}

    def write_npz_zeroed(self, names_in_csv_order, path: Path) -> dict:
        """An all-zeros npz of the right shape. The wiring check's control arm: same CLI, same
        widened encoder, same code path, only the VALUES differ -- so a difference in predictions
        can only come from the values having reached the model."""
        arrays = [np.zeros_like(self.matrices[n]) for n in names_in_csv_order]
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(path, *arrays)
        return {"path": str(path), "n_datapoints": len(arrays), "zeroed": True}


# ------------------------------------------------------------------ the widened encoder
def build_widened_chemeleon(k: int, path: Path) -> dict:
    """Write a local foundation file: CheMeleon with `k` zero-initialised extra atom-feature
    columns. Asserts at-init numerical identity with CheMeleon before saving."""
    import torch
    from rdkit import Chem

    from chemprop.data import BatchMolGraph
    from chemprop.featurizers import MultiHotAtomFeaturizer, SimpleMoleculeMolGraphFeaturizer
    from chemprop.models import MPNN
    from chemprop.models.utils import save_model
    from chemprop.nn import BondMessagePassing, MeanAggregation, RegressionFFN

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


if __name__ == "__main__":
    import json

    src = AtomFeatureSource()
    print(f"per-atom matrices: {len(src.matrices)} compounds, "
          f"{sum(m.shape[0] for m in src.matrices.values())} heavy atoms, "
          f"{K_ATOM_FEATURES} features {ATOM_FEATURE_COLS}")
    info = build_widened_chemeleon(K_ATOM_FEATURES, OUT / "chemeleon_widened.pt")
    print(json.dumps(info, indent=2))
