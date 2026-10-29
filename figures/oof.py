"""Out-of-fold predictions for the shipped arm, assembled from the three notebooks that hold them.

``plain`` is ``chemprop_chemeleoninit``, the arm submitted as ``10c`` and the reference every other
arm is compared against. Its 25 cluster-disjoint folds are **split across three notebooks**: repeats
0 and 1 are notebook 34's, repeat 2 is notebook 35's, repeats 3 and 4 are notebook 37's. Reading
notebook 37's ``predictions/`` alone would silently give 10 folds.

Nothing is refitted or re-scored. The prediction files are read, concatenated and joined to the
curated labels; the only arithmetic is a mean over the five repeats and a signed residual.
"""

from __future__ import annotations

from functools import lru_cache

import pandas as pd

from figures.style.base import repo_root
from figures.style.palette import ISOFORMS

ARM = "plain"
#: notebook directory -> the repeats whose prediction files it holds.
REPEAT_SOURCE = {"34_butina_5x5": (0, 1), "35_butina_repeat3": (2,),
                 "37_butina_5x5_complete": (3, 4)}


@lru_cache(maxsize=1)
def oof_long() -> pd.DataFrame:
    """One row per (isoform, compound): measured, predicted (mean over 5 repeats), residual."""
    root = repo_root()
    frames = []
    for directory, repeats in REPEAT_SOURCE.items():
        d = root / "outputs" / directory / "predictions"
        for r in repeats:
            for f in sorted(d.glob(f"{ARM}__r{r}f*.csv")):
                t = pd.read_csv(f, float_precision="round_trip")
                t["repeat"] = r
                frames.append(t)
    preds = pd.concat(frames, ignore_index=True)
    assert preds.repeat.nunique() == 5, f"{preds.repeat.nunique()} repeats"

    curated = pd.read_csv(root / "data" / "processed" / "train_inhibition_curated.csv",
                          float_precision="round_trip")
    rows = []
    for iso in ISOFORMS:
        col = f"{iso}_pIC50_direct_inhibition"
        p = (preds.groupby("inchikey")[col].mean().rename("predicted").reset_index())
        # The raw SMILES column, not canonical_smiles: notebook 43 established the two disagree
        # on 1,443 of 5,655 compounds in atom order, and anything drawing a structure wants the
        # molecule as the challenge supplied it.
        m = curated[["inchikey", "Molecule_Name", "SMILES", col]].rename(columns={col: "measured"})
        j = m.merge(p, on="inchikey", how="inner").dropna(subset=["measured"])
        j["isoform"] = iso
        j["residual"] = j.predicted - j.measured
        rows.append(j)
    out = pd.concat(rows, ignore_index=True)
    assert out.groupby("isoform").size().to_dict() == {
        "CYP1A2": 1412, "CYP2C9": 1285, "CYP2D6": 1493, "CYP3A4": 2335}, \
        out.groupby("isoform").size().to_dict()
    return out
