"""The valid (arm, isoform) joins between cross-validation and the blind board.

Shared by the parity and rank figures. Built here rather than in either, so both answer the same
question over the same points.

THE JOIN, AND WHY IT IS SMALL
-----------------------------
A point needs a submission whose column came from ONE of notebook 37's eleven cross-validated arms.
From ``outputs/55_candidate_comparison/inventory.csv``, which traced all sixty submitted columns to
models on ``Spearman == 1.0``:

* ``plain``    <- 10c            (also NB16, whose CYP1A2/CYP3A4 reproduce 10c's board values)
* ``aid``      <- NB30
* ``mae_only`` <- NB47-mae
* **eight arms were never submitted**: ``deadzone``, ``single``, and all six tree arms.

NB32-A and NB36-widened are **not** the ``deadzone`` arm. Their model is ``deadzone+aid`` -- the
dead-zone target on the AID architecture -- which is not one of the eleven, and notebook 32
deliberately never reconfirmed that combination at 25 folds. Every other submission is a
per-isoform assembly or an ensemble of several arms.

``04b`` is excluded: it is the same recipe as ``plain`` but trained in the older ``cyp-admet``
environment (torch 2.2.2 against 2.13.0), a different draw. Admitting it would put two board
readings against one cross-validated mean and double-count the arm.

CORRECTION FILTER -- AND WHY IT APPLIES TO ST-RAE BUT NOT TO RANK
-----------------------------------------------------------------
Only an **uncorrected** column is a model-versus-model read on ST-RAE, because a placement or
spread correction moves ST-RAE directly. 10c and NB30 are raw on all four isoforms; NB47-mae
corrected CYP1A2 and CYP2D6 and left CYP2C9/CYP3A4 raw. That leaves **ten** points.

Spearman is **affine-invariant**, and every correction this project deploys is affine, so for a
*ranking* comparison no filter is needed and all **twelve** (3 arms x 4 isoforms) points are valid.
That is why the two figures have different n.
"""

from __future__ import annotations

import pandas as pd

from figures.board import wide
from figures.data import load
from figures.style.palette import ISOFORMS

#: submission -> notebook 37 arm, for submissions that are exactly one arm.
ARM_OF = {"10c": "plain", "NB30": "aid", "NB47-mae": "mae_only"}

#: Isoforms whose column carried a post-hoc correction, per inventory.csv's correction_layers.
CORRECTED = {"10c": frozenset(), "NB30": frozenset(),
             "NB47-mae": frozenset({"CYP1A2", "CYP2D6"})}

EXCLUDED = {
    "04b": "same recipe as `plain` but the older cyp-admet environment; a different training draw",
    "NB16": "maps to `plain`, but its uncorrected columns reproduce 10c's board values exactly",
}


def _cv(metric_mean: str, metric_sd: str) -> pd.DataFrame:
    cv = load("37_butina_5x5_complete", "arm_summary_25fold.csv")
    assert set(cv.n_folds) == {25}, sorted(set(cv.n_folds))
    return cv.rename(columns={"isoform": "isoform", "arm": "arm",
                              metric_mean: "cv", metric_sd: "cv_sd"})[["isoform", "arm", "cv", "cv_sd"]]


def strae_pairs() -> pd.DataFrame:
    """Ten (arm, isoform) points: uncorrected board ST-RAE against 25-fold CV ST-RAE."""
    board = wide("ST-RAE")
    cv = _cv("ST-RAE_mean", "ST-RAE_std")
    rows = []
    for sub, arm in ARM_OF.items():
        for iso in ISOFORMS:
            if iso in CORRECTED[sub]:
                continue
            c = cv[(cv.arm == arm) & (cv.isoform == iso)]
            assert len(c) == 1, f"{arm}/{iso}: {len(c)} CV rows"
            rows.append({"arm": arm, "isoform": iso, "submission": sub,
                         "board": float(board.loc[sub, iso]),
                         "cv": float(c.cv.iat[0]), "cv_sd": float(c.cv_sd.iat[0])})
    df = pd.DataFrame(rows)
    assert len(df) == 10, f"expected 10 valid points, got {len(df)}"
    return df


def spearman_pairs() -> pd.DataFrame:
    """Twelve (arm, isoform) points: board Spearman against 25-fold CV Spearman.

    No correction filter: Spearman is affine-invariant, so a corrected column's ranking is the
    uncorrected model's ranking. Reuses notebook 55's own join, with 04b's four rows dropped.
    """
    df = load("55_candidate_comparison", "ordering_board_vs_cv.csv")
    df = df[df.submission != "04b"].copy()
    df = df.rename(columns={"board_spearman": "board", "cv_spearman_25fold": "cv", "cv_sd": "cv_sd"})
    assert set(df.arm) == set(ARM_OF.values()), sorted(set(df.arm))
    assert len(df) == 12, f"expected 12 points, got {len(df)}"
    # Cross-check the CV side against notebook 37's own summary rather than trusting nb55's copy.
    cv = _cv("Spearman_R_mean", "Spearman_R_std")
    for _, r in df.iterrows():
        ref = cv[(cv.arm == r.arm) & (cv.isoform == r.isoform)].cv.iat[0]
        assert abs(ref - r.cv) < 1e-6, f"{r.arm}/{r.isoform}: nb55 {r.cv} vs nb37 {ref}"
    return df[["submission", "arm", "isoform", "board", "cv", "cv_sd"]].reset_index(drop=True)
