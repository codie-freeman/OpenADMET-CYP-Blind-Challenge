"""Observed against predicted, out of fold. The standard figure the set did not have.

Method and provenance live here; the caption is ``CAPTION`` below, which ``build.py`` writes to
``build/captions.md``.

Predictions are the shipped arm's (``plain`` = ``chemprop_chemeleoninit``, submitted as ``10c``),
pooled over all 25 cluster-disjoint folds via ``figures.oof``. Measured values are the curated
labels. The three statistics printed in each panel are computed on the POOLED predictions, which is
a different quantity from notebook 37's mean-of-per-fold scores and is stated as such below --
pooling and averaging give different numbers and neither is wrong.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr

from figures.oof import oof_long
from figures.style.base import Target, figure, figure_title, save
from figures.style.palette import ISOFORM_COLOURS, ISOFORMS, REFERENCE, TEXT

SLUG = "observed_vs_predicted"
SIZES = {"thesis": (6.3, 4.6), "slides": (9.0, 5.4), "web": (8.0, 5.2)}
BAND = 1.0          # +/- one log unit

CAPTION = """
Measured against out-of-fold predicted pIC50 for the shipped model
(`chemprop_chemeleoninit`, the arm submitted as `10c`), pooled across all 25 cluster-disjoint folds,
one panel per isoform. The dashed line is y = x and the shaded band is +/- one log unit.

**The two useful readings are the slope and the spread.** Every panel is visibly flatter than y = x
— predictions are compressed toward the mean, which is regression to the mean and is the same
compression every submitted column shows against the training label SD. And the band tells you how
much of that matters: a prediction inside +/- 1 log unit is, for most screening purposes, the right
answer about a compound.

The printed percentage is the share of compounds predicted within one log unit: **79 / 92 / 83 /
86 per cent**, so on CYP2C9 more than nine in ten predictions are close enough to be useful for
screening even though its R-squared is 0.455.

**CYP3A4 is the best-predicted isoform and CYP2D6 the worst**, by every statistic printed, and the
gap is large: R-squared 0.611 against 0.169. CYP3A4 has nearly twice CYP2D6's training labels, and
CYP2D6 additionally uses a different assay from the other three.

The statistics are computed on the pooled predictions. Notebook 37 reports the mean of the 25
per-fold scores instead, which for this arm gives an ST-RAE of 0.86 / 0.66 / 1.00 / 0.54 and a
Spearman of 0.50 / 0.63 / 0.40 / 0.76 — pooling and averaging are different operations and the
difference is not an error in either.
"""


def variant(target: Target) -> str:
    return "2 x 2 grid (4 isoforms)"


def stats(g):
    resid = g.predicted - g.measured
    r2 = 1.0 - (resid ** 2).sum() / ((g.measured - g.measured.mean()) ** 2).sum()
    return np.abs(resid).mean(), r2, spearmanr(g.predicted, g.measured)[0]


def draw(target: Target):
    d = oof_long()
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        axes = fig.subplots(2, 2, sharex=True, sharey=True)
        lo = min(d.measured.min(), d.predicted.min()) - 0.2
        hi = max(d.measured.max(), d.predicted.max()) + 0.2
        for ax, iso in zip(axes.ravel(), ISOFORMS):
            g = d[d.isoform == iso]
            colour = ISOFORM_COLOURS[iso]
            ax.fill_between([lo, hi], [lo - BAND, hi - BAND], [lo + BAND, hi + BAND],
                            color=REFERENCE, alpha=0.25, linewidth=0, zorder=1)
            ax.plot([lo, hi], [lo, hi], color=REFERENCE, linestyle="--",
                    linewidth=target.base_size * 0.07, zorder=2)
            ax.scatter(g.predicted, g.measured, s=target.base_size ** 2 * 0.05,
                       color=colour, alpha=0.30, edgecolor="none", zorder=3)
            mae, r2, rho = stats(g)
            inside = int((np.abs(g.predicted - g.measured) <= BAND).sum())
            ax.set_title(iso, color=colour)
            ax.set_title(f"n = {len(g):,}", loc="right", fontsize="small", color=TEXT)
            # ASCII "+/-", not mathtext "$\pm$": mathtext.fontset is dejavusans, so a single
            # mathtext glyph pulls DejaVu into an otherwise all-Arial PDF. Caught by the
            # font-embedding check, which is what it is for.
            ax.text(0.03, 0.97, f"MAE {mae:.3f}   R2 {r2:.3f}\n"
                                f"Spearman {rho:.3f}\n{inside / len(g):.0%} within +/-1 log",
                    transform=ax.transAxes, va="top", ha="left", fontsize="small", color=TEXT)
            ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
            ax.grid(True, zorder=0); ax.set_axisbelow(True)

        fig.supxlabel("Out-of-fold predicted pIC50")
        fig.supylabel("Measured pIC50")
        figure_title(fig, target, "Predictions are compressed toward the mean")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
