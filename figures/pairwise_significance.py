"""Which arm pairs the Ash protocol actually separates. The evidence behind the CLD letters.

Method and provenance live here; the caption is ``CAPTION`` below, which ``build.py`` writes to
``build/captions.md``.

Read straight from ``outputs/37_butina_5x5_complete/figure9_pairwise_heatmaps.csv``, which notebook
37 wrote: one row per (isoform, row_arm, col_arm) with the Benjamini-Hochberg-corrected p-value and
the mean difference. **No test is recomputed here**; the figure renders a stored verdict.

Two states per cell -- separated or not -- so there is no colour scale and no colourbar. A
continuous ramp would invite reading a p-value off a shade, which is not a thing a p-value supports.
"""

from __future__ import annotations

from figures.data import load
from figures.style.base import Target, figure, figure_title, save
from figures.style.palette import ISOFORM_COLOURS, ISOFORMS, NEUTRAL, REFERENCE, TEXT, arm_label

SLUG = "pairwise_significance"
SIZES = {"thesis": (6.3, 6.4), "slides": (9.0, 5.0), "web": (8.0, 7.0)}
ALPHA = 0.05

CAPTION = """
Every pair of the eleven modelling arms, per isoform, after the Ash protocol's omnibus test and
Benjamini-Hochberg correction. A cell is **filled where the pair is significantly different** at
alpha = 0.05 and **empty where it is not**; two states only, so no colour scale and no colourbar.
Lower triangle, arms in the same order as the compact-letter-display forest.

The compact letter display reports the *verdict* of these tests and none of the evidence. This is
the evidence: **48, 51, 36 and 49 of the 55 pairs are separated** on CYP1A2, CYP2C9, CYP2D6 and
CYP3A4. Empty cells are where the protocol could not tell two arms apart at 25 folds, and they are
not scattered at random — they cluster among arms that share a representation, and CYP2D6 has by
far the most of them, which is the same isoform every other diagnostic in this work finds hardest.

Read it alongside the fold-spread figure: that shows how much the 25 folds disagree within an arm,
which is what decides whether a pair can be separated at all.
"""


def variant(target: Target) -> str:
    return "2 x 2 grid (4 isoforms, lower triangle)"


def table():
    t = load("37_butina_5x5_complete", "figure9_pairwise_heatmaps.csv")
    assert set(t.comparison) == set(ISOFORMS), sorted(set(t.comparison))
    for iso, g in t.groupby("comparison"):
        assert len(g) == 121, f"{iso}: {len(g)} cells, expected 11 x 11"
    return t


def arm_order():
    """The forest plot's own order, so the two figures can be read against each other."""
    forest = load("37_butina_5x5_complete", "figure7_forest_ST-RAE.csv")
    f = forest[forest.comparison == "MA"] if "MA" in set(forest.comparison) else forest
    return list(f.sort_values("mean").arm)


def draw(target: Target):
    t = table()
    order = arm_order()
    width, height = SIZES[target.name]
    n = len(order)

    with figure(target, width, height) as fig:
        axes = fig.subplots(2, 2)
        for ax, iso in zip(axes.ravel(), ISOFORMS):
            m = t[t.comparison == iso].pivot(index="row_arm", columns="col_arm", values="p_BH")
            m = m.loc[order, order]
            n_sig = 0
            for i in range(n):
                for j in range(i):
                    sig = float(m.iat[i, j]) < ALPHA
                    n_sig += sig
                    ax.add_patch(__import__("matplotlib").patches.Rectangle(
                        (j - 0.42, i - 0.42), 0.84, 0.84,
                        facecolor=ISOFORM_COLOURS[iso] if sig else "none",
                        edgecolor=REFERENCE, linewidth=0.4, zorder=3))
            ax.set_title(f"{iso}   {n_sig} of {n * (n - 1) // 2} separated", color=ISOFORM_COLOURS[iso])
            ax.set_xlim(-0.6, n - 1.4)
            ax.set_ylim(n - 0.4, 0.4)
            ax.set_xticks(range(n - 1))
            ax.set_xticklabels([arm_label(a) for a in order[:-1]], rotation=90, fontsize="small")
            ax.set_yticks(range(1, n))
            ax.set_yticklabels([arm_label(a) for a in order[1:]], fontsize="small")
            ax.tick_params(length=0)
            for sp in ax.spines.values():
                sp.set_visible(False)

        fig.supxlabel("filled = significantly different after Benjamini-Hochberg correction; "
                      "empty = not separated", fontsize="small", color=TEXT)
        figure_title(fig, target, "What the letters are standing in for")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
