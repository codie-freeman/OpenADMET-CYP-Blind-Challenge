"""The live half against the interim full set, on identical columns.

Method and provenance live here; the caption is ``CAPTION`` below, which ``build.py``
writes to ``build/captions.md`` alongside the renders.
"""

from __future__ import annotations

from figures.data import load
from figures.style.base import Target, figure, figure_title, save
from figures.style.palette import ISOFORM_COLOURS, TEXT

CAPTION = """
The only full-set reading in the record, against the live half, for the same submitted columns. The
interim board of 2026-09-21 scored an accidental byte-identical re-send of 27b, so this is
like-for-like: the same file, scored on half the blind set and on all of it.

**The three isoforms built by hit expansion all got easier on the full set; CYP2D6, the one isoform
excluded from it, got harder -- and by more than the other three moved combined.** Solid lines are
hit-expanded isoforms, dashed is not. Spearman moves too (not shown here; see the record), so the
halves differ in ordering difficulty as well as in placement, and no affine correction addresses
that.

This is the only generalisation estimate the record contains, and it is why the parity figure draws
no vertical error bars: a board reading is exact for a fixed file, so the uncertainty is here, not
in the measurement.

Source: ``outputs/56_cyp1a2_identification/half_vs_full.csv``, written by notebook 56 and matching
the record's own half-to-full table.
"""

SLUG = "half_vs_full_set"
SIZES = {"thesis": (6.3, 3.2), "slides": (9.0, 4.6), "web": (8.0, 3.8)}

PANELS = [("ST_half", "ST_full", "Board ST-RAE", True),
          ("R2_half", "R2_full", "Board R-squared", False)]


def variant(target: Target) -> str:
    return "two panels (ST-RAE, R2)"


def table():
    t = load("56_cyp1a2_identification", "half_vs_full.csv", index_col=0)
    t.index.name = "isoform"
    assert set(t.index) == set(ISOFORM_COLOURS), sorted(t.index)
    return t


def draw(target: Target):
    t = table()
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        axes = fig.subplots(1, 2)
        for ax, (lo_col, hi_col, label, lower_better) in zip(axes, PANELS):
            for iso, r in t.iterrows():
                colour = ISOFORM_COLOURS[iso]
                a, b = float(r[lo_col]), float(r[hi_col])
                ax.plot([0, 1], [a, b], color=colour, marker="o",
                        linewidth=target.base_size * 0.13,
                        markersize=target.base_size * 0.38,
                        linestyle="-" if bool(r.hit_expanded) else "--", zorder=3)
                ax.text(-0.04, a, f"{a:.4f}", ha="right", va="center",
                        fontsize="small", color=TEXT)
                ax.text(1.04, b, f"{b:.4f}  {iso}", ha="left", va="center",
                        fontsize="small", color=TEXT)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["live half", "interim\nfull set"])
            ax.set_xlim(-0.62, 2.25)      # room for the value labels
            ax.set_title(f"{label}{' (lower better)' if lower_better else ' (higher better)'}")
            ax.yaxis.grid(True, zorder=0)
            ax.set_axisbelow(True)
        fig.supxlabel("solid = built by hit expansion, dashed = not (CYP2D6)",
                      fontsize="small", color=TEXT)
        figure_title(fig, target, "The two halves disagree, and not in the same direction")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
