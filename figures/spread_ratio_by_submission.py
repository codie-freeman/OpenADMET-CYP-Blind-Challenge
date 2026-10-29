"""Prediction spread against the training label SD, per submission. The compression story.

Method and provenance live here; the caption is ``CAPTION`` below.

**COVERAGE GAP, stated on the figure as well as here.** The fullest stored table of spread ratios
is notebook 36's ``all_submissions_ratio_vs_board.csv`` and it stops at NB32-A: **11 of the 17
scored submissions**. The six missing are NB36-widened, NB40-ensemble, NB47-mae, NB48-ensemble,
NB56-identify and NB56b-s050 -- which includes both deliberately widened submissions, so the figure
shows the compression but not the whole of the correction that answered it. Computing the missing
six would mean taking the SD of six submitted CSVs, which is a new quantity and belongs upstream,
not in a figure script.
"""

from __future__ import annotations

from figures.board import chronological
from figures.data import load
from figures.style.base import Target, figure, figure_title, save
from figures.style.palette import ISOFORM_COLOURS, ISOFORMS, REFERENCE, TEXT

SLUG = "spread_ratio_by_submission"
SIZES = {"thesis": (6.3, 3.4), "slides": (9.0, 4.8), "web": (8.0, 4.0)}

CAPTION = """
Each submission's prediction spread as a fraction of that isoform's training-label SD, in submission
order. A ratio of 1.0 means the predictions are as spread out as the labels; everything below it is
compression.

**Every submission is compressed on every isoform, and CYP2D6 worst.** This is the regression to the
mean the observed-against-predicted figure shows directly, measured on the submitted columns. It
matters because a compressed column is penalised twice by the scored metric: the predictions sit too
close to the mean, and the soft-threshold band they have to land inside is narrowest exactly where
the compression bites hardest.

The step up on CYP2D6 between NB13-calib and NB16 is the board-solved placement correction being
deployed for the first time, and it is the largest single improvement in this project's record.

**This figure is incomplete and the gap matters.** The stored table covers 11 of the 17 scored
submissions and stops before NB36-widened, so the two submissions that were *deliberately widened*
to address exactly this — NB36-widened and NB47-mae — are not on it. The compression is shown; the
remedy mostly is not. Regenerate from a complete table before using this as a standalone claim.
"""


def variant(target: Target) -> str:
    t = table()
    return f"single panel ({t.submission.nunique()} of 17 submissions -- coverage gap)"


def table():
    t = load("36_spread_corrected_nb32a", "all_submissions_ratio_vs_board.csv")
    order = [s for s in chronological().submission.drop_duplicates() if s in set(t.submission)]
    t = t[t.submission.isin(order)].copy()
    t["order"] = t.submission.map({s: i for i, s in enumerate(order)})
    return t.sort_values("order")


def draw(target: Target):
    t = table()
    order = list(t.drop_duplicates("submission").submission)
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        ax = fig.add_subplot(111)
        ax.axhline(1.0, color=REFERENCE, linestyle="--",
                   linewidth=target.base_size * 0.08, zorder=1)
        for iso in ISOFORMS:
            g = t[t.isoform == iso].sort_values("order")
            ax.plot(g.order, g.ratio, marker="o", color=ISOFORM_COLOURS[iso],
                    linewidth=target.base_size * 0.11,
                    markersize=target.base_size * 0.32, zorder=3)
            ax.annotate(iso, (g.order.iloc[-1] + 0.18, g.ratio.iloc[-1]),
                        va="center", ha="left", fontsize="small", color=TEXT)

        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=55, ha="right")
        ax.set_ylabel("prediction SD / training-label SD")
        ax.set_xlim(-0.5, len(order) + 1.6)
        ax.set_title("1.0 = as spread out as the labels", loc="right",
                     fontsize="small", color=TEXT)
        ax.yaxis.grid(True, zorder=0); ax.set_axisbelow(True)
        fig.supxlabel("covers 11 of 17 scored submissions: the stored table stops at NB32-A",
                      fontsize="small", color=TEXT)
        figure_title(fig, target, "Every submitted column is compressed")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
