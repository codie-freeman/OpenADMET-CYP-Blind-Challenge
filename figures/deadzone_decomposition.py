"""How much of the dead-zone gain is the loss function alone. Read from notebook 37.

Every number is read from ``outputs/37_butina_5x5_complete/mae_only_decomposition.csv``, which
notebook 37 computed at 25 folds; nothing is recomputed.
"""

from __future__ import annotations

from figures.data import load
from figures.style.base import Target, figure, figure_title, save
from figures.style.palette import ISOFORM_COLOURS, ISOFORMS, REFERENCE, TEXT

SLUG = "deadzone_decomposition"
SIZES = {"thesis": (6.3, 2.8), "slides": (9.0, 4.4), "web": (8.0, 3.4)}
FOLD_SET = "25-fold"

CAPTION = """
The dead-zone arm changes two things at once against the plain baseline — the training loss (mean
squared to mean absolute) and the target (the point estimate to that estimate clipped into its
credible band). A third arm, `mae_only`, changes **only** the loss and keeps the literal label, so
the gain splits cleanly. Bars show the share attributable to the loss switch alone, at 25
cluster-disjoint folds; the remainder is what clipping adds on top.

**The split is strongly isoform-specific and it inverts across the set.** On CYP3A4 the loss switch
accounts for 5 per cent and clipping does essentially all the work; on CYP2D6 it accounts for 66 per
cent and clipping adds the smaller part. CYP1A2 and CYP2C9 sit between, at 36 and 14 per cent.

That matters for what the dead-zone result can be claimed to show. `mae_only` carries **none** of
the fold-crosstalk exposure the clipped target does — its target is the curated label, not an
out-of-fold model's prediction — so wherever the loss-alone share is large, the gain is on firmer
ground. CYP2D6's is the largest share and therefore the best-supported, and it is also the isoform
whose out-of-fold rankings the board has contradicted twice.
"""


def variant(target: Target) -> str:
    return f"single-panel ({FOLD_SET})"


def table():
    t = load("37_butina_5x5_complete", "mae_only_decomposition.csv")
    t = t[t.fold_set == FOLD_SET].set_index("isoform").loc[list(ISOFORMS)]
    assert t.deadzone_beats_plain.all(), "dead-zone does not beat plain on some isoform"
    return t


def draw(target: Target):
    t = table()
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        ax = fig.add_subplot(111)
        x = range(len(t))
        share = t.loss_alone_share.to_numpy() * 100
        ax.bar(x, share, color=[ISOFORM_COLOURS[i] for i in t.index], width=0.6, zorder=3)
        ax.bar(x, 100 - share, bottom=share, color=REFERENCE, alpha=0.55,
               width=0.6, zorder=2)

        for i, (iso, r) in enumerate(t.iterrows()):
            # Below ~12% the bar is too short to hold its own label, so it goes just above it.
            inside = share[i] > 12
            ax.text(i, share[i] / 2 if inside else share[i] + 3,
                    f"{share[i]:.0f}%", ha="center",
                    va="center" if inside else "bottom",
                    color="white" if inside else TEXT,
                    fontweight="bold", fontsize="small")
            ax.text(i, share[i] + (100 - share[i]) / 2, f"{100 - share[i]:.0f}%",
                    ha="center", va="center", color=TEXT, fontsize="small")
            ax.text(i, 104, f"total gain {r.total_gain_plain_minus_deadzone:.3f}",
                    ha="center", va="bottom", fontsize="small", color=TEXT)

        ax.set_xticks(list(x)); ax.set_xticklabels(t.index)
        ax.set_ylim(0, 122)
        ax.set_ylabel("Share of the dead-zone gain (%)")
        ax.yaxis.grid(True, zorder=0); ax.set_axisbelow(True)
        fig.supxlabel("coloured = loss switch alone (MSE to MAE)   "
                      "grey = clipping into the credible band, on top",
                      fontsize="small", color=TEXT)
        figure_title(fig, target, "What the dead-zone gain is actually made of")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
