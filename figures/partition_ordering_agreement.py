"""Does the split scheme change which model wins? Almost never.

Method and provenance live here; the caption is ``CAPTION`` below. Every value is read from
notebook 37's ``ordering_agreement.csv`` and ``partition_agreement.csv``; nothing is recomputed.
"""

from __future__ import annotations

from matplotlib.lines import Line2D

from figures.data import load
from figures.style.base import Target, figure, figure_title, save
from figures.style.palette import ISOFORM_COLOURS, ISOFORMS, NEUTRAL, REFERENCE, TEXT

SLUG = "partition_ordering_agreement"
SIZES = {"thesis": (6.3, 3.8), "slides": (9.0, 5.0), "web": (8.0, 4.4)}

CAPTION = """
Whether the choice of cross-validation split changes the answer. Each point is one (arm, isoform)
cell: its ST-RAE difference from `plain` under the random partition against the same difference
under the cluster-disjoint one, with the identity as reference. A point in the lower-left or
upper-right quadrant is an arm both partitions agree about; a point in an off-diagonal quadrant is
one they disagree about in *sign*.

**Hollow points are cells where neither partition resolved a sign** — the difference is smaller
than its own standard error, so there is no direction to agree or disagree about, and they are
excluded from the count rather than scored either way. There are 15 of them against 25 resolved.

Rank correlation of the two orderings, which is a single number per isoform and so is stated here
rather than drawn, is 0.988 (CYP1A2), 0.964 (CYP2C9), 0.842 (CYP2D6) and 0.988 (CYP3A4), and the
best arm is the same under both partitions on all four.

**They agree almost completely.** Of the 25 cells where both partitions resolve a sign, **none
flips**.

That matters because an earlier five-fold run (notebook 33) found the AID arm's sign flipping on
three of four isoforms between partitions, which if real would have undermined every comparison in
this work. Raising the cluster-disjoint design from 5 folds to 25 removed the flip: it was a
low-power artefact of one cluster-to-fold draw, not a property of the split scheme. The count went
1 of 25 at ten folds, to 0 of 25 at fifteen, to 0 of 25 here.

**What this does not license** is calling either partition interchangeable with the other. The
agreement is about *ordering*; absolute levels differ, and notebook 34 found the cluster-disjoint
numbers were produced under a thread setting that changes floating-point reduction order, so
cross-partition level comparisons carry a caveat that ordering comparisons do not.
"""


def variant(target: Target) -> str:
    return "single panel (paired differences)"


def draw(target: Target):
    paired = load("37_butina_5x5_complete", "partition_agreement.csv")
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        # One panel. The ordering-agreement bars that used to sit beside it were four numbers
        # (0.988, 0.964, 0.842, 0.988) drawn as a chart -- a sentence, which the caption now
        # carries instead.
        ax = fig.add_subplot(111)
        both = paired.dropna(subset=["butina_diff", "random_diff"])
        lim = max(both.butina_diff.abs().max(), both.random_diff.abs().max()) * 1.18
        ax.plot([-lim, lim], [-lim, lim], color=REFERENCE, linestyle="--",
                linewidth=target.base_size * 0.07, zorder=1)
        ax.axhline(0, color=REFERENCE, linewidth=target.base_size * 0.05, zorder=1)
        ax.axvline(0, color=REFERENCE, linewidth=target.base_size * 0.05, zorder=1)

        resolved = paired.sign_agrees == True                      # noqa: E712
        for iso in ISOFORMS:
            g = both[(both.isoform == iso)]
            r = g[g.index.isin(paired[resolved].index)]
            u = g[~g.index.isin(paired[resolved].index)]
            ax.scatter(r.random_diff, r.butina_diff, s=target.base_size ** 2 * 0.16,
                       color=ISOFORM_COLOURS[iso], alpha=0.85, edgecolor="none", zorder=4)
            ax.scatter(u.random_diff, u.butina_diff, s=target.base_size ** 2 * 0.16,
                       facecolors="none", edgecolors=ISOFORM_COLOURS[iso],
                       linewidths=target.base_size * 0.06, zorder=3)

        n_res, n_flip = int(resolved.sum()), int((paired.sign_agrees == False).sum())  # noqa: E712
        ax.text(0.03, 0.97, f"{n_flip} of {n_res} sign-resolved cells flip",
                transform=ax.transAxes, va="top", ha="left", fontsize="small", color=TEXT)

        handles = [Line2D([], [], linestyle="none", marker="o", color=ISOFORM_COLOURS[i],
                          markersize=target.base_size * 0.40, label=i) for i in ISOFORMS]
        handles += [Line2D([], [], linestyle="none", marker="o", markerfacecolor="none",
                           markeredgecolor=NEUTRAL, markersize=target.base_size * 0.40,
                           label="sign not resolved")]
        ax.legend(handles=handles, loc="lower right", frameon=False, ncol=2,
                  handletextpad=0.3, columnspacing=1.0, labelcolor=TEXT)

        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
        ax.set_xlabel("random partition: arm ST-RAE minus plain")
        ax.set_ylabel("cluster-disjoint: arm ST-RAE minus plain")
        ax.grid(True, zorder=0); ax.set_axisbelow(True)
        figure_title(fig, target, "The split scheme does not change which model wins")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
