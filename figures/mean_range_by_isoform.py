"""Every submitted prediction mean, per isoform, against identification status.

Method and provenance live here; the caption is ``CAPTION`` below, which ``build.py``
writes to ``build/captions.md`` alongside the renders.

ON THE OTHER THREE ISOFORMS, AND WHAT THIS FIGURE DOES NOT CLAIM
---------------------------------------------------------------
CYP3A4's largest internal gap is also wide (0.413, 57% of its range), but it is a
single outlier rather than two clusters -- twelve submissions between 4.618 and
4.928 plus NB13-calib alone at 5.341 -- so only CYP2D6 is bimodal in the sense the
caption uses. CYP2C9 and CYP1A2 are continuous (largest gap 47% and 38% of a range
an order of magnitude smaller).

Nor does bimodality alone account for all four of CYP2D6's equations. Two of them
(10c, NB13-calib) sit in the upper cluster and two (NB16, NB27-widened) in the
lower, but those two share a mean to six decimals and differ in spread -- so the
clusters supply the mean-shifted pair that breaks the mirror root, which is the
condition the other three isoforms failed, and the fourth equation is a spread
change at an unchanged mean.

THE STATUS SHOWN HERE IS PRE-PROBE, DELIBERATELY
------------------------------------------------
Only CYP2D6 is marked identified. That was true at 28 October 2026 and is the
diagnostic that motivated the probe, which is why the figure shows it rather than
today's state. Notebook 56b identified the other three on 29 October
(``outputs/56b_cyp1a2_corrected/solved_populations.csv``: unique admissible root
and recovered rho for CYP1A2, CYP2C9 and CYP3A4). ``docs/findings.md`` finding 2
still carries the pre-probe table and a Consequence paragraph that notebook 56b
falsified; correcting it is a separate job.

THIS FIGURE WILL BE STALE BEFORE THE CHALLENGE CLOSES
-----------------------------------------------------
``submitted_means.csv`` holds 13 submissions. There are now **seventeen** real
sends (``docs/leaderboard_submissions.md``); the file omits NB10, NB40-ensemble,
the notebook 56 probe, and notebook 56b's s = 0.5 candidate (sent 2026-10-29
14:18 UTC, score pending). Submissions close 3 November 2026.
**Regenerate this figure from a complete record after that date** -- every row
gains points, and CYP1A2's range must grow, since the probe moved it by design.
"""

from __future__ import annotations

from figures.data import load
from figures.style.base import Target, figure, figure_title, limit_ticks, save, wilkinson_bins
from figures.style.palette import ISOFORM_COLOURS, ISOFORMS, TEXT

CAPTION = """
Each of the thirteen submissions on record at 28 October 2026 -- immediately before
the identification probe was sent -- plotted as its own submitted prediction mean,
by isoform, on a shared pIC50 axis. Rows are ordered by the range that mean has
covered, annotated at the right of each rule. Identifiability is a property of this
submission history, not of the chemistry: a blind population is solvable only from
submissions that are genuinely distinct affine variants of one ordering, and a
mean-PRESERVING variant adds an equation without breaking the sign ambiguity, so
what counts is whether the mean has actually moved.

**CYP2D6's 1.478 range is two clusters, not a spread.** Eight submissions sit
within 0.006 pIC50 of each other at 3.156, five within 0.125 at 4.509-4.634, and
the 1.347 gap between them is 91% of the whole range. The split is not
chronological drift but whether a submission carries the board-solved placement
correction: NB30, sent after three corrected submissions, sits in the upper cluster
because its CYP2D6 column is raw. That jump is what gave CYP2D6 mean-shifted
variants, and so four distinct equations, where CYP1A2, CYP2C9 and CYP3A4 yielded
two apiece -- their repeated submissions being byte-identical resubmissions or
mean-preserving rescales. A deliberate +0.25 pIC50 shift on those three (notebook
56) supplied the missing third equation and identified all three (notebook 56b).
Statuses shown are pre-probe.
"""

SLUG = "mean_range_by_isoform"

SIZES = {
    "thesis": (6.3, 3.4),   # taller than the bar version: CYP2D6 stacks eight deep
    "slides": (9.0, 4.4),
    "web": (8.0, 4.0),
}

# Pre-probe identification status, as at 28 October 2026.
# Source: docs/findings.md finding 2 -- distinct affine equations per isoform
# (CYP2D6 four, the other three two each). A verdict, not a computation, so it is
# not derivable from this figure's own source file. PRE-PROBE: superseded for the
# three unidentified isoforms by outputs/56b_cyp1a2_corrected/solved_populations.csv
# on 29 October. See the module docstring.
IDENTIFIED_PRE_PROBE = {
    "CYP1A2": False,
    "CYP2C9": False,
    "CYP2D6": True,
    "CYP3A4": False,
}


def variant(target: Target) -> str:
    return "single-panel (all four isoforms)"


def submitted_means():
    """The thirteen submitted prediction means per isoform, as submitted."""
    # The index header is blank in this file, hence index_col=0.
    return load("56_cyp1a2_identification", "submitted_means.csv",
                index_col=0)[list(ISOFORMS)]


def draw(target: Target):
    means = submitted_means()
    spread = (means.max() - means.min()).sort_values(ascending=False)
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        ax = fig.add_subplot(111)

        isoforms = list(spread.index)            # widest range first, so CYP2D6 leads
        lo_all, hi_all = means.to_numpy().min(), means.to_numpy().max()
        span = hi_all - lo_all

        # Wilkinson dot plot: coincident submissions stack upward so the COUNT is visible.
        # Eight of CYP2D6's thirteen share a mean to six decimals, and rendering them as one dot
        # -- which alpha alone does -- undercuts the whole argument, since the two clusters and
        # their sizes are the finding.
        bin_width = span * 0.012
        dot = target.base_size * 0.34
        # One stack level, in row units. The tallest stack is eight, so this keeps a full stack
        # inside 8 x 0.075 = 0.6 of the 1.0 row pitch and clear of the row above.
        step = 0.075

        for i, iso in enumerate(isoforms):
            values = means[iso].to_numpy()
            colour = ISOFORM_COLOURS[iso]
            centre, level = wilkinson_bins(values, bin_width)

            # NO connecting rule. It would draw a line straight across the 1.347 gap between
            # CYP2D6's two clusters, which is the single thing this figure exists to show.
            # Each stack is CENTRED on its row baseline rather than growing upward from it, so
            # a tall stack reads as belonging to its own row and not to the one above.
            offsets = [(lv - level.max() / 2) * step for lv in level]
            ax.plot(centre, [i + o for o in offsets],
                    linestyle="none", marker="o", markersize=dot, color=colour,
                    alpha=0.9, markeredgecolor="white",
                    markeredgewidth=target.base_size * 0.035, zorder=3)

            tallest = int(level.max()) + 1
            if tallest > 1:
                k = int(level.argmax())
                ax.text(centre[k] + span * 0.013, i - (tallest / 2 + 0.2) * step,
                        f"x{tallest}", ha="left", va="center",
                        fontsize="small", color=TEXT)

            ax.text(values.max() + span * 0.030, i, f"range {spread[iso]:.3f}",
                    va="center", ha="left", fontsize="small", color=TEXT)

        ax.set_xlim(lo_all - span * 0.05, hi_all + span * 0.30)   # right: range annotations
        limit_ticks(ax, lo_all, hi_all)
        ax.set_xlabel("Submitted prediction mean (pIC50)")
        ax.xaxis.grid(True, zorder=0)
        ax.set_axisbelow(True)
        ax.set_ylim(len(isoforms) - 0.55, -0.55)

        # Status reads as one unit with the isoform name, on a second line beneath
        # it. Drawn as a secondary axis rather than free text so the layout engine
        # measures it and reserves the margin; free text outside the axes is
        # invisible to constrained layout and gets clipped.
        ax.set_yticks(range(len(isoforms)))
        ax.set_yticklabels(isoforms)
        status_axis = ax.secondary_yaxis("left")
        status_axis.set_yticks([i + 0.27 for i in range(len(isoforms))])
        status_axis.set_yticklabels(
            ["identified" if IDENTIFIED_PRE_PROBE[i] else "not identified"
             for i in isoforms],
            fontsize="small", color=TEXT,
        )
        status_axis.tick_params(length=0)
        status_axis.spines["left"].set_visible(False)

        figure_title(fig, target, "Identifiability is a property of the submission history")
        for i, iso in enumerate(isoforms):
            label = ax.get_yticklabels()[i]
            if IDENTIFIED_PRE_PROBE[iso]:
                label.set_fontweight("bold")

        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
