"""Cross-validated ST-RAE against blind-board ST-RAE. The LEVEL transferred; see the rank figure.

Method and provenance live here; the caption is ``CAPTION`` below, written to
``build/captions.md`` by ``build.py``. The join and the rules that make it small are documented in
``figures/cv_board_pairs.py``.

Neither this chart nor its rank companion appears in any of the comparator papers, which make the
cross-validation-to-board link in prose only.
"""

from __future__ import annotations

from matplotlib.lines import Line2D

from figures.cv_board_pairs import strae_pairs
from figures.style.base import Target, figure, figure_title, save
from figures.style.palette import ISOFORM_COLOURS, ISOFORMS, REFERENCE, TEXT, arm_label

#: Arm is encoded by MARKER SHAPE, isoform by colour, so no point needs a text label.
ARM_MARKER = {"plain": "o", "aid": "s", "mae_only": "^"}

SLUG = "cv_board_parity"
SIZES = {"thesis": (6.3, 4.4), "slides": (9.0, 5.0), "web": (8.0, 4.8)}

CAPTION = """
Each point is one (modelling arm, isoform) pair for which an **uncorrected** board column exists:
25-fold cluster-disjoint cross-validated ST-RAE against blind-board ST-RAE, with the identity as
reference. **Marker shape is the arm, colour is the isoform**, both in one legend, so the data area
carries no text. Horizontal bars are the fold-to-fold standard deviation across the
25 folds.

**Seven of the ten points sit BELOW the diagonal: most arms scored better on the board than they
did out of fold.** The level transferred, and it transferred in this project's favour. The likely
reason is in the construction of the test set — the blind 750 are hit-expanded catalogue analogues
of known actives, so they sit closer to the training compounds than a deliberately
cluster-disjoint held-out fold does, and a cluster-disjoint fold is therefore the harder test.
The three exceptions are both CYP2D6 points and `aid` on CYP3A4.

**So the cross-validation-to-board failure on this project is not one of level. It is one of
ordering**, which this chart cannot show and its rank companion can.

**There are ten points, and that is the honest number.** Only three of the eleven cross-validated
arms were ever submitted (`plain` as 10c, `aid` as NB30, `mae_only` as NB47-mae); the other eight —
including `deadzone`, which holds the compact-letter-display `a` on all four isoforms — have no
board reading at all. Of the twelve cells those three arms span, two carried a post-hoc correction
(NB47-mae's CYP1A2 and CYP2D6) and are excluded, because a placement or spread correction moves
ST-RAE directly and the comparison would stop being model against model. `04b` is also excluded: it
maps to the same arm as `plain` but is a different training draw in an older environment, and one
cross-validated mean cannot rank two draws of itself. `NB16` adds nothing, because its uncorrected
columns reproduce 10c's board values exactly.

**There are no vertical error bars, deliberately.** A board reading is exact for a fixed file: the
board returns identical scores for identical bytes, confirmed six times on this record by
byte-identical re-sends. There is no measurement interval to draw. The uncertainty that does exist
is in generalisation, and it lives in the half-set-against-full-set figure, which is the only such
estimate the record contains.
"""


def variant(target: Target) -> str:
    return f"single-panel ({len(strae_pairs())} uncorrected pairs)"


def draw(target: Target):
    df = strae_pairs()
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        ax = fig.add_subplot(111)

        lo = min(df.board.min(), (df.cv - df.cv_sd).min())
        hi = max(df.board.max(), (df.cv + df.cv_sd).max())
        pad = (hi - lo) * 0.08
        lims = (lo - pad, hi + pad)
        ax.plot(lims, lims, color=REFERENCE, linestyle="--",
                linewidth=target.base_size * 0.07, zorder=1)

        for _, r in df.iterrows():
            colour = ISOFORM_COLOURS[r.isoform]
            # Rule 4: the bar takes its own mark's colour at reduced weight -- thinner than the
            # marker edge, alpha 0.6, no caps -- so it reads as uncertainty rather than as data.
            ax.errorbar(r.cv, r.board, xerr=r.cv_sd, fmt="none",
                        ecolor=colour, elinewidth=target.base_size * 0.045,
                        alpha=0.6, capsize=0, zorder=3)
            ax.plot(r.cv, r.board, marker=ARM_MARKER[r.arm], color=colour,
                    markersize=target.base_size * 0.46, markeredgecolor="none", zorder=4)

        below = int((df.board < df.cv).sum())
        ax.text(0.03, 0.97, f"{below} of {len(df)} below the line: better on the board",
                transform=ax.transAxes, va="top", ha="left",
                fontsize="small", color=TEXT)

        # Extra room on the right so the rightmost arm label is not clipped.
        ax.set_xlim(lims[0], lims[1] + (hi - lo) * 0.22)
        ax.set_ylim(*lims)
        ax.set_xlabel("25-fold cluster-disjoint cross-validated ST-RAE")
        ax.set_ylabel("Blind-board ST-RAE")
        ax.grid(True, zorder=0); ax.set_axisbelow(True)

        # ONE legend carrying both encodings, so no point needs a label: shape is the arm,
        # colour is the isoform. This is what lets the data area hold no text at all.
        handles = [Line2D([], [], linestyle="none", marker=ARM_MARKER[a], color=TEXT,
                          markersize=target.base_size * 0.42, label=arm_label(a))
                   for a in ARM_MARKER]
        handles += [Line2D([], [], linestyle="none", marker="o", color=ISOFORM_COLOURS[i],
                           markersize=target.base_size * 0.42,
                           label=f"{i} ({int((df.isoform == i).sum())})")
                    for i in ISOFORMS]
        ax.legend(handles=handles, loc="lower right", frameon=False, ncol=2,
                  handletextpad=0.3, borderpad=0.2, columnspacing=1.0, labelcolor=TEXT)

        figure_title(fig, target, "The level transferred, and in our favour")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
