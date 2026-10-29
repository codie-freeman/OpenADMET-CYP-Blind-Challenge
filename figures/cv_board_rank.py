"""THE LEAD FIGURE OF THE PAIR. Cross-validated rank against board rank, as a slopegraph.

Method and provenance live here; the caption is ``CAPTION`` below, which ``build.py`` writes to
``build/captions.md``. The join and the rules that make it small are in ``figures/cv_board_pairs.py``.

Structured so that no text sits in the data area: arm names in the left margin, ranks on both
axes, lines between. Crossing lines are the finding and are deliberately not sorted away.
"""

from __future__ import annotations

from figures.cv_board_pairs import spearman_pairs
from figures.style.base import Target, figure, figure_title, save
from figures.style.palette import ISOFORM_COLOURS, ISOFORMS, TEXT, arm_label

CAPTION = """
**The lead figure of the cross-validation-to-board pair.** Its companion shows that the *level*
transferred, and transferred favourably; this asks the question the validation protocol exists to
answer — whether it picked the right model — and the answer is that half the time it did not.

Twelve (arm, isoform) pairs, ranked by 25-fold cluster-disjoint out-of-fold Spearman on the left
and by board Spearman on the right, 1 best. One line per pair, coloured by isoform. **A flat line
is agreement; a crossing is disagreement, and the crossings are the result.** Within an isoform the
orderings agree on **two of four**, and **CYP3A4 is fully reversed** — the arm cross-validation
ranked best is the one the board ranked worst of the three. That count is computed from these same
twelve points and reproduces `docs/findings.md` finding 3 exactly.

A pooled rank correlation over all twelve is 0.81, but it overstates the agreement badly: Spearman
differs far more between isoforms than between arms within one, so pooling mostly measures that
both sources agree CYP3A4 ranks above CYP2D6, which is not the question.

**This is on Spearman, not ST-RAE, for two independent reasons.** ST-RAE is not comparable across
isoforms, so a pooled rank would largely sort by isoform. And Spearman is affine-invariant, so no
post-hoc correction can move it and all twelve pairs are valid model-versus-model reads with no
correction filter — which is why this figure has twelve points where the parity figure has ten.

**Coverage limitation.** Three of the eleven cross-validated arms were ever submitted. `deadzone`,
which holds the compact-letter-display `a` on all four isoforms, is not among them, so this is
measured over the arms that happened to be sent rather than the ones cross-validation ranked best.
"""

SLUG = "cv_board_rank"
SIZES = {"thesis": (6.3, 4.6), "slides": (9.0, 5.0), "web": (8.0, 5.0)}


def variant(target: Target) -> str:
    return f"slopegraph ({len(spearman_pairs())} pairs)"


def ranked():
    df = spearman_pairs().copy()
    df["cv_rank"] = df.cv.rank(ascending=False)
    df["board_rank"] = df.board.rank(ascending=False)
    return df


def agreement(df):
    """Within-isoform ordering agreement — the question the protocol exists to answer."""
    agree, reversed_ = 0, []
    for iso, g in df.groupby("isoform"):
        cv = list(g.sort_values("cv", ascending=False).arm)
        bd = list(g.sort_values("board", ascending=False).arm)
        agree += cv == bd
        if cv == bd[::-1]:
            reversed_.append(iso)
    return agree, reversed_


def draw(target: Target):
    df = ranked()
    n = len(df)
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        ax = fig.add_subplot(111)
        for _, r in df.iterrows():
            colour = ISOFORM_COLOURS[r.isoform]
            ax.plot([0, 1], [r.cv_rank, r.board_rank], color=colour,
                    linewidth=target.base_size * 0.14, alpha=0.85, zorder=3)
            ax.plot([0, 1], [r.cv_rank, r.board_rank], linestyle="none", marker="o",
                    color=colour, markersize=target.base_size * 0.38, zorder=4)

        # Arm names in the LEFT MARGIN, as y tick labels: outside the data area entirely, so
        # nothing can collide with a line. Rank is read off both axes.
        left = df.sort_values("cv_rank")
        ax.set_yticks(left.cv_rank)
        ax.set_yticklabels([arm_label(a) for a in left.arm])
        right = ax.secondary_yaxis("right")
        right.set_yticks(df.sort_values("board_rank").board_rank)
        right.set_yticklabels([f"{int(v)}" for v in sorted(df.board_rank)])
        right.tick_params(length=0)
        right.set_ylabel("Board rank (1 = best)")

        ax.set_ylim(n + 0.6, 0.4)
        ax.set_xlim(-0.06, 1.06)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["cross-validated rank", "board rank"])
        ax.set_ylabel("Cross-validated rank (1 = best out of fold)")
        for side in ("top", "right", "bottom"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="x", length=0)

        # One legend, isoform only: it is the only thing colour encodes here.
        for iso in ISOFORMS:
            ax.plot([], [], color=ISOFORM_COLOURS[iso], linewidth=target.base_size * 0.14,
                    label=iso)
        ax.legend(loc="center left", bbox_to_anchor=(1.12, 0.5), frameon=False,
                  handlelength=1.4, borderpad=0.2, labelcolor=TEXT)
        figure_title(fig, target, "Validation picked the right model half the time")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
