"""Are the compounds every model gets wrong the ones furthest from training? No.

Method and provenance live here; the caption is ``CAPTION`` below, which ``build.py``
writes to ``build/captions.md`` alongside the renders.
"""

from __future__ import annotations

from figures.data import load
from figures.style.base import Target, figure, figure_title, save
from figures.style.palette import ISOFORM_COLOURS, ISOFORMS, NEUTRAL, TEXT

CAPTION = """
For each isoform, the distribution of nearest-neighbour Tanimoto similarity to the training set
(top row) directly above the distribution of soft-threshold error (bottom row), on a shared isoform
facet, split into the consensus-hard and consensus-easy deciles.

**The bottom row separates by construction -- that is what "hard" means -- and the top row does
not separate at all.** Notebook 25 measured the hard-versus-easy difference in nearest-neighbour
similarity at a rank-biserial |r| of 0.02 to 0.09, significant on no isoform after
Benjamini-Hochberg correction. So the compounds every one of eleven configurations gets wrong are
structurally **ordinary** relative to training, not distant from it, uniformly across all four
isoforms. A retrieved-neighbour corpus aimed at "unseen chemistry" is aimed at the wrong population.

**The easy decile's soft-threshold error is almost all exactly zero**, which is why the bottom row
is a spike at the origin plus a tail: a prediction that falls inside its compound's credible interval
scores zero under this metric by definition. The two deciles have identical n, so the histograms are
counts rather than densities and are directly comparable.

Hard and easy are the worst and best deciles by mean soft-threshold error across all eleven
configurations, which is notebook 25's primary definition (its plain-absolute-error sets differ
materially -- Jaccard overlap 0.06 to 0.48).

Source: ``outputs/25_hard_compound_analysis/hard_compound_table.csv``.
"""

SLUG = "hard_compound_similarity"
SIZES = {"thesis": (6.3, 4.0), "slides": (9.0, 5.0), "web": (8.0, 4.6)}

ROWS = [("loo_nn_similarity", "NN similarity\nto training"),
        ("mean_soft_threshold_error", "Soft-threshold\nerror")]


def variant(target: Target) -> str:
    return "2 x 4 grid (similarity over error, per isoform)"


def table():
    t = load("25_hard_compound_analysis", "hard_compound_table.csv")
    assert set(t.isoform) == set(ISOFORMS)
    assert t.loo_nn_similarity.notna().all()
    return t


def draw(target: Target):
    t = table()
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        axes = fig.subplots(2, 4, sharex="row")
        for col, iso in enumerate(ISOFORMS):
            sub = t[t.isoform == iso]
            hard = sub[sub.hard_soft.astype(bool)]
            easy = sub[sub.easy_soft.astype(bool)]
            colour = ISOFORM_COLOURS[iso]

            for row, (field, _) in enumerate(ROWS):
                ax = axes[row][col]
                # Counts, not density: the two deciles have IDENTICAL n by construction
                # (they are the worst and best tenth of the same population), so counts are
                # directly comparable and density only rescales the easy decile's zero spike.
                bins = 24
                ax.hist(easy[field], bins=bins, color=NEUTRAL, alpha=0.30, zorder=2)
                ax.hist(hard[field], bins=bins, color=colour, alpha=0.65, zorder=3)
                for data, c in ((easy[field], NEUTRAL), (hard[field], colour)):
                    ax.axvline(data.median(), color=c, linestyle="--",
                               linewidth=target.base_size * 0.07, zorder=4)
                ax.set_yticks([])
                if row == 0:
                    ax.set_title(iso, color=colour, pad=target.base_size * 1.9)
                    # Above the axes, not inside: inside it lands on the histogram on CYP2D6
                    # and CYP3A4, whose mass sits to the right.
                    # The right title SLOT: a title may sit outside the axes box, free text
                    # above it may not.
                    # Under the x label rather than in the right title slot, which on a narrow
                    # four-column grid ran straight into the isoform name beside it.
                    ax.set_xlabel(f"Tanimoto\nmedian hard {hard[field].median():.2f} / "
                                  f"easy {easy[field].median():.2f}", fontsize="small")
                if col == 0:
                    # A row TITLE, not a y-axis label: the y axis is a count and the x axis is
                    # Tanimoto, so labelling the y axis "NN similarity" would misread. Drawn as a
                    # FIGURE artist because it deliberately sits outside the axes.
                    fig.text(-0.34, 0.5, ROWS[row][1], transform=ax.transAxes,
                             rotation=90, ha="center", va="center",
                             fontweight="bold", color=TEXT)

        for col in range(4):
            axes[1][col].set_xlabel("pIC50", fontsize="small")
        fig.supxlabel("coloured = consensus-hard decile, grey = consensus-easy decile; "
                      "dashed = median", fontsize="small", color=TEXT)
        fig.get_layout_engine().set(w_pad=0.10)
        figure_title(fig, target, "Hard compounds are not the distant ones")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
