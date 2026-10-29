"""Residual against nearest-neighbour similarity. The project's original hypothesis, finally drawn.

Method and provenance live here; the caption is ``CAPTION`` below. Both quantities are stored
columns of notebook 25's ``hard_compound_table.csv``; the per-isoform correlation is read from its
own ``hard_vs_easy_tests.csv`` where available and otherwise stated as unavailable rather than
recomputed.
"""

from __future__ import annotations

from figures.data import load
from figures.style.base import Target, figure, figure_title, save
from figures.style.palette import ISOFORM_COLOURS, ISOFORMS, REFERENCE, TEXT

SLUG = "residual_vs_similarity"
SIZES = {"thesis": (6.3, 3.4), "slides": (9.0, 4.8), "web": (8.0, 4.0)}

CAPTION = """
Mean absolute error across the eleven cross-validated configurations against leave-one-out
nearest-neighbour Tanimoto similarity to the rest of the training set, in decile bins with a 95 per
cent interval on each bin mean, one panel per isoform. **All four panels share one x-range**, which
they must for the comparison to be readable.

**This was the project's founding hypothesis.** The expectation was that the compounds the models
get wrong are the ones furthest from anything they were trained on — that error should fall as
similarity rises. **Binning changed the answer, and the figure reports what it shows rather than what was expected.**

On **CYP2C9 and CYP2D6 there is no trend**: first and last decile differ by 0.03 and 0.003 mean
absolute error, with wholly overlapping intervals. That is the flat null.

On **CYP1A2 error RISES with similarity** — 0.499 in the least-similar decile to 0.753 in the
most-similar, intervals not overlapping — which is the *opposite* of the founding hypothesis. On
**CYP3A4 it falls**, 0.648 to 0.439, also non-overlapping, which is the hypothesis direction. So two
isoforms show real monotone structure in opposite directions and two show none.

**This does not contradict notebook 25's stored result, because the two measure different things.**
Its rank-biserial |r| of 0.02 to 0.09, significant on no isoform, compares the *consensus-hard
decile against the consensus-easy decile* — the extremes of the error distribution. This compares
mean error across deciles of *similarity*. A population can show no hard-versus-easy separation and
still carry a trend across similarity, and on two isoforms it does.

What survives unchanged is the conclusion that mattered: nearest-neighbour distance does not
identify the compounds the pool gets wrong, since the effect is absent on half the isoforms and
points the wrong way on a quarter. A retrieved-neighbour corpus aimed at unseen chemistry would
have been aimed at CYP3A4 alone, and the wrong way round on CYP1A2.

The consequence was concrete: a retrieved-neighbour pretraining corpus, aimed at covering "unseen
chemistry", was planned and then not built. What *is* hard turned out to be potency extremes rather
than structural distance — the structure-grid figure shows that directly.
"""


def variant(target: Target) -> str:
    return "1 x 4 facets (one per isoform)"


def table():
    t = load("25_hard_compound_analysis", "hard_compound_table.csv")
    assert t.loo_nn_similarity.notna().all()
    return t


def effect_sizes():
    """Notebook 25's own stored rank-biserial effect size for this property, if present."""
    tests = load("25_hard_compound_analysis", "hard_vs_easy_tests.csv")
    prop = [c for c in tests.columns if "prop" in c.lower() or "variable" in c.lower()]
    if not prop:
        return {}
    col = prop[0]
    sub = tests[tests[col].astype(str).str.contains("nn_similarity", case=False, na=False)]
    eff = [c for c in tests.columns if "biserial" in c.lower() or c.lower() in ("r", "effect")]
    if sub.empty or not eff:
        return {}
    return dict(zip(sub.isoform, sub[eff[0]]))


def binned(t, n_bins: int = 10):
    """Mean absolute error per within-isoform NN-similarity decile, with a 95% interval.

    A scatter of 1,300 to 2,300 points per isoform renders as one blob, which is the wrong picture
    of a null: a flat line with overlapping intervals says "no relationship" at a glance and a blob
    says nothing at all. Binning is a reshape of stored values; the interval is the standard error
    of the mean within each bin.
    """
    import numpy as np
    import pandas as pd
    rows = []
    for iso, g in t.groupby("isoform"):
        q = pd.qcut(g.loo_nn_similarity, n_bins, duplicates="drop")
        for interval, gg in g.groupby(q, observed=True):
            m = gg.mean_abs_error
            rows.append({"isoform": iso, "centre": interval.mid, "n": len(gg),
                         "mean": m.mean(), "stderr": m.std(ddof=1) / np.sqrt(len(gg))})
    return pd.DataFrame(rows)


def draw(target: Target):
    t = table()
    b = binned(t)
    eff = effect_sizes()
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        axes = fig.subplots(1, 4, sharey=True, sharex=True)      # SHARED x: see the caption
        for ax, iso in zip(axes, ISOFORMS):
            g = b[b.isoform == iso].sort_values("centre")
            colour = ISOFORM_COLOURS[iso]
            overall = t.loc[t.isoform == iso, "mean_abs_error"].mean()
            ax.axhline(overall, color=REFERENCE, linestyle="--",
                       linewidth=target.base_size * 0.07, zorder=2)
            # g["stderr"], not g.stderr: a column named "sem" shadowed DataFrame.sem and
            # attribute access returned the method. Named to avoid the collision.
            ax.errorbar(g.centre, g["mean"], yerr=1.96 * g["stderr"], fmt="o-",
                        color=colour, ecolor=colour, alpha=0.9,
                        elinewidth=target.base_size * 0.05, capsize=0,
                        linewidth=target.base_size * 0.10,
                        markersize=target.base_size * 0.34, zorder=3)
            ax.set_title(iso, color=colour)
            if iso in eff:
                ax.set_title(f"r = {eff[iso]:+.2f}", loc="right", fontsize="small", color=TEXT)
            ax.set_xlabel("NN Tanimoto")
            ax.grid(True, zorder=0)
            ax.set_axisbelow(True)
        axes[0].set_ylabel("Mean absolute error\n(11 configurations)")
        fig.supxlabel("decile bins, 95% interval on the bin mean; "
                      "dashed = that isoform's overall mean error",
                      fontsize="small", color=TEXT)
        figure_title(fig, target, "Hard compounds are not the distant ones")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
