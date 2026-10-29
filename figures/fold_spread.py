"""The 25 fold values behind every mean and interval in the forest plot.

Method and provenance live here; the caption is ``CAPTION`` below, which ``build.py`` writes to
``build/captions.md``.

Per-fold scores are **read**, not recomputed, and they are split across three notebooks: repeats 0
and 1 are notebook 34's, repeat 2 is notebook 35's, repeats 3 and 4 are notebook 37's, and
``mae_only``'s repeats 0 to 2 are **all** notebook 35's because that arm was added later. Taking
notebook 37's own ``run_scores.csv`` alone would silently give 10 folds labelled as 25.
"""

from __future__ import annotations

import pandas as pd

from figures.data import load
from figures.pairwise_significance import arm_order
from figures.style.base import Target, figure, figure_title, save
from figures.style.palette import ISOFORM_COLOURS, ISOFORMS, NEUTRAL, REFERENCE, TEXT, arm_label

SLUG = "fold_spread"
SIZES = {"thesis": (6.3, 6.2), "slides": (9.0, 5.0), "web": (8.0, 6.8)}

#: Which notebook holds which repeats. See the module docstring.
REPEAT_SOURCE = {"34_butina_5x5": (0, 1), "35_butina_repeat3": (2,),
                 "37_butina_5x5_complete": (3, 4)}

CAPTION = """
The 25 cluster-disjoint fold scores behind every point in the compact-letter-display forest, by
isoform. Each strip is one arm; each mark is one fold; the heavy line is the arm's mean. Arms are
in the forest's own order, best at top, with the best arm accented.

The forest plot shows a mean and an interval. **This shows what they are made of, and the folds
disagree by far more than the arms do.** On every isoform the within-arm spread across folds is
wider than most of the gaps between arm means — which is why the protocol needs 25 *paired* folds
rather than a comparison of means, and why pairs that look far apart on the forest can still fail
to separate in the pairwise matrix.

**What decides whether a pair separates is the ratio of the two, not the spread alone**, and
checking that corrects an obvious-sounding reading of this figure. CYP2D6 has the fewest separated
pairs (36 of 55) but **not** the widest folds — its median within-arm fold SD is 0.048, the second
narrowest of the four. What is narrow on CYP2D6 is the gap between the arms: their means span only
0.195 ST-RAE, so the arm-range-to-fold-spread ratio is 4.1, the lowest of the four. CYP3A4 has the
tightest folds (0.037) and a wide arm range (0.332), giving 8.9 and 49 separated pairs. CYP2C9 has
the *widest* folds of all (0.057) and still separates the most (51), because its arms are spread
wider still (0.393).

So CYP2D6 is hard to resolve because its models are close together, not because its folds are
noisy.
"""


def variant(target: Target) -> str:
    return "2 x 2 grid (11 arms x 25 folds, per isoform)"


def per_fold() -> pd.DataFrame:
    parts = []
    for directory, repeats in REPEAT_SOURCE.items():
        t = load(directory, "run_scores.csv")
        parts.append(t[t.repeat.isin(repeats)])
    # mae_only's repeats 0-2 are all notebook 35's: it was added after notebook 34 ran.
    extra = load("35_butina_repeat3", "run_scores.csv")
    parts.append(extra[(extra.arm == "mae_only") & extra.repeat.isin((0, 1))])

    s = pd.concat(parts, ignore_index=True)
    s = s[s.Endpoint.str.startswith("CYP")].copy()
    s["isoform"] = s.Endpoint.str.split("_").str[0]
    cells = s.groupby(["isoform", "arm"]).apply(
        lambda g: g[["repeat", "fold"]].drop_duplicates().shape[0], include_groups=False)
    assert (cells == 25).all(), f"not 25 folds everywhere: {sorted(cells.unique())}"
    return s


def draw(target: Target):
    s = per_fold()
    order = arm_order()
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        axes = fig.subplots(2, 2, sharey=True)
        for ax, iso in zip(axes.ravel(), ISOFORMS):
            sub = s[s.isoform == iso]
            accent = ISOFORM_COLOURS[iso]
            for i, arm in enumerate(order):
                v = sub.loc[sub.arm == arm, "ST-RAE"].to_numpy()
                best = i == 0
                colour = accent if best else NEUTRAL
                ax.plot(v, [i] * len(v), linestyle="none", marker="o",
                        markersize=target.base_size * 0.26, color=colour,
                        alpha=0.45, markeredgecolor="none", zorder=3)
                ax.plot([v.mean()], [i], marker="|", color=colour,
                        markersize=target.base_size * 1.1,
                        markeredgewidth=target.base_size * 0.16, zorder=4)
            ax.set_title(iso, color=accent)
            ax.set_ylim(len(order) - 0.4, -0.6)
            ax.xaxis.grid(True, zorder=0)
            ax.set_axisbelow(True)

        for ax in axes[:, 0]:
            ax.set_yticks(range(len(order)))
            ax.set_yticklabels([arm_label(a) for a in order])
        for ax in axes[1, :]:
            ax.set_xlabel("ST-RAE")
        fig.supxlabel("one mark per fold, 25 per arm; the heavy tick is the arm's mean",
                      fontsize="small", color=TEXT)
        figure_title(fig, target, "The folds disagree more than the arms do")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
