"""Models agree more about their predictions than about their errors -- the quantity that matters.

Method and provenance live here; the caption is ``CAPTION`` below, which ``build.py``
writes to ``build/captions.md`` alongside the renders.
"""

from __future__ import annotations

import pandas as pd

from figures.data import load
from figures.style.base import Target, figure, figure_title, save
from figures.style.palette import ISOFORM_COLOURS, NEUTRAL, REFERENCE, TEXT

CAPTION = """
Every pair among the eleven-configuration pool, on the two isoforms for which both correlations were
saved: raw out-of-fold prediction correlation on the x axis against residual (prediction minus
truth) error correlation on the y, with the identity as reference.

**Each point is one PAIR of models. Most sit above the line, and on CYP1A2 every one of the 55
does -- but on CYP3A4 only 37 of 55, so the effect is strong rather than universal and is reported
as such.** CYP1A2 also splits into two visible bands, at error correlation around 0.80-0.88 and
around 0.88-0.96; the upper band is the pairs that share a molecular representation and the lower
band the pairs that cross between graph and tabular features, which is the same graph-versus-tabular
diversity axis notebook 11b found in this pool. Two models can agree
closely on their raw outputs simply because both track the same true values; what determines whether
combining them helps is whether their *errors* decorrelate (Dietterich, 2000). On CYP1A2 the gap is
stark: prediction correlation spans 0.32 to 0.83 while error correlation never falls below 0.79 --
shared ignorance rather than useful diversity, which is the pattern behind this project's
repeated finding that candidate ensemble members sit at or above the pool-typical error-correlation
band rather than below it. On CYP3A4 the two
quantities overlap substantially (prediction 0.62 to 0.93, error 0.67 to 0.94), and a third of pairs
fall on or below the line.

**ONLY TWO OF FOUR ISOFORMS ARE SHOWN, AND THAT IS A DATA LIMITATION, NOT A CHOICE.** Error
correlation was saved for all four (notebook 24), but **raw-prediction correlation exists only for
CYP1A2 and CYP3A4** (notebook 11b). Producing it for CYP2C9 or CYP2D6 would mean recomputing from
``outputs/11_caruana_prep/``, which this layer does not do. Where both exist they are the same pool
on the same partition, so they are comparable.

Sources: ``outputs/24_nonnegative_stacking/error_correlation_long.csv`` and
``outputs/11b_caruana_selection/correlation_matrix_{CYP1A2,CYP3A4}_raw_prediction.csv``.
"""

SLUG = "error_vs_prediction_correlation"
SIZES = {"thesis": (6.3, 3.7), "slides": (9.0, 4.8), "web": (8.0, 4.3)}

#: Raw-prediction correlation was only ever saved for these two.
ISOFORMS_WITH_BOTH = ("CYP1A2", "CYP3A4")


def variant(target: Target) -> str:
    return f"two panels ({', '.join(ISOFORMS_WITH_BOTH)} -- the only isoforms with both)"


def pairs() -> pd.DataFrame:
    err = load("24_nonnegative_stacking", "error_correlation_long.csv")
    out = []
    for iso in ISOFORMS_WITH_BOTH:
        raw = load("11b_caruana_selection", f"correlation_matrix_{iso}_raw_prediction.csv",
                   index_col=0)
        e = err[err.isoform == iso]
        wide = e.pivot(index="config_1", columns="config_2", values="error_correlation")
        configs = sorted(set(raw.index) & set(wide.index))
        assert len(configs) == 11, f"{iso}: {len(configs)} shared configs"
        for i, a in enumerate(configs):
            for b in configs[i + 1:]:
                out.append({"isoform": iso, "a": a, "b": b,
                            "prediction": float(raw.loc[a, b]),
                            "error": float(wide.loc[a, b])})
    df = pd.DataFrame(out)
    assert len(df) == len(ISOFORMS_WITH_BOTH) * 55, f"{len(df)} pairs"
    return df


def draw(target: Target):
    df = pairs()
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        # NOT sharey: each panel has its own range and an equal aspect, and sharing the y axis
        # with a fixed aspect squashes one panel and overflows the other.
        axes = fig.subplots(1, len(ISOFORMS_WITH_BOTH))
        fig.get_layout_engine().set(h_pad=0.08)       # keep the panel titles off the top edge
        for ax, iso in zip(axes, ISOFORMS_WITH_BOTH):
            sub = df[df.isoform == iso]
            colour = ISOFORM_COLOURS[iso]
            lo = min(sub.prediction.min(), sub.error.min()) - 0.03
            hi = max(sub.prediction.max(), sub.error.max()) + 0.03
            ax.plot([lo, hi], [lo, hi], color=REFERENCE, linestyle="--",
                    linewidth=target.base_size * 0.08, zorder=1)
            # Rule 1: one series per panel, so the marks take the isoform colour.
            ax.scatter(sub.prediction, sub.error, s=target.base_size ** 2 * 0.18,
                       color=colour, alpha=0.6, edgecolor="none", zorder=3)
            above = int((sub.error > sub.prediction).sum())
            # Below the diagonal is empty by construction wherever most pairs sit above it.
            ax.text(0.96, 0.06, f"{above} of {len(sub)} above the line",
                    transform=ax.transAxes, va="bottom", ha="right",
                    fontsize="small", color=TEXT,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.6))
            ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
            ax.set_title(iso, color=colour)
            ax.set_xlabel("Raw prediction correlation\n(one point = one model pair)",
                          fontsize="small")
            ax.grid(True, zorder=0); ax.set_axisbelow(True)
        for ax in axes:
            ax.set_ylabel("Error (residual) correlation")
        figure_title(fig, target, "Models share errors more than they share predictions")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
