"""Credible-band width against measured potency. The precondition for the dead-zone target.

Both quantities are read per compound from notebook 25's saved table; the per-isoform Spearman is
read from notebook 21's own ``interval_summary.csv`` rather than recomputed here.
"""

from __future__ import annotations

from figures.data import load
from figures.style.base import Target, figure, figure_title, save
from figures.style.palette import ISOFORM_COLOURS, ISOFORMS, NEUTRAL, TEXT

SLUG = "interval_width_vs_potency"
SIZES = {"thesis": (6.3, 3.6), "slides": (9.0, 4.8), "web": (8.0, 4.2)}

CAPTION = """
Each training compound's credible-interval width against its measured pIC50, by isoform, with the
rank correlation annotated. **Width falls steeply with potency on every isoform (Spearman −0.56 to
−0.93): weak compounds carry wide bands and potent ones carry narrow bands.**

This one fact explains two things that otherwise arrive unmotivated.

**It makes the scored metric asymmetric.** ST-RAE charges nothing while a prediction sits inside
its compound's credible band and charges from the nearest edge once it leaves. A wide band at the
inactive end is therefore cheap to over-predict into and a narrow band at the potent end is
expensive to miss — so at raw placement, total under-prediction error exceeds total
over-prediction error by 1.65 to 2.79 times across the isoforms, and the ST-RAE-optimal uniform
offset sits 0.06 to 0.12 pIC50 units *above* the R-squared-optimal one.

**And it is what the dead-zone training target exploits.** Fitting to the prediction clipped into
the band, under absolute-error loss, stops spending model capacity on precision the assay never
had — which is why that target is uniquely best on all four isoforms under cross-validation.

Sources: per-compound widths from ``outputs/25_hard_compound_analysis/hard_compound_table.csv``;
correlations from ``outputs/21_strae_offset_curve/interval_summary.csv``, computed there, not here.
"""


def variant(target: Target) -> str:
    return "1 x 4 facets (one per isoform)"


def table():
    t = load("25_hard_compound_analysis", "hard_compound_table.csv")
    stats = load("21_strae_offset_curve", "interval_summary.csv").set_index("isoform")
    assert (stats.n_conf_high_lt_conf_low == 0).all(), "a credible interval is inverted"
    return t, stats


def draw(target: Target):
    t, stats = table()
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        axes = fig.subplots(1, 4, sharey=True)
        for ax, iso in zip(axes, ISOFORMS):
            sub = t[t.isoform == iso].drop_duplicates("inchikey")
            colour = ISOFORM_COLOURS[iso]
            # Rule 1: one series per panel, so the marks take the isoform colour.
            ax.scatter(sub.y_true, sub.interval_width, s=target.base_size ** 2 * 0.06,
                       color=colour, alpha=0.3, edgecolor="none", zorder=3)
            rho = float(stats.loc[iso, "spearman_width_vs_y_true"])
            ax.set_title(iso, color=colour)
            ax.text(0.95, 0.95, f"rho {rho:+.2f}", transform=ax.transAxes,
                    ha="right", va="top", fontsize="small", color=TEXT,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.6))
            ax.set_xlabel("pIC50")
            ax.grid(True, zorder=0); ax.set_axisbelow(True)
        axes[0].set_ylabel("Credible-interval width (pIC50)")
        figure_title(fig, target, "Weak compounds carry wide bands; potent ones do not")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
