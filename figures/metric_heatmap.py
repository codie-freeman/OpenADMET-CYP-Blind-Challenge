"""Every board metric, every submission, every isoform. The annotated-heatmap genre figure.

Method and provenance live here; the caption is ``CAPTION`` below, which ``build.py`` writes to
``build/captions.md``.

Three panels, submissions as rows in true chronological order, isoforms as columns, the value
printed in every cell. Four things this gets right that an earlier version did not:

1. **Normalised per COLUMN, not per panel.** The claim is that ST-RAE and R-squared travel down the
   rows while Spearman does not. Normalising per panel let CYP2D6's much larger between-isoform
   range dominate the ramp, so three of four columns saturated and the very contrast the figure
   exists to show was invisible.
2. **ST-RAE's ramp is flipped**, so the dark end is the better end in all three panels. Otherwise
   dark means good in one panel and bad in the next.
3. **A colourbar**, per the project's sequential-colour rule. The earlier version had none.
4. **Conditional text colour**, so a value on a dark cell stays readable.

The ramp is the project's single sequential hue, which is deliberately not an isoform colour.
"""

from __future__ import annotations

import numpy as np

from figures.board import wide
from figures.style.base import Target, figure, figure_title, save
from figures.style.palette import ISOFORMS, TEXT, sequential

SLUG = "metric_heatmap"
SIZES = {"thesis": (6.3, 6.6), "slides": (9.0, 5.0), "web": (8.0, 7.2)}

#: (metric, lower_is_better)
PANELS = [("ST-RAE", True), ("R²", False), ("Spearman", False)]

CAPTION = """
Every scored submission against every isoform, for the three board metrics that matter, in true
submission order. **Each column is normalised within itself**, and each ramp is oriented so the
dark end is the better end — so a column darkening downward means that isoform improved over time,
whatever its absolute level.

**The contrast between the panels is the result.** ST-RAE and R-squared darken substantially down
their columns; Spearman barely changes at all. Almost every submission after 27b is an affine
correction of an earlier one, and placement and spread are exactly what an affine correction moves,
while Spearman is invariant under every such transform. Runs of identical Spearman values down a
column — CYP1A2's 0.7375 across NB16, NB19, NB27-widened and 27b — are not a rendering artefact:
those submissions share one ranking.

**Read the colour against the span printed under each column header, not on its own.** Per-column
normalisation stretches each column's own range across the full ramp, which is what makes the
within-column trend legible — but it also means a column spanning 0.046 looks exactly as vividly
contrasted as one spanning 0.349. CYP1A2's Spearman column runs 0.710 to 0.756, a range of **0.046**;
its ST-RAE column runs 0.563 to 0.913, a range of **0.349**, seven times larger. Without the printed
spans the colour would say the opposite of the claim.

Pooling across a panel instead was tried and is worse: CYP2D6's range swamps the other three, three
of four columns saturate to one shade, and the within-column trend disappears entirely.
"""


def variant(target: Target) -> str:
    return f"{len(PANELS)} panels, per-column normalisation"


def draw(target: Target):
    tables = {m: wide(m) for m, _ in PANELS}
    rows = list(tables["ST-RAE"].index)
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        axes = fig.subplots(1, len(PANELS), sharey=True)
        for ax, (metric, lower_better) in zip(axes, PANELS):
            t = tables[metric]
            v = t.to_numpy(dtype=float)

            # Per COLUMN: each isoform's own min and max set its own ramp.
            lo = np.nanmin(v, axis=0, keepdims=True)
            hi = np.nanmax(v, axis=0, keepdims=True)
            frac = (v - lo) / np.where(hi - lo == 0, 1.0, hi - lo)
            # Dark = better in EVERY panel: flip the fraction where lower is better.
            good = 1.0 - frac if lower_better else frac

            cmap = sequential()
            im = ax.imshow(good, cmap=cmap, aspect="auto", vmin=0.0, vmax=1.0)
            for i in range(v.shape[0]):
                for j in range(v.shape[1]):
                    ax.text(j, i, f"{v[i, j]:.3f}", ha="center", va="center", fontsize="small",
                            color="white" if good[i, j] > 0.55 else TEXT)

            # Per-column normalisation stretches a column's own range across the whole ramp, so
            # a column spanning 0.046 looks as vividly contrasted as one spanning 0.349 and
            # visually contradicts the claim. The actual span is printed under each header.
            spans = (hi - lo).ravel()
            ax.set_xticks(range(len(ISOFORMS)))
            # The number alone: "span 0.349" is wider than a one-column cell at four columns
            # across, and the word goes in the axis label once instead.
            ax.set_xticklabels([f"{i.replace('CYP', '')}\n{sp:.3f}"
                                for i, sp in zip(ISOFORMS, spans)], fontsize="small")
            ax.set_xlabel("isoform / span of this column", fontsize="small")
            ax.set_title(metric, color=TEXT)
            ax.set_xticks(np.arange(-0.5, len(ISOFORMS)), minor=True)
            ax.set_yticks(np.arange(-0.5, len(rows)), minor=True)
            ax.grid(which="minor", color="white", linewidth=0.8)
            ax.tick_params(which="minor", length=0)
            for sp in ax.spines.values():
                sp.set_visible(False)

        axes[0].set_yticks(range(len(rows)))
        axes[0].set_yticklabels(rows)

        cb = fig.colorbar(im, ax=axes, orientation="horizontal", fraction=0.04, pad=0.04,
                          ticks=[0.0, 1.0])
        cb.ax.set_xticklabels(["worst in column", "best in column"], fontsize="small")
        cb.outline.set_visible(False)
        cb.ax.tick_params(length=0)
        figure_title(fig, target, "Calibration moves the metrics; ranking does not")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
