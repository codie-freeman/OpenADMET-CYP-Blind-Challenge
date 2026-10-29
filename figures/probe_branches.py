"""The identification probe: predicted R-squared per mirror branch against the translation applied.

Method and provenance live here; the caption is ``CAPTION`` below, which ``build.py``
writes to ``build/captions.md`` alongside the renders.

THE IDENTITY, AND WHY THIS FIGURE EVALUATES IT
----------------------------------------------
Notebook 56 solved each branch's parameters and stored them, but stored only two points per branch
-- R-squared at a translation of 0 and at +0.25. R-squared is **quadratic** in the translation, so
joining two points with a straight line would misrepresent the shape the figure exists to show.
This figure therefore evaluates the identity the notebook already solved, over a grid:

    R2(d) = 2*rho*k - k^2 - b(d)^2,    b(d) = (pred_mean + d - ybar) / sd_y,   k = pred_sd / sd_y

with ``rho``, ``ybar``, ``sd_y``, ``k`` and ``b(0)`` all read from
``outputs/56_cyp1a2_identification/cyp1a2_roots.csv``. Nothing is re-solved. As verification, the
two values notebook 56 stored for each branch in ``branch_predictions.csv`` are overplotted as open
marks and must lie on the curves. See ``figures/README.md`` for where this sits relative to the
rule that figures read rather than compute.
"""

from __future__ import annotations

import numpy as np

from figures.data import load
from figures.style.base import Target, figure, figure_title, place_labels, save
from figures.style.palette import ISOFORM_COLOURS, REFERENCE, TEXT

CAPTION = """
Before the probe, CYP1A2's blind population was consistent with **four** candidate solutions in two
mirror pairs, because its two distinct submitted equations shared a prediction mean and a
mean-preserving variant does not break the sign ambiguity. Each curve is one candidate's predicted
board R-squared as a function of a uniform translation of the submitted column; the mirror pairs
separate under translation, which is the whole reason a shifted column identifies the system. The
+0.25 pIC50 probe (notebook 56) is marked, with each branch's predicted reading at it, and the
**actual board reading of 0.1154** overlaid.

**It lands on the `ybar = 4.3007, sd_y = 1.4912` branch** -- predicted 0.1107, observed 0.1154, a
difference of 0.0047 against the 0.02 bar notebook 56 pre-registered. So the LOW-centred pair is
correct and the hit-expansion sign prior, which predicted HIGH, is falsified; notebook 56b then
solved all three remaining isoforms exactly from that reading.
"""

SLUG = "probe_branches"
SIZES = {"thesis": (6.3, 3.6), "slides": (9.0, 5.0), "web": (8.0, 4.2)}

ISOFORM = "CYP1A2"
PROBE_SHIFT = 0.25
OBSERVED_R2 = 0.1154        # docs/leaderboard_submissions.md, NB56-identify, CYP1A2 R2


def variant(target: Target) -> str:
    return f"single-panel ({ISOFORM}, {len(branches())} mirror branches)"


def branches():
    b = load("56_cyp1a2_identification", "cyp1a2_roots.csv")
    return b[b.isoform == ISOFORM].reset_index(drop=True)


def stored_points():
    s = load("56_cyp1a2_identification", "branch_predictions.csv")
    return s[s.isoform == ISOFORM].reset_index(drop=True)


def draw(target: Target):
    br, stored = branches(), stored_points()
    colour = ISOFORM_COLOURS[ISOFORM]
    width, height = SIZES[target.name]

    d = np.linspace(-0.6, 0.6, 241)

    with figure(target, width, height) as fig:
        ax = fig.add_subplot(111)
        ax.axvline(0, color=REFERENCE, linewidth=target.base_size * 0.06, zorder=1)
        ax.axvline(PROBE_SHIFT, color=REFERENCE, linestyle="--",
                   linewidth=target.base_size * 0.08, zorder=2)

        dashes = [(), (), (1, 1.6), (4, 1.6)]      # two HIGH branches get distinct patterns
        hi_seen = 0
        curves: list = []
        label_pts: list = []
        for _, r in br.iterrows():
            # b(d) = b(0) + d / sd_y, since a uniform translation moves only the prediction mean.
            b_d = r.b + d / r.sd_y
            r2 = 2 * r.rho_bvn * r.k - r.k ** 2 - b_d ** 2
            high = r.b < 0                      # the mirror partner above the training mean
            style = {}
            if high:
                style["dashes"] = dashes[2 + hi_seen]
                hi_seen += 1
            ax.plot(d, r2, color=colour, alpha=0.95 if not high else 0.55,
                    linewidth=target.base_size * 0.12, zorder=3, **style)
            at_probe = 2 * r.rho_bvn * r.k - r.k ** 2 - (r.b + PROBE_SHIFT / r.sd_y) ** 2
            curves.append((d, r2))
            label_pts.append((PROBE_SHIFT, at_probe,
                              f"ybar {r.ybar:.4f}, sd {r.sd_y:.3f}: {at_probe:.4f}"))

        # Verification: notebook 56's own two stored values per branch must sit on the curves.
        ax.scatter(np.zeros(len(stored)), stored.R2_now, facecolors="none",
                   edgecolors=TEXT, s=target.base_size ** 2 * 0.5, zorder=4)
        ax.scatter(np.full(len(stored), PROBE_SHIFT), stored.R2_at_probe, facecolors="none",
                   edgecolors=TEXT, s=target.base_size ** 2 * 0.5, zorder=4)

        ax.plot([PROBE_SHIFT], [OBSERVED_R2], marker="*", color=TEXT,
                markersize=target.base_size * 1.3, zorder=6)

        # All five annotations -- four branch readings and the observed value -- are placed
        # together with every curve in `avoid`, so none lands on a line or on another label.
        label_pts.append((PROBE_SHIFT, OBSERVED_R2, f"observed {OBSERVED_R2:.4f}"))
        place_labels(ax, [p[0] for p in label_pts], [p[1] for p in label_pts],
                     [p[2] for p in label_pts], colours=TEXT, avoid=curves,
                     radius_pt=target.base_size * 0.95)

        ax.set_xlabel(f"Uniform translation applied to the submitted {ISOFORM} column (pIC50)")
        ax.set_ylabel("Predicted board R-squared")
        ax.set_xlim(d.min() - 0.12, d.max() + 0.22)
        ax.grid(True, zorder=0); ax.set_axisbelow(True)
        # A FIGURE artist: at 16 pt this legend line is wider than the axes, and it is a note
        # about the figure rather than a label on part of it.
        fig.supxlabel("solid = LOW-centred, dashed = HIGH-centred; "
                      "open marks = notebook 56's stored values",
                      fontsize="small", color=TEXT)
        figure_title(fig, target, "One shifted column identifies the blind population")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
