"""Board ST-RAE against board R-squared, per isoform, with the corrected series picked out.

Method and provenance live here; the caption is ``CAPTION`` below, which ``build.py``
writes to ``build/captions.md`` alongside the renders.

ON THE TREND LINE
-----------------
A least-squares line over the points the figure already displays is a presentation aid, not a
result: it asserts nothing the scatter does not already show, and no number from it is reported.
See ``figures/README.md`` for where that sits relative to the rule that figures read rather than
compute.
"""

from __future__ import annotations

import numpy as np

from figures.board import wide
from matplotlib.lines import Line2D

from figures.style.base import Target, figure, figure_title, place_labels, save
from figures.style.palette import ISOFORM_COLOURS, NEUTRAL, REFERENCE, TEXT

CAPTION = """
Every submission as a point, board ST-RAE against board R-squared, for **the two isoforms that
have a correction chain to show**. CYP2D6 and CYP3A4 are omitted rather than drawn empty: on
CYP2D6 all three chain submissions carry an identical column and on CYP3A4 two of the three
coincide, so the chain there is a single mark. The two are
tightly related because ST-RAE and R-squared both respond to placement and spread, which is what
an affine correction moves; the grey line is a least-squares trend over the plotted points.
The corrected series -- NB32-A to NB36-widened to NB56b-s050, each one an affine correction of the
one before on at least one isoform -- is drawn connected, so the bend is visible where correcting
stopped buying R-squared. The s = 1.0 candidate is built and unsent, which is the gap at the end of
the CYP1A2 and CYP2C9 chains. **On CYP2D6 the chain collapses to a single mark, because all three
submissions carry an identical CYP2D6 column** -- that series only ever corrected CYP1A2 and
CYP2C9.

**A darker grey mark is several submissions at the same point, not a different series.** Board
scores are identical for identical bytes, so every re-send and every submission that carries an
isoform's column over unchanged lands exactly on top of its predecessor: the dark mark on CYP1A2
is five submissions (10c, NB16, NB19, NB27-widened and 27b), and CYP2C9 has one of four and one of
two.

**The bare number at each chain endpoint is its leave-one-out residual from the trend**: refit
without that point, then measure it against the line the rest of the history sets. That divergence
is what the figure exists to show. It is +0.050 on CYP1A2 -- the endpoint sits well off a trend it
played no part in setting -- and +0.016 on CYP2C9. The in-sample residual would read +0.023 on
CYP1A2, because that endpoint has high leverage (board R-squared 0.51 against a cluster at 0.11 to
0.31) and drags the line towards itself.

For reference, the record carries an externally supplied field fit, macro
``ST-RAE = 0.695 - 0.210 x R2``, stated as valid for R2 in [0.40, 0.50]. It is not drawn here: it
is a macro relation from outside this project and these panels are per isoform.
"""

SLUG = "strae_vs_r2"
SIZES = {"thesis": (6.3, 3.6), "slides": (9.0, 5.0), "web": (8.0, 4.2)}

#: Only these two isoforms have a chain worth drawing. On CYP2D6 all three submissions carry an
#: IDENTICAL column, and on CYP3A4 two of the three coincide, so a "chain" there is one or two
#: marks on top of each other -- a panel showing nothing the caption cannot say in a clause.
PANEL_ISOFORMS = ("CYP1A2", "CYP2C9")

#: Each is an affine correction of its predecessor on at least one isoform.
CORRECTED_SERIES = ["NB32-A", "NB36-widened", "NB56b-s050"]


def variant(target: Target) -> str:
    return f"two panels ({', '.join(PANEL_ISOFORMS)} -- the only isoforms with a chain)"


def draw(target: Target):
    st, r2 = wide("ST-RAE"), wide("R²")
    assert list(st.index) == list(r2.index)
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        axes = fig.subplots(1, len(PANEL_ISOFORMS))
        for ax, iso in zip(axes, PANEL_ISOFORMS):
            x, y = r2[iso].to_numpy(float), st[iso].to_numpy(float)
            colour = ISOFORM_COLOURS[iso]

            m, c = np.polyfit(x, y, 1)                       # trend line over plotted points only
            xs = np.array([x.min(), x.max()])
            ax.plot(xs, m * xs + c, color=REFERENCE, linewidth=target.base_size * 0.07,
                    zorder=1)
            # LEAVE-ONE-OUT residual for the chain's endpoint: refit without that point, then
            # measure it against the trend the rest of the history establishes. Not the in-sample
            # residual -- CYP1A2's endpoint is a high-leverage point (R-squared 0.51 against a
            # cluster at 0.11-0.31), so including it drags the line towards itself and understates
            # the divergence by more than half (+0.023 in sample against +0.050 left out).
            def loo_residual(k: int) -> float:
                keep = [i for i in range(len(x)) if i != k]
                mk, ck = np.polyfit(x[keep], y[keep], 1)
                return float(y[k] - (mk * x[k] + ck))

            # Rule 2: faceted by isoform, so background submissions are NEUTRAL and only the
            # emphasised series -- the correction chain -- takes the isoform colour.
            ax.scatter(x, y, s=target.base_size ** 2 * 0.22, color=NEUTRAL, alpha=0.5,
                       edgecolor="none", zorder=3)

            idx = [list(st.index).index(s) for s in CORRECTED_SERIES]
            ax.plot(x[idx], y[idx], color=colour, linewidth=target.base_size * 0.13,
                    marker="o", markersize=target.base_size * 0.40, zorder=4)
            # Label only where the chain actually moves. On CYP2D6 all three submissions carry
            # the IDENTICAL column (0.7502) and on CYP3A4 two of three do, so the chain collapses
            # to a single mark there -- which is itself the finding: this series only ever
            # corrected CYP1A2 and CYP2C9.
            # Deduplicate on the board's own 4-decimal precision: where two submissions in the
            # chain carry an identical column, they are one mark and get one label.
            keep, seen = [], set()
            for k, name in zip(idx, CORRECTED_SERIES):
                key = (round(x[k], 4), round(y[k], 4))
                if key in seen:
                    continue
                seen.add(key)
                keep.append((k, name))
            # Surviving labels can still sit on top of each other where the chain barely moves,
            # so their text is spread vertically while the marks stay put.
            # Short labels, placed per label rather than at a fixed offset, with the trend line
            # avoided. The residual is the bare number; what it means is in the caption, and the
            # chain-collapse explanation moved there entirely rather than sitting over the cluster
            # it describes.
            place_labels(ax, [x[k] for k, _ in keep], [y[k] for k, _ in keep],
                         [f"{nm}  {loo_residual(k):+.3f}" if k == idx[-1] else nm
                          for k, nm in keep],
                         colours=TEXT, avoid=[(xs, m * xs + c)],
                         radius_pt=target.base_size * 0.80)

            pad = (x.max() - x.min())
            ax.set_xlim(x.min() - pad * 0.10, x.max() + pad * 0.30)
            ax.set_title(iso, color=colour)
            ax.grid(True, zorder=0)
            ax.set_axisbelow(True)
        for ax in axes:
            ax.set_xlabel("Board R-squared")
        axes[0].set_ylabel("Board ST-RAE")

        # The figure had no legend, and three of its four mark types were undecodable: the neutral
        # background points, the darker final submission among them, the chain, and the fit line.
        handles = [
            Line2D([], [], linestyle="none", marker="o", color=NEUTRAL, alpha=0.5,
                   markersize=target.base_size * 0.40,
                   label="other submissions"),
            Line2D([], [], linestyle="-", marker="o", color=ISOFORM_COLOURS[PANEL_ISOFORMS[0]],
                   markersize=target.base_size * 0.40,
                   label="correction chain"),
            Line2D([], [], linestyle="-", color=REFERENCE,
                   label="least-squares trend"),
        ]
        # "outside lower center" is the one placement constrained layout actually reserves room
        # for; an anchored fig.legend below the axes renders off the page.
        fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False,
                   handlelength=1.6, columnspacing=1.4, labelcolor=TEXT, fontsize="small")
        figure_title(fig, target, "ST-RAE tracks R-squared, until the corrected chain leaves it")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
