"""The normal-theory rho proxy is near-exact for well-ranked isoforms and badly wrong for CYP2D6.

Every placement correction this project deployed needs a Pearson rho it cannot observe, and takes
one of two substitutes from the board's Spearman: Spearman used directly (notebook 16), or the
bivariate-normal relation ``rho = 2 * sin(pi * rho_S / 6)`` (notebook 38). Both were unverifiable
until an exactly solved population produced a RECOVERED rho to compare them against.

THE IDENTITY, AND THE VERIFICATION THAT LICENSES EVALUATING IT HERE
-------------------------------------------------------------------
``rho_bvn = 2 * sin(pi * rho_S / 6)`` is a closed form, not a fit. Notebook 55 stored the whole
comparison for **CYP2D6** (``outputs/55_superseded_preserved/rho_recovered_vs_proxy.csv``:
recovered 0.560016, board Spearman 0.4, proxy 0.4158234, bias 0.1441926), and that stored row is
**overplotted as an open mark and asserted against the evaluated identity** before anything else is
drawn. The other three isoforms' recovered rho come from notebook 56b's exact three-equation solves
(``outputs/56b_cyp1a2_corrected/solved_populations.csv``) and have no stored proxy, so the identity
is evaluated for them. See ``figures/README.md``.

**The file named in the brief, ``outputs/56b_cyp1a2_corrected/rho_proxy_bias.csv``, does not exist**
-- the data is split across the two files above, and this module joins them.

Each recovered rho belongs to a specific submitted ORDERING, so each is paired with that ordering's
own board Spearman: NB36-widened for the three notebook 56b solved, and the plain ordering (10c)
for CYP2D6, which is the ordering notebook 55 solved.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from figures.board import metrics
from figures.data import load
from figures.style.base import Target, figure, figure_title, save, spread_labels
from figures.style.palette import ISOFORM_COLOURS, REFERENCE, TEXT

SLUG = "rho_proxy_bias"
SIZES = {"thesis": (6.3, 3.4), "slides": (9.0, 4.8), "web": (8.0, 4.0)}

#: Which submitted ordering each recovered rho was solved from.
ORDERING_OF = {"CYP1A2": "NB36-widened", "CYP2C9": "NB36-widened",
               "CYP3A4": "NB36-widened", "CYP2D6": "10c"}

CAPTION = """
The bias in the Pearson rho substitute every deployed placement correction depends on, against the
board Spearman of the ordering it was taken from. Positive means the proxy **understates** the true
correlation. Four points, one per isoform; the open mark is the one value a notebook had already
stored, overplotted as a check on the identity evaluated for the other three.

**Above Spearman 0.75 the normal-theory proxy is near-exact — it *overstates* by 0.002 to 0.015 —
and at CYP2D6's 0.40 it *understates* by 0.144, an order of magnitude larger and in the opposite
direction.** The sign flip is worth keeping: the brief version of this result, that the proxy
"undershoots", is true only of CYP2D6. On the three isoforms that rank well the error is both tiny
and the other way.

So the substitution was harmless where models rank well and materially wrong on the one isoform
that does not — which is also the isoform this project corrected most often. Because rho sets the
shrinkage slope, understating it means every deployed CYP2D6 correction **shrank too hard**;
notebook 26 reached the same sign independently from an out-of-fold spread sweep.

The practical reading is that a proxy validated only where models rank well tells you nothing about
where they rank badly, and badly-ranked is exactly the case a placement correction is for.
"""


def variant(target: Target) -> str:
    return "single-panel (four isoforms)"


def _bvn(spearman):
    """rho = 2 sin(pi rho_S / 6) -- the bivariate-normal relation, notebook 38's refinement."""
    return 2.0 * np.sin(np.pi * np.asarray(spearman, dtype=float) / 6.0)


def table() -> pd.DataFrame:
    solved = load("56b_cyp1a2_corrected", "solved_populations.csv").set_index("isoform")
    stored = load("55_superseded_preserved", "rho_recovered_vs_proxy.csv").set_index("isoform")
    sp = metrics().query("metric == 'Spearman'").pivot(
        index="submission", columns="isoform", values="value")

    rows = []
    for iso, ordering in ORDERING_OF.items():
        recovered = (float(stored.loc[iso, "rho_recovered"]) if iso in stored.index
                     else float(solved.loc[iso, "rho"]))
        s = float(sp.loc[ordering, iso])
        rows.append({"isoform": iso, "ordering": ordering, "board_spearman": s,
                     "recovered": recovered, "proxy": float(_bvn(s)),
                     "stored": iso in stored.index})
    df = pd.DataFrame(rows)
    df["bias"] = df.proxy - df.recovered

    # The stored row is the verification: the evaluated identity must reproduce it.
    for _, r in df[df.stored].iterrows():
        for col, want in (("board_spearman", "board_spearman"), ("proxy", "rho_bvn_proxy"),
                          ("recovered", "rho_recovered")):
            assert abs(r[col] - float(stored.loc[r.isoform, want])) < 1e-6, \
                f"{r.isoform}/{col}: {r[col]} against stored {stored.loc[r.isoform, want]}"
        assert abs(-r.bias - float(stored.loc[r.isoform, "bias_vs_bvn"])) < 1e-6
    return df.sort_values("board_spearman").reset_index(drop=True)


def draw(target: Target):
    df = table()
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        ax = fig.add_subplot(111)
        ax.axhline(0, color=REFERENCE, linewidth=target.base_size * 0.08, zorder=1)
        ax.plot(df.board_spearman, -df.bias, color=REFERENCE,
                linewidth=target.base_size * 0.06, zorder=2)

        for _, r in df.iterrows():
            ax.plot(r.board_spearman, -r.bias, marker="o", color=ISOFORM_COLOURS[r.isoform],
                    markerfacecolor="white" if r.stored else ISOFORM_COLOURS[r.isoform],
                    markeredgewidth=target.base_size * 0.09,
                    markersize=target.base_size * (0.60 if r.stored else 0.48), zorder=4)

        # Three of the four sit within 0.015 of zero and within 0.05 of each other on x, so their
        # labels are stacked to the right of the cluster rather than placed at each mark.
        hi, lo = (-df.bias).max(), (-df.bias).min()
        label_y = spread_labels((-df.bias).tolist(), (hi - lo) * 0.17)
        # Short labels: the sign convention is on the y axis and the proxy-versus-recovered
        # numbers are in the caption, so repeating them here only overflows the axes.
        for (_, r), ly in zip(df.iterrows(), label_y):
            ax.annotate(f"{r.isoform}  {-r.bias:+.3f}",
                        (r.board_spearman + 0.014, ly), ha="left", va="center",
                        fontsize="small", color=TEXT)

        ax.set_xlabel("Board Spearman of the ordering the rho was solved from")
        ax.set_ylabel("Proxy error (positive = understates recovered rho)")
        ax.set_xlim(df.board_spearman.min() - 0.05, df.board_spearman.max() + 0.18)
        ax.set_ylim(-0.055, (-df.bias).max() * 1.35)
        ax.grid(True, zorder=0); ax.set_axisbelow(True)
        ax.text(0.02, 0.04, "open mark = the one value a notebook had already stored,\n"
                            "overplotted as a check on the evaluated identity",
                transform=ax.transAxes, ha="left", va="bottom",
                fontsize="small", color=TEXT)
        figure_title(fig, target, "The rho proxy fails exactly where it is needed")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
