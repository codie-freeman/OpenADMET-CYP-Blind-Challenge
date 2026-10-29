"""The full error-correlation matrix: shared ignorance across the whole model pool.

Method and provenance live here; the caption is ``CAPTION`` below. Read from notebook 24's
``error_correlation_long.csv``; nothing is recomputed. Sequential ramp with a colourbar, per the
project's colour rule -- magnitude is a different job from identity, so it does not use the isoform
palette.
"""

from __future__ import annotations

import numpy as np

from figures.data import load
from figures.style.base import Target, figure, figure_title, save
from figures.style.palette import ISOFORM_COLOURS, ISOFORMS, TEXT, sequential

SLUG = "error_correlation_matrix"
SIZES = {"thesis": (6.3, 6.0), "slides": (9.0, 5.0), "web": (8.0, 6.6)}

CAPTION = """
Residual (prediction minus truth) correlation between every pair of the eleven cross-validated
configurations, per isoform, lower triangle. Darker is more correlated.

**The pool shares its errors.** Correlations run roughly 0.6 to 0.97 and most of the matrix is
dark — these models mostly get the same compounds wrong, which is why nine successive attempts at
building a usefully diverse ensemble on this project came back pool-typical. What determines
whether combining models helps is error decorrelation, not prediction agreement (Dietterich, 2000),
and this is the quantity that matters shown in full rather than the two-isoform scatter.

**The structure is a representation structure, not a model-family one.** The pale cells cluster
where a graph-based model meets a descriptor-based one; within the tabular block, swapping random
forest for LightGBM for XGBoost barely changes which compounds are missed. That is the same finding
notebook 11b reached from the other direction — the real diversity axis is graph against tabular,
not the choice of head on top of one representation.
"""


def variant(target: Target) -> str:
    return "2 x 2 grid (4 isoforms, lower triangle)"


def matrices():
    t = load("24_nonnegative_stacking", "error_correlation_long.csv")
    assert set(t.isoform) == set(ISOFORMS), sorted(set(t.isoform))
    out = {}
    for iso, g in t.groupby("isoform"):
        m = g.pivot(index="config_1", columns="config_2", values="error_correlation")
        m = m.combine_first(m.T)                      # the file stores one triangle
        order = m.mean().sort_values(ascending=False).index.tolist()
        out[iso] = m.loc[order, order]
    return out


def short(name: str) -> str:
    return (name.replace("chemprop_chemeleoninit", "chemprop/CheMeleon")
                .replace("chemprop_randominit", "chemprop/random")
                .replace("chemeleon__", "CheMeleon/").replace("ecfp4_narrow__", "ECFP4/")
                .replace("mordred_pca__", "Mordred/").replace("lightgbm", "LGBM")
                .replace("xgboost", "XGB").replace("rf", "RF"))


def draw(target: Target):
    mats = matrices()
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        axes = fig.subplots(2, 2)
        vmin = min(float(np.nanmin(m.to_numpy())) for m in mats.values())
        for ax, iso in zip(axes.ravel(), ISOFORMS):
            m = mats[iso]
            v = m.to_numpy(dtype=float).copy()
            v[np.triu_indices_from(v)] = np.nan            # lower triangle only
            im = ax.imshow(v, cmap=sequential(), aspect="equal", vmin=vmin, vmax=1.0)
            ax.set_title(iso, color=ISOFORM_COLOURS[iso])
            n = len(m)
            ax.set_xticks(range(n - 1))
            ax.set_xticklabels([short(c) for c in m.index[:-1]], rotation=90, fontsize="small")
            ax.set_yticks(range(1, n))
            ax.set_yticklabels([short(c) for c in m.index[1:]], fontsize="small")
            ax.set_xlim(-0.5, n - 1.5); ax.set_ylim(n - 0.5, 0.5)
            ax.tick_params(length=0)
            for sp in ax.spines.values():
                sp.set_visible(False)

        cb = fig.colorbar(im, ax=axes, orientation="horizontal", fraction=0.035, pad=0.03)
        cb.set_label("residual correlation", color=TEXT, fontsize="small")
        cb.outline.set_visible(False)
        figure_title(fig, target, "The pool shares its errors")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
