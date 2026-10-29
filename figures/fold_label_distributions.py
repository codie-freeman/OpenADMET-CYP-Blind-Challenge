"""pIC50 distribution per cluster-disjoint fold, faceted by isoform. Appendix.

Method and provenance live here; the caption is ``CAPTION`` below, which ``build.py``
writes to ``build/captions.md`` alongside the renders.
"""

from __future__ import annotations

import pandas as pd

from figures.data import load
from figures.style.base import repo_root, Target, figure, figure_title, save
from figures.style.palette import ISOFORM_COLOURS, ISOFORMS, NEUTRAL, REFERENCE, TEXT

CAPTION = """
Measured pIC50 in each of the 25 cluster-disjoint held-out folds, by isoform. Ash et al. recommend
showing this whenever a non-random split is used, because a cluster-disjoint partition can
concentrate activity in a fold in a way a random split cannot, and a fold whose label distribution
is unrepresentative will produce a fold score that is not comparable with the others.

Folds are ``repeat x fold`` from ``data/folds/cv_folds_butina_5x5.csv`` -- the frozen partition every
arm in the model comparison was trained and scored on -- joined to the curated labels by InChIKey.
Boxes are the interquartile range with the median marked; whiskers 1.5 x IQR; the dashed rule is
each isoform's own overall median. Per-isoform fold membership varies because 26 to 48 per cent of
compounds carry a label for any given isoform.

Note the fold columns in that file are named ``butina_repeat_0`` to ``butina_repeat_4``, not
``repeat_0`` to ``repeat_4``.
"""

SLUG = "fold_label_distributions"
SIZES = {"thesis": (6.3, 6.6), "slides": (9.0, 5.0), "web": (8.0, 7.2)}

N_REPEATS, N_FOLDS = 5, 5


def variant(target: Target) -> str:
    return f"4 x 1 facets ({N_REPEATS * N_FOLDS} folds each)"


def long_labels() -> pd.DataFrame:
    """One row per (fold, isoform, compound) with a measured label. A join, not a computation."""
    root = repo_root()
    folds = pd.read_csv(root / "data" / "folds" / "cv_folds_butina_5x5.csv",
                        float_precision="round_trip")
    curated = pd.read_csv(root / "data" / "processed" / "train_inhibition_curated.csv",
                          float_precision="round_trip")
    cols = [f"butina_repeat_{r}" for r in range(N_REPEATS)]
    assert set(cols) <= set(folds.columns), sorted(folds.columns)

    merged = folds.merge(curated, on="inchikey", how="inner", suffixes=("", "_c"))
    assert len(merged) == len(folds), f"{len(merged)} joined vs {len(folds)} fold rows"

    rows = []
    for r, col in enumerate(cols):
        for iso in ISOFORMS:
            y = merged[f"{iso}_pIC50_direct_inhibition"]
            m = y.notna()
            rows.append(pd.DataFrame({"repeat": r, "fold": merged.loc[m, col].astype(int),
                                      "isoform": iso, "pIC50": y[m]}))
    out = pd.concat(rows, ignore_index=True)
    out["cell"] = out.repeat.astype(str) + "." + out.fold.astype(str)
    assert out.cell.nunique() == N_REPEATS * N_FOLDS, out.cell.nunique()
    return out


def draw(target: Target):
    df = long_labels()
    cells = sorted(df.cell.unique(), key=lambda c: tuple(int(x) for x in c.split(".")))
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        axes = fig.subplots(len(ISOFORMS), 1, sharex=True)
        for ax, iso in zip(axes, ISOFORMS):
            sub = df[df.isoform == iso]
            colour = ISOFORM_COLOURS[iso]
            data = [sub.loc[sub.cell == c, "pIC50"].to_numpy() for c in cells]
            bp = ax.boxplot(data, positions=range(len(cells)), widths=0.62,
                            patch_artist=True, showfliers=False, zorder=3)
            # Rule 1: one series per panel, so the marks take the isoform colour.
            for box in bp["boxes"]:
                box.set(facecolor=colour, alpha=0.35, edgecolor=colour,
                        linewidth=target.base_size * 0.06)
            for part in ("whiskers", "caps", "medians"):
                for art in bp[part]:
                    art.set(color=colour, linewidth=target.base_size * 0.07)

            ax.axhline(sub.pIC50.median(), color=REFERENCE, linestyle="--",
                       linewidth=target.base_size * 0.07, zorder=1)
            ax.set_ylabel(iso, color=colour)
            # In the RIGHT TITLE SLOT, outside the data area: inside the axes this sat behind
            # the whiskers in every panel. "n" would read as a sample size and is not one --
            # every labelled compound appears once per repeat -- so it says what it is.
            n_compounds = len(sub) // N_REPEATS
            ax.set_title(f"{n_compounds:,} compounds x {N_REPEATS} repeats   "
                         f"median {sub.pIC50.median():.2f}",
                         loc="right", fontsize="small", color=TEXT)
            ax.yaxis.grid(True, zorder=0); ax.set_axisbelow(True)

        axes[-1].set_xticks(range(len(cells)))
        axes[-1].set_xticklabels(cells, rotation=90, fontsize="small")
        fig.supxlabel("repeat . fold (cluster-disjoint partition)")
        fig.supylabel("Measured pIC50")
        figure_title(fig, target, "Label distribution is even across the 25 folds")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
