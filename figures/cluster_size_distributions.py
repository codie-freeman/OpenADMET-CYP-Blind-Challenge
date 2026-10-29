"""Scaffold clusters against Butina clusters. The evidence for the split design.

Method and provenance live here; the caption is ``CAPTION`` below.

Butina cluster sizes are notebook 33's stored ``cluster_size_distribution.csv``. Scaffold cluster
sizes are a COUNT of notebook 25's stored ``scaffold`` column over the 4,905 unique compounds --
a count of a stored column, not a re-derivation, and it reproduces notebook 02's 96.4 per cent
singleton figure exactly (4,363 of 4,526).
"""

from __future__ import annotations

import pandas as pd

from figures.data import load
from figures.style.base import Target, figure, figure_title, save
from figures.style.palette import CONTRAST, NEUTRAL, TEXT

SLUG = "cluster_size_distributions"
SIZES = {"thesis": (6.3, 2.8), "slides": (9.0, 4.4), "web": (8.0, 3.4)}

CAPTION = """
How the 4,905 training compounds group, two ways. **Bemis-Murcko scaffolds: 4,526 groups, 96.4 per
cent of them containing exactly one compound.** Butina clusters at Morgan radius 3 and Tanimoto
distance 0.65: 2,617 groups, 56.0 per cent singleton.

This pair of numbers is why the split design is what it is. A scaffold-based split was the obvious
first choice and is not available: with 96 per cent of scaffolds appearing once, scaffold groups
carry almost no information about which compounds resemble each other, and a scaffold split
degenerates to a random one. Butina clustering on fingerprint similarity does group this library —
**70.1 per cent of compounds sit in a non-singleton cluster** — which is what makes a genuinely
cluster-disjoint partition possible, and that partition is what every model comparison in this work
is run on.

Both panels share one x-range, cut at the largest Butina cluster; the scaffold panel's sparse tail
beyond it is reported as a count rather than stretching both axes to reach it.

Note the two percentages count different things and are easy to conflate: 96.4 and 56.0 are shares
of *groups* that are singletons, while 70.1 is the share of *compounds* in a non-singleton group.
An earlier version of this project's own notes quoted the cluster figure as "44 per cent of
compounds", which was the share of clusters, not of compounds.
"""


def variant(target: Target) -> str:
    return "two panels (scaffold, Butina)"


def distributions():
    butina = load("33_butina_split", "cluster_size_distribution.csv")
    hard = load("25_hard_compound_analysis", "hard_compound_table.csv")
    scaff = (hard.drop_duplicates("inchikey").scaffold.fillna("").value_counts()
                 .value_counts().rename("n_clusters").rename_axis("cluster_size")
                 .reset_index().sort_values("cluster_size"))
    return {"Bemis-Murcko scaffold": scaff, "Butina (r=3, distance 0.65)": butina}


def draw(target: Target):
    dists = distributions()
    width, height = SIZES[target.name]

    # The shared limit covers the Butina distribution in full; the scaffold tail beyond it holds
    # a handful of groups and is reported as a count instead.
    shared_max = int(max(d.cluster_size.max() for name, d in dists.items() if "Butina" in name))

    with figure(target, width, height) as fig:
        axes = fig.subplots(1, 2, sharey=True, sharex=True)
        for ax, (name, d) in zip(axes, dists.items()):
            n_clusters = int(d.n_clusters.sum())
            n_compounds = int((d.cluster_size * d.n_clusters).sum())
            singleton = int(d.loc[d.cluster_size == 1, "n_clusters"].iloc[0])
            in_multi = n_compounds - singleton
            # NEUTRAL and CONTRAST, not isoform colours: this is scaffold-against-Butina, not
            # a statement about any isoform.
            colour = CONTRAST if "Butina" in name else NEUTRAL
            ax.bar(d.cluster_size, d.n_clusters, color=colour, width=0.8, zorder=3)
            ax.set_yscale("log")
            ax.set_title(name, color=TEXT)
            ax.set_title(f"{singleton / n_clusters:.1%} singleton", loc="right",
                         fontsize="small", color=TEXT)
            ax.text(0.97, 0.88, f"{n_clusters:,} groups\n"
                                f"{in_multi / n_compounds:.1%} of compounds\nin a group of 2+",
                    transform=ax.transAxes, ha="right", va="top", fontsize="small", color=TEXT)
            # One shared x-range so the two distributions can actually be compared. The scaffold
            # panel's tail runs to 50 and is almost empty; it is cut at the shared limit and the
            # count beyond is stated rather than stretching both axes to accommodate it.
            ax.set_xlim(0.3, shared_max + 0.7)
            beyond = int(d.loc[d.cluster_size > shared_max, "n_clusters"].sum())
            if beyond:
                ax.text(0.97, 0.70, f"+{beyond} group(s) larger\nthan {shared_max}",
                        transform=ax.transAxes, ha="right", va="top",
                        fontsize="small", color=TEXT)
            ax.set_xlabel("compounds in the group")
            ax.grid(True, zorder=0); ax.set_axisbelow(True)
        axes[0].set_ylabel("number of groups (log)")
        figure_title(fig, target, "Why the split is Butina and not scaffold")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
