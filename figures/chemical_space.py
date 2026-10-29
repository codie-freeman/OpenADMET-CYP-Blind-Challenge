"""Chemical space: training against the blind 750, in two descriptors.

Method and provenance live here; the caption is ``CAPTION`` below.

Deliberately a plain two-descriptor scatter rather than a t-SNE or UMAP projection. With 96.4 per
cent of Bemis-Murcko scaffolds appearing exactly once in this library, a neighbour-embedding has
very little local structure to preserve and its axes would carry no interpretable units, whereas
these two do.

Both descriptors are read: molecular weight from notebook 41's ``per_compound_audit.csv`` (which
also carries the train/blind flag) and Crippen logP from notebook 17's ``charge_logp_features.csv``,
joined by InChIKey. Both cover all 5,655 compounds.
"""

from __future__ import annotations

from figures.data import load
from figures.style.base import Target, figure, figure_title, save
from figures.style.palette import CONTRAST, NEUTRAL, TEXT

SLUG = "chemical_space"
SIZES = {"thesis": (6.3, 3.4), "slides": (9.0, 4.8), "web": (8.0, 4.0)}

CAPTION = """
The training set and the blind test set in two interpretable descriptors, molecular weight against
Crippen logP, with marginal distributions.

**The blind 750 sit inside the training distribution, not beside it.** They were built by taking
the most potent known hits for three of the four isoforms and purchasing catalogue near-neighbours
of them, so they are drawn from the same chemical space by construction — which is the single most
important fact for reading every cross-validation result in this work. A deliberately
cluster-disjoint held-out fold is a *harder* test than the blind set, and that is why most
cross-validated scores transferred to the board in this project's favour rather than against it.

The molecular-weight axis is clipped at the 99.5th percentile; the excluded count is stated on the
axis label. Uncapped, a handful of outliers above 700 Da stretched the axis to 1,000 and compressed
the real cloud — which sits between about 200 and 500 Da — into the left third of the panel.

Plotted as two real descriptors rather than a neighbour embedding: with 96.4 per cent of scaffolds
appearing once, there is little local neighbourhood structure for a t-SNE or UMAP to preserve, and
its axes would mean nothing.
"""


def variant(target: Target) -> str:
    return "scatter with marginals (train vs blind)"


def table():
    audit = load("41_dataset_audit_3d", "per_compound_audit.csv")[["inchikey", "set", "mol_wt"]]
    logp = load("17_protonation_charge_check", "charge_logp_features.csv")
    col = next(c for c in logp.columns if "logp" in c.lower())
    j = audit.merge(logp[["inchikey", col]].rename(columns={col: "logp"}),
                    on="inchikey", how="inner").dropna(subset=["mol_wt", "logp"])
    assert set(j["set"].unique()) == {"train", "blind"}, sorted(j["set"].unique())
    return j


def draw(target: Target):
    t = table()
    # CONTRAST, not an isoform colour: this comparison is train-against-blind, not about
    # isoform identity, so borrowing CYP2D6's orange would assert something the figure is not.
    groups = {"train": NEUTRAL, "blind": CONTRAST}
    # The real cloud sits between about 200 and 500 Da; a handful of outliers stretched the axis
    # to 1000 and squeezed everything into the left third. Clipped, with the excluded count said.
    mw_cap = float(t.mol_wt.quantile(0.995))
    n_excluded = int((t.mol_wt > mw_cap).sum())
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        gs = fig.add_gridspec(2, 2, width_ratios=(5, 1), height_ratios=(1, 4),
                              wspace=0.04, hspace=0.04)
        ax = fig.add_subplot(gs[1, 0])
        axx = fig.add_subplot(gs[0, 0], sharex=ax)
        axy = fig.add_subplot(gs[1, 1], sharey=ax)

        for name, colour in groups.items():
            g = t[t["set"] == name]
            ax.scatter(g.mol_wt, g.logp, s=target.base_size ** 2 * (0.04 if name == "train" else 0.07),
                       color=colour, alpha=0.22 if name == "train" else 0.55,
                       edgecolor="none", zorder=3 if name == "train" else 4,
                       label=f"{name} (n = {len(g):,})")
            axx.hist(g.mol_wt, bins=50, density=True, color=colour,
                     alpha=0.45, zorder=3, histtype="stepfilled")
            axy.hist(g.logp, bins=50, density=True, color=colour, alpha=0.45,
                     orientation="horizontal", zorder=3, histtype="stepfilled")

        ax.set_xlim(right=mw_cap)
        axx.set_xlim(right=mw_cap)
        ax.set_xlabel(f"Molecular weight (Da), axis clipped at {mw_cap:.0f}; "
                      f"{n_excluded} compounds above it not shown")
        ax.set_ylabel("Crippen logP")
        ax.legend(loc="lower right", frameon=False, labelcolor=TEXT, handletextpad=0.3)
        ax.grid(True, zorder=0); ax.set_axisbelow(True)
        for a in (axx, axy):
            a.set_axis_off()
        figure_title(fig, target, "The blind set sits inside the training space")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
