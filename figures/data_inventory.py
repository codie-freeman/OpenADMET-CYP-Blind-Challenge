"""What the models were trained on and what they were scored on. Opening figure.

Counts only -- a count is not a computation. Training labels from
``data/processed/train_inhibition_curated.csv``, cross-checked against notebook 23's own saved
per-isoform counts; the blind set from ``data/processed/test_blinded_curated.csv``.
"""

from __future__ import annotations

import pandas as pd

from figures.data import load
from figures.style.base import repo_root, Target, figure, figure_title, save
from figures.style.palette import ISOFORM_COLOURS, ISOFORMS, REFERENCE, TEXT

SLUG = "data_inventory"
SIZES = {"thesis": (6.3, 2.8), "slides": (9.0, 4.4), "web": (8.0, 3.4)}

CAPTION = """
The data this project had. **4,905 unique training compounds, but only 6,525 measurements across
the four isoforms** — so each isoform was trained on between 1,285 and 2,335 labels, not on the
whole set, and every compound carries a label for a median of one isoform rather than four. The
blind test set is 750 compounds, scored on all four.

Two consequences run through everything that follows. **The label matrix is sparse** — between 26
and 48 per cent of training compounds are labelled for any given isoform — which is why the models
are multitask, why a change to one isoform's labels moves the others' predictions, and why
single-isoform screens are harder to power than their compound counts suggest. And **CYP3A4 has
nearly twice CYP2C9's labels**, which is part of why it is consistently the best-predicted isoform
and CYP2D6, despite having more labels than CYP2C9, the worst.
"""


def counts() -> pd.DataFrame:
    root = repo_root()
    train = pd.read_csv(root / "data" / "processed" / "train_inhibition_curated.csv",
                        float_precision="round_trip")
    blind = pd.read_csv(root / "data" / "processed" / "test_blinded_curated.csv",
                        float_precision="round_trip")
    rows = [{"isoform": i,
             "train": int(train[f"{i}_pIC50_direct_inhibition"].notna().sum()),
             "blind": len(blind)} for i in ISOFORMS]
    df = pd.DataFrame(rows)
    df["share"] = df.train / len(train)

    # Cross-check against notebook 23's own saved counts rather than trusting this read.
    # NOT n_already_scored_same_isoform, which is a per-source-file overlap count and differs
    # (CYP3A4 reads 1,805 there). existing_scored_pool is the isoform's own scored-label count and
    # is constant across that table's three source files.
    ref = load("23_compound_pool_overlap", "pool_overlap_counts.csv")
    for _, r in df.iterrows():
        want = ref.loc[ref.isoform == r.isoform, "existing_scored_pool"].unique()
        assert len(want) == 1 and int(want[0]) == r.train, \
            f"{r.isoform}: {r.train} here vs {want} in notebook 23"
    df.attrs["n_train_compounds"] = len(train)
    df.attrs["n_blind"] = len(blind)
    return df


def variant(target: Target) -> str:
    return "single-panel (training labels per isoform against the blind set)"


def draw(target: Target):
    df = counts()
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        ax = fig.add_subplot(111)
        y = range(len(df))
        ax.barh(y, df.train, color=[ISOFORM_COLOURS[i] for i in df.isoform],
                height=0.6, zorder=3)
        for i, r in df.iterrows():
            ax.text(r.train + df.train.max() * 0.015, i,
                    f"{r.train:,}  ({r.share:.0%} of the {df.attrs['n_train_compounds']:,} compounds)",
                    va="center", ha="left", fontsize="small", color=TEXT)

        ax.axvline(df.attrs["n_blind"], color=REFERENCE, linestyle="--",
                   linewidth=target.base_size * 0.09, zorder=4)
        ax.text(df.attrs["n_blind"], -0.85, f" blind set: {df.attrs['n_blind']} compounds,\n"
                                            " all four isoforms",
                va="top", ha="left", fontsize="small", color=TEXT)

        ax.set_yticks(list(y)); ax.set_yticklabels(df.isoform)
        ax.set_ylim(len(df) - 0.4, -1.1)
        ax.set_xlim(0, df.train.max() * 1.95)
        ax.set_xlabel("Compounds with a measured pIC50")
        ax.xaxis.grid(True, zorder=0)
        ax.set_axisbelow(True)
        figure_title(fig, target, "Sparse labels, four isoforms, 750 blind compounds")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
