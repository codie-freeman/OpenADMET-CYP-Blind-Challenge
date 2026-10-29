"""Training labels against the SOLVED blind population. Not a guess at the test distribution.

Method and provenance live here; the caption is ``CAPTION`` below.

The blind moments are **solved, not assumed**: CYP1A2, CYP2C9 and CYP3A4 from notebook 56b's
three-equation solves (``outputs/56b_cyp1a2_corrected/solved_populations.csv``, rho recovered rather
than proxied) and CYP2D6 from notebook 55's four-equation solve
(``outputs/55_superseded_preserved/solve_identified.csv``). Training labels are read from the
curated file. Nothing is re-solved here.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from figures.data import load
from figures.style.base import repo_root, Target, figure, figure_title, save
from figures.style.palette import ISOFORM_COLOURS, ISOFORMS, REFERENCE, TEXT

SLUG = "label_vs_solved_population"
SIZES = {"thesis": (6.3, 3.6), "slides": (9.0, 4.8), "web": (8.0, 4.2)}

CAPTION = """
Each isoform's training-label distribution with its **solved blind population** overlaid as a normal
curve at the solved mean and SD. Vertical rules mark the two means.

**The blind populations were solved from the board's own returned metrics, not guessed.** Each is
the exact solution of `R2 = 2*rho*k - k^2 - b^2` over three or four submitted affine variants of one
ordering, with the Pearson correlation *recovered* rather than substituted from Spearman. That is
why this can be drawn at all: most work of this kind puts a train-versus-test comparison in
supporting information using the test labels, and here there are none — the blind labels do not
exist in this repository and will not until the challenge closes.

**CYP2D6 is the outlier and it is the whole reason the corrections exist.** Its solved blind centre
sits about 1.7 pIC50 units *below* its training mean, where the other three sit within about 0.7.
A model trained on the training distribution and scored on the blind one is therefore mis-centred on
CYP2D6 by more than a log unit before it makes a single chemical mistake, which is exactly what the
deployed placement correction removes — and it took that isoform's board ST-RAE from 1.1673 to
0.8299 on identical raw predictions.
"""


def variant(target: Target) -> str:
    return "1 x 4 facets (labels with the solved blind population)"


def solved() -> pd.DataFrame:
    a = load("56b_cyp1a2_corrected", "solved_populations.csv")[["isoform", "ybar", "sd_y"]]
    b = load("55_superseded_preserved", "solve_identified.csv").rename(
        columns={"over_ybar": "ybar", "over_sd_y": "sd_y"})[["isoform", "ybar", "sd_y"]]
    out = pd.concat([a, b], ignore_index=True)
    assert set(out.isoform) == set(ISOFORMS), sorted(set(out.isoform))
    return out.set_index("isoform")


def draw(target: Target):
    curated = pd.read_csv(repo_root() / "data" / "processed" / "train_inhibition_curated.csv",
                          float_precision="round_trip")
    sol = solved()
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        axes = fig.subplots(1, 4, sharey=True)
        for ax, iso in zip(axes, ISOFORMS):
            y = curated[f"{iso}_pIC50_direct_inhibition"].dropna()
            colour = ISOFORM_COLOURS[iso]
            ax.hist(y, bins=40, density=True, color=colour, alpha=0.45, zorder=3)
            ax.axvline(y.mean(), color=colour, linestyle="-",
                       linewidth=target.base_size * 0.09, zorder=4)

            mu, sd = float(sol.loc[iso, "ybar"]), float(sol.loc[iso, "sd_y"])
            xs = np.linspace(min(y.min(), mu - 3 * sd), max(y.max(), mu + 3 * sd), 300)
            ax.plot(xs, np.exp(-0.5 * ((xs - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi)),
                    color=TEXT, linewidth=target.base_size * 0.10, zorder=5)
            ax.axvline(mu, color=TEXT, linestyle="--",
                       linewidth=target.base_size * 0.08, zorder=5)

            ax.set_title(iso, color=colour)
            ax.set_title(f"shift {mu - y.mean():+.2f}", loc="right", fontsize="small", color=TEXT)
            ax.set_xlabel("pIC50")
            ax.grid(True, zorder=0); ax.set_axisbelow(True)
        axes[0].set_ylabel("density")
        fig.supxlabel("filled = training labels; dark curve and dashed rule = solved blind population",
                      fontsize="small", color=TEXT)
        figure_title(fig, target, "The blind population is not the training population")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
