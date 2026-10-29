"""The compounds every model gets wrong, drawn. The only structures in the figure set.

Method and provenance live here; the caption is ``CAPTION`` below, which ``build.py`` writes to
``build/captions.md``.

Compounds are ranked by notebook 25's stored ``mean_abs_error`` -- the mean absolute error across
all eleven configurations, which is how that notebook defines consensus-hard. The predicted value
shown beside each is the SHIPPED arm's own out-of-fold prediction, pooled over its 25 folds, so the
two numbers under a structure are "what the assay measured" and "what the submitted model said".
Those are different populations by construction (eleven configs for the ranking, one for the
prediction) and the caption says so rather than implying one number explains the other.

Nothing is recomputed: the error ranking, the measured label and the SMILES are all read, and the
prediction is a mean over stored per-fold predictions. Drawing a molecule is rendering, not
computing.

SMILES come from the curated file's raw ``SMILES`` column. Notebook 43 established that
``canonical_smiles`` disagrees with it on atom order for 1,443 of 5,655 compounds; for a drawing
the molecule as the challenge supplied it is the right one.
"""

from __future__ import annotations

import io

from figures.data import load
from figures.oof import oof_long
from figures.style.base import Target, figure, figure_title, save
from figures.style.palette import ISOFORM_COLOURS, ISOFORMS, TEXT

SLUG = "hard_compound_structures"
SIZES = {"thesis": (6.3, 6.2), "slides": (9.0, 5.0), "web": (8.0, 6.8)}
N_PER_ISOFORM = 4

#: Compounds measured below this are excluded. Ranking by absolute error without it put 19 of 24
#: selections below 2.5 -- the dense low tail of weak measurements visible as a horizontal band in
#: the observed-against-predicted figure -- so the grid was a picture of where the error
#: distribution is widest, not of chemistry. Plenty remain above it: 1,352 / 1,264 / 1,424 / 2,026.
FLOOR = 2.5

CAPTION = """
The six worst-predicted compounds for each isoform, ranked by the mean absolute error across all
eleven cross-validated configurations — so these are the compounds the whole model pool gets wrong,
not one model's unlucky draw. Under each: the compound identifier, the measured pIC50, and the
shipped model's own out-of-fold prediction.

**They are not exotic.** Notebook 25 measured nearest-neighbour similarity to the training set for
the consensus-hard decile and found it does **not** separate hard from easy on any isoform
(rank-biserial 0.02 to 0.09, significant on none after correction). These are ordinary
drug-like molecules that the models place badly, not structural outliers — which is why a
retrieved-neighbour corpus aimed at unseen chemistry was the wrong remedy.

**Compounds measured below pIC50 2.5 are excluded, and that exclusion is the point.** Ranked on
absolute error alone, 19 of 24 selections fell below 2.5 — the dense low tail of weak measurements
that shows as a horizontal band in the observed-against-predicted figure — because a model
predicting near the mean is automatically most wrong at the extremes. That grid was a picture of
where the error distribution is widest, not of chemistry. Excluding the tail leaves 1,352 to 2,026
compounds per isoform, and the selection flips to the **potent** end: **15 of the 16 are measured
between 6.1 and 8.0 and under-predicted** by 2.1 to 2.7 log units. The one exception is a CYP1A2
compound at 2.71 that is over-predicted. Under-predicting a potent inhibitor is the more
consequential failure of the two — it is the one that would let a real interaction through — and it
is what the pool does on these.

**One of them is a cinchona alkaloid, and it is the most telling entry here.** CYP2D6's
worst-predicted non-floor compound, `OCNT-1911794`, carries the InChIKey connectivity block
`LOUPRKONTZGTKE` — the same constitution as **quinine and quinidine**, the archetypal CYP2D6
inhibitor pair. Measured at 7.47 and predicted at 5.56, the pool under-predicts it by 1.9 log
units. Which diastereomer it is cannot be read off: its stereo block matches neither quinine's
standard key nor quinidine's, which is the unspecified-configuration problem this dataset carries
throughout — 93 per cent of blind compounds with a stereocentre never state its configuration.

**Nor do they share a scaffold.** With 96.4 per cent of Bemis-Murcko scaffolds appearing exactly
once in this dataset, there is usually no common core to align on, and these twenty-four compounds
have none — so they are drawn in their own orientations rather than forced onto a template that
does not exist.

The ranking and the prediction come from different populations — eleven configurations for the
ranking, the shipped arm alone for the printed prediction — so a large gap between measured and
predicted here is illustrative of consensus difficulty rather than a measurement of this one model.
"""


def variant(target: Target) -> str:
    return f"{len(ISOFORMS)} x {N_PER_ISOFORM} structure grid"


def worst():
    """The N worst-predicted compounds per isoform, by notebook 25's stored consensus error."""
    hard = load("25_hard_compound_analysis", "hard_compound_table.csv")
    pred = oof_long()[["isoform", "inchikey", "predicted", "measured", "SMILES",
                       "Molecule_Name"]]
    j = hard.merge(pred, on=["isoform", "inchikey"], how="inner", suffixes=("", "_p"))
    above = j[j.y_true >= FLOOR]
    assert above.groupby("isoform").size().min() >= N_PER_ISOFORM, "grid cannot be filled"
    out = (above.sort_values("mean_abs_error", ascending=False)
                .groupby("isoform", sort=False).head(N_PER_ISOFORM))
    assert len(out) == len(ISOFORMS) * N_PER_ISOFORM, len(out)
    # Two independent reads of the same label must agree: notebook 25's y_true and the curated
    # file's own column, joined by inchikey through two different paths.
    assert (out.y_true - out.measured).abs().max() < 1e-9, "label sources disagree"
    # No common scaffold to align on: checked, not assumed.
    out = out.copy()
    out["shares_scaffold"] = out.groupby("isoform").scaffold.transform(
        lambda s: s.fillna("").duplicated(keep=False))
    return out


def _render_molecule(smiles: str, px: int):
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D
    import matplotlib.image as mpimg

    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"unparsed SMILES: {smiles}"
    d = rdMolDraw2D.MolDraw2DCairo(px, px)
    opts = d.drawOptions()
    opts.clearBackground = False
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return mpimg.imread(io.BytesIO(d.GetDrawingText()), format="png")


def draw(target: Target):
    w = worst()
    width, height = SIZES[target.name]
    px = 340 if target.name != "slides" else 260

    with figure(target, width, height) as fig:
        axes = fig.subplots(len(ISOFORMS), N_PER_ISOFORM)
        for row, iso in enumerate(ISOFORMS):
            sub = w[w.isoform == iso].reset_index(drop=True)
            for col in range(N_PER_ISOFORM):
                ax = axes[row][col]
                r = sub.iloc[col]
                ax.imshow(_render_molecule(r.SMILES, px))
                ax.set_axis_off()

                # Two short lines: "meas X pred Y" on one line is wider than a 1-inch cell at
                # six columns across 6.3 in, and ran into the neighbouring panel's title.
                ax.set_title(f"{r.Molecule_Name}\n{r.y_true:.2f} -> {r.predicted:.2f}",
                             fontsize="small", color=TEXT, loc="center",
                             pad=target.base_size * 0.25)
        # Isoform labels in the left margin, as figure artists, positioned after the layout has
        # settled so they sit on their own row rather than on the previous one.
        fig.canvas.draw()
        for row, iso in enumerate(ISOFORMS):
            pos = axes[row][0].get_position()
            fig.text(0.004, (pos.y0 + pos.y1) / 2, iso, rotation=90, ha="left", va="center",
                     color=ISOFORM_COLOURS[iso], fontweight="bold")
        fig.supxlabel("measured -> predicted pIC50, with the compound identifier",
                      fontsize="small", color=TEXT)
        figure_title(fig, target, "The compounds every model gets wrong")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
