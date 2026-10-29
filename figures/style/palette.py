"""Colour and label definitions shared by every figure in this layer.

The palette is Okabe-Ito, chosen for colourblind safety. It is defined once, here.
No figure module defines a colour of its own.

Greyscale legibility is not a requirement for these targets, so colour carries
meaning rather than merely decorating.
"""

# --- isoform colours -------------------------------------------------------
# Okabe-Ito. One colour per isoform, used consistently across every figure so a
# reader who learns the mapping in one figure keeps it in the next.
ISOFORM_COLOURS = {
    "CYP1A2": "#0072B2",   # blue
    "CYP2C9": "#009E73",   # bluish green
    "CYP2D6": "#D55E00",   # vermillion
    "CYP3A4": "#CC79A7",   # reddish purple
}

#: Canonical isoform order: the order they appear in every result table in this repo.
ISOFORMS = tuple(ISOFORM_COLOURS)

# --- four jobs, four treatments --------------------------------------------
# Colour does different work in different figures, and branching from the
# isoform palette is legitimate where the JOB is different -- not where the
# figure merely looks crowded. The four jobs:
#
#   CATEGORICAL (which isoform)  -> the four-colour set below. IDENTITY ONLY.
#   SEQUENTIAL (magnitude)       -> ONE single-hue ramp, `sequential()`, the
#                                   same hue everywhere it appears, ALWAYS with
#                                   a colourbar. Deliberately not drawn from the
#                                   isoform palette, so a shade can never be
#                                   mistaken for an isoform.
#   BINARY STATE (yes/no)        -> filled and empty. No colour, no colourbar.
#                                   `pairwise_significance` is the reference.
#   EMPHASIS within a panel      -> the isoform colour as an accent against
#                                   NEUTRAL.
#
# Text is always TEXT; panel titles are always isoform-coloured.
#
# Within the categorical job, which marks take colour:
#
#   1. ONE SERIES IN THE PANEL -- marks TAKE the isoform colour. This holds
#      whether every isoform shares one set of axes or each has its own panel:
#      what matters is that nothing else in the panel needs distinguishing.
#      (fold_label_distributions, interval_width_vs_potency,
#       error_vs_prediction_correlation, and every multi-isoform axis.)
#   2. SEVERAL SERIES IN THE PANEL -- marks go NEUTRAL and exactly ONE accent
#      series takes the isoform colour, because colour is then the only thing
#      separating that series from the rest.
#      (model_comparison_cld: the best arm. hard_compound_similarity: the hard
#       decile. strae_vs_r2: the correction chain.)
#      Panel titles stay isoform-coloured either way.
#   3. ALL TEXT IS `TEXT`. The only coloured text in this project is a panel
#      title. No tinted annotations, point labels or series labels.
#   4. ERROR BARS take their own mark's colour at reduced weight: thinner than
#      the marker edge, alpha about 0.6, small or no caps.
#   5. REFERENCE GEOMETRY is always `REFERENCE` -- never coloured, and never
#      heavier than the data it is a reference for.

#: Every annotation, axis label and tick. Also set as `text.color` in each style file.
TEXT = "#333333"
#: Marks where isoform is not the variable being encoded.
NEUTRAL = "#595959"
#: Diagonals, fit lines, leader lines, any reference geometry.
REFERENCE = "#B0B0B0"
#: Grid lines.
GRID = "#E0E0E0"

#: One accent for a BINARY comparison that is not about isoform identity --
#: train against blind, scaffold against Butina. Deliberately outside the
#: four-colour isoform set. A dark violet: its minimum RGB distance to any of the
#: four isoform colours is 0.48, against 0.30 for the dark green first tried,
#: which sat close enough to CYP2C9 to be mistaken for it.
CONTRAST = "#762A83"

#: Deprecated alias for REFERENCE. Kept so nothing breaks mid-edit; do not use.
NEUTRAL_LIGHT = REFERENCE


# --- modelling-arm display names -------------------------------------------
# The raw arm names run to 22 characters (``ecfp4_narrow__lightgbm``), which
# consumes most of a 3.15 in thesis panel. These abbreviations are what the
# figures display; any caption using them must point a reader at the methods,
# where the mapping below belongs verbatim.
#
#   RAW NAME                  DISPLAYED AS     WHAT IT IS
#   ------------------------  ---------------  ----------------------------------
#   plain                     plain            chemprop_chemeleoninit: Chemprop
#                                              multitask, CheMeleon encoder, four
#                                              regression heads, MSE loss. The
#                                              reference arm; submitted as 10c.
#   deadzone                  dead-zone        as plain, but fitted to
#                                              clip(OOF prediction, conf_low,
#                                              conf_high) under MAE loss.
#   mae_only                  MAE loss         as plain with the loss switched to
#                                              MAE and the raw label kept. The
#                                              leak-free half of dead-zone;
#                                              submitted as NB47-mae.
#   aid                       AID heads        as plain plus eight auxiliary heads
#                                              from PubChem AID 1851. Submitted as
#                                              NB30.
#   single                    single-task      one model per isoform rather than
#                                              one multitask model.
#   ecfp4_narrow__rf          ECFP4 / RF       random forest on the ecfp4_narrow
#                                              feature block.
#   ecfp4_narrow__lightgbm    ECFP4 / LGBM     LightGBM, same block.
#   ecfp4_narrow__xgboost     ECFP4 / XGB      XGBoost, same block.
#   mordred_pca__rf           Mordred / RF     random forest on PCA-reduced
#                                              Mordred 2D descriptors.
#   mordred_pca__lightgbm     Mordred / LGBM   LightGBM, same block.
#   mordred_pca__xgboost      Mordred / XGB    XGBoost, same block.
ARM_LABELS = {
    "plain": "plain",
    "deadzone": "dead-zone",
    "mae_only": "MAE loss",
    "aid": "AID heads",
    "single": "single-task",
    "ecfp4_narrow__rf": "ECFP4 / RF",
    "ecfp4_narrow__lightgbm": "ECFP4 / LGBM",
    "ecfp4_narrow__xgboost": "ECFP4 / XGB",
    "mordred_pca__rf": "Mordred / RF",
    "mordred_pca__lightgbm": "Mordred / LGBM",
    "mordred_pca__xgboost": "Mordred / XGB",
}


def arm_label(arm: str) -> str:
    """Display name for a modelling arm, falling back to the raw name."""
    return ARM_LABELS.get(arm, arm)


# --- sequential ramp for heatmaps ------------------------------------------
# Derived from the palette rather than introducing a new hue: white through
# CYP1A2's blue. Heatmap panels orient it so the DARK end is always "better",
# which means reversing it for metrics where lower is better.
#: The one sequential hue, used wherever magnitude is encoded. Deliberately a
#: desaturated slate rather than any isoform colour, so a shade on a ramp can
#: never be read as an isoform identity.
SEQUENTIAL_STOPS = ["#FFFFFF", "#DCE3EA", "#A9BBCC", "#6D88A6", "#3C5C80", "#1F3550"]


def sequential(reverse: bool = False):
    """The project's single sequential ramp. Always pair it with a colourbar.

    ``reverse=True`` puts the dark end at LOW values, which is what a
    lower-is-better metric needs so that "dark = good" holds in every panel of a
    multi-metric figure.
    """
    from matplotlib.colors import LinearSegmentedColormap
    stops = SEQUENTIAL_STOPS[::-1] if reverse else SEQUENTIAL_STOPS
    return LinearSegmentedColormap.from_list("seq", stops)
