"""Freeze notebook 54's pre-registration BEFORE the `interval_2d6` arm is trained.

This script REFUSES TO RUN if any `interval_2d6` prediction already exists on disk, so the
thresholds below provably cannot have been tuned to the arm's own result.

Every quantity here is a function of things that already existed: notebook 37's 25-fold paired
tests and the per-fold score tables of notebooks 34/35/37. None involves the new arm.

THE COMPARISON. One arm, `interval_2d6`, against notebook 37's existing `mae_only`, paired on the
25 cells of the frozen cluster-disjoint partition. The two differ in exactly one thing: CYP2D6's
label column carries the assay bound on compounds that were screened, failed promotion, and show
no inhibition. Nothing else in the table moves, and the SCORED population is unchanged -- interval
rows are never scored, on any fold, because the ground truth used for scoring is the unmodified
curated file where those cells are NaN.

THE BARS. Notebook 37's own pooled-SE procedure: 2 x the mean paired standard error across its ten
arms against `plain`, per isoform -- never 2 x the SE of the comparison under test, which would
just restate the significance threshold. ST-RAE comes from `paired_tests_vs_plain_25fold.csv`'s own
`se` column; that file has no Spearman rows, so Spearman's bars are recomputed from the per-fold
`Spearman_R` values the same way (notebook 45's method).

THE ASYMMETRY, recorded before the arm runs. The figure this arm is testing (CYP2D6 rho
0.378 -> 0.476) was measured by talindrew on their V2 surface, which contains the interval
compounds. V1 does not: it asks the harsher question of whether teaching the model where the low
end sits improves rank order among compounds that ALL have curves. The labelled CYP2D6 range runs
1.95 to 7.53 and the low end is where compression bites, so a gain is plausible -- but the
inference is one-sided. A positive result on V1 is strong evidence; a null is ambiguous, because
V1 cannot see the population the intervention is about.
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "outputs" / "54_interval_2d6_arm"
PRED_DIR = OUT / "predictions"
NB34, NB35, NB37 = (REPO / "outputs" / d for d in
                    ["34_butina_5x5", "35_butina_repeat3", "37_butina_5x5_complete"])
ISO = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]

# ---- refuse to run if the arm has already produced anything.
existing = sorted(PRED_DIR.glob("interval_2d6__*.csv")) if PRED_DIR.exists() else []
if existing:
    sys.exit(f"REFUSING TO FREEZE: {len(existing)} interval_2d6 prediction file(s) already exist "
             f"(e.g. {existing[0].name}). The pre-registration must precede the arm.")
OUT.mkdir(parents=True, exist_ok=True)

# ---- ST-RAE bars, from notebook 37's own paired-test table.
pt = pd.read_csv(NB37 / "paired_tests_vs_plain_25fold.csv")
assert sorted(pt.metric.unique()) == ["MAE", "ST-RAE"], "paired_tests metrics changed"
assert pt.arm.nunique() == 10, f"expected 10 arms vs plain, found {pt.arm.nunique()}"
assert (pt.n == 25).all(), "a paired test is not on 25 folds"
strae_bar = (pt[pt.metric == "ST-RAE"].groupby("isoform").se.mean() * 2).round(6).to_dict()

# ---- Spearman bars, recomputed from per-fold scores the same way (nb45's method).
rs = pd.concat([pd.read_csv(d / "run_scores.csv") for d in [NB34, NB35, NB37]], ignore_index=True)
rs = rs[rs.Endpoint != "MA"].copy()
rs["isoform"] = rs.Endpoint.str.split("_").str[0]
rows = []
for iso, g in rs.groupby("isoform"):
    wide = g.pivot_table(index=["repeat", "fold"], columns="arm", values="Spearman_R")
    for arm in [a for a in wide.columns if a != "plain"]:
        d = (wide[arm] - wide["plain"]).dropna()
        assert len(d) == 25, f"{iso}/{arm}: {len(d)} paired folds, expected 25"
        rows.append({"isoform": iso, "arm": arm, "se": d.std(ddof=1) / np.sqrt(len(d))})
sp = pd.DataFrame(rows)
assert sp.arm.nunique() == 10, "Spearman bars not built from the same ten arms"
spearman_bar = (sp.groupby("isoform").se.mean() * 2).round(6).to_dict()

# ---- assert the bars against the values the task states, to 3dp / 4dp as quoted.
EXPECT_ST = {"CYP1A2": 0.023, "CYP2C9": 0.021, "CYP2D6": 0.028, "CYP3A4": 0.016}
EXPECT_SP = {"CYP1A2": 0.0183, "CYP2C9": 0.0160, "CYP2D6": 0.0212, "CYP3A4": 0.0100}
for i in ISO:
    assert abs(round(strae_bar[i], 3) - EXPECT_ST[i]) < 1e-9, \
        f"ST-RAE bar {i}: recomputed {strae_bar[i]:.6f}, task states {EXPECT_ST[i]}"
    assert abs(round(spearman_bar[i], 4) - EXPECT_SP[i]) < 1e-9, \
        f"Spearman bar {i}: recomputed {spearman_bar[i]:.6f}, task states {EXPECT_SP[i]}"

PREREG = {
 "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
 "refuses_to_run_if": "any outputs/54_interval_2d6_arm/predictions/interval_2d6__*.csv exists",
 "interval_2d6_predictions_present_at_freeze": 0,
 "comparison": {
   "arm": "interval_2d6",
   "comparator": "mae_only, notebook 37's existing 25-fold arm -- nothing retrained for it",
   "paired_on": "the 25 (repeat, fold) cells of data/folds/cv_folds_butina_5x5.csv",
   "single_variable": ("CYP2D6's label column only: screened + unpromoted + no-inhibition readout "
                       "compounds carry the assay bound instead of NaN. All four scored columns on "
                       "curve-bearing rows, and the other three isoforms' columns entirely, are "
                       "bit-identical to mae_only's."),
   "scored_population": ("UNCHANGED from mae_only -- curve-bearing held-out compounds only. The "
                         "ground truth passed to the scorer is the unmodified curated file, in "
                         "which every interval compound is NaN on CYP2D6 and so is dropped by the "
                         "scorer for that endpoint. No interval window reaches the scoring side."),
 },
 "primary_metric": "Spearman",
 "co_primary_metric": "ST-RAE",
 "materiality_bars": {"ST-RAE": strae_bar, "Spearman": spearman_bar},
 "bar_procedure": ("notebook 37's pooled-SE rule: 2 x the mean paired SE across its ten arms "
                   "against plain, per isoform. NOT 2 x the SE of this comparison."),
 "verdict_rule": ("per isoform and metric: SIGNIFICANT if p_BH < 0.05 under Ash's parametric "
                  "route; MATERIAL additionally requires |mean paired difference| >= the bar. "
                  "A result that is significant but under bar is reported as DETECTED NOT "
                  "MATERIAL, which is this project's standing language."),
 "directional_prediction": ("CYP2D6 should gain. It carries both the strongest promotion filter "
                            "(AUC of screen activity for promotion = 0.9998) and the strongest "
                            "over-prediction on the population that filter removed (66.4% of "
                            "confirmed non-inhibitors predicted above a bound they are known to "
                            "be below) -- notebook 53 Parts 1 and 3.4."),
 "inference_asymmetry_recorded_before_the_run": (
   "The target figure (CYP2D6 rho 0.378 -> 0.476) was measured on talindrew's V2 surface, which "
   "CONTAINS the interval compounds. V1 does not. V1 asks whether teaching the model where the "
   "low end sits improves rank order among compounds that all have curves -- a harsher question. "
   "The labelled CYP2D6 range is 1.95 to 7.53 and the low end is where compression bites, so a "
   "gain is plausible. But the inference is ONE-SIDED: a positive result on V1 is strong evidence "
   "for the training claim; a NULL is AMBIGUOUS, because V1 cannot see the population the "
   "intervention is about. This is written here, before the arm runs, rather than concluded "
   "afterwards."),
 "built_in_control": (
   "CYP1A2, CYP2C9 and CYP3A4 labels are untouched, so ANY movement on those three measures "
   "Chemprop's multitask loss reweighting by valid-label count (notebooks 06 and 20), not "
   "chemistry. CYP2D6's labelled count goes 1,493 -> 4,141 (x2.77), so movement is expected. "
   "These three cells are reported under that heading and are NOT evidence about selection bias "
   "in either direction."),
 "secondary_outcome_prediction_spread": (
   "NOT a gate, measured and reported. 2,648 of the arm's 4,141 CYP2D6 training labels (64%) will "
   "be the identical constant 4.3054. Under MAE loss the conditional median is what is optimised, "
   "so the arm is at real risk of compressing its CYP2D6 prediction spread toward that point. "
   "This project's own standing association (notebook 27 Section 1: Spearman -0.822 between a "
   "submitted column's spread ratio and its board ST-RAE, observational) makes compression a "
   "known risk, so the arm's out-of-fold CYP2D6 prediction SD is compared against mae_only's and "
   "reported whichever way it falls."),
 "stated_limitation_of_the_target": (
   "A point target at the bound UNDER-corrects: a compound genuinely at pIC50 2.0 is pulled only "
   "to 4.3054, not to 2.0. The error is bounded, one-directional and conservative. Accepted in "
   "exchange for needing no donor model and no fitted window -- see decisions_taken.json."),
 "analysis": ("Ash's parametric route, Levene measured per cell, RM-ANOVA per metric at "
              "alpha = 0.05/2, post-hoc where the omnibus clears, CLD. With two arms the omnibus "
              "reduces to the paired test, as notebook 29 recorded for the same situation. "
              "BH across the four isoforms within each metric: two families of four. 10-fold and "
              "15-fold reported alongside 25-fold."),
 "provenance": {
   "fold_file_sha256": hashlib.sha256(
       (REPO / "data" / "folds" / "cv_folds_butina_5x5.csv").read_bytes()).hexdigest(),
   "nb37_paired_tests_sha256": hashlib.sha256(
       (NB37 / "paired_tests_vs_plain_25fold.csv").read_bytes()).hexdigest(),
 },
}
(OUT / "prereg.json").write_text(json.dumps(PREREG, indent=1) + "\n")
print("FROZEN -> outputs/54_interval_2d6_arm/prereg.json")
print(f"  interval_2d6 predictions present at freeze: 0 (asserted)")
print(f"  ST-RAE   bars: " + "  ".join(f"{i}={strae_bar[i]:.4f}" for i in ISO))
print(f"  Spearman bars: " + "  ".join(f"{i}={spearman_bar[i]:.4f}" for i in ISO))
print("  both bar sets reproduce the values stated in the task, asserted.")
