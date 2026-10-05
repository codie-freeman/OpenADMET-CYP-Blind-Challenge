# Notebooks

All notebooks run in numbered order and load only the frozen artifacts earlier notebooks wrote —
no notebook recomputes another's output. This file gives one short paragraph per notebook, grouped
by purpose, each tagged **COMPLETE**, **SUPERSEDED** (by what), or **ONGOING**. For the full
reasoning behind any result summarized here, read `CLAUDE.md`'s "Notebook status log" — that
section is the authoritative, detailed record; this file is the map. For any leaderboard figure,
`docs/leaderboard_submissions.md` is canonical — figures below are described, not restated, except
where a table is the clearest way to show a notebook's own internal (non-leaderboard) comparison.

Every score quoted anywhere in this repository — CV or leaderboard — is a half-set measurement
where it concerns the blind leaderboard: the live board scores only half the blind test set, split
by chemical series, with the full set scored only at the intermediate reveal (24 Sep 2026) and at
the competition's close (3 Nov 2026).

## Data (00–02)

**[`00_schema_audit.ipynb`](00_schema_audit.ipynb)** — COMPLETE. Structural audit of all five raw
HuggingFace files. Key finding: `TRAIN_TDI.csv`'s `direct_inhibition` columns are byte-identical to
`TRAIN_inhibition.csv` (a copy, not an independent measurement), so `TRAIN_TDI.csv`,
`TRAIN_Emax.csv`, and `single-concentration-TRAIN.csv` are set aside from the Direct Inhibition
track, leaving `TRAIN_inhibition.csv` (4,905 rows) and `TEST-BLINDED.csv` (750 rows) in scope. Also
traces the source of the `_conf_high`/`_conf_low`/`_std` columns via OpenADMET's own scoring code
and a Discord confirmation.

**[`01_data_curation.ipynb`](01_data_curation.ipynb)** — COMPLETE. Canonicalizes SMILES, generates
InChIKeys (`src/features.py`), and checks salts, parse failures, duplicates, and train/test
leakage (all by InChIKey). All checks clean: 0 salts, 0 parse failures, 0 duplicates, 0 leakage.
Writes `train_inhibition_curated.csv` and `test_blinded_curated.csv`.

**[`02_chemical_space_exploration.ipynb`](02_chemical_space_exploration.ipynb)** — COMPLETE.
Feasibility-checks a cluster-based CV split. Four findings: low pairwise redundancy; the blind set
is chemically *closer* to train than train is to itself; SALI activity-cliff prevalence is
literature-low but reference-range-plausible at looser thresholds; 96.4% of Bemis–Murcko scaffolds
are singletons. Conclusion: a scaffold-exact cluster split would functionally resemble random
splitting for most of the data. Also notes CYP2D6 has the lowest per-isoform NN similarity and did
not seed the real test set's construction the way the other three isoforms did — relevant later to
CYP2D6's leaderboard result.

## Framework (03–05)

**[`03_features_and_fold_split.ipynb`](03_features_and_fold_split.ipynb)** — COMPLETE, frozen.
Freezes the two artifacts everything downstream depends on: a 5×5 repeated random CV fold
assignment (`data/folds/cv_folds.csv`, seed 42) and three tabular feature sets (ECFP4 +
physicochemical descriptors, full RDKit 2D, CheMeleon embeddings). Never touched again after this
point.

**[`03b_log2fc_pretrained_encoder.ipynb`](03b_log2fc_pretrained_encoder.ipynb)** — ONGOING, paused.
Reports on a random-init Chemprop encoder pretrained on the 4,376-compound log2fc primary screen
(trained by `scripts/train_log2fc_encoder.py`) and frozen as a feature extractor. Not yet ablated
into any downstream model; paused pending a leakage check
(`scripts/train_log2fc_encoder_blind.py` exists but has not been run to completion).

**[`04a_baseline_screen.ipynb`](04a_baseline_screen.ipynb)** — SUPERSEDED by `05` (for config
selection purposes). Single-fold screen across 21 configs (naive baselines, 9 tabular, 10
Chemprop). CheMeleon-init multitask Chemprop won clearly on every macro metric (0.700 ST-RAE vs.
0.774 next-best). Its winning config fed `04b`'s submission before `05`'s properly-powered
comparison existed.

**[`04b_final_submission.ipynb`](04b_final_submission.ipynb)** — COMPLETE, ACTIVE reference. First
live submission — a single, unensembled `chemprop_chemeleoninit`, full-data retrain (via
`scripts/train_final_submission_multitask.py`), config chosen by `04a`. Every later submission in
this project is compared back against this one.

**[`04c_submission_diagnostics.ipynb`](04c_submission_diagnostics.ipynb)** — COMPLETE. Diagnoses
`04b`'s CYP2D6 underperformance. Rules out systematic bias, CI width, NN similarity, and
activity-cliff proximity as the cause; the leading hypothesis (CYP2D6 not hit-expansion-seeded the
way the other three isoforms were) was independently confirmed via the OpenADMET Discord.

**[`05_cv_comparison.ipynb`](05_cv_comparison.ipynb)** — COMPLETE, ACTIVE, foundational. The full
5×5 repeated CV (25 folds) comparison across all 11 real configs plus naive baselines, using Ash et
al. (2025)'s statistical protocol. This is the authoritative CV table for the whole project — pull
figures from `outputs/05_cv_comparison/`, not from any display subset.

**[`05b_cluster_cv_comparison.ipynb`](05b_cluster_cv_comparison.ipynb)** — COMPLETE, inconclusive
by design. Tests a Tanimoto/Butina-clustered alternative 5×5 CV split (`chemprop_chemeleoninit`
only) against whether a cluster-aware split predicts the real CV-to-blind gap better than the
original random split. Result: mixed, small-magnitude, no clean win either way — CYP2D6 does not
stand alone here, it pairs with CYP3A4. Nothing downstream depends on this result.

**[`05c_cluster_cv_leakage_sensitivity.ipynb`](05c_cluster_cv_leakage_sensitivity.ipynb)** —
COMPLETE. Tests whether the non-`chemprop_chemeleoninit` configs used in later Caruana ensembles
swing more under cluster-CV than random-CV — a second candidate mechanism for the CYP1A2/CYP3A4
blind regression seen later in notebook 12. Real signal found: CYP3A4's six configs shift 2–4×
`chemprop_chemeleoninit`'s own reference shift, all in the leakage-consistent direction. CYP1A2
shows no meaningful signal. Descriptive only.

## CYP2D6 investigation (06–09)

**[`06_outlier_check.ipynb`](06_outlier_check.ipynb)** — COMPLETE, flagged ambiguous. CYP2D6
residual- and CI-width-based outlier exclusion. Residual exclusion helped three configs in CV, but
also found that excluding CYP2D6 labels shifts `chemprop_chemeleoninit`'s predictions on the
*other*, untouched isoforms too — an unresolved, leading hypothesis of Chemprop's multitask loss
reweighting by valid-label count (later given strong, though not CYP2D6-conclusive, cross-isoform
support by notebook `20`). Its specific outputs feed only notebook `08`'s ensembles, which are
themselves superseded — flagged rather than forced into "superseded" or "core."

**[`07_weighting_tuning.ipynb`](07_weighting_tuning.ipynb)** — COMPLETE, flagged ambiguous (same
reasoning as `06`). CYP2D6 sample weighting (negative for both tested configs) and
hyperparameter tuning (helped `chemeleon__rf`, adopted; not significant for `chemeleon__lightgbm`).

**[`08_ensemble_selection.ipynb`](08_ensemble_selection.ipynb)** — SUPERSEDED by `11b`'s Caruana
selection. Per-isoform ensemble selection via `05`'s compact-letter-display CV tiers and simple
averaging. Adopted ensembles for CYP1A2/CYP2D6/CYP3A4; its per-isoform CV results are still valid
data, reused as input to notebook `11`'s wider candidate pool.

**[`09_placement_recalibration.ipynb`](09_placement_recalibration.ipynb)** — COMPLETE, negative
result. R² = 2ρk − k² − b² placement/recalibration diagnostic (SuperCowPowers' decomposition). The
placement-error share of the R² gap sits under 10% for all four isoforms — no correction adopted;
all four isoforms use identity recalibration. A real finding: ranking quality, not placement,
dominates the CV-to-blind gap at this point in the project.

## Ensembling (10–12)

**[`10_final_retrain_predict.ipynb`](10_final_retrain_predict.ipynb)** — COMPLETE, flagged
ambiguous. Second live submission, using `08`'s simple-average ensembles, full retrain. Scored
substantially worse than `04b` on every isoform, including the two isoforms `06`'s exclusion never
touched — a real, deliberate submission whose recipe has since been replaced by `11b`/`12`'s
Caruana selection.

**[`10b_blind_regression_investigation.ipynb`](10b_blind_regression_investigation.ipynb)** —
COMPLETE, root cause not definitively isolated. Investigates `10`'s regression. The
spread-compression hypothesis weakens under the correct metric (ST-RAE); a Caruana et al.
(2004)-style selection-bias mechanism is identified as a plausible contributor for CYP1A2/CYP3A4
specifically (not CYP2D6, which was already broken pre-ensembling).

**[`10c_control_submission.ipynb`](10c_control_submission.ipynb)** — COMPLETE, ACTIVE reference.
Control submission: plain, unensembled `chemprop_chemeleoninit`, full retrain under
`cyp-admet-v2`, no CYP2D6 exclusion — isolates the environment migration from `08`'s
ensembling/exclusion as the cause of `10`'s regression. Scored close to `04b`, resolving the
question: the `cyp-admet-v2` migration is **not** the cause of `10`'s regression.

**[`11_caruana_prep.ipynb`](11_caruana_prep.ipynb)** — COMPLETE. Assembles one pooled
(repeat, fold, inchikey)-level out-of-fold table per isoform across all 25 CV folds and all 11 real
configs from `05` — the candidate library the Caruana selection in `11b` runs against. No
retraining.

**[`11b_caruana_selection.ipynb`](11b_caruana_selection.ipynb)** — COMPLETE, ACTIVE, current
ensemble-selection method. Caruana bagged ensemble selection (`src/ensemble/caruana.py`) per
isoform against `05`'s full pool, scored on real ST-RAE, then capped at each isoform's empirical
CV-minimizing member count after a marginal-value-trajectory check found the uncapped ensemble
measurably worse in CV for two isoforms. A later correlation-matrix diagnostic (residual, not raw
prediction, correlation) confirmed the CYP1A2/CYP3A4 candidate pools are more redundant than
diverse.

**[`12_caruana_retrain_predict.ipynb`](12_caruana_retrain_predict.ipynb)** — COMPLETE. Third live
submission (NB12), using `11b`'s capped Caruana-weighted ensembles, full retrain. CYP2C9 reached
its best value on record here (later reused unmodified by notebook `19`).

**[`12b_caruana_vs_control_comparison.ipynb`](12b_caruana_vs_control_comparison.ipynb)** —
COMPLETE, descriptive only. Compares NB12 vs. `10c` prediction spread. Capping genuinely widened
blind-set prediction spread for all four isoforms versus the uncapped version, but did not fully
close the gap to `10c`'s own spread — a real, partial mitigation, not a fix.

## Calibration (13–16, 19)

**[`13_aid1851_blind_population_calibration.ipynb`](13_aid1851_blind_population_calibration.ipynb)**
— COMPLETE. (Notebook number reused 2026-09-12; an earlier, unrelated "13" — a placement-correction
retest — was built, found not adopted, and deleted; see `CLAUDE.md` for that history.) Fetches and
curates PubChem AID 1851 (an external NCATS qHTS cytochrome panel) as a population-moment estimate
only, never as training data, per two real challenge entrants' reported gains from external
population recalibration. Applies floored-mean calibration and submits it uniformly to all four
isoforms as **NB13-calib**, the fifth live submission — worse than both `04b` and `10c` on every
isoform.

**[`14_censored_mle_population_calibration.ipynb`](14_censored_mle_population_calibration.ipynb)**
— COMPLETE, negative result. (Number also reused; an earlier, unrelated "14" — a CYP2D6-revert
candidate — was built, never submitted, and deleted.) Re-estimates AID 1851's population moments
via censored (Tobit-style) MLE instead of `13`'s floored-mean substitution — moves CYP2D6's
estimated population center much closer to an external precedent. The actual point of the
notebook, an out-of-fold acid test, is a clean negative result: neither calibrated version beats
raw OOF ST-RAE for any isoform.

**[`15_analog_holdout_comparison.ipynb`](15_analog_holdout_comparison.ipynb)** — COMPLETE, decisive
negative. Builds a single train/test split approximating how the real blind set was actually
constructed (top-25 hits + Tanimoto-nearest analogs per isoform, CYP2D6 excluded from
hit-expansion per its confirmed construction). Result: this split's CV-to-blind gap is 2.3–9.3×
*larger* than the random-CV or cluster-CV gap — the opposite of the hoped-for outcome. Limited
going forward to error-driven investigation, not model selection.

**[`16_board_solved_population_calibration.ipynb`](16_board_solved_population_calibration.ipynb)**
— COMPLETE, SUPERSEDED for CYP2C9 by `19` (CYP1A2/CYP2D6/CYP3A4 remain the current recipe). Solves
each isoform's blind-population mean/spread algebraically from this project's own four earlier
published board metrics (R²/MAE/Spearman), then recenters `10c`'s raw predictions for
CYP2C9/CYP2D6. An addendum overrides CYP2D6 to a lower-root solution matching an external
precedent. Submitted as **NB16**, the sixth live submission — the best result at the time, a real
win for CYP2D6, a regression for CYP2C9. See `docs/leaderboard_submissions.md` for the figures.

**[`19_cyp2c9_revert_submission_candidate.ipynb`](19_cyp2c9_revert_submission_candidate.ipynb)** —
COMPLETE, ACTIVE, current best. Pure recombination, no retraining: CYP1A2/CYP2D6/CYP3A4 taken
byte-identical from `16`'s submitted file; CYP2C9 reverted to `12`'s own uncalibrated,
capped-Caruana-ensemble column, since `16`'s recentring had regressed it. Submitted as **NB19**,
the seventh live submission and the best result to date — see `docs/leaderboard_submissions.md`
for the full comparison.

## Screens (17–18, 20–23)

**[`17_protonation_charge_check.ipynb`](17_protonation_charge_check.ipynb)** — COMPLETE. Computes
physiological-pH net formal charge (Dimorphite-DL) and a logP proxy for every compound, then checks
per-isoform whether charge carries signal against pIC50 — a gate for whether a descriptor-augmented
Chemprop run is worth building. Mechanism confirmed for CYP1A2/CYP2C9/CYP2D6 (CYP3A4 fails the
bar); a confound check against lipophilicity shows all three surviving signals are largely real,
not logP artifacts, with CYP2D6 the cleanest case.

**[`18_charge_augmented_chemprop_screen.ipynb`](18_charge_augmented_chemprop_screen.ipynb)** —
COMPLETE, negative result. Single-fold, three-seed screen adding `17`'s net-charge column as a
Chemprop extra descriptor. Solo ST-RAE ties on all four isoforms; the charge arm is no more
decorrelated from the existing candidate pool than the pool already is from itself. No isoform
advances to a full 25-fold comparison on this evidence.

**[`20_residual_exclusion_screen.ipynb`](20_residual_exclusion_screen.ipynb)** — COMPLETE, mixed.
Single-fold, three-seed isolated test of `06`'s CYP2D6 residual-exclusion finding, never
previously tested on its own, with a random-masking control arm to separate a data-quality effect
from Chemprop's suspected multitask loss-reweighting-by-valid-label-count behavior. CYP2D6's own
score is NULL (unresolved either way, likely underpowered), but the cross-isoform result gives the
strongest evidence yet for the reweighting mechanism: two disjoint masked sets produce nearly
indistinguishable divergence on the three untouched isoforms.

**[`21_strae_offset_curve.ipynb`](21_strae_offset_curve.ipynb)** — COMPLETE, report only. Derives
the ST-RAE-optimal uniform additive offset per isoform directly from pooled OOF predictions and
training credible-interval widths. Finds a real but small (0.06–0.12 pIC50 unit) gap between the
ST-RAE-optimal and R²-optimal placement on every isoform — confirms the scoring-asymmetry mechanism
but concludes it is a much smaller lever than blind-population-centre placement (`16`'s mechanism).
No offset applied to any submission.

**[`22_auxiliary_heads_screen.ipynb`](22_auxiliary_heads_screen.ipynb)** — COMPLETE, mixed. Tests
whether auxiliary Chemprop prediction heads from this project's own already-released TDI/Emax files
add real information (unlike every prior screen, which only rearranged existing information).
CYP3A4 (the only isoform with a non-trivial pool-expansion ratio) improves as predicted; CYP2C9
improves unpredicted; CYP1A2 regresses under the combined TDI+Emax arm; CYP2D6 ties. CYP2C9 and
CYP3A4 are candidates for a full 25-fold CV comparison — not yet run.

**[`23_compound_pool_overlap.ipynb`](23_compound_pool_overlap.ipynb)** — COMPLETE, report only.
Read-only counting exercise: quantifies how many new compounds each of this project's own unused
challenge files would add per isoform. The large single-concentration file contributes zero new
compounds to any isoform (a strict subset of the curated pool); a small correction to notebook
`22`'s Emax count is found and reported (1 compound, immaterial). Compound-pool expansion from this
project's own released files is, on this evidence, largely exhausted outside `22`'s existing
CYP3A4/TDI finding.

## Ensemble and calibration diagnostics (24–27, 27b)

**[`24_nonnegative_stacking.ipynb`](24_nonnegative_stacking.ipynb)** — COMPLETE, report only.
Refits the ensemble stacker on the same 11-config pooled OOF predictions (`11`'s pool) with
non-negative combiners (NNLS, non-negative ridge), using a genuine nested outer-fold design —
unlike `11b`'s in-sample Caruana fit. Finds no coefficient-cancellation problem in this project's
own pool (unlike a cited real entrant's diagnosed failure mode on their own pool); every fitted
combiner beats the best single config, and non-negative combiners are never worse than unconstrained
ones. No adoption recommendation made.

**[`25_hard_compound_analysis.ipynb`](25_hard_compound_analysis.ipynb)** — COMPLETE, diagnostic
only. Characterizes the compounds every one of the 11 pooled configs gets wrong, to test whether a
planned retrieved-neighbour-corpus pretraining step targets the right population. Finding: consensus
-hard compounds are structurally *ordinary* relative to training — nearest-neighbour similarity does
NOT separate hard from easy on any isoform — so a retrieved corpus for "unseen chemistry" is aimed
at the wrong problem as currently conceived; the real, isoform-dependent hard population is instead
potent, narrow-interval, and lipophilic/rigid compounds.

**[`26_spread_sweep.ipynb`](26_spread_sweep.ipynb)** — COMPLETE, report only. Two analyses:
confirms consensus-hard compounds (from `25`) are systematically under-predicted (a clean
regression-to-the-mean signature), and sweeps a pure spread multiplier against OOF predictions,
finding the ST-RAE-optimal multiplier sits between the theoretical rho-shrinkage value and 1.0 on
every isoform — rho-shrinkage over-shrinks relative to what ST-RAE itself rewards. A separate
CYP2D6-only sweep against the deployed, blind-targeted column disagrees with the raw-OOF sweep and
is reported side by side, not reconciled. No correction applied to any submission.

**[`27_calibration.ipynb`](27_calibration.ipynb)** — ONGOING, the project's standing calibration
notebook: future placement/spread/board-metric work appends new sections here rather than spawning
new numbered notebooks, **unless a real submission already lives in this notebook** — see `27b`'s
entry below for why that carve-out now exists. **Section 1**: tests an external entrant's claim that
board ST-RAE is rank-ordered by submitted-column spread compression; finds a robust pooled
association across this project's own submission history (Spearman −0.82), though confounded with
ensembling/calibration differences across those submissions. **Section 2** (2026-09-16): widens
CYP2D6's spread to 0.85× of its training-label SD on top of `19`'s submitted file — sent as
**NB27-widened**, the eighth live submission (macro ST-RAE 0.6155, since superseded by `27b`'s
0.5982 as the best result to date); a
post-hoc note appended after the fact corrects the section's own "no action taken" closing text once
the real board result came back. A third section (solving NB30's CYP2D6 population and building a
corrected candidate) was briefly appended here on 2026-09-18 and then moved out the same day into
its own notebook, `27b`, after appending to and re-executing this notebook caused real problems —
see `27b`'s entry immediately below.

**[`27b_aid_cyp2d6_corrected.ipynb`](27b_aid_cyp2d6_corrected.ipynb)** — COMPLETE, standalone. **Not
a further section of `27`**, deliberately: `27` already carries a real, already-scored submission
(`NB27-widened`, Section 2), and appending further work there and re-executing top-to-bottom caused
real problems, so this notebook is fully self-contained instead — it reads nothing from `27`, only
files already on disk. Solves the blind CYP2D6 population from `30`'s own published board metrics
(reusing `16`'s solve verbatim, not reimplementing it), cross-checks the result against `16`'s
independent solve from four earlier submissions (close agreement, 0.18%/0.05% apart — real
convergence evidence), then applies both a placement correction (using the AID model's own OOF rho,
computed from `29`) and the same 0.85× widening `27` Section 2 used to `30`'s CYP2D6 column.
Assembles a mixed-recipe candidate (`10c`'s CYP1A2/CYP3A4 + `19`'s CYP2C9 + the corrected CYP2D6
column), validated PASS, and **sent 2026-09-19 12:11 UTC** — the best submission to date (macro
ST-RAE 0.5982). See `docs/leaderboard_submissions.md` for the full NB27-widened, NB30, and 27b
figures.

## AID 1851 as training signal (28–30)

**[`28_external_data.ipynb`](28_external_data.ipynb)** — COMPLETE, single-fold screen. Gate-then-
screen test of PubChem AID 1851 (already fetched by `13`) used as **training signal via auxiliary
Chemprop heads** — a different role from `13`/`14`'s calibration-target use. A pre-registered gate
found the potency readout's apparent signal was mostly a censoring-floor artifact; the efficacy
readout (`max_inhibition`) was chosen instead, user-confirmed after the notebook's own frozen rule
required stopping to ask. Result on this one fold: CYP1A2 and CYP2C9 improve; CYP2D6 and CYP3A4 tie.
CYP1A2 (much the stronger) and CYP2C9 flagged as candidates for a full 25-fold CV confirmation.

**[`29_cv_confirmation.ipynb`](29_cv_confirmation.ipynb)** — COMPLETE. Full, honest 5×5 repeated-CV
confirmation of `28`'s screen (50 real Chemprop runs, paired significance tests, Benjamini-Hochberg
corrected), run fully autonomously. `28`'s single-fold verdict reproduces for only 2 of 4 isoforms:
CYP1A2 and CYP2D6 are significantly better (p_BH = 0.0017, 0.00014); CYP2C9 and CYP3A4 are not
different from baseline — CYP2C9's earlier single-fold gain is identified as a selection artifact of
that one fold. Ensemble-diversity check: still no, error correlation sits at each isoform's own
fold-to-fold noise floor, not below it.

**[`30_aid_full_retrain.ipynb`](30_aid_full_retrain.ipynb)** — COMPLETE. Full-data retrain of `29`'s
confirmed recipe; submission candidate validated and, per the user's own review, submitted as
**NB30**, the ninth live submission. Flagged before submission that this model's raw predictions are
more compressed than `10c`'s on all four isoforms. **Scored worse than `10c` on all four isoforms**
(macro ST-RAE 0.8160) — the third instance on this project of a CV result failing to transfer to the
real blind board, and the most rigorously confirmed of the three. See
`docs/leaderboard_submissions.md`'s NB30 section for the full comparison and the resulting
mixed-recipe correction built in `27b`.

## Dead-zone training target (31–32)

**[`31_deadzone.ipynb`](31_deadzone.ipynb)** — COMPLETE, single-fold screen. Tests fitting to
`clip(oof_prediction, conf_low, conf_high)` under absolute-error loss instead of the point estimate
under squared error — the training target sits inside the credible band whenever the model's own
out-of-fold prediction already does, removing the incentive to chase noise within it. Built the
clipped target from `chemprop_chemeleoninit`'s pooled out-of-fold predictions only (never a model's
own in-fold prediction — two guards confirm this, including a hit-rate ceiling check). **The first
unanimous, clean win in this project's screening history**: the DEADZONE arm improves ST-RAE on all
four isoforms (4.5–6× each isoform's own baseline seed spread) and Spearman rises on all four too —
the finding that matters most, since an affine placement/spread correction provably cannot move
rank order. A post-hoc MAE_ONLY diagnostic arm decomposes the gain into a loss-function-alone
component and a clipping-on-top component: clipping does real, independent work on three of four
isoforms; CYP2D6's gain is mostly loss-function-alone, carrying `21`'s standing OOF-unreliability
caveat for that isoform. A stated fold-crosstalk caveat means the screen's *magnitude* is plausibly
optimistic (every training target came from an OOF donor model that had itself seen this screen's
held-out fold) but does not apply to a future blind submission and does not affect the *direction*
of the result. Combined with the AID architecture in a DEADZONE+AID diagnostic arm, the gain
extends further still — never itself proposed as an adoption candidate at this stage.

**[`32_deadzone_aid_retrain.ipynb`](32_deadzone_aid_retrain.ipynb)** — COMPLETE, two candidates
prepared and validated. Full-data retrain combining `29`'s CV-confirmed AID auxiliary-head
architecture with `31`'s screened dead-zone target — explicitly **not** itself confirmed at 25-fold
CV for this specific combination, a deliberate, reasoned deviation (cost, `31`'s own
fold-crosstalk caveat, and this project's own two-for-two prior record of CV-confirmed results
reversing on the board). CYP2D6 corrected using this model's own OOF rho, placed onto `16`'s
board-validated blind-population target, then widened to 0.85× the CYP2D6 training-label SD — the
same recipe `27` Section 2 and `27b` used. Two candidates built: **Candidate A** (all four columns
from this model, CYP2D6 corrected) and **Candidate B** (`10c`'s CYP1A2/CYP3A4 + NB19's CYP2C9 +
this model's corrected CYP2D6, the conservative single-variable change from `27b`). Both validated
PASS; both `gradio_client` cells fully commented out. **Candidate A was submitted as NB32-A, the
eleventh live submission** — worse than the appropriate per-isoform reference on all four isoforms
(macro ST-RAE 0.6488 vs. 27b's 0.5982), a fourth instance of a CV-favorable result failing to
transfer to the board. Candidate B was never submitted; its board value is nonetheless computable
because it shares Candidate A's CYP2D6 column byte-for-byte — worked out at macro ST-RAE 0.6066,
still worse than 27b, confirming the submission slot was correctly not spent on it. See
`docs/leaderboard_submissions.md`'s NB32-A section for the full comparison, the Spearman
decomposition, and the Candidate B deduction.

## Cluster-disjoint cross-validation (33–35, 37)

**[`33_butina_cluster_split.ipynb`](33_butina_cluster_split.ipynb)** — COMPLETE, comparison only.
Builds a cluster-disjoint Butina partition (2048-bit Morgan r=3, Tanimoto distance 0.65 — Pat
Walters' Discord-confirmed CYP3A4 recipe) and compares it against the frozen random partition,
motivated by entrants placed above this project claiming random splits flatter every model. **The
honest prior did not hold**: this recipe gives only 56.0% singleton clusters on this data, not the
93–96% that `02`'s scaffold analysis and an external entrant's figure implied. Two findings cut
against the motivating hypothesis: Butina held-out sets are *less* similar to their training
portion (median 0.365) than random ones are (0.438), both further from the blind set's own profile
(0.587) — so "Butina looks more like the blind set" is ruled out as the mechanism; and a
5,000-shuffle Monte Carlo pairing check finds the deliberate activity-contrast signature on
CYP1A2/CYP2C9/CYP3A4 (p≈0 each) but **a clean null on CYP2D6 (p=0.995)** — independent, internal
confirmation that CYP2D6 was not built by the same hit-expansion recipe. The consequential result
is that the AID-vs-PLAIN verdict flips sign on three of four isoforms between partitions. Three
implementation bugs were hit and documented rather than silently patched.

**[`34_butina_5x5_repeated_cv.ipynb`](34_butina_5x5_repeated_cv.ipynb)** — COMPLETE, evidence
only. Raises `33` to a repeated design: 2 repeats × 5 folds × 10 arms = 100 units, 10.71 h. **The
sanity check failed on first execution and was right to** — a real environment bug, not noise:
copying `33`'s `OMP_NUM_THREADS=1` changed BLAS/OMP floating-point reduction order and put the
retrained baseline 0.0378 ST-RAE off `05`'s stored scores, with bit-identical val_loss for two
epochs before diverging. Unset, it reproduces `05` exactly. The consequence for `33` is stated
rather than left implicit: its runs were produced under that setting, so its *absolute* Butina
numbers are not on the same footing as `05`/`29`/`31`'s. **`33`'s AID sign-flip does not survive a
repeated design** — 1 of 25 sign-resolved cells flips, 24 agree. The question TREES and SINGLE
existed to answer is answered cleanly: tree models on raw fingerprint/descriptor features sit
below their isoform's own re-seed floor on all 24 cells, this project's **first genuinely
decorrelated ensemble candidate**, while single-task models decorrelate not at all. A unit error in
`33`'s own prose is corrected (44% is the fraction of *clusters* that are non-singleton; the
compound share is 70.1%). The timing gate's 2-repeat reduction is reported as a limitation of the
rule — recomputed from realised costs it would have allowed 3.

**[`35_butina_repeat3.ipynb`](35_butina_repeat3.ipynb)** — COMPLETE, evidence only. Adds repeat 2
plus a new MAE_ONLY decomposition arm (65 units, 7.56 h), taking the design to 15 folds. Adding a
repeat is not re-deciding `34`'s gate, and the notebook makes that checkable: the fold assignment
was frozen before any result existed and the extra folds are added uniformly to every arm. **The
third repeat changed nothing directional** — 6 of 40 cells changed verdict, none a reversal, the
largest shift in any arm's mean 0.0126. One of `34`'s headline findings does weaken: CYP2D6/`aid`
is no longer resolved against BH. MAE_ONLY separates the dead-zone gain into a loss-function
component and a clipping-on-top component, and the split is isoform-specific (CYP2C9 82%
clipping; CYP2D6 71% loss-alone). Two real bugs found and recorded, one of which — a join that
silently emptied the blend pool — had passed a mid-run smoke test only because the arm it was
written for had not finished yet.

**[`37_butina_5x5_complete.ipynb`](37_butina_5x5_complete.ipynb)** — COMPLETE, evidence only.
Completes the design (repeats 3–4, 110 units, 12.15 h) to the **25 samples Ash et al. specify**,
with all three quoted sentences verified verbatim against the PDF at runtime. The caveat is stated
prominently: reaching 25 satisfies Ash's *count*, not their independence argument, since a repeat
here re-randomises cluster-to-fold assignment over one fixed clustering. **Completing the 5×5 moved
almost nothing** — 1 of 40 verdict cells changed, 0 sign reversals. The consequential result is
that **AID resolves on only one of four isoforms (CYP1A2), and CYP2D6 — the isoform `29` found
strongest — does not reproduce on cluster-disjoint folds**, the first matched 25-against-25
comparison this project has had. DEADZONE is significantly better than PLAIN on all four and
uniquely best under `05`'s own CLD protocol. The split question is now unanimous: **0 of 25
sign-resolved cells flip between partitions**. Its finding that CYP2D6's blend collapses to a
single member is **superseded by `39`** — see that entry.

## Spread correction and placement diagnostics (36, 38)

**[`36_spread_corrected_nb32a.ipynb`](36_spread_corrected_nb32a.ipynb)** — COMPLETE, submitted.
Applies a mean-preserving spread correction to `32`'s candidate A. **The brief's central premise
was wrong and is corrected in Part 0**: that candidate was not awaiting a check — it had already
been sent as NB32-A and scored 0.6488. The correction improved the work, turning the compression
question from a forecast into something answerable from NB32-A's own per-isoform board scores. Two
real discrepancies in `32`'s saved diagnostics are reported (a mixed `ddof` convention, and a
CYP2D6 row that describes the raw model rather than the candidate's column). Targets are each
isoform's best-scoring board column; all three judged isoforms clear the frozen 10% bar. **An
unprompted finding is the most consequential thing in it**: this candidate *cannot* beat `27b` even
if the widening worked perfectly, because its CYP2D6 column passes through unchanged and that gap
is a ranking decline no affine correction reaches. Five defects in its own verification code were
found in review and fixed. **Submitted as NB36-widened, macro ST-RAE 0.627** — it confirmed all
five pre-registered predictions, the first clean sweep on this project.

**[`38_cyp3a4_placement_diagnostic.ipynb`](38_cyp3a4_placement_diagnostic.ipynb)** — COMPLETE,
diagnosis only. Asks whether CYP3A4's residual gap after `36` is a placement offset. **Answer: no
— it is predominantly ranking, which no affine correction can reach.** The discriminating
comparison needs no model: after `36`'s widening this column's CYP3A4 SD is identical to `04b`'s to
2.2e-16, because `36`'s target rule selected `04b`'s own ratio, so **scale is eliminated by
construction** and what remains is Spearman 0.8020 against 0.8329. A joint R² decomposition across
all on-record submissions — using the bivariate-normal relation rather than `16`'s Spearman-as-
Pearson substitution, which removes an infeasibility that affected 5 of 11 rows — attributes
`04b`'s advantage to correlation (+0.0427) with placement working slightly in the *submitted*
column's favour (−0.0057). The OOF cross-check **disagrees in sign**, reproducing `21`'s
OOF-versus-blind population mismatch on a second isoform. When the missing board row was later
supplied, the population solve passed a genuine held-out test: it had predicted this column's R² at
0.6515 before the value existed in the repo, against a measured 0.6563 — inside the fit's own RMS
residual.

## Ensemble diversity and the single-object blend (39–40)

**[`39_ambiguity_and_shared_blend.ipynb`](39_ambiguity_and_shared_blend.ipynb)** — COMPLETE,
analysis only. Implements Krogh & Vedelsby's ambiguity decomposition (quotations verified against
the PDF; the identity `E = Ē − A` checked to 4.4e-16) to ask whether a blend's diversity survives
into the blind set. **It does not collapse — it grows**: `A_blind/A_oof` is 1.28–1.50 across
isoforms. Retrospectively, `08`'s ensembles — the ones NB10 lost 0.116 macro with — were *also* not
collapsed (1.11–1.79), so **vanished diversity is not the mechanism behind that failure** and it
still needs its explanation elsewhere. Part B finds a single shared member set costs +0.0003 macro,
and **one shared weight vector is 0.0027 macro *better* than the unconstrained per-isoform blend**,
driven almost entirely by CYP2D6 — the per-isoform search was overfitting its own folds, which
supersedes `37`'s reading of CYP2D6's collapsed blend. A later-appended Part D records the sharpest
limitation: the ambiguity check is **one-directional**, able to close a route but not to open one.

**[`40_single_object_ensemble.ipynb`](40_single_object_ensemble.ipynb)** — COMPLETE,
prepare-and-validate. Builds `39`'s shared-weight object on full data (`deadzone` 0.625 / `aid`
0.25 / `ecfp4_narrow__xgboost` 0.125, re-derived rather than quoted) and prepares submission
candidates. CYP3A4's widening is *predicted* rather than asserted, from the R² identity: it helps
iff ρ > 0.676, which every plausible value exceeds. Part 8 corrects its own earlier diagnosis —
the xgboost member is **not** seed-unstable; at these hyperparameters XGBoost is deterministic
given its data, and the entire instability is the 15% early-stopping split draw. The one piece of
genuinely adverse evidence is reported as such: the blind ambiguity ratio came back **mixed and
partly adverse** for this specific three-member object, below 1.0 on three of four isoforms under a
scale-free normalisation. **An unintended submission fired from this notebook on 2026-10-02** — the
first non-deliberate submission in this project's history, caused by the comment-gate convention
being unsafe against an editor. The notebook now uses an explicit `SEND = False` flag with the
printed message computed from it, so the cell cannot contradict itself.

## Quantum-chemistry descriptors (41–43)

**[`41_dataset_audit_3d.ipynb`](41_dataset_audit_3d.ipynb)** — COMPLETE, audit only. Replaces the
assumptions behind a 3D/QM descriptor plan with measurements over all 5,655 compounds. **The size
tail is not the problem it was expected to be**: only 3 compounds exceed 50 heavy atoms and none of
them is in the blind set, so cost is dominated by about ten named compounds rather than by an
isoform. Nine elements are present and nothing exotic, so GFN2's parameterisation covers
everything. **Chirality splits the verdict**: canonicalisation preserves stereochemistry exactly
(0 differing compounds across four sources), but only 9.4% of training and **2.4% of blind**
compounds have every centre specified, and no blind compound has defined double-bond
stereochemistry — so a descriptor claiming to separate enantiomers would have to invent a
configuration for most scored compounds. A basic-nitrogen SMARTS, validated 10/10 against
reference compounds, covers 18.9%/11.6% of train/blind but reaches **37.4% on CYP2D6 against
12.7–13.8% elsewhere**. Redundancy against existing feature sets is **nil**: `mordred_pca`'s block
is exactly Mordred's `ignore_3D=True` set and contains 0 of its 213 3D-exclusive descriptors.

**[`42_orca_timing_run.ipynb`](42_orca_timing_run.ipynb)** — COMPLETE, measurement only. The first
notebook here to run real quantum chemistry: 43 ORCA jobs over 15 deliberately chosen compounds,
to choose the production route on evidence. **Three of its own brief's environment premises were
wrong and were corrected by measurement** — `! XTB2` invokes a legacy external driver that is not
installed (the native keyword is `NATIVE-GFN2-XTB`); MPI-parallel DFT fails at every `nprocs ≥ 2`
under the briefed keyword line, which turns RIJCOSX on, while adding `NORI` both avoids the broken
path and runs 3.2× faster for ≤0.0019 eV of frontier-orbital difference; and **ORCA exits 0 on
error termination**, so the exit code is trusted nowhere afterwards. Verdict: full-set DFT is
affordable. Zero SCF convergence failures, including on all three chemotypes `41` flagged as
risks. A Discord source is transcribed into `docs/stereochemistry_provenance.md` with its status
recorded as `CONSISTENT_WITH_LOCAL_DATA` rather than `VERIFIED`, since a chat message has no
hashable artefact — and the six compounds it names check out 6/6 against this project's own data.

**[`43_coordinating_atom_candidates.ipynb`](43_coordinating_atom_candidates.ipynb)** — COMPLETE,
candidate definitions and evidence only; no pattern chosen and no feature file written. Part A puts
the ORCA parse on record (see `scripts/README.md` and `outputs/qm_parsed/README.md`). Parts B–E
measure three SMARTS candidates for the atom coordinating the haem iron. **A premise in its own
brief is corrected**: Correia & Hollenberg's "neutral ferric P450 species" qualifies the *enzyme's*
redox state, not the ligand, and the passage concerns anions — so the QM run's neutral-species
choice stands but is re-attributed to this project rather than to the source. The same sources
force a reframing: the chapter specifies **no amine substitution class** and puts aliphatic and
aromatic nitrogen on equal footing, so the *faithful* reading is the wider pattern, not the narrow
one. **The central result is a bind**: the conservative candidate (identical to `41`'s strict
pattern on all 5,655 compounds) reproduces its 2.84× CYP2D6 enrichment but is empty for 88.4% of
the scored blind set, while the faithful candidate reaches 79% coverage with no isoform spread at
all and whose presence flag stops separating labels on two isoforms. The tie-break is close to a
coin flip (the two proposed rules disagree on 43–48% of multi-match compounds, and the choice of
charge *scheme* is a third axis), and **the presence flag alone is significant on all four
isoforms**, reproducing `17`'s sign pattern from one bit — evidence bearing on whether the arm is
worth building at all. All three decisions are left open.
