# Standing findings

Results that are **independent of the notebook that produced them** and must survive it. Each
entry records the finding, the method, the provenance, and what it does and does not license.

A finding earns a place here when it is still true after the notebook that found it is gone —
because it is a property of the data, the board's returned metrics, or the project's own record,
rather than of a particular analysis. Notebook-specific results stay in CLAUDE.md's status log.

---

## 1. The blind CYP2D6 Pearson ρ is 0.5600, and both a-priori estimates undershoot it

**Finding.** The Pearson correlation between this project's `chemprop_chemeleoninit` CYP2D6
predictions and the blind truth is **ρ = 0.5600**, recovered with **no Spearman-for-Pearson
substitution anywhere**. The two estimates in use both undershoot it, in the same direction:

| quantity | value | error |
|---|---:|---:|
| **ρ recovered** (4 equations, no proxy) | **0.5600** | — |
| board Spearman, used directly as ρ by notebook 16 | 0.4000 | **−0.1600** |
| ρ = 2·sin(πρ_S/6), notebook 38's normal-theory conversion | 0.4158 | **−0.1442** |

Three-subset range 0.5502–0.6156. The normal-theory conversion is an improvement on using Spearman
raw, but closes only about a tenth of the gap.

**Method.** Submissions whose prediction columns are affine transforms of one another share a
single Pearson ρ against the blind truth. Each **distinct** member of such an ordering group
contributes one equation

> R² = 2ρk − k² − b²,  k = σ_p/σ_y,  b = (p̄ − ȳ)/σ_y

in the three unknowns (ȳ, σ_y, ρ). Three distinct members identify the system exactly. CYP2D6 has
four: `10c`, `NB13-calib`, `NB16`, `NB27-widened` — four genuinely different affine maps of one
ranking. All four exact three-subsets solve to rms ≈ 1e-16; the overdetermined four-equation fit
gives **ȳ = 3.0515, σ_y = 1.5397, ρ = 0.5600** at rms 6.4e-3.

**The sign ambiguity is broken.** The b² term is even in (p̄ − ȳ), which is why notebook 16, working
from one submission at a time, could not recover the sign of the population offset. With four
variants sitting at *different* prediction means the per-equation mirrors no longer coincide, and
the solver finds **exactly one root**.

**Consequence.** ρ sets the slope of every ρ-shrunk placement correction this project has deployed
(NB16 → NB19 → NB27-widened → 27b → NB32-A → NB47-mae → NB48-ensemble). An understated ρ shrinks
too hard. This is the same direction notebook 26 reached independently from an out-of-fold spread
sweep, where the ST-RAE-optimal multiplier sits strictly above ρ on every isoform — two routes,
same sign.

**What it does not license.** It does **not** invalidate the deployed corrections. The same solve
puts the blind CYP2D6 population at mean 3.0515 / SD 1.5397 against NB16's deployed target of
3.1558 / 1.5114 — agreeing to ~0.10 and ~0.03 pIC50 despite the ρ error, because ρ enters the
*slope* and the population enters the *target*, so the two errors do not compound. It is also one
isoform, and the one notebook 21 established cannot be ranked on out-of-fold ST-RAE.

**ρ is a property of a model's ordering, not of the isoform.** 0.5600 belongs to the
`chemprop_chemeleoninit` CYP2D6 ranking. The *population* (ȳ, σ_y) is model-independent and
transfers; ρ does not. Any other model's ρ must be implied from its own board R² against that
population, not copied.

**Provenance.** Derived 2026-10-22 from the fourteen fully-scored submissions' published R² plus
each submitted CSV's own mean and SD (`ddof=1`). Backing data preserved at
`outputs/55_superseded_preserved/{rho_recovered_vs_proxy,solve_three_subsets,solve_identified}.csv`.
Originating notebook superseded and deleted; the finding does not depend on it.

---

## 2. Identification status per isoform — count distinct equations, not group members

**Finding.** Counting how many *distinct* affine variants of one ordering have been submitted:

| isoform | distinct equations | status |
|---|---:|---|
| **CYP2D6** | **4** | **IDENTIFIED** — (ȳ, σ_y, ρ) solvable with no proxy |
| CYP1A2 | 2 | unidentified — 2 equations, 3 unknowns |
| CYP2C9 | 2 | unidentified |
| CYP3A4 | 2 | unidentified |

**A byte-identical resubmission contributes no new equation.** It carries the same mean, the same
SD and the same board R², so it restates an equation already in the system. Counted as *members*,
every isoform's largest group has 5 or 6 and all four would look identified. Counted as *distinct
equations*, only CYP2D6 is.

**Sharpened 2026-10-27: being a distinct equation is necessary but NOT sufficient. A
mean-PRESERVING variant adds an equation without breaking the sign ambiguity.** The b² term is even
in (p̄ − ȳ), so the mirror root ȳ′ = 2p̄ − ȳ is killed only by a member sitting at a *different
prediction mean*. Two members that share a mean share their mirror, and it survives into the
solved system. CYP1A2 is the worked case: its two distinct equations are NB32-A and NB36-widened,
and notebook 36's widening was **mean-preserving** (both columns sit at mean 5.0478), so solving
the pair returns **four roots** — ȳ = 4.2543 / 4.3007 / 5.7950 / 5.8413 — in two mirror pairs, and
notebook 16's own ambiguity ratio against the training mean 4.9554 is 1.26 and 1.28, *undetermined*
under its own 1.5 bar. σ_y is not pinned either (1.49 against 2.64 at the same ρ).

The structural cause is visible in one number per isoform — how far this project's own submitted
prediction mean has ever moved: **CYP2D6 1.4778 pIC50** (identified), CYP3A4 0.7233, CYP2C9 0.3731,
**CYP1A2 0.0892** (four-fold ambiguous). Identifiability is a property of the submission history,
not of the chemistry.

**Why CYP2D6 and not the others.** It is the isoform this project corrected most often *and with
the most different transforms* — raw, population calibration, placement recentring, spread
widening. CYP1A2's and CYP3A4's six-member group is `10c` resubmitted byte-identically five times
plus one calibration: two equations.

**Consequence.** Placement is computable from the board's own returned R² for CYP2D6 and not for
the other three. Any correction proposed for CYP1A2, CYP2C9 or CYP3A4 on a solved-population basis
rests on an assumed ρ, not a recovered one.

**The cheapest way to identify another isoform** is to submit a variant of an ordering already on
the board whose **prediction mean differs** — not merely a different affine map, and not a
resubmission. A mean shift is what breaks the mirror (see the sharpening above); a pure rescale
does not. With two mean-preserving equations already in hand, one mean-shifted third column makes
the system three equations in three unknowns and returns a unique root, with ρ *recovered* rather
than proxied — the same standard as CYP2D6's solve in finding 1.

**Provenance.** Grouping verified on `Spearman == 1.0` to floating point, then each within-group
member verified genuinely **affine** to its base (max |residual| ≈ 1e-14 across all groups), which
is what licenses one shared ρ per group. Backing data at
`outputs/55_superseded_preserved/{ordering_groups,affine_constants}.csv`.

---

## 3. Seventh transfer observation — the board and cross-validation disagree on model *ordering*

**This is the strongest finding on this project's record.** It belongs in the main results.

**Finding.** Ranking submitted orderings on board Spearman and on 25-fold cluster-disjoint
out-of-fold Spearman, for the four orderings whose architecture exists as an arm, the two sources
**disagree on 2 of 4 isoforms, and CYP3A4 is fully reversed**:

| Isoform | 25-fold CV says | The board says | |
|---|---|---|---|
| CYP1A2 | NB30 > NB47-mae > 10c | NB30 > NB47-mae > 10c | agree |
| CYP2C9 | NB30 > NB47-mae | NB47-mae > NB30 | **disagree** |
| CYP2D6 | NB30 > NB47-mae > 10c | NB30 > NB47-mae > 10c | agree |
| CYP3A4 | NB30 > NB47-mae > 10c | 10c > NB47-mae > NB30 | **disagree, fully reversed** |

**Why it is the strongest of the seven.** Transfer failures 1–6 on this project are all **ST-RAE**
failures, and ST-RAE mixes model quality with placement — so "the model was fine, the calibration
was not" was always an available explanation, and on several occasions it was the right one.
**Spearman is affine-invariant on both sides.** No placement or spread correction can move it, in
cross-validation or on the board. So this disagreement cannot be explained by calibration: the two
estimates disagree about the models themselves.

**It was predicted.** Notebook 40 wrote down, *in advance*, that the out-of-fold-to-board Spearman
mapping is not monotone on this project's record, citing `aid`'s out-of-fold CYP3A4 Spearman
exceeding `plain`'s while NB30's board Spearman (0.7651) sits below `10c`'s (0.8177). Notebook 48
then hit exactly that. This is the same effect measured systematically at 25 folds.

**Consequence.** Cross-validated ordering quality is not a reliable guide to board ordering quality
on this dataset, on at least half the isoforms, and the failure is not fixable by better
calibration. Model selection on out-of-fold Spearman carries the same risk this project has already
paid for six times on ST-RAE.

**Coverage limitation — the comparison does not cover the arms the cross-validation ranked
highest.** Eleven arms were compared at 25 cluster-disjoint folds; **three have ever been
submitted** — `plain`, `aid` and `mae_only`. The other eight (`deadzone`, `single`, and all six
tree arms) have no board reading at all, and **`deadzone` is among them: the arm that holds the
compact-letter-display `a` on all four isoforms and on the macro, uniquely and tied with nothing.**
So every CV-to-board comparison on this project is made over the three arms that happened to be
submitted, which are not the three the cross-validation ranked highest, and the disagreement above
is measured on the arms that were available rather than on the ones under test. The one submission
built on the dead-zone target, NB32-A, applies it to a **different architecture** (the AID
auxiliary-head model) and so does not fill the gap. Whether cross-validated ordering transfers for
the arm it most favours is untested and, with submissions closing 2026-11-03, likely to stay that
way.

**Arms and directories.** `plain` ← `10c`, `aid` ← NB30, `mae_only` ← NB47-mae. `plain`'s repeats 0
and 1 are in `outputs/34_butina_5x5/`, repeat 2 in `outputs/35_butina_repeat3/`, repeats 3 and 4 in
`outputs/37_butina_5x5_complete/`; `mae_only`'s repeats 0–2 are **all** notebook 35's. All arms
verified at exactly 25 distinct (repeat, fold) cells. `04b` is excluded from the comparison: it maps
to the same arm as `10c` but is a different training draw forming a different board ordering, and
one arm cannot rank two draws of itself.

**Provenance.** Backing data at
`outputs/55_superseded_preserved/{ordering_rank_agreement,ordering_rank_cv_crosscheck}.csv`.

---

## 4. Pooling across ordering groups costs 0.005 pIC50 — delavega's 0.5 is not a pooling artefact

**Finding.** Solving the blind population by pooling prediction vectors across ordering groups with
Spearman standing in for Pearson — the external entrant **delavega**'s route — was reproduced on
this project's own record and decomposed against the identified solve:

| route | ȳ (CYP2D6) |
|---|---:|
| (a) one ordering group, ρ **free** (identified) | 3.0515 |
| (b) one ordering group, ρ = proxy | 3.1231 |
| (c) all groups pooled, ρ = proxy (delavega's route) | 3.1284 |

| component | effect |
|---|---:|
| **proxy** (b − a) | +0.0717 pIC50 |
| **pooling** (c − b) | **+0.0053 pIC50** |
| total (c − a) | +0.0769 pIC50 |

**Pooling across ordering groups costs essentially nothing**, and the total disagreement between
the two routes is **about 6× smaller** than the ~0.5 pIC50 delavega reports. **So delavega's 0.5
pIC50 disagreement is not reproduced as a pooling artefact, and remains unexplained.**

**Where a 0.5-scale disagreement does live on this record.** It is the identified solve's own
**three-subset spread**, 0.479 pIC50 — i.e. *which three submissions you solve from* matters about
90× more than whether you pool orderings or proxy ρ. The dominant uncertainty is the conditioning
of the system, not either of the two things that distinguish the methods.

**Consequence, stated against the identified solve's own interest.** Grouping by ordering does
**not** earn its place for estimating the population — it costs equations (it is why three of four
isoforms are unidentified) and buys a shift of 0.005 pIC50 over simply pooling. It earns its place
for the one thing pooling cannot do: recovering ρ without a proxy (finding 1).

**Independent corroboration of the population machinery.** The pooled CYP3A4 solve gives mean
**4.7780**, SD **1.4150**, against notebook 38's own out-of-sample-validated solve of **4.7778** /
**1.3971** — agreeing to **0.0002** and 0.018 by a different route on different inputs.

**Status of delavega's figures.** `delavega` appears nowhere in this repository (checked directly).
Every figure attributed to that source is **external and not independently verifiable here**,
recorded in the same spirit as the `ST-RAE = 0.695 − 0.210·R²` field fit and
`docs/stereochemistry_provenance.md`'s `CONSISTENT_WITH_LOCAL_DATA` status. What is local is the
reproduction of the *method* on this project's own record.

**Provenance.** Backing data at
`outputs/55_superseded_preserved/{disagreement_decomposition,pooled_population_proxy}.csv`.

---

## 5. A one-variable controlled send pair: rank order does not move, and ST-RAE disagrees with R²

**This is the only controlled experiment this project has run on the live board.** Every other
pair of submissions differs in model, feature block, correction layers or several at once.
NB56b-s050 and NB56b-s100 differ in **exactly one scalar** — the shrinkage fraction `s` applied to
the affine placement/scale correction on CYP1A2 and CYP2C9 — with the same base file
(NB36-widened), the same solved blind populations, the same recovered ρ, the same code path, sent
twelve hours apart. Both are half-set measurements.

| | s = 0.5 (NB56b-s050) | s = 1.0 (NB56b-s100) | Δ |
|---|---:|---:|---:|
| macro ST-RAE | **0.5861** | 0.5933 | **+0.0072** |
| macro MAE | 0.6973 | **0.6864** | −0.0109 |
| macro R² | 0.5249 | **0.5366** | +0.0117 |
| macro Spearman | 0.6894 | 0.6894 | **0.0000** |
| macro Kendall | 0.5118 | 0.5118 | **0.0000** |

### 5a. The lever moved only the two isoforms it touched

CYP2D6 and CYP3A4 were carried byte-identical into both files. **All ten of their cells are exactly
equal across the two submissions** — not equal to four significant figures, equal at every digit
the board prints:

| Isoform | ST-RAE | MAE | R² | Spearman | Kendall |
|---|---:|---:|---:|---:|---:|
| CYP2D6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| CYP3A4 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |

Verified at source as well as on the board: max |difference| between the two submitted CSVs is
`0.00e+00` on both columns. This is the **seventh** confirmation on this record that the board
returns identical scores for identical input, and the first where it serves as a within-submission
control rather than a resend check.

The two isoforms that did move:

| Isoform | ST-RAE | MAE | R² | Spearman | Kendall |
|---|---:|---:|---:|---:|---:|
| CYP1A2 | **+0.0136** | −0.0413 | +0.0397 | 0.0000 | 0.0000 |
| CYP2C9 | **+0.0154** | −0.0022 | +0.0072 | 0.0000 | 0.0000 |

### 5b. Rank order did not move — measured, not argued from algebra

**Not one of the eight rank-metric cells differs.** Spearman and Kendall are identical on all four
isoforms and at macro level.

This is the project's central claim about its own method — an affine correction changes where
predictions sit and how far apart they are, and cannot change which compound ranks above which —
and it has until now rested on algebra plus resend checks. Here it is measured on predictions that
genuinely changed: CYP1A2 and CYP2C9 were transformed by two *different* affine maps and returned
the *same* ordering. On CYP2D6 and CYP3A4 the result is trivial, the files being identical; the
load-bearing evidence is the two corrected isoforms.

It is also a **three-way** agreement, not a two-way one. Both submissions reproduce NB36-widened's
own CYP1A2 Spearman (0.7565) and CYP2C9 Spearman (0.7554) exactly, so the base and both corrected
descendants agree across three separate sends.

**What this licenses.** Any claim that a placement or spread correction improved *ranking* is
false by construction, and the board now says so directly. It follows that the ordering quality of
every affine descendant in this project's record is a property of the model underneath it — which
is what makes finding 3's transfer disagreement a statement about models rather than calibration.

**What it does not license.** It says nothing about whether the correction improved the submission
overall, which is 5c, and nothing about rank order under a *non*-affine change — the two ensembles
(NB40-ensemble, NB48-ensemble) did move Spearman, as they must.

### 5c. ST-RAE and R² disagree about which submission is better

On both isoforms that moved, **ST-RAE got worse while R² got better**:

| | CYP1A2 | CYP2C9 | macro |
|---|---:|---:|---:|
| ST-RAE (lower better) | +0.0136 **worse** | +0.0154 **worse** | +0.0072 **worse** |
| R² (higher better) | +0.0397 **better** | +0.0072 **better** | +0.0117 **better** |
| MAE (lower better) | −0.0413 **better** | −0.0022 **better** | −0.0109 **better** |

**ST-RAE is the odd one out: MAE moves with R², not with ST-RAE.** That is a stronger statement
than a two-metric divergence, and it is the one worth carrying — fully correcting the placement
improved two of the three error metrics and worsened the third.

Stated as leaderboard records rather than deltas: **s = 1.0 takes the record on six cells** — macro
MAE and R², CYP1A2 MAE and R², CYP2C9 MAE and R² — and **every one of them from s = 0.5**, which
keeps the ST-RAE record. The same pair of files holds the best readings this project has achieved
on two metrics and on the third.

**Why it happens, mechanically.** ST-RAE scores a prediction as zero error anywhere inside the
compound's credible interval and as distance to the nearest edge outside it. Those intervals are
wide at low pIC50 and narrow at high (Spearman −0.56 to −0.93 against pIC50, notebook 21), so
over-prediction is cheap and under-prediction is punished. A correction that shifts predictions
**down** onto a low-centred blind population — which is what every correction here does, all three
isoforms being low-centred — moves mass from the cheap side to the expensive side. R² and MAE have
no such asymmetry and score the same move as an improvement. The optimum therefore sits at
different places for the two metrics, and notebook 21 measured that gap directly on this project's
own out-of-fold predictions: the ST-RAE optimum sits **+0.06 to +0.12 pIC50 above** the R² optimum
on all four isoforms. **This pair is the first time that gap has been observed on the real board
rather than out of fold**, and its direction matches.

**Provenance of the mechanism claim.** The asymmetry and the +0.06–0.12 figure are this project's
own (notebook 21, `outputs/21_strae_offset_curve/`). The framing was prompted by an external
entrant's claim; that source is **not present in this repository and is not independently
verifiable here**, and is recorded in the same spirit as the `ST-RAE = 0.695 − 0.210·R²` field fit
and `docs/stereochemistry_provenance.md`'s `CONSISTENT_WITH_LOCAL_DATA` status. What is local is
the measurement.

### 5d. s = 0.5 needs no cherry-picking — on ST-RAE, and only on ST-RAE

**s = 0.5 wins or ties on every isoform on ST-RAE**: it wins CYP1A2 (0.5631 vs 0.5767) and CYP2C9
(0.5408 vs 0.5562) and ties the two untouched ones. There is no isoform on which s = 1.0 is the
better ST-RAE choice, so selecting s = 0.5 requires no per-isoform mixing and keeps the submission
one model plus one scalar.

**That dominance is metric-conditional, and the condition is the whole reason it holds.** On R² the
comparison reverses on all four (two wins and two ties for s = 1.0); on MAE likewise. s = 0.5 is
the right final submission **conditional on ST-RAE being the metric that decides** — not because it
is the better-calibrated file. On the evidence, s = 1.0 is the better-calibrated file.

**The condition is assumed rather than confirmed.** ST-RAE is first in the vendored evaluator's
`ACTIVITY_METRICS` and the metric its docstring uses as the exemplar `primary_metric`, and this
project has treated it as governing throughout; but no default `primary_metric` is set in the
vendored code, the live backend's sort key is not in this repository, and `CLAUDE.md` carries the
standing caveat that the vendored evaluator has never been diffed against the live backend. The
recorded board ranks cannot settle it either — they are not monotone in macro ST-RAE across dates,
which an actively-resubmitting field explains and which makes ranks useless for cross-date
comparison. **So finding 5d is conditional, and the condition is the one thing here that could make
the final-submission choice wrong without any of 5a–5c being wrong.**

**What the pair does not settle.** `s` was sampled at two points, so nothing here locates the
optimum; the ST-RAE optimum lies at or below 0.5 and the R² optimum at or above 1.0, and neither is
bracketed. Both readings are half-set, and the half-to-full drift measured on CYP2D6 (−0.2789 or
+0.4989 pIC50, `outputs/56_cyp1a2_identification/half_to_full_drift.json`) is larger than the
0.0072 macro gap being used to choose between them — so the full-set reveal could reverse this
choice. And neither submission tests the base model: three of NB36-widened's four columns are
unchanged in s = 0.5.

**Provenance.** Board readings in `docs/leaderboard_submissions.md` (NB56b-s050 and NB56b-s100
sections, both rows in the main table with their sha256). Files
`outputs/56b_cyp1a2_corrected/submission_candidate_s050.csv` (sha256 `a1da00e2…`) and
`submission_candidate_s100.csv` (sha256 `b2c53f13…`), with server receipts for both sends.

---

## Preserved artefacts from a superseded notebook

`outputs/55_superseded_preserved/` holds the backing CSVs for findings 1–4 and `bugs_found.json`,
a record of five bugs found and fixed in the superseded notebook's own code. The two worth
carrying forward are prose rather than code: both were sentences written before the numbers
existed, contradicted by the numbers when they arrived, and both happened to flatter the method
being used. They were caught by reading executed output against the prose, not by any assertion.

The notebook itself and its `outputs/55_candidate_audit/` directory were deleted on 2026-10-22: it
answered a question about the correction machinery when what was needed was a comparison of
candidates. The four findings above are the part that stands independently.
