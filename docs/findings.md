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

## Preserved artefacts from a superseded notebook

`outputs/55_superseded_preserved/` holds the backing CSVs for findings 1–4 and `bugs_found.json`,
a record of five bugs found and fixed in the superseded notebook's own code. The two worth
carrying forward are prose rather than code: both were sentences written before the numbers
existed, contradicted by the numbers when they arrived, and both happened to flatter the method
being used. They were caught by reading executed output against the prose, not by any assertion.

The notebook itself and its `outputs/55_candidate_audit/` directory were deleted on 2026-10-22: it
answered a question about the correction machinery when what was needed was a comparison of
candidates. The four findings above are the part that stands independently.
