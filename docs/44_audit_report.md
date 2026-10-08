# Audit of notebook 44 — correctness of the pre-registered negative

Read-only audit, 2026-10-08. Nothing retrained, no notebook or script modified, nothing written
outside this file. Every figure below was recomputed from the per-fold saved files rather than read
from the notebook's narrative.

---

## Verdict

**The headline negative survives.** Every primary figure reproduces exactly from
`outputs/44_qm_molecular_block/run_scores.csv`; the effect-size bars reproduce exactly from
notebook 37's own file, which is dated **2026-10-02**, five days before this work and untouched; no
pre-registered arm carries a forbidden column; all 19 arms are complete at 25 folds and every primary
comparison is paired on the same 25 cells in the same order; the decorrelation is confirmed
**residual**, not raw-prediction; the chemprop sanity gate is bit-exact at 0.0 on all four isoforms.
There is no case in which the notebook's prose diverges from its saved files — which is worth stating
explicitly, since this project has caught four such divergences before.

**The single weakest link is the freeze artefact, and it is weaker than the notebook admits.**
`prereg.json`'s mtime is `2026-10-07 23:18:30`, which post-dates **459 of 476** prediction files,
**all 228** tier-1 prediction files, and **all 44** feature-block CSVs. Nothing is committed to git, so
there is **zero independent commit-order evidence**. The only support for its contents being the
originals is (a) the `written_at_utc` field inside the file, which is self-reported, and (b) an
in-session comparison against values printed earlier in the same session. Both are exactly the kind
of evidence this audit was told not to trust.

Two things materially limit that damage, and both are independently checkable:

- **The bars cannot have been tuned to results.** They recompute to the digit as 2× each isoform's
  mean paired ST-RAE SE from `outputs/37_butina_5x5_complete/paired_tests_vs_plain_25fold.csv`, mtime
  **2026-10-02 09:43**. A file five days older than the experiment cannot encode its outcome.
- **The rebuild did not alter the training data.** Refitting four arm/fold cells from the *current*
  feature CSVs reproduces predictions written at 21:44–22:05, before the 23:18 rebuild, to
  **8.882e-16** — machine epsilon from the CSV round-trip. The arms were trained on exactly these bytes.

So what is lost is the *provenance* of the criteria text, not the data and not the bars. If the
criteria themselves matter to a reader — and for a pre-registered negative they do — that gap cannot
be closed from what was saved.

### What invalidates the headline: nothing found.
### What is a detail worth noting: four items, below (A2, A3, A4, B1).

---

## Item-by-item

### 1. Was the pre-registration frozen before results? — **FAIL (artefact), PASS (substance)**

| evidence | result |
|---|---|
| `prereg.json` mtime | 2026-10-07 23:18:30 |
| prediction files older than it | **459 / 476** |
| tier-1 prediction files older than it | **228 / 228** |
| feature CSVs older than it | **44 / 44** |
| tracked in git | **0 files** — no commit-order evidence exists |
| contents vs what the plan promised | all 11 promised elements present (see below) |

Contents check — all **PASS**: `RDKIT_MATCH18` frozen by name (18 cols); `RDKIT_MATCH3_ELEC` =
`[fr_aniline, FractionCSP3, NumAromaticRings]`; `RDKIT_MATCH3_NEW3` =
`[NumHeteroatoms, BalabanJ, PEOE_VSA8]`; 20 draw seeds = 44000–44019; per-isoform bars present;
`criterion_I` with I.1/I.2/I.3; `criterion_II`; the pre-declared negative; `bh_family_rule`; the
9-arm display subset; the model-selection rule.

`decisions_taken.json` **D4** records the rewrite and states what is not recoverable. That is the
right disclosure, but the notebook's own Part 0 prints `prereg.json`'s self-reported timestamp
without flagging that it post-dates training — a reader of the notebook alone would not learn this.

### 2. Was the model type selected on the control block only? — **PASS**

- The filter is `rs["arm"].str.startswith("rdkit2d__")`. Verified by enumeration that this matches
  `rdkit2d__{model}` and **excludes** `rdkit2d_qm__`, `rdkit2d_3d__`, `rdkit2d_new3__`,
  `rdkit2d_elec__` and `qm_full__`. Provably control-only.
- `model_selection.json` was written **21:38:57**, inside the single unattended `--tier1` process
  (21:23–21:41), between phase A and phase B. No human could intervene mid-run.
- Recorded scores: ridge 0.7559 / rf 0.7899 / lightgbm 0.8202 / xgboost 0.8626.

**Qualification worth stating:** QM-arm tier-1 predictions were written **21:32–21:36**, i.e. they
existed on disk *before* the 21:38:57 selection. Blindness therefore rests on the selection being a
deterministic code path in an unattended process, not on the ordering. That is a sound basis, but
had anyone intervened there would be no saved evidence either way.

### 3. Were the BH corrections applied within the specified families? — **MIXED**

| family the plan specified | what ran |
|---|---|
| I.1 / I.2 / I.3 each own family of 4 | **PASS** — one `bh()` call per test inside the `for tid` loop, 4 p-values each |
| `plain_qm` own family of 4 | **PASS** — one `bh()` call, 4 p-values |
| notebook 37's `p_BH` read, never recomputed | **PASS** — the notebook never opens `paired_tests_vs_plain_25fold.csv` at all; the bars arrive via `prereg.json`. Its only read of 34/35/37 `run_scores.csv` pulls per-fold ST-RAE for a *new* `plain` vs `plain_qm` test, not 37's verdicts |
| secondary arms share one family across arm × isoform | **FAIL — never ran.** There are exactly **two** `bh()` call sites in the entire notebook. No secondary family was constructed |

The notebook contains only two BH adjustments. Everything outside the three primaries and `plain_qm`
is reported on **raw p-values**.

### 4. Were the bars derived from notebook 37's saved output? — **PASS, exactly**

| isoform | my recomputed SE | `prereg` SE | 2×SE | `prereg` bar |
|---|---|---|---|---|
| CYP1A2 | 0.0117 | 0.0117 | 0.023 | 0.023 |
| CYP2C9 | 0.0104 | 0.0104 | 0.021 | 0.021 |
| CYP2D6 | 0.0138 | 0.0138 | 0.028 | 0.028 |
| CYP3A4 | 0.0080 | 0.0080 | 0.016 | 0.016 |

10 arms per isoform, mean of the `se` column. Source file mtime 2026-10-02. The "30% of its bar"
framing is correct.

---

## The cancellation finding — **post-hoc in its testing, and the notebook does not say so**

This is finding **A2**, and it is the one I would weight down most after the freeze.

- `prereg.json` declares `QM_ELEC`, `RDKIT2D_ELEC` and `RDKIT_MATCH3_ELEC` as **blocks**, and
  `RDKIT2D_ELEC` is in the confirm spec. So the arm was planned in advance — revision 2's R2.1 added
  it specifically to make an I.1 PASS attributable.
- But **no electronic arm appears anywhere in `criterion_I`**. Verified: the strings `rdkit2d_elec`,
  `qm_elec`, `RDKIT2D_ELEC`, `QM_ELEC` are all absent from it. Declared tests are I.1, I.2, I.3, the
  I.4 sensitivity, I.5 and I.6. There is **no criterion, no BH family and no stated test** for the
  electronic/shape split.
- The cancellation table carries columns `elec, elec_p, elec_folds, shape3d, shape_p, shape_folds, …`
  and **no `p_BH`**. The significance quoted in the headline is **raw and uncorrected, borrowed from
  no declared family**.

So the honest status is: *pre-registered arms, post-hoc comparison*. The notebook presents it as
Part 5, the lead finding, with no post-hoc fence — while Part 8's group-(c) probe *is* properly
fenced. That inconsistency is the problem, not the finding itself.

**How much it matters, by isoform:** CYP2C9's halves are p = 7.2e-06 and 2.5e-04 and would survive
any plausible correction, so the cancellation claim there is robust to the missing family. CYP3A4's
are p = 0.084 and 0.054 — they would survive no correction at all, and the notebook calls CYP3A4
"the same shape". The CYP3A4 half of that claim should be read as illustrative only, which the
notebook does say ("consistent with … rather than independent evidence for it") but says after
having led with "on CYP2C9 **and CYP3A4**".

---

## Recomputed headline figures — all agree with the notebook

Primary tests, independently recomputed; `notebook agrees` compares `mean_diff` and `p_BH` to
`primary_tests.csv` at 1e-9:

| test | isoform | mean_diff | se | p_BH | dz | folds | bar | % of bar | agrees |
|---|---|---|---|---|---|---|---|---|---|
| I.1 | CYP1A2 | **−0.00690** | 0.00204 | **0.00983** | −0.677 | 19/25 | 0.023 | **30%** | yes |
| I.1 | CYP2C9 | −0.00027 | 0.00140 | 0.85032 | −0.038 | 13/25 | 0.021 | 1% | yes |
| I.1 | CYP2D6 | −0.00391 | 0.00225 | 0.18947 | −0.348 | 16/25 | 0.028 | 14% | yes |
| I.1 | CYP3A4 | +0.00056 | 0.00115 | 0.83887 | +0.098 | 11/25 | 0.016 | 4% | yes |
| I.2 | CYP1A2 | −0.00398 | 0.00197 | 0.07319 | −0.404 | 20/25 | 0.023 | 17% | yes |
| I.2 | CYP2C9 | **−0.00512** | 0.00119 | **0.00100** | −0.859 | 20/25 | 0.021 | **24%** | yes |
| I.2 | CYP2D6 | −0.00014 | 0.00118 | 0.90746 | −0.023 | 13/25 | 0.028 | 0% | yes |
| I.2 | CYP3A4 | −0.00169 | 0.00084 | 0.07319 | −0.404 | 19/25 | 0.016 | 11% | yes |
| I.3 | all four | −0.00385 … +0.00093 | — | 0.235–0.787 | ≤0.40 | — | — | ≤17% | yes |

- **Largest favourable primary effect: I.1 / CYP1A2 at 0.00690 against a 0.023 bar = 30%.** It is
  also the largest in magnitude across all 12 cells. The largest *standardised* effect is
  I.2 / CYP2C9 at dz = −0.859. 0 of 12 cells clear a bar.
- **Nested blend:** 11-arm column reproduces notebook 37's published `raw_nested` at **0.0e+00** on
  all four isoforms. Improvements +0.00366 / +0.00094 / 0.00000 / 0.00000 against required
  0.0252 / 0.0225 / 0.0280 / 0.0160. All fail.
- **`plain_qm`:** −0.01104 / −0.00267 / −0.00364 / −0.00366, all `p_BH` 0.8236, |dz| ≤ 0.209,
  15/14/14/12 of 25. Matches `plain_qm_test.csv`. Its residual correlation with `plain` is
  0.9499 / 0.9337 / 0.9565 / 0.9281 — above the re-seed floor on all four, i.e. not decorrelated.
- **Residual vs raw:** the saved `decorrelation_new_arms.csv` matches the **residual**
  (prediction − true) computation on all 16 spot-checked cells and matches the raw-prediction
  computation on **none**. The distinction is load-bearing: `qm_3d__ridge` on CYP2D6 reads 0.845
  residual against **0.048** raw. CLAUDE.md's rule (Dietterich 2000) is satisfied.

The notebook's own methods-result cell states "the largest effect is 0.0069 against a bar of 0.023"
and fractions of 30% and 24%, which match. No divergence found anywhere.

---

## Did the plan's verification gates run?

| gate | result |
|---|---|
| chemprop sanity vs notebook 05's stored `repeat0_fold0` | **PASS** — `abs_diff` 0.0 on all four isoforms, bit-exact; `sanity_verdict.json` `pass: true` |
| fold-file sha256 `52f39b1a…` | **PASS** — file hashes to the expected value; asserted in `Ctx.__init__`, which the notebook calls, and the notebook ran error-free |
| `OMP_NUM_THREADS` unset and asserted | **PASS** — asserted in both the driver (line 49) and the notebook's Part 0 |
| `KMP_DUPLICATE_LIB_OK` reported, not asserted | **PASS** — and the printed output demonstrates the documented behaviour: `None` at kernel start, `'True'` after imports |
| scope check at 0 violations | **PASS** — 4 tests (the plan said two; four were run), 16 protected paths, 0/0/0/0 |
| cache-hit re-run → byte-identical `run_timings.csv` | **FAIL — never ran.** Both confirm logs show **0 cache-hit lines** against 450 and 25 real-run lines. No full replay was performed, so neither the cache-skip contract nor the `run_timings.csv` anti-clobber guard was exercised for the confirm run |
| no arm missing folds | **PASS** — 19 arms, all 25 units, 25 distinct (repeat, fold) cells |
| every primary comparison paired on the same 25 cells | **PASS** — identical cells in identical order for all four increment arms × four isoforms; one seed per cell |

---

## Also checked

- **Forbidden columns — PASS.** None of the 35 pre-registered blocks carries a group (b)/(c)/(d)
  column, verified both from the `prereg.json` manifest and from the CSV headers on disk. **`n_basis`
  appears in no pre-registered block.** The 9 `PH_*` post-hoc blocks do carry group (c) columns, by
  design and correctly labelled.
- **`decisions_taken.json` — 6 decisions.** Two touch the primaries: **D1** (M* selection, which
  determines the model every primary uses — audited clean at item 2) and **D6** (the Ipc check, which
  bears on the comparator; its conclusion that the clip *helps* `rdkit2d` on 4/4 isoforms means any
  residual scaling effect works against the QM arms, not for them). **D2/D3** are post-hoc and fenced.
  **D4** is the prereg rewrite. **D5** is an attribution fix.
- **Execution record — PASS.** 25 code cells, execution counts contiguous 1…25, 0 errors, 0 cells
  unexecuted. The final state is a single clean in-order run. The notebook is pure reporting — every
  decision that could be contaminated (M*, bars, criteria, arm list) was fixed by the scripts before
  the notebook existed — so out-of-order execution could not have let a later result inform an
  earlier decision even in principle.

### Pre-registered analyses missing from the notebook — **A3**

Four declared items never appear:

| declared | status |
|---|---|
| **I.4** sensitivity, `qm_full` vs `rdkit_match18` | arms ran; comparison **absent** from the notebook |
| **I.5** Spearman co-reported on the primaries | **never computed at all.** `Spearman_R` appears nowhere in the notebook's code; `arm_folds` defaults to ST-RAE |
| **I.6** model-class robustness on the primaries | **absent.** The increment arms were run at ridge only; `QM_FULL` *standalone* was run at four classes, which is a different comparison |
| DFT-vs-GFN2 cost attribution (a *named* secondary question in `prereg.json`) | arms ran; comparison **absent** |

**I.5 is the substantive one.** It was pre-registered precisely because Spearman is the one CV
quantity whose board counterpart no affine correction can move. Its absence means the notebook cannot
say whether the QM block moved *ranking* while failing to move ST-RAE — a question its own
pre-registration raised and which the saved per-fold files can still answer.

---

## What could not be determined from what was saved

1. **The original `prereg.json`'s contents.** The first builder run's output went to the terminal, not
   a log; nothing is in git. The in-session content comparison is the only evidence, and it is not
   independent. The bars and the feature data *are* independently confirmed; the criteria text is not.
2. **Whether anyone intervened between tier-1 phase A and the M\* selection.** Procedurally impossible
   within a single unattended process, and the code is provably control-only, but there is no saved
   record of non-intervention as such.
3. **Whether the cache-skip contract holds for the confirm run.** It was never exercised, so its
   correctness is asserted by code inspection only.

---

# Remediation, 2026-10-08

Four items from the audit, all run on per-fold predictions already on disk. **No model retrained and
no prediction recomputed.** Notebooks 32 and 40 were never opened. `OMP_NUM_THREADS` and
`KMP_DUPLICATE_LIB_OK` left unset. Notebook 44 goes from 43 cells to **54 (30 code), 0 errors,
execution counts contiguous 1-30, scope check 0 violations on all four tests.**

## Blocking verification, before anything was computed

The saved per-fold predictions are `outputs/44_qm_molecular_block/predictions/{arm}__r{R}f{F}.csv`,
**475 files** = 19 arms x 25 cells (plus one sanity file), each carrying `Molecule_Name`, `inchikey`
and all four `*_pIC50_direct_inhibition` columns. The fold assignment inside them was checked against
`data/folds/cv_folds_butina_5x5.csv`, which hashes to
`52f39b1a1cea9200b58dce0a061d4a21ce94d9d58ec77067cf6508b0f4b3e87e` as expected: for 5 arms x 25
cells the compound set in each file is exactly the partition's own membership for that cell,
**0 mismatches**. The repeated-measures pairing is valid and the work proceeded.

`Spearman_R` turned out to be **already saved per fold** in `run_scores.csv`, as the vendored
scorer's own bootstrap mean. I.5 therefore reads the same column from the same source the ST-RAE
cells read, which keeps it on identical footing rather than introducing a second estimator.

## Item 1: I.5 Spearman, computed (new Part 6b)

Pre-registered in specification, **run late with the ST-RAE result already known**, and labelled
that way in the cell output. Sign convention stated explicitly in the output and deliberately not
inherited from ST-RAE: Spearman is higher-is-better, so positive = QM arm better.

| test | isoform | Spearman diff | p_BH | dz | folds better | bar | % of bar |
|---|---|---|---|---|---|---|---|
| I.1 | CYP1A2 | **+0.00575** | **0.0252** | +0.54 | 19/25 | 0.0183 | 31% |
| I.1 | CYP2C9 | +0.00183 | 0.3167 | +0.24 | 16/25 | 0.0160 | 11% |
| I.1 | CYP2D6 | **+0.01123** | **0.0055** | +0.72 | 18/25 | 0.0212 | 53% |
| I.1 | CYP3A4 | −0.00039 | 0.6513 | −0.09 | 9/25 | 0.0100 | — |
| I.2 | CYP1A2 | **+0.00659** | **0.0005** | +0.92 | 20/25 | 0.0183 | 36% |
| I.2 | CYP2C9 | **+0.00463** | **0.0046** | +0.65 | 20/25 | 0.0160 | 29% |
| I.2 | CYP2D6 | **+0.00351** | **0.0023** | +0.74 | 21/25 | 0.0212 | 17% |
| I.2 | CYP3A4 | +0.00106 | 0.0576 | +0.40 | 18/25 | 0.0100 | 11% |
| I.3 | CYP1A2 | **+0.00644** | **0.0008** | +0.88 | 21/25 | 0.0183 | 35% |
| I.3 | CYP2C9 | **−0.00301** | **0.0240** | −0.51 | 10/25 | 0.0160 | (a loss) |
| I.3 | CYP2D6 | **+0.00361** | **0.0018** | +0.76 | 19/25 | 0.0212 | 17% |
| I.3 | CYP3A4 | −0.00004 | 0.9141 | −0.02 | 15/25 | 0.0100 | — |

**Spearman shows BH-significant ranking gains on 7 of 12 cells, against ST-RAE's 2 of 12.** The two
metrics genuinely disagree about how consistent the signal is, and the notebook reports that rather
than reconciling it. **No cell is material.** The largest gain is I.1/CYP2D6 at +0.01123, 53% of its
bar. One BH-significant *decline* is reported too: I.3 on CYP2C9, −0.00301.

**The bar derivation was the one real judgement call, and the first attempt was unsound.** Deriving
2 x the paired SE of *the comparison being tested* gives 0.0012 to 0.0034 and would call **7 of 12
cells material**. That is circular: the increment arms differ from `rdkit2d` by 18 columns in 235, so
that SE measures the precision of a near-identical pair and makes almost any detectable difference
"material" by construction. The pre-registered ST-RAE bar was **not** built that way: it is 2 x the
mean paired SE across notebook 37's **ten different arms** against plain, an independent noise scale.
Rebuilding the Spearman bar the same way gives **0.0183 / 0.0160 / 0.0212 / 0.0100**, which is 5.4x
to 8.3x larger and calls **none** material. Both derivations and the resulting flip are printed in
the notebook, because the verdict turns on the choice and a reader is entitled to see that.

`plain_qm` on Spearman: **+0.01523 on CYP1A2**, p_BH 0.0517, 18/25 folds, **83% of its bar** - the
closest anything in this notebook comes to materiality on either metric, and named as the one cell a
future design might target. CYP2C9/CYP2D6/CYP3A4 are +0.0046/+0.0051/+0.0034, none significant.

**Bearing on notebook 38: none.** That notebook established CYP3A4's remaining gap is ranking, not
placement: Spearman 0.8020 against 04b's 0.8329, a gap of **0.0309**. The QM block moves CYP3A4
ranking by at most **+0.00106**, roughly 30x too small, and two of the three tests point the wrong
way on that isoform. The QM block does not address that gap.

## Item 2: Part 5 fenced

Heading changed to "**Part 5 - POST-HOC**", with a fence block on Part 8's template recording that
`prereg.json` declares the electronic arm as a *block* and puts it in the confirm spec, but that **no
electronic arm appears anywhere in `criterion_I`**: no criterion, no pre-declared BH family, no
stated test. The decomposition was formulated after the block-level null was read.

**CYP3A4 is dropped from the headline claim.** Its halves are only marginal on raw p (0.084, 0.054)
and neither survives correction; the electronic half reaches p_BH 0.0565 in the declared family. The
cell now says the earlier "CYP2C9 and CYP3A4" framing was not supportable, and retains CYP3A4 as
illustration of the pattern rather than evidence for it. **CYP2D6 is also dropped** - see item 4.
The claim now rests on CYP2C9 alone, where the electronic half is p_BH **1.44e-05** and the shape
half p_BH 0.0010 from the I.2 primary family. The notebook states plainly that the two halves sit in
**different families**.

## Item 3: prereg rewrite disclosed in the notebook

A markdown limitation block now sits immediately after the title, before the setup cell. It states
the mtime, the 459/476, 228/228 and 44/44 counts, the absence of git history, and that every "frozen
before any arm ran" claim rests on a self-reported field plus a session-internal comparison. It then
gives the two independent mitigations - bars reproducing from a file dated 2026-10-02, training data
reproducing to 8.882e-16 - and states plainly that **the provenance of the criteria text itself is
not recoverable**. Written as a limitation to weigh, not a defence.

## Item 4: secondary BH family run (new Part 4b)

The family is enumerated from the plan's own Part 3.2 secondary table plus revision 2's electronic
increment: **9 comparisons x 4 isoforms = 36 tests, BH across all 36**. `prereg.json` codifies only
I.4 explicitly, so family *membership* is a remediation decision and is labelled as one; what is not
a judgement call is that correcting within some declared family beats quoting raw p.

It is placed **before** Part 5 because Part 5's electronic half is a member of it. That required
lifting `rs`, `arm_folds`, `paired_test` and `bh` out of Part 5's first cell into a shared helpers
cell, which is an improvement in structure regardless.

**28 of 36 significant on raw p, 27 of 36 after BH. Exactly one conclusion changes:**
`rdkit2d_elec` vs `rdkit2d` on **CYP2D6**, p_raw 0.0440 to p_BH 0.0565. That is the CYP2D6 electronic
half of Part 5, which is why Part 5 now rests on CYP2C9 alone. Everything else is robust, so the
original raw-p reporting overstated nothing except that one cell.

Three secondary cells clear their ST-RAE bar: `qm_3d` vs `rdkit_match3_new3` on CYP2C9 (−0.0285) and
CYP3A4 (−0.0433), and the `hirshfeld_min` sensitivity on CYP3A4 (−0.0174). All three are
weak-against-weak: `qm_3d` and `qm_full` sit at 0.93-0.97 macro against the `rdkit2d` comparator at
0.772, so beating a 3-column cheap control materially does not make either competitive. prereg I.4
says exactly that, and it holds.

**Cost attribution, now adjusted.** The original run called `qm_dft3` vs `qm_gfn23` a tie on macro.
Per isoform it is **mixed-sign and mostly significant**: DFT better on CYP2C9 (−0.0097, p_BH 2.6e-04)
and CYP3A4 (−0.0108, 1.4e-04), **GFN2 better on CYP2D6** (+0.0248, p_BH < 1e-5), CYP1A2 not
significant. Every effect is under its bar and both arms sit at ~0.995 macro, i.e. naive-predictor
level, so neither level of theory carries usable standalone signal. The informative version of the
question remains Part 2's feature-level increment test (+0.14 to +0.23 out-of-sample R2).

## Did the conclusion move?

**No.** Criterion I still fails on all three primary tests, Criterion II still fails, and
Criterion III - the pre-declared negative - remains the recorded conclusion. I.5 adds a real,
consistent, BH-significant **ranking** signal that ST-RAE did not show, on 7 of 12 cells rather than
2, and it is immaterial on an independently derived bar. So the negative is now tested on **both**
metrics rather than one, and it holds on both, with the ranking side closer to the bar than the
ST-RAE side ever came (53% and 83% against 30%).

The honest residual: I.5 narrows the ranking question rather than closing it. A reader who derives
the Spearman bar from the tested comparison reaches the opposite verdict on 7 cells, and `plain_qm`
on CYP1A2 sits at 83% of its bar with p_BH 0.0517. Neither fact rescues the block, and neither is
buried.
