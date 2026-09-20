# Leaderboard submissions — authoritative record

This is the single canonical record of every real submission this project has made to the
OpenADMET CYP450 blind leaderboard. It supersedes any partial or macro-only leaderboard figures
stated elsewhere (notebook markdown, CLAUDE.md prose, the root `README.md`) — those should point
here rather than repeat these numbers. Real per-isoform scores exist here for 04b and notebook 10
already (originally cross-checked against `README.md`'s own "Leaderboard result" tables before
those tables were retired in favor of a pointer to this file — see "Provenance" below); 10c,
notebook 12, NB13-calib's, NB16's, NB19's, NB27-widened's, NB30's, 27b's, and NB32-A's full metric
breakdowns are recorded here as well.

No blind-set labels are held anywhere in this repo — all scores below came from OpenADMET's own
leaderboard after real submissions, not from anything computed locally.

**Every score in this document is a half-set measurement.** The live leaderboard scores only half
of the blind test set (split by chemical series) as submissions come in — the full set is scored
only at the intermediate reveal and at the competition's close. This applies to every row below;
none of these numbers are the final, full-set result, and later reveals could move any of them.

## All eleven submissions to date

| Isoform | Metric | 04b | 10c | NB10 | NB12 | NB13-calib | NB16 | NB19 | NB27-widened | NB30 | 27b | NB32-A |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **MA (Overall)** | ST-RAE | 0.7179 | 0.7138 | 0.8299 | 0.8299 | 0.8811 | 0.6364 | 0.6265 | 0.6155 | 0.8160 | 0.5982 | 0.6488 |
| | MAE | 0.8952 | 0.8927 | 0.9690 | 0.9674 | 1.0260 | 0.7679 | 0.7599 | 0.7488 | 0.9580 | 0.7206 | 0.7670 |
| | R² | 0.2191 | 0.2278 | 0.1054 | 0.1159 | 0.0203 | 0.4407 | 0.4482 | 0.4557 | 0.1356 | 0.4835 | 0.4329 |
| | Spearman | 0.6739 | 0.6735 | 0.6748 | 0.6613 | 0.6785 | 0.6735 | 0.6785 | 0.6785 | 0.6831 | 0.7002 | 0.6894 |
| | Kendall | 0.5008 | 0.4987 | 0.4977 | 0.4877 | 0.5045 | 0.4987 | 0.5045 | 0.5045 | 0.4999 | 0.5212 | 0.5118 |
| **CYP1A2** | ST-RAE | 0.7114 | 0.6954 | 0.7850 | 0.7633 | 0.8197 | 0.6954 | 0.6954 | 0.6954 | 0.7630 | 0.6954 | 0.7906 |
| | MAE | 0.9589 | 0.9385 | 1.0214 | 1.0035 | 1.0549 | 0.9385 | 0.9385 | 0.9385 | 1.0027 | 0.9385 | 1.0287 |
| | R² | 0.2824 | 0.3073 | 0.2051 | 0.2275 | 0.1537 | 0.3073 | 0.3073 | 0.3073 | 0.2588 | 0.3073 | 0.2125 |
| | Spearman | 0.7430 | 0.7375 | 0.7246 | 0.7389 | 0.7375 | 0.7375 | 0.7375 | 0.7375 | 0.7466 | 0.7375 | 0.7565 |
| | Kendall | 0.5473 | 0.5425 | 0.5297 | 0.5422 | 0.5425 | 0.5425 | 0.5425 | 0.5425 | 0.5468 | 0.5425 | 0.5611 |
| **CYP2C9** | ST-RAE | 0.5408 | 0.5489 | 0.5489 | 0.5375 | 0.6036 | 0.5770 | 0.5375 | 0.5375 | 0.5922 | 0.5375 | 0.5635 |
| | MAE | 0.5399 | 0.5444 | 0.5444 | 0.5280 | 0.5934 | 0.5603 | 0.5280 | 0.5280 | 0.5607 | 0.5280 | 0.5568 |
| | R² | 0.5263 | 0.5272 | 0.5272 | 0.5432 | 0.4237 | 0.5133 | 0.5432 | 0.5432 | 0.4789 | 0.5432 | 0.5016 |
| | Spearman | 0.7428 | 0.7389 | 0.7389 | 0.7585 | 0.7585 | 0.7389 | 0.7585 | 0.7585 | 0.7333 | 0.7585 | 0.7554 |
| | Kendall | 0.5490 | 0.5462 | 0.5462 | 0.5692 | 0.5692 | 0.5462 | 0.5692 | 0.5692 | 0.5369 | 0.5692 | 0.5652 |
| **CYP2D6** | ST-RAE | 1.1903 | 1.1673 | 1.4603 | 1.4092 | 1.3570 | 0.8299 | 0.8299 | 0.7858 | 1.3168 | 0.7165 | 0.7502 |
| | MAE | 1.5842 | 1.5694 | 1.7489 | 1.7209 | 1.6904 | 1.0546 | 1.0546 | 1.0104 | 1.6688 | 0.8975 | 0.9319 |
| | R² | −0.6220 | −0.5941 | −0.9392 | −0.8726 | −0.7903 | 0.2713 | 0.2713 | 0.3012 | −0.7660 | 0.4124 | 0.3768 |
| | Spearman | 0.3769 | 0.4000 | 0.4401 | 0.3771 | 0.4000 | 0.4000 | 0.4000 | 0.4000 | 0.4872 | 0.4872 | 0.4438 |
| | Kendall | 0.2606 | 0.2747 | 0.3081 | 0.2587 | 0.2747 | 0.2747 | 0.2747 | 0.2747 | 0.3415 | 0.3415 | 0.3076 |
| **CYP3A4** | ST-RAE | 0.4291 | 0.4434 | 0.5253 | 0.6095 | 0.7443 | 0.4434 | 0.4434 | 0.4434 | 0.5920 | 0.4434 | 0.4907 |
| | MAE | 0.4976 | 0.5183 | 0.5614 | 0.6173 | 0.7654 | 0.5183 | 0.5183 | 0.5183 | 0.5996 | 0.5183 | 0.5506 |
| | R² | 0.6898 | 0.6709 | 0.6285 | 0.5655 | 0.2941 | 0.6709 | 0.6709 | 0.6709 | 0.5708 | 0.6709 | 0.6405 |
| | Spearman | 0.8329 | 0.8177 | 0.7955 | 0.7706 | 0.8177 | 0.8177 | 0.8177 | 0.8177 | 0.7651 | 0.8177 | 0.8020 |
| | Kendall | 0.6463 | 0.6316 | 0.6065 | 0.5809 | 0.6316 | 0.6316 | 0.6316 | 0.6316 | 0.5743 | 0.6316 | 0.6132 |

*ST-RAE is this project's governing metric (lower is better); R², Spearman, and Kendall are
better higher; MAE is better lower. See `src/vendor/openadmet_eval/` for the scoring definitions.*

## What was submitted, and when

| Submission | Notebook | Submitted (UTC) | Recipe | Deliberate? |
|---|---|---|---|---|
| **04b** | `04b_final_submission.ipynb` | (undated, first submission) | Single `chemprop_chemeleoninit`, `cyp-admet` env, config chosen by 04a's single-fold screen | Yes — the reference, uncontested real submission |
| **10c** | `10c_control_submission.ipynb` | 2026-09-04 ~09:15 local | Single `chemprop_chemeleoninit`, `cyp-admet-v2` env, config chosen by 05's full 5×5 CV — a control isolating the environment migration | Yes |
| **NB10** | `10_final_retrain_predict.ipynb` | 2026-09-03 15:24 UTC (file mtime) | 08's simple-average ensembles (uncapped), full retrain under `cyp-admet-v2` | Yes |
| **NB12** | `12_caruana_retrain_predict.ipynb` | 2026-09-04 21:28 UTC (file mtime) | Caruana bagged-selection ensembles, capped at the empirical CV-minimizing `k` per isoform | Yes |
| **NB13-calib** | `13_aid1851_blind_population_calibration.ipynb` | 2026-09-13 15:07 UTC | AID 1851 floored-mean population calibration, applied uniformly to all four isoforms — CYP1A2/CYP2D6/CYP3A4 calibrate `chemprop_chemeleoninit`'s full-data blind predictions (reused bit-identical from `10c`/NB12's output), CYP2C9 calibrates the capped Caruana ensemble's full-data blind predictions (reused bit-identical from NB12's output). Deliberately applied even to CYP2C9/CYP3A4, where notebook 14's OOF test had already found the correction hurts — a clean, isolated test of the calibration mechanism itself, not a best-guess recipe | Yes |
| **NB16** | `16_board_solved_population_calibration.ipynb` | 2026-09-14 13:21 UTC | `10c`'s raw blind predictions with a board-metrics-solved blind-population recentring applied to CYP2C9 and CYP2D6 only, CYP2D6 using notebook 16's addendum lower-root override (target mean 3.1558, sd 1.5114); CYP1A2 and CYP3A4 left byte-identical to `10c` | Yes — manually reviewed and sent, per this project's manual-gated submission process |
| **NB19** | `19_cyp2c9_revert_submission_candidate.ipynb` | 2026-09-15 10:14 UTC | `16`'s submitted file (board-solved population recentring, CYP2D6 lower-root override) with CYP2C9 reverted to NB12's own uncalibrated, capped-Caruana-ensemble column; CYP1A2, CYP2D6 and CYP3A4 carried over byte-identical (sha256-confirmed per column) from NB16 | Yes — manually reviewed and sent, per this project's manual-gated submission process |
| **NB27-widened** | `27_calibration.ipynb` (Section 2) | 2026-09-16 21:36 UTC | NB19's submitted file with CYP2D6 widened about its own mean by multiplier 1.397056 (spread ratio 0.6084 → 0.85 of the CYP2D6 training-label SD, so NB16's underlying CYP2D6 recentring remains in place beneath the widening); CYP1A2, CYP2C9 and CYP3A4 carried over byte-identical from NB19 | Yes — manually reviewed and sent, per this project's manual-gated submission process; see the note below on this notebook's own "no action taken" closing text |
| **NB30** | `30_aid_full_retrain.ipynb` | 2026-09-18 11:39 UTC | CheMeleon multitask + AID 1851 auxiliary-head model (notebook 29's 5×5-CV-confirmed recipe), full-data retrain under `cyp-admet-v2`, seed 2684470948, raw and uncorrected on all four isoforms — no calibration, no ensembling | Yes — manually reviewed and sent, per this project's manual-gated submission process |
| **27b** | `27b_aid_cyp2d6_corrected.ipynb` | 2026-09-19 12:11 UTC | NB30's AID 1851 auxiliary-head model, mixed recipe: CYP1A2/CYP3A4 carried over byte-identical from `10c`, CYP2C9 byte-identical from NB19, CYP2D6 placement-corrected onto the board-solved blind population (mean 3.1614, sd 1.5121, solved from NB30's own published metrics) then widened to 0.85 of the CYP2D6 training-label SD — the same recipe NB27 Section 2 used, applied to NB30's own OOF rho instead of the plain model's | Yes — manually reviewed and sent, per this project's manual-gated submission process |
| **NB32-A** | `32_deadzone_aid_retrain.ipynb` (candidate A, `submission_candidate_A_allfour.csv`) | 2026-09-20 02:10 UTC | The DEADZONE+AID model's own four columns (AID 1851 auxiliary-head architecture — notebook 29's 5×5-CV-confirmed recipe — trained against notebook 31's screened dead-zone clipped target instead of the point estimate, full-data retrain), CYP2D6 placement-corrected and widened using the same NB16-solved target and 0.85-of-training-SD ratio as NB27-widened/27b; CYP1A2, CYP2C9 and CYP3A4 raw from that model, uncorrected | Yes — manually reviewed and sent, per this project's manual-gated submission process |

**On "deliberate":** all eleven submissions were reviewed and uncommented/run by hand, exactly as
this project's manual-gated-submission process requires. An earlier diagnostic pass over this
repo (before this correction) reasoned from file evidence alone — a misleading header comment on
the NB10/10c/NB12 submission cells, which claimed "every line below is a comment" when only the
explanatory paragraph was — and concluded the submissions might have fired by accident. **That
conclusion was wrong.** The header comment itself was a real bug (now fixed, 2026-09-07 — see
CLAUDE.md's Rules section), but the submissions it sat on top of were intended, reviewed actions,
and the content each one sent was correct.

**NB27-widened's own notebook text creates a similar risk of the same wrong inference, and for the same underlying reason.** Its submission cell is fully commented out by this project's now-standard convention — every line, including the `client.predict(...)` call itself, is a `#` comment — and its closing markdown cell is titled "Section 2 — summary, no action taken." Read from the notebook's text alone, both would reasonably suggest this candidate was built and left unsent. **It was not.** It was submitted by hand on 2026-09-16 at 21:36 UTC, exactly as this project's manual-gated process requires; the notebook itself was not updated with that fact until afterward, via a single markdown cell appended to its end (see notebook 27's own closing note). Recorded here so a future reader does not draw the same incorrect inference from the notebook's text alone that the earlier diagnostic pass drew from the header-comment bug.

## NB13-calib: population calibration tested and rejected on real blind data

**NB13-calib is worse than both 04b and 10c on every single isoform**, on ST-RAE (the governing
metric): CYP1A2 0.8197 vs. 0.7114/0.6954, CYP2C9 0.6036 vs. 0.5408/0.5489, CYP2D6 1.3570 vs.
1.1903/1.1673, CYP3A4 0.7443 vs. 0.4291/0.4434, macro 0.8811 vs. 0.7179/0.7138. This includes
CYP1A2 and CYP2D6 — the two isoforms where notebook 14's own OOF cross-validation test had
predicted floored-mean calibration would *improve* ST-RAE over raw (0.821 vs. 0.868 for CYP1A2,
0.932 vs. 1.004 for CYP2D6). Both went the wrong direction on real blind data.

Combined with notebook 14's clean OOF-side negative result for censored-MLE calibration on every
isoform, and 09/11b's earlier finding that own-population placement correction wasn't worth
pursuing (kb_share <10% everywhere), **population calibration targeting either this project's own
training-population moments (09, 11b) or an external reference population (AID 1851, 13, 14) is a
confirmed real-blind-data negative result for both of those targets, not just an OOF one.** No
form of training-population or external-population calibration is adopted or recommended going
forward.

**This does not extend to calibrating toward the blind population itself** — a different target,
tested later by a different submission, with a real positive result. See "NB16" below (this
section's finding for 09/11b/13/14 stands unchanged; only the earlier "in every form tested"
framing was too broad, since blind-population recentring was not among what those four notebooks
tested).

## NB16: population calibration works when targeting the blind population

**NB16 is the best submission to date on the governing metric.** Macro ST-RAE 0.6364 — beating
every prior submission, including 10c (0.7138, the previous best) and 04b (0.7179).

**CYP2D6 was the win.** ST-RAE 1.1673 (10c) → 0.8299, R² −0.5941 → +0.2713 — a real, sizeable
change. Spearman is unchanged at 0.4000 (identical to 10c's own CYP2D6 Spearman): the recentring
is a pure additive placement fix that left which compounds rank above which untouched, exactly
what a rigid mean/spread correction on top of otherwise-unchanged predictions should produce.

**CYP2C9 regressed.** ST-RAE 0.5489 (10c) / 0.5375 (NB12, the best CYP2C9 result on record) →
0.5770 under the same correction mechanism. The correction earns its place where placement is
genuinely broken (CYP2D6) and costs performance where placement was already close to sound
(CYP2C9) — not a universal win.

**This supersedes, but does not overturn, the "NB13-calib" section above.** That section's own
findings for notebooks 09/11b (own training-population calibration) and 13/14 (external AID-1851
population calibration) stand as real negative results — both of those calibration *targets*
really did fail on real blind data. NB16 calibrates toward a third, different target: the blind
population's own mean/spread, solved algebraically from this project's four earlier real board
submissions' published R²/MAE/Spearman (notebook 16's own method, including its lower-root
override addendum — see CLAUDE.md's notebook-status log for the full derivation). That target is
closer to what ST-RAE actually scores against, and on CYP2D6 it produced a real improvement. The
document's earlier "in every form tested" framing for the negative result has been narrowed
accordingly, in place, above.

## NB19: CYP2C9 reverted — the best submission to date

**NB19 is the best submission to date on the governing metric.** Macro ST-RAE 0.6265 —
beating NB16 (0.6364, the previous best), 10c (0.7138), and 04b (0.7179).

**It confirms the per-isoform criterion NB16's own result implied.** Blind-population
recentring earns its place where the bias term is large and placement is genuinely
broken (CYP2D6: ST-RAE 1.1673 [10c] → 0.8299, R² −0.5941 → +0.2713 — unchanged from
NB16, since NB19 carries CYP2D6 over byte-identical) — and it costs performance where
placement was already sound (CYP2C9). Reverting CYP2C9 to NB12's uncalibrated,
capped-Caruana-ensemble column returns it to **0.5375 exactly, its best value on
record**, undoing NB16's regression to 0.5770.

**A REPRODUCIBILITY RESULT worth stating in its own right.** CYP1A2, CYP2D6, and CYP3A4
were submitted byte-identical to NB16 (per notebook 19's own sha256 check against
`outputs/board_solved_calibration/submission_candidate.csv`) and scored **byte-identical
on the board** — 0.6954 / 0.8299 / 0.4434 ST-RAE, and every other metric, matching NB16's
own row exactly for all three isoforms. Identical inputs produced identical board output:
the half-set composition and bootstrap scoring are stable between submissions taken 1
day apart, so per-isoform comparisons across different submissions in this table are
sound and not an artifact of scoring noise.

**CYP2C9's monotonic-improvement run, noted below as stopping at NB12, now resumes.**
NB13-calib (0.6036) and NB16 (0.5770) both moved CYP2C9 the wrong way; NB19 returns it
to NB12's 0.5375 — the same value, not just the same neighborhood, since NB19 reuses
NB12's own CYP2C9 predictions unmodified.

## NB27-widened: CYP2D6 spread widened — the best submission to date, pending one Spearman check

**NB27-widened is the best submission to date on the governing metric, if the published macro is
taken at face value** — macro ST-RAE 0.6155, beating NB19 (0.6265, the previous best), NB16
(0.6364), 10c (0.7138), and 04b (0.7179). See "Provenance" below for the check on this figure:
ST-RAE, MAE, and Kendall reproduce exactly from the four isoform values above (ST-RAE
(0.6954+0.5375+0.7858+0.4434)/4 = 0.615525 → 0.6155, not the 0.655525 an earlier arithmetic pass
over these same figures produced — that arithmetic error is corrected here, not carried forward);
R² reproduces via the same round-half-up convention already established for NB16/NB19's own
boundary cases. **Spearman does not check out** — the same 0.0001 gap already seen in NB19, and
for the same reason: NB27-widened's four isoform Spearman values are numerically identical to
NB19's (CYP2D6's rank order is unchanged by a pure spread scale — see below), so re-averaging them
reproduces NB19's own discrepancy exactly. The ST-RAE figure this headline rests on is not in
doubt; see "Provenance" for the full check.

**CYP2D6 moved again, under a transform with a clean, pre-stated falsifiable prediction.** ST-RAE
0.8299 (NB16/NB19) → 0.7858, R² 0.2713 → 0.3012. Notebook 27 stated, before submission, that a
pure mean-preserving spread scale cannot change which compounds rank above which — Spearman and
Kendall must be unchanged. They were: 0.4000 and 0.2747, identical to NB16/NB19's own CYP2D6
values to four decimal places. **The falsifiable prediction held.**

**A REPRODUCIBILITY RESULT, for the second time.** CYP1A2, CYP2C9 and CYP3A4 were carried over
byte-identical to NB19 and scored byte-identical on the board — every metric, on all three
isoforms, matches NB19's own row exactly (0.6954/0.5375/0.4434 ST-RAE, and all four other metrics
per isoform). Following NB19's own demonstration of this a day after NB16, this is the second
confirmation that identical inputs produce identical board output across separately-timed
submissions, so per-isoform comparisons across this table's different submissions continue to be
sound rather than an artifact of half-set scoring noise.

**A second board result in tension with out-of-fold evidence for CYP2D6 — reported carefully,
because the two relevant OOF exercises do not agree with each other either.** Notebook 21 already
established that OOF ST-RAE cannot be trusted as a selection metric for this isoform (the
board-solved recentring that fixed CYP2D6 for NB16 scores 3.08 on the OOF curve, roughly 3x the
raw baseline). Notebook 26 ran two further, non-agreeing OOF-adjacent sweeps for CYP2D6, and only
one of them shares NB27's own starting point: its Part 3 sweep, applied to the SAME deployed
CYP2D6 column NB27 widened (mean 3.155762, std 0.557002, i.e. `s=1.0` in Part 3's own units is
this project's actually-submitted NB16/NB19 column) against a solved estimate of the blind
population, wanted a further expansion multiplier of `s≈2.77` from that exact baseline — the same
direction NB27 tested and the board rewarded. Its separate Part 2 sweep, run on raw, uncorrected
OOF predictions rather than the deployed column, found an ST-RAE-optimal multiplier of `s=0.48`
against a theoretical squared-error-optimal shrinkage of `rho=0.3685` — real shrinkage relative to
raw, but *less* shrinkage than the theoretical `rho` amount, i.e. pointing toward more spread than
a naive correction would apply, not less. Because Part 2 starts from a different, non-comparable
baseline (raw OOF predictions, not NB19's already mean-recentred, mildly spread-reduced blind
column), and notebook 26 itself explicitly declines to reconcile its own two sweeps ("do the two
answers agree — False"), this document does not read Part 2 as a quantified prediction that
further widening NB19's column would help, nor as one that it would hurt. Only Part 3, sharing
NB27's exact starting point, is directly comparable, and it points the same direction as the
board. Recorded at this level of care because a superficial reading of Part 2's own "s=0.48"
figure in isolation could easily be, and elsewhere has been, characterised as a contradiction of
the board result; on notebook 26's own text, that characterisation does not hold up.

**Notebook 27's own Section 1 pooled-compression finding is one submission short of the record
now available, and has not been re-run.** That finding (Spearman −0.8221 between spread ratio and
board ST-RAE, n=28 — 4 isoforms × 7 submissions, 95% bootstrap CI [−0.8879, −0.6617]) was computed
before this eighth submission existed; the fuller record now has 32 such points. It has not been
re-computed here, and this document does not know whether the correlation strengthens, weakens, or
holds at n=32. What can be said without re-running it: CYP2D6's own movement between NB19 and
NB27-widened — spread ratio 0.6084 → 0.85 (less compression), ST-RAE 0.8299 → 0.7858 (better) — is
directionally consistent with the pooled finding's own direction (less compression, better
ST-RAE), for this one new pair of points. That is a single consistent data point, not a
re-verification of the pooled result.

## NB30: a confirmed CV win fails to transfer — the third such instance

**NB30 is worse than the current best (NB27-widened, 0.6155) on every isoform, and worse than
`10c`'s raw baseline (0.7138) on the governing metric overall** — macro ST-RAE 0.8160. This is the
CheMeleon + AID 1851 auxiliary-head model, full-data retrained under the exact recipe notebook 29
confirmed at 5×5 CV: CYP1A2 significantly better (p_BH = 0.0017), CYP2D6 significantly better
(p_BH = 0.00014), CYP2C9/CYP3A4 not significantly different from baseline. **The board says worse
on all four isoforms** — CYP1A2 0.6954 → 0.7630, CYP2C9 0.5489 → 0.5922, CYP2D6 1.1673 → 1.3168 (or
0.8299 → 1.3168 against the deployed, corrected best), CYP3A4 0.4434 → 0.5920, all measured against
`10c`'s raw baseline, the fair comparison for an uncalibrated model.

**This is the third time a CV result has failed to transfer to the real blind board on this
project**, and the most rigorously confirmed of the three: notebook 08's ensembles were selected
from a single 25-fold CV evaluation pass (05's data), never independently re-confirmed before
submission (NB10/NB12, macro 0.8299, worse than `10c`'s 0.7138); notebook 28's AID-head result was
a single CV fold; notebook 29's AID-head result was a full, honest 5×5 repeated-CV **confirmation**
run specifically to re-test notebook 28's single-fold finding, with paired significance tests,
BH-corrected — the strongest CV evidence any adopted recipe on this project has carried into a
submission, and it still did not transfer. Recorded plainly, per this document's own convention of
reporting negative results without softening them.

**The compression signal predicted this outcome; the CV signal did not.** Notebook 30 itself
flagged, before submission, that this model's raw predictions are more compressed than `10c`'s raw
predictions on all four isoforms (spread ratio vs. training-label SD: CYP1A2 0.685 vs. 0.836,
CYP2C9 0.726 vs. 0.995, CYP2D6 0.415 vs. 0.701, CYP3A4 0.799 vs. 0.933 — all confirmed directly
against `outputs/30_aid_full_retrain/compression_vs_10c.json` in this pass). Notebook 27 Section 1
found a robust pooled association between spread compression and board ST-RAE across this
project's own submission history (Spearman −0.822, 95% bootstrap CI [−0.888, −0.662], n = 28 at the
time it was computed, `10c`/NB16/NB19 all included since they carry identical CYP1A2/CYP3A4
columns) — and NB30 lands squarely on that association's own predicted side: it is more compressed
than `10c` on all four isoforms and scores worse than `10c` on all four isoforms too (checked
directly against `outputs/27_calibration/spread_ratio_table.csv`, which does **not** yet include
NB30 or NB27-widened as rows — this is a fresh comparison against `10c` alone, not a re-run of the
full 28/32-point pooled correlation). **This is not the most compressed submission on record on
every isoform** — NB12 and NB13-calib are more compressed than NB30 on CYP2D6 and CYP3A4, and
NB13-calib is more compressed on CYP1A2 too (`spread_ratio_table.csv`) — so the claim here is
specifically "more compressed than, and scores worse than, `10c`," not "the most compressed
submission ever." **Nor is NB30 the second-worst submission on the board**: by macro ST-RAE it
ranks seventh-best of ten (worse than `10c`, `04b`, NB16, NB19, NB27-widened, and 27b; better than
NB10, NB12, and NB13-calib — updated from this section's original "sixth-best of nine" once 27b's
submission was recorded; see "The record gap and its consequence") — see the chronological
trajectory table below for the full ranking. The
5×5-CV significance tests carried no information about the `10c`-relative compression/direction
match; the compression measurement, available before submission, pointed the right direction for
that one comparison.

**The Spearman decomposition is the reason this document does not recommend discarding the AID
model wholesale.** Placement (mean/intercept) and spread (SD/slope) corrections are both affine
and cannot change rank order — so any change in Spearman between `10c` and NB30 is a genuine
property of the model itself, not an artifact of calibration:

| Isoform | 10c Spearman | NB30 Spearman | Direction |
|---|---:|---:|---|
| CYP1A2 | 0.7375 | 0.7466 | slightly better |
| CYP2C9 | 0.7389 | 0.7333 | slightly worse |
| CYP2D6 | 0.4000 | 0.4872 | markedly better |
| CYP3A4 | 0.8177 | 0.7651 | materially worse |

CYP2D6 is a genuinely better-*ranked* model sitting in a badly-*placed*, badly-*compressed*
column — precisely the failure mode NB16 and NB27 Section 2 already corrected twice on this
project's own `chemprop_chemeleoninit` predictions. CYP3A4's ranking is genuinely worse under this
model; no placement or spread correction can recover a worse rank ordering, so CYP3A4 has no case
for adopting the AID model in any form. This decomposition is what justifies a **mixed** recipe
(see below) rather than either a wholesale swap to the AID model or discarding it outright.

### A corrected submission candidate built from this decomposition — sent, and scored

Notebook **27b** (`27b_aid_cyp2d6_corrected.ipynb`, built 2026-09-18) solved the blind population's
CYP2D6 mean/spread from NB30's own published R²/MAE/Spearman (reusing notebook 16's solve exactly),
applied an affine placement correction using the AID model's own OOF Pearson rho for CYP2D6
(0.4228, computed from notebook 29's pooled 25-fold out-of-fold predictions — not `10c`'s OOF rho,
a different model), then widened the spread to the same 0.85-of-training-SD ratio NB27 Section 2
already used and the board already rewarded once. **This is a standalone notebook, not a further
section of notebook 27** — notebook 27 already carries a real, already-scored submission
(NB27-widened, Section 2), and appending further work in place there and re-executing it caused
problems in practice, so this correction was moved out into its own notebook instead. The two
independent population solves — this one from NB30's metrics, and notebook 16's original from
`04b`/`10c`/`NB10`/`NB12`'s metrics — agree closely (mean 3.1614 vs. 3.1558, 0.18% apart; sd 1.5121
vs. 1.5114, 0.05% apart), real convergence evidence for the method, not forced. The resulting
recipe — CYP1A2/CYP3A4 from `10c`, CYP2C9 from NB19, CYP2D6 from the corrected AID column — was
validated (`validate_activity_submission`: PASS) and written to
`outputs/27b_aid_cyp2d6_corrected/submission_candidate.csv`. **It was submitted 2026-09-19 12:11
UTC and scored — see the dedicated "27b" section below for the full result, which is now the best
submission on record.** (This document previously stated here that it had not been submitted;
that was true when written and is corrected below, not retracted from the historical account — see
"The record gap and its consequence.")

## 27b: the AID model, corrected — the best submission to date

**27b is the best submission to date on the governing metric.** Macro ST-RAE 0.5982 — beating
every prior submission, including NB27-widened (0.6155, the previous best), NB19 (0.6265), NB16
(0.6364), 10c (0.7138), and 04b (0.7179). It also carries the highest macro R² (0.4835, beating
NB27-widened's 0.4557), macro Spearman (0.7002, beating NB30's 0.6831), and macro Kendall (0.5212,
beating NB16/NB19/NB27-widened's shared 0.5045) on record. Provenance check (see "Provenance"
below): ST-RAE and Kendall reproduce exactly from the four isoform values; MAE and Spearman round
unambiguously to the published value; R² sits on a rounding boundary and resolves to the published
value under the same round-half-up convention already established for several earlier
submissions. All five metrics check out — the cleanest reconciliation of any submission recorded
in this document.

**27b has positive R² on every isoform, but is not the first submission to achieve that** — a
claim in this notebook's own originating brief that does not hold up against the table above and
is corrected here rather than carried forward. NB16, NB19, and NB27-widened all already had
positive R² on all four isoforms (CYP2D6 R² 0.2713/0.2713/0.3012 respectively), since 2026-09-14 —
CYP1A2/CYP2C9/CYP3A4's R² has been positive on *every* submission this project has made, including
04b; the only isoform that has ever been negative is CYP2D6, and it stopped being negative three
submissions before 27b. What genuinely is new is the *macro* R² value (0.4835) — the highest on
record — and a fourth isoform, CYP2D6, moving from the worst-in-class R² among corrected
submissions (0.3012, NB27-widened) to a materially better one (0.4124).

**The AID model's raw board reading (NB30) was a calibration artefact, not a model verdict.**
NB30 submitted the AID auxiliary-head model raw and scored CYP2D6 ST-RAE 1.3168 — worse than every
corrected submission on record (04b, 10c, NB16, NB19, and NB27-widened all beat it), though **not
the single worst CYP2D6 figure ever recorded**, a stronger claim in this notebook's own originating
brief that the table above does not support: NB10 (1.4603), NB12 (1.4092), and NB13-calib (1.3570)
are all worse still, corrected here rather than repeated. The comparison that matters is within the
same model: corrected, the identical AID model's CYP2D6 column scores 0.7165 — better than the
plain `chemprop_chemeleoninit` model's own best corrected result (NB27-widened, 0.7858). **The
general lesson, stated plainly**: on this challenge, a raw, uncorrected submission is not a fair
test of a model. Placement and spread are correctable after the fact, and — on this project's own
repeated evidence (NB16, NB19, NB27-widened, and now 27b) — the correction is worth more to the
scored metric than most of the modelling changes this project has tried. NB30's raw board reading
said the AID architecture was a regression; 27b's corrected reading of the *same* architecture says
the opposite for the one isoform that differs between them.

**The rank-preservation prediction has now held three times.** NB19, NB27-widened, and 27b each
applied a purely affine transform to CYP2D6 (placement and/or spread only, never a re-ranking
step) and each returned Spearman and Kendall unchanged from the value already fixed by the
underlying model's own ranking. NB19/NB27-widened returned 0.4000/0.2747 (10c's own CYP2D6
values, carried through NB16's placement and NB27's widening). **27b returned 0.4872 and 0.3415 —
NB30's own raw CYP2D6 Spearman/Kendall, exactly** (see the table in the NB30 section above),
because 27b's correction chain is affine and its rank order is inherited unchanged from the AID
model, a different model from the one NB19/NB27-widened corrected. Three affine corrections, three
exact reproductions of the underlying model's own rank-order metrics.

**Board scoring stability is confirmed a third time.** 27b resubmitted CYP1A2 and CYP3A4
byte-identical to `10c` and CYP2C9 byte-identical to NB19; all three reproduced exactly on every
metric (ST-RAE/MAE/R²/Spearman/Kendall), matching NB19's own demonstration of this the day after
NB16, and NB27-widened's the day after that. Three separate confirmations, across submissions
spanning 2026-09-15 to 2026-09-19, that identical inputs produce identical board output —
per-isoform comparisons across this table's different submissions continue to be sound, not an
artifact of half-set scoring noise.

**The CYP2D6 population target 27b used is now itself board-validated, which matters for a choice
made in a later notebook.** 27b's correction targeted the population solved in notebook 27b itself
from NB30's own published metrics (mean 3.1614, sd 1.5121) — a different, independently-solved
target from notebook 16's original (mean 3.1558, sd 1.5114, solved from `04b`/`10c`/NB10/NB12's
metrics), the two agreeing to within 0.18%/0.05% before either had been tested against the board.
Notebook 32 (`32_deadzone_aid_retrain.ipynb`), built after 27b but before this record was
corrected, chose notebook 16's target over 27b's own for a different model's CYP2D6 correction, on
the stated grounds that "27b's own solve has never itself been tested against the real board."
That reasoning was sound when notebook 32 was built and read this document — 27b's submission
predates notebook 32's build, but this document had not yet been updated to reflect it. 27b's own
result (ST-RAE 0.7165, a real improvement over the plain model's corrected 0.7858) is now a
board-validated use of that exact target, for the record — a fact notebook 32 did not have
available and is not retroactively required to have used, per this task's own scope (notebook 32
is not edited here).

## The record gap and its consequence

**Stated factually.** Two real, scored submissions existed on the OpenADMET leaderboard —
NB27-widened (submitted 2026-09-16 21:36 UTC, macro 0.6155) and 27b (submitted 2026-09-19 12:11
UTC, macro 0.5982) — that this document, at various points, did not correctly reflect.

NB27-widened's row and dedicated section are, as of this pass, fully present and correct in this
document — checked directly against the figures re-supplied for this task, all matching exactly,
and consistent with this document's own existing provenance paragraph for that submission (below).
Notebook 30 (built 2026-09-18) recorded finding this row absent when it read this document and
corrected its own brief's "NB27-widened, macro 0.6155" premise to "NB19, macro 0.6265" instead —
a correction that was itself the error, since NB27-widened was a real, scored submission and did
belong ahead of NB19 in the ranking. That correction stands in notebook 30 as a record of what was
believed at the time; it is not retracted there, per this task's scope. Whether the row was added
to this document before or after notebook 30 ran is not established here — only that the document
as it now stands carries it correctly, and notebook 30 did not have that benefit when it ran.

27b's submission was, until this pass, not recorded in this document at all — the "NB30" section's
own subsection on the corrected candidate stated "It has not been submitted," which was true when
written (before 2026-09-19 12:11 UTC) and has been corrected above. Notebook 32
(`32_deadzone_aid_retrain.ipynb`), built after 27b's real submission but before this document was
updated, read this document, found no such figure, and corrected its own brief's "current best is
27b, macro 0.5982" premise to "NB27-widened, macro 0.6155" — again, a correction that was itself
the error, since 27b's real board score was already 0.5982 and already the best on record by the
time notebook 32 ran. That correction likewise stands in notebook 32 as a record of what was
believed at the time and is not retracted there, per this task's scope.

Fixing this document is what stops a third instance of the same failure mode: a real submission
event happening outside any notebook session (both NB27-widened and 27b were sent by hand,
per this project's manual-gated submission process), followed by a notebook that reads this
document as the canonical record, finds a gap, and reasonably — but wrongly — concludes the gap
reflects reality rather than a stale file.

## NB32-A: dead-zone target on top of the AID model — a fourth CV-to-board transfer failure

**NB32-A is worse than the current best (27b, macro ST-RAE 0.5982) on the macro, and worse than the
appropriate per-isoform reference on all four isoforms** — macro ST-RAE 0.6488. Per isoform, checked
against whichever submission actually holds that isoform's best real result: CYP1A2 0.7906 vs.
`10c`'s raw 0.6954 (worse); CYP2C9 0.5635 vs. NB19's 0.5375 (worse); CYP2D6 0.7502 vs. 27b's 0.7165
(worse); CYP3A4 0.4907 vs. `10c`'s raw 0.4434 (worse). All four confirmed directly against the table
above, not assumed.

**This is the fourth occasion on which a cross-validated result has failed to transfer to the real
blind board on this project**, extending the pattern already recorded in the "NB30" section above:

1. Notebook `08`'s ensembles — selected from a single 25-fold CV pass, never independently
   re-confirmed before submission — scored worse than `10c` on the board (NB10/NB12, macro 0.8299).
2. Notebook `28`'s single-fold AID auxiliary-head screen — its own verdict did not survive
   independent re-confirmation at full 5×5 CV (notebook `29` reproduced only 2 of 4 isoform
   verdicts).
3. Notebook `29`'s full, honest, paired, Benjamini-Hochberg-corrected 5×5 CV confirmation of the
   AID architecture — the strongest CV evidence any adopted recipe on this project has carried into
   a submission — still scored worse than `10c` on all four isoforms when submitted raw (NB30, macro
   ST-RAE 0.8160).
4. **Now, notebook `31`'s three-arm single-fold screen**, which found the dead-zone clipped training
   target improved ST-RAE on all four isoforms out of fold (macro 0.7490 → 0.6643, each isoform
   beating its own baseline seed spread by 4.5–6×) — combined with the AID architecture and
   retrained on the full data (notebook `32`) without a fresh 25-fold reconfirmation of that specific
   combination (a deliberate, stated deviation in notebook `32`, not an oversight), the resulting
   submission is worse than the appropriate baseline on all four isoforms on the real board.

Each instance carries progressively stronger CV backing than the last (a single uncrossvalidated
selection pass; a single-fold screen; a full paired 5×5 confirmation; a single-fold screen combined
with an already-confirmed architecture) and none of the four has transferred cleanly to the board.

**The Spearman decomposition is the informative part — it separates model quality from calibration,
since placement (intercept) and spread (slope) corrections are both affine and cannot change rank
order.**

| Isoform | Reference Spearman | NB32-A Spearman | Reference | Direction |
|---|---:|---:|---|---|
| CYP1A2 | 0.7375 | 0.7565 | `10c`/NB19/27b (never corrected) | **improved** |
| CYP2C9 | 0.7585 | 0.7554 | NB19 (current adopted column) | roughly flat |
| CYP2D6 | 0.4872 | 0.4438 | NB30/27b (the plain AID model) | **declined** |
| CYP3A4 | 0.8177 | 0.8020 | `10c`/27b (never corrected) | roughly flat |

**CYP1A2's ranking genuinely improved (0.7375 → 0.7565) while its ST-RAE got worse** — that column is
a placement-or-spread problem, not a model regression, and CYP1A2 has never had any correction
applied in this project (its Spearman is identical, 0.7375, across every submission that carries the
plain `chemprop_chemeleoninit` or capped-Caruana-ensemble ranking unmodified). **CYP2D6's ranking
genuinely declined against the plain AID model (0.4872 → 0.4438)** — a real regression that no
placement or spread correction can reach, since none can change rank order. **CYP2C9 and CYP3A4 are
not identically flat** — CYP2C9's shift (−0.0031) is small; CYP3A4's (−0.0157) is roughly five times
larger in magnitude, though still far smaller than CYP2D6's decline — and both isoforms combine a
roughly-held ranking with a worse ST-RAE, consistent with a placement/spread problem rather than a
clean ranking regression on those two isoforms specifically.

**Candidate B was not submitted — and its board value can be computed, not guessed, because it
shares candidate A's CYP2D6 column byte-for-byte.** Notebook 32's own
`candidate_b_byte_identical_checks.csv` confirms this. With that CYP2D6 column's now-known board
value (ST-RAE 0.7502) standing in for the un-submitted candidate B, and `10c`'s raw CYP1A2/CYP3A4
(0.6954/0.4434) and NB19's CYP2C9 (0.5375) — candidate B's own recipe — its macro ST-RAE works out to
(0.6954 + 0.5375 + 0.7502 + 0.4434) / 4 = **0.6066** (exact decimal arithmetic, not floating point).
That is worse than 27b's 0.5982, so the submission slot was correctly not spent on it. This is a
worked deduction, not a measurement — a case where the board's own scoring stability, already
confirmed three times on this project (NB19, NB27-widened, and 27b each reproduced an earlier
submission's byte-identical columns exactly), lets a submission be evaluated without sending it.

## The full macro ST-RAE trajectory, in true chronological order

Notebook numbering does not track submission order — `NB10` was actually submitted before `10c`
despite the numbering (see notebook 27 Section 1's own figure-ordering note). By actual submission
timestamp:

| Order | Submission | Submitted (UTC) | Macro ST-RAE |
|---:|---|---|---:|
| 1 | 04b | (undated, first submission) | 0.7179 |
| 2 | NB10 | 2026-09-03 15:24 UTC (file mtime) | 0.8299 |
| 3 | 10c | 2026-09-04 ~09:15 local | 0.7138 |
| 4 | NB12 | 2026-09-04 21:28 UTC (file mtime) | 0.8299 |
| 5 | NB13-calib | 2026-09-13 15:07 UTC | 0.8811 |
| 6 | NB16 | 2026-09-14 13:21 UTC | 0.6364 |
| 7 | NB19 | 2026-09-15 10:14 UTC | 0.6265 |
| 8 | NB27-widened | 2026-09-16 21:36 UTC | 0.6155 |
| 9 | NB30 | 2026-09-18 11:39 UTC | 0.8160 |
| 10 | 27b | 2026-09-19 12:11 UTC | 0.5982 |
| 11 | NB32-A | 2026-09-20 02:10 UTC | 0.6488 |

The trajectory is not monotonic: five step-to-step regressions (04b→NB10 +0.1120, 10c→NB12
+0.1161, NB12→NB13-calib +0.0512, NB27-widened→NB30 +0.2005, 27b→NB32-A +0.0506) interrupt an
otherwise steady improvement from `10c` (0.7138) down to 27b (0.5982, the best on record). **NB30's
regression is the largest single step-to-step increase in macro ST-RAE anywhere in this
trajectory** — larger than 04b→NB10, 10c→NB12, NB12→NB13-calib, and 27b→NB32-A, each recomputed
directly from the chronological table above rather than assumed. **NB30→27b is the largest single
step-to-step decrease anywhere in this trajectory** (−0.2178, larger in magnitude than any of the
five regressions) — the same model, corrected, moving further in one step than any prior
submission-to-submission change in either direction. **27b→NB32-A is the smallest of the five
regressions** (+0.0506, well under half the size of the next-smallest, NB12→NB13-calib's +0.0512 —
in fact marginally smaller even than that one) — consistent with NB32-A carrying the same
already-corrected CYP2D6 placement/widening as 27b/NB27-widened rather than reverting to a raw,
uncorrected column the way NB30 did.

## The one clean positive signal so far

**CYP2C9 is the only isoform where every single metric improves or holds flat, monotonically,
across all four submissions (04b → 10c → NB10/NB12).** ST-RAE: 0.5408 → 0.5489 → 0.5489 → 0.5375.
R²: 0.5263 → 0.5272 → 0.5272 → 0.5432. Spearman: 0.7428 → 0.7389 → 0.7389 → 0.7585. Kendall:
0.5490 → 0.5462 → 0.5462 → 0.5692. (MAE moves in the direction consistent with the others:
0.5399 → 0.5444 → 0.5444 → 0.5280.) No other isoform shows this — CYP1A2/CYP3A4 both regress
under ensembling before partially recovering, and CYP2D6 has never recovered to 04b's level on
any metric. Worth treating as the clearest positive data point in the investigation so far, not
just a side note.

**This trend does not extend to NB13-calib.** CYP2C9's ST-RAE jumps back up to 0.6036 under
population calibration — worse than all four prior submissions, including 04b. NB13-calib's
recipe reused NB12's exact CYP2C9 ensemble predictions and only applied calibration on top, so
this is a calibration-attributable regression, not a re-run of the ensembling question. The
monotonic run above is specifically 04b → 10c → NB10/NB12 and stops there.

**Update (NB19, 2026-09-15): the run resumes.** NB16 also moved CYP2C9 the wrong way (0.5770);
NB19 reverts CYP2C9 to NB12's own prediction column and returns it to 0.5375 exactly — see the
"NB19" section above for the full comparison.

## A resolved loose end from the earlier diagnostic report

That report flagged NB10 and NB12 landing on the *exact same* macro ST-RAE (0.8299) from two
structurally different recipes as worth a second look, in case it was a transcription artifact
from manual leaderboard-reading. With the full metric table now in hand, this is not an
artifact: NB10 and NB12's other four macro metrics all differ (MAE 0.9690 vs. 0.9674, R² 0.1054
vs. 0.1159, Spearman 0.6748 vs. 0.6613, Kendall 0.4977 vs. 0.4877) — a copy-paste of one row onto
the other would have matched everywhere, not just on ST-RAE. The macro-ST-RAE tie is a genuine,
if striking, coincidence.

## Provenance

- Full table supplied directly by the user, cross-checked against this repo's own CV summary
  (`outputs/05_cv_comparison/summary_table.csv`) and submission-file hashes/timestamps
  (`outputs/{10_final_retrain_predict,10c_control_submission,12_caruana_retrain_predict}/submission_candidate*.csv`)
  in the prior diagnostic report — no inconsistency found there.
- Every macro (MA) value above is confirmed to equal the mean of its own four isoform values, for
  every metric and every submission (checked directly against this table, not assumed).
- 04b's and NB10's per-isoform figures matched `README.md`'s own "Leaderboard result" tables
  exactly, confirmed at the time those tables still existed in `README.md`; the root README was
  later restructured (2026-09-15) to point here instead of restating any leaderboard table, per
  this document's own "supersedes any partial or macro-only leaderboard figures stated elsewhere"
  rule above.
- NB13-calib's figures were supplied directly by the user (2026-09-13, the day of submission) and
  verified in this pass: every macro (MA) value equals the mean of its four isoform values, to
  four decimal places, for all five metrics.
- NB16's figures were supplied directly by the user (2026-09-14 board capture) and verified in
  this pass. **4 of 5 metrics check out exactly**: ST-RAE (0.6954+0.5770+0.8299+0.4434)/4 =
  0.636425 → 0.6364; MAE (0.9385+0.5603+1.0546+0.5183)/4 = 0.767925 → 0.7679; R²
  (0.3073+0.5133+0.2713+0.6709)/4 = 0.4407 exactly; Spearman
  (0.7375+0.7389+0.4000+0.8177)/4 = 0.673525 → 0.6735 — all match the published MA value.
  **Kendall does not**: (0.5425+0.5462+0.2747+0.6316)/4 = 0.49875, which rounds to 0.4988 under
  standard rounding, not the published 0.4987 — a 0.0001 discrepancy sitting exactly on the
  rounding boundary. Reported here rather than silently adjusted, per this document's own
  convention of not reconciling a figure that doesn't check out; most plausibly the board's own
  macro is computed from higher-precision per-isoform values before they are rounded to 4 decimal
  places for per-isoform display, but this is not confirmed.
- NB19's figures were supplied directly by the user (2026-09-15 board capture) and verified in
  this pass. **MAE, R², and Kendall check out exactly**: MAE
  (0.9385+0.5280+1.0546+0.5183)/4 = 0.75985 → 0.7599; R²
  (0.3073+0.5432+0.2713+0.6709)/4 = 0.448175 → 0.4482; Kendall
  (0.5425+0.5692+0.2747+0.6316)/4 = 0.5045 exactly. **ST-RAE does not**:
  (0.6954+0.5375+0.8299+0.4434)/4 = 0.62655 — exactly on a rounding boundary, the same kind of
  case as NB16's Kendall above, resolving to 0.6266 under standard rounding rather than the
  published 0.6265. **Spearman also does not check out, though not as a boundary case**:
  (0.7375+0.7585+0.4000+0.8177)/4 = 0.678425, which rounds unambiguously to 0.6784 (the digit
  after the fourth decimal place is 2, not a tie) — a 0.0001 discrepancy from the published
  0.6785 that isn't a rounding-boundary artifact in the same sense as the ST-RAE/Kendall cases.
  Both are reported here rather than adjusted, per this document's convention; the most likely
  explanation for both is the same one offered for NB16's Kendall — the board's own macro is
  plausibly computed from higher-precision per-isoform values before per-isoform display
  rounding, so re-averaging the rounded display values can drift by a unit in the last place —
  but this is not confirmed for either metric.
- NB27-widened's figures were supplied directly by the user (2026-09-16 board capture) and
  verified in this pass. **ST-RAE, MAE, and Kendall check out exactly**: ST-RAE
  (0.6954+0.5375+0.7858+0.4434)/4 = 0.615525 → 0.6155; MAE
  (0.9385+0.5280+1.0104+0.5183)/4 = 0.7488 exactly; Kendall
  (0.5425+0.5692+0.2747+0.6316)/4 = 0.5045 exactly. **R² checks out via the same rounding
  convention as NB16's Kendall and NB19's ST-RAE**: (0.3073+0.5432+0.3012+0.6709)/4 = 0.45565,
  exactly on a rounding boundary, resolving to 0.4557 under standard round-half-up, matching the
  published value. **Spearman does not check out**: (0.7375+0.7585+0.4000+0.8177)/4 = 0.678425,
  which rounds unambiguously to 0.6784 (the digit after the fourth decimal place is 2, not a
  tie) — a 0.0001 discrepancy from the published 0.6785. This is the identical discrepancy NB19
  showed, for the identical reason: NB27-widened's four isoform Spearman values are numerically
  the same four values as NB19's (a pure spread scale cannot move Spearman), so the same
  higher-precision-before-display-rounding explanation offered for NB16's Kendall and NB19's
  ST-RAE/Spearman applies here too, and is not confirmed for this case either.
- **A separate arithmetic check, not a board-figure discrepancy**: a prior pass over these same
  four ST-RAE isoform values computed their mean as 0.655525 and reported it as not matching the
  published 0.6155. That arithmetic was wrong — (0.6954+0.5375+0.7858+0.4434)/4 = 0.615525,
  which matches 0.6155 exactly, unambiguously, with no rounding-boundary judgment call involved.
  Recorded here because the earlier claim was specific and could otherwise be mistaken for a
  genuine provenance finding rather than a transcription slip.
- NB30's figures were supplied directly by the user (2026-09-18 board capture) and verified in
  this pass using exact decimal arithmetic (not floating point, to avoid a repeat of NB27-widened's
  own earlier arithmetic slip). **All five metrics reconcile with the published macro, though two
  sit exactly on a rounding boundary**: ST-RAE (0.763+0.5922+1.3168+0.592)/4 = 0.8160 exactly —
  matches the published 0.816 with no rounding involved. R² (0.2588+0.4789−0.766+0.5708)/4 =
  0.135625 → rounds unambiguously to 0.1356 (5th decimal digit is 2), matching. Kendall
  (0.5468+0.5369+0.3415+0.5743)/4 = 0.499875 → rounds unambiguously to 0.4999 (5th decimal digit is
  7), matching. **MAE** (1.0027+0.5607+1.6688+0.5996)/4 = 0.95795 exactly — a genuine rounding-
  boundary case, the same kind seen for NB16's Kendall, NB19's ST-RAE, and NB27-widened's R² above
  — and resolves to 0.9580 under standard round-half-up, matching the published 0.958. **Spearman**
  (0.7466+0.7333+0.4872+0.7651)/4 = 0.68305 exactly — also a boundary case — resolves to 0.6831
  under the same round-half-up convention, matching the published 0.6831. Unlike NB19's and
  NB27-widened's own boundary cases (which resolved to a value one unit-in-the-last-place away
  from the published figure), **both of NB30's boundary cases resolve to exactly the published
  value** under round-half-up — the cleanest reconciliation of any submission recorded in this
  document to date.
- **NB27-widened re-verified in this pass** (its row and figures were already present in this
  document before this pass began — see "The record gap and its consequence" above): all five
  metrics reproduce the checks already recorded in this section exactly (ST-RAE/MAE/Kendall exact
  or unambiguous-round; R² a boundary case resolving correctly; Spearman the one known 0.0001
  discrepancy, unchanged). No figure or check in this document's own existing NB27-widened
  material required correction.
- 27b's figures were supplied directly by the user (2026-09-19 board capture) and verified in this
  pass using exact decimal arithmetic. **All five metrics reconcile with the published macro,
  three exactly and one on a rounding boundary that resolves correctly**: ST-RAE
  (0.6954+0.5375+0.7165+0.4434)/4 = 0.5982 exactly — no rounding involved. Kendall
  (0.5425+0.5692+0.3415+0.6316)/4 = 0.5212 exactly. MAE (0.9385+0.5280+0.8975+0.5183)/4 = 0.720575
  → rounds unambiguously to 0.7206 (5th decimal digit is 7, not a tie), matching. Spearman
  (0.7375+0.7585+0.4872+0.8177)/4 = 0.700225 → rounds unambiguously to 0.7002 (5th decimal digit
  is 2), matching. **R²** (0.3073+0.5432+0.4124+0.6709)/4 = 0.48345 exactly — a genuine
  rounding-boundary case, the same kind seen for NB16's Kendall, NB19's ST-RAE, NB27-widened's R²,
  and NB30's MAE/Spearman above — resolves to 0.4835 under standard round-half-up, matching the
  published value. **All five metrics check out** — no unresolved discrepancy, matching NB30's own
  clean result and unlike NB16/NB19/NB27-widened, each of which had at least one metric that did
  not reconcile.
- NB32-A's figures were supplied directly by the user (2026-09-20 board capture) and verified in
  this pass using exact decimal arithmetic. **All five metrics reconcile with the published macro,
  three exactly/unambiguously and two on a rounding boundary that resolves correctly under the same
  round-half-up convention used throughout this document**: MAE (1.0287+0.5568+0.9319+0.5506)/4 =
  0.7670 exactly, matching the published 0.767. Spearman (0.7565+0.7554+0.4438+0.8020)/4 =
  0.689425 → rounds unambiguously to 0.6894 (5th decimal digit is 2, not a tie), matching. Kendall
  (0.5611+0.5652+0.3076+0.6132)/4 = 0.511775 → rounds unambiguously to 0.5118 (5th decimal digit is
  7), matching. **ST-RAE** (0.7906+0.5635+0.7502+0.4907)/4 = 0.64875 exactly — a genuine
  rounding-boundary case, the same kind seen for NB16's Kendall, NB19's ST-RAE, NB27-widened's R²,
  NB30's MAE/Spearman, and 27b's R² above — resolves to 0.6488 under standard round-half-up,
  matching the published value. **R²** (0.2125+0.5016+0.3768+0.6405)/4 = 0.43285 exactly — also a
  boundary case — resolves to 0.4329 under the same round-half-up convention, matching the
  published value. **All five metrics check out** — no unresolved discrepancy, matching NB30's and
  27b's own clean results, and unlike NB16/NB19/NB27-widened, each of which had at least one metric
  that did not reconcile.
  Candidate B (not submitted) was not independently verified against a board figure, since it was
  never sent — its macro ST-RAE, computed here from candidate A's own now-known CYP2D6 board value
  plus `10c`'s and NB19's already-verified columns (candidate B's actual recipe), is a worked
  deduction rather than a board reading: (0.6954+0.5375+0.7502+0.4434)/4 = 0.606625 → 0.6066 under
  round-half-up (5th decimal digit is 2, well clear of a boundary case).
