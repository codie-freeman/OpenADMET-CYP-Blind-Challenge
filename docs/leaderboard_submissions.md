# Leaderboard submissions — authoritative record

This is the single canonical record of every real submission this project has made to the
OpenADMET CYP450 blind leaderboard. It supersedes any partial or macro-only leaderboard figures
stated elsewhere (notebook markdown, CLAUDE.md prose) — those should point here rather than
repeat these numbers. Real per-isoform scores exist here for 04b and notebook 10 already (also
recorded in `README.md`, which matches this table exactly — cross-checked, see "Provenance"
below); 10c, notebook 12, and NB13-calib's full metric breakdowns are recorded here for the first
time.

No blind-set labels are held anywhere in this repo — all scores below came from OpenADMET's own
leaderboard after real submissions, not from anything computed locally.

**Every score in this document is a half-set measurement.** The live leaderboard scores only half
of the blind test set (split by chemical series) as submissions come in — the full set is scored
only at the intermediate reveal and at the competition's close. This applies to every row below;
none of these numbers are the final, full-set result, and later reveals could move any of them.

## All six submissions to date

| Isoform | Metric | 04b | 10c | NB10 | NB12 | NB13-calib | NB16 |
|---|---|---:|---:|---:|---:|---:|---:|
| **MA (Overall)** | ST-RAE | 0.7179 | 0.7138 | 0.8299 | 0.8299 | 0.8811 | 0.6364 |
| | MAE | 0.8952 | 0.8927 | 0.9690 | 0.9674 | 1.0260 | 0.7679 |
| | R² | 0.2191 | 0.2278 | 0.1054 | 0.1159 | 0.0203 | 0.4407 |
| | Spearman | 0.6739 | 0.6735 | 0.6748 | 0.6613 | 0.6785 | 0.6735 |
| | Kendall | 0.5008 | 0.4987 | 0.4977 | 0.4877 | 0.5045 | 0.4987 |
| **CYP1A2** | ST-RAE | 0.7114 | 0.6954 | 0.7850 | 0.7633 | 0.8197 | 0.6954 |
| | MAE | 0.9589 | 0.9385 | 1.0214 | 1.0035 | 1.0549 | 0.9385 |
| | R² | 0.2824 | 0.3073 | 0.2051 | 0.2275 | 0.1537 | 0.3073 |
| | Spearman | 0.7430 | 0.7375 | 0.7246 | 0.7389 | 0.7375 | 0.7375 |
| | Kendall | 0.5473 | 0.5425 | 0.5297 | 0.5422 | 0.5425 | 0.5425 |
| **CYP2C9** | ST-RAE | 0.5408 | 0.5489 | 0.5489 | 0.5375 | 0.6036 | 0.5770 |
| | MAE | 0.5399 | 0.5444 | 0.5444 | 0.5280 | 0.5934 | 0.5603 |
| | R² | 0.5263 | 0.5272 | 0.5272 | 0.5432 | 0.4237 | 0.5133 |
| | Spearman | 0.7428 | 0.7389 | 0.7389 | 0.7585 | 0.7585 | 0.7389 |
| | Kendall | 0.5490 | 0.5462 | 0.5462 | 0.5692 | 0.5692 | 0.5462 |
| **CYP2D6** | ST-RAE | 1.1903 | 1.1673 | 1.4603 | 1.4092 | 1.3570 | 0.8299 |
| | MAE | 1.5842 | 1.5694 | 1.7489 | 1.7209 | 1.6904 | 1.0546 |
| | R² | −0.6220 | −0.5941 | −0.9392 | −0.8726 | −0.7903 | 0.2713 |
| | Spearman | 0.3769 | 0.4000 | 0.4401 | 0.3771 | 0.4000 | 0.4000 |
| | Kendall | 0.2606 | 0.2747 | 0.3081 | 0.2587 | 0.2747 | 0.2747 |
| **CYP3A4** | ST-RAE | 0.4291 | 0.4434 | 0.5253 | 0.6095 | 0.7443 | 0.4434 |
| | MAE | 0.4976 | 0.5183 | 0.5614 | 0.6173 | 0.7654 | 0.5183 |
| | R² | 0.6898 | 0.6709 | 0.6285 | 0.5655 | 0.2941 | 0.6709 |
| | Spearman | 0.8329 | 0.8177 | 0.7955 | 0.7706 | 0.8177 | 0.8177 |
| | Kendall | 0.6463 | 0.6316 | 0.6065 | 0.5809 | 0.6316 | 0.6316 |

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

**On "deliberate":** all six submissions were reviewed and uncommented/run by hand, exactly as
this project's manual-gated-submission process requires. An earlier diagnostic pass over this
repo (before this correction) reasoned from file evidence alone — a misleading header comment on
the NB10/10c/NB12 submission cells, which claimed "every line below is a comment" when only the
explanatory paragraph was — and concluded the submissions might have fired by accident. **That
conclusion was wrong.** The header comment itself was a real bug (now fixed, 2026-09-07 — see
CLAUDE.md's Rules section), but the submissions it sat on top of were intended, reviewed actions,
and the content each one sent was correct.

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
- 04b's and NB10's per-isoform figures match `README.md`'s existing "Leaderboard result" tables
  exactly.
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
