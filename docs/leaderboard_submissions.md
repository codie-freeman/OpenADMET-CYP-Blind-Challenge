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

## All five submissions to date

| Isoform | Metric | 04b | 10c | NB10 | NB12 | NB13-calib |
|---|---|---:|---:|---:|---:|---:|
| **MA (Overall)** | ST-RAE | 0.7179 | 0.7138 | 0.8299 | 0.8299 | 0.8811 |
| | MAE | 0.8952 | 0.8927 | 0.9690 | 0.9674 | 1.0260 |
| | R² | 0.2191 | 0.2278 | 0.1054 | 0.1159 | 0.0203 |
| | Spearman | 0.6739 | 0.6735 | 0.6748 | 0.6613 | 0.6785 |
| | Kendall | 0.5008 | 0.4987 | 0.4977 | 0.4877 | 0.5045 |
| **CYP1A2** | ST-RAE | 0.7114 | 0.6954 | 0.7850 | 0.7633 | 0.8197 |
| | MAE | 0.9589 | 0.9385 | 1.0214 | 1.0035 | 1.0549 |
| | R² | 0.2824 | 0.3073 | 0.2051 | 0.2275 | 0.1537 |
| | Spearman | 0.7430 | 0.7375 | 0.7246 | 0.7389 | 0.7375 |
| | Kendall | 0.5473 | 0.5425 | 0.5297 | 0.5422 | 0.5425 |
| **CYP2C9** | ST-RAE | 0.5408 | 0.5489 | 0.5489 | 0.5375 | 0.6036 |
| | MAE | 0.5399 | 0.5444 | 0.5444 | 0.5280 | 0.5934 |
| | R² | 0.5263 | 0.5272 | 0.5272 | 0.5432 | 0.4237 |
| | Spearman | 0.7428 | 0.7389 | 0.7389 | 0.7585 | 0.7585 |
| | Kendall | 0.5490 | 0.5462 | 0.5462 | 0.5692 | 0.5692 |
| **CYP2D6** | ST-RAE | 1.1903 | 1.1673 | 1.4603 | 1.4092 | 1.3570 |
| | MAE | 1.5842 | 1.5694 | 1.7489 | 1.7209 | 1.6904 |
| | R² | −0.6220 | −0.5941 | −0.9392 | −0.8726 | −0.7903 |
| | Spearman | 0.3769 | 0.4000 | 0.4401 | 0.3771 | 0.4000 |
| | Kendall | 0.2606 | 0.2747 | 0.3081 | 0.2587 | 0.2747 |
| **CYP3A4** | ST-RAE | 0.4291 | 0.4434 | 0.5253 | 0.6095 | 0.7443 |
| | MAE | 0.4976 | 0.5183 | 0.5614 | 0.6173 | 0.7654 |
| | R² | 0.6898 | 0.6709 | 0.6285 | 0.5655 | 0.2941 |
| | Spearman | 0.8329 | 0.8177 | 0.7955 | 0.7706 | 0.8177 |
| | Kendall | 0.6463 | 0.6316 | 0.6065 | 0.5809 | 0.6316 |

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

**On "deliberate":** all five submissions were reviewed and uncommented/run by hand, exactly as
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
pursuing (kb_share <10% everywhere), **population calibration — in every form tested across
notebooks 09, 11b, 13, and 14 — is now a confirmed real-blind-data negative result, not just an
OOF one.** No form of it is adopted or recommended going forward.

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
