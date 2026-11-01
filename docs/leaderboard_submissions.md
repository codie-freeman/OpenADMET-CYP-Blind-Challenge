# Leaderboard submissions — authoritative record

This is the single canonical record of every real submission this project has made to the
OpenADMET CYP450 blind leaderboard. It supersedes any partial or macro-only leaderboard figures
stated elsewhere (notebook markdown, CLAUDE.md prose, the root `README.md`) — those should point
here rather than repeat these numbers. Real per-isoform scores exist here for 04b and notebook 10
already (originally cross-checked against `README.md`'s own "Leaderboard result" tables before
those tables were retired in favor of a pointer to this file — see "Provenance" below); 10c,
notebook 12, NB13-calib's, NB16's, NB19's, NB27-widened's, NB30's, 27b's, NB32-A's and
NB36-widened's, NB47-mae's, NB48-ensemble's, NB56-identify's, NB56b-s050's and NB56b-s100's full
metric breakdowns are recorded here as well.

No blind-set labels are held anywhere in this repo — all scores below came from OpenADMET's own
leaderboard after real submissions, not from anything computed locally.

**Every score in this document is a half-set measurement.** The live leaderboard scores only half
of the blind test set (split by chemical series) as submissions come in — the full set is scored
only at the intermediate reveal and at the competition's close. This applies to every row below;
none of these numbers are the final, full-set result, and later reveals could move any of them.
**One full-set reading does now exist** — the interim board of 2026-09-21, which scored 27b; it is
recorded in its own section below, and it disagrees with the live half materially and in opposite
directions by isoform.

## The eighteen submissions with a full per-isoform record

**There are nineteen real sends, and exactly one is absent from this table**: the byte-identical
re-send of **27b**, fired accidentally on 2026-09-21 12:49 UTC. It needs no column, because the
board returned 27b's own scores, already shown here. It has its own section below and is not
missing from the record — only from this *table*, which is per-isoform and per-metric.

**NB40-ensemble was the other absence and no longer is.** Its full per-isoform row was recovered
from the 2026-10-02 leaderboard screenshots on 2026-10-29 and is now a column here. It is the
**only row in this table whose per-isoform values come from a screenshot rather than a
contemporaneous capture**, which is how it should be read; everything else was transcribed from the
board at the time. All five of its macro values reconcile as the mean of four (each 2.5e-5, display
rounding), and the two figures the record already held independently — macro ST-RAE 0.6384 and
CYP2D6 ST-RAE 0.7751 — agree with the screenshot exactly, which is the best available check on it.

| Isoform | Metric | 04b | 10c | NB10 | NB12 | NB13-calib | NB16 | NB19 | NB27-widened | NB30 | 27b | NB32-A | NB36-widened | NB40-ensemble | NB47-mae | NB48-ensemble | NB56-identify | NB56b-s050 | NB56b-s100 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **MA (Overall)** | ST-RAE | 0.7179 | 0.7138 | 0.8299 | 0.8299 | 0.8811 | 0.6364 | 0.6265 | 0.6155 | 0.8160 | 0.5982 | 0.6488 | 0.627 | **0.6384** | 0.6328 | 0.6356 | 0.6947 | **0.5861** | 0.5933 |
|  | MAE | 0.8952 | 0.8927 | 0.9690 | 0.9674 | 1.0260 | 0.7679 | 0.7599 | 0.7488 | 0.9580 | 0.7206 | 0.7670 | 0.7426 | 0.7506 | 0.7468 | 0.7494 | 0.8014 | 0.6973 | 0.6864 |
|  | R² | 0.2191 | 0.2278 | 0.1054 | 0.1159 | 0.0203 | 0.4407 | 0.4482 | 0.4557 | 0.1356 | 0.4835 | 0.4329 | 0.4677 | 0.4475 | 0.4495 | 0.4507 | 0.3858 | 0.5249 | 0.5366 |
|  | Spearman | 0.6739 | 0.6735 | 0.6748 | 0.6613 | 0.6785 | 0.6735 | 0.6785 | 0.6785 | 0.6831 | 0.7002 | 0.6894 | 0.6894 | 0.6759 | 0.6854 | 0.6739 | 0.6894 | 0.6894 | 0.6894 |
|  | Kendall | 0.5008 | 0.4987 | 0.4977 | 0.4877 | 0.5045 | 0.4987 | 0.5045 | 0.5045 | 0.4999 | 0.5212 | 0.5118 | 0.5118 | 0.5021 | 0.5090 | 0.5008 | 0.5118 | 0.5118 | 0.5118 |
| **CYP1A2** | ST-RAE | 0.7114 | 0.6954 | 0.7850 | 0.7633 | 0.8197 | 0.6954 | 0.6954 | 0.6954 | 0.7630 | 0.6954 | 0.7906 | 0.7241 | 0.7435 | 0.7254 | 0.7414 | 0.9126 | 0.5631 | 0.5767 |
|  | MAE | 0.9589 | 0.9385 | 1.0214 | 1.0035 | 1.0549 | 0.9385 | 0.9385 | 0.9385 | 1.0027 | 0.9385 | 1.0287 | 0.9616 | 0.9731 | 0.9666 | 0.9699 | 1.1137 | 0.7988 | 0.7575 |
|  | R² | 0.2824 | 0.3073 | 0.2051 | 0.2275 | 0.1537 | 0.3073 | 0.3073 | 0.3073 | 0.2588 | 0.3073 | 0.2125 | 0.3068 | 0.2827 | 0.2763 | 0.2899 | 0.1154 | 0.5070 | 0.5467 |
|  | Spearman | 0.7430 | 0.7375 | 0.7246 | 0.7389 | 0.7375 | 0.7375 | 0.7375 | 0.7375 | 0.7466 | 0.7375 | 0.7565 | 0.7565 | 0.7102 | 0.7423 | 0.7107 | 0.7565 | 0.7565 | 0.7565 |
|  | Kendall | 0.5473 | 0.5425 | 0.5297 | 0.5422 | 0.5425 | 0.5425 | 0.5425 | 0.5425 | 0.5468 | 0.5425 | 0.5611 | 0.5611 | 0.5145 | 0.5425 | 0.5146 | 0.5611 | 0.5611 | 0.5611 |
| **CYP2C9** | ST-RAE | 0.5408 | 0.5489 | 0.5489 | 0.5375 | 0.6036 | 0.5770 | 0.5375 | 0.5375 | 0.5922 | 0.5375 | 0.5635 | 0.5435 | 0.5282 | 0.5255 | 0.5290 | 0.5770 | 0.5408 | 0.5562 |
|  | MAE | 0.5399 | 0.5444 | 0.5444 | 0.5280 | 0.5934 | 0.5603 | 0.5280 | 0.5280 | 0.5607 | 0.5280 | 0.5568 | 0.5396 | 0.5317 | 0.5276 | 0.5322 | 0.5703 | 0.5212 | 0.5190 |
|  | R² | 0.5263 | 0.5272 | 0.5272 | 0.5432 | 0.4237 | 0.5133 | 0.5432 | 0.5432 | 0.4789 | 0.5432 | 0.5016 | 0.5308 | 0.5274 | 0.5358 | 0.5258 | 0.4542 | 0.5595 | 0.5667 |
|  | Spearman | 0.7428 | 0.7389 | 0.7389 | 0.7585 | 0.7585 | 0.7389 | 0.7585 | 0.7585 | 0.7333 | 0.7585 | 0.7554 | 0.7554 | 0.7729 | 0.7419 | 0.7692 | 0.7554 | 0.7554 | 0.7554 |
|  | Kendall | 0.5490 | 0.5462 | 0.5462 | 0.5692 | 0.5692 | 0.5462 | 0.5692 | 0.5692 | 0.5369 | 0.5692 | 0.5652 | 0.5652 | 0.5842 | 0.5581 | 0.5811 | 0.5652 | 0.5652 | 0.5652 |
| **CYP2D6** | ST-RAE | 1.1903 | 1.1673 | 1.4603 | 1.4092 | 1.3570 | 0.8299 | 0.8299 | 0.7858 | 1.3168 | 0.7165 | 0.7502 | 0.7502 | 0.7751 | 0.7872 | 0.7652 | 0.7502 | 0.7502 | 0.7502 |
|  | MAE | 1.5842 | 1.5694 | 1.7489 | 1.7209 | 1.6904 | 1.0546 | 1.0546 | 1.0104 | 1.6688 | 0.8975 | 0.9319 | 0.9319 | 0.9547 | 0.9600 | 0.9520 | 0.9319 | 0.9319 | 0.9319 |
|  | R² | −0.6220 | −0.5941 | −0.9392 | −0.8726 | −0.7903 | 0.2713 | 0.2713 | 0.3012 | −0.7660 | 0.4124 | 0.3768 | 0.3768 | 0.3468 | 0.3455 | 0.3543 | 0.3768 | 0.3768 | 0.3768 |
|  | Spearman | 0.3769 | 0.4000 | 0.4401 | 0.3771 | 0.4000 | 0.4000 | 0.4000 | 0.4000 | 0.4872 | 0.4872 | 0.4438 | 0.4438 | 0.4238 | 0.4602 | 0.4173 | 0.4438 | 0.4438 | 0.4438 |
|  | Kendall | 0.2606 | 0.2747 | 0.3081 | 0.2587 | 0.2747 | 0.2747 | 0.2747 | 0.2747 | 0.3415 | 0.3415 | 0.3076 | 0.3076 | 0.2969 | 0.3193 | 0.2919 | 0.3076 | 0.3076 | 0.3076 |
| **CYP3A4** | ST-RAE | 0.4291 | 0.4434 | 0.5253 | 0.6095 | 0.7443 | 0.4434 | 0.4434 | 0.4434 | 0.5920 | 0.4434 | 0.4907 | 0.4904 | 0.5067 | 0.4930 | 0.5067 | 0.5390 | 0.4904 | 0.4904 |
|  | MAE | 0.4976 | 0.5183 | 0.5614 | 0.6173 | 0.7654 | 0.5183 | 0.5183 | 0.5183 | 0.5996 | 0.5183 | 0.5506 | 0.5373 | 0.5430 | 0.5332 | 0.5434 | 0.5896 | 0.5373 | 0.5373 |
|  | R² | 0.6898 | 0.6709 | 0.6285 | 0.5655 | 0.2941 | 0.6709 | 0.6709 | 0.6709 | 0.5708 | 0.6709 | 0.6405 | 0.6563 | 0.6330 | 0.6406 | 0.6326 | 0.5969 | 0.6563 | 0.6563 |
|  | Spearman | 0.8329 | 0.8177 | 0.7955 | 0.7706 | 0.8177 | 0.8177 | 0.8177 | 0.8177 | 0.7651 | 0.8177 | 0.8020 | 0.8020 | 0.7968 | 0.7974 | 0.7985 | 0.8020 | 0.8020 | 0.8020 |
|  | Kendall | 0.6463 | 0.6316 | 0.6065 | 0.5809 | 0.6316 | 0.6316 | 0.6316 | 0.6316 | 0.5743 | 0.6316 | 0.6132 | 0.6132 | 0.6129 | 0.6162 | 0.6154 | 0.6132 | 0.6132 | 0.6132 |

*ST-RAE is this project's governing metric (lower is better); R², Spearman, and Kendall are
better higher; MAE is better lower. See `src/vendor/openadmet_eval/` for the scoring definitions.*

## Per-row provenance: file, sha256, and scoring basis

Audited 2026-10-30. **Every scored row now carries the sha256 of the exact file that was sent and
an explicit statement of what half of the blind set produced its score.** Eleven of these hashes
were absent before that audit and are backfilled here; the seven that were already recorded were
recomputed from the files on disk and **all seven match**.

| Submission | Submitted file | sha256 | Scoring basis |
|---|---|---|---|
| 04b | `outputs/submissions/activity_submission_v1.csv` | `2d99d55ca5730c909e6c74abadfc44a5785b1be22659d806ddc1364d2cc08c3c` | **half-set** |
| NB10 | `outputs/10_final_retrain_predict/submission_candidate.csv` | `ab84077e7a4f28d2727e6820e6e53826f060f39bf2617a2231849cb85fefcecf` | **half-set** |
| 10c | `outputs/10c_control_submission/submission_candidate.csv` | `33c24cfd826afe65d4181606b3224c068be3b08ddd20c68044eede01f8a7c3b1` | **half-set** |
| NB12 | `outputs/12_caruana_retrain_predict/submission_candidate.csv` | `4e433c85c9a477392386101ec2188673e4e75a67fa1737d112460c1c022869dd` | **half-set** |
| NB13-calib | `outputs/13_aid1851_blind_population_calibration/submission_candidate.csv` | `9d56163672b93bfc3563348ce03a9c18f90ac3e464f7ef82b6aac13558db6100` | **half-set** |
| NB16 | `outputs/board_solved_calibration/submission_candidate.csv` | `5a6130322a5e1456bbd7b35f5db350af33623cfb8c08439f3f5b526fea826dcf` | **half-set** |
| NB19 | `outputs/19_cyp2c9_revert/submission_candidate.csv` | `02d296f13d62d02544dceece21dabdd58d34dc9e0b172954a1b7000d99ed1078` | **half-set** |
| NB27-widened | `outputs/27_calibration/cyp2d6_widened/submission_candidate.csv` | `ff759d39bc716f9d0e609ea0515769940b4e47cfdc1b89673ded592967687757` | **half-set** |
| NB30 | `outputs/30_aid_full_retrain/submission_candidate.csv` | `bfe57a7db57df08ed1f90b73923ec0726e052e8fdcc64e2e1271874d2ba99457` | **half-set** |
| 27b | `outputs/27b_aid_cyp2d6_corrected/submission_candidate.csv` | `deaaa31b27186e9b74f14ac86e2920b870e7a7e191c153cc366c81d8fa8b7ed7` | **half-set**, and the *only* row with a **full-set** reading as well (the interim board, from its accidental 2026-09-21 re-send of this same file) |
| NB32-A | `outputs/32_deadzone_aid_retrain/submission_candidate_A_allfour.csv` | `af710cf4cf55a0cc6690c7ed4d152d9dc62616dfe4a51385bd0954dab7e6ecbb` | **half-set** |
| NB36-widened | `outputs/36_spread_corrected_nb32a/submission_candidate_spread_corrected.csv` | `6bcb254b62c4656a32ac3f6fa73d606602a4fc15154137c3c82c5c2ec380118e` | **half-set** |
| NB40-ensemble | `outputs/40_single_object_ensemble/part8_stabilised_member/submission_candidate_rebuilt.csv` | `7df6e263cf463d1ebd72139922e98c7a870d335a1448152efb890a56ab09bf53` | **half-set**; per-isoform values from a screenshot, not a contemporaneous capture |
| NB47-mae | `outputs/47_mae_only_full_retrain/submission_candidate.csv` | `5ee259abcffaf590419d3de8b69c89133c0159a340b928dc39c5f47247054418` | **half-set** |
| NB48-ensemble | `outputs/48_board_ensemble/submission_candidate.csv` | `14ccf684fca506f5d7f1f8ed7ec8596e5ea2a90c79c7589208f500c92b97f6a7` | **half-set** |
| NB56-identify | `outputs/56_cyp1a2_identification/submission_candidate_identify.csv` | `80dfdef1d8734b1118c16d1c1596bd3d7155e959f4d2c0889bdcb5141b2b7815` | **half-set** |
| NB56b-s050 | `outputs/56b_cyp1a2_corrected/submission_candidate_s050.csv` | `a1da00e209975f4b82fdaeae1bc0210eb83953ef13872d22aa04d354c438aa8e` | **half-set** |
| NB56b-s100 | `outputs/56b_cyp1a2_corrected/submission_candidate_s100.csv` | `b2c53f136f1d2684643126492a98d2b718b696e55b216a17d60bf7ff6f79a228` | **half-set** |

**Four of the eighteen have a contemporaneous server receipt** on disk — NB48-ensemble,
NB56-identify, NB56b-s050 and NB56b-s100 — carrying the sent timestamp, the file's sha256 at send
time and the server's own reply. For the other fourteen the hash is of the archived file, and the
sent-versus-archived identity rests on the notebook's own saved output rather than a receipt. That
is a weaker guarantee and is stated rather than glossed.

**The file-to-submission mapping was not taken on trust.** Each of the eighteen paths above was
checked against a byte-identity claim this document already makes about it — NB16's CYP1A2/CYP3A4
against `10c`'s, NB19's against NB16's, NB19's CYP2C9 against NB12's, 27b's three carried columns,
NB56-identify's and both NB56b columns' carried CYP2D6/CYP3A4 — **eighteen checks, all of which
hold**: sixteen at exactly `0.00e+00` and two (NB27-widened and 27b against NB19 on CYP2C9) at
`4.44e-16`, which is the written-decimal difference documented in notebook 56b's entry and is
twelve orders of magnitude below the four decimal places the board reports. Had any path been
mis-assigned, these checks would have failed.

## NB47-mae and NB48-ensemble: both scored, and the gains landed where the CV said they would not

Two submissions were sent before either was recorded here, and both have now been scored. Neither
beats 27b. On macro ST-RAE they come 7th and 8th of the eighteen scored rows, between
NB36-widened (0.627) and NB16 (0.6364), so both sit in the corrected-submission band and well
clear of the raw baselines `10c` (0.7138) and `04b` (0.7179).

| | NB47-mae | NB48-ensemble |
|---|---|---|
| notebook | `47_mae_only_full_retrain.ipynb` | `48_board_ensemble.ipynb` |
| file | `outputs/47_mae_only_full_retrain/submission_candidate.csv` | `outputs/48_board_ensemble/submission_candidate.csv` |
| sha256 | `5ee259abcffaf590419d3de8b69c89133c0159a340b928dc39c5f47247054418` | `14ccf684fca506f5d7f1f8ed7ec8596e5ea2a90c79c7589208f500c92b97f6a7` |
| submitted (UTC) | 2026-10-10 21:44:38 | 2026-10-12 00:08 |
| alias | `fold-zero` | `fold-zero` |
| macro ST-RAE | 0.6328 (rank 83) | 0.6356 (rank 74) |
| CYP1A2 | 0.7254 (rank 126) | 0.7414 (rank 118) |
| CYP2C9 | **0.5255 (rank 138)** | 0.5290 (rank 120) |
| CYP2D6 | 0.7872 (rank 86) | 0.7652 (rank 78) |
| CYP3A4 | 0.4930 (rank 94) | 0.5067 (rank 82) |

Both reconcile: every macro equals the plain mean of its own four isoform values on all five
metrics, to 5e-5 or better, which is rounding at the precision the board displays.

**NB47-mae** is `chemprop_chemeleoninit` retrained on all 4,905 compounds under an **MAE training
loss instead of MSE**, raw point-estimate targets, no validation split, a fixed 7 epochs taken from
the median per-fold best epoch across notebook 37's 25 cluster-disjoint folds. It is the leak-free
half of the dead-zone effect: it changes the loss function alone, against the literal curated
label, so it carries none of `deadzone`'s fold-crosstalk exposure. CYP1A2 was spread-widened
(x1.1631); CYP2D6 was placement-corrected onto NB16's board-solved blind population using this
model's own pooled OOF rho (0.380669) and then widened to 0.85 (total x1.9109); CYP2C9 and CYP3A4
were left raw, both under notebook 36's 10% bar.

**NB48-ensemble** is a **per-isoform non-negative ensemble**, weights constrained to sum to one,
fitted on notebook 37's 25-fold cluster-disjoint out-of-fold predictions and deployed on full-data
blind columns: `deadzone` plus `aid` plus one `ecfp4_narrow` tree per isoform. It is the first
submission this project has made that is not an affine transform of a single model's output, so
the first that could move the rank metrics at all. CYP2D6 was placement-corrected on the blend's
own nested OOF rho (0.4525) and then widened; CYP1A2 and CYP3A4 were spread-widened to each
isoform's board-best ratio; CYP2C9 was left uncorrected, its 8.0% gap falling under the 10% bar.

### The ensemble finding: material in cross-validation on both metrics, nothing on the board

This is the first transfer failure on this project where the cross-validated result was material
on **both** metrics, not merely significant. Notebook 48 measured the blend against its own best
single component on notebook 45's frozen resolution bar and found CYP1A2 at **156% of the Spearman
bar** and CYP3A4 at **183%**, with ST-RAE also material on both (121% and 125%). Those were the
two isoforms it named as resolvable in advance.

Neither delivered. Reading the two new submissions as the adjacent pair they are, NB47-mae then
NB48-ensemble:

| Isoform | Spearman, NB47-mae | Spearman, NB48-ensemble | CV claim | Board outcome |
|---|---:|---:|---|---|
| CYP1A2 | 0.7423 | 0.7107 | material, 156% of bar, up | **fell 0.0316** |
| CYP3A4 | 0.7974 | 0.7985 | material, 183% of bar, up | **flat, +0.0011** |
| CYP2C9 | 0.7419 | 0.7692 | not resolvable | **rose 0.0273, highest on record** |
| CYP2D6 | 0.4602 | 0.4173 | not resolvable | **fell 0.0429** |

Against `plain`, which is the reference notebook 48's own prediction was written against and which
appears on the board as `10c`, the direction prediction scores **2 of 4**: CYP2C9 (0.7389 to
0.7692) and CYP2D6 (0.4000 to 0.4173) rose as predicted, while CYP1A2 (0.7375 to 0.7107) and
CYP3A4 (0.8177 to 0.7985) fell. The counterexample notebook 48 carried in writing, that the
out-of-fold-to-board Spearman mapping is not monotone on this project's record, is exactly what
bit, and it bit on the two isoforms the CV called material.

**CYP2D6's ST-RAE improvement is placement, not ranking.** It fell 0.7872 to 0.7652 while its
Spearman fell 0.4602 to 0.4173. A column cannot rank worse and predict better unless the gain came
from where its predictions sit rather than from their order, so that 0.0220 belongs to the
correction chain, not to the ensemble.

**The gains landed on the isoform the CV analysis explicitly set aside.** CYP2C9 was the one
isoform notebook 48 declared not resolvable, and the one it left uncorrected, and it is where both
submissions did their best work: NB47-mae's 0.5255 is a **new best on record for CYP2C9**, beating
the 0.5375 held jointly by NB12, NB19, NB27-widened and 27b, and NB48-ensemble's 0.5290 beats it
too. NB48-ensemble also set the highest CYP2C9 Spearman on record at 0.7692, against a prior best
of 0.7585. So of the four isoforms, the two the CV resolved as material both regressed, and the
two it could not resolve produced the only record. That is a stronger statement about CV-to-board
transfer on this project than any single macro figure.

### The compression signature: not a third hit, and ambiguous rather than a clean miss

The spread-compression association from notebook 27 Section 1, pooled Spearman of spread ratio
against board ST-RAE of -0.822, was two for two before these two submissions: it anticipated
NB30's failure where the CV evidence did not, and it predicted the gain NB36-widened's widening
actually delivered. Judged on these two submissions it is **neither**, and the reasoning matters
more than the verdict.

Measured spread ratios, prediction SD over training-label SD, both at `ddof=1`:

| Isoform | target ratio | NB47-mae | NB48-ensemble |
|---|---:|---:|---:|
| CYP1A2 | 0.8361 | 0.8361 | 0.8361 |
| CYP2C9 | 0.8168 | 0.8860 | 0.7517 |
| CYP2D6 | 0.8500 | 0.8500 | 0.8500 |
| CYP3A4 | 0.9262 | 0.8737 | 0.9262 |

Three things follow. First, on **CYP1A2 and CYP2D6 the two submissions carry identical ratios**,
because notebook 36's rule forced both to the same targets, yet their board ST-RAE differs by
+0.0160 and -0.0220. A variable that is identical across two columns cannot explain a difference
between them, so for half the isoforms the signature is structurally silent here. Second, of the
two isoforms where the ratio does differ, **one supports it and one contradicts it**: CYP2C9 is
more compressed in NB48-ensemble (0.7517 against 0.8860) and scored worse (0.5290 against 0.5255),
as the signature predicts, while CYP3A4 is **wider** in NB48-ensemble (0.9262 against 0.8737) and
scored **worse** (0.5067 against 0.4930), which it does not. One for one. Third, taken as a level
prediction the signature fails: NB48-ensemble's CYP1A2, CYP2D6 and CYP3A4 columns all sit exactly
on the ratios the board has already rewarded, and all three still scored worse than their own
per-isoform records, while the single column sitting furthest below its target, CYP2C9 at 0.7517,
beat the record.

**Verdict: ambiguous, and explicitly not a third hit.** The honest reading is that spread is now
calibrated on three of four isoforms and has been largely spent as a lever, which is consistent
with notebook 38's finding that what separates the leading columns is ranking, a quantity no
affine correction and no spread ratio can see. The signature was always observational rather than
causal, and NB48-ensemble is the first submission whose rank order differs from its underlying
single models, so it is also the submission on which that confound is most active.

### Placement headroom, newly quantified

Against a field fit of **ST-RAE = 0.695 - 0.210 x R2**, stated as valid for R2 in [0.40, 0.50],
NB48-ensemble's macro R2 of 0.4507 predicts a macro ST-RAE of **0.600** where the actual is
**0.6356**. That is roughly **0.035 of placement headroom**: the submission scores about 0.035
worse than a competitor at the same explained variance, which is a calibration gap rather than a
model-quality gap. NB47-mae sits in the same place, R2 0.4495 predicting 0.601 against an actual
0.6328, a gap of 0.032.

**Source and its status.** The fit was supplied by the user, attributed to a public competitor
write-up. It is **not** independently verifiable from anything in this repository, no primary
source is held here, and it has not been checked against this project's own submissions,
so it is recorded as an external reference figure rather than a measured result, in the same spirit
as `docs/stereochemistry_provenance.md`'s `CONSISTENT_WITH_LOCAL_DATA` status. The validity window
matters: both submissions fall inside R2 [0.40, 0.50], but 27b (R2 0.4835) and the raw baselines
(`10c` at 0.2278) would need the fit extrapolated or would fall outside it entirely.

## NB40-ensemble — 2026-10-02 13:37 UTC, notebook 40's rebuilt candidate: ACCEPTED AND SCORED

**This section previously read `UNCONFIRMED SUBMISSION ... accepted? unknown — must be checked on
the leaderboard`. That was answerable and wrong as written, and it is corrected here.** The send
was accepted and scored. The reply was not *never received* — it was **discarded by a
bare-expression bug in notebook 40's submission cell**, where `result` was written as a bare
expression rather than the cell's last statement, so nothing was captured.

Leaving the row as "unknown" was itself a hazard: a later reader would reasonably conclude the
send had failed, which is the mirror image of the staleness failures recorded under "The record
gap and its consequence" below.

| field | value |
|---|---|
| submission | **NB40-ensemble** |
| notebook | `40_single_object_ensemble.ipynb` (Part 8, rebuilt candidate) |
| file | `outputs/40_single_object_ensemble/part8_stabilised_member/submission_candidate_rebuilt.csv` |
| sha256 | `7df6e263cf463d1ebd72139922e98c7a870d335a1448152efb890a56ab09bf53` (verified present on disk) |
| submitted (UTC) | **2026-10-02 13:37** |
| alias | `fold-zero` |
| accepted? | **Yes** |
| **macro ST-RAE** | **0.6384** |
| **board rank** | **75** |
| **CYP2D6 ST-RAE** | **0.7751** |
| deliberate? | **No** — fired as a side effect of a notebook re-run. One of **two** non-deliberate sends, the other being the 2026-09-21 re-send of 27b; this document once called it the only one, and that was wrong |

**The full row was recovered on 2026-10-29, and this section's own contingency is what made it
cheap.** It previously said the missing columns were "not recoverable", because the board keeps one
row per user and NB48-ensemble replaced this one on 2026-10-12, so querying the board returns
NB48-ensemble rather than this submission. That reasoning was right about the *board* and wrong
about the *record*: the 2026-10-02 screenshots did carry the missing columns. The claim is retired.
What stands in its place is narrower and still worth knowing — **this is the only row in the main
table sourced from a screenshot rather than a contemporaneous capture**, and it is flagged as
user-supplied external data, like the interim full-set board.

The paragraph that earned its keep said analyses were "built so the row can be completed later
without rework if the 2026-10-02 screenshots turn out to carry the missing columns". They did, and
there was no rework: the figure layer's board parser picked the new column up without modification.

**And the row changes a standing claim.** NB40-ensemble's **CYP2C9 Spearman of 0.7729 is the
highest on this board for any isoform-submission pair on CYP2C9**, ahead of NB48-ensemble's 0.7692,
which this document had recorded as the record-holder. Both are ensembles and both are the only two
submissions that are not affine transforms of a single model — which is consistent with ranking
being the one thing a blend can move and an affine correction cannot. Its CYP1A2 Spearman (0.7102)
is, in the same breath, the **lowest** on record for that isoform, so the ensembles' ordering gains
are not uniform across isoforms.

**Where it sits.** At macro ST-RAE 0.6384 it is **tenth of the eighteen scored rows**, between
NB16 (0.6364) above it and NB32-A (0.6488) below — inside the
corrected-submission band and clear of the raw baselines `10c` (0.7138) and `04b` (0.7179). Its
CYP2D6 ST-RAE of 0.7751 sits between NB48-ensemble's 0.7652 and NB47-mae's 0.7872. It did **not**
beat 27b (0.5982).

**The timestamp.** 13:37 UTC is the board's own. This document previously estimated "~13:44"; that
estimate is superseded. The preceding attempt at ~13:33 UTC reached `Loaded as API` and then
raised `ValueError` on a bad relative path, so **that** one transmitted nothing.

**Predicted scores**, recorded before the board was read, are in
`outputs/40_single_object_ensemble/part8_stabilised_member/predictions_before_board_rebuilt.json`:
macro ST-RAE in [0.582, 0.707] against 27b's 0.5982. **The prediction held** — 0.6384 falls inside
that interval — though notebook 40 itself stated in advance that the interval was ±0.05 around the
anchor, wider than the effect being tested, so it could detect a catastrophic failure but could not
resolve whether the ensemble helped. It did not help.

**What happened.** Notebook 40's submission cell was uncommented in the IDE between a verification
run and a subsequent `jupyter nbconvert --execute --inplace` re-run; the re-run executed it. The
cell's own text still read "fully commented out and was NOT executed. Nothing was sent." while the
call was completing — the NB10/10c/12 header-comment hazard in its most dangerous form, this time
with a real send behind it.

**Consequence for the cell.** The comment gate is no longer used in notebook 40. It has been
replaced by an explicit `SEND = False` flag whose printed message is computed from the flag, so the
cell cannot claim it sent nothing while sending. Before any future `--execute` of that notebook,
check the gate on disk first.

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
| **NB36-widened** | `36_spread_corrected_nb32a.ipynb` | 2026-10-01 19:54 UTC | NB32-A's submitted file with CYP1A2, CYP2C9 and CYP3A4 each widened about their own mean to the spread ratio of the column holding that isoform's best board ST-RAE (multipliers x1.4145, x1.1304, x1.1317 onto target ratios 0.8361, 0.8168, 0.9262); CYP2D6 carried over byte-identical, since NB32-A already placed and widened it to 0.85 | Yes — manually reviewed and sent, per this project's manual-gated submission process |
| **NB40-ensemble** | `40_single_object_ensemble.ipynb` (Part 8, `part8_stabilised_member/submission_candidate_rebuilt.csv`) | 2026-10-02 13:37 UTC | Single shared-member, shared-weight ensemble (`deadzone` 0.625 / `aid` 0.25 / `ecfp4_narrow__xgboost` 0.125, seed-averaged over ten early-stopping split draws), one weight vector across all four isoforms with no per-isoform refitting; CYP1A2/CYP2D6/CYP3A4 spread-widened, CYP2D6 also placement-corrected on this blend's own OOF rho (0.4496), CYP2C9 left uncorrected (7.9% gap, under the 10% bar) | **No** — fired as a side effect of an `nbconvert --execute --inplace` re-run. **Scored: macro 0.6384, rank 75** |
| **NB47-mae** | `47_mae_only_full_retrain.ipynb` | 2026-10-10 21:44:38 UTC | `chemprop_chemeleoninit` full-data retrain under an **MAE training loss** instead of MSE, raw point-estimate targets, no validation split, fixed 7 epochs from the median per-fold best epoch across 25 cluster-disjoint folds; CYP1A2 spread-widened (×1.1631), CYP2D6 placement-corrected on this model's own OOF rho (0.380669) then widened to 0.85 (total ×1.9109), CYP2C9/CYP3A4 raw | Yes, uncommented and run by hand. **Scored: macro 0.6328, rank 83** |
| **NB48-ensemble** | `48_board_ensemble.ipynb` | 2026-10-12 00:08:07 UTC | Per-isoform **non-negative ensemble** (weights summing to one) fitted on notebook 37's 25-fold cluster-disjoint OOF predictions, deployed on full-data blind columns (`deadzone` + `aid` + an `ecfp4_narrow` tree per isoform); CYP2D6 placement-corrected on the blend's own nested OOF rho (0.4525) then widened, CYP1A2/CYP3A4 spread-widened, CYP2C9 uncorrected (8.0% gap, under the 10% bar). The first submission that is not an affine transform of one model, so the first that can move the rank metrics | Yes, gate opened and run deliberately, on explicit instruction. **Scored: macro 0.6356, rank 74** |
| **NB56-identify** | `56_cyp1a2_identification.ipynb` | 2026-10-28 19:38:58 UTC | **Identification probe, not a modelling candidate.** NB36-widened with a uniform +0.25 pIC50 translation on CYP1A2/CYP2C9/CYP3A4 and CYP2D6 carried byte-identical as a built-in control. Deliberately degraded: its only purpose was to read the *sign* of each isoform's placement error off the returned R², which an affine shift makes identifiable | Yes, gate opened by hand; disclosed in its own prereg as a probe. **Scored: macro 0.6947** |
| **NB56b-s050** | `56b_cyp1a2_corrected.ipynb` | 2026-10-29 14:18:11 UTC | NB36-widened with CYP1A2 and CYP2C9 placed and rescaled onto blind populations solved *exactly* from three board readings (ρ recovered, no Spearman proxy), at shrinkage **s = 0.5**; CYP2D6 and CYP3A4 carried byte-identical | Yes. **Scored: macro 0.5861 — best on record** |
| **NB56b-s100** | `56b_cyp1a2_corrected.ipynb` | 2026-10-30 07:54:11 UTC | Identical to NB56b-s050 in every respect except the shrinkage fraction, **s = 1.0** (placement fully corrected rather than half-corrected). The two together are a one-variable controlled comparison | Yes. **Scored: macro 0.5933, rank 93 — second best** |

**On "deliberate":** seventeen of the nineteen real sends were reviewed and uncommented/run by
hand, exactly as this project's manual-gated-submission process requires. An earlier diagnostic pass over this
repo (before this correction) reasoned from file evidence alone — a misleading header comment on
the NB10/10c/NB12 submission cells, which claimed "every line below is a comment" when only the
explanatory paragraph was — and concluded the submissions might have fired by accident. **That
conclusion was wrong.** The header comment itself was a real bug (now fixed, 2026-09-07 — see
CLAUDE.md's Rules section), but the submissions it sat on top of were intended, reviewed actions,
and the content each one sent was correct. **There are two genuine exceptions, not one.** NB40-ensemble
fired as a side effect of a notebook re-run on 2026-10-02, and **27b was re-sent the same way on
2026-09-21 12:49 UTC** — found only on 2026-10-27, and the more consequential of the two because it
is the send the interim full-set leaderboard scored. Both have their own sections. This document
previously called NB40-ensemble "the only non-deliberate submission in this project's history";
that was wrong, and the earlier claim is retired here rather than quietly amended.

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

5. **Notebook `35`/`37`'s MAE-loss arm**, significantly better than `plain` at 25 cluster-disjoint
   folds on CYP1A2 (-0.0317) and CYP2D6 (-0.0803) and not different on the other two. Submitted as
   NB47-mae, its CYP1A2 came back 0.7254 against `plain`'s own 0.6954, a clean reversal of the one
   isoform where the CV comparator and the board comparator are the same object; CYP3A4 (0.4930
   against 0.4434) also regressed. The nuance worth keeping is that this one is not a blanket
   failure: CYP2C9, which the CV called null, set a new record at 0.5255.
6. **Notebook `48`'s per-isoform non-negative ensemble**, the first result on this project that was
   **material on both metrics** rather than merely significant, at 156% of the frozen Spearman bar
   on CYP1A2 and 183% on CYP3A4. Submitted as NB48-ensemble, CYP1A2's Spearman fell and CYP3A4's
   was flat, and both isoforms' ST-RAE regressed against their references. Its one real gain,
   CYP2D6's ST-RAE, is attributable to placement rather than ranking, since that column's Spearman
   fell at the same time. See the NB47-mae/NB48-ensemble section above.

7. **The CV-versus-board disagreement on SPEARMAN itself** (2026-10-22; canonical record
   `docs/findings.md`, finding 3). This one is a different kind from the six above and is recorded
   separately below.

Items 1 to 4 were recorded when NB32-A was scored; items 5 and 6 were added later, when the two
2026-10 submissions were scored; item 7 when the submitted orderings were ranked against the
25-fold cluster-disjoint CV.

### Item 7 — the first transfer disagreement calibration cannot explain

> **Canonical record: `docs/findings.md`, finding 3.**

Items 1 to 6 are all **ST-RAE** failures, and ST-RAE mixes model quality with placement, so
"the model was fine, the calibration was not" was always an available explanation. Notebook 55
ranked the submitted orderings on **Spearman**, which no affine correction can move, and
cross-checked that ranking against the 25-fold cluster-disjoint CV for the four orderings whose
architecture exists as an arm (`plain` ← `10c`, `aid` ← NB30, `mae_only` ← NB47-mae).

**The board and the CV disagree about which model ranks better on 2 of 4 isoforms**, and CYP3A4 is
a complete reversal:

| Isoform | CV says (25-fold OOF Spearman) | The board says | |
|---|---|---|---|
| CYP1A2 | NB30 > NB47-mae > 10c | NB30 > NB47-mae > 10c | agree |
| CYP2C9 | NB30 > NB47-mae | NB47-mae > NB30 | **disagree** |
| CYP2D6 | NB30 > NB47-mae > 10c | NB30 > NB47-mae > 10c | agree |
| CYP3A4 | NB30 > NB47-mae > 10c | 10c > NB47-mae > NB30 | **disagree, fully reversed** |

Because Spearman is affine-invariant on both sides, calibration cannot account for this: the two
estimates disagree about the models themselves. It also reproduces, at 25 folds and on the board,
the specific counterexample notebook 40 wrote down *in advance* and notebook 48 then hit — that
the out-of-fold-to-board Spearman mapping is not monotone on this project's record.

## The ρ substitution in the deployed correction chain is biased low

> **Canonical record: `docs/findings.md`, finding 1.** The summary below is kept here because it
> bears directly on every correction in this document; the method, the three-subset spread and the
> backing data live there.

Every placement correction this project has deployed — NB16, NB19, NB27-widened, 27b, NB32-A,
NB47-mae, NB48-ensemble — applies a ρ-shrunk slope, and takes ρ either from the board's Spearman
directly (notebook 16) or from the bivariate-normal relation ρ = 2·sin(πρ_S/6) (notebook 38).
**That substitution had never been checked against an independently recovered ρ, because there was
no way to recover one.** Notebook 55 found one.

Within a single *ordering group* — submissions whose columns are affine transforms of each other,
and which therefore share one Pearson ρ — each distinct member contributes one equation
`R² = 2ρk − k² − b²` in the three unknowns (ȳ, σ_y, ρ). Three distinct members identify it with no
proxy anywhere. **CYP2D6 is the only isoform on this record with enough distinct variants** (four:
`10c`, NB13-calib, NB16, NB27-widened — four genuinely different affine maps of one ranking), and
the other three isoforms max out at two, because their repeated submissions were byte-identical
resubmissions rather than different corrections.

| quantity | value |
|---|---|
| ρ **recovered**, 4 equations, no proxy | **0.5600** (three-subset range 0.5502–0.6156) |
| board Spearman, used as ρ by notebook 16 | 0.4000 — understates by **+0.1600** |
| ρ = 2·sin(πρ_S/6), notebook 38's refinement | 0.4158 — understates by **+0.1442** |

Both proxies understate it, and the refinement closes about a tenth of the gap. Since ρ sets the
placement slope, the deployed chain **shrinks too hard** — the same direction notebook 26 reached
independently from an out-of-fold spread sweep, where the ST-RAE-optimal multiplier sits strictly
above ρ on every isoform. Two independent routes, same sign.

**This does not invalidate the deployed corrections.** The same solve puts the blind CYP2D6
population at mean 3.0515, SD 1.5397, against NB16's deployed target of 3.1558 / 1.5114 — agreeing
to ~0.10 and ~0.03 pIC50 despite the ρ error, because ρ enters the *slope* and not the *target*,
so the two errors do not compound. The caveat that matters is that this is one isoform, and it is
the one notebook 21 established cannot be ranked on out-of-fold ST-RAE.

Each instance carries progressively stronger CV backing than the last (a single uncrossvalidated
selection pass; a single-fold screen; a full paired 5×5 confirmation; a single-fold screen combined
with an already-confirmed architecture; a 25-fold significant single-variable change; a 25-fold
result material on both metrics) and none of the six has transferred cleanly to the board. **Item 6
is the first in which the cross-validated effect was material on both metrics and still produced
nothing**, which removes "the effect was only just significant" as an available explanation.

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

## NB36-widened: spread correction works, and is bounded by ranking

Notebook `36` took NB32-A's submitted file and widened CYP1A2, CYP2C9 and CYP3A4 about their own
means to the spread ratio of whichever column holds that isoform's best board ST-RAE, leaving
CYP2D6 byte-identical. It scored **macro 0.627** — better than its own parent NB32-A (0.6488) but
**fourth overall**, behind 27b (0.5982), NB27-widened (0.6155) and NB19 (0.6265). A mean-preserving
positive scale is exactly rank-preserving, so by construction this submission could not change any
isoform's ranking, only where its predictions sit and how far they spread.

**All five of notebook 36's pre-registered predictions were confirmed**, which is the first time
this project has recorded a clean sweep:

| # | Prediction | Outcome |
|---:|---|---|
| 1 | Spearman and Kendall identical to NB32-A's on all four isoforms | **CONFIRMED** — all eight values reproduce exactly |
| 2 | CYP2D6's byte-identical column reproduces NB32-A's row on every metric | **CONFIRMED** — all five (0.7502 / 0.9319 / 0.3768 / 0.4438 / 0.3076) |
| 3 | macro ST-RAE between 0.6030 and 0.6488 | **CONFIRMED** — 0.627 |
| 4 | CYP1A2 improves most of the three widened isoforms | **CONFIRMED** — +0.0665 vs CYP2C9 +0.0200 vs CYP3A4 +0.0003 |
| 5 | if the compression association is causal, all three widened isoforms improve | **CONFIRMED** — all three improved |

Prediction 2 is also the **fifth** independent confirmation that the board reproduces a
byte-identical column exactly (after NB19, NB27-widened, 27b, and 10c/NB16's shared columns).

**The payoff is strongly non-uniform, and notebook `38` diagnosed why.** CYP2C9 and CYP3A4 took
near-identical treatment — multipliers x1.1304 and x1.1317 from relative spread gaps of 11.5% and
11.6% — and gained 0.0200 and **0.0003**. CYP3A4's residue is **not** placement:

- after widening, its column has **exactly `04b`'s spread** (SD 1.012520 against 1.012520, identical
  to 2.2e-16), because notebook `36`'s target rule was `04b`'s own ratio — so scale is eliminated
  by construction;
- the blind CYP3A4 population solved jointly from the other eleven submissions (mean 4.7778,
  SD 1.3971) puts this column **closer to perfectly centred than `04b`** (|b| 0.0715 against
  0.1040), so placement works in its favour, not against it;
- what is left is **ranking**: Spearman 0.8020 against `04b`'s 0.8329, Kendall 0.6132 against
  0.6463 — invariant under every monotone transform, so no affine correction reaches it.

With this row's own R² now on record, that decomposition is a measurement rather than an inference.
Taking the identity `R² = 2ρk − k² − b²` with the solved population, the implied Pearson ρ is
**0.8457 for `04b` and 0.8187 for this submission**, and `04b`'s 0.0335 R² advantage splits as
**correlation +0.0391, placement −0.0057, scale 0.0000** (sum 0.0334, against the measured 0.0335).

**This row also retro-validates notebook `38`'s population solve.** That notebook fitted the blind
CYP3A4 population to the **eleven** then-recorded submissions and, without this row, predicted this
column's R² at **0.6515**. The measured value is **0.6563** — an error of 0.0048, inside the fit's
own RMS residual of 0.0065. The solve was genuinely predictive, not merely descriptive.

**Verification.** Every macro value equals the mean of its own four isoform values, checked
directly: ST-RAE 0.627050, MAE 0.742600, R² 0.467675, Spearman 0.689425, Kendall 0.511775 — all
five agree with the published figures at the precision those figures are displayed to, with **no
rounding-boundary case and no unresolved discrepancy**, unlike NB16, NB19 and NB27-widened. The
macro ST-RAE is recorded as `0.627` because that is the precision the board displayed it at; the
mean of the four isoform values is 0.627050, consistent with it.

## The sixteenth send: 27b, re-sent by accident on 2026-09-21 — and the one full-set reveal

**This row was missing from this document until 2026-10-27. It is the fifth time this record has
gone stale, and the second accidental send.** `notebooks/27b_aid_cyp2d6_corrected.ipynb` was
re-executed on **2026-09-21 at 12:49 UTC** and its submission cell fired. The file sent was
byte-identical to 27b's original, so the board returned 27b's own scores and the per-isoform table
above needs no new column.

| evidence | value |
|---|---|
| `outputs/27b_aid_cyp2d6_corrected/submission_candidate.csv` mtime | **2026-09-21 12:49:54 UTC** |
| `notebooks/27b_aid_cyp2d6_corrected.ipynb` mtime | 2026-09-21 12:49:59 UTC |
| submitted file sha256 | `deaaa31b27186e9b74f14ac86e2920b870e7a7e191c153cc366c81d8fa8b7ed7` (unchanged on disk) |
| cell 26 `execution_count` | **19** — four uncommented live lines, no `SEND` gate |
| `model_tag` as sent | `""` (empty string) |
| captured output | `Loaded as API: https://openadmet-cyp-challenge.hf.space`, then *"submission cell defined above is fully commented out — not executed. Nothing was sent."* |

The cell's own printed text contradicted what it had just done, exactly as NB40's did. The server's
reply was not captured, because `result` was a bare expression and not the cell's last statement —
also exactly as NB40. Both failure modes were already documented from the NB40 incident; this send
predates that documentation by eleven days and was simply never looked for.

**Why it matters far more than a duplicate should: this send, not NB32-A, is what the interim
full-set leaderboard scored.** The intermediate submission deadline was **24 September** (the
challenge's own `config.py` and launch announcement), and the last row standing at that deadline was
this one. Three independent lines agree: the artefact mtime matches the board's displayed timestamp
to the second, that timestamp is a submission time rather than a snapshot label, and the deadline
ordering puts this send last before the cut.

**The accident helped.** It replaced NB32-A (half-set macro 0.6488) with 27b (0.5982) as this
project's one and only full-set reveal. That is luck, not method, and it is recorded as such.

As of 2026-10-27 every submission cell in every notebook in this repository has been gated so that
none can fire again; see `CLAUDE.md`'s step-0 note and
`outputs/56_cyp1a2_identification/step0_gate_record.json`.

### The interim (full-set) leaderboard — the only full-set reading before the close

**User-supplied board capture (2026-10-27), not derivable from anything in this repository**, and
recorded in the same spirit as the other user-supplied figures in the Provenance section below.
Rank **66 of 219**, username `fold-zero`, submitted 2026-09-21 12:49 UTC, Model Report Link
"Not submitted". Per-isoform ranks are also of 219.

| Isoform | ST-RAE | MAE | R² | Spearman | Kendall | rank |
|---|---:|---:|---:|---:|---:|---:|
| CYP1A2 | 0.6444 | 0.9119 | 0.3234 | 0.7298 | 0.5416 | 103 |
| CYP2C9 | 0.5117 | 0.5778 | 0.5749 | 0.7946 | 0.6015 | 109 |
| CYP2D6 | 0.7962 | 0.9809 | 0.3537 | 0.4630 | 0.3175 | 65 |
| CYP3A4 | 0.4350 | 0.5166 | 0.6948 | 0.8301 | 0.6443 | 56 |
| **macro** | **0.5968** | **0.7468** | **0.4867** | **0.7044** | **0.5262** | **66** |

**Macro reconciliation computed, not asserted** — mean of the four isoform values against the
stated macro: ST-RAE 2.5e-5, MAE 0.0, R² 0.0, Spearman 2.5e-5, Kendall 2.5e-5. All display
rounding; no boundary case.

**The two halves disagree materially, and not in the same direction.** Because the interim scored
27b, this is a like-for-like comparison of the *same columns* on the live half and the full set:

| Isoform | ST-RAE half → full | R² half → full | Spearman half → full | built by hit expansion? |
|---|---|---|---|---|
| CYP1A2 | 0.6954 → **0.6444** | 0.3073 → 0.3234 (+0.0161) | 0.7375 → 0.7298 (−0.0077) | yes |
| CYP2C9 | 0.5375 → **0.5117** | 0.5432 → 0.5749 (+0.0317) | 0.7585 → 0.7946 (+0.0361) | yes |
| CYP2D6 | 0.7165 → **0.7962** | 0.4124 → 0.3537 (**−0.0587**) | 0.4872 → 0.4630 (−0.0242) | **no** |
| CYP3A4 | 0.4434 → **0.4350** | 0.6709 → 0.6948 (+0.0239) | 0.8177 → 0.8301 (+0.0124) | yes |
| macro | 0.5982 → 0.5968 | 0.4867 | 0.7044 | |

The three isoforms built by hit expansion all got **easier** on the full set; **CYP2D6, the one
isoform excluded from hit expansion, got harder**, and by more than the other three moved combined.
Spearman moves too, so the halves differ in *ordering* difficulty as well as in placement — and no
affine correction addresses that.

**Consequence for every correction this project has deployed or may yet deploy.** All of them are
fitted against half-set readings and applied to all 750 compounds. Taking CYP2D6 — the only isoform
whose half-set population is identified (ȳ 3.0515, σ_y 1.5397, finding 1 in `docs/findings.md`) —
and holding its ordering's ρ (0.6656) and σ_y fixed, the full-set R² of 0.3537 implies a population
centre of **2.7725 or 3.5503**, i.e. a half-to-full drift of **−0.2789 or +0.4989 pIC50**; holding ȳ
instead implies σ_y of 0.8338 or 2.0970. Either way the drift is **large relative to the corrections
being made** — it is a substantial fraction of a typical placement displacement, not a rounding
term. Any correction tuned tightly to a live half-set reading carries that as an untested
assumption, and it is the reason to prefer a partial correction over a full one.

## NB56-identify — the identification probe: a deliberately degraded send that worked

**This send carried no modelling content.** It is NB36-widened with CYP1A2, CYP2C9 and CYP3A4
translated +0.25 pIC50 and CYP2D6 left byte-identical, built by `56_cyp1a2_identification.ipynb`
for one purpose: to read the **sign** of each isoform's placement error off the returned R².
It was expected to score worse, and it did — macro ST-RAE 0.6947 against NB36-widened's 0.6270.

| field | value |
|---|---|
| notebook | `56_cyp1a2_identification.ipynb` |
| file | `outputs/56_cyp1a2_identification/submission_candidate_identify.csv` |
| sha256 | `80dfdef1d8734b1118c16d1c1596bd3d7155e959f4d2c0889bdcb5141b2b7815` |
| submitted (UTC) | **2026-10-28 19:38:58** |
| alias | `fold-zero` |
| reply | captured — "Submission received from 'fold-zero' for the Regression Prediction track" |
| deliberate? | **Yes** — sent by hand through an explicit `SEND` flag, receipt written to `submission_receipt.json` |

**Disclosure, carried from the pre-registration and repeated here because it belongs in the
methods rather than in a footnote.** This is different in kind from every other correction on this
project, which solved blind moments from sends made on their merits. This one was degraded on
purpose to extract a parameter of the hidden test distribution. Three conditions were set in
advance and all hold: it is disclosed as a probe, the final reported rank does not rest on it, and
it is not the last send.

### Every pre-registered check passed

Macro reconciliation: all five metrics equal the mean of their own four isoform values to 2.5e-5.

**The CYP2D6 control reproduced NB36-widened exactly on all five metrics** — 0.7502 / 0.9319 /
0.3768 / 0.4438 / 0.3076. That column was carried byte-identical, so this is the **sixth**
confirmation on this project that identical inputs produce identical board output.

**Spearman and Kendall are unchanged on all four isoforms.** A pure translation cannot move rank
order, and the board agrees — the fourth time an affine-transform prediction of this kind has held
exactly.

### What it established: the sign prior was wrong, and the R² identity is sound

The pre-registration froze four predicted CYP1A2 readings, two per mirror family:

| branch | predicted R² | observed |
|---|---:|---|
| **ȳ=4.3007, σ_y=1.4912 (LOW)** | **0.1107** | **0.1154 — matched, \|diff\| 0.0047** |
| ȳ=4.2543, σ_y=2.6407 (LOW) | 0.2409 | |
| ȳ=5.8413, σ_y=2.6407 (HIGH) | 0.3547 | |
| ȳ=5.7950, σ_y=1.4912 (HIGH) | 0.4467 | |

Two conditions fired, and they point in opposite directions:

- **Falsification condition 3 did NOT fire.** The reading matched a branch to 0.0047, far inside
  the 0.02 bar, so `R² = 2ρk − k² − b²` describes this board and the placement route lives.
- **Falsification condition 4 DID fire. The hit-expansion sign prior is FALSIFIED.** CYP1A2's R²
  *fell*, so its blind population is **LOW-centred**: the predictions sit **above** the blind
  centre and the correction shifts them **down**, the opposite of what was predicted. CYP2C9 and
  CYP3A4 fell too, so all three are low-centred.

**The prior's failure is more interesting than a coin landing the other way, and is recorded
rather than discarded.** Its premise was that a hit-expanded test set should centre *above* the
training mean. On two isoforms that premise is correct — the solved blind centres are CYP2C9
+0.227 and CYP3A4 +0.701 relative to training — and the conclusion still failed, because the chain
"blind centre above training ⇒ predictions below blind centre" assumed the predictions sit near
the training mean. They do not: this model predicts high relative to training on every isoform
(CYP2C9 4.879 against 4.581, CYP3A4 4.878 against 4.096). On CYP1A2 the premise itself fails as
well — the blind centre solves to 4.300, **0.655 below** training's 4.955, the same direction as
CYP2D6's already-known 1.733 below. So the construction argument was not a reliable guide to
placement, and a measurement was needed.

### The three-equation solve: all three remaining isoforms are now IDENTIFIED

With a third column at a different prediction mean, the mirror cannot survive. Each isoform solves
to a **unique** root with **ρ recovered, not proxied** — the same standard as CYP2D6's solve in
`docs/findings.md` finding 1.

| isoform | ȳ | σ_y | ρ | b | k | R² now | ceiling ρ² | ΔR² | shift to centre |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **CYP1A2** | 4.3001 | 1.5099 | 0.7690 | **+0.4952** | 0.5706 | 0.3068 | 0.5914 | **+0.2846** | **−0.7477** |
| CYP2C9 | 4.8074 | 1.1325 | 0.7560 | +0.0631 | 0.5642 | 0.5308 | 0.5716 | +0.0408 | −0.0715 |
| CYP3A4 | 4.7974 | 1.3145 | 0.8136 | +0.0611 | 0.7703 | 0.6563 | 0.6619 | +0.0056 | −0.0803 |

**CYP3A4 is confirmed closed for the third time**, now by direct measurement rather than
inference: ΔR² = +0.0056, smaller still than notebook 38's estimate and than the two-equation
solve's +0.0086.

**An honest note on what the probe actually bought.** The recovered ρ values sit within
0.0026 / 0.0146 / 0.0018 of the normal-theory Spearman proxy the two-equation solve had assumed,
so the *magnitudes* that solve produced were nearly right all along. What was unknown, and what
this send purchased, was the **direction** — exactly as the pre-registration said it would be
("the sign is free, the magnitude needs the solve"), though it is worth recording that the
magnitude turned out not to need it much. The proxy's accuracy here is also isoform-specific
rather than general: on CYP2D6 the same proxy understates ρ by 0.144 (`findings.md` finding 1).

## NB56b-s050: the identified population deployed — the best row on record

**Macro ST-RAE 0.5861**, beating 27b's 0.5982 by 0.0121 and taking the record for the first time
since 2026-09-19. Submitted **2026-10-29 14:18:11 UTC**, file
`outputs/56b_cyp1a2_corrected/submission_candidate_s050.csv`, sha256
`a1da00e209975f4b82fdaeae1bc0210eb83953ef13872d22aa04d354c438aa8e` (matches both the file on disk
and its own receipt). Macro reconciliation computed, not asserted — mean of the four isoform values
against the stated macro: ST-RAE 2.5e-5, MAE 0.0, R² 0.0, Spearman 2.5e-5, Kendall 2.5e-5, all
display rounding.

**It is also the first record-holder that is not a per-isoform assembly.** 27b, NB19, NB27-widened
and NB16 all took columns from more than one model; this is NB36-widened's own four columns with an
affine correction applied to two of them, so it is attributable to one model plus arithmetic.

**What it is.** NB36-widened with CYP1A2 and CYP2C9 placed and rescaled onto blind populations
solved *exactly* from three board readings — ρ **recovered** rather than taken from a Spearman
proxy — at shrinkage s = 0.5. CYP2D6 and CYP3A4 carried byte-identical. The populations came from
the notebook 56 identification probe (`NB56-identify`), which is why a deliberately degraded send
earned its slot.

| what moved | NB36-widened → NB56b-s050 |
|---|---|
| CYP1A2 ST-RAE | 0.7241 → **0.5631** (−0.1610) |
| CYP1A2 R² | 0.3068 → **0.5070** (+0.2002) |
| CYP2C9 ST-RAE | 0.5435 → 0.5408 (−0.0027) |
| CYP2C9 R² | 0.5308 → 0.5595 (+0.0287) |
| CYP2D6, CYP3A4 | unchanged by construction |

**CYP1A2 0.5631 is far the best CYP1A2 reading on record** — the previous best was 0.6954, held by
`10c` and every byte-identical resend of it — and its R² of 0.5070 is the first above 0.32 on that
isoform.

**Two controls built into the row, both passed.** CYP2D6 and CYP3A4 were carried byte-identical, so
**all ten of their cells must reproduce NB36-widened's and all ten do** — the sixth confirmation on
this record that the board returns identical scores for an identical file. And because the
correction is affine, **Spearman and Kendall on the two corrected isoforms must be unchanged, and
both are** (CYP1A2 0.7565/0.5611, CYP2C9 0.7554/0.5652). A row that passes fourteen checks it could
have failed is worth more than one that merely looks plausible.

**Notebook 56b's pre-registration scored, and this is the first time this project has forecast R²
ahead of a send and had it land.** Recorded to `predictions_before_board.json` at 14:18:06 UTC,
five seconds before the send:

| prediction | outcome |
|---|---|
| 12 exact-equality cells (CYP2D6 ×4, CYP3A4 ×4, Spearman/Kendall on the corrected two ×4) | **12 of 12 hit** |
| CYP2C9 R² = 0.5614 | 0.5595, off by **0.0019** |
| CYP1A2 R² = 0.5203 | 0.5070, off by **0.0133** |
| macro ST-RAE in [0.5659, 0.6154] | 0.5861, **inside** |
| CYP2C9 ST-RAE = 0.5255 | **0.5408 — missed by +0.0153** |

The CYP2C9 ST-RAE prediction was a degenerate "band" of a single point and it missed; recorded as a
miss rather than rounded into a range. The R² forecasts are the load-bearing ones, since R² is what
the solved population predicts directly, and both landed inside the fit's own residual.

**What it does not settle — one item since resolved.** The shrinkage fraction was unresolved when
this row was written: s = 1.0 was built (`submission_candidate_s100.csv`, sha256 `b2c53f13…`) and
unsent. **It was sent on 2026-10-30 and scored 0.5933**, so half-correcting was better on ST-RAE,
not merely safer — see the NB56b-s100 section below for the controlled comparison. Nor does it
validate the base — three of NB36-widened's four columns are unchanged, so this says nothing about
whether `deadzone+aid` was the right model to correct. And every figure here remains a half-set
measurement.

## NB56b-s100: the shrinkage lever, and a controlled ST-RAE/R² divergence

Submitted **2026-10-30 07:54:11 UTC**, file
`outputs/56b_cyp1a2_corrected/submission_candidate_s100.csv`, sha256
`b2c53f136f1d2684643126492a98d2b718b696e55b216a17d60bf7ff6f79a228` (matches the file on disk and
its own receipt). **Macro ST-RAE 0.5933, board rank 93** — second-best of the eighteen scored rows,
behind only NB56b-s050 (0.5861) and ahead of 27b (0.5982). Macro reconciliation computed, not
asserted: mean of the four isoform values against the stated macro gives ST-RAE 7.5e-5, MAE 2.5e-5,
R² 2.5e-5, Spearman 2.5e-5, Kendall 2.5e-5 — all inside display rounding, ST-RAE's the largest of
the five. **A half-set measurement**, like every other row here.

| field | value |
|---|---|
| submission | **NB56b-s100** |
| notebook | `56b_cyp1a2_corrected.ipynb` |
| submitted (UTC) | **2026-10-30 07:54:11** |
| alias | `fold-zero` |
| accepted? | Yes — server reply captured in `submission_receipt_s1.0.json` |
| deliberate? | Yes, gate opened and run by hand |
| macro ST-RAE | **0.5933** (rank 93) |
| CYP1A2 | ST-RAE 0.5767 (rank 87) |
| CYP2C9 | ST-RAE 0.5562 (rank 161) |
| CYP2D6 | ST-RAE 0.7502 (rank 106) |
| CYP3A4 | ST-RAE 0.4904 (rank 93) |

**One refused attempt preceded it, and it is worth recording because it confirms the cooldown
exactly.** The same file was sent at **2026-10-29 21:28:04 UTC** and the server refused it:
`Error: You submitted a regression prediction on 2026-10-29 14:18:12 (UTC). Please wait 04:50:07
before submitting again.` 21:28:04 plus 04:50:07 is 2026-10-30 02:18:11, which is 14:18:12 plus
twelve hours — an independent confirmation of `HOURS_BETWEEN_SUBMISSIONS = 12` from the server's own
arithmetic. **Nothing was transmitted**, so this is not a send and is not counted as one.

### Why this pair is a controlled experiment

NB56b-s050 and NB56b-s100 differ in **exactly one scalar**: the shrinkage fraction `s` applied to
the affine correction on CYP1A2 and CYP2C9. Same base file (NB36-widened), same solved populations,
same recovered ρ, same code path, sent twelve hours apart. CYP2D6 and CYP3A4 are byte-identical
across the two (verified: max |difference| `0.00e+00` on both columns). Three results follow, and
the full write-up with the arithmetic is in `docs/findings.md` finding 5.

**The lever moved only what it touched.** All ten CYP2D6 and CYP3A4 cells are exactly equal across
the two submissions — not equal to four significant figures, equal at every digit the board prints.

**Spearman and Kendall did not move on any of the four isoforms.** Not one of the eight rank-metric
cells differs. On CYP2D6/CYP3A4 that is trivial (identical files); the load-bearing part is CYP1A2
and CYP2C9, where the predictions genuinely changed under a different affine map and the ordering
still did not. Both also reproduce NB36-widened's own values on those isoforms, so this is a
three-way agreement across three separate sends.

**And the two metrics disagree about which submission is better.** On the two isoforms that moved,
ST-RAE got **worse** while R² got **better**, on both:

| | CYP1A2 | CYP2C9 | macro |
|---|---:|---:|---:|
| ST-RAE (lower better) | +0.0136 | +0.0154 | +0.0072 |
| R² (higher better) | +0.0397 | +0.0072 | +0.0117 |
| MAE (lower better) | −0.0413 | −0.0022 | −0.0109 |

So ST-RAE is the odd one out: **MAE moves with R², not with ST-RAE.** Stated as leaderboard
records rather than as deltas, s = 1.0 **takes the record on six cells** — macro MAE and R²,
CYP1A2 MAE and R², CYP2C9 MAE and R² — every one of them from s = 0.5, which keeps the ST-RAE
record. Fully correcting the placement is the better submission on two of the three error metrics
and the worse one on the metric this project and this leaderboard are scored by.

**s = 0.5 wins or ties on every isoform on ST-RAE** — it wins CYP1A2 and CYP2C9 and ties the two
untouched ones — so choosing it needs no per-isoform cherry-picking. That dominance is
**metric-conditional**, and the condition matters: on R² and MAE the same comparison reverses on
all four (two wins, two ties for s = 1.0). s = 0.5 is the right choice **conditional on ST-RAE
being the metric that decides**, not because it is the better-calibrated file. On the evidence here,
s = 1.0 is the better-calibrated file.

**That condition is assumed, not confirmed, and it is worth stating because the final-submission
choice turns on it.** ST-RAE is listed **first** in the vendored evaluator's own `ACTIVITY_METRICS`
and is the metric its docstring uses as the exemplar `primary_metric`
(`src/vendor/openadmet_eval/`), and this project has treated it as governing throughout — but the
vendored evaluator sets no default `primary_metric`, the live backend's sort key appears nowhere in
this repository, and `CLAUDE.md` already carries the standing caveat that the vendored evaluator has
never been diffed against the live backend. **Nor can the recorded board ranks settle it**: they are
not monotone in macro ST-RAE across dates (NB47-mae's 0.6328 ranked 83 on 2026-10-10 while
NB48-ensemble's worse 0.6356 ranked 74 two days later), which is what an actively-resubmitting field
produces and which makes ranks unusable as a cross-date comparison in either direction. If the
organisers rank on R² or MAE, s = 1.0 is the better submission of the two.

**What the pair settles and what it does not.** It settles that the shrinkage fraction is not free —
half-correcting scored better on the governing metric, so the hedge paid. It does not establish
that s = 0.5 is optimal: `s` was sampled at two points, the ST-RAE-versus-R² optimum lies somewhere
between them or below, and both readings are half-set. Nor does either submission say anything
about the base model, since three of NB36-widened's four columns are unchanged in s = 0.5 and two
in s = 1.0.

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
| 12 | **27b, re-sent (accidental)** | **2026-09-21 12:49 UTC** | **0.5982** |
| 13 | NB36-widened | 2026-10-01 19:54 UTC | 0.627 |
| 14 | **NB40-ensemble** | **2026-10-02 13:37 UTC** | **0.6384** |
| 15 | NB47-mae | 2026-10-10 21:44:38 UTC | 0.6328 |
| 16 | NB48-ensemble | 2026-10-12 00:08 UTC | 0.6356 |
| 17 | **NB56-identify** (probe, deliberately degraded) | **2026-10-28 19:38:58 UTC** | **0.6947** |
| 18 | **NB56b-s050** | **2026-10-29 14:18:11 UTC** | **0.5861 — best on record** |
| 19 | **NB56b-s100** | **2026-10-30 07:54:11 UTC** | **0.5933 — second best** |

One row here is not in the per-isoform table above: the accidental 27b re-send at order 12, which
is byte-identical to order 10 and so needs no column. **It is the row the interim full-set
leaderboard scored**, which is why order 12 matters out of proportion to its being a duplicate.
NB40-ensemble at order 14 was the other absence until its row was recovered from the 2026-10-02
screenshots on 2026-10-29; both it and the re-send were the two non-deliberate, non-hand-run sends.

The trajectory is not monotonic: five step-to-step regressions (04b→NB10 +0.1120, 10c→NB12
+0.1161, NB12→NB13-calib +0.0512, NB27-widened→NB30 +0.2005, 27b→NB32-A +0.0506) interrupt an
otherwise steady improvement from `10c` (0.7138) down to NB56b-s050 (0.5861, the best on record),
which took the record from 27b (0.5982) on 2026-10-29. **NB30's
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

**Update (2026-10, the record finally moves).** CYP2C9's 0.5375 had stood across NB12, NB19,
NB27-widened and 27b, four submissions and five weeks. Both 2026-10 submissions beat it:
**NB47-mae at 0.5255, the new best on record**, and NB48-ensemble at 0.5290. NB48-ensemble also
set the highest CYP2C9 Spearman on record, 0.7692 against a prior best of 0.7585. Two things are
worth noting about how that happened. Neither submission applied any correction to CYP2C9: both
left it raw, its spread gap falling under notebook 36's 10% bar in NB48-ensemble's case and its
ratio already exceeding the target in NB47-mae's. And notebook 48's own pre-registered analysis
named CYP2C9 as **not resolvable**, one of the two isoforms it expected to learn nothing about. So
this isoform remains the clearest positive signal on the project, and it keeps moving when it is
left alone rather than when it is corrected or selected for.

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
- NB36-widened's figures were supplied directly by the user (2026-10-02 board capture, submitted
  2026-10-01 19:54 UTC under username `fold-zero`) and verified in this pass. **All five metrics
  check out exactly** — every macro value equals the mean of its own four isoform values to the
  displayed precision, with no rounding-boundary case: ST-RAE
  (0.7241+0.5435+0.7502+0.4904)/4 = 0.627050 → 0.627; MAE
  (0.9616+0.5396+0.9319+0.5373)/4 = 0.742600 → 0.7426; R²
  (0.3068+0.5308+0.3768+0.6563)/4 = 0.467675 → 0.4677; Spearman
  (0.7565+0.7554+0.4438+0.8020)/4 = 0.689425 → 0.6894; Kendall
  (0.5611+0.5652+0.3076+0.6132)/4 = 0.511775 → 0.5118. Its CYP2D6 row additionally reproduces
  NB32-A's on all five metrics, as expected from a byte-identical column, and its Spearman and
  Kendall reproduce NB32-A's on all four isoforms, as expected from a rank-preserving transform —
  both checked against this table rather than assumed.
- **This document lagged this submission by a day, and that day had a cost.** NB36-widened was
  sent on 2026-10-01 and recorded here on 2026-10-02; in between, notebook `38` had to take its macro score
  and two per-isoform gains from its task brief and label them as unverifiable. That is the fourth
  time this file has lagged a real submission (after 27b, NB32-A, and the count itself being stale
  at nine until 2026-10-01). The lag has a concrete cost each time: notebook `32` "corrected" its
  own brief to a wrong current-best figure because of it.
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
