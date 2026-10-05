# outputs/

Per-run predictions, params, metrics tables, and submission files. Almost every subdirectory is
named after the notebook that wrote it — read that notebook (or its entry in `notebooks/README.md`
/ `CLAUDE.md`) for what the numbers mean. This file is an index only; it does not enumerate
individual files.

| Directory | Written by | Holds |
|---|---|---|
| `04a_baseline_screen/` | `04a` | Single-fold 21-config screen results. |
| `05_cv_comparison/` | `05` | The authoritative 5×5 CV comparison — manifest, per-fold predictions/scores, summary table. Canonical source for any CV figure in this project. |
| `05b_cluster_cv_comparison/` | `05b` | Cluster-aware CV split (`cluster_cv_folds.csv`) and its `chemprop_chemeleoninit` results. |
| `05c_cluster_cv_leakage_sensitivity/` | `05c` | Leakage-sensitivity check results, 6 configs. |
| `06_outlier_check/` | `06` | CYP2D6 outlier-flagging criteria and flagged-compound lists. |
| `07_weighting_tuning/` | `07` | CYP2D6 weighting/tuning results. |
| `08_ensemble_selection/` | `08` | Superseded per-isoform simple-average ensemble selections. |
| `09_placement_recalibration/` | `09` | Own-population placement/recalibration diagnostic (negative result). |
| `10_final_retrain_predict/` | `10` (+`scripts/10_final_retrain_predict.py`) | Second live submission's training/prediction artifacts. |
| `10b_blind_regression_investigation/` | `10b` | Investigation into `10`'s regression. |
| `10c_control_submission/` | `10c` | Control submission's artifacts (env-migration isolation). |
| `11_caruana_prep/` | `11` | Pooled OOF prep tables feeding Caruana selection. |
| `11b_caruana_selection/` | `11b` | Caruana ensemble weights (capped, canonical, and uncapped), recalibration params, correlation diagnostics. |
| `12_caruana_retrain_predict/` | `12` (+`scripts/12_caruana_retrain_predict.py`) | Third live submission's training/prediction artifacts, including the CYP2C9 column notebook `19` later reuses. |
| `12b_caruana_vs_control_comparison/` | `12b` | NB12-vs-`10c` spread/agreement comparison. |
| `13_aid1851_blind_population_calibration/` | `13` | Fifth live submission (NB13-calib) build artifacts. |
| `15_analog_holdout_comparison/` | `15` (+`scripts/15_run_analog_holdout.py`) | Analog-holdout split results. |
| `16_board_solved_calibration_figures/` | `16` | Figures for the board-metrics-solved population estimate; the submission candidate itself is in `board_solved_calibration/`. |
| `17_protonation_charge_check/` | `17` | Charge/logP features, decision criteria, partial-correlation verdict. |
| `18_charge_augmented_chemprop_screen/` | `18` | Charge-augmented Chemprop screen results and decorrelation summary. |
| `19_cyp2c9_revert/` | `19` | Seventh live submission (NB19) build artifacts — see `docs/leaderboard_submissions.md` for the current best. |
| `20_residual_exclusion_screen/` | `20` | Residual-exclusion screen with random-masking control. |
| `21_strae_offset_curve/` | `21` | ST-RAE-optimal offset sweep and figures. |
| `22_auxiliary_heads_screen/` | `22` | TDI/Emax auxiliary-head screen results, prereg, figures. |
| `23_compound_pool_overlap/` | `23` | Compound-pool overlap counts and figures. |
| `24_nonnegative_stacking/` | `24` | Non-negative stacker coefficients and combiner scores, with the honest nested outer-loop protocol. |
| `25_hard_compound_analysis/` | `25` | Per-compound difficulty table — which compounds every config gets wrong, and what they have in common. |
| `26_spread_sweep/` | `26` | Pure spread-multiplier sweeps and the joint offset x spread grids. |
| `27_calibration/` | `27` | **The standing calibration notebook's outputs.** Spread-ratio table, widening multipliers, and NB27-widened's own candidate (submitted, macro 0.6155). |
| `27b_aid_cyp2d6_corrected/` | `27b` | The AID model's CYP2D6 column placement-corrected and widened — **the best submission on record (macro 0.5982)**. |
| `28_external_data/` | `28` | AID 1851 auxiliary-head single-fold screen: gate criteria, run scores, figures. |
| `29_cv_confirmation/` | `29` | The 5x5 CV confirmation of `28` — 50 real runs, paired tests, per-fold predictions. |
| `30_aid_full_retrain/` | `30` | NB30's full-data retrain artifacts and its pre-board falsifiable predictions. |
| `31_deadzone/` | `31` | Dead-zone target screen, including the post-hoc MAE_ONLY decomposition arm. |
| `32_deadzone_aid_retrain/` | `32` | NB32-A's build artifacts (submitted) and candidate B (never sent). |
| `33_butina_split/` | `33` | First cluster-disjoint partition and its 5-fold diagnostic comparison. |
| `34_butina_5x5/` | `34` | The repeated cluster-disjoint design's fold file, run scores, and OOF tables. |
| `35_butina_repeat3/` | `35` | Repeat 2 plus the MAE_ONLY arm; 15-fold analysis. |
| `36_spread_corrected_nb32a/` | `36` | NB36-widened's build artifacts (submitted, macro 0.627). |
| `37_butina_5x5_complete/` | `37` | The completed 25-sample design, Ash quote verification, and `05`-protocol figures. |
| `38_cyp3a4_placement_diagnostic/` | `38` | CYP3A4 R2 decomposition and the board-solved population fit. |
| `39_ambiguity_and_shared_blend/` | `39` | Krogh & Vedelsby ambiguity tables and the shared-blend cost comparison. |
| `40_single_object_ensemble/` | `40` | The shared-weight blend's full-data artifacts and submission candidates. |
| `41_dataset_audit_3d/` | `41` | Per-compound 3D-readiness audit — size, elements, stereochemistry, basic nitrogen. Reused widely by `42`/`43`. |
| `42_orca_timing/` | `42` | ORCA timing measurements, cost fits, and core-scaling curves. The raw job tree is gitignored. |
| `43_coordinating_atom/` | `43` | Candidate SMARTS coverage/overlap/tie-break tables, molecule panels, and the three open decisions. |
| `board_solved_calibration/` | `16` | Sixth live submission (NB16) candidate, board-metrics-solved recentring. |
| `calibration_blind_population/` | `13`, `14` | External-population calibration params (floored-mean and censored-MLE). |
| `final_submission/` | `04b` (`scripts/train_final_submission_multitask.py`) | First live submission's training/prediction artifacts. |
| `mordred_pca/` | `scripts/generate_mordred_pca_features.py` | Mordred descriptor NaN report and PCA explained-variance curve. |
| `population_moments/` | `13`, `14`, `16` | Per-isoform population-moment estimates (AID 1851 floored/censored-MLE, board-solved). |
| `qm_descriptors/` | `scripts/run_qm_descriptors.py` | The production QM run: per-compound descriptors and Hirshfeld charges for all 5,655 compounds. Has its own README. The 8.8 GB per-compound ORCA job tree is gitignored. |
| `qm_parsed/` | `scripts/parse_orca_output.py` | Per-atom, per-bond and molecule-level quantities parsed out of that ORCA tree. Has its own README. |
| `submissions/` | `04b` | `activity_submission_v1.csv`, the first submitted file. |

Notebooks `03`, `03b`, `04c`, and `04a`'s two supporting screen scripts write into shared or
already-listed directories above rather than a directory of their own; see `notebooks/README.md`
for what each notebook produces.
