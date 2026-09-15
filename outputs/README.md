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
| `17_protonation_charge_check/` | `17` | Charge/logP features, decision criteria, partial-correlation verdict. |
| `18_charge_augmented_chemprop_screen/` | `18` | Charge-augmented Chemprop screen results and decorrelation summary. |
| `19_cyp2c9_revert/` | `19` | Seventh live submission (NB19, current best) build artifacts. |
| `20_residual_exclusion_screen/` | `20` | Residual-exclusion screen with random-masking control. |
| `21_strae_offset_curve/` | `21` | ST-RAE-optimal offset sweep and figures. |
| `22_auxiliary_heads_screen/` | `22` | TDI/Emax auxiliary-head screen results, prereg, figures. |
| `23_compound_pool_overlap/` | `23` | Compound-pool overlap counts and figures. |
| `board_solved_calibration/` | `16` | Sixth live submission (NB16) candidate, board-metrics-solved recentring. |
| `calibration_blind_population/` | `13`, `14` | External-population calibration params (floored-mean and censored-MLE). |
| `final_submission/` | `04b` (`scripts/train_final_submission_multitask.py`) | First live submission's training/prediction artifacts. |
| `mordred_pca/` | `scripts/generate_mordred_pca_features.py` | Mordred descriptor NaN report and PCA explained-variance curve. |
| `population_moments/` | `13`, `14`, `16` | Per-isoform population-moment estimates (AID 1851 floored/censored-MLE, board-solved). |
| `submissions/` | `04b` | `activity_submission_v1.csv`, the first submitted file. |

Notebooks `03`, `03b`, `04c`, and `04a`'s two supporting screen scripts write into shared or
already-listed directories above rather than a directory of their own; see `notebooks/README.md`
for what each notebook produces.
