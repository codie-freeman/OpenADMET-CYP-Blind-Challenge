# OpenADMET CYP450 Inhibition Blind Challenge — Direct Inhibition Track

Solo entry to the [OpenADMET](https://openadmet.org) CYP450 Inhibition Blind Challenge, **Direct Inhibition regression track only**. BSc Pharmaceutical Chemistry dissertation (module **CH3PRO**, University of Reading), supervised by **Dr Mauricio Cafiero**.

**Timeline:** opened 17 Aug 2026 · intermediate deadline 24 Sep 2026 · final close 3 Nov 2026.

**Task:** predict pIC50 (direct inhibition) for four CYP450 isoforms — CYP1A2, CYP2C9, CYP2D6, CYP3A4 — on a 750-compound blinded test set, from a 4,905-compound labeled training set.

## Repository structure

```
notebooks/   analysis and modelling notebooks, run in numbered order (00 → 04b)
src/         shared library code (src/features.py) + vendored third-party code (src/vendor/)
scripts/     standalone long-running training scripts, launched by hand
data/
  raw/       untouched HuggingFace release (gitignored except PROVENANCE.md)
  processed/ curated/derived artifacts from notebooks 01 and 03 (gitignored, regenerated)
  folds/     frozen 5x5 repeated-CV fold assignment (data/folds/cv_folds.csv)
outputs/     per-run predictions, params, metrics tables, submission file(s)
logs/        stdout/progress logs from scripts/ runs
models/      frozen model checkpoints from scripts/ (gitignored)
.cc-history/ Claude Code session transcripts (AI-usage disclosure trail)
```

`data/raw/` and `data/processed/` are gitignored — see `data/raw/PROVENANCE.md` for download instructions and the pinned HuggingFace revision. A fresh clone must regenerate both by running the notebooks in order (see "Reproducing this work").

## Notebook status overview

A scannable index of all 25 notebooks currently in `notebooks/`, verified against what's actually
on disk. This is a high-level index only — it sits alongside, and does not replace, the detailed
per-notebook prose below or CLAUDE.md's own notebook-status log (the authoritative source for the
full reasoning behind each entry; read it for detail this table deliberately omits).

Status categories: **Core pipeline** (data/model work the project's real submissions depend on),
**Diagnostic — negative result** (tested something, found it doesn't help or doesn't resolve the
question, documented and closed), **Diagnostic — in progress** (not yet complete or not yet acted
on), **Superseded** (replaced by a later notebook), and **Flagged — ambiguous** (doesn't cleanly
fit one bucket — see the note below the table rather than treating the label as a guess).

| Notebook | Purpose | Status |
|---|---|---|
| [`00_schema_audit.ipynb`](notebooks/00_schema_audit.ipynb) | Structural audit of raw HuggingFace files; scopes the track to `TRAIN_inhibition.csv`/`TEST-BLINDED.csv` | Core pipeline |
| [`01_data_curation.ipynb`](notebooks/01_data_curation.ipynb) | SMILES canonicalization, InChIKey generation, duplicate/leakage checks | Core pipeline |
| [`02_chemical_space_exploration.ipynb`](notebooks/02_chemical_space_exploration.ipynb) | Chemical-space and cluster-split feasibility exploration | Core pipeline |
| [`03_features_and_fold_split.ipynb`](notebooks/03_features_and_fold_split.ipynb) | Freezes the 5×5 repeated-CV fold assignment and tabular feature sets (never touched again) | Core pipeline |
| [`03b_log2fc_pretrained_encoder.ipynb`](notebooks/03b_log2fc_pretrained_encoder.ipynb) | log2fc-pretrained Chemprop encoder embeddings, reporting only | Diagnostic — in progress (not yet ablated downstream; paused pending a leakage check) |
| [`04a_baseline_screen.ipynb`](notebooks/04a_baseline_screen.ipynb) | Single-fold 21-config screen that picked the CheMeleon-init recipe | Superseded (by 05's full 5×5 CV comparison) |
| [`04b_final_submission.ipynb`](notebooks/04b_final_submission.ipynb) | First live submission — single, unensembled `chemprop_chemeleoninit` | Core pipeline (reference baseline every later submission is compared against) |
| [`04c_submission_diagnostics.ipynb`](notebooks/04c_submission_diagnostics.ipynb) | Diagnoses 04b's CYP2D6 underperformance; confirms cause via OpenADMET Discord | Diagnostic — negative result |
| [`05_cv_comparison.ipynb`](notebooks/05_cv_comparison.ipynb) | Full 5×5 CV comparison, all 11 real configs — the authoritative CV table | Core pipeline |
| [`05b_cluster_cv_comparison.ipynb`](notebooks/05b_cluster_cv_comparison.ipynb) | Cluster-aware (Butina) CV split vs. random split, `chemprop_chemeleoninit` only | Diagnostic — negative result (mixed/inconclusive, no adoption) |
| [`05c_cluster_cv_leakage_sensitivity.ipynb`](notebooks/05c_cluster_cv_leakage_sensitivity.ipynb) | Tests whether non-chemprop configs swing more under cluster CV (leakage-sensitivity check) | Diagnostic — negative result (real CYP3A4 signal found, descriptive only, no adoption) |
| [`06_outlier_check.ipynb`](notebooks/06_outlier_check.ipynb) | CYP2D6 residual/CI-width outlier exclusion | Flagged — ambiguous (see note) |
| [`07_weighting_tuning.ipynb`](notebooks/07_weighting_tuning.ipynb) | CYP2D6 sample weighting + hyperparameter tuning, tabular configs | Flagged — ambiguous (see note) |
| [`08_ensemble_selection.ipynb`](notebooks/08_ensemble_selection.ipynb) | Per-isoform ensemble selection via CV tiers, simple averaging | Superseded (by 11b's Caruana selection) |
| [`09_placement_recalibration.ipynb`](notebooks/09_placement_recalibration.ipynb) | R²-decomposition placement/recalibration diagnostic | Diagnostic — negative result |
| [`10_final_retrain_predict.ipynb`](notebooks/10_final_retrain_predict.ipynb) | Second live submission, using 08's simple-average ensembles | Flagged — ambiguous (see note) |
| [`10b_blind_regression_investigation.ipynb`](notebooks/10b_blind_regression_investigation.ipynb) | Investigates notebook 10's blind-score regression | Diagnostic — negative result (root cause not definitively isolated) |
| [`10c_control_submission.ipynb`](notebooks/10c_control_submission.ipynb) | Control submission isolating environment migration from ensembling | Core pipeline (active reference point; reused directly by later notebooks) |
| [`11_caruana_prep.ipynb`](notebooks/11_caruana_prep.ipynb) | Pooled OOF prep table across all 25 CV folds, all 11 configs | Core pipeline |
| [`11b_caruana_selection.ipynb`](notebooks/11b_caruana_selection.ipynb) | Caruana bagged ensemble selection, capped — current ensemble-selection method | Core pipeline |
| [`12_caruana_retrain_predict.ipynb`](notebooks/12_caruana_retrain_predict.ipynb) | Third live submission, using 11b's capped Caruana ensembles | Core pipeline |
| [`12b_caruana_vs_control_comparison.ipynb`](notebooks/12b_caruana_vs_control_comparison.ipynb) | Compares notebook 12 vs. 10c prediction spread/agreement | Diagnostic — negative result (capping helped but didn't fully close the spread-compression gap) |
| [`13_aid1851_blind_population_calibration.ipynb`](notebooks/13_aid1851_blind_population_calibration.ipynb) | AID 1851 external population calibration; fifth live submission | Diagnostic — negative result (worse than 04b/10c on every isoform on real blind data) |
| [`14_censored_mle_population_calibration.ipynb`](notebooks/14_censored_mle_population_calibration.ipynb) | Censored-MLE re-estimate of AID 1851 population moments | Diagnostic — negative result (clean OOF negative for every isoform) |
| [`15_analog_holdout_comparison.ipynb`](notebooks/15_analog_holdout_comparison.ipynb) | Analog-holdout split approximating the real blind-set construction | Diagnostic — negative result (gap 2.3–9.3x larger than CV gap, opposite of hoped-for outcome) |

**On the three "Flagged — ambiguous" rows:** these don't cleanly fit one status without a judgment
call the project hasn't explicitly made, so they're flagged here rather than guessed.
- **06** and **07** found real, adopted-at-the-time CV improvements (CYP2D6 outlier exclusion;
  `chemeleon__rf` tuning), but both fed only into notebook 08's ensembles, which are now
  superseded — 11b's candidate pool deliberately excludes 06's residual-excluded and 07's tuned
  configs (an isolation test, not an oversight), so neither notebook's specific output is part of
  the currently active recipe. Whether they should be relabeled "superseded" or kept as
  standalone diagnostics is an open call.
- **10** was a real, deliberate live submission (not a dry run), which argues for "core pipeline"
  the way 04b/10c/12 are treated here — but its recipe (08's simple-average ensembles) has since
  been replaced by 11b/12's Caruana selection, which argues for "superseded" the way 08 is
  treated. Kept out of both buckets rather than forced into either.

## Notebooks

Run in order; each loads only the frozen artifacts the previous ones wrote.

**[`00_schema_audit.ipynb`](notebooks/00_schema_audit.ipynb)** — Structural audit of all five raw HuggingFace files. Key finding: `TRAIN_TDI.csv`'s `direct_inhibition` columns are **byte-identical** to `TRAIN_inhibition.csv` — a copy, not an independent measurement — so `TRAIN_TDI.csv`, `TRAIN_Emax.csv`, and `single-concentration-TRAIN.csv` are set aside from this track, leaving `TRAIN_inhibition.csv` (4,905 rows) and `TEST-BLINDED.csv` (750 rows) in scope. Also resolves the source of `_conf_high`/`_conf_low`/`_std`, confirmed via OpenADMET's scoring code and a Discord confirmation from OpenADMET staff.

**[`01_data_curation.ipynb`](notebooks/01_data_curation.ipynb)** — Canonicalizes SMILES, generates InChIKeys (`src/features.py`), checks salts, parse failures, duplicates, and train/test leakage (all by InChIKey). All checks clean: 0 salts, 0 parse failures, 0 duplicates, 0 leakage. Writes `train_inhibition_curated.csv` and `test_blinded_curated.csv`.

**[`02_chemical_space_exploration.ipynb`](notebooks/02_chemical_space_exploration.ipynb)** — Feasibility-checks a planned cluster-based CV split. Four findings: (1) low pairwise redundancy (mean Tanimoto ~0.13–0.14 per isoform); (2) the blinded test set is *closer* to train (median max-Tanimoto 0.587) than train is to itself (0.450, Mann-Whitney p≈1.3×10⁻²²⁵) — the opposite of "held out to be maximally novel"; (3) SALI activity-cliff prevalence is undersampled at the literature threshold (Tanimoto≥0.9: 0–4 pairs per isoform) but well above the PXR reference range (4–12%) at better-sampled thresholds (14–45%); (4) 96.4% of Bemis-Murcko scaffolds are singletons, covering 89.0% of compounds. Conclusion: a scaffold-exact cluster split would functionally resemble random splitting for most of the data; whichever split is used should measure its actual achieved shift rather than assume one. Also notes CYP2D6 has the lowest per-isoform NN similarity (mean 0.375 vs. 0.403–0.488 for the others) and didn't seed the real test set's construction the way the other three isoforms did — both later relevant to CYP2D6's leaderboard result (see "Initial baseline submission").

**[`03_features_and_fold_split.ipynb`](notebooks/03_features_and_fold_split.ipynb)** — Freezes two artifacts everything downstream depends on: a 5×5 repeated random CV fold assignment (`data/folds/cv_folds.csv`, `RepeatedKFold`, seed 42, verified as 5 genuinely distinct shuffles) and three tabular feature sets (ECFP4+9 physicochemical descriptors, full RDKit 2D descriptors, CheMeleon embeddings). Rejects CYP2D6-specific CV stratification — its per-fold representation is already balanced under plain random splitting (4.3% relative std, in line with the other isoforms). Reconciles the single-concentration "~1,500 hit compounds" figure against a directly computed 3,326 (a 2.2x discrepancy, reported not forced); finds the single-concentration population is almost entirely already pIC50-labeled, ruling out "restrict log2fc pretraining to unlabeled compounds" for notebook 03b.

**[`03b_log2fc_pretrained_encoder.ipynb`](notebooks/03b_log2fc_pretrained_encoder.ipynb)** — Reporting notebook; training happens in `scripts/train_log2fc_encoder.py` (a random-init Chemprop encoder trained on the 4,376-compound log2fc primary screen, then frozen as a feature extractor — not CheMeleon continued-pretraining). Loads that script's output, confirms shapes/NaN-freedom, and plots the training loss (still decreasing, not plateaued, at epoch 50; 639s/10.7min total). No validation split — a deliberate leakage-avoidance choice, since the encoder trains on data overlapping almost completely with the pIC50-labeled set. One candidate feature option among four; usefulness is notebook 04's job.

**[`04a_baseline_screen.ipynb`](notebooks/04a_baseline_screen.ipynb)** — Single-fold screen across 21 configs: naive mean/median floor, 9 tabular configs (3 feature sets × RF/XGBoost/LightGBM), 10 Chemprop configs (4 single-task + 1 multitask, random-init and CheMeleon-init variants). Scored with the vendored OpenADMET bootstrap evaluation code. TabICLv2 dropped — a single `.predict()` call tried to allocate 36–50GB against this machine's 36GB RAM. Result: **CheMeleon-init multitask Chemprop wins clearly on every macro metric** (0.700 ST-RAE vs. 0.774 next-best) — see "Initial baseline submission" below.

**[`04b_final_submission.ipynb`](notebooks/04b_final_submission.ipynb)** — Validates and (on manual, deliberate action) submits the activity-track prediction file. The model is trained by `scripts/train_final_submission_multitask.py`, retraining the screened CheMeleon-init multitask recipe on all 4,905 labeled compounds (15% internal validation slice for early stopping only, independent of `cv_folds.csv`), predicting on the real 750-compound blind test set. Format validated against the tutorial repo's own `validate_activity_submission` (vendored). Training: 21 epochs (patience=5), ≈297 min (~5hrs), CPU-only.

## Shared infrastructure

**`src/features.py`** — the one shared module all fingerprint/similarity/split logic goes through:
- `canonicalize_smiles` / `smiles_to_inchikey` / `add_canonical_smiles_and_inchikey` — SMILES canonicalization and InChIKey generation.
- `ecfp4_fingerprints`, `pairwise_tanimoto_matrix`, `nearest_neighbor_similarity` — ECFP4 (radius 2, chirality off by default) and Tanimoto similarity.
- `compute_sali` — SALI activity-cliff scoring (Guha & Van Drie 2008), per-isoform only, never unioned.
- `isoform_structural_descriptors` — a narrow 9-descriptor set (Kiani & Jabeen 2019), including two labeled proxies (`logd_proxy`, `vsa_acc_proxy`) where RDKit has no native equivalent.
- `rdkit_2d_descriptors` — full RDKit 2D set, matching the official OpenADMET baseline recipe.
- `bemis_murcko_scaffold` — scaffold extraction (notebook 02).
- `chemeleon_embeddings` / `frozen_encoder_embeddings` — frozen forward-pass feature extraction, shared batching loop reused by `scripts/train_log2fc_encoder.py`.
- `assign_screen_split` / `assign_final_submission_split` — the two non-CV split functions, called identically from notebooks and scripts for a byte-identical split without freezing an intermediate file.

`src/scoring.py` and `src/chemprop_screen.py` hold project code built *on top of* the vendored evaluation code (a macro ST-RAE helper, and shared Chemprop training/logging orchestration) — the vendored files themselves are never edited.

**`src/vendor/`** — unmodified third-party code, each with its own `PROVENANCE.md`:
- **`openadmet_eval/`** — OpenADMET's official evaluation code, pulled byte-for-byte from `OpenADMET/CYP-Challenge-Tutorial` and verified against GitHub's blob SHAs. Used for every bootstrapped metric table (ST-RAE, MAE, R², Spearman's ρ, Kendall's τ) so scores match the real leaderboard exactly. Caveat: `utils.bootstrap_sampling` uses a single fixed, globally-shared seed under `@lru_cache` — fine for cross-submission comparability, but must be wrapped before use inside a per-fold CV loop.
- **`validation/`** — the tutorial's own `validate_activity_submission` (plus `tdi_validation.py`, vendored only because `__init__.py` imports it), used unmodified before submission.

## Initial baseline submission

**`chemprop_multitask_chemeleoninit`** — a single multitask Chemprop model (all four isoforms, one shared encoder) with a CheMeleon-initialized message-passing encoder — beat every alternative on every macro metric in the single-fold screen:

| Config | ST-RAE | MAE | R² | Spearman's ρ |
|---|---|---|---|---|
| **chemprop_multitask_chemeleoninit** | **0.700** | **0.574** | **0.352** | **0.583** |
| tabular_baseline (ECFP4+descriptors) / RF | 0.774 | 0.603 | 0.292 | 0.518 |
| tabular_baseline / LightGBM | 0.783 | 0.609 | 0.278 | 0.495 |
| tabular_rdkit2d (full RDKit 2D) / RF | 0.787 | 0.613 | 0.271 | 0.516 |
| chemeleon embeddings / RF | 0.793 | 0.625 | 0.262 | 0.527 |
| chemprop_multitask_randominit | 0.926 | 0.708 | 0.082 | 0.403 |
| naive median | 0.988 | 0.738 | −0.025 | 0.000 |
| naive mean | 1.007 | 0.745 | −0.007 | 0.000 |

It beat its own random-init counterpart by the largest margin of any of the five paired Chemprop comparisons (0.700 vs. 0.926 — the four single-task pairs range only 0.003–0.084), evidence CheMeleon's pretrained weights are doing real work, not just the multitask architecture — at a real cost (28–83x slower per model than random-init on this hardware). TabICLv2 was excluded from this screen entirely (hardware infeasibility, see 04a).

This config was retrained by `scripts/train_final_submission_multitask.py` on the full 4,905-compound labeled set to produce `outputs/submissions/activity_submission_v1.csv`.

### Leaderboard result

Submitted under the alias **`fold-zero`** (repo is currently private; revisit the alias if it's ever made public):

| Track | Rank | ST-RAE | MAE | R² | Spearman's ρ | Kendall's τ |
|---|---|---|---|---|---|---|
| Overall (macro) | 21 | 0.7179 | 0.8952 | 0.2191 | 0.6739 | 0.5008 |
| CYP1A2 | 22 | 0.7114 | 0.9589 | 0.2824 | 0.7430 | 0.5473 |
| CYP2C9 | 18 | 0.5408 | 0.5399 | 0.5263 | 0.7428 | 0.5490 |
| CYP2D6 | 23 | 1.1903 | 1.5842 | −0.6220 | 0.3769 | 0.2606 |
| CYP3A4 | 3 | 0.4291 | 0.4976 | 0.6898 | 0.8329 | 0.6463 |

One trained model from one training run — not the 5×5 repeated-CV comparison the baseline screen was built to precede.

## Final retrain submission (notebook 10)

`notebooks/10_final_retrain_predict.ipynb` (backed by `scripts/10_final_retrain_predict.py`) retrained every distinct model config referenced by notebook 08's per-isoform winning ensemble/single-best on the full 4,905-compound labeled set and combined them per isoform (CYP1A2: 4-way average, CYP2D6/CYP3A4: 2-way average, CYP2C9: single model, unaveraged) — intended to replace `04b`'s result as the live leaderboard entry.

### Leaderboard result

Submitted under the alias **`fold-zero`**, 2026-09-03 15:26 UTC:

| Isoform | Rank | ST-RAE | MAE | R² | Spearman's ρ | Kendall's τ |
|---|---|---|---|---|---|---|
| Overall (macro) | 70 | 0.8299 | 0.9690 | 0.1054 | 0.6748 | 0.4977 |
| CYP1A2 | 76 | 0.7850 | 1.0214 | 0.2051 | 0.7246 | 0.5297 |
| CYP2C9 | 56 | 0.5489 | 0.5444 | 0.5272 | 0.7389 | 0.5462 |
| CYP2D6 | 87 | 1.4603 | 1.7489 | −0.9392 | 0.4401 | 0.3081 |
| CYP3A4 | 46 | 0.5253 | 0.5614 | 0.6285 | 0.7955 | 0.6065 |

Scored substantially worse than `04b` on every isoform, on both ST-RAE and R². A read-only pipeline audit (6 checks: training volume, training-pool lineage, merge/row-alignment, prediction-scale sanity, individual-compound tracing, cross-script comparison against `04b`) found no pipeline bug — see `notebooks/10b_blind_regression_investigation.ipynb` for the follow-up investigation into *why* the decline happened (a spread-compression hypothesis, tested against both R² and the project's actual governing metric, ST-RAE; verdict: partially supported at best, not a clean or sufficient explanation).

## Reproducing this work

**Environment:** conda env `cyp-admet-v2`, Python 3.11, Apple Silicon (native arm64, no Rosetta), built from [`environment.yml`](environment.yml) at the repo root via `conda env create -f environment.yml`. The original `cyp-admet` environment (x86_64/Rosetta) still exists on this machine as a fallback but is no longer the primary environment. A secondary `pip freeze` reference is provided at [`requirements.txt`](requirements.txt) — not an independent install path (some packages, RDKit in particular, were installed via specific channels/pip index ordering that a plain `pip install -r requirements.txt` into an empty environment may not replicate correctly); `environment.yml` is the authoritative spec.

**Hardware constraint:** RDKit and PyTorch each bundle their own OpenMP runtime, and running them multi-threaded in the same process segfaults on this machine — `OMP_NUM_THREADS=1` and `KMP_DUPLICATE_LIB_OK=TRUE` must be set before RDKit/PyTorch/Chemprop are first imported (every notebook/script sets this at the top). Chemprop/CheMeleon training remains CPU-only regardless of environment, now for a precisely identified reason: MPS's `scatter_reduce` operation has no deterministic kernel, and Chemprop enables PyTorch Lightning's `deterministic=True` mode whenever a pytorch seed is set — `chemprop/cli/train.py` (v2.3.1), lines 1798–1803: `deterministic = True` if `args.pytorch_seed is not None` else `False`, passed straight through to `pl.Trainer(..., deterministic=deterministic)` at line 1938 — and this project's convention is to always set a seed. This is a different limitation than previously documented here: the op is no longer literally unimplemented on MPS, it's implemented but non-deterministic, and Chemprop's CLI exposes no way to relax that check. Other computation — tabular models (RF/XGBoost/LightGBM), RDKit feature generation, general torch use — now runs natively on Apple Silicon under `cyp-admet-v2` rather than under Rosetta emulation: measured directly on an identical small task (RandomForestRegressor fit, 500-row subset of `tabular_rdkit2d` features, CYP1A2 endpoint, 3 runs each), `cyp-admet-v2` completed in a median 1.12s vs. `cyp-admet`'s median 1.58s — about 1.4x faster on this task (a single small-fit data point, not a full-pipeline projection).

**Run order:**
1. `notebooks/00_schema_audit.ipynb` → `01_data_curation.ipynb` → `02_chemical_space_exploration.ipynb` → `03_features_and_fold_split.ipynb`, in order (each loads the previous one's frozen output).
2. `scripts/train_log2fc_encoder.py`, launched by hand in the background (**not** from a notebook or agentic session — this project's convention for any long-running training):