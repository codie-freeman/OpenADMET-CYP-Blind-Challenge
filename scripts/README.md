# scripts/

Standalone scripts, launched by hand from a terminal (`caffeinate -i nohup python scripts/... &`),
never from a notebook or an agentic session — this project's convention for anything that is
long-running (minutes to hours of real training) or needs to keep running unattended in the
background. Notebooks read these scripts' saved output from disk; they never invoke them directly.

## Baseline screen (notebook `04a`)

- **`run_baseline_screen_chemprop.py`** — trains the 5 random-init Chemprop configs for `04a`'s
  single-fold screen (4 single-task + 1 multitask).
- **`run_baseline_screen_chemprop_chemeleon.py`** — trains the CheMeleon-initialized counterpart of
  the same 5 configs, alongside (not replacing) the random-init run, for direct comparison.
- **`train_log2fc_encoder.py`** — trains a fresh, random-init Chemprop encoder on the log2fc
  primary-screen data (notebook `03b`), freezes it, and extracts embeddings.
- **`train_log2fc_encoder_blind.py`** — retrains that same encoder with `04a`'s screen-test
  compounds excluded from pretraining, as a leakage check on `04a`'s suspiciously strong
  `log2fc_pretrained_*` tabular results. Exists but has not been run to completion.

## Final submission (notebook `04b`)

- **`train_final_submission_multitask.py`** — retrains `04a`'s winning config
  (`chemprop_multitask_chemeleoninit`) on the full 4,905-compound labeled set and predicts on the
  real blind test set. Long-running (~5 hours, CPU-only).

## Full CV comparison (notebook `05` and its variants)

- **`generate_5x5_cv_manifest.py`** — writes the static, never-mutated manifest
  (`outputs/05_cv_comparison/manifest.csv`) of all (config, repeat, fold) rows the 5×5 comparison
  needs to fill in.
- **`run_5x5_cv_comparison.py`** — the resumable driver that fills in that manifest: for each row
  not already done (checked purely by whether its own output files exist on disk), trains,
  predicts, and scores it. Multiple instances can run concurrently against the same manifest with
  no shared mutable state. Long-running (hours).
- **`generate_mordred_pca_features.py`** — computes Mordred 2D descriptors + PCA reduction for one
  of `05`'s tabular feature arms, replacing full RDKit2D there.
- **`05b_run_cluster_cv.py`** — notebook `05b`'s driver: refits `chemprop_chemeleoninit` only, over
  the 25 fold-partitions of the cluster-aware split built by the notebook itself. Mirrors
  `run_5x5_cv_comparison.py`'s own seeding/resumability conventions.
- **`05c_run_cluster_cv_randominit.py`** — notebook `05c`'s driver: the same pattern as
  `05b_run_cluster_cv.py`, for `chemprop_randominit` only (the one config in `05c`'s 6-config check
  too slow to fit directly in the notebook).

## CYP2D6 investigation (notebooks `06`–`07`)

- **`cyp2d6_outlier_check.py`** — computes notebook `06`'s two outlier-flagging criteria
  (CV-residual and CI-width) from already-computed `05` artifacts; no new training.
- **`cyp2d6_weighting_tuning.py`** — notebook `07`'s CYP2D6 sample-weighting and
  hyperparameter-tuning runs, built on `06`'s residual-excluded training pool.

## Ensembling and submissions (notebooks `10`, `12`)

- **`10_final_retrain_predict.py`** — retrains every distinct config referenced by notebook `08`'s
  per-isoform ensemble winners on the full training set and predicts on the blind set (second live
  submission, NB10).
- **`12_caruana_retrain_predict.py`** — retrains every config with non-zero Caruana weight
  (`11b`'s selection) that isn't already reusable from `10`'s output, and predicts on the blind set
  (third live submission, NB12).

## External calibration data (notebook `13`)

- **`13_fetch_aid1851.py`** — fetches PubChem AID 1851 (NCATS qHTS cytochrome panel) in full,
  batched under PubChem's 10,000-SID PUG REST cap, resumable via on-disk per-batch CSV caching.
  Used strictly as an external population-moment estimate, never as training/fine-tuning data.

## Analog-holdout split (notebook `15`)

- **`15_run_analog_holdout.py`** — trains `chemprop_chemeleoninit` once on the single train/test
  partition notebook `15` builds (`data/folds/analog_holdout_split.csv`), structurally like one
  fold of `05b_run_cluster_cv.py`.

All scripts import `src/features.py`, `src/scoring.py`, and (where relevant) `src/chemprop_screen.py`
or `src/cv_bootstrap.py` rather than reimplementing any of that logic locally — see `src/README.md`.
