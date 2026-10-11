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

## Quantum-chemistry descriptors (notebooks `41`–`46`)

- **`run_qm_descriptors.py`** — the production QM run: all 5,655 compounds through notebook `42`'s
  measured pipeline (RDKit ETKDGv3+MMFF conformer → GFN2-xTB optimisation → `B3LYP D3 def2-SV(P)
  NORI Hirshfeld` single point), 11,310 ORCA jobs, **zero failures**. 281.6 core-hours in 47 h of
  wall time across four stop-start sessions at ~100% parallel efficiency on six single-core
  workers. Restartable per job on a content-based done-check (ORCA exits 0 on error, so the exit
  code is trusted nowhere); SIGTERM drains cleanly. Writes `outputs/qm_descriptors/`. See that
  directory's own README for the three corrections it makes to notebook `42`'s projections.
- **`qm_supervisor.sh`** — a polling supervisor that relaunches `run_qm_descriptors.py` if it dies
  or hits an `--hours` cap, stops itself on completion, and caps relaunches so a crash-loop cannot
  churn. Off-switch is the sentinel file `outputs/qm_descriptors/KEEP_RUNNING`. Added beyond the
  original brief and tested by deliberately stopping a healthy run rather than trusted on first
  use.
- **`parse_orca_output.py`** — extracts per-atom, per-bond and molecule-level quantities from the
  completed run into three tidy tables (`outputs/qm_parsed/`). Prefers each job's `.property.txt`
  over its `.out` wherever a quantity exists in both — typed indexed arrays rather than
  fixed-width text — which confines the one real parser trap to a single block; that block is then
  read by fixed character position, because a two-character element symbol runs into the shell
  label (`7 Cls       :`) and whitespace splitting silently drops those atoms on 815 compounds.
  Eight hard assertions per compound, including an atom-order check against RDKit and a
  bond-subset check that catches same-element permutations the element check cannot see. Parses
  the full set in well under a minute; restartable per compound via sharded writes.

### The QM modelling arms (notebooks `44`–`46`)

- **`44_build_qm_feature_blocks.py`**, **`44_redundancy_triage.py`**, **`44_run_qm_arms.py`**,
  **`44_posthoc_group_c_probe.py`** — notebook `44`'s molecule-level arm: feature blocks and its
  pre-registration, the cheap-recoverability triage, the training driver (tier-1 model screen plus
  the 25-fold confirmation), and the post-hoc probe that added the excluded group-(c) columns back.
- **`46_build_feature_blocks.py`** — notebook `46`'s six feature blocks and its pre-registration,
  written **before any training**: the script reads no fold file and no prediction, so the criteria
  cannot be tuned to a result. Its baseline is notebook `45`'s own axis A (rdkit2d's 217 plus
  ecfp4_narrow's 9), which is the block the descriptors were screened against — deliberately not
  notebook `44`'s `RDKIT2D`, which is the 217 alone, and the reason `44`'s arms could not simply be
  reused despite sharing the partition. The halogen block ships its element-conditioned columns as
  NaN with presence flags rather than pre-filled, because the fill has to be fitted per fold.
- **`46_run_arms.py`** — notebook `46`'s training driver. Self-contained by choice rather than
  importing `44_run_qm_arms.py`: every path that script writes to is a module-level constant inside
  `outputs/44_qm_molecular_block/`, which is protected, and monkeypatching those globals risks
  missing one. Functions are copied verbatim with attribution, following the precedent by which
  notebooks 35 and 37 each copied 34's harness. Two things are new — the loader tolerates NaN for
  blocks the pre-registration declares may carry it, and the per-fold median fill for those columns
  is fitted on that fold's training portion only.

### The per-atom QM arm (notebook `52`)

- **`52_build_atom_ordering.py`** — builds and verifies the per-compound ORCA-to-chemprop atom
  permutation, which is the whole risk in notebook `52`. The QM run featurised every compound as
  `Chem.AddHs(Chem.MolFromSmiles(row.SMILES))` — the **raw** `SMILES` column — while chemprop is
  given `-s canonical_smiles`, and `canonical_smiles` reproduces ORCA's element sequence on only
  4,212 of 5,655 compounds. Both molecules are routed through one `Chem.MolToSmiles` call and
  matched on `_smilesAtomOutputOrder`, so RDKit's own tie-breaking resolves symmetry identically on
  both sides. Every check is a hard failure rather than a warning, including notebook `43`'s
  element-sequence and bond-subset assertions, the same element check again against
  `qm_atoms.csv` itself, and a **negative control** that corrupts a verified permutation and
  requires the checks to reject it — without which the assertions could be vacuous.
- **`52_atom_arrays.py`** — the `.npz` per-atom arrays and the widened CheMeleon encoder. Two
  things here are corrections rather than preferences, both documented in the module docstring.
  chemprop 2.3.1 **refuses** `--atom-features-path` when CheMeleon is named, because its published
  `W_i` is `(2048, 86)` against a 72-wide V2 featurizer; the refusal is guarded by
  `if not local_foundation:`, so this builder writes a local foundation file carrying CheMeleon's
  own weights widened by *k* zero columns at index `d_v`, asserted numerically identical to
  CheMeleon at initialisation. And on the `--from-foundation` path chemprop computes an
  atom-feature scaler, logs it, and never attaches the `GraphTransform` to the loaded encoder,
  while `normalize_inputs` mutates the training arrays in place — so training would be scaled and
  `predict` unscaled. The scaling is therefore done here, from each fold's `screen_inner_train`
  atoms only, with `--no-atom-feature-scaling` passed.
- **`52_run_qm_atom_arms.py`** — notebook `52`'s training driver: `--sanity` (reproduce notebook
  `05`'s stored `chemprop_chemeleoninit__repeat0_fold0`), `--wiring` (two fits differing only in
  the `.npz` contents, real against all-zeros, as the gate on whether the CV is worth running),
  `--confirm` (the frozen 25 cluster-disjoint folds) and `--final` (full-data fit). The harness is
  notebook `44`'s, itself notebook `37`'s, with the differences the arm forces.

All scripts import `src/features.py`, `src/scoring.py`, and (where relevant) `src/chemprop_screen.py`
or `src/cv_bootstrap.py` rather than reimplementing any of that logic locally — see `src/README.md`.
`parse_orca_output.py` additionally reuses `run_qm_descriptors.py`'s own `read_orca_xyz()` rather
than writing a second reader for the same file.
