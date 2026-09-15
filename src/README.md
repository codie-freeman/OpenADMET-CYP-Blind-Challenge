# src/

Shared library code. Per `CLAUDE.md`'s rule that all fingerprint/similarity/split logic goes
through one shared module, no notebook or script recomputes these operations independently.

## `features.py`

The one shared module for canonicalization, fingerprints, similarity, and splitting:

- `canonicalize_smiles`, `smiles_to_inchikey`, `add_canonical_smiles_and_inchikey` — SMILES
  canonicalization and InChIKey generation. Used by notebook `01` (curation) and every later
  notebook that needs to identify a compound by InChIKey (leakage checks, overlap counting in
  `22`/`23`, etc.).
- `ecfp4_fingerprints`, `pairwise_tanimoto_matrix`, `nearest_neighbor_similarity` — ECFP4 (radius
  2, chirality off by default) and Tanimoto similarity. Used by notebook `02` (chemical space),
  `04c`, `05b`/`05c` (cluster splitting), and `15` (analog-holdout construction).
- `compute_sali` — SALI activity-cliff scoring (Guha & Van Drie 2008), per-isoform only, never
  unioned. Used by notebook `02`.
- `isoform_structural_descriptors` — a narrow 9-descriptor set (Kiani & Jabeen 2019). Used by
  notebook `04a`'s tabular baseline screen.
- `rdkit_2d_descriptors` — full RDKit 2D descriptor set, matching the official OpenADMET baseline
  recipe. Used by notebook `03`.
- `mordred_2d_descriptors` — Mordred 2D descriptors, feeding `scripts/generate_mordred_pca_features.py`.
- `chemeleon_embeddings`, `frozen_encoder_embeddings` — frozen forward-pass feature extraction,
  shared batching loop. Used by notebook `03` and `scripts/train_log2fc_encoder.py`.
- `butina_clusters` — Butina clustering, added for notebook `05b`'s cluster-aware CV split.
- `bemis_murcko_scaffold` — scaffold extraction, used by notebook `02`.
- `assign_screen_split`, `assign_final_submission_split` — the two non-CV split functions, called
  identically from notebooks and scripts for a byte-identical split without freezing an
  intermediate file. Used by `04a`'s screen and the final-submission training scripts.

## `scoring.py`

Project code built on top of the vendored evaluation code — a macro ST-RAE helper
(`get_ma_st_rae`). The vendored files themselves (`src/vendor/`) are never edited; anything
project-specific goes here instead. Used by every notebook and script that scores a model
(`05`–`23`, all training scripts).

## `calibration.py`

Calibration of out-of-fold predictions onto an *external*, independent population's mean/std
(`fit_blind_population_calibration`, `apply_calibration`) — distinct from the *own*-population
placement corrections tried in `09`/`11b`. Built for, and used by, notebooks `13` and `14` (AID
1851 external-population calibration). Not used by `16` (which solves for the blind population's
own moments directly from board metrics, a different mechanism entirely).

## `chemprop_screen.py`

Shared Chemprop CLI orchestration — training-CSV/predict-CSV builders, subprocess launch and
streamed progress logging, prediction verification. Factored out so the per-epoch logging
mechanics aren't duplicated across scripts. Used by `scripts/run_baseline_screen_chemprop.py` and
`scripts/run_baseline_screen_chemprop_chemeleon.py` (notebook `04a`), and by notebooks `05c`
(`05c_run_cluster_cv_randominit.py`) and `20` for their own Chemprop screens. Notebook `18`
extended this pattern locally (an extra-descriptor-column hook) rather than modifying this module,
since the module's existing builders have no hook for an extra column.

## `cv_bootstrap.py`

Per-fold bootstrap-seed wrapper (`per_fold_bootstrap_seed`) for the 5×5 CV comparison
(`scripts/run_5x5_cv_comparison.py`). Exists because the vendored evaluator's own
`bootstrap_sampling` uses a single fixed, globally-shared seed under `@lru_cache` — fine for
cross-submission comparability, but wrong to call directly inside a per-fold CV loop. Per
`src/vendor/openadmet_eval/README.md`'s own convention ("if a change is needed, note it in a
wrapper module instead"), the fix lives here, not in the vendored file.

## `ensemble/`

- `caruana.py` — Caruana bagged ensemble selection (Caruana, Niculescu-Mizil, Crew & Ksikes 2004,
  ICML'04, Sections 2.1–2.3): sorted top-1 init, greedy selection-with-replacement with early
  stopping, bagged over random subsets of the candidate library. Used by notebook `11b` (selection)
  and reused by notebook `12` (retrain against the selected weights).

## `vendor/`

**Unmodified upstream code. Do not edit these files** — if a change is needed, add it to a wrapper
module above instead (`scoring.py`, `cv_bootstrap.py`), per `src/vendor/openadmet_eval/README.md`'s
own stated convention (that file already exists; this section does not replace it).

- `openadmet_eval/` — OpenADMET's official evaluation code, pulled byte-for-byte from
  `OpenADMET/CYP-Challenge-Tutorial` and verified against GitHub's blob SHAs at the time it was
  vendored. Used for every bootstrapped metric table (ST-RAE, MAE, R², Spearman, Kendall) so CV
  scores are computed the same way as the real leaderboard. **Caveat, not yet resolved**: this
  code's provenance against the live leaderboard backend has never been directly re-diffed against
  the tutorial repo's current HEAD — see `CLAUDE.md`'s Rules section.
- `validation/` — the tutorial's own `validate_activity_submission` (plus `tdi_validation.py`,
  vendored only because `src/vendor/validation/__init__.py` imports it), run unmodified before
  every submission.
