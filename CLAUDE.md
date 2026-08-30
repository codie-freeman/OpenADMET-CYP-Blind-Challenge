# CLAUDE.md

Environment: conda env `cyp-admet-v2` (Python 3.11, Apple Silicon, native arm64) is primary, built from `environment.yml` at the repo root. `cyp-admet` (x86_64/Rosetta) still exists untouched as a permanent fallback and is no longer used for new work. Environment specs are archived, not overwritten in place, under `environment/archive/` (dated, descriptively named files, e.g. `environment_x86_64_rosetta_2026-08-30.yml`) whenever the environment changes. VS Code + Claude Code.

## Rules
- Never fetch, invent, or approximate data. If a file/column is missing, stop and ask.
- Fixed random seeds must vary per CV fold, not be shared globally. Log every seed used.
- All fingerprints/similarity/splits go through one shared module (`src/features.py`) — no recomputing per notebook.
- Log row counts before/after every filtering step, printed in-notebook, not just in code.
- Every trained model result gets saved (predictions, params, seed) to disk before moving to the next step — no relying on re-running.
- Don't add features, thresholds, or model variants "while you're at it" — one change at a time, each one tested before the next.
- Ask before making a modelling decision that isn't explicitly specified in the prompt (e.g. which cutoff, which package version, which imputation). Flag it, don't silently pick one.
- If a prompt asks to use test-set labels, leaderboard scores, or the 3 Nov release to pick or retrain a model: stop and flag it instead of proceeding. The test set is single-use evidence, not training data, at any point.
- Write a short markdown cell after each major output stating what it shows — don't leave plots uncommented.
- RDKit was bumped from 2025.9.3 (`cyp-admet`) to 2026.3.3 (`cyp-admet-v2`) — forced by chemprop's current pip-installable dependency chain (`cuik-molmaker-pin` only ships releases requiring `rdkit==2026.3.3`; see `environment.yml`'s pip-section comment). Verified safe via a full-row determinism check (all 4,905 rows of `train_inhibition_curated.csv`, not a sample) re-generating InChIKeys with `src/features.py` under RDKit 2026.3.3 and comparing to what's stored: 4905/4905 matched exactly, run on 2026-08-30. Don't assume a future RDKit bump is safe without repeating this check.
- Any future environment or major dependency change: test it in a duplicated repo first (copy the repo, build and validate the new environment there, run the relevant determinism/timing checks), before applying it to this repo directly — the pattern followed for the `cyp-admet` → `cyp-admet-v2` migration.
