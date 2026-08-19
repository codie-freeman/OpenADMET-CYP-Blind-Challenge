# CLAUDE.md

Environment: conda env `cyp-admet`, Python 3.11, Apple Silicon (MPS). VS Code + Claude Code.

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
