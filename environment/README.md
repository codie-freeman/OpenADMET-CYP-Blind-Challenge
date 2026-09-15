# environment/

This directory holds archived environment specs. The live spec is `environment.yml` at the repo
root, not here — see the root `README.md`'s Quick Start for how to build it.

## `cyp-admet` vs. `cyp-admet-v2`

`cyp-admet-v2` (Python 3.11, native arm64) is the primary environment for all current work.
`cyp-admet` (x86_64, running under Rosetta emulation) still exists untouched on this machine as a
permanent fallback and is no longer used for new work.

**Why the migration happened, and what it bought:** `cyp-admet` ran every computation — RDKit,
PyTorch, Chemprop, tabular models — under x86_64/Rosetta emulation on Apple Silicon hardware.
`cyp-admet-v2` runs the same stack natively. The measured effect on Chemprop training specifically
is large: an isolated same-day timing check (`scripts/run_5x5_cv_comparison.py`'s own header
comment, 2026-08-30) found roughly a **13.2× speedup already present just from the arm64-native
migration itself** (one fold, two epochs, one trial — not yet validated across repeats or full
epoch counts), on top of which leaving `OMP_NUM_THREADS` unset (rather than this project's usual
`=1`) gave a further ~2.07× on top of that. A separate, smaller, more controlled micro-benchmark
(a single `RandomForestRegressor` fit, 500-row subset, 3 runs each) found a more modest ~1.4×
speedup (1.12s vs. 1.58s median) — a legitimate data point, but a single small-fit case, not a
full-pipeline projection. Both are real; they measure different things (a whole-pipeline
bottleneck vs. a cheap tabular fit) and shouldn't be conflated.

**The path gotcha:** `cyp-admet-v2` lives under `~/miniconda3-arm64` (a separate, native-arm64
conda installation), not the Intel `/opt/anaconda3` this machine's default `conda` points at — so
plain `conda activate cyp-admet-v2` fails unless that installation's own `conda` is on `PATH`
first. Every long-running script in this project that launches itself in the background invokes
the environment's Python by full path instead of relying on activation
(e.g. `/Users/codiefreeman/miniconda3-arm64/envs/cyp-admet-v2/bin/python`, see
`scripts/05b_run_cluster_cv.py` and similar) — this is deliberate, not an oversight, and the
pattern to follow for any new background-launched script.

**MPS remains unusable for Chemprop, and this is now precisely diagnosed, not a guess.** MPS's
`scatter_reduce` operation has no deterministic kernel, and Chemprop enables PyTorch Lightning's
`deterministic=True` mode whenever a PyTorch seed is set (`chemprop/cli/train.py` v2.3.1, lines
1798–1803, passed straight through to `pl.Trainer(..., deterministic=deterministic)` at line
1938) — and this project's convention is to always set a seed (`CLAUDE.md`'s Rules). The op is not
literally unimplemented on MPS; it's implemented but non-deterministic, and Chemprop's CLI exposes
no way to relax the determinism check. CPU-only Chemprop/CheMeleon training is therefore the
correct conclusion here, not a workaround pending a fix.

**Hardware constraint independent of the above:** RDKit and PyTorch each bundle their own OpenMP
runtime, and running them multi-threaded in the same process segfaults on this machine —
`OMP_NUM_THREADS=1` and `KMP_DUPLICATE_LIB_OK=TRUE` must be set before RDKit/PyTorch/Chemprop are
first imported, in every notebook and script (the one deliberate, documented exception being
`scripts/run_5x5_cv_comparison.py`, which leaves `OMP_NUM_THREADS` unset for the measured speedup
above).

## `archive/`

Environment specs are archived here, never overwritten in place, whenever the environment changes
— dated, descriptively named files:

- `environment_x86_64_rosetta_2026-08-30.yml` — the original `cyp-admet` spec, archived at the
  point of migrating to `cyp-admet-v2`.
- `environment_cyp-admet-v2_2026-09-14_pre-dimorphite-dl.yml` — `cyp-admet-v2`'s spec immediately
  before notebook `17` installed `dimorphite-dl`, which silently downgrades the pinned
  `rdkit==2026.3.3` to a `<2026` release as a hard dependency — fixed via a forced reinstall of the
  pinned RDKit version afterward (verified both packages still work together); see `environment.yml`'s
  own pip-section comments and `CLAUDE.md`'s notebook `17` entry for the full story.

Any future environment or major dependency change should be tested in a duplicated repo first
(copy the repo, build and validate the new environment there, rerun the relevant
determinism/timing checks) before being applied here directly — the pattern this migration itself
followed.
