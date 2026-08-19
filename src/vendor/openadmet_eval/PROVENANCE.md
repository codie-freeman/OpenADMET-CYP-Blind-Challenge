# Provenance

Source repo: https://github.com/OpenADMET/CYP-Challenge-Tutorial
Source directory: `evaluation/`
Commit hash: `858ae63ce79934113bccdb7fc65467de5f7b1935`
  (`main` branch HEAD at time of pull — "Merge pull request #1 from OpenADMET/evaluation-code")
Date pulled: 2026-08-19
License: Apache License 2.0 (as declared by the source repository)

## Files vendored (unmodified, byte-for-byte)

| File | Size (bytes) | Blob SHA (GitHub) |
|---|---|---|
| `config.py` | 2691 | `f48ccec5f754cca4e77a8738f937f9af4d1722a3` |
| `custom_scoring_functions.py` | 14174 | `cfc1136de48454f83be2738c3a70ab946161bac5` |
| `evaluate_predictions.py` | 20808 | `d9682805a7362a5b57487685d628f9c5136607a2` |
| `utils.py` | 1414 | `54942eba80418003efb19069ae357cf1c4f92321` |

Fetched from `raw.githubusercontent.com` pinned to the commit hash above, and verified
byte-identical (`cmp`) against the GitHub API's reported file sizes/blob SHAs before
being placed here.

To re-fetch this exact snapshot later:

```
https://raw.githubusercontent.com/OpenADMET/CYP-Challenge-Tutorial/858ae63ce79934113bccdb7fc65467de5f7b1935/evaluation/<filename>
```

## Pre-vendoring review notes

Before vendoring `evaluate_predictions.py` and `utils.py`, they were read in full to
check for CV-breaking assumptions. Summary (see project conversation history for full
detail):

- **No hardcoded test-set-size (750) assumption in the code itself.** `config.py` does
  define `ACTIVITY_DATASET_SIZE = 750`, but neither `evaluate_predictions.py` nor
  `utils.py` imports or references it — it's inert for our purposes. All functions
  operate on the actual shape of whatever DataFrame is passed in, so they run fine on
  any fold size.
- **`utils.bootstrap_sampling` uses a single fixed, globally-shared seed
  (`BOOTSTRAP_SEED = 0`), wrapped in `@lru_cache(maxsize=3)`.** This is intentional in
  the original context (ensures every leaderboard submission is bootstrapped on
  identical samples for fair comparison), but it conflicts with this project's rule
  that seeds must vary per CV fold and be logged. Every same-sized CV fold calling this
  function unmodified would get bit-identical bootstrap resampling. **Do not call this
  function directly inside a CV loop — wrap it (see repo `README.md` in this folder)
  if per-fold bootstrap variation is needed.**
- **Dependencies:** `evaluate_predictions.py` depends on `config.py` and `utils.py`
  (relative imports), and `config.py` depends on `custom_scoring_functions.py` — so all
  four vendored files are mutually required as a set. External packages:
  `numpy`, `pandas`, `scipy` (`kendalltau`, `spearmanr`), `scikit-learn` (all present in
  `cyp-admet`), and `loguru` (used only by `evaluate_predictions.py` — **not installed
  in `cyp-admet` as of this vendoring**; that module will fail to import until it is).
  `custom_scoring_functions.py` itself has no dependency beyond `numpy`/`pandas`.
- **Scope:** these files score a single prediction set against ground truth (with
  bootstrap resampling for a metric's mean/std), not a full multi-submission leaderboard
  workflow — no file I/O, no looping over submissions. Reusable directly for scoring one
  set of predictions (e.g. one CV fold), modulo the seed caveat above.
