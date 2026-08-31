"""Generate the static manifest for the 5x5 (25-fold) CV model comparison --
`outputs/05_cv_comparison/manifest.csv`, one row per (config, repeat, fold) =
12 configs x 5 repeats x 5 folds = 300 rows.

STATIC PLAN, NEVER MUTATED BY THE RUNS: unlike a typical job-queue manifest, this file
is written once here and never written back to by `scripts/run_5x5_cv_comparison.py`.
"Done" is determined purely by checking whether that row's own output files already
exist on disk (see that script's `is_done`) -- not by any shared status column. This is
deliberate: it lets multiple driver processes (different `--family`/`--shard` values,
per the parallelization plan) run concurrently against the same manifest with no
shared mutable state and therefore no file-locking/race-condition concerns.

SEEDS: CLAUDE.md requires fixed seeds to vary per CV fold, not be shared globally, and
every seed used to be logged. One 32-bit seed is generated per (repeat, fold) pair --
NOT per (config, repeat, fold) -- via `np.random.SeedSequence(CV_SEED_BASE).spawn(25)`,
enumerated in a fixed (repeat 0..4, fold 0..4) order. That one seed is then reused by
every config evaluated on that (repeat, fold)'s held-out compounds for: the inner
train/val split (`assign_screen_split`'s own `seed` arg), every model's own random
state (RF/XGBoost/LightGbM `random_state`, Chemprop's `--data-seed`/`--pytorch-seed`),
and the bootstrap-scoring seed (`src.cv_bootstrap.per_fold_bootstrap_seed`). Sharing
one seed across configs *within* a (repeat, fold) is intentional, not a shortcut: the
independence CLAUDE.md's rule protects is *across* the 25 (repeat, fold) replicates for
a given model (the axis the later ANOVA/Tukey step compares) -- reusing the same seed
across configs *at* a given (repeat, fold) instead means every config is evaluated
against literally the same held-out compounds via the same bootstrap resampling
pattern, i.e. a paired comparison, which is a stricter, more standard design for this
kind of model comparison, not a weaker one.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

OUT = REPO_ROOT / "outputs" / "05_cv_comparison"
MANIFEST_PATH = OUT / "manifest.csv"

CV_SEED_BASE = 42  # matches this project's established global-SEED precedent (04a, final submission)
N_REPEATS = 5
N_FOLDS = 5

FEATURE_SETS = ["ecfp4_narrow", "chemeleon", "mordred_pca"]
ALGORITHMS = ["rf", "xgboost", "lightgbm"]
TABULAR_CONFIGS = [f"{feat}__{algo}" for feat in FEATURE_SETS for algo in ALGORITHMS]

CONFIGS = [
    ("naive_mean", "naive"),
    *[(c, "tabular") for c in TABULAR_CONFIGS],
    ("chemprop_randominit", "chemprop_randominit"),
    ("chemprop_chemeleoninit", "chemprop_chemeleoninit"),
]


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    seed_sequences = np.random.SeedSequence(CV_SEED_BASE).spawn(N_REPEATS * N_FOLDS)
    fold_seeds = {}
    i = 0
    for repeat in range(N_REPEATS):
        for fold in range(N_FOLDS):
            # generate_state returns a numpy uint32 array; take a plain Python int
            fold_seeds[(repeat, fold)] = int(seed_sequences[i].generate_state(1)[0])
            i += 1

    rows = []
    for config, family in CONFIGS:
        for repeat in range(N_REPEATS):
            for fold in range(N_FOLDS):
                rows.append(
                    {
                        "config": config,
                        "family": family,
                        "repeat": repeat,
                        "fold": fold,
                        "repeat_col": f"repeat_{repeat}",
                        "seed": fold_seeds[(repeat, fold)],
                    }
                )

    manifest = pd.DataFrame(rows)
    assert len(manifest) == len(CONFIGS) * N_REPEATS * N_FOLDS == 300, (
        f"expected 300 rows, got {len(manifest)}"
    )
    assert manifest.groupby(["repeat", "fold"])["seed"].nunique().eq(1).all(), (
        "seed did not come out constant within a (repeat, fold) pair"
    )
    assert manifest["seed"].nunique() == 25, "expected exactly 25 distinct fold seeds"

    manifest.to_csv(MANIFEST_PATH, index=False)
    print(f"wrote {MANIFEST_PATH}: {manifest.shape}")
    print(f"configs ({len(CONFIGS)}): {[c for c, _ in CONFIGS]}")
    print("fold seeds (repeat, fold) -> seed:")
    for (repeat, fold), seed in sorted(fold_seeds.items()):
        print(f"  (repeat={repeat}, fold={fold}) -> {seed}")


if __name__ == "__main__":
    main()
