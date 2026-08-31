"""Per-fold bootstrap-seed wrapper for the 5x5 CV comparison
(`scripts/run_5x5_cv_comparison.py`) -- never edits the vendored file itself
(`src/vendor/openadmet_eval/README.md`'s own convention: "if a change is needed, note
it in a wrapper module instead"). Directly implements the fix flagged in that file's
own `PROVENANCE.md`: "Do not call [bootstrap_sampling] directly inside a CV loop --
wrap it ... if per-fold bootstrap variation is needed."

`openadmet_eval.utils.bootstrap_sampling(size, n_repeats)` is `@lru_cache`d and reads
a module-level `BOOTSTRAP_SEED = 0` global -- fixed, not a function argument. Two
things must happen together for a fold to actually get a fresh sample: the module
attribute must be set to that fold's seed, AND the cache must be cleared -- the cache
key is `(size, n_repeats)` only, so a result already cached for a given held-out-fold
size is returned unchanged even after the seed changes, and many folds/configs in this
comparison share the same size. Confirmed directly (not assumed): clearing the cache
alone, with the seed left at its default, still returns byte-identical draws call to
call -- the seed change is what matters, the cache clear is what makes it take effect.
"""

import contextlib

from src.vendor.openadmet_eval import utils as _openadmet_utils


@contextlib.contextmanager
def per_fold_bootstrap_seed(seed: int):
    """Force every `bootstrap_sampling` call made inside this block to draw fresh
    samples seeded by `seed`, then restore the vendored module's default seed (and
    clear the cache again) on exit so a later default-seed caller elsewhere in the
    process isn't left pointing at a stale non-default result.
    """
    original_seed = _openadmet_utils.BOOTSTRAP_SEED
    _openadmet_utils.BOOTSTRAP_SEED = seed
    _openadmet_utils.bootstrap_sampling.cache_clear()
    try:
        yield
    finally:
        _openadmet_utils.BOOTSTRAP_SEED = original_seed
        _openadmet_utils.bootstrap_sampling.cache_clear()
