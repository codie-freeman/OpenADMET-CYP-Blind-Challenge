"""Caruana bagged ensemble selection (Roadmap Step 3 rebuild, notebook 11).

Implements Caruana, Niculescu-Mizil, Crew & Ksikes (2004), "Ensemble Selection
from Libraries of Models," ICML'04, Sections 2.1-2.3:

- 2.1 Ensemble Selection: greedy forward selection **with replacement** over a
  fixed library of already-computed out-of-fold (OOF) predictions -- at each
  step, add whichever library member most improves the running ensemble's
  score on a fixed hillclimbing set, stopping as soon as no library member
  improves it further (rather than running a fixed number of steps and
  picking the best point in hindsight).
- 2.2 Selecting From Libraries of Models: **sorted ensemble initialization**
  -- start the greedy loop from the `n_init` individually best-scoring
  library members (rather than an empty ensemble), which the paper found
  reduces overfitting to the hillclimbing set from a purely greedy start.
- 2.3 Bagged Ensemble Selection: repeat the whole sorted-init + greedy
  procedure `n_bags` times, each run restricted to a random `bag_frac`
  subset **of the candidate library** (not the data -- this project's task
  brief is explicit that bagging is over library subsets), then combine the
  `n_bags` runs' selected multisets by selection frequency into final
  per-config weights. This is the paper's mitigation for a known failure
  mode of plain (2.1+2.2) selection: picking the best-scoring combination out
  of several similar candidates using the same data used to score them is an
  overfitting risk in its own right (see this project's CLAUDE.md and
  notebook 08/10/10b history) -- bagging over the library reduces variance in
  which models end up selected.

Provenance note: the task brief that requested this module named an external
reference implementation, `jeremycheminf/openadmet_scripts`'s
`src/cyp_submission/caruana.py`, as "already reviewed in this project." A
repo-wide search (filesystem, `git log --all`, all branches, this session's
other working directories) turned up no local trace of that file -- it was
not actually available in this project to consult, despite the brief's
description. It does exist externally, at a fuller path than the brief gave
(`CYP_Challenge/src/cyp_submission/caruana.py`, not the bare
`src/cyp_submission/caruana.py` named in the brief) -- the same
`jeremycheminf/openadmet_scripts` repo notebook 09 already cites (and
verified via the GitHub API) for a sibling file, `calibration.py`. Fetched
directly from GitHub this session (raw.githubusercontent.com) and reviewed as
the working-implementation-pattern reference the brief intended. Two design
choices below follow it deliberately: **early stopping** the per-bag greedy
loop as soon as no library member strictly improves the running score
(rather than running a fixed `max_iter` and picking the best-scoring prefix
in hindsight), and this project's own **default hyperparameter values**
(`n_bags=20, bag_frac=0.5, n_init=1, max_iter=50, seed=42`, set explicitly by
the caller in notebook 11b -- never defaulted in this module itself) match
that file's own defaults, now that they're a reviewed community reference
for this specific challenge rather than an arbitrary guess. Scoring, input
validation, the per-config dict return, and per-bag logging are this
project's own additions -- not present in the reference, and not copied from
it; this module is an independent implementation, not a port.

Scoring uses this project's actual ST-RAE metric --
``src.vendor.openadmet_eval.custom_scoring_functions.rae_soft_threshold_absolute_error``,
the same function that backs "ST-RAE" in
``src.vendor.openadmet_eval.config.ACTIVITY_METRICS`` -- not a generic error
function. It is called directly (no bootstrap resampling) as the greedy
search's hillclimbing objective: bootstrapping is this project's mechanism
for *significance testing* (``src.cv_bootstrap.per_fold_bootstrap_seed``),
not something a hillclimbing objective needs, and calling the plain metric
function directly keeps the many candidate evaluations a selection run
performs (``n_bags * max_iter * library_subset_size``) cheap. The full
bootstrap pipeline is used downstream, once, in notebook 11b to report the
final selected ensemble's fold-level ST-RAE for the paired comparison against
notebook 08 -- see that notebook, not this module.

Every hyperparameter (n_bags, bag_frac, n_init, max_iter, seed) is a
required argument with no default, per this project's "log every seed used,
no silent defaults" rule -- callers must state them explicitly, and this
module logs them (and every bag's outcome) via its module logger before and
during a selection run.
"""

from __future__ import annotations

import logging
from typing import Sequence

import numpy as np

from src.vendor.openadmet_eval.custom_scoring_functions import (
    rae_soft_threshold_absolute_error,
)

logger = logging.getLogger(__name__)


def _st_rae(
    y_pred: np.ndarray, y_true: np.ndarray, conf_low: np.ndarray, conf_high: np.ndarray
) -> float:
    """This project's ST-RAE metric, called directly (no bootstrap)."""
    return rae_soft_threshold_absolute_error(
        y_true, y_pred, y_true_upper=conf_high, y_true_lower=conf_low
    )


def _select_one_bag(
    bag_index: int,
    bag_seed: int,
    oof: np.ndarray,
    y_true: np.ndarray,
    conf_low: np.ndarray,
    conf_high: np.ndarray,
    config_names: Sequence[str],
    bag_frac: float,
    n_init: int,
    max_iter: int,
) -> list[str]:
    """One bag of Caruana Sections 2.1+2.2: sorted-top-`n_init` initialization
    over a random `bag_frac` subset of the library, then greedy
    selection-with-replacement, stopping as soon as no subset member improves
    the running average -- or after `max_iter` additions, whichever comes
    first. Only the *library* is subsetted per bag; `oof`/`y_true`/`conf_low`/
    `conf_high` are always the full hillclimbing set.
    """
    rng = np.random.default_rng(bag_seed)
    n_configs = len(config_names)
    bag_size = max(1, round(bag_frac * n_configs))
    bag_n_init = min(n_init, bag_size)  # cap init to this bag's own (possibly small) subset
    subset_idx = np.sort(rng.choice(n_configs, size=bag_size, replace=False))

    # Sorted top-n_init initialization (Section 2.2): score every subset member
    # alone on the hillclimbing set, keep the individually-best `bag_n_init`.
    individual_scores = np.array(
        [_st_rae(oof[:, idx], y_true, conf_low, conf_high) for idx in subset_idx]
    )
    init_order = np.argsort(individual_scores)  # ascending ST-RAE = best first
    selected = list(subset_idx[init_order[:bag_n_init]])

    running_sum = oof[:, selected].sum(axis=1)
    current_score = _st_rae(running_sum / len(selected), y_true, conf_low, conf_high)

    # Greedy selection with replacement (Section 2.1): at each step, add whichever
    # library member (repeats allowed) most improves the running average; stop the
    # instant nothing does, rather than forcing exactly max_iter additions.
    n_added = 0
    for _ in range(max_iter):
        candidate_sums = running_sum[:, None] + oof[:, subset_idx]  # (n_samples, bag_size)
        candidate_preds = candidate_sums / (len(selected) + 1)
        candidate_scores = np.array(
            [_st_rae(candidate_preds[:, j], y_true, conf_low, conf_high) for j in range(bag_size)]
        )
        best_j = int(np.argmin(candidate_scores))
        best_score = float(candidate_scores[best_j])
        if best_score >= current_score:
            break  # no library member improves on the current ensemble -- stop growing it
        best_add = int(subset_idx[best_j])
        selected.append(best_add)
        running_sum = running_sum + oof[:, best_add]
        current_score = best_score
        n_added += 1

    final_selected = [config_names[i] for i in selected]
    logger.info(
        "  bag %d (seed=%d): library subset=%s -> selected=%s "
        "(final ST-RAE=%.4f, %d greedy addition(s) beyond the %d-model init, cap was %d)",
        bag_index,
        bag_seed,
        [config_names[i] for i in subset_idx],
        final_selected,
        current_score,
        n_added,
        bag_n_init,
        max_iter,
    )
    return final_selected


def caruana_bagged_selection(
    oof: np.ndarray,
    y_true: np.ndarray,
    conf_low: np.ndarray,
    conf_high: np.ndarray,
    config_names: Sequence[str],
    n_bags: int,
    bag_frac: float,
    n_init: int,
    max_iter: int,
    seed: int,
) -> dict[str, float]:
    """Run Caruana bagged ensemble selection and return per-config weights.

    Args:
        oof (np.ndarray): OOF prediction matrix, shape (n_samples, n_configs) --
            one column per candidate model config, one row per (out-of-fold)
            prediction.
        y_true (np.ndarray): True values, shape (n_samples,), aligned to `oof`.
        conf_low (np.ndarray): Lower confidence bound on `y_true`, shape
            (n_samples,) -- used by this project's ST-RAE metric
            (soft-thresholded RAE), not a per-model prediction interval.
        conf_high (np.ndarray): Upper confidence bound on `y_true`, counterpart
            to `conf_low`.
        config_names (Sequence[str]): Name of each `oof` column, length
            n_configs.
        n_bags (int): Number of bagging iterations (Section 2.3). Required, no
            default.
        bag_frac (float): Fraction of the library sampled (without replacement)
            per bag. Required, no default.
        n_init (int): Number of individually-best library members used to
            initialize each bag's ensemble (Section 2.2), capped down to that
            bag's own subset size if `bag_frac` makes it smaller than n_init.
            Required, no default.
        max_iter (int): Maximum number of greedy selection-with-replacement
            rounds run per bag, beyond the `n_init` initialization (Section
            2.1) -- an upper bound; a bag's greedy loop stops earlier as soon
            as no library member improves the running score. Required, no
            default.
        seed (int): Master seed. Per-bag seeds are drawn from this master seed
            and logged individually -- never a single global seed reused across
            bags, per this project's seeding rule.

    Returns:
        dict[str, float]: One entry per `config_names` entry (including configs
            that received zero weight), summing to 1.0. A config's weight is its
            selection frequency across all bags: the number of times it appears
            in a bag's selected multiset, divided by the total number of
            selections across all bags.

    Raises:
        ValueError: On shape mismatches or NaNs anywhere in the inputs -- this
            function never silently drops or imputes incomplete rows.
        RuntimeError: If every bag selects an empty ensemble (should not happen
            for n_init >= 1, but checked explicitly rather than returning an
            all-zero/undefined weight dict).

    """
    oof = np.asarray(oof, dtype=float)
    y_true = np.asarray(y_true, dtype=float)
    conf_low = np.asarray(conf_low, dtype=float)
    conf_high = np.asarray(conf_high, dtype=float)

    n_samples, n_configs = oof.shape
    if len(config_names) != n_configs:
        raise ValueError(
            f"config_names has {len(config_names)} entries but oof has {n_configs} columns"
        )
    for name, arr in [("y_true", y_true), ("conf_low", conf_low), ("conf_high", conf_high)]:
        if arr.shape != (n_samples,):
            raise ValueError(f"{name} has shape {arr.shape}, expected ({n_samples},)")
    if (
        np.isnan(oof).any()
        or np.isnan(y_true).any()
        or np.isnan(conf_low).any()
        or np.isnan(conf_high).any()
    ):
        raise ValueError(
            "NaNs found in caruana_bagged_selection inputs -- refusing to select "
            "on incomplete data"
        )
    if not (0 < bag_frac <= 1):
        raise ValueError(f"bag_frac must be in (0, 1], got {bag_frac}")
    if n_init < 1:
        raise ValueError(f"n_init must be >= 1, got {n_init}")

    hyperparameters = {
        "n_bags": n_bags,
        "bag_frac": bag_frac,
        "n_init": n_init,
        "max_iter": max_iter,
        "seed": seed,
    }
    logger.info(
        "Caruana bagged selection starting: n_samples=%d, n_configs=%d, configs=%s, "
        "hyperparameters=%s",
        n_samples,
        n_configs,
        list(config_names),
        hyperparameters,
    )

    rng = np.random.default_rng(seed)
    bag_seeds = rng.integers(0, 2**31 - 1, size=n_bags)
    logger.info("bag seeds (drawn from master seed=%d): %s", seed, bag_seeds.tolist())

    selection_counts = {name: 0 for name in config_names}
    for b in range(n_bags):
        bag_selected = _select_one_bag(
            b,
            int(bag_seeds[b]),
            oof,
            y_true,
            conf_low,
            conf_high,
            config_names,
            bag_frac,
            n_init,
            max_iter,
        )
        for name in bag_selected:
            selection_counts[name] += 1

    total_selections = sum(selection_counts.values())
    if total_selections == 0:
        raise RuntimeError(
            "Caruana selection produced an empty ensemble across all bags -- no "
            "weights to return"
        )
    weights = {name: count / total_selections for name, count in selection_counts.items()}
    logger.info(
        "final weights (selection frequency across %d bags, %d total selections): %s",
        n_bags,
        total_selections,
        {k: round(v, 4) for k, v in weights.items() if v > 0},
    )
    return weights
