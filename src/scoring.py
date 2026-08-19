"""Scoring helpers built on top of the vendored OpenADMET evaluation code.

See `src/vendor/openadmet_eval/README.md` and `PROVENANCE.md` — the vendored files
themselves are never edited; anything project-specific goes here instead.
"""

import numpy as np

from src.vendor.openadmet_eval.config import REGRESSION_ENDPOINTS
from src.vendor.openadmet_eval.custom_scoring_functions import (
    rae_soft_threshold_absolute_error,
)


def get_ma_st_rae(
    y_true_dict: dict[str, np.ndarray],
    y_pred_dict: dict[str, np.ndarray],
    conf_bounds_dict: dict[str, dict[str, np.ndarray]],
) -> float:
    """Macro-averaged soft-threshold RAE (MA ST-RAE) across the four regression endpoints.

    Computes `rae_soft_threshold_absolute_error` independently for each endpoint in
    `config.REGRESSION_ENDPOINTS` (CYP1A2/CYP2C9/CYP2D6/CYP3A4 direct-inhibition pIC50),
    then macro-averages the four per-endpoint scores with a plain arithmetic mean —
    matching how `evaluate_predictions.compute_macro_bootstrap_results` computes the
    `MACRO_ENDPOINT_LABEL` ("MA") row (plain mean, not a Fisher z-transform; see that
    function's docstring for why).

    Each endpoint is expected to carry its own arrays, since compounds are not tested
    against every CYP isoform — the four endpoints generally have different support
    (different compounds, different lengths) and must not be concatenated together.

    Args:
        y_true_dict: `{endpoint: y_true}` for each endpoint in `REGRESSION_ENDPOINTS`.
        y_pred_dict: `{endpoint: y_pred}`, aligned elementwise with `y_true_dict[endpoint]`.
        conf_bounds_dict: `{endpoint: {"lower": y_true_lower, "upper": y_true_upper}}`,
            aligned elementwise with `y_true_dict[endpoint]` — the credible-interval
            bounds consumed by the soft-thresholded RAE metric (e.g. this project's
            `_conf_low` / `_conf_high` columns).

    Returns:
        The macro-averaged ST-RAE, as a single float.

    Raises:
        KeyError: If any endpoint in `REGRESSION_ENDPOINTS` is missing from
            `y_true_dict`, `y_pred_dict`, or `conf_bounds_dict` (or from a bounds
            dict's `"lower"`/`"upper"` keys).

    """
    missing = {
        endpoint
        for endpoint in REGRESSION_ENDPOINTS
        if endpoint not in y_true_dict or endpoint not in y_pred_dict or endpoint not in conf_bounds_dict
    }
    if missing:
        raise KeyError(
            f"Missing endpoint(s) {sorted(missing)} — get_ma_st_rae requires all of "
            f"{REGRESSION_ENDPOINTS} in y_true_dict, y_pred_dict, and conf_bounds_dict."
        )

    per_endpoint_scores = []
    for endpoint in REGRESSION_ENDPOINTS:
        bounds = conf_bounds_dict[endpoint]
        score = rae_soft_threshold_absolute_error(
            np.asarray(y_true_dict[endpoint]),
            np.asarray(y_pred_dict[endpoint]),
            y_true_upper=np.asarray(bounds["upper"]),
            y_true_lower=np.asarray(bounds["lower"]),
        )
        per_endpoint_scores.append(score)

    return float(np.mean(per_endpoint_scores))
