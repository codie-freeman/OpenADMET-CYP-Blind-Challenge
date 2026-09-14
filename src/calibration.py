"""Calibration of OOF predictions onto an external, independent population's
mean/std -- distinct from the OWN-population placement corrections already
tried in `outputs/09_placement_recalibration` and
`outputs/11b_caruana_selection`.

Those two checks both fit `slope, intercept = np.polyfit(y_pred, y_true)`
against the *same* population that produced `y_pred` (this project's own
training/OOF population) -- an honest ranking-quality diagnostic, but one
that can never correct a population-level mean/std mismatch between this
project's training distribution and the true blind-test distribution,
because it has no information about that distribution. Both found the
scale+bias share of the R^2 gap under the 10% pursuit threshold and adopted
identity (slope=1, intercept=0) for every isoform.

The functions here fit a *different* correction: they still use this
project's own pooled OOF predictions to estimate `rho` (how well the model's
raw predictions track truth, on the only ground truth available), but they
target the mean/std of an external reference population instead of this
project's own training population -- e.g. PubChem AID 1851 (notebook 13),
used purely as a population-moment estimate, never as training/fine-tuning
data. This is the same mechanism two real entrants on this challenge
(SuperCowPowers, jeremycheminf) report using to correct a placement mismatch
between this project's training population and the true blind population.

`fit_blind_population_calibration` derives a shrunk variance-matching
transform (slope scaled by `rho`, not a full 1:1 variance match) so that
applying it moves a model's raw predictions toward the target population's
mean/std in proportion to how much the model's own predictions are actually
correlated with truth. `apply_calibration` applies the resulting affine
transform to new predictions.
"""

import numpy as np


def fit_blind_population_calibration(oof_pred, oof_true, target_mean, target_std) -> dict:
    """Fit slope/intercept mapping `oof_pred` onto an external population's moments.

    `rho` is the Pearson correlation between `oof_pred` and `oof_true` --
    computed on this project's own OOF data, since that is the only ground
    truth available. `target_mean`/`target_std` come from an external
    reference population (e.g. AID 1851), never from `oof_true` itself.
    """
    oof_pred = np.asarray(oof_pred, dtype=float)
    oof_true = np.asarray(oof_true, dtype=float)

    rho = float(np.corrcoef(oof_pred, oof_true)[0, 1])
    own_mean = float(np.mean(oof_pred))
    own_std = float(np.std(oof_pred, ddof=0))

    slope = rho * (target_std / own_std)
    intercept = target_mean - slope * own_mean

    return dict(
        slope=slope,
        intercept=intercept,
        rho=rho,
        own_mean=own_mean,
        own_std=own_std,
        target_mean=float(target_mean),
        target_std=float(target_std),
        target_source="aid1851",
    )


def apply_calibration(raw_pred, params: dict) -> np.ndarray:
    """Apply an affine transform fit by `fit_blind_population_calibration`."""
    return params["slope"] * np.asarray(raw_pred, dtype=float) + params["intercept"]
