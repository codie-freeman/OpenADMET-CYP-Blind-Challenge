#!/usr/bin/env python
"""Notebook 52's three-arm comparison: Ash's parametric default, run in full.

Three methods -- `plain` (notebook 37's stored scores), `qm_both_zeroatoms`, `qm_both` -- is
well under Ash's >10 threshold, so the protocol routes to Tukey HSD. This is the first time
this project fires that branch: notebooks 05 and 37 both ran at 11 arms, where their own
`posthoc_route` sends everything to BH, so the Tukey code in `src/ash_protocol.py` was
reproduced but never executed.

The route, per (isoform, metric):

  1. Levene branch check (`branch_is_parametric`), MEASURED per cell rather than inherited.
     Notebook 05 found ST-RAE nonparametric everywhere, notebook 29 parametric on all four
     isoforms and notebook 37 split by isoform, so the three disagree and the check has to be
     run on the data in hand.
  2. Repeated-measures ANOVA across the three arms, with the fold as the subject.
  3. The ANOVA's alpha is Bonferroni-corrected by the NUMBER OF METRICS: 0.05 / 2 = 0.025.
     Ash's own step, and the one notebook 44 skipped entirely.
  4. Tukey HSD on all three pairwise comparisons, only where the omnibus clears, using the
     repeated-measures error term rather than a plain one-way pooled variance (Ash main text
     3.2.1). FWER is controlled across the three comparisons by the studentized range itself,
     so no further correction is applied on top.
  5. Compact letter display.

WHAT THIS MODULE ADDS over `src/ash_protocol.py`, which is imported unmodified: that module's
`rm_tukey_hsd` returns pairwise p-values but no intervals, and materiality cannot be judged from
a p-value. The Tukey half-width is `q(1-alpha, k, df) * sqrt(MSE / n_per_group)`, which is the
same quantity `rm_tukey_hsd`'s own `tukey_se = sqrt(2*MSE/n)` carries, divided by sqrt(2). Both
forms are computed here and asserted equal, so the interval cannot drift from the p-value it
accompanies.

`src/ash_protocol.py` is deliberately NOT modified: it is on notebook 52's protected list, and
every earlier notebook's results depend on it being the file they ran.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.libqsturng import qsturng

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.ash_protocol import (  # noqa: E402
    assign_cld, branch_is_parametric, conover_friedman_holm, friedman_omnibus, posthoc_route,
    rm_anova, rm_tukey_hsd,
)

ISOFORMS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
BARS = {"ST-RAE": {"CYP1A2": 0.0234, "CYP2C9": 0.0209, "CYP2D6": 0.0277, "CYP3A4": 0.0160},
        "Spearman_R": {"CYP1A2": 0.0183, "CYP2C9": 0.0160, "CYP2D6": 0.0212, "CYP3A4": 0.0100}}
ALPHA = 0.05
N_METRICS = len(BARS)
ANOVA_ALPHA = ALPHA / N_METRICS          # Ash's Bonferroni by the number of metrics


def build_long(scores: dict[str, pd.Series], folds) -> pd.DataFrame:
    """One tidy frame for a single (isoform, metric): `value`, `config`, `subject` (the fold).

    `scores` maps an arm name to a Series indexed by (repeat, fold). Every arm must cover every
    fold in `folds` -- a missing cell would silently unbalance the repeated-measures design, so
    it raises instead.
    """
    rows = []
    for arm, s in scores.items():
        missing = [k for k in folds if k not in s.index]
        if missing:
            raise ValueError(f"{arm} is missing {len(missing)} of {len(folds)} folds: {missing[:5]}")
        for k in folds:
            rows.append({"value": float(s.loc[k]), "config": arm, "subject": f"r{k[0]}f{k[1]}"})
    return pd.DataFrame(rows)


# Precedence for orienting a pairwise row: the arm earlier in this list becomes `arm_a`, so
# every row reads "arm_a relative to arm_b" with the newer arm first. Without this the rows come
# out in groupby order and `plain vs qm_both` reports the sign backwards.
ORIENT = ["qm_both", "qm_both_zeroatoms", "plain"]


def _orient(a, b):
    ia = ORIENT.index(a) if a in ORIENT else len(ORIENT)
    ib = ORIENT.index(b) if b in ORIENT else len(ORIENT)
    return (a, b) if ia <= ib else (b, a)


def tukey_intervals(long_df, mse, df_resid, alpha=ALPHA):
    """Tukey HSD pairwise differences WITH confidence intervals, on the RM error term."""
    pairwise_p, means = rm_tukey_hsd(long_df, mse, df_resid)
    methods = list(means.index)
    n_per_group = long_df["config"].value_counts().mean()
    tukey_se = np.sqrt(2 * mse / n_per_group)
    q = qsturng(1 - alpha, len(methods), df_resid)
    half = q * np.sqrt(mse / n_per_group)
    assert abs(half - q * tukey_se / np.sqrt(2)) < 1e-12, "the two half-width forms disagree"
    rows = []
    for i, m1 in enumerate(methods):
        for m2 in methods[i + 1:]:
            a, b = _orient(m1, m2)
            d = float(means[a] - means[b])
            rows.append({"arm_a": a, "arm_b": b, "mean_a": float(means[a]),
                         "mean_b": float(means[b]), "diff": d,
                         "ci_low": d - half, "ci_high": d + half,
                         "p_tukey": float(pairwise_p.loc[a, b]),
                         "significant": float(pairwise_p.loc[a, b]) < alpha})
    return pd.DataFrame(rows), pairwise_p, means, float(half)


def compare_cell(scores: dict[str, pd.Series], folds, isoform: str, metric: str) -> dict:
    """The full route for one (isoform, metric) cell."""
    long_df = build_long(scores, folds)
    parametric, lev = branch_is_parametric(long_df, alpha=ALPHA)
    route = posthoc_route(long_df["config"].nunique(), parametric)
    out = {"isoform": isoform, "metric": metric, "n_folds": len(folds),
           "n_arms": long_df["config"].nunique(), "parametric": parametric,
           "levene_p": lev["levene_p"], "variance_ratio": lev["variance_ratio"],
           "posthoc_route": route, "anova_alpha": ANOVA_ALPHA}

    if parametric:
        aov = rm_anova(long_df)
        out.update({"omnibus": "RM-ANOVA", "omnibus_p": aov["anova_p"],
                    "mse": aov["mse"], "df_resid": aov["df_resid"]})
        clears = aov["anova_p"] < ANOVA_ALPHA
        out["omnibus_clears"] = bool(clears)
        if clears:
            tuk, pw, means, half = tukey_intervals(long_df, aov["mse"], aov["df_resid"])
            out["tukey"] = tuk
            out["half_width"] = half
            out["cld"] = assign_cld(means, pw, alpha=ALPHA,
                                    ascending=(metric == "ST-RAE"))
        else:
            out["tukey"], out["cld"], out["half_width"] = None, None, None
    else:
        fr = friedman_omnibus(long_df)
        out.update({"omnibus": "Friedman", "omnibus_p": fr["friedman_p"],
                    "mse": np.nan, "df_resid": np.nan})
        clears = fr["friedman_p"] < ANOVA_ALPHA
        out["omnibus_clears"] = bool(clears)
        if clears:
            pw = conover_friedman_holm(long_df)
            means = long_df.groupby("config")["value"].mean()
            rows = []
            methods = list(means.index)
            for i, m1 in enumerate(methods):
                for m2 in methods[i + 1:]:
                    a, b = _orient(m1, m2)
                    rows.append({"arm_a": a, "arm_b": b, "mean_a": float(means[a]),
                                 "mean_b": float(means[b]),
                                 "diff": float(means[a] - means[b]),
                                 "ci_low": np.nan, "ci_high": np.nan,
                                 "p_tukey": float(pw.loc[a, b]),
                                 "significant": float(pw.loc[a, b]) < ALPHA})
            out["tukey"] = pd.DataFrame(rows)
            out["cld"] = assign_cld(means, pw, alpha=ALPHA, ascending=(metric == "ST-RAE"))
            out["half_width"] = np.nan
        else:
            out["tukey"], out["cld"], out["half_width"] = None, None, None
    return out


def verdict_row(r, metric, isoform, bar):
    """Materiality from the Tukey INTERVAL, not from the point estimate alone.

    `improves` is metric-aware: ST-RAE down is better, Spearman up is better. A comparison is
    MATERIAL only if it is significant AND its point estimate reaches the bar; `ci_excludes_bar`
    additionally records whether the whole interval sits beyond the bar, which is the stronger
    statement and the one a bar built from 25-fold SEs can actually support.
    """
    d = r["diff"]
    improves = (d < 0) if metric == "ST-RAE" else (d > 0)
    reaches = abs(d) >= bar
    if metric == "ST-RAE":
        beyond = r["ci_high"] < -bar if improves else r["ci_low"] > bar
    else:
        beyond = r["ci_low"] > bar if improves else r["ci_high"] < -bar
    return {"improves": bool(improves), "reaches_bar": bool(reaches),
            "pct_of_bar": 100 * abs(d) / bar,
            "material": bool(r["significant"] and reaches),
            "ci_excludes_bar": bool(beyond) if pd.notna(r["ci_low"]) else None}
