"""Ash et al. (2025) method-comparison protocol: the shared implementation.

Ash, Wognum, Rodriguez-Perez et al., "Practically Significant Method Comparison Protocols for
Machine Learning in Small Molecule Drug Discovery", J. Chem. Inf. Model. 2025, 65, 9398-9411,
DOI 10.1021/acs.jcim.5c01609.

WHY THIS MODULE EXISTS. These functions were written in notebook 05, copied verbatim into
notebook 37, and would have been copied a third time by notebook 46. They are lifted here so
there is one implementation, per CLAUDE.md's own "goes through one shared module" rule. The
bodies are unchanged from notebook 37 cell 116, which is itself notebook 05 cell 7 verbatim;
`branch_is_parametric` and `posthoc_route` are new, and only name rules those notebooks already
applied inline.

THE ROUTE, from the paper's own Section 3.2.2 (main text, verified against the PDF at sha256
7c2caa48c19daff58a768f7893cc4521e1f29fc7f3d19d44d57a0cb3ea6e4d2a):

    "(1) Perform repeated measures ANOVA for each metric. (2) Perform Bonferroni correction on
    the ANOVA significance level. Divide the significance level by the number of metrics to
    obtain an adjusted significance level. (3) For all metrics with a significant ANOVA test,
    perform a Tukey HSD posthoc test."

and, on restricting the family to one-against-all rather than all pairs:

    "In the cases where only the comparisons of one method against all others are of interest,
    then less tests need to be corrected for than all pairwise comparisons as in Tukey HSD.
    However, since Tukey HSD has good statistical power compared to other procedures, we
    recommend this for general use. This gives researchers a 'license to fish', providing
    protection against any pairwise comparison that might be performed."

A PROVENANCE LIMIT, recorded because it bears on which branch to trust. Only the 14-page main
text is available in this project; the Supporting Information is not in the repository or in the
reference library. The main text carries the parametric route above, and says of Tukey that it
"has good statistical power compared to other procedures" and controls the family-wise error
rate "regardless of the number of tests performed". It does NOT contain the nonparametric
workflow (`Friedman` and `Conover` both appear zero times), the >10-method BH exception, or the
9x variance tolerance -- notebooks 05 and 37 cite those to SI Sections A.2 and C.1, which cannot
be checked here. `rm_tukey_hsd` is therefore the branch whose justification is verifiable, and a
comparison that can stay at or under `LARGE_COMPARISON_THRESHOLD` methods should.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pingouin as pg
import scikit_posthocs as sp
from scipy.stats import levene, ttest_rel
from statsmodels.stats.libqsturng import psturng
from statsmodels.stats.multitest import multipletests

# Ash SI Section A.2 (not verifiable here -- see the module docstring).
LARGE_COMPARISON_THRESHOLD = 10
# Ash SI Section C.1's tolerance for equally-sized groups (likewise not verifiable here).
VARIANCE_RATIO_TOLERANCE = 9


def levene_check(long_df):
    groups = [g["value"].to_numpy() for _, g in long_df.groupby("config")]
    stat, pvalue = levene(*groups)
    variances = {name: float(np.var(g["value"].to_numpy(), ddof=1))
                 for name, g in long_df.groupby("config")}
    min_var, max_var = min(variances.values()), max(variances.values())
    ratio = float("inf") if min_var == 0 else max_var / min_var
    return {"levene_stat": float(stat), "levene_p": float(pvalue), "variance_ratio": ratio,
            "zero_variance_group": min_var == 0}


def rm_anova(long_df):
    """Repeated-measures ANOVA via pingouin. Row 0 = the `config` main effect, row 1 = its
    within-subject Error term -- confirmed directly against this project's installed pingouin
    0.6.1 (`detailed=True`) before relying on this indexing, not assumed.
    """
    aov = pg.rm_anova(dv="value", within="config", subject="subject", data=long_df, detailed=True)
    return {"anova_p": float(aov.loc[0, "p_unc"]), "mse": float(aov.loc[1, "MS"]),
            "df_resid": float(aov.loc[1, "DF"])}


def rm_tukey_hsd(long_df, mse, df_resid):
    """Pairwise Tukey HSD using the repeated-measures ANOVA's own error term rather than a plain
    one-way pooled variance -- Ash main text Section 3.2.1. Translated from the paper's own
    public reference implementation (polaris-hub/polaris-method-comparison,
    ADME_example/model_comparison.py::rm_tukey_hsd). Not exercised here: with 11 arms notebook
    05's own >10-group rule routes post-hoc testing to BH.
    """
    means = long_df.groupby("config")["value"].mean()
    methods = list(means.index)
    n_groups = len(methods)
    n_per_group = long_df["config"].value_counts().mean()
    tukey_se = np.sqrt(2 * mse / n_per_group)
    pairwise_p = pd.DataFrame(index=methods, columns=methods, data=1.0)
    for i, m1 in enumerate(methods):
        for m2 in methods[i + 1:]:
            diff = means[m1] - means[m2]
            studentized_range = abs(diff) / tukey_se
            p_adj = psturng(studentized_range * np.sqrt(2), n_groups, df_resid)
            p_adj = float(p_adj[0]) if isinstance(p_adj, np.ndarray) else float(p_adj)
            pairwise_p.loc[m1, m2] = pairwise_p.loc[m2, m1] = p_adj
    return pairwise_p.astype(float), means


def friedman_omnibus(long_df):
    res = pg.friedman(long_df, dv="value", within="config", subject="subject")
    return {"friedman_p": float(res["p_unc"].iloc[0])}


def _posthoc_conover_friedman(wide, p_adjust):
    # scikit-posthocs raises a benign RuntimeWarning from its own diagonal (self-vs-self)
    # computation -- confirmed in notebook 05 not to leak a NaN/Inf into the off-diagonal
    # p-values it returns, so it is suppressed rather than left to clutter the output.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning,
                                message="divide by zero encountered")
        return sp.posthoc_conover_friedman(wide, p_adjust=p_adjust)


def conover_friedman_holm(long_df):
    wide = long_df.pivot(index="subject", columns="config", values="value")
    pairwise_p = _posthoc_conover_friedman(wide, p_adjust="holm")
    return pairwise_p.astype(float), wide.mean()


def bh_pairwise(long_df, parametric):
    """BH-adjusted pairwise p-values, REPLACING Tukey HSD's / Conover-Friedman's own correction
    -- for >10-group comparisons. Uses the same underlying pairwise test as whichever branch the
    comparison is on (paired t for parametric, Conover-Friedman's raw p for nonparametric), just
    correcting with BH instead of the studentized range / Holm-Bonferroni.
    """
    wide = long_df.pivot(index="subject", columns="config", values="value")
    methods = list(wide.columns)
    pairs = [(m1, m2) for i, m1 in enumerate(methods) for m2 in methods[i + 1:]]
    if parametric:
        raw_p = [float(ttest_rel(wide[m1], wide[m2]).pvalue) for m1, m2 in pairs]
    else:
        raw = _posthoc_conover_friedman(wide, p_adjust=None)
        raw_p = [float(raw.loc[m1, m2]) for m1, m2 in pairs]
    _, adj_p, _, _ = multipletests(raw_p, method="fdr_bh")
    pairwise_p = pd.DataFrame(index=methods, columns=methods, data=1.0)
    for (m1, m2), p in zip(pairs, adj_p):
        pairwise_p.loc[m1, m2] = pairwise_p.loc[m2, m1] = float(p)
    return pairwise_p.astype(float), wide.mean()


def assign_cld(means, pairwise_p, alpha=0.05, ascending=True):
    """Compact letter display via a contiguous-run sweep on the mean-sorted order (Piepho
    2004-style insert/absorb). Correctness property: any two groups sharing a letter are
    confirmed by `pairwise_p` not to be significantly different. Validated against Ash SI
    Table 1's own worked example below before use.
    """
    order = sorted(means.index, key=lambda g: means[g], reverse=not ascending)
    n = len(order)

    def all_same(i, j):
        return all(
            pairwise_p.loc[order[a], order[b]] >= alpha
            for a in range(i, j + 1) for b in range(a + 1, j + 1)
        )

    runs = []
    i = 0
    while i < n:
        j = i
        while j + 1 < n and all_same(i, j + 1):
            j += 1
        runs.append((i, j))
        i += 1

    kept = sorted(set(r for r in runs
                      if not any(r != o and o[0] <= r[0] and r[1] <= o[1] for o in runs)))

    letters = {g: "" for g in order}
    for k, (i, j) in enumerate(kept):
        letter = chr(ord("a") + k)
        for idx in range(i, j + 1):
            letters[order[idx]] += letter
    return letters


# validate assign_cld against Ash SI Table 1's own worked example before trusting it below
def validate_assign_cld():
    """Notebook 05's own validation of `assign_cld` against Ash SI Table 1's worked example,
    verbatim -- but as a function. In notebooks 05 and 37 this ran at cell scope; as a module it
    would run on every import and print, so callers invoke it explicitly and notebook 46 does.
    Returns the letter assignment so a caller can show it."""
    means = pd.Series({"chemprop_st": 0.30, "chemprop_mt": 0.30, "lgbm_morgan": 0.37})
    pw = pd.DataFrame(1.0, index=means.index, columns=means.index)
    pw.loc["chemprop_st", "lgbm_morgan"] = pw.loc["lgbm_morgan", "chemprop_st"] = 0.001
    pw.loc["chemprop_mt", "lgbm_morgan"] = pw.loc["lgbm_morgan", "chemprop_mt"] = 0.001
    letters = assign_cld(means, pw, ascending=True)
    assert letters["chemprop_st"] == letters["chemprop_mt"] != letters["lgbm_morgan"], letters
    return letters


# ---------------------------------------------------------------- new in notebook 46
def branch_is_parametric(long_df, alpha=0.05):
    """Notebook 37's own branch rule, named rather than left inline: homogeneous iff Levene's p
    exceeds `alpha` AND the max/min variance ratio is inside `VARIANCE_RATIO_TOLERANCE`.

    Returns `(is_parametric, levene_result_dict)`. Measured per (comparison, metric), never
    inherited -- notebook 05 found ST-RAE nonparametric everywhere, notebook 29 parametric on all
    four isoforms, and notebook 37 split by isoform, so the three disagree and the check has to
    be run on the data in hand.
    """
    res = levene_check(long_df)
    homogeneous = (res["levene_p"] > alpha) and (res["variance_ratio"] <= VARIANCE_RATIO_TOLERANCE)
    return homogeneous, res


def posthoc_route(n_groups, parametric):
    """Which post-hoc Ash's protocol selects. Separated out so a design can assert which branch
    it expects to fire -- notebooks 05 and 37 both ran at >10 groups and so never exercised
    Tukey, and a design that intends to use it should check rather than assume."""
    if n_groups > LARGE_COMPARISON_THRESHOLD:
        return "bh"
    return "tukey" if parametric else "conover_holm"
