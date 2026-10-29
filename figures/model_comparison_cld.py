"""Eleven modelling arms compared on 25-fold cluster-disjoint cross-validated ST-RAE.

Method and provenance live here; the caption is ``CAPTION`` below, which ``build.py``
writes to ``build/captions.md`` alongside the renders.
"""

from __future__ import annotations

from figures.data import load
from functools import lru_cache
from figures.style.base import Target, figure, figure_title, limit_ticks, save
from figures.style.palette import GRID, ISOFORM_COLOURS, ISOFORMS, NEUTRAL, REFERENCE, TEXT, arm_label

CAPTION = """
Twenty-five-fold cluster-disjoint cross-validated ST-RAE (lower is better) for
eleven modelling arms, by isoform; point and 95% confidence interval over 25 folds,
arms sorted best-first. Letters are compact letter displays: arms sharing a letter
are not significantly different within that isoform. The accented arm is each
panel's best. Panels do not share an x-scale -- each one's ST-RAE range is given on
its title line -- and the letters are within-isoform by construction, so arms must
not be compared across panels. Arm-name abbreviations are defined in the methods.

**The dead-zone target is uniquely best on all four isoforms here, and it has never
been tested on the blind board.** No submission carries the ``deadzone`` arm alone;
of these eleven arms only three (``plain``, ``aid``, ``mae_only``) were ever
submitted, and the arm holding the letter ``a`` in every panel is not among them.

The nearest board evidence applies the same clipped target to a DIFFERENT
architecture: NB32-A is the dead-zone target on the AID auxiliary-head model, a
combination notebook 32 deliberately never reconfirmed at 25 folds. It left CYP1A2,
CYP2C9 and CYP3A4 uncorrected, so those three read directly against 10c, the
submission of the ``plain`` arm, and board ST-RAE rose on all three: CYP1A2
0.6954 -> 0.7906, CYP2C9 0.5489 -> 0.5635, CYP3A4 0.4434 -> 0.4907. **That is
evidence against the target on that architecture, not against the arm shown here.**
The losses also differ in kind: CYP1A2's Spearman *improved* (0.7375 -> 0.7565) and
a later mean-preserving widening recovered most of its ST-RAE (0.7906 -> 0.7241,
NB36-widened), so that loss was largely placement and spread; CYP3A4's Spearman
declined (0.8177 -> 0.8020) and the same widening recovered nothing
(0.4907 -> 0.4904), so that loss is ranking, which no affine correction can reach.
NB32-A is the fourth of six recorded cross-validation-to-board transfer failures.

Provenance of every board figure quoted above: ``docs/leaderboard_submissions.md``
-- the per-isoform metric table (lines 34-56), the NB32-A section (line 637 on) and
the Spearman decomposition table (line 759 on). None is quoted from memory.

The caption is long because the figure is misleading without its second paragraph.
It splits cleanly into a short caption plus body text.

SOURCE
------
``outputs/37_butina_5x5_complete/figure7_forest_ST-RAE.csv`` -- mean, 95% CI, CLD
letter and tier per (comparison, arm), already aggregated over all 25 folds by
notebook 37. Nothing is recomputed here.

Note that the *per-fold* 25-fold scores are NOT all in notebook 37: its own
``run_scores.csv`` holds repeats 3 and 4 only, with repeats 0-1 in
``outputs/34_butina_5x5/`` and repeat 2 plus all of ``mae_only`` in
``outputs/35_butina_repeat3/``. This figure needs none of that, which is what keeps
it genuinely read-only.

MA (the macro average) is present in the source file and is dropped here: it is not
an isoform.
"""

SLUG = "model_comparison_cld"

#: Eleven rows at 16 pt need roughly 3.4 in of panel height, so a 2x2 grid needs
#: about 6.8 in -- against a 5.0 in slide content area. Rather than overflow, or
#: silently shrink the type below the target's own base size, the slides target
#: renders the one isoform a talk is about.
SLIDES_ISOFORM = "CYP2D6"

SIZES = {
    "thesis": (6.3, 6.0),
    "slides": (9.0, 5.0),
    "web": (8.0, 6.4),
}


def panels_for(target: Target) -> list[str]:
    if target.name == "slides":
        return [SLIDES_ISOFORM]
    return list(ISOFORMS)


def variant(target: Target) -> str:
    shown = panels_for(target)
    if len(shown) == 1:
        return f"single-panel ({shown[0]})"
    return f"2 x 2 grid ({len(shown)} isoforms)"


@lru_cache(maxsize=1)
def omnibus() -> dict:
    """Per-isoform omnibus test and p-value for ST-RAE, as notebook 37 actually ran it.

    The branch is NOT uniform: notebook 37 measured variance homogeneity per cell with Levene and
    routed accordingly, which sent **CYP2D6's ST-RAE down the nonparametric branch** while the other
    three stayed parametric. Labelling all four "ANOVA" would misstate the method, so the test name
    is read from the record rather than assumed.
    """
    t = load("37_butina_5x5_complete", "omnibus_tests.csv")
    t = t[(t.metric == "ST-RAE") & t.comparison.isin(ISOFORMS)]
    assert len(t) == 4, f"{len(t)} ST-RAE omnibus rows"
    assert t.significant.all(), "an omnibus did not clear its adjusted alpha"
    names = {"parametric": "RM-ANOVA", "nonparametric": "Friedman"}
    return {r.comparison: (names[r.branch], float(r.omnibus_p), int(r.n_groups))
            for _, r in t.iterrows()}


def table():
    """Per-arm aggregated ST-RAE, isoforms only, with the source's claims checked."""
    df = load("37_butina_5x5_complete", "figure7_forest_ST-RAE.csv")
    df = df[df["comparison"].isin(ISOFORMS)].copy()

    # The caption says "25 folds" and "the accented arm is each panel's best".
    # Both are claims about this file, so they are checked rather than trusted.
    folds = sorted(df["n_folds"].unique())
    if folds != [25]:
        raise ValueError(f"expected every row at n_folds=25, found {folds}")
    for iso, group in df.groupby("comparison"):
        if len(group) != 11:
            raise ValueError(f"{iso}: expected 11 arms, found {len(group)}")
        best = group.loc[group["tier"] == "best", "arm"].tolist()
        if len(best) != 1:
            raise ValueError(f"{iso}: expected exactly one best arm, found {best}")
        if group.loc[group["arm"] == best[0], "mean"].iat[0] != group["mean"].min():
            raise ValueError(f"{iso}: tier=='best' arm {best[0]} is not the lowest mean")
    return df


def _panel(ax, group, isoform: str, target: Target, show_xlabel: bool) -> None:
    # Sorted per panel, best at top. A single shared order across panels was tried and abandoned:
    # at eleven arms the compact letter displays stop running in sequence down the panel and
    # become a jumble, which costs more than cross-panel scanning gains.
    group = group.sort_values("mean").reset_index(drop=True)
    y = list(range(len(group)))

    lower = (group["mean"] - group["ci_low"]).to_numpy()
    upper = (group["ci_high"] - group["mean"]).to_numpy()
    accent = ISOFORM_COLOURS[isoform]
    is_best = (group["tier"] == "best").to_numpy()

    capsize = target.base_size * 0.22
    for i, best in enumerate(is_best):
        colour = accent if best else NEUTRAL
        ax.errorbar(
            group["mean"].iat[i], i,
            xerr=[[lower[i]], [upper[i]]],
            fmt="o", color=colour,
            ecolor=colour if best else REFERENCE,
            elinewidth=target.base_size * 0.09,
            capsize=capsize, capthick=target.base_size * 0.09,
            markersize=target.base_size * (0.52 if best else 0.42),
            zorder=4 if best else 3,
        )

    # The CLD letters sit just clear of the longest interval in THIS panel, at a
    # small offset past max(ci_high), rather than at a far-right position shared
    # across panels -- which left a wide gap after the error bars wherever a
    # panel's widest interval was not its rightmost point.
    lo = group["ci_low"].min()
    hi = group["ci_high"].max()
    span = hi - lo
    letter_x = hi + span * 0.05
    ax.set_xlim(lo - span * 0.06, letter_x + span * 0.13)   # room for two letters
    limit_ticks(ax, lo, hi)

    # The best arm's CI as two dashed rules spanning the panel: the tied group is whatever
    # overlaps the band, readable without decoding the letters.
    best = group[group.tier == "best"].iloc[0]
    # On CYP3A4 the two bounds sit about 0.015 apart, close enough that at full weight they read as
    # one thick rule under the marker. Lighten and thin them in proportion to how close they are.
    tight = (best.ci_high - best.ci_low) / span < 0.06
    for edge in (best.ci_low, best.ci_high):
        ax.axvline(edge, color=accent, linestyle="--",
                   alpha=0.30 if tight else 0.55,
                   linewidth=target.base_size * (0.045 if tight else 0.07), zorder=2)

    for i, letter in enumerate(group["CLD"]):
        ax.text(letter_x, i, letter,
                va="center", ha="left", fontsize="small", color=TEXT,
                fontweight="bold" if is_best[i] else "normal")

    ax.set_yticks(y)
    ax.set_yticklabels([arm_label(a) for a in group["arm"]])
    for i, best in enumerate(is_best):
        if best:
            ax.get_yticklabels()[i].set_fontweight("bold")

    ax.set_ylim(len(group) - 0.5, -0.5)      # shared arm order, best-pooled at the top

    # Each panel has its own x-scale, and a reader compares positions across
    # panels long before reaching the caption. The range goes on the title line,
    # where it is read before any position comparison -- an in-axes note was tried
    # and collides with the bottom row's marker. The numerals are NEUTRAL rather
    # than the isoform colour: coloured numerals cost legibility and add nothing
    # the title's own colour is not already saying.
    test, pval, _ = omnibus()[isoform]
    ax.set_title(isoform, color=accent)
    # The right title SLOT, not free text above the axes: a title is allowed to sit outside the
    # axes box, free text is not, and the check in Part 4 treats them differently for that reason.
    ax.set_title(f"{test} p < 0.001   {lo:.2f}-{hi:.2f}", loc="right",
                 fontsize="small", color=TEXT)
    ax.xaxis.grid(True, color=GRID, zorder=0)
    ax.set_axisbelow(True)
    if show_xlabel:
        ax.set_xlabel("ST-RAE")


def draw(target: Target):
    df = table()
    shown = panels_for(target)
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        if len(shown) == 1:
            axes = [fig.add_subplot(111)]
        else:
            axes = fig.subplots(2, 2).ravel().tolist()

        for ax, isoform in zip(axes, shown):
            _panel(ax, df[df["comparison"] == isoform], isoform, target,
                   show_xlabel=len(shown) == 1)

        if len(shown) > 1:
            fig.supxlabel(
                "ST-RAE, 25-fold cluster-disjoint cross-validation (lower is better)",
                fontsize="medium",
            )
        figure_title(fig, target, "The dead-zone target wins every panel, and was never submitted")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
