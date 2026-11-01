"""Macro and per-isoform board ST-RAE across every submission, in true chronological order.

Method and provenance live here; the caption is ``CAPTION`` below, written to
``build/captions.md`` by ``build.py``.

ORDER IS NOT TABLE ORDER. The metric table in ``docs/leaderboard_submissions.md`` lists ``NB10``
(2026-09-03) *after* ``10c`` (2026-09-04), so its column order is not chronological. Order here
comes from that document's own timestamped table, via ``figures.board.chronological``.

Since NB40-ensemble's per-isoform row was recovered from the 2026-10-02 screenshots, every send
except the byte-identical re-send has a full set of cells, so there is no missing data to handle.

Every value is parsed from the record, never retyped.
"""

from __future__ import annotations

import pandas as pd

from figures.board import MACRO, chronological, metrics
from figures.style.base import Target, figure, figure_title, save, spread_labels
from figures.style.palette import ISOFORM_COLOURS, ISOFORMS, NEUTRAL, REFERENCE, TEXT

SLUG = "leaderboard_trajectory"
SIZES = {"thesis": (6.3, 3.8), "slides": (9.0, 5.0), "web": (8.0, 4.4)}

INTERIM_ORDER = 12          # the accidental 27b re-send, which the interim full set scored

CAPTION = """
Board ST-RAE (lower is better) for every real submission, in submission order, with the macro in
grey and heavier and the four isoforms in their own colours. Order is taken from the record's own
timestamped table, not from the metric table's column order, which is not chronological.

**The open markers at order 12 are the one byte-identical re-send of 27b**, fired accidentally by
a notebook re-run; every other mark is solid. Its values repeat order 10's exactly because the board returns the same
scores for the same bytes; that repetition is the point, not a plotting artefact, and it is the
sixth confirmation of board determinism on this record. It is also the send the interim full-set
reveal scored, which the labelled dashed rule marks.

The trajectory is not monotonic, and it does **not** end on its own best point. The best row on
record is NB56b-s050 at macro 0.5861, at order 18 — the first record-holder that is not a
per-isoform assembly. Order 19 is its full-strength twin, NB56b-s100, which scored 0.5933: the two
differ in one scalar, and sending both is what established that half strength was the right
setting rather than an assumption.
"""


def variant(target: Target) -> str:
    return "single-panel (macro + four isoforms)"


def series() -> pd.DataFrame:
    chrono = chronological()
    m = metrics()
    wide = m[m.metric == "ST-RAE"].pivot(index="submission", columns="isoform", values="value")
    rows = []
    for _, r in chrono.iterrows():
        row = {"order": r.order, "submission": r.submission, "accidental": r.accidental,
               MACRO: r.macro_strae}
        for iso in ISOFORMS:
            row[iso] = float(wide.loc[r.submission, iso])
        rows.append(row)
    out = pd.DataFrame(rows)
    # Every send now has a complete row -- NB40's was recovered on 2026-10-29.
    assert out[list(ISOFORMS)].notna().all().all(), "a submission is missing isoform cells"
    return out


def draw(target: Target):
    s = series()
    width, height = SIZES[target.name]

    with figure(target, width, height) as fig:
        ax = fig.add_subplot(111)
        x = s.order.to_numpy()
        resend = s.accidental.to_numpy()

        named = [(MACRO, "macro", NEUTRAL)] + [(i, i, ISOFORM_COLOURS[i]) for i in ISOFORMS]
        for key, _, colour in named:
            macro = key == MACRO
            y = s[key].astype(float)
            ax.plot(x, y, color=colour,
                    linewidth=target.base_size * (0.20 if macro else 0.09),
                    zorder=4 if macro else 3)
            # The re-send is drawn hollow: a real send whose values repeat by construction.
            ax.plot(x[~resend], y[~resend], linestyle="none", marker="o", color=colour,
                    markersize=target.base_size * (0.42 if macro else 0.30), zorder=5)
            ax.plot(x[resend], y[resend], linestyle="none", marker="o",
                    markerfacecolor="white", markeredgecolor=colour,
                    markeredgewidth=target.base_size * 0.07,
                    markersize=target.base_size * (0.46 if macro else 0.34), zorder=5)

        # The interim-reveal rule, with its label restored: without it the dashed line is
        # unexplained, and it marks the only full-set reading in the record.
        ax.axvline(INTERIM_ORDER, color=REFERENCE, linestyle="--",
                   linewidth=target.base_size * 0.08, zorder=1)
        ax.text(INTERIM_ORDER - 0.16, ax.get_ylim()[1], "interim full-set reveal ",
                rotation=90, va="top", ha="right", fontsize="small", color=TEXT)

        # In the right title slot, outside the data area entirely. Every in-axes position tried
        # sat on one series line or another -- the final five columns are crowded.
        best = s.loc[s[MACRO].idxmin()]
        ax.set_title(f"best macro: {best.submission} at {best[MACRO]:.4f}", loc="right",
                     fontsize="small", color=TEXT)

        ax.set_xticks(x)
        ax.set_xticklabels([f"{r.submission}{' (re-send)' if r.accidental else ''}"
                            for _, r in s.iterrows()], rotation=55, ha="right")
        ax.set_ylabel("Board ST-RAE")
        ax.yaxis.grid(True, zorder=0)
        ax.set_axisbelow(True)

        last = s.iloc[-1]
        lo, hi = ax.get_ylim()
        ys = spread_labels([float(last[k]) for k, _, _ in named], (hi - lo) * 0.062)
        for (key, text, colour), ly in zip(named, ys):
            ax.text(x.max() + 0.4, ly, text, va="center", ha="left", fontsize="small",
                    color=TEXT, fontweight="bold" if key == MACRO else "normal")
        ax.set_xlim(x.min() - 0.5, x.max() + 3.0)
        figure_title(fig, target, "Every submission, in the order it was sent")
        return fig


def render(target: Target):
    return save(draw(target), SLUG, target)
