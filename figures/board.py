"""The board record, parsed once.

`docs/leaderboard_submissions.md` is the authoritative record and is markdown, not a CSV, so this
is the one place in the figure layer that parses it. Every board figure reads through here rather
than re-implementing the parse, and no figure retypes a board number.

Two properties of the document that a naive parse gets wrong:

* **Negative R² uses a Unicode minus (U+2212), not a hyphen**, so `float()` raises on it.
* **The metric table's column order is NOT chronological.** `NB10` (2026-09-03) is listed after
  `10c` (2026-09-04). Submission order comes from the document's own chronological table, which
  states the timestamps; `chronological()` returns that order.

`NB40-ensemble` is deliberately absent from the metric table: only its macro ST-RAE, rank and
CYP2D6 ST-RAE were ever captured, and the board keeps one row per user, so the rest is
unrecoverable. It appears in `chronological()` with a macro value and no per-isoform cells.
"""

from __future__ import annotations

import re
from functools import lru_cache

import pandas as pd

from figures.style.base import repo_root
from figures.style.palette import ISOFORMS

MACRO = "MA (Overall)"
METRICS = ("ST-RAE", "MAE", "R²", "Spearman", "Kendall")

_ROW = re.compile(
    r"\|\s*(?:\*\*(?P<iso>.+?)\*\*)?\s*\|\s*(?P<met>ST-RAE|MAE|R²|Spearman|Kendall)\s*\|(?P<vals>.+)\|\s*$"
)
_CHRONO = re.compile(r"^\|\s*(?P<n>\d+)\s*\|\s*(?P<sub>.+?)\s*\|\s*(?P<when>.+?)\s*\|\s*(?P<macro>.+?)\s*\|\s*$")


def _doc() -> list[str]:
    return (repo_root() / "docs" / "leaderboard_submissions.md").read_text().splitlines()


def _clean(cell: str) -> float:
    return float(cell.strip().strip("*").replace("−", "-"))


@lru_cache(maxsize=1)
def metrics() -> pd.DataFrame:
    """Long table: submission, isoform, metric, value. Isoform includes the macro row."""
    doc = _doc()
    header = next(l for l in doc if l.startswith("| Isoform | Metric | 04b |"))
    subs = [c.strip() for c in header.strip("|").split("|")][2:]

    out, cur = [], None
    for line in doc[doc.index(header) + 2:]:
        if not line.startswith("|"):
            break
        m = _ROW.match(line)
        if not m:
            raise ValueError(f"unparsed metric-table row: {line[:70]}")
        if m.group("iso"):
            cur = m.group("iso")
        vals = m.group("vals").split("|")
        if len(vals) != len(subs):
            raise ValueError(f"{cur}/{m.group('met')}: {len(vals)} values for {len(subs)} columns")
        for sub, cell in zip(subs, vals):
            out.append({"submission": sub, "isoform": cur, "metric": m.group("met"),
                        "value": _clean(cell)})

    df = pd.DataFrame(out)
    # The document's own claims about its table, checked rather than trusted.
    assert set(df.isoform) == set(ISOFORMS) | {MACRO}, sorted(set(df.isoform))
    assert set(df.metric) == set(METRICS)
    assert len(df) == len(subs) * 5 * 5
    return df


@lru_cache(maxsize=1)
def chronological() -> pd.DataFrame:
    """Submission order as the record states it: order, submission, when, macro ST-RAE.

    Read from the document's own chronological table, NOT from metric-table column order, which is
    not chronological. Accidental re-sends are flagged rather than dropped.
    """
    doc = _doc()
    start = doc.index("| Order | Submission | Submitted (UTC) | Macro ST-RAE |")
    out = []
    for line in doc[start + 2:]:
        if not line.startswith("|"):
            break
        m = _CHRONO.match(line)
        if not m:
            raise ValueError(f"unparsed chronological row: {line[:70]}")
        name = m.group("sub").replace("**", "")
        # Strip the parenthetical BEFORE splitting on the comma: "NB56-identify (probe,
        # deliberately degraded)" has its comma inside the parentheses, so splitting first
        # leaves "NB56-identify (probe" and the name no longer matches the metric table.
        bare = re.sub(r"\s*\([^)]*\)", "", name).split(",")[0].strip()
        out.append({
            "order": int(m.group("n")),
            "submission": bare,
            "label": name,
            "when": m.group("when").replace("**", ""),
            "macro_strae": _clean(m.group("macro").split("—")[0]),
            "accidental": "accidental" in name,
            "in_metric_table": False,
        })
    df = pd.DataFrame(out).sort_values("order").reset_index(drop=True)
    df["in_metric_table"] = df.submission.isin(set(metrics().submission))
    assert df.order.tolist() == list(range(1, len(df) + 1))
    return df


def wide(metric: str, chrono: bool = True, isoforms_only: bool = True) -> pd.DataFrame:
    """One metric as submissions x isoforms, one row per submission, in submission order.

    The accidental 27b re-send (order 12) is dropped: it is byte-identical to order 10 and the
    board returned the same scores, so as a *row* it is a duplicate. It is a real send and
    :func:`chronological` keeps it, because a trajectory over time should show it.
    """
    df = metrics()
    w = df[df.metric == metric].pivot(index="submission", columns="isoform", values="value")
    w = w[list(ISOFORMS) + ([] if isoforms_only else [MACRO])]
    if chrono:
        c = chronological()
        order = [s for s in c[~c.accidental].submission if s in w.index]
        missing = set(w.index) - set(order)
        assert not missing, f"in the metric table but not the chronological table: {missing}"
        assert len(order) == len(w), f"{len(order)} ordered vs {len(w)} rows"
        w = w.loc[order]
    return w
