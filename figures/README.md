# `figures/`

A rendering layer. It reads the result CSVs the notebooks already wrote and renders them to three
targets (thesis PDF, slides PDF, web PNG) from one shared style definition. It is separate from the
notebooks on purpose: several of those are 12–17 hour runs and are the record of what was run.

```
python figures/build.py                         # every figure, every target
python figures/build.py --figure <slug>         # one figure
python figures/build.py --target thesis         # one target
```

Output lands in `build/`, which is gitignored: every file is regeneratable from tracked CSVs.

## Rules

1. **Figures read; they do not compute.** A number shown in a figure is a number a notebook already
   computed and saved. No figure refits a model, re-scores a prediction, or re-runs a test.
2. **No figure defines a colour.** Everything comes from `style/palette.py`.
3. **No figure imports from a notebook, opens one, or executes one.**
4. **Values are printed on the mark** wherever they fit without crowding, so a figure doubles as its
   own results table.
5. **Filenames are descriptive slugs**, never figure numbers — numbering belongs to the document.
6. Each figure script opens with a docstring holding its draft caption.

## The one exception to rule 1, and what qualifies

> A figure may evaluate a closed-form identity that a notebook has already solved, provided (a) the
> identity is written in the figure's docstring, (b) the notebook's own stored values are
> overplotted as verification, and (c) the source notebook is named. Anything else is computed
> upstream and read as a CSV.

Two further things look like computation and are not, because neither asserts anything the figure
does not already display:

* **A presentation aid derived from the plotted points themselves** — a trend line, an axis range, a
  rank, a count of points above a reference line. It restates the scatter; it reports no new
  quantity. `strae_vs_r2.py`'s least-squares line is this.
* **A join or a reshape** of stored tables, including deriving a range as `max − min` from stored
  values.

### Current exceptions, named

| figure | what it does | why it qualifies |
|---|---|---|
| `probe_branches.py` | evaluates `R² = 2ρk − k² − b(d)²` over a grid of translations | the identity is in the docstring; notebook 56's own `branch_predictions.csv` values are overplotted as a check; source named |
| `mean_range_by_isoform.py` | carries a four-entry pre-probe identification-status literal | a *verdict* from `docs/findings.md` finding 2, not derivable from the figure's own source file, and cited in place |
| `strae_vs_r2.py` | least-squares trend line | a presentation aid over the plotted points, per the clause above |

A new exception goes in this table, in the same change that introduces it.

## Layout

| path | role |
|---|---|
| `style/palette.py` | isoform colours (Okabe-Ito), neutrals, the sequential ramp, arm display names |
| `style/base.py` | `Target`/`TARGETS`, the `figure()` context manager, `save()`, shared helpers |
| `style/{thesis,slides,web}.mplstyle` | one style file per target |
| `data.py` | `load(notebook_dir, filename)` — the only door to `outputs/` |
| `board.py` | parses `docs/leaderboard_submissions.md` once; every board figure reads through it |
| `cv_board_pairs.py` | the cross-validation-to-board join, and the rules that make it small |
| `build.py` | renders everything; prints the variant each target got |

**Saving happens inside the target's style**, not merely alongside it: a `Text` artist created
without an explicit family stores the generic name `serif` and resolves it against
`rcParams["font.serif"]` at draw time, and `savefig` is a draw. `save()` re-enters the style so that
mistake cannot be made per-figure.

**`font.serif` must not name `Nimbus Roman` or `TeX Gyre Termes`**: neither is installed here and
matplotlib falls back to DejaVu *Sans* silently, which looks like a working serif chain and is not
one. The chain is Times New Roman → STIX Two Text → DejaVu Serif.
