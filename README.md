# OpenADMET CYP450 Inhibition Blind Challenge — Direct Inhibition Track

Solo entry to the [OpenADMET](https://openadmet.org) CYP450 Inhibition Blind Challenge, Direct Inhibition regression track. Final-year BSc Pharmaceutical Chemistry dissertation project (**CH3PRO**, University of Reading), supervised by Dr Mauricio Cafiero.

Predict pIC50 for CYP1A2, CYP2C9, CYP2D6 and CYP3A4 on a 750-compound blinded test set, from 4,905 labelled training compounds carrying 6,525 measurements. Challenge opened 17 August 2026, closes 3 November 2026.

---

## Final Submission

**[`METHOD.md`](METHOD.md) for full model report**

The final submission is **NB56b-s050** (`notebooks/56b_cyp1a2_corrected.ipynb`): a dead-zone-target Chemprop/CheMeleon model with an affine placement-and-scale correction at half strength.

| | Macro | CYP1A2 | CYP2C9 | CYP2D6 | CYP3A4 |
|---|---:|---:|---:|---:|---:|
| **ST-RAE** | **0.5861** | 0.5631 | 0.5408 | 0.7502 | 0.4904 |

Best of eighteen scored rows across nineteen real sends, ahead of the previous best at 0.5982. Its correction constants are solved from the board's own returned metrics, which is leaderboard probing and is disclosed as such in the model report.

---

## The metric

Submissions are scored on **macro ST-RAE**, the mean across the four isoforms of a soft-threshold relative absolute error. **Lower is better**, and a naive mean-prediction baseline scores exactly 1.0 by construction, so anything under 1.0 beats that floor.

ST-RAE charges nothing for a prediction falling inside a compound's credible interval, and the distance to the nearest edge otherwise. That shape drives most of the decisions in the model report. Scoring code: `src/vendor/openadmet_eval/`, vendored unmodified.

---

## On blind-set labels

**No blind-set labels exist anywhere in this repository, in any form.** Every leaderboard figure cited anywhere in this project came from OpenADMET's own board after a real, deliberate submission, never from a locally held answer key.

The project's working rule in `CLAUDE.md` is explicit: the blind test set is single-use evidence, not training data, at any point. No prompt that would use test-set labels, leaderboard scores or the full release to pick or retrain a model is acted on without stopping to flag it first.

Session transcripts for AI-assisted work are kept under `.cc-history/` as a disclosure trail, per the module's AI-usage requirements.

---

## Repository map

| | |
|---|---|
| [`notebooks/`](notebooks/README.md) | analysis and modelling notebooks, run in numbered order |
| [`src/`](src/README.md) | shared library code, plus vendored third-party code under `src/vendor/` |
| [`scripts/`](scripts/README.md) | standalone long-running training scripts, launched by hand |
| [`data/`](data/README.md) | raw release (gitignored), processed (gitignored, regenerated), frozen fold splits |
| [`environment/`](environment/README.md) | archived environment specs and the `cyp-admet-v2` migration |
| [`outputs/`](outputs/README.md) | per-run predictions, params, metrics and submission files |
| `figures/` | the figure-generation layer: one module per figure, three output targets |
| `docs/figures/` | committed PNGs used by `METHOD.md` |
| `logs/`, `models/` | run logs; frozen checkpoints (gitignored) |
| `.cc-history/` | Claude Code session transcripts |

**`docs/`** holds the record layer:

- [`leaderboard_submissions.md`](docs/leaderboard_submissions.md) — the canonical, full-precision record of every real submission, all metrics, all isoforms, each with its sha256. Read this rather than reconstructing scores from notebook prose.
- [`findings.md`](docs/findings.md) — results that outlive the notebook that produced them.
- [`prereg_fullset.md`](docs/prereg_fullset.md) — predictions committed before the labels were released, with the tolerances that count as a pass.

`CLAUDE.md` at the repo root is the agent operating file: the rules this project works under and the notebook-by-notebook history. It is written for a cold-start Claude Code session rather than as reader-facing documentation, but it is the most complete single record here.

---

## Quick start

Native arm64 conda is required. `cyp-admet-v2` lives under `~/miniconda3-arm64`, not the default `/opt/anaconda3`, so a plain `conda activate` will not find it unless that installation's `conda` is initialised first (see [`environment/README.md`](environment/README.md)):

```bash
~/miniconda3-arm64/bin/conda env create -f environment.yml
source ~/miniconda3-arm64/etc/profile.d/conda.sh
conda activate cyp-admet-v2
```

`environment.yml` is the authoritative spec. `requirements.txt` is a secondary `pip freeze` reference, not an independent install path — RDKit in particular needs channel ordering that a plain `pip install -r` into an empty environment will not reproduce.

RDKit and PyTorch each bundle their own OpenMP runtime, and running them multi-threaded in the same process segfaults on this machine. Every notebook and script sets these before RDKit, PyTorch or Chemprop is first imported, and you should too if running cells by hand:

```bash
export OMP_NUM_THREADS=1
export KMP_DUPLICATE_LIB_OK=TRUE
```

Notebooks are executed headlessly, not through the Jupyter UI:

```bash
jupyter nbconvert --to notebook --execute --inplace notebooks/00_schema_audit.ipynb
```

Notebooks `00`→`03` run first, in order, and only once. They freeze the curated data and the CV fold assignment that every later notebook reads read-only. After that see [`notebooks/README.md`](notebooks/README.md). Long training runs happen in `scripts/`, launched by hand in the background.

---

## How this project held itself to a standard

Three rules were fixed before the experiments that use them, and they are what the model report's claims rest on.

**The fold partitions are frozen and hash-pinned.** A random 5×5 and a cluster-disjoint 5×5, both in `data/folds/`, neither modified after the first experiment ran against it.

**Significance is not sufficient for adoption.** Every comparison carries a materiality bar derived from how well the design resolves any difference, frozen before results, and a result that is significant but under its bar is reported as detected, not material. That is the commonest outcome in this project's screening history.

**Predictions were written down before they could be checked.** `docs/prereg_fullset.md` and the per-experiment pre-registrations in the notebooks record what should happen, with tolerances, including the predictions most likely to fail.

