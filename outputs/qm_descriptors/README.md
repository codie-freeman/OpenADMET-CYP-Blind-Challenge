# QM descriptor run — COMPLETE

Produced by `scripts/run_qm_descriptors.py`. **Finished 2026-10-05 17:21 UTC**:
all **5,655 compounds** (4,905 curated train + 750 blind), **11,310 ORCA jobs, zero failures**.

Final log line: `stopped: queue exhausted | submitted 3050 | done 5655/5655 | failed 0 | 19.31 h elapsed`

This run produced a descriptor table and nothing else. **No 3D/QM descriptor has been added to any
feature set, no model has been trained on these, no ablation has been run, and no submission was
built, validated or sent.** Whether this block carries signal is a separate, later question.

## Files

| file | contents |
|---|---|
| `descriptors.csv` | 5,655 rows x 35 columns — one row per compound |
| `hirshfeld_charges.csv` | 252,237 rows — one per atom (`atom_index`, `element`, `charge`, `spin`) |
| `status.json` | final state; rewritten in place during the run, so it holds only the last snapshot |
| `decisions_taken.json` | decisions taken without explicit approval, rewritten at every launch |
| `restart_test.json` | the stop-start test record, including what was **not** verified |
| `completion_timeline.csv` | 5,655 rows — per-compound completion timestamp, duration, basis, gap, cumulative count. Extracted from the run logs because those are gitignored |
| `orca/{compound}/` | per-compound `.inp`, full `.out`, `.property.txt`, geometries. **8.8 GB, gitignored.** ORCA scratch (`.gbw`/`.densities`/`.tmp`) pruned after each compound's parse assertions passed; every `.out` kept in full so SCF convergence stays inspectable |
| `../../logs/qm_full_run.log` | full run log, millisecond-timestamped. **Gitignored** — the timestamps it carries are preserved in `completion_timeline.csv` (see Drift) |
| `../../logs/qm_supervisor.log` | supervisor activity |

## Method

RDKit ETKDGv3 conformer ensemble (10 embeddings, seed 42, MMFF94 minimised, lowest kept)
→ ORCA `NATIVE-GFN2-XTB Opt` geometry optimisation
→ ORCA `B3LYP D3 def2-SV(P) NORI Hirshfeld` single point.

That is qcMol's published level (B3LYP-D3/def2-SV(P)//GFN2-xTB), so their validation of the
workflow applies. `NORI` removes the RI approximation — see notebook 42 for why it is required
here. Multiplicity 1 throughout (notebook 41 established the population contains zero radicals);
formal charge read per compound from RDKit, never defaulted.

Concurrency was a pool of 6 concurrent **single-core** jobs, not one wide job. `%pal nprocs 1`
invokes no MPI, which eliminates the `MPI_ERR_ARG` ceiling notebook 42 measured at 936–1,018 basis
functions — `OCNT-2328942` (1,018 basis) completed normally with no fallback.

## Cost

| | |
|---|---|
| ORCA invocations | 11,310 |
| Total ORCA compute | **281.6 core-h** (11.7 core-days) |
| GFN2 geometry optimisation | 18.3 core-h — **6.5%**, mean 11.7 s |
| DFT single point | 263.3 core-h — **93.5%**, mean 167.6 s |
| RDKit conformer stage | 0.62 h, mean 0.40 s |
| Active wall time | **47.03 h** over 4 sessions (49.05 h calendar, 2.01 h paused) |
| Pool speedup | **5.99x on 6 workers — ~100% efficiency** |
| Per compound | mean 179.8 s, median 167.5 s, min 5.1 s, max 1328.0 s |
| Basis functions | 2,181,424 total; 80–1,018, mean 386 |

Refitted cost model over all 5,655 completions, as two stages:

    total_s  =  0.4 (conformer)  +  10.0 (GFN2 floor)  +  2.835e-05 * basis^2.61     MAPE 13.3%

## Validation

All checks passed; none is a self-report from `status.json`.

- **Completeness** — 5,655 unique compounds, 0 duplicates, every row `status == OK`. Train/blind
  checked against `train_inhibition_curated.csv` and `test_blinded_curated.csv` directly:
  **4,905 + 750, 0 missing, 0 extra.**
- **Convergence** — 0 DFT non-convergence, 0 GFN2 optimisation failures, 0 SCF-trouble flags across
  all 11,310 jobs. SCF mean 10.1 cycles, max 13.
- **No missing values** in any of the 14 key descriptor columns.
- **Physical sanity** — 0 HOMO ≥ LUMO violations, 0 non-positive gaps, gap range 1.883–7.777 eV
  (median 4.873), D3 dispersion stabilising in all cases, all total energies negative.
- **Charge conservation (the strongest check)** — Hirshfeld charges sum to each molecule's formal
  charge: **max |residual| 0.001094 e, mean 0.000143 e, 0 compounds above 0.01 e.** That is
  numerical-integration noise, and it confirms both the Hirshfeld parser and that per-compound
  charge reached ORCA rather than defaulting to neutral. The population holds 5,634 neutral,
  11 cations and 10 anions; all balance.
- **Hirshfeld coverage** — 252,237 atom rows, **0 mismatches** against `total_atoms_with_h`.

## Three corrections to notebook 42

A 15-compound timing study extrapolated three things wrongly. This run measures all three directly.

1. **The stage-dominance reversal is wrong.** Notebook 42 projected GFN2 at **66%** of population
   cost and DFT at 32%. Measured: **GFN2 6.5%, DFT 93.5%.** GFN2 *is* nearly size-independent as it
   said, but at a mean of **11.7 s, not 70 s** — most likely because notebook 42 timed it at
   `nprocs=6`, where MPI setup overhead dominates a job that cheap. DFT dominates at every size
   present. Notebook 42's figure 1, its suptitle and its caption all assert the reversal and should
   be read with this correction attached.
2. **The 166.8 h projection was 3.55x pessimistic** — actual 47.03 h. Partly the `nprocs=1`
   throughput gain notebook 42 measured itself (2.54x), partly its cost model being fitted on 15
   deliberately top-weighted compounds.
3. **The MPI ceiling never arose**, by construction — see Method.

## Drift: measured, and absent

Notebook 42 listed thermal and contention degradation over sustained multi-day load as explicitly
unquantified, with notebook 28's 10.6x overnight degradation on this machine as the live precedent.

Raw hourly throughput rose from 54 mol/h in the first full hour to a peak of 221, but that is mostly
the largest-first queue cheapening — the mean basis count of the compounds completing each hour fell
661 → 232 — so an uncontrolled reading is meaningless. Two readings, which agree:

- **Size-controlled regression** over all 5,655 compounds (`log(total_s) ~ log(basis) + active_hours`,
  session pauses removed) — drift coefficient **−0.0027 per active hour** (SE 0.00033, t = −8.1),
  i.e. **−11.9% over the full 47 h. The machine got slightly *faster* under sustained load.**
- **Banded median comparison** of early vs late halves — **0.91x** on the one basis band with enough
  overlap to support it. Largest-first ordering separates the halves almost completely, which is why
  the regression is the primary reading and the banding only a cross-check.

No degradation of any kind over 47 h of six-way load.

This is reconstructable **only** from `logs/qm_full_run.log`, whose `OK` lines carry millisecond
timestamps with compound name, duration and basis count. `descriptors.csv` has durations but no
wall-clock timestamp, and `status.json` is rewritten in place — hence the log being opened in append
mode so a same-day relaunch could not truncate it.

The repo gitignores every `.log` under `logs/`, **so that log is not preserved in version control.** The
per-compound completion times are therefore extracted into `completion_timeline.csv`, which is
tracked, covers all 5,655 compounds with 0 basis mismatches against `descriptors.csv`, and
reproduces the drift result on its own.

## Stop-start

Verified **six times**, and once at per-stage granularity, which is stronger than the per-compound
check originally specified: a hard kill mid-DFT left `OCNT-0496415` with a truncated job B, and the
restart skipped its conformer and GFN2 stages (mtimes unchanged) while re-running only the
truncated DFT stage.

The done-check is **content-based, never exit code or file existence** — notebook 42 established
ORCA exits 0 on error termination. A job counts as done only with `****ORCA TERMINATED NORMALLY****`
present, the not-fully-converged marker absent, and a `FINAL SINGLE POINT ENERGY` line.

Graceful SIGTERM was exercised twice: dispatch stops, in-flight jobs finish and are written, then
exit. `--hours` terminates through the same machinery, so a capped run's clean shutdown depends on
that path.

## Known gaps

1. **`descriptors.parquet` was never written** — neither `pyarrow` nor `fastparquet` is installed in
   `cyp-admet-v2`, so the write raises and is caught (`final parquet write skipped: Unable to find a
   usable engine`). Non-fatal and deliberate; `descriptors.csv` carries everything. To produce it:
   `conda install pyarrow`, then
   `python -c "import pandas as pd; pd.read_csv('outputs/qm_descriptors/descriptors.csv').to_parquet('outputs/qm_descriptors/descriptors.parquet', index=False)"`
2. **The specified pre-run test was aborted part-way**, on user instruction, and the clean
   all-5-cache-skip observation was never made (the user later waived it). `restart_test.json`
   records what was and was not verified.
3. **`scripts/qm_supervisor.sh` was added beyond the brief** to carry the run past a `--hours` cap
   the user lifted mid-campaign. It was tested by deliberately stopping a healthy run rather than
   trusted on first use. `run_qm_descriptors.py` itself was never modified.

## Re-running

Everything already on disk cache-skips, so a re-run is a near no-op and safe:

    python scripts/run_qm_descriptors.py --jobs 6

To regenerate from scratch, delete `outputs/qm_descriptors/orca/` first — the cache is keyed on
those `.out` files, not on `descriptors.csv`.
