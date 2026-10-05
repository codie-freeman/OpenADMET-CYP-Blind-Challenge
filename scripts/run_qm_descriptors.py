"""Production QM descriptor run over all 5,655 compounds (4,905 curated train + 750 blind).

Per compound: RDKit ETKDGv3 conformer ensemble (10 confs, MMFF94 minimised, lowest kept)
-> ORCA `! NATIVE-GFN2-XTB Opt` geometry optimisation -> ORCA
`! B3LYP D3 def2-SV(P) NORI Hirshfeld` single point on that geometry. Exactly the pipeline
notebook 42 timed, with the same keywords and the same output parsers, carried over verbatim.

This is a SCRIPT, not a notebook, because it runs across many stop-start sessions over days --
this project's convention is that `notebooks/` report and `scripts/` compute.

USAGE
    cd /path/to/OpenADMET-CYP-Blind-Challenge
    caffeinate -i nohup /Users/codiefreeman/miniconda3-arm64/envs/cyp-admet-v2/bin/python \\
        scripts/run_qm_descriptors.py --jobs 6 --hours 40 \\
        > logs/qm_$(date +%F).log 2>&1 &
    # then, in the morning:
    cat outputs/qm_descriptors/status.json

    # a 5-compound smoke test first, which is what was actually run before any long session:
    python scripts/run_qm_descriptors.py --jobs 2 --limit 5

ENVIRONMENT FACTS, every one established by notebook 42 running ORCA rather than by reading
documentation. Each is re-asserted at startup by `preflight()` so a changed install fails loudly
instead of producing wrong numbers.

  * The geometry keyword is `! NATIVE-GFN2-XTB Opt`. `! XTB2 Opt` FAILS on this installation --
    it is the legacy external-driver keyword and shells out to `otool_xtb` or a standalone `xtb`,
    neither of which is present. No xtb is installed and none is needed.
  * The DFT line needs `NORI`. Without it ORCA turns RIJCOSX on by default and MPI-parallel runs
    abort with MPI_ERR_ARG inside the RI/J Cholesky decomposition. `NORI` *removes* the RI
    approximation, so the single point computes exact Coulomb and exchange; notebook 42 measured
    the effect on HOMO/LUMO/gap at <= 0.0019 eV.
  * ORCA EXITS 0 ON ERROR TERMINATION. The exit code is never used as a success signal anywhere
    in this script. A job counts as done only when its .out contains `****ORCA TERMINATED
    NORMALLY****` AND a `FINAL SINGLE POINT ENERGY` line that is NOT the
    `(Wavefunction not fully converged!)` variant.
  * Formal charge is read per compound from RDKit, never hardcoded -- 21 training compounds are
    non-neutral. Multiplicity is 1 throughout, which notebook 41 licensed by establishing there
    are zero radicals in either set; the script asserts that rather than assuming it.

CONCURRENCY
Notebook 42 measured 6-core parallel DFT at 72% efficiency (4.33x on 6 cores) and found that
`floor(11 cores / nprocs) / t(nprocs)` peaks at nprocs=1: eleven concurrent single-core jobs
finish the set 2.54x faster than one 6-core job at a time. This script therefore dispatches a
POOL of `--jobs` concurrent jobs at `%pal nprocs 1 end`, not one wide job at a time.

Default `--jobs 6`: the same thermal load as the previous approach, roughly double the throughput,
and 6 x `%maxcore 3000` = 18 GB against 36 GB installed. `--jobs 11` would be faster still but
needs 33 GB and has not been measured; notebook 42 flagged it as reported-not-recommended.

A free consequence of running single-core: `%pal nprocs 1` invokes no MPI at all, so the
MPI_ERR_ARG ceiling notebook 42 measured (parallel DFT fails somewhere between 936 and 1,018
basis functions) CANNOT arise. Exactly one of the 5,655 compounds sits above that ceiling. The
serial-fallback machinery is kept as a safety net but is expected never to fire.

STOP-START
The run is designed to be killed and resumed indefinitely.
  * Cache-skip is per JOB, not per compound: stage 1, job A and job B are each checked
    independently, so an interruption between A and B costs nothing on resume.
  * Every result is written to disk the moment it completes. Nothing is buffered.
  * `--hours H` stops cleanly BETWEEN jobs once the budget is spent, so no calculation is ever
    killed mid-flight.
  * SIGINT/SIGTERM stop new dispatch, let in-flight jobs finish, write state and exit 0.
  * The queue is sorted LARGEST FIRST by an exact basis-function count, so the expensive
    compounds get the early hours of a session and a bad cost estimate surfaces immediately
    rather than on the last day.

READ-ONLY: `data/`, `src/`, `notebooks/`, and `outputs/` for notebooks 30-42. This script writes
only under `outputs/qm_descriptors/`. `notebooks/32_deadzone_aid_retrain.ipynb` and
`notebooks/40_single_object_ensemble.ipynb` are never opened -- both carry live submission code.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import signal
import subprocess
import sys
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

# ---------------------------------------------------------------- configuration

ORCA_BIN = Path("/Users/codiefreeman/.local/bin/orca")
GFN2_KEYWORDS = "NATIVE-GFN2-XTB Opt"
DFT_KEYWORDS = "B3LYP D3 def2-SV(P) NORI Hirshfeld"
MAXCORE_MB = 3000
N_CONFORMERS = 10
CONFORMER_SEED = 42
MULTIPLICITY = 1

CURATED = REPO_ROOT / "data" / "processed" / "train_inhibition_curated.csv"
BLIND = REPO_ROOT / "data" / "processed" / "test_blinded_curated.csv"

OUT = REPO_ROOT / "outputs" / "qm_descriptors"
ORCA_DIR = OUT / "orca"
DESCRIPTORS_CSV = OUT / "descriptors.csv"
DESCRIPTORS_PARQUET = OUT / "descriptors.parquet"
HIRSHFELD_CSV = OUT / "hirshfeld_charges.csv"
FAILURES_CSV = OUT / "failures.csv"
STATUS_JSON = OUT / "status.json"
DECISIONS_JSON = OUT / "decisions_taken.json"

TERMINATED_OK = "****ORCA TERMINATED NORMALLY****"
NOT_CONVERGED = "(Wavefunction not fully converged!)"

# Contracted def2-SV(P) basis functions per element. NOT a regression and NOT taken from a table:
# solved by least squares against notebook 42's 15 measured compounds (exact additivity, max
# residual 0.000000) and completed for Br and P with two dedicated single-point jobs, since
# neither element appeared in that selection with enough independence to be determined. Verified
# to reproduce all 15 measured counts exactly and to cover every element present in the
# population. Br (32) exceeds I (26) because iodine carries an ECP in def2 while bromine is
# all-electron -- the ordering is physical, not a transcription slip.
BASIS_PER_ELEMENT = {"H": 2, "C": 14, "N": 14, "O": 14, "F": 14,
                     "P": 18, "S": 18, "Cl": 18, "Br": 32, "I": 26}

# Scratch ORCA leaves behind that is not needed once the output is parsed. Every .out is kept IN
# FULL so SCF behaviour stays inspectable. On notebook 42's 15 compounds this cut 455 MB to 68 MB;
# across 5,655 compounds the wavefunction files alone would be well over 100 GB.
KEEP_SUFFIXES = (".inp", ".out", ".xyz", ".property.txt")

_LOCK = threading.Lock()
_STOP = threading.Event()
DECISIONS: list[dict] = []
log = logging.getLogger("qm")


def record_decision(decision: str, reason: str, alternatives: str, authority: str) -> None:
    DECISIONS.append({"decision": decision, "reason": reason, "alternatives": alternatives,
                      "authority": authority, "at_utc": datetime.now(timezone.utc).isoformat()})
    DECISIONS_JSON.write_text(json.dumps(DECISIONS, indent=2))


# ---------------------------------------------------------------- ORCA output parsers
# Carried over verbatim from notebook 42, where every anchor below was validated against real
# ORCA 6.1.1 output before anything depended on it. Each raises rather than returning NaN, so a
# format surprise surfaces immediately instead of producing a silently blank column.


class OrcaParseError(RuntimeError):
    pass


def final_energy(txt: str) -> tuple[float, bool]:
    hits = re.findall(r"^FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)(.*)$", txt, re.M)
    if not hits:
        raise OrcaParseError("no FINAL SINGLE POINT ENERGY line")
    energy, trailer = hits[-1]
    return float(energy), ("not fully converged" in trailer.lower())


def n_basis_functions(txt: str) -> int:
    # Anchored to line start so it cannot match "   # of basis functions in Aux-J   ...   71".
    m = re.search(r"^Number of basis functions\s+\.\.\.\s+(\d+)\s*$", txt, re.M)
    if not m:
        raise OrcaParseError("no 'Number of basis functions' line")
    return int(m.group(1))


def scf_cycles(txt: str):
    hits = re.findall(r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES", txt)
    return int(hits[-1]) if hits else None


def scf_trouble(txt: str) -> list[str]:
    # "Aborting the run." is deliberately absent: ORCA's optimiser prints it for a RECOVERABLE
    # augmented-Hessian event and then retries the step in Cartesians and finishes normally.
    # Notebook 42 shipped that false positive once and caught it.
    flags = []
    for pat, label in [(r"SCF NOT CONVERGED", "SCF NOT CONVERGED"),
                       (r"Wavefunction not fully converged", "wavefunction not fully converged"),
                       (r"SCF CONVERGENCE PROBLEM", "SCF convergence problem")]:
        if re.search(pat, txt, re.I):
            flags.append(label)
    return flags


def orbital_energies(txt: str) -> tuple[float, float, float]:
    blocks = txt.split("ORBITAL ENERGIES")
    if len(blocks) < 2:
        raise OrcaParseError("no ORBITAL ENERGIES block")
    occ, vir = [], []
    for line in blocks[-1].split("\n"):
        m = re.match(r"^\s*(\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", line)
        if m:
            (occ if float(m.group(2)) > 0.5 else vir).append(float(m.group(4)))
    if not occ or not vir:
        raise OrcaParseError(f"ORBITAL ENERGIES found but occ={len(occ)} vir={len(vir)}")
    homo, lumo = max(occ), min(vir)
    return homo, lumo, lumo - homo


def hirshfeld_charges(txt: str) -> list[dict]:
    if "HIRSHFELD ANALYSIS" not in txt:
        raise OrcaParseError("no HIRSHFELD ANALYSIS block")
    rows, started = [], False
    for line in txt.split("HIRSHFELD ANALYSIS")[-1].split("\n"):
        m = re.match(r"^\s*(\d+)\s+([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", line)
        if m:
            started = True
            rows.append({"atom_index": int(m.group(1)), "element": m.group(2),
                         "charge": float(m.group(3)), "spin": float(m.group(4))})
        elif started and line.strip().startswith("TOTAL"):
            break
    if not rows:
        raise OrcaParseError("HIRSHFELD ANALYSIS found but no charge rows parsed")
    return rows


def dipole_debye(txt: str) -> float:
    m = re.findall(r"^Magnitude \(Debye\)\s+:\s+(-?\d+\.\d+)\s*$", txt, re.M)
    if not m:
        raise OrcaParseError("no 'Magnitude (Debye)' line")
    return float(m[-1])


def dispersion_correction(txt: str):
    m = re.findall(r"^Dispersion correction\s+(-?\d+\.\d+)\s*$", txt, re.M)
    return float(m[-1]) if m else None


def orca_wall_seconds(txt: str) -> float:
    m = re.search(r"TOTAL RUN TIME:\s+(\d+) days (\d+) hours (\d+) minutes (\d+) seconds "
                  r"(\d+) msec", txt)
    if not m:
        raise OrcaParseError("no TOTAL RUN TIME line")
    d, h, mi, s, ms = (int(g) for g in m.groups())
    return d * 86400 + h * 3600 + mi * 60 + s + ms / 1000.0


def optimisation_converged(txt: str) -> bool:
    return "*** OPTIMIZATION RUN DONE ***" in txt


def read_orca_xyz(path: Path) -> list[str]:
    lines = Path(path).read_text().rstrip("\n").split("\n")
    n = int(lines[0].strip())
    body = lines[2:2 + n]
    if len(body) != n:
        raise OrcaParseError(f"{path}: expected {n} atoms, got {len(body)}")
    return body


# ---------------------------------------------------------------- job state


def job_is_done(out_path: Path) -> bool:
    """A job counts as done only on all three conditions, per the brief and notebook 42.

    ORCA exits 0 on error termination, so neither the exit code nor the mere existence of a .out
    means anything. A truncated .out from a killed job, or one whose SCF did not converge, is
    treated as not-done and re-run.
    """
    if not out_path.exists():
        return False
    try:
        txt = out_path.read_text(errors="replace")
    except OSError:
        return False
    if TERMINATED_OK not in txt:
        return False
    if NOT_CONVERGED in txt:
        return False
    return bool(re.search(r"^FINAL SINGLE POINT ENERGY\s+-?\d+\.\d+", txt, re.M))


def write_inp(path: Path, keywords: str, nprocs: int, charge: int, geom: list[str]) -> None:
    body = [f"! {keywords}"]
    if nprocs > 1:
        body.append(f"%pal nprocs {nprocs} end")
    body.append(f"%maxcore {MAXCORE_MB}")
    body.append(f"* xyz {charge} {MULTIPLICITY}")
    body.extend(geom)
    body.append("*")
    path.write_text("\n".join(body) + "\n")


def _kill_group(proc: subprocess.Popen, tag: str) -> None:
    # ORCA's children survive a plain terminate() and would hold cores for the rest of the run.
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        log.error("%s: TIMEOUT, sent SIGTERM to the process group", tag)
        time.sleep(15)
        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        pass


def run_orca(tag: str, keywords: str, geom: list[str], charge: int, workdir: Path,
             nprocs: int, timeout_s: float) -> tuple[Path, float]:
    workdir.mkdir(parents=True, exist_ok=True)
    inp, out = workdir / f"{tag}.inp", workdir / f"{tag}.out"
    write_inp(inp, keywords, nprocs, charge, geom)
    t0 = time.time()
    # cwd is the job directory so ORCA resolves its input and drops all scratch there;
    # start_new_session gives the job its own process group so a timeout can kill its children.
    proc = subprocess.Popen([str(ORCA_BIN), inp.name], cwd=str(workdir),
                            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, bufsize=1,
                            start_new_session=True)
    timer = threading.Timer(timeout_s, _kill_group, args=(proc, tag))
    timer.start()
    try:
        with open(out, "w") as fh:
            for line in proc.stdout:
                fh.write(line)
        proc.wait()
    finally:
        timer.cancel()
    return out, time.time() - t0


def prune_scratch(workdir: Path) -> int:
    freed = 0
    for f in workdir.iterdir():
        if f.is_file() and not any(f.name.endswith(k) for k in KEEP_SUFFIXES):
            try:
                freed += f.stat().st_size
                f.unlink()
            except OSError:
                pass
    return freed


# ---------------------------------------------------------------- the per-compound pipeline


def estimate_basis(smiles: str) -> int:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    return sum(BASIS_PER_ELEMENT[a.GetSymbol()] for a in mol.GetAtoms())


def stage1_conformer(name: str, smiles: str, n_atoms_expected: int, charge_expected: int,
                     workdir: Path) -> tuple[list[str], float, bool]:
    xyz = workdir / f"{name}__stage1.xyz"
    if xyz.exists():
        try:
            body = read_orca_xyz(xyz)
            if len(body) == n_atoms_expected:
                return body, 0.0, True
        except (OrcaParseError, ValueError, OSError):
            pass
    t0 = time.time()
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = CONFORMER_SEED
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=N_CONFORMERS, params=params)
    if not cids:
        raise RuntimeError("ETKDG embedded no conformers")
    res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)
    best = min(range(len(res)), key=lambda i: res[i][1])
    if mol.GetNumAtoms() != n_atoms_expected:
        raise RuntimeError(f"{mol.GetNumAtoms()} atoms, expected {n_atoms_expected}")
    if Chem.GetFormalCharge(mol) != charge_expected:
        raise RuntimeError(f"charge {Chem.GetFormalCharge(mol)}, expected {charge_expected}")
    conf = mol.GetConformer(cids[best])
    body = []
    for i, at in enumerate(mol.GetAtoms()):
        p = conf.GetAtomPosition(i)
        body.append(f"  {at.GetSymbol():2s} {p.x:18.12f} {p.y:18.12f} {p.z:18.12f}")
    workdir.mkdir(parents=True, exist_ok=True)
    xyz.write_text(f"{len(body)}\nETKDGv3(seed={CONFORMER_SEED}) + MMFF94, conformer {best} of "
                   f"{len(cids)}, E={res[best][1]:.4f}\n" + "\n".join(body) + "\n")
    return body, time.time() - t0, False


def run_one_compound(row, timeout_a: float, timeout_b: float) -> dict:
    """Stage 1 -> job A -> job B for one compound, each stage independently cache-skipped.

    Returns a result dict. Never raises: a failure is returned as a row so the pool continues.
    """
    name, smiles = row.Molecule_Name, row.SMILES
    charge, n_atoms = int(row.formal_charge), int(row.total_atoms_with_h)
    wd = ORCA_DIR / name
    rec = {"Molecule_Name": name, "set": row.set, "inchikey": row.inchikey,
           "heavy_atoms": int(row.heavy_atoms), "total_atoms_with_h": n_atoms,
           "formal_charge": charge, "n_basis_estimated": int(row.n_basis_est),
           "status": "OK", "error": "", "cached": False}
    t_start = time.time()
    try:
        geom1, t1, cached1 = stage1_conformer(name, smiles, n_atoms, charge, wd)
        rec["stage1_s"] = t1

        # ---- job A: GFN2 geometry optimisation
        tag_a = f"{name}__A_gfn2opt"
        out_a = wd / f"{tag_a}.out"
        if job_is_done(out_a):
            rec["jobA_s"], cached_a = 0.0, True
        else:
            out_a, wall_a = run_orca(tag_a, GFN2_KEYWORDS, geom1, charge, wd, 1, timeout_a)
            rec["jobA_s"], cached_a = wall_a, False
            if not job_is_done(out_a):
                raise RuntimeError("job A did not complete (see .out)")
        txt_a = out_a.read_text(errors="replace")
        rec["gfn2_opt_converged"] = optimisation_converged(txt_a)
        h, l, g = orbital_energies(txt_a)
        rec.update({"gfn2_homo_ev": h, "gfn2_lumo_ev": l, "gfn2_gap_ev": g,
                    "gfn2_energy_eh": final_energy(txt_a)[0],
                    "gfn2_orca_wall_s": orca_wall_seconds(txt_a)})
        geom_opt = read_orca_xyz(out_a.with_suffix(".xyz"))

        # ---- job B: DFT single point on job A's optimised geometry
        tag_b = f"{name}__B_dftsp"
        out_b = wd / f"{tag_b}.out"
        if job_is_done(out_b):
            rec["jobB_s"], cached_b = 0.0, True
        else:
            out_b, wall_b = run_orca(tag_b, DFT_KEYWORDS, geom_opt, charge, wd, 1, timeout_b)
            rec["jobB_s"], cached_b = wall_b, False
            if not job_is_done(out_b):
                raise RuntimeError("job B did not complete (see .out)")
        txt_b = out_b.read_text(errors="replace")
        e, warn = final_energy(txt_b)
        h, l, g = orbital_energies(txt_b)
        ch = hirshfeld_charges(txt_b)
        if len(ch) != n_atoms:
            raise RuntimeError(f"{len(ch)} Hirshfeld rows for {n_atoms} atoms")
        if abs(sum(c["charge"] for c in ch) - charge) > 5e-3:
            raise RuntimeError("Hirshfeld charges do not sum to the formal charge")
        rec.update({
            "dft_energy_eh": e, "dft_not_fully_converged": warn,
            "dft_homo_ev": h, "dft_lumo_ev": l, "dft_gap_ev": g,
            "dft_dipole_debye": dipole_debye(txt_b),
            "dft_dispersion_eh": dispersion_correction(txt_b),
            "n_basis": n_basis_functions(txt_b),
            "dft_scf_cycles": scf_cycles(txt_b),
            "dft_scf_trouble": "|".join(scf_trouble(txt_b)),
            "dft_orca_wall_s": orca_wall_seconds(txt_b),
            "hirshfeld_min": min(c["charge"] for c in ch),
            "hirshfeld_max": max(c["charge"] for c in ch),
            "hirshfeld_abs_mean": float(np.mean([abs(c["charge"]) for c in ch])),
        })
        rec["cached"] = cached1 and cached_a and cached_b
        rec["freed_bytes"] = prune_scratch(wd)
        rec["hirshfeld"] = [{**c, "Molecule_Name": name} for c in ch]
    except Exception as exc:  # noqa: BLE001 -- a failure must not stop the pool
        rec["status"] = "FAILED"
        rec["error"] = f"{type(exc).__name__}: {exc}"
    rec["total_s"] = time.time() - t_start
    return rec


# ---------------------------------------------------------------- persistence


def append_rows(path: Path, rows: list[dict], key: str = "Molecule_Name") -> None:
    """Append to CSV, de-duplicating on `key`. Written on every completion, never buffered."""
    if not rows:
        return
    new = pd.DataFrame(rows)
    if path.exists():
        try:
            old = pd.read_csv(path)
            new = pd.concat([old[~old[key].isin(new[key])], new], ignore_index=True)
        except (pd.errors.EmptyDataError, OSError):
            pass
    tmp = path.with_suffix(path.suffix + ".tmp")
    new.to_csv(tmp, index=False)
    tmp.replace(path)          # atomic: a kill mid-write cannot corrupt the real file


def write_status(total: int, done_names: set, failed: list[dict], durations: list[float],
                 started_at: float, args) -> dict:
    remaining = total - len(done_names)
    mean_s = float(np.mean(durations)) if durations else None
    # Projection assumes the pool stays full, which it does except at the very end of the queue.
    proj = None
    if mean_s and remaining:
        proj_s = remaining * mean_s / max(1, args.jobs)
        proj = (datetime.now(timezone.utc) + timedelta(seconds=proj_s)).isoformat()
    status = {
        "updated_at_utc": datetime.now(timezone.utc).isoformat(),
        "total_compounds": total,
        "done": len(done_names),
        "remaining": remaining,
        "failed": len(failed),
        "pct_done": round(100.0 * len(done_names) / total, 2) if total else None,
        "jobs": args.jobs,
        "mean_seconds_per_compound": round(mean_s, 1) if mean_s else None,
        "completed_this_session": len(durations),
        "session_elapsed_h": round((time.time() - started_at) / 3600.0, 3),
        "hours_budget": args.hours,
        "projected_completion_utc": proj,
        "projected_remaining_h": round(remaining * mean_s / max(1, args.jobs) / 3600.0, 2)
                                 if mean_s and remaining else None,
        "failed_compounds": [f["Molecule_Name"] for f in failed][:50],
        "dft_keywords": DFT_KEYWORDS,
        "gfn2_keywords": GFN2_KEYWORDS,
        "nprocs_per_job": 1,
    }
    tmp = STATUS_JSON.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(status, indent=2))
    tmp.replace(STATUS_JSON)
    return status


# ---------------------------------------------------------------- startup checks


def preflight() -> dict:
    """Re-assert every environment fact this script depends on. Fails loudly."""
    env = {}
    if not ORCA_BIN.exists():
        raise SystemExit(f"ORCA launcher not found at {ORCA_BIN}")
    real = Path.home() / "orca" / "orca"
    if not real.exists():
        raise SystemExit(f"ORCA binary not found at {real}")
    ver = subprocess.run(["strings", str(real)], capture_output=True, text=True).stdout
    m = re.search(r"Program Version\s+([0-9]+\.[0-9]+\.[0-9]+)", ver)
    env["orca_version"] = m.group(1) if m else None
    if env["orca_version"] != "6.1.1":
        raise SystemExit(f"expected ORCA 6.1.1, found {env['orca_version']} -- notebook 42's "
                         f"keyword and parser findings were established against 6.1.1")
    import rdkit
    env["rdkit"] = rdkit.__version__
    if rdkit.__version__ != "2026.03.3":
        raise SystemExit(f"expected the pinned rdkit 2026.3.3, found {rdkit.__version__}")
    # The legacy keyword must still be absent, or the environment has changed under us.
    if (Path.home() / "orca" / "otool_xtb").exists():
        log.warning("otool_xtb now exists -- the XTB2 keyword may work again, but this script "
                    "deliberately keeps NATIVE-GFN2-XTB, which notebook 42 validated.")
    env["physical_cores"] = int(subprocess.run(["sysctl", "-n", "hw.physicalcpu"],
                                               capture_output=True, text=True).stdout.strip())
    return env


def load_population(limit: int | None) -> pd.DataFrame:
    train = pd.read_csv(CURATED)[["Molecule_Name", "SMILES", "inchikey"]].assign(set="train")
    blind = pd.read_csv(BLIND)[["Molecule_Name", "SMILES", "inchikey"]].assign(set="blind")
    pop = pd.concat([train, blind], ignore_index=True)
    if len(pop) != 5655:
        raise SystemExit(f"expected 5,655 compounds, got {len(pop)}")
    if not pop.Molecule_Name.is_unique:
        raise SystemExit("Molecule_Name is not unique across the two sets")

    heavy, atoms, charges, radicals, basis = [], [], [], [], []
    for s in pop.SMILES:
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            raise SystemExit(f"SMILES failed to parse: {s}")
        molh = Chem.AddHs(mol)
        heavy.append(mol.GetNumHeavyAtoms())
        atoms.append(molh.GetNumAtoms())
        charges.append(Chem.GetFormalCharge(mol))
        radicals.append(sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()))
        basis.append(sum(BASIS_PER_ELEMENT[a.GetSymbol()] for a in molh.GetAtoms()))
    pop["heavy_atoms"], pop["total_atoms_with_h"] = heavy, atoms
    pop["formal_charge"], pop["n_basis_est"] = charges, basis
    # Multiplicity 1 is hardcoded, so the absence of radicals is asserted, not assumed.
    if max(radicals) != 0:
        raise SystemExit("a radical-bearing compound is present -- multiplicity 1 is not safe")

    # LARGEST FIRST: the expensive compounds get the early hours of a session, so a bad cost
    # estimate shows up on day one rather than on the last day.
    pop = pop.sort_values(["n_basis_est", "Molecule_Name"], ascending=[False, True])
    pop = pop.reset_index(drop=True)
    return pop.head(limit) if limit else pop


def already_done(pop: pd.DataFrame) -> set:
    """Compounds whose job A and job B are both already complete on disk."""
    done = set()
    for name in pop.Molecule_Name:
        wd = ORCA_DIR / name
        if job_is_done(wd / f"{name}__A_gfn2opt.out") and job_is_done(wd / f"{name}__B_dftsp.out"):
            done.add(name)
    return done


# ---------------------------------------------------------------- main


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=6,
                    help="concurrent single-core ORCA jobs (default 6; 6 x 3 GB = 18 GB)")
    ap.add_argument("--hours", type=float, default=None,
                    help="stop cleanly BETWEEN jobs once this much wall time has elapsed")
    ap.add_argument("--limit", type=int, default=None,
                    help="process only the first N compounds of the largest-first queue")
    ap.add_argument("--timeout-a-min", type=float, default=90.0,
                    help="per-job timeout for the GFN2 optimisation")
    ap.add_argument("--timeout-b-min", type=float, default=180.0,
                    help="per-job timeout for the DFT single point")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    ORCA_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                        stream=sys.stdout)
    RDLogger.DisableLog("rdApp.info")

    started_at = time.time()
    env = preflight()
    log.info("ORCA %s | rdkit %s | %d physical cores | --jobs %d | nprocs 1 per job",
             env["orca_version"], env["rdkit"], env["physical_cores"], args.jobs)
    if args.jobs > env["physical_cores"]:
        log.warning("--jobs %d exceeds %d physical cores; jobs will contend",
                    args.jobs, env["physical_cores"])
    log.info("memory commitment: %d x %d MB = %.1f GB",
             args.jobs, MAXCORE_MB, args.jobs * MAXCORE_MB / 1024)

    record_decision(
        decision="Dispatched a pool of single-core ORCA jobs rather than one wide job at a time",
        reason="Notebook 42 measured 6-core DFT at 72%% parallel efficiency and a throughput peak "
               "at nprocs=1. A free consequence is that %%pal nprocs 1 invokes no MPI at all, so "
               "the MPI_ERR_ARG ceiling it measured between 936 and 1,018 basis functions cannot "
               "arise -- exactly one of the 5,655 compounds sits above it.",
        alternatives="One 6-core job at a time (2.54x slower by notebook 42's own measurement, and "
                     "exposed to the MPI defect on the largest compounds).",
        authority="the brief, which specifies the pool design")
    record_decision(
        decision="Estimated basis functions exactly from element composition for the largest-first "
                 "queue, rather than regressing on atom count",
        reason="Basis functions are exactly additive over atoms. The per-element table was solved "
               "by least squares against notebook 42's 15 measured counts (max residual 0.000000) "
               "and completed for Br and P by two dedicated single-point jobs, since neither "
               "element was determined by that selection. It reproduces all 15 exactly and covers "
               "every element in the population, so the queue order is exact rather than fitted.",
        alternatives="Use notebook 42's basis-vs-atoms regression (R2 0.79, so it would mis-order "
                     "the queue); or order by heavy atoms (ignores that iodine contributes 26 "
                     "basis functions against carbon's 14).",
        authority="none -- the brief asks for a basis-function estimate without specifying one")

    pop = load_population(args.limit)
    done = already_done(pop)
    todo = pop[~pop.Molecule_Name.isin(done)]
    log.info("population %d | already complete %d | to run %d", len(pop), len(done), len(todo))
    log.info("queue is largest-first: basis %d down to %d",
             int(pop.n_basis_est.max()), int(pop.n_basis_est.min()))

    failed: list[dict] = []
    if FAILURES_CSV.exists():
        try:
            failed = pd.read_csv(FAILURES_CSV).to_dict("records")
        except (pd.errors.EmptyDataError, OSError):
            failed = []
    durations: list[float] = []
    write_status(len(pop), done, failed, durations, started_at, args)

    def _signal(signum, _frame):
        if not _STOP.is_set():
            log.warning("signal %d received -- no new jobs will start; in-flight jobs will "
                        "finish and state will be written", signum)
            _STOP.set()

    signal.signal(signal.SIGINT, _signal)
    signal.signal(signal.SIGTERM, _signal)

    deadline = started_at + args.hours * 3600 if args.hours else None
    timeout_a, timeout_b = args.timeout_a_min * 60, args.timeout_b_min * 60
    rows = list(todo.itertuples(index=False))
    submitted = 0
    stop_reason = "queue exhausted"

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {}
        it = iter(rows)

        def _submit_next() -> bool:
            nonlocal submitted
            if _STOP.is_set():
                return False
            if deadline and time.time() >= deadline:
                return False
            try:
                r = next(it)
            except StopIteration:
                return False
            futures[pool.submit(run_one_compound, r, timeout_a, timeout_b)] = r.Molecule_Name
            submitted += 1
            return True

        for _ in range(args.jobs):
            _submit_next()

        while futures:
            for fut in as_completed(list(futures)):
                name = futures.pop(fut)
                try:
                    rec = fut.result()
                except Exception as exc:  # noqa: BLE001
                    rec = {"Molecule_Name": name, "status": "FAILED",
                           "error": f"worker crashed: {type(exc).__name__}: {exc}",
                           "total_s": 0.0}
                with _LOCK:
                    hirsh = rec.pop("hirshfeld", None)
                    if rec["status"] == "OK":
                        append_rows(DESCRIPTORS_CSV, [rec])
                        if hirsh:
                            append_rows(HIRSHFELD_CSV, hirsh, key="Molecule_Name")
                        done.add(name)
                        if not rec.get("cached"):
                            durations.append(rec["total_s"])
                        log.info("%-14s OK   %6.1fs  basis=%s  gap=%.3f eV  (%d/%d)",
                                 name, rec["total_s"], rec.get("n_basis"),
                                 rec.get("dft_gap_ev", float("nan")), len(done), len(pop))
                    else:
                        failed = [f for f in failed if f["Molecule_Name"] != name]
                        failed.append({"Molecule_Name": name, "error": rec["error"],
                                       "at_utc": datetime.now(timezone.utc).isoformat()})
                        append_rows(FAILURES_CSV, failed)
                        log.error("%-14s FAILED  %s", name, rec["error"])
                    status = write_status(len(pop), done, failed, durations, started_at, args)
                    if durations and len(durations) % 25 == 0:
                        try:
                            pd.read_csv(DESCRIPTORS_CSV).to_parquet(DESCRIPTORS_PARQUET,
                                                                    index=False)
                        except Exception as exc:  # noqa: BLE001 -- parquet is a convenience
                            log.warning("parquet write skipped: %s", exc)
                if not _submit_next():
                    if _STOP.is_set():
                        stop_reason = "signal"
                    elif deadline and time.time() >= deadline:
                        stop_reason = "hours budget spent"
                break  # re-enter as_completed over the updated future set

    try:
        if DESCRIPTORS_CSV.exists():
            pd.read_csv(DESCRIPTORS_CSV).to_parquet(DESCRIPTORS_PARQUET, index=False)
    except Exception as exc:  # noqa: BLE001
        log.warning("final parquet write skipped: %s", exc)
    status = write_status(len(pop), done, failed, durations, started_at, args)
    log.info("stopped: %s | submitted %d | done %d/%d | failed %d | %.2f h elapsed",
             stop_reason, submitted, status["done"], status["total_compounds"],
             status["failed"], status["session_elapsed_h"])
    if status["projected_remaining_h"]:
        log.info("projected %.1f h of compute remaining at --jobs %d",
                 status["projected_remaining_h"], args.jobs)
    return 0


if __name__ == "__main__":
    sys.exit(main())
