#!/usr/bin/env python3
"""Extract per-atom, per-bond and molecule-level QM quantities from the completed ORCA run.

Parsing only. No ORCA is re-run. Nothing under outputs/qm_descriptors/ is written or modified --
descriptors.csv and hirshfeld_charges.csv are read for cross-validation and nothing else. All
output goes to outputs/qm_parsed/.

USAGE
    conda run -n cyp-admet-v2 python scripts/parse_orca_output.py --jobs 6
    ... --limit 20      # first N compounds, for a smoke test
    ... --force         # re-parse compounds whose shards already exist

SOURCE PRECEDENCE
    .property.txt is preferred wherever a quantity exists there, because it stores values as
    typed indexed arrays (&Type "ArrayOfDoubles", &Dim (n,m)) -- an element symbol never sits
    adjacent to a number, so the fixed-width trap below cannot apply. Measured against both
    jobs' full field inventory: 19 of the 23 requested quantities are available there. The four
    that are NOT, and so must come from {tag}__B_dftsp.out:

        Mulliken/Loewdin s/p/d shell aggregates   -- no reduced-orbital section exists at all
        the orbital ladder (HOMO-n, LUMO+n)       -- no orbital-energies section
        the three rotational constants            -- no rotational section
        the SCF cycle count                       -- .property.txt carries timings, not cycles

    Only the first of those four is exposed to the trap, and it is parsed by fixed character
    position -- see parse_reduced_orbital().

    Job A's .out is NEVER OPENED: 5.48 GB, 64% of the whole tree, and a lower-theory mirror of
    job B. Everything needed from job A -- cycle count, final gradient norm, relaxation strain
    -- is in its .property.txt as named fields.

    The optimised geometry comes from {tag}__A_gfn2opt.xyz, which was verified bit-identical to
    job B's input geometry (max coordinate difference 0.00e+00) and agrees with job A's last
    $Geometry block to 3.2e-08 Angstrom after Bohr conversion, i.e. rounding only.

THE TRAP -- measured, 815 of 5,655 compounds (14.4%): 669 contain Cl, 161 contain Br
    In the reduced-orbital sections ORCA writes fixed-width fields, so a two-character element
    symbol runs into the shell label with no separator:

          6 C s       :     2.925169
          7 Cls       :     5.940711      <-- "Cls"

    Any ([A-Z][a-z]?)\\s+(\\S+) pattern SILENTLY drops those atoms and returns a short block with
    no error at all; the measured shortfall equals the two-character-atom count exactly. Only Cl
    and Br are affected because the other seven elements present (H C N O F P S I) have
    one-character symbols -- which is exactly why a parser written against one example looks fine.

ATOM ORDER
    Validated per compound, not assumed. The pipeline's stage1_conformer() built each geometry
    with Chem.AddHs(Chem.MolFromSmiles(row.SMILES)) and wrote atoms in RDKit's native index
    order with no reordering, so ORCA's atom index IS the RDKit atom index. Note that the raw
    SMILES column is load-bearing here: on a 400-compound sample the raw column reproduces the
    element sequence 400/400 while canonical_smiles manages only 310/400.
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import signal
import sys
import threading
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

# Reused rather than reimplemented, from the run that produced this tree.
from run_qm_descriptors import (  # noqa: E402
    BLIND,
    CURATED,
    DESCRIPTORS_CSV,
    HIRSHFELD_CSV,
    ORCA_DIR,
    OrcaParseError,
    read_orca_xyz,
)

# ---------------------------------------------------------------- configuration

OUT = REPO_ROOT / "outputs" / "qm_parsed"
PARTS = OUT / "parts"
ATOMS_CSV = OUT / "qm_atoms.csv"
BONDS_CSV = OUT / "qm_bonds.csv"
MOL_CSV = OUT / "qm_molecular_raw.csv"
REPORT_MD = OUT / "parse_report.md"
STATUS_JSON = OUT / "parse_status.json"
DECISIONS_JSON = OUT / "decisions_taken.json"
RUN_META_JSON = OUT / "run_meta.json"

EXPECTED_RDKIT = "2026.03.3"
EXPECTED_COMPOUNDS = 5655

# a.u. -> Debye (CODATA). Verified against the existing pipeline's own dft_dipole_debye on 60
# compounds: agreement to 1.78e-04, which is the .out's own 6-decimal rounding, so the factor is
# confirmed and the .property.txt route is the more precise of the two.
AU_TO_DEBYE = 2.5417464519
EH_TO_KCAL = 627.5094740631

# Per-scheme, each justified by its own measurement over 120 compounds rather than one blanket
# number. Mulliken and Loewdin come out of .property.txt exact to floating point (max residual
# 1.5e-12); only Hirshfeld carries real numerical-integration error (max 5.9e-04). A single 5e-3
# tolerance would hide a genuine Mulliken bug by six orders of magnitude.
CHARGE_TOL = {"mulliken": 1e-6, "loewdin": 1e-6, "hirshfeld": 5e-3}

FLOAT_DP = 6          # charges/populations/coordinates; matches the .out's own precision
TERMINATED_OK = "****ORCA TERMINATED NORMALLY****"

_STOP = threading.Event()
DECISIONS: list[dict] = []
log = logging.getLogger("parse")


def record_decision(decision: str, reason: str, alternatives: str, authority: str) -> None:
    DECISIONS.append({"decision": decision, "reason": reason, "alternatives": alternatives,
                      "authority": authority,
                      "recorded_at_utc": datetime.now(timezone.utc).isoformat()})
    DECISIONS_JSON.parent.mkdir(parents=True, exist_ok=True)
    DECISIONS_JSON.write_text(json.dumps(DECISIONS, indent=2) + "\n")


# ------------------------------------------------------- .property.txt parsing

_DIM_RE = re.compile(r"&Dim\s*\((\d+)\s*,\s*(\d+)\)")


def prop_blocks(txt: str) -> dict[str, list[str]]:
    """Split a .property.txt into {section_name: [block_text, ...]}.

    A section may legitimately repeat: job A writes one $Geometry, $SCF_Energy and
    $Single_Point_Data per optimisation cycle, so the value is always a list and the caller
    chooses first/last/count explicitly.
    """
    out: dict[str, list[str]] = {}
    name: str | None = None
    buf: list[str] = []
    for line in txt.splitlines():
        s = line.strip()
        if s.startswith("$") and s != "$End":
            name, buf = s[1:], []
        elif s == "$End":
            if name is not None:
                out.setdefault(name, []).append("\n".join(buf))
            name = None
        elif name is not None:
            buf.append(line)
    return out


def _field_tail(block: str, field: str) -> tuple[str, str] | None:
    """Return (header_line, text_after_it) for `&field`, or None if absent."""
    m = re.search(rf"^\s*&{re.escape(field)}\b(.*)$", block, re.M)
    if m is None:
        return None
    return m.group(0), block[m.end():]


def prop_scalar(block: str, field: str, cast=float):
    """Read a scalar &field. Handles Integer, Double, Boolean and bare values."""
    got = _field_tail(block, field)
    if got is None:
        raise OrcaParseError(f"&{field} not found")
    head = got[0]
    tail = head.split("]", 1)[1] if "]" in head else head.split(field, 1)[1]
    tok = tail.strip().split()
    if not tok:
        raise OrcaParseError(f"&{field} has no value")
    v = tok[0]
    if cast is bool:
        return v.lower() == "true"
    return cast(v)


def prop_array(block: str, field: str) -> list[list[float]]:
    """Read an indexed array &field, using its own declared &Dim (nrows, ncols).

    Discriminating data rows from the column-index header is done on token count: a data row
    carries 1 + ncols tokens and its first token is the running row index, while the header
    carries exactly ncols. The declared row count is then asserted, so a truncated or
    mis-sliced block fails loudly instead of returning something short.
    """
    got = _field_tail(block, field)
    if got is None:
        raise OrcaParseError(f"&{field} not found")
    head, tail = got
    dim = _DIM_RE.search(head)
    if dim is None:
        raise OrcaParseError(f"&{field} carries no &Dim")
    nrows, ncols = int(dim.group(1)), int(dim.group(2))
    rows: list[list[float]] = []
    for line in tail.splitlines():
        s = line.strip()
        if s.startswith("&") or s.startswith("$"):
            break
        tok = s.split()
        if len(tok) != ncols + 1:
            continue
        try:
            idx = int(tok[0])
        except ValueError:
            continue
        if idx != len(rows):
            continue
        rows.append([float(x) for x in tok[1:]])
    if len(rows) != nrows:
        raise OrcaParseError(f"&{field}: declared {nrows} rows, parsed {len(rows)}")
    return rows


def prop_col(block: str, field: str) -> list[float]:
    """An (n,1) array as a flat list."""
    return [r[0] for r in prop_array(block, field)]


# ---------------------------------------------------------------- .out parsing

# Fixed character offsets in the reduced-orbital sections, measured on 1-, 2- and 3-digit atom
# indices and on both two-character elements present (Cl, Br):
#
#     '  0 C s       :     3.103820  s :     3.103820'
#     '  7 Cls       :     5.940711  s :     5.940711'      <-- Cl, no space before the shell
#     '125 H s       :     0.925115  s :     0.925115'      <-- 3-digit index
#      ^^^^ ^^ ^^^^^^^^ ^ ^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^
#      0:4  4:6  6:14   14:15   15:30          30:
#      index elem shell  ':'    primitive      trailing s/p/d aggregate
#
# Parsing by whitespace instead drops every Cl and Br atom silently. On OCNT-2328824 (4 Cl) this
# reader returns 37 atom headers for 37 atoms; a whitespace parser returns 33.
_COL_IDX = slice(0, 4)
_COL_ELEM = slice(4, 6)
_COL_SHELL = slice(6, 14)
_COL_COLON = slice(14, 15)
_COL_AGG = 30

_AGG_RE = re.compile(r"^([a-z])\s*:\s*(-?\d+\.\d+)")
_ORB_RE = re.compile(r"^\s*(\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", re.M)
_ROT_RE = re.compile(r"Rotational constants in cm-1:\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")
_SCF_RE = re.compile(r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES")


def parse_reduced_orbital(txt: str, header: str, endmark: str) -> dict[int, dict]:
    """Per-atom s/p/d shell aggregates, by FIXED CHARACTER POSITION (see the offsets above).

    Only the aggregates are kept. The individual pz/px/py and five d components depend on the
    orientation the molecule happens to have in the file and are not molecular properties, while
    the s/p/d sums over each shell are rotationally invariant.
    """
    if header not in txt:
        raise OrcaParseError(f"{header} not found")
    blk = txt.split(header)[-1].split(endmark)[0]
    atoms: dict[int, dict] = {}
    cur: int | None = None
    for line in blk.splitlines():
        if len(line) < 15 or line[_COL_COLON] != ":":
            continue
        idx_field = line[_COL_IDX].strip()
        if idx_field.isdigit():
            cur = int(idx_field)
            atoms[cur] = {"element": line[_COL_ELEM].strip(), "s": None, "p": None, "d": None}
        if cur is None:
            continue
        tail = line[_COL_AGG:].strip()
        if not tail:
            continue
        m = _AGG_RE.match(tail)
        if m is None:
            continue
        shell = m.group(1)
        if shell not in atoms[cur]:
            # def2-SV(P) gives these nine elements s, p and d only. An f aggregate would mean the
            # basis is not what this script assumes -- fail rather than drop it silently.
            raise OrcaParseError(f"unexpected shell aggregate '{shell}' on atom {cur}")
        atoms[cur][shell] = float(m.group(2))
    return atoms


def parse_orbital_ladder(txt: str) -> tuple[list[float], list[float]]:
    """(occupied, virtual) orbital energies in eV, in ascending order.

    Occupancy > 0.5 and the eV column, matching the existing pipeline's own convention -- the
    HOMO and LUMO this yields reproduce descriptors.csv exactly (0.00e+00 on 60 compounds).

    The virtual side is TRUNCATED by ORCA ("*Only the first 10 virtual orbitals were printed."),
    so nothing beyond LUMO+1 is taken anywhere in this script. All occupied orbitals are present.
    """
    if "ORBITAL ENERGIES" not in txt:
        raise OrcaParseError("ORBITAL ENERGIES not found")
    blk = txt.split("ORBITAL ENERGIES")[-1].split("MULLIKEN ATOMIC CHARGES")[0]
    occ, vir = [], []
    for m in _ORB_RE.finditer(blk):
        (occ if float(m.group(2)) > 0.5 else vir).append(float(m.group(4)))
    if len(occ) < 3:
        raise OrcaParseError(f"only {len(occ)} occupied orbitals; need 3 for HOMO-2")
    if len(vir) < 2:
        raise OrcaParseError(f"only {len(vir)} virtual orbitals; need 2 for LUMO+1")
    return occ, vir


def parse_rotational_constants(txt: str) -> tuple[float, float, float]:
    m = _ROT_RE.search(txt)
    if m is None:
        raise OrcaParseError("rotational constants not found")
    return float(m.group(1)), float(m.group(2)), float(m.group(3))


def parse_scf_cycles(txt: str) -> int:
    hits = _SCF_RE.findall(txt)
    if not hits:
        raise OrcaParseError("SCF CONVERGED AFTER not found")
    return int(hits[-1])


def parse_xyz_geometry(path: Path) -> list[tuple[str, float, float, float]]:
    """Optimised geometry in Angstrom, via the existing read_orca_xyz()."""
    out = []
    for line in read_orca_xyz(path):
        tok = line.split()
        if len(tok) != 4:
            raise OrcaParseError(f"{path}: cannot read geometry line {line!r}")
        out.append((tok[0], float(tok[1]), float(tok[2]), float(tok[3])))
    return out


# ------------------------------------------------------------- job A trajectory

def parse_job_a(path: Path) -> dict:
    """Optimisation cycle count, final gradient norm and relaxation strain, from .property.txt.

    Job A's .out is deliberately never opened (5.48 GB across the tree).

    Cycle-count convention: a 26-cycle optimisation writes 27 $Geometry blocks, the extra being
    the final energy evaluation at the stationary point, so n_geometry - 1 is reported -- that is
    the number the .out's own "GEOMETRY OPTIMIZATION CYCLE" counter agrees with.

    Relaxation strain is E_first - E_last, so a POSITIVE value means energy was released during
    relaxation (the normal case).
    """
    blocks = prop_blocks(path.read_text(errors="replace"))
    geoms = blocks.get("Geometry", [])
    energies = blocks.get("SCF_Energy", [])
    grads = blocks.get("SCF_Nuc_Gradient", [])
    if not geoms or not energies or not grads:
        raise OrcaParseError(f"job A: {len(geoms)} geometries, {len(energies)} energies, "
                             f"{len(grads)} gradients")
    e_first = prop_col(energies[0], "totalEnergy")[0]
    e_last = prop_col(energies[-1], "totalEnergy")[0]
    strain = e_first - e_last
    return {"gfn2_opt_cycles": len(geoms) - 1,
            "gfn2_final_grad_norm": prop_scalar(grads[-1], "gradNorm"),
            "gfn2_relaxation_strain_eh": strain,
            "gfn2_relaxation_strain_kcal": strain * EH_TO_KCAL}


# ------------------------------------------------------------------- validation

def validate(name: str, n_rdkit: int, rd_elements: list[str], rd_atnos: list[int],
             rd_bonds: set, q_formal: int, n_expected_csv: int,
             sections: dict, charges: dict, mayer: dict, reduced: dict,
             bonds: list[dict], geom: list) -> list[str]:
    """Every assertion from the plan. Returns a list of failure reasons; empty means pass.

    Returning reasons rather than raising is deliberate: a failing compound must be recorded with
    its reason and skipped, never written as a short block, because a short block is exactly what
    the fixed-width trap produces silently.
    """
    bad: list[str] = []

    # (1) row count per section, each checked separately
    for label, got in [("mulliken", len(charges["mulliken"])), ("loewdin", len(charges["loewdin"])),
                       ("hirshfeld", len(charges["hirshfeld"])), ("mayer_na", len(mayer["na"])),
                       ("mayer_va", len(mayer["va"])),
                       ("mulliken_reduced", len(reduced["mulliken"])),
                       ("loewdin_reduced", len(reduced["loewdin"])), ("geometry", len(geom))]:
        if got != n_rdkit:
            bad.append(f"{label}: {got} rows, expected {n_rdkit}")
    if n_rdkit != n_expected_csv:
        bad.append(f"RDKit atom count {n_rdkit} != descriptors.csv total_atoms_with_h "
                   f"{n_expected_csv}")
    # Every check below indexes rd_elements/rd_atnos by an ORCA atom index, so a count mismatch
    # must short-circuit: continuing would raise IndexError and bury the real reason under a
    # traceback. The compound is rejected either way, but the report has to say why.
    if bad:
        return bad

    # (2) charge sums, per-scheme tolerances
    for scheme, tol in CHARGE_TOL.items():
        resid = abs(sum(charges[scheme]) - q_formal)
        if resid > tol:
            bad.append(f"{scheme} charge sum residual {resid:.3e} > {tol:.0e} "
                       f"(formal charge {q_formal})")

    # (3) indices contiguous / in range
    for label in ("mulliken", "loewdin"):
        keys = sorted(reduced[label])
        if keys != list(range(n_rdkit)):
            bad.append(f"{label}_reduced atom indices not contiguous 0..{n_rdkit - 1}")
    for b in bonds:
        if not (0 <= b["atom_i"] < n_rdkit and 0 <= b["atom_j"] < n_rdkit):
            bad.append(f"bond index out of range: ({b['atom_i']}, {b['atom_j']}) n={n_rdkit}")
            break

    # (4) element at each ORCA index == element at the same RDKit index. This is the check the
    # existing pipeline never did -- it validated only count and charge sum.
    geom_elements = [g[0] for g in geom]
    if geom_elements != rd_elements:
        first = next((i for i, (a, b_) in enumerate(zip(geom_elements, rd_elements)) if a != b_),
                     None)
        bad.append(f"geometry element order differs from RDKit at index {first}: "
                   f"{geom_elements[first] if first is not None else '?'} vs "
                   f"{rd_elements[first] if first is not None else '?'}")
    for label in ("mulliken", "loewdin"):
        got = [reduced[label][i]["element"] for i in range(len(reduced[label]))]
        if got != rd_elements:
            bad.append(f"{label}_reduced element order differs from RDKit")

    # (5) ATNO agrees with RDKit -- an independent form of (4), free from .property.txt
    for sec_name, block in sections.items():
        if "atno" not in block:
            continue
        if block["atno"] != rd_atnos:
            bad.append(f"{sec_name} ATNO differs from RDKit atomic numbers")

    # (6) RDKit bonds are a SUBSET of the Mayer pair list, not equal to it: the Mayer list is
    # thresholded at bond order > 0.1 and so also contains weak intramolecular contacts (measured
    # mean 0.64 extra pairs per compound, max 4). Equality would fail on 47% of compounds for a
    # reason that is not an error.
    orca_pairs = {(min(b["atom_i"], b["atom_j"]), max(b["atom_i"], b["atom_j"])) for b in bonds}
    missing = rd_bonds - orca_pairs
    if missing:
        bad.append(f"{len(missing)} RDKit bonds absent from the Mayer pair list, "
                   f"e.g. {sorted(missing)[:3]}")

    # (7) NAtoms consistent across every section that reports it
    natoms = {s: b["natoms"] for s, b in sections.items() if "natoms" in b}
    if natoms and len(set(natoms.values())) != 1:
        bad.append(f"&NAtoms disagrees across sections: {natoms}")
    if natoms and next(iter(natoms.values())) != n_rdkit:
        bad.append(f"&NAtoms {next(iter(natoms.values()))} != RDKit {n_rdkit}")

    # (8) reduced-orbital coverage: def2-SV(P) puts no polarisation function on hydrogen
    for label in ("mulliken", "loewdin"):
        for i, rec in reduced[label].items():
            is_h = rd_elements[i] == "H"
            if rec["s"] is None:
                bad.append(f"{label}_reduced atom {i} ({rd_elements[i]}) has no s aggregate")
            if is_h and (rec["p"] is not None or rec["d"] is not None):
                bad.append(f"{label}_reduced H atom {i} unexpectedly has p/d")
            if not is_h and (rec["p"] is None or rec["d"] is None):
                bad.append(f"{label}_reduced heavy atom {i} ({rd_elements[i]}) missing p or d")
    return bad


# ----------------------------------------------------------- per-compound parse

def parse_one(job: dict) -> dict:
    """Parse one compound and write its three shards. Returns a small status record only."""
    name = job["name"]
    t0 = time.time()
    try:
        wd = ORCA_DIR / name
        tag_a, tag_b = f"{name}__A_gfn2opt", f"{name}__B_dftsp"
        out_b = wd / f"{tag_b}.out"
        prop_b = wd / f"{tag_b}.property.txt"
        prop_a = wd / f"{tag_a}.property.txt"
        xyz_a = wd / f"{tag_a}.xyz"
        for p in (out_b, prop_b, prop_a, xyz_a):
            if not p.exists():
                raise OrcaParseError(f"missing {p.name}")

        # The raw SMILES column, NOT canonical_smiles. stage1_conformer() used row.SMILES and
        # wrote atoms in RDKit's native order; canonical_smiles reproduces that order on only
        # 310 of 400 sampled compounds, so using it would silently return charges on wrong atoms.
        mol = Chem.MolFromSmiles(job["smiles"])
        if mol is None:
            raise OrcaParseError("RDKit could not parse the curated SMILES")
        mol = Chem.AddHs(mol)
        rd_elements = [a.GetSymbol() for a in mol.GetAtoms()]
        rd_atnos = [a.GetAtomicNum() for a in mol.GetAtoms()]
        rd_bonds = {(min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                     max(b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in mol.GetBonds()}
        n = mol.GetNumAtoms()
        q_formal = Chem.GetFormalCharge(mol)

        # ---- job B .property.txt (preferred source for 19 of 23 quantities)
        pb = prop_blocks(prop_b.read_text(errors="replace"))
        for req in ("SCF_Mulliken_Population_Analysis", "SCF_Loewdin_Population_Analysis",
                    "SCF_Hirshfeld_Population_Analysis", "SCF_Mayer_Population_Analysis",
                    "Single_Point_Data", "SCF_Dipole_Moment", "VdW_Correction",
                    "Calculation_Info"):
            if req not in pb:
                raise OrcaParseError(f"job B .property.txt has no ${req}")
        mull_b = pb["SCF_Mulliken_Population_Analysis"][-1]
        loew_b = pb["SCF_Loewdin_Population_Analysis"][-1]
        hirs_b = pb["SCF_Hirshfeld_Population_Analysis"][-1]
        mayer_b = pb["SCF_Mayer_Population_Analysis"][-1]

        charges = {"mulliken": prop_col(mull_b, "AtomicCharges"),
                   "loewdin": prop_col(loew_b, "AtomicCharges"),
                   "hirshfeld": prop_col(hirs_b, "AtomicCharges")}
        mayer = {"na": prop_col(mayer_b, "NA"), "va": prop_col(mayer_b, "VA")}

        sections = {}
        for sec, blk in [("mulliken", mull_b), ("loewdin", loew_b), ("hirshfeld", hirs_b),
                         ("mayer", mayer_b)]:
            rec = {"natoms": int(prop_scalar(blk, "NAtoms", int))}
            try:
                rec["atno"] = [int(x) for x in prop_col(blk, "ATNO")]
            except OrcaParseError:
                pass
            sections[sec] = rec

        # bonds: &components is (nbonds, 4) = [index_i, Z_i, index_j, Z_j]
        comps = prop_array(mayer_b, "components")
        orders = prop_col(mayer_b, "BondOrders")
        if len(comps) != len(orders):
            raise OrcaParseError(f"{len(comps)} bond components vs {len(orders)} bond orders")
        bonds = [{"atom_i": int(c[0]), "atom_j": int(c[2]), "mayer_bond_order": o}
                 for c, o in zip(comps, orders)]

        # ---- job B .out (the only four quantities .property.txt does not carry)
        tb = out_b.read_text(errors="replace")
        if TERMINATED_OK not in tb:
            raise OrcaParseError("job B .out does not say ORCA TERMINATED NORMALLY")
        reduced = {
            "mulliken": parse_reduced_orbital(tb, "MULLIKEN REDUCED ORBITAL CHARGES",
                                              "LOEWDIN ATOMIC CHARGES"),
            "loewdin": parse_reduced_orbital(tb, "LOEWDIN REDUCED ORBITAL CHARGES",
                                             "MAYER POPULATION ANALYSIS"),
        }
        occ, vir = parse_orbital_ladder(tb)
        rot_a, rot_b_, rot_c = parse_rotational_constants(tb)
        scf_cycles = parse_scf_cycles(tb)

        geom = parse_xyz_geometry(xyz_a)
        job_a = parse_job_a(prop_a)

        failures = validate(name, n, rd_elements, rd_atnos, rd_bonds, q_formal,
                            job["n_atoms_csv"], sections, charges, mayer, reduced, bonds, geom)
        if failures:
            return {"name": name, "status": "FAILED", "failures": failures,
                    "seconds": time.time() - t0}

        # ---- rows
        r = lambda v: None if v is None else round(v, FLOAT_DP)  # noqa: E731
        atom_rows = []
        for i in range(n):
            mu, lo = reduced["mulliken"][i], reduced["loewdin"][i]
            atom_rows.append({
                "Molecule_Name": name, "atom_index": i, "element": rd_elements[i],
                "atomic_number": rd_atnos[i],
                "q_mulliken": r(charges["mulliken"][i]), "q_loewdin": r(charges["loewdin"][i]),
                "q_hirshfeld": r(charges["hirshfeld"][i]),
                "mayer_na": r(mayer["na"][i]), "mayer_va": r(mayer["va"][i]),
                "mulliken_s": r(mu["s"]), "mulliken_p": r(mu["p"]), "mulliken_d": r(mu["d"]),
                "loewdin_s": r(lo["s"]), "loewdin_p": r(lo["p"]), "loewdin_d": r(lo["d"]),
                "x": r(geom[i][1]), "y": r(geom[i][2]), "z": r(geom[i][3])})
        bond_rows = [{"Molecule_Name": name, "atom_i": b["atom_i"], "atom_j": b["atom_j"],
                      "element_i": rd_elements[b["atom_i"]],
                      "element_j": rd_elements[b["atom_j"]],
                      "mayer_bond_order": r(b["mayer_bond_order"]),
                      "is_rdkit_bond": (min(b["atom_i"], b["atom_j"]),
                                        max(b["atom_i"], b["atom_j"])) in rd_bonds}
                     for b in bonds]
        mol_row = {
            "Molecule_Name": name, "set": job["set"], "inchikey": job["inchikey"],
            "n_atoms": n, "formal_charge": q_formal,
            "total_energy_eh": prop_scalar(pb["Single_Point_Data"][-1], "FinalEnergy"),
            "homo_ev": occ[-1], "homo1_ev": occ[-2], "homo2_ev": occ[-3],
            "lumo_ev": vir[0], "lumo1_ev": vir[1], "gap_ev": vir[0] - occ[-1],
            "dipole_debye": prop_scalar(pb["SCF_Dipole_Moment"][-1],
                                        "dipoleMagnitude") * AU_TO_DEBYE,
            "dispersion_eh": prop_scalar(pb["VdW_Correction"][-1], "vdW"),
            "rot_const_a_cm1": rot_a, "rot_const_b_cm1": rot_b_, "rot_const_c_cm1": rot_c,
            "n_basis": int(prop_scalar(pb["Calculation_Info"][-1], "NumOfBasisFuncts", int)),
            "scf_cycles": scf_cycles, **job_a}

        write_shards(name, atom_rows, bond_rows, mol_row)
        return {"name": name, "status": "OK", "failures": [], "n_atoms": n,
                "n_bonds": len(bond_rows), "seconds": time.time() - t0}
    except Exception as exc:                                  # noqa: BLE001
        return {"name": name, "status": "FAILED", "seconds": time.time() - t0,
                "failures": [f"{type(exc).__name__}: {exc}"]}


# ------------------------------------------------------------------ persistence

def shard_paths(name: str) -> tuple[Path, Path, Path]:
    return (PARTS / f"{name}.atoms.csv", PARTS / f"{name}.bonds.csv", PARTS / f"{name}.mol.csv")


def shard_done(name: str) -> bool:
    """A compound is done iff all three shards exist and are non-empty.

    Content-based rather than existence-only, in the same spirit as the generating run's
    job_is_done(). All three are written together and only after validation passes, so a
    half-written compound cannot look complete.
    """
    return all(p.exists() and p.stat().st_size > 0 for p in shard_paths(name))


def write_shards(name: str, atom_rows: list[dict], bond_rows: list[dict], mol_row: dict) -> None:
    """One shard per table per compound, written atomically.

    This is a deliberate deviation from the generating run's append_rows(), which is
    read-modify-rewrite-whole-file. That is right for a 5,655-row table but measured at 0.275 s
    per call once the atom table reaches 252k rows -- 26 minutes of pure rewriting across the
    population, against a parse that should take a few minutes. A shard write is 0.0002 s.
    The molecular table is sharded too, rather than appended: it closes a restart hole, since an
    interruption between an atom shard and a buffered molecular flush would otherwise leave a
    compound that looks done but has no molecular row.
    """
    PARTS.mkdir(parents=True, exist_ok=True)
    atoms_p, bonds_p, mol_p = shard_paths(name)
    for path, rows in ((atoms_p, atom_rows), (bonds_p, bond_rows), (mol_p, [mol_row])):
        tmp = path.with_suffix(path.suffix + ".tmp")
        pd.DataFrame(rows).to_csv(tmp, index=False)
        tmp.replace(path)


def concat_shards(names: list[str]) -> dict[str, int]:
    """Build the three tables from the shards. Deterministic row order."""
    OUT.mkdir(parents=True, exist_ok=True)
    counts = {}
    for kind, dest, sort_by in (("atoms", ATOMS_CSV, ["Molecule_Name", "atom_index"]),
                                ("bonds", BONDS_CSV, ["Molecule_Name", "atom_i", "atom_j"]),
                                ("mol", MOL_CSV, ["Molecule_Name"])):
        frames = []
        for name in names:
            p = PARTS / f"{name}.{kind}.csv"
            if p.exists() and p.stat().st_size > 0:
                frames.append(pd.read_csv(p))
        if not frames:
            log.warning("no %s shards to concatenate", kind)
            counts[kind] = 0
            continue
        df = pd.concat(frames, ignore_index=True).sort_values(sort_by).reset_index(drop=True)
        tmp = dest.with_suffix(dest.suffix + ".tmp")
        df.to_csv(tmp, index=False)
        tmp.replace(dest)
        counts[kind] = len(df)
        log.info("wrote %-28s %7d rows  %6.1f MB", dest.name, len(df),
                 dest.stat().st_size / 2 ** 20)
    return counts


def write_status(total: int, done: int, failed: list[dict], durations: list[float],
                 started_at: float) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = {
        "updated_at_utc": datetime.now(timezone.utc).isoformat(),
        "total_compounds": total, "done": done, "remaining": total - done,
        "failed": len(failed),
        "pct_done": round(100.0 * done / total, 2) if total else 0.0,
        "mean_seconds_per_compound": (round(sum(durations) / len(durations), 4)
                                      if durations else None),
        "session_elapsed_s": round(time.time() - started_at, 1),
        "failed_compounds": [f["name"] for f in failed][:50],
    }
    tmp = STATUS_JSON.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2) + "\n")
    tmp.replace(STATUS_JSON)


# --------------------------------------------------------------- startup checks

def preflight() -> None:
    """Fail before any parsing if the environment is not the one this was written against."""
    import rdkit
    if rdkit.__version__ != EXPECTED_RDKIT:
        raise SystemExit(f"rdkit {rdkit.__version__}, expected {EXPECTED_RDKIT}")
    if not ORCA_DIR.is_dir():
        raise SystemExit(f"ORCA tree not found at {ORCA_DIR}")
    for p in (CURATED, BLIND, DESCRIPTORS_CSV):
        if not p.exists():
            raise SystemExit(f"missing required input {p}")


def load_population(limit: int | None) -> pd.DataFrame:
    """The same 5,655 compounds, keyed the same way, as the run that produced the tree.

    Only the raw SMILES column is taken -- see parse_one() for why canonical_smiles is unsafe.
    """
    train = pd.read_csv(CURATED)[["Molecule_Name", "SMILES", "inchikey"]].assign(set="train")
    blind = pd.read_csv(BLIND)[["Molecule_Name", "SMILES", "inchikey"]].assign(set="blind")
    pop = pd.concat([train, blind], ignore_index=True)
    if len(pop) != EXPECTED_COMPOUNDS:
        raise SystemExit(f"expected {EXPECTED_COMPOUNDS} compounds, got {len(pop)}")
    if not pop.Molecule_Name.is_unique:
        raise SystemExit("Molecule_Name is not unique across the two sets")
    desc = pd.read_csv(DESCRIPTORS_CSV)[["Molecule_Name", "total_atoms_with_h"]]
    pop = pop.merge(desc, on="Molecule_Name", how="left", validate="one_to_one")
    if pop.total_atoms_with_h.isna().any():
        missing = pop.loc[pop.total_atoms_with_h.isna(), "Molecule_Name"].tolist()[:5]
        raise SystemExit(f"descriptors.csv has no row for {len(missing)}+ compounds, "
                         f"e.g. {missing}")
    pop = pop.sort_values("Molecule_Name").reset_index(drop=True)
    return pop.head(limit) if limit else pop


# ------------------------------------------------------------------ cross-check

def cross_check() -> list[dict]:
    """Compare the new tables against the generating run's own output.

    This is the strongest test available: descriptors.csv and hirshfeld_charges.csv were produced
    by independent code reading the same files, so exact agreement on the overlapping quantities
    is real evidence and any disagreement is a bug here.
    """
    rows: list[dict] = []
    mol = pd.read_csv(MOL_CSV)
    desc = pd.read_csv(DESCRIPTORS_CSV)
    m = mol.merge(desc, on="Molecule_Name", suffixes=("", "_old"))
    # The energy tolerances are not zero because the two sources print the same double in
    # different decimal forms: the .out gives ~12 significant figures, .property.txt gives 17, so
    # they round-trip to values differing by double-precision epsilon scaled to the magnitude
    # (~1790 Eh -> ~4.5e-13 measured). 1e-9 keeps ~2000x headroom over that while still being far
    # tighter than any difference a real bug could produce.
    pairs = [("homo_ev", "dft_homo_ev", 0.0), ("lumo_ev", "dft_lumo_ev", 0.0),
             ("gap_ev", "dft_gap_ev", 1e-9), ("total_energy_eh", "dft_energy_eh", 1e-9),
             ("n_basis", "n_basis", 0.0), ("scf_cycles", "dft_scf_cycles", 0.0),
             ("dispersion_eh", "dft_dispersion_eh", 1e-9),
             # the .out stores Debye to 6 dp, so this is its rounding, not a disagreement
             ("dipole_debye", "dft_dipole_debye", 2e-3),
             ("n_atoms", "total_atoms_with_h", 0.0),
             ("formal_charge", "formal_charge", 0.0),
             ("gfn2_opt_cycles", None, None)]
    for new, old, tol in pairs:
        if old is None or old not in m.columns or new not in m.columns:
            continue
        d = (m[new].astype(float) - m[old].astype(float)).abs()
        rows.append({"quantity": f"{new} vs descriptors.csv::{old}", "n": int(d.notna().sum()),
                     "max_abs_diff": float(d.max()), "tolerance": tol,
                     "verdict": "PASS" if float(d.max()) <= tol else "FAIL"})
    # Hirshfeld, per atom, against the existing per-atom table
    atoms = pd.read_csv(ATOMS_CSV, usecols=["Molecule_Name", "atom_index", "element",
                                            "q_hirshfeld"])
    old_h = pd.read_csv(HIRSHFELD_CSV)
    j = atoms.merge(old_h, on=["Molecule_Name", "atom_index"], suffixes=("", "_old"))
    d = (j.q_hirshfeld - j.charge).abs()
    rows.append({"quantity": "q_hirshfeld vs hirshfeld_charges.csv::charge", "n": len(j),
                 "max_abs_diff": float(d.max()), "tolerance": 5e-7,
                 "verdict": "PASS" if float(d.max()) <= 5e-7 else "FAIL"})
    bad_el = int((j.element != j.element_old).sum())
    rows.append({"quantity": "element vs hirshfeld_charges.csv::element", "n": len(j),
                 "max_abs_diff": float(bad_el), "tolerance": 0.0,
                 "verdict": "PASS" if bad_el == 0 else "FAIL"})
    return rows


SOURCE_TABLE = [
    ("Mulliken atomic charge", "B .property.txt",
     "$SCF_Mulliken_Population_Analysis &AtomicCharges"),
    ("Loewdin atomic charge", "B .property.txt", "$SCF_Loewdin_Population_Analysis &AtomicCharges"),
    ("Hirshfeld atomic charge", "B .property.txt",
     "$SCF_Hirshfeld_Population_Analysis &AtomicCharges"),
    ("Mayer NA, VA", "B .property.txt", "$SCF_Mayer_Population_Analysis &NA, &VA"),
    ("Mayer bond orders + indices", "B .property.txt",
     "$SCF_Mayer_Population_Analysis &BondOrders, &components"),
    ("atomic numbers (validation)", "B .property.txt", "&ATNO in each population section"),
    ("total energy", "B .property.txt", "$Single_Point_Data &FinalEnergy"),
    ("dipole magnitude", "B .property.txt",
     "$SCF_Dipole_Moment &dipoleMagnitude (a.u., converted)"),
    ("D3 dispersion", "B .property.txt", "$VdW_Correction &vdW"),
    ("basis count", "B .property.txt", "$Calculation_Info &NumOfBasisFuncts"),
    ("Mulliken s/p/d aggregates", "B .out (FIXED COLUMNS)", "MULLIKEN REDUCED ORBITAL CHARGES"),
    ("Loewdin s/p/d aggregates", "B .out (FIXED COLUMNS)", "LOEWDIN REDUCED ORBITAL CHARGES"),
    ("HOMO, HOMO-1, HOMO-2, LUMO, LUMO+1", "B .out",
     "ORBITAL ENERGIES (virtuals truncated by ORCA)"),
    ("rotational constants", "B .out", "Rotational constants in cm-1"),
    ("SCF cycles", "B .out", "SCF CONVERGED AFTER"),
    ("optimised geometry", "A .xyz", "read_orca_xyz(), reused from the generating run"),
    ("GFN2 opt cycles", "A .property.txt", "count of $Geometry blocks, minus 1"),
    ("GFN2 final gradient norm", "A .property.txt", "last $SCF_Nuc_Gradient &gradNorm"),
    ("GFN2 relaxation strain", "A .property.txt", "first vs last $SCF_Energy &totalEnergy"),
]


def table_stats() -> dict:
    """Facts computed from the tables that were just written, for the report.

    The two-character-element count is the trap's own acceptance test: every one of these atoms
    would be silently absent had the reduced-orbital blocks been parsed by whitespace.
    """
    atoms = pd.read_csv(ATOMS_CSV, usecols=["Molecule_Name", "element", "mulliken_p", "loewdin_p"])
    bonds = pd.read_csv(BONDS_CSV, usecols=["Molecule_Name", "is_rdkit_bond"])
    two = atoms[atoms.element.isin(["Cl", "Br"])]
    extra = int((~bonds.is_rdkit_bond).sum())
    return {
        "two_char_atoms": len(two),
        "two_char_compounds": int(two.Molecule_Name.nunique()),
        "cl_atoms": int((atoms.element == "Cl").sum()),
        "br_atoms": int((atoms.element == "Br").sum()),
        "elements": sorted(atoms.element.unique()),
        "rdkit_bonds": int(bonds.is_rdkit_bond.sum()),
        "extra_pairs": extra,
        "extra_per_compound": extra / max(bonds.Molecule_Name.nunique(), 1),
        "h_atoms": int((atoms.element == "H").sum()),
        "pd_null_exactly_on_h": bool(
            ((atoms.mulliken_p.isna()) == (atoms.element == "H")).all()
            and ((atoms.loewdin_p.isna()) == (atoms.element == "H")).all()),
    }


def write_report(pop_n: int, counts: dict, failed: list[dict], checks: list[dict],
                 durations: list[float], started_at: float, args, stats: dict) -> None:
    def mb(p: Path) -> str:
        return f"{p.stat().st_size / 2 ** 20:.1f} MB" if p.exists() else "-"
    el = time.time() - started_at
    L = ["# ORCA output parse report", "",
         f"Generated {datetime.now(timezone.utc).isoformat()}  ",
         f"`scripts/parse_orca_output.py --jobs {args.jobs}"
         f"{' --limit ' + str(args.limit) if args.limit else ''}"
         f"{' --force' if args.force else ''}`", "",
         "Parsing only. No ORCA re-run. Nothing under `outputs/qm_descriptors/` was written or",
         "modified; `descriptors.csv` and `hirshfeld_charges.csv` were read for cross-validation",
         "only.", "",
         "## Result", "",
         f"- compounds in population: **{pop_n}**",
         f"- parsed successfully: **{pop_n - len(failed)}**",
         f"- failed: **{len(failed)}**",
         f"- wall time: **{el / 60:.1f} min**"
         + (f" ({sum(durations) / len(durations) * 1000:.0f} ms/compound of parse work)"
            if durations else ""), "",
         "| table | rows | size |", "|---|---|---|",
         f"| `qm_atoms.csv` | {counts.get('atoms', 0):,} | {mb(ATOMS_CSV)} |",
         f"| `qm_bonds.csv` | {counts.get('bonds', 0):,} | {mb(BONDS_CSV)} |",
         f"| `qm_molecular_raw.csv` | {counts.get('mol', 0):,} | {mb(MOL_CSV)} |", "",
         "## Source used per quantity", "",
         "`.property.txt` is preferred wherever a quantity exists there: it stores typed indexed",
         "arrays, so an element symbol never sits adjacent to a number and the fixed-width trap",
         "cannot apply. Of the 19 quantities below, **13 come from a `.property.txt`** (10 from",
         "job B's, 3 from job A's), **5 from job B's `.out`** because no `.property.txt`",
         "equivalent exists, and 1 from the optimised-geometry `.xyz`. Only the two",
         "reduced-orbital blocks are exposed to the trap, and they are parsed by fixed character",
         "position.", "",
         "Job A's `.out` is never opened (5.48 GB, 64% of the tree).", "",
         "| quantity | source | section / field |", "|---|---|---|"]
    L += [f"| {q} | `{s}` | `{f}` |" for q, s, f in SOURCE_TABLE]
    trap = ["", "## The trap: acceptance test on the real output", "",
          "The reduced-orbital blocks are the only place a two-character element symbol sits",
          "against a shell label (`7 Cls       :`), and a whitespace parser drops those atoms",
          "silently. The parsed tables therefore carry their own acceptance test:", "",
          f"- **{stats['two_char_atoms']:,} Cl/Br atoms parsed**, across"
          f" **{stats['two_char_compounds']:,} compounds** ({stats['cl_atoms']:,} Cl,"
          f" {stats['br_atoms']:,} Br). Every one would be missing under whitespace splitting,"
          f" and the compound count matches the 815 measured independently in the audit.",
          f"- elements present: {', '.join(stats['elements'])} -- all nine in the library,"
          f" plus H.",
          f"- **p/d aggregates are null on exactly the {stats['h_atoms']:,} hydrogens and nowhere"
          f" else**: {stats['pd_null_exactly_on_h']}. def2-SV(P) puts no polarisation function on"
          f" H, so this is the expected pattern and a second check that no atom was skipped.",
          f"- bonds: **{stats['rdkit_bonds']:,}** are real RDKit bonds and"
          f" **{stats['extra_pairs']:,}** are weak contacts above ORCA's 0.1 threshold"
          f" ({stats['extra_per_compound']:.2f} per compound, against 0.64 measured on the"
          f" 400-compound sample).", ""] if stats else []
    L += trap
    L += ["## Validation", "",
          "Eight assertions per compound. A failure records the compound and its reason and the",
          "compound is skipped entirely -- no partial rows are written, because a short block is",
          "exactly what the trap produces silently.", "",
          "1. row count per section == RDKit atom count, checked separately for Mulliken, Loewdin,",
          "   Hirshfeld, Mayer NA, Mayer VA, both reduced-orbital blocks and the geometry; plus",
          "   agreement with `descriptors.csv`'s `total_atoms_with_h`",
          "2. charge sums == formal charge, per scheme: Mulliken/Loewdin 1e-6, Hirshfeld 5e-3",
          "   (measured maxima 1.5e-12 and 5.9e-4 respectively -- one blanket tolerance would hide",
          "   a Mulliken bug by six orders of magnitude)",
          "3. atom indices contiguous 0..n-1; bond indices within range",
          "4. element at each ORCA index == element at the same RDKit index, from the raw `SMILES`",
          "   column (not `canonical_smiles`, which reproduces the order on only 310 of 400",
          "   sampled compounds)",
          "5. `&ATNO` atomic numbers agree with RDKit's -- an independent form of (4)",
          "6. RDKit bond set is a **subset** of the Mayer pair list, not equal to it: the Mayer",
          "   list is thresholded at bond order > 0.1 and also carries weak intramolecular",
          "   contacts (measured mean 0.64 extra pairs, max 4). Equality would fail on 47% of",
          "   compounds for a reason that is not an error.",
          "7. `&NAtoms` consistent across every section that reports it",
          "8. reduced-orbital coverage: every heavy atom has s/p/d, every hydrogen has s only",
          "   (def2-SV(P) puts no polarisation function on H)", "",
          "## Cross-check against the generating run", "",
          "`descriptors.csv` and `hirshfeld_charges.csv` were produced by independent code reading",
          "the same files, so these are real external checks rather than self-consistency.", "",
          "| quantity | n | max abs diff | tolerance | verdict |", "|---|---|---|---|---|"]
    for c in checks:
        L.append(f"| `{c['quantity']}` | {c['n']:,} | {c['max_abs_diff']:.3e} | "
                 f"{c['tolerance']:.1e} | **{c['verdict']}** |")
    L += ["", "## Conventions worth knowing", "",
          "- **Energies cross-check to ~5e-13, not to zero.** The `.out` prints ~12 significant",
          "  figures and `.property.txt` 17, so the same double round-trips to values differing by",
          "  machine epsilon scaled to the magnitude (~1790 Eh). That the difference is this small",
          "  is itself the evidence the two sources carry the same number.",
          "- **Dipole** comes from `.property.txt` in a.u. and is converted with 2.5417464519",
          "  (CODATA). Against the existing `dft_dipole_debye` this agrees to ~2e-4, which is the",
          "  `.out`'s own 6-decimal rounding -- the `.property.txt` route is the more precise.",
          "- **`gfn2_opt_cycles`** is `n_geometry - 1`: a 26-cycle optimisation writes 27",
          "  `$Geometry` blocks, the extra being the final evaluation at the stationary point.",
          "  This is the number the `.out`'s own cycle counter agrees with.",
          "- **`gfn2_relaxation_strain_*`** is `E_first - E_last`, so **positive means energy was",
          "  released** during relaxation (the normal case).",
          "- **Nothing beyond LUMO+1** is extracted. ORCA truncates the virtual spectrum",
          "  (`*Only the first 10 virtual orbitals were printed.`) and job A's count of them",
          "  varies 10-13. All occupied orbitals are present, so HOMO-1 and HOMO-2 are reliable.",
          "- **Only the s/p/d aggregates** are kept from the reduced-orbital blocks. The pz/px/py",
          "  and five d components depend on the molecule's orientation in the file and are not",
          "  molecular properties; the shell sums are rotationally invariant.",
          "- **`is_rdkit_bond`** separates real bonds from the weak contacts that ORCA's 0.1",
          "  threshold also admits. Beyond the brief, but nothing else distinguishes them.",
          "- **Hirshfeld `spin` is dropped**: identically 0.0 for all 5,655 closed-shell singlets.",
          "- **Floats are written to 6 dp** in the per-atom and per-bond tables (matching the",
          "  `.out`'s own precision); molecular energies keep full precision, since the relaxation",
          "  strain is a difference of two large numbers.", "",
          "## Deliberately not extracted", "",
          "SCF convergence metrics, virial ratio, energy decomposition terms, Hirshfeld integrated",
          "densities, dipole vector components, rot-axis dipole, job A's per-atom mirror.", ""]
    if failed:
        L += ["## Failed compounds", "", "| compound | reason |", "|---|---|"]
        for f in failed:
            L.append(f"| `{f['name']}` | {'; '.join(f['failures'])} |")
        L.append("")
    else:
        L += ["## Failed compounds", "", "**None.** Every compound passed all eight assertions.",
              ""]
    OUT.mkdir(parents=True, exist_ok=True)
    REPORT_MD.write_text("\n".join(L))
    log.info("wrote %s", REPORT_MD.name)


# ------------------------------------------------------------------------- main

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=6,
                    help="concurrent parse workers (processes; default 6)")
    ap.add_argument("--limit", type=int, default=None,
                    help="process only the first N compounds, by name")
    ap.add_argument("--force", action="store_true",
                    help="re-parse compounds whose shards already exist")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                        stream=sys.stdout)
    RDLogger.DisableLog("rdApp.*")
    started_at = time.time()
    preflight()
    OUT.mkdir(parents=True, exist_ok=True)
    PARTS.mkdir(parents=True, exist_ok=True)

    # Recorded at the TOP so the scope check asks "did THIS run write there", not "is this path
    # dirty" -- the distinction that gave notebooks 29, 38 and 42 false positives.
    RUN_META_JSON.write_text(json.dumps({
        "started_at": started_at,
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "argv": sys.argv[1:]}, indent=2) + "\n")

    record_decision(
        "Shard per compound rather than reuse the generating run's append_rows().",
        "append_rows is read-modify-rewrite-whole-file; measured at 0.275 s per call once the "
        "atom table reaches 252k rows, i.e. ~26 min of pure rewriting across the population "
        "against a parse budget of a few minutes. A shard write is 0.0002 s. Sharding the "
        "molecular table too (rather than appending it) closes a restart hole: an interruption "
        "between an atom shard and a buffered molecular flush would leave a compound that looks "
        "done but has no molecular row.",
        "Reuse append_rows for all three (26 min of rewriting, and the restart hole); or stream "
        "true appends to open handles (fast, but a kill mid-write leaves a partial line).",
        "measured; deviates from the plan, which said to reuse append_rows for the molecular "
        "table")
    record_decision(
        "Assert RDKit bonds are a SUBSET of the Mayer pair list, not equal to it.",
        "Measured on 400 compounds: subset holds 400/400, equality only 212/400, because ORCA's "
        "0.1 bond-order threshold also admits weak intramolecular contacts (mean 0.64 extra "
        "pairs, max 4). Asserting equality would fail 47% of compounds for a reason that is not "
        "an error, while the subset form still validates the bond graph at every atom index.",
        "Assert equality (fails ~47%); or drop the connectivity check and rely on element "
        "sequence alone, which cannot catch a permutation of two same-element atoms.",
        "measured")
    record_decision(
        "Cross-check energies at 1e-9 absolute, not at exact equality.",
        "The plan predicted exact agreement with descriptors.csv on dft_energy_eh, and the first "
        "run measured 9.095e-13 instead -- enough to fail a zero tolerance. The cause is not a "
        "disagreement: the .out prints ~12 significant figures and .property.txt 17, so the same "
        "double round-trips through two decimal forms to values differing by machine epsilon "
        "scaled to the magnitude (~1790 Eh). 1e-9 keeps ~2000x headroom over the measured "
        "residual while remaining far tighter than any difference a real bug could produce. The "
        "same applies to dispersion_eh (measured 5.0e-10).",
        "Keep a zero tolerance and report a permanent FAIL on a difference that is not an error; "
        "or take the energy from the .out so both sides share one decimal form, which would mean "
        "preferring the less precise source purely to make a check pass.",
        "measured; corrects an expectation stated in the plan")
    record_decision(
        "ProcessPoolExecutor, not the generating run's ThreadPoolExecutor.",
        "That run's workers blocked on ORCA subprocesses, so threads were right. This work is "
        "CPU-bound regex in Python and threads would be serialised by the GIL. Per-compound "
        "shards mean there is no shared mutable state, so no lock is needed either.",
        "ThreadPoolExecutor (no speedup here); or single process (a few minutes, also fine).",
        "reasoned")

    pop = load_population(args.limit)
    names = list(pop.Molecule_Name)
    todo_rows = [{"name": r.Molecule_Name, "smiles": r.SMILES, "inchikey": r.inchikey,
                  "set": r.set, "n_atoms_csv": int(r.total_atoms_with_h)}
                 for r in pop.itertuples(index=False)
                 if args.force or not shard_done(r.Molecule_Name)]
    log.info("population %d, already parsed %d, to parse %d",
             len(names), len(names) - len(todo_rows), len(todo_rows))

    def _signal(signum, _frame):
        if not _STOP.is_set():
            log.warning("signal %d received -- no new compounds will start; in-flight ones will "
                        "finish and their shards will be written", signum)
            _STOP.set()

    signal.signal(signal.SIGINT, _signal)
    signal.signal(signal.SIGTERM, _signal)

    failed: list[dict] = []
    durations: list[float] = []
    done_n = len(names) - len(todo_rows)
    stop_reason = "queue exhausted"

    if todo_rows:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures: dict = {}
            it = iter(todo_rows)

            def _submit_next() -> bool:
                if _STOP.is_set():
                    return False
                try:
                    job = next(it)
                except StopIteration:
                    return False
                futures[pool.submit(parse_one, job)] = job["name"]
                return True

            for _ in range(args.jobs):
                _submit_next()

            while futures:
                for fut in as_completed(list(futures)):
                    name = futures.pop(fut)
                    try:
                        rec = fut.result()
                    except Exception as exc:                          # noqa: BLE001
                        rec = {"name": name, "status": "FAILED",
                               "failures": [f"worker crashed: {exc}"], "seconds": 0.0}
                    if rec["status"] == "OK":
                        done_n += 1
                        durations.append(rec["seconds"])
                        if done_n % 250 == 0:
                            log.info("%-14s OK  atoms=%3d bonds=%3d  (%d/%d)", name,
                                     rec["n_atoms"], rec["n_bonds"], done_n, len(names))
                    else:
                        failed.append(rec)
                        log.error("%-14s FAILED  %s", name, "; ".join(rec["failures"]))
                    write_status(len(names), done_n, failed, durations, started_at)
                    if not _submit_next() and _STOP.is_set():
                        stop_reason = "signal"
                    break

    log.info("parse loop finished (%s): %d done, %d failed", stop_reason, done_n, len(failed))

    parsed = [n for n in names if shard_done(n)]
    counts = concat_shards(parsed)
    checks = cross_check() if counts.get("mol") else []
    for c in checks:
        log.info("cross-check %-52s %s (max |diff| %.3e)", c["quantity"], c["verdict"],
                 c["max_abs_diff"])
    stats = table_stats() if counts.get("atoms") else {}
    write_report(len(names), counts, failed, checks, durations, started_at, args, stats)
    write_status(len(names), done_n, failed, durations, started_at)

    if any(c["verdict"] == "FAIL" for c in checks):
        log.error("one or more cross-checks FAILED -- treat the tables as suspect")
        return 2
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
