# ORCA output parse — COMPLETE

Written by `scripts/parse_orca_output.py`. Extracts per-atom, per-bond and molecule-level
quantities from the completed QM run (`outputs/qm_descriptors/orca/`, 11,310 ORCA jobs over 5,655
compounds) into three tidy tables. **5,655 compounds parsed, 0 failures**, in well under a minute.

Parsing only: no ORCA was re-run, and nothing under `outputs/qm_descriptors/` is written or
modified — `descriptors.csv` and `hirshfeld_charges.csv` are read for cross-validation and nothing
else. Notebook `43` Part A puts this on record in reader-facing form.

## Files

| File | Rows | Holds |
|---|---|---|
| `qm_atoms.csv` | 252,237 | One row per atom: Mulliken/Loewdin/Hirshfeld charges, Mayer NA and VA, Mulliken and Loewdin s/p/d shell aggregates, and the optimised geometry (x, y, z in Angstrom). |
| `qm_bonds.csv` | 269,256 | One row per printed Mayer pair: both atom indices, both elements, the bond order, and `is_rdkit_bond`. |
| `qm_molecular_raw.csv` | 5,655 | Total energy, HOMO/−1/−2, LUMO/+1, gap, dipole, D3 dispersion, three rotational constants, basis count, SCF cycles, plus GFN2 optimisation cycles, final gradient norm and relaxation strain. |
| `parse_report.md` | — | Source used per quantity, the validation set, the cross-checks, and every failed compound (there are none). |
| `parse_status.json`, `decisions_taken.json`, `run_meta.json` | — | Run state, the four decisions taken, and the run-start timestamp. |

`parts/` holds one shard per compound per table (16,965 files) and is **gitignored** — it exists
only so the parse is restartable per compound, and the three tables are concatenated from it at the
end. **It is safe to delete once the tables are built, and was deleted after this run**: it is pure
intermediate, and a fresh parse rebuilds it in under a minute. Deleting it is worth knowing about
because 16,965 files is enough to make an editor's file watcher visibly unhappy even though git
ignores them. The script recreates the directory on its next run and simply re-parses every
compound.

The three tables are tracked because the 8.8 GB ORCA tree they derive from is not, which makes them
the only durable form of this data.

Atom indices are **0-based and are RDKit's own**, which is what makes it valid to join a
SMARTS-derived atom index straight into `qm_atoms.csv` — see Validation below.

## Source precedence

The parse prefers each job's `.property.txt` over its `.out` wherever a quantity exists in both.
That is the opposite of the obvious approach, and it was the most consequential design choice:
`.property.txt` stores values as typed indexed arrays (`&Type "ArrayOfDoubles", &Dim (n,m)`) in
which an element symbol never sits adjacent to a number.

Measured against both jobs' full field inventory: **13 of 19 quantities come from a
`.property.txt`** (10 from job B's, 3 from job A's), **5 from job B's `.out`** because no equivalent
section exists there at all — the two reduced-orbital blocks, the orbital ladder, the rotational
constants and the SCF cycle count — and 1 from the optimised-geometry `.xyz`. The effect is that
the trap below can only reach two of the nineteen.

**Job A's `.out` is never opened**: 5.48 GB, 64% of the whole tree, and a lower-theory mirror of
job B. Everything needed from job A is in its `.property.txt` as named fields.

## The trap

In the reduced-orbital sections ORCA writes fixed-width fields, so a two-character element symbol
runs into the shell label with no separator:

```
  6 C s       :     2.925169
  7 Cls       :     5.940711      <-- "Cls"
```

Any `([A-Z][a-z]?)\s+(\S+)` pattern **silently drops those atoms** and returns a short block with
no error at all. Of the nine elements in this library only Cl and Br have two-character symbols —
C, N, O, S, P, F, I and H are all safe, which is exactly why a parser written against one example
looks fine. The affected population is **815 compounds, 14.4%**.

That block is therefore parsed by fixed character position (`[0:4]` index, `[4:6]` element,
`[6:14]` shell, `[15:30]` value, `[30:]` aggregate), and the parsed table carries its own
acceptance test: **924 Cl/Br atoms across exactly 815 compounds**, matching the count notebook `41`
derived independently from SMILES. Two different methods on two different files.

## Validation

Eight assertions per compound; a failure records the compound and its reason and writes **no rows
at all** for it, because a short block is exactly what the trap produces silently. All 5,655
passed.

The two worth knowing, because each was checked against a deliberate negative control rather than
trusted:

- **Atom order.** ORCA's order came from `Chem.AddHs(Chem.MolFromSmiles(row.SMILES))` with no
  reordering, so RDKit's atom index *is* ORCA's. The raw `SMILES` column reproduces it on all
  5,655 compounds; **`canonical_smiles`, which sits in the same CSV, manages 4,212 — it is wrong on
  1,443 compounds, a quarter of the library.** Anything joining an atom index into these tables
  must read the raw column.
- **Bond topology.** Element sequence cannot see a permutation of two atoms of the *same* element,
  so RDKit's bond set is also asserted to be a **subset** of the Mayer pair list — subset, not
  equality, because ORCA's 0.1 threshold also admits weak contacts (mean 0.64 extra pairs per
  compound). Negative control: swapping two carbons leaves the element sequence identical and the
  bond-subset check still fails, with 7 bonds missing.

Charge-sum tolerances are **per scheme**, each justified by measurement rather than one blanket
number: Mulliken and Loewdin come out of `.property.txt` exact to floating point (max residual
1.5e-12, tolerance 1e-6), while only Hirshfeld carries real numerical-integration error (max
5.9e-4, tolerance 5e-3). A single 5e-3 tolerance would hide a genuine Mulliken bug by six orders of
magnitude.

## Cross-check against the generating run

`descriptors.csv` and `hirshfeld_charges.csv` were produced by independent code reading the same
ORCA files, so these are real external checks rather than self-consistency. **All 12 pass, 9 at
exactly `0.000e+00`** — HOMO, LUMO, gap, basis count, SCF cycles, atom count, formal charge, and
all 252,237 Hirshfeld charges and element symbols.

The three that are not exactly zero are explained and none is a disagreement: `dft_energy_eh` at
9.095e-13 and `dispersion_eh` at 5.0e-10 are machine epsilon at ~1790 Eh between the `.out`'s ~12
significant figures and `.property.txt`'s 17, and `dipole_debye` at 5.9e-04 is the `.out`'s own
6-decimal rounding. That the energies differ by *this little* is itself the evidence the two
sources carry the same number.

## Known limits

- **The orbital ladder is truncated by ORCA** — both jobs print `*Only the first 10 virtual
  orbitals were printed.`, and job A's virtual count varies 10–13. All occupied orbitals are
  present, so HOMO−n is reliable, but **nothing beyond LUMO+1 is extracted** and nothing beyond it
  is recoverable without re-running.
- **Only the s/p/d aggregates** are kept from the reduced-orbital blocks. The individual pz/px/py
  and five d components depend on the orientation the molecule happens to have in the file and are
  not molecular properties; the shell sums are rotationally invariant.
- **Deliberately not extracted**: SCF convergence metrics, virial ratio, energy decomposition
  terms, Hirshfeld integrated densities, dipole vector components, rot-axis dipole, and job A's
  per-atom mirror. Quadrupole, polarizability and thermochemistry are **not available at all** —
  they were never requested in the keyword line, so they would need a re-run, not a re-parse.
- Floats are written to 6 decimal places in the per-atom and per-bond tables, matching the `.out`'s
  own precision; molecular energies keep full precision, since the relaxation strain is a
  difference of two large numbers.

## Re-running

```bash
conda run -n cyp-admet-v2 python scripts/parse_orca_output.py --jobs 6
```

Restartable: a compound is done iff all three of its shards exist and are non-empty, so an
interruption loses at most the compounds in flight. `--force` re-parses shards that already exist;
`--limit N` takes the first N compounds. SIGTERM drains cleanly rather than killing in-flight work.
Exit code is 0 on success, 1 if any compound failed, 2 if a cross-check failed.
