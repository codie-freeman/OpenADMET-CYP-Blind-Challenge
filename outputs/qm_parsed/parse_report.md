# ORCA output parse report

Generated 2026-10-05T19:18:00.975060+00:00  
`scripts/parse_orca_output.py --jobs 6`

Parsing only. No ORCA re-run. Nothing under `outputs/qm_descriptors/` was written or
modified; `descriptors.csv` and `hirshfeld_charges.csv` were read for cross-validation
only.

## Result

- compounds in population: **5655**
- parsed successfully: **5655**
- failed: **0**
- wall time: **0.2 min**

| table | rows | size |
|---|---|---|
| `qm_atoms.csv` | 252,237 | 31.9 MB |
| `qm_bonds.csv` | 269,256 | 9.3 MB |
| `qm_molecular_raw.csv` | 5,655 | 1.4 MB |

## Source used per quantity

`.property.txt` is preferred wherever a quantity exists there: it stores typed indexed
arrays, so an element symbol never sits adjacent to a number and the fixed-width trap
cannot apply. Of the 19 quantities below, **13 come from a `.property.txt`** (10 from
job B's, 3 from job A's), **5 from job B's `.out`** because no `.property.txt`
equivalent exists, and 1 from the optimised-geometry `.xyz`. Only the two
reduced-orbital blocks are exposed to the trap, and they are parsed by fixed character
position.

Job A's `.out` is never opened (5.48 GB, 64% of the tree).

| quantity | source | section / field |
|---|---|---|
| Mulliken atomic charge | `B .property.txt` | `$SCF_Mulliken_Population_Analysis &AtomicCharges` |
| Loewdin atomic charge | `B .property.txt` | `$SCF_Loewdin_Population_Analysis &AtomicCharges` |
| Hirshfeld atomic charge | `B .property.txt` | `$SCF_Hirshfeld_Population_Analysis &AtomicCharges` |
| Mayer NA, VA | `B .property.txt` | `$SCF_Mayer_Population_Analysis &NA, &VA` |
| Mayer bond orders + indices | `B .property.txt` | `$SCF_Mayer_Population_Analysis &BondOrders, &components` |
| atomic numbers (validation) | `B .property.txt` | `&ATNO in each population section` |
| total energy | `B .property.txt` | `$Single_Point_Data &FinalEnergy` |
| dipole magnitude | `B .property.txt` | `$SCF_Dipole_Moment &dipoleMagnitude (a.u., converted)` |
| D3 dispersion | `B .property.txt` | `$VdW_Correction &vdW` |
| basis count | `B .property.txt` | `$Calculation_Info &NumOfBasisFuncts` |
| Mulliken s/p/d aggregates | `B .out (FIXED COLUMNS)` | `MULLIKEN REDUCED ORBITAL CHARGES` |
| Loewdin s/p/d aggregates | `B .out (FIXED COLUMNS)` | `LOEWDIN REDUCED ORBITAL CHARGES` |
| HOMO, HOMO-1, HOMO-2, LUMO, LUMO+1 | `B .out` | `ORBITAL ENERGIES (virtuals truncated by ORCA)` |
| rotational constants | `B .out` | `Rotational constants in cm-1` |
| SCF cycles | `B .out` | `SCF CONVERGED AFTER` |
| optimised geometry | `A .xyz` | `read_orca_xyz(), reused from the generating run` |
| GFN2 opt cycles | `A .property.txt` | `count of $Geometry blocks, minus 1` |
| GFN2 final gradient norm | `A .property.txt` | `last $SCF_Nuc_Gradient &gradNorm` |
| GFN2 relaxation strain | `A .property.txt` | `first vs last $SCF_Energy &totalEnergy` |

## The trap: acceptance test on the real output

The reduced-orbital blocks are the only place a two-character element symbol sits
against a shell label (`7 Cls       :`), and a whitespace parser drops those atoms
silently. The parsed tables therefore carry their own acceptance test:

- **924 Cl/Br atoms parsed**, across **815 compounds** (761 Cl, 163 Br). Every one would be missing under whitespace splitting, and the compound count matches the 815 measured independently in the audit.
- elements present: Br, C, Cl, F, H, I, N, O, P, S -- all nine in the library, plus H.
- **p/d aggregates are null on exactly the 113,596 hydrogens and nowhere else**: True. def2-SV(P) puts no polarisation function on H, so this is the expected pattern and a second check that no atom was skipped.
- bonds: **265,364** are real RDKit bonds and **3,892** are weak contacts above ORCA's 0.1 threshold (0.69 per compound, against 0.64 measured on the 400-compound sample).

## Validation

Eight assertions per compound. A failure records the compound and its reason and the
compound is skipped entirely -- no partial rows are written, because a short block is
exactly what the trap produces silently.

1. row count per section == RDKit atom count, checked separately for Mulliken, Loewdin,
   Hirshfeld, Mayer NA, Mayer VA, both reduced-orbital blocks and the geometry; plus
   agreement with `descriptors.csv`'s `total_atoms_with_h`
2. charge sums == formal charge, per scheme: Mulliken/Loewdin 1e-6, Hirshfeld 5e-3
   (measured maxima 1.5e-12 and 5.9e-4 respectively -- one blanket tolerance would hide
   a Mulliken bug by six orders of magnitude)
3. atom indices contiguous 0..n-1; bond indices within range
4. element at each ORCA index == element at the same RDKit index, from the raw `SMILES`
   column (not `canonical_smiles`, which reproduces the order on only 310 of 400
   sampled compounds)
5. `&ATNO` atomic numbers agree with RDKit's -- an independent form of (4)
6. RDKit bond set is a **subset** of the Mayer pair list, not equal to it: the Mayer
   list is thresholded at bond order > 0.1 and also carries weak intramolecular
   contacts (measured mean 0.64 extra pairs, max 4). Equality would fail on 47% of
   compounds for a reason that is not an error.
7. `&NAtoms` consistent across every section that reports it
8. reduced-orbital coverage: every heavy atom has s/p/d, every hydrogen has s only
   (def2-SV(P) puts no polarisation function on H)

## Cross-check against the generating run

`descriptors.csv` and `hirshfeld_charges.csv` were produced by independent code reading
the same files, so these are real external checks rather than self-consistency.

| quantity | n | max abs diff | tolerance | verdict |
|---|---|---|---|---|
| `homo_ev vs descriptors.csv::dft_homo_ev` | 5,655 | 0.000e+00 | 0.0e+00 | **PASS** |
| `lumo_ev vs descriptors.csv::dft_lumo_ev` | 5,655 | 0.000e+00 | 0.0e+00 | **PASS** |
| `gap_ev vs descriptors.csv::dft_gap_ev` | 5,655 | 0.000e+00 | 1.0e-09 | **PASS** |
| `total_energy_eh vs descriptors.csv::dft_energy_eh` | 5,655 | 9.095e-13 | 1.0e-09 | **PASS** |
| `n_basis vs descriptors.csv::n_basis` | 5,655 | 0.000e+00 | 0.0e+00 | **PASS** |
| `scf_cycles vs descriptors.csv::dft_scf_cycles` | 5,655 | 0.000e+00 | 0.0e+00 | **PASS** |
| `dispersion_eh vs descriptors.csv::dft_dispersion_eh` | 5,655 | 4.999e-10 | 1.0e-09 | **PASS** |
| `dipole_debye vs descriptors.csv::dft_dipole_debye` | 5,655 | 5.908e-04 | 2.0e-03 | **PASS** |
| `n_atoms vs descriptors.csv::total_atoms_with_h` | 5,655 | 0.000e+00 | 0.0e+00 | **PASS** |
| `formal_charge vs descriptors.csv::formal_charge` | 5,655 | 0.000e+00 | 0.0e+00 | **PASS** |
| `q_hirshfeld vs hirshfeld_charges.csv::charge` | 252,237 | 0.000e+00 | 5.0e-07 | **PASS** |
| `element vs hirshfeld_charges.csv::element` | 252,237 | 0.000e+00 | 0.0e+00 | **PASS** |

## Conventions worth knowing

- **Energies cross-check to ~5e-13, not to zero.** The `.out` prints ~12 significant
  figures and `.property.txt` 17, so the same double round-trips to values differing by
  machine epsilon scaled to the magnitude (~1790 Eh). That the difference is this small
  is itself the evidence the two sources carry the same number.
- **Dipole** comes from `.property.txt` in a.u. and is converted with 2.5417464519
  (CODATA). Against the existing `dft_dipole_debye` this agrees to ~2e-4, which is the
  `.out`'s own 6-decimal rounding -- the `.property.txt` route is the more precise.
- **`gfn2_opt_cycles`** is `n_geometry - 1`: a 26-cycle optimisation writes 27
  `$Geometry` blocks, the extra being the final evaluation at the stationary point.
  This is the number the `.out`'s own cycle counter agrees with.
- **`gfn2_relaxation_strain_*`** is `E_first - E_last`, so **positive means energy was
  released** during relaxation (the normal case).
- **Nothing beyond LUMO+1** is extracted. ORCA truncates the virtual spectrum
  (`*Only the first 10 virtual orbitals were printed.`) and job A's count of them
  varies 10-13. All occupied orbitals are present, so HOMO-1 and HOMO-2 are reliable.
- **Only the s/p/d aggregates** are kept from the reduced-orbital blocks. The pz/px/py
  and five d components depend on the molecule's orientation in the file and are not
  molecular properties; the shell sums are rotationally invariant.
- **`is_rdkit_bond`** separates real bonds from the weak contacts that ORCA's 0.1
  threshold also admits. Beyond the brief, but nothing else distinguishes them.
- **Hirshfeld `spin` is dropped**: identically 0.0 for all 5,655 closed-shell singlets.
- **Floats are written to 6 dp** in the per-atom and per-bond tables (matching the
  `.out`'s own precision); molecular energies keep full precision, since the relaxation
  strain is a difference of two large numbers.

## Deliberately not extracted

SCF convergence metrics, virial ratio, energy decomposition terms, Hirshfeld integrated
densities, dipole vector components, rot-axis dipole, job A's per-atom mirror.

## Failed compounds

**None.** Every compound passed all eight assertions.
