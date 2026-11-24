# Gaussian input builder

**One page that opens an .xyz file, lets you pick constrained atoms or set up a TS search, and writes a ready Gaussian input file without typos.**

## IN (first version)

- One HTML file, opens with a double-click, works offline.
- Button to open an `.xyz` file. On first open it already shows a made-up sample molecule.
- Atoms shown as a numbered list (number, element, x, y, z), so atom numbers are easy to read.
- Five calculation types:
  - **Optimisation:** `Opt=Tight Freq` (regular, no constraints).
  - **Constrained optimisation:** `Opt=ModRedundant Freq` with frozen bonds (B), angles (A) and dihedrals (D), marked with `F`.
  - **TS search:** `Opt=(TS,CalcFC,NoEigenTest) Freq`, optionally with constraints too (adds `ModRedundant`).
  - **Relaxed scan:** `Opt=ModRedundant` with one or more `S` lines (number of steps, step size), plus optional frozen (`F`) coordinates. No `Freq` by default.
  - **Single point:** energy only, no `Opt`/`Freq`.
- Convergence choice for every optimisation type: Loose · Normal · Tight. Default: Tight for regular optimisation, Normal for constrained opt, TS search and scan.
- Add constraints by typing atom numbers (2 for a bond, 3 for an angle, 4 for a dihedral); the tool shows the current value (Å or °) so you can see you picked the right atoms. For a scan it also shows the final value after all steps.
- Settings with defaults, so you only change what you need:
  - Method / basis for opt, TS and scan: `wB97XD/6-31G**`
  - Method / basis for single point: `wB97XD/6-311++G(d,p)`
  - Solvent for all: `SCRF=(PCM,Solvent=Acetonitrile)` (can be changed or switched off)
  - Charge / multiplicity: `0 1`
  - Resources: `%nprocshared=16`, `%mem=5GB`, `%chk=<file name>.chk`
  - Everything can be overwritten by typing; button "Save as my defaults" (remembered in this browser only).
- Live preview of the `.gjf` file, plus "Copy" and "Download .gjf" buttons.
- Short hints next to each field so colleagues and students understand it.

## OUT (not now)

- Running Gaussian.
- 3D molecule viewer.
- Cluster job scripts (SLURM/PBS) — good candidate for later.
- ORCA — maybe later.

## Sample data

- Made-up small molecule, e.g. a simple reaction guess structure with ~8 atoms, stored as `sample.xyz` in this folder and built into the page.

## Done when

1. For the sample molecule, the generated file matches what you would write by hand (route line, title, charge/multiplicity, coordinates, blank lines, constraint lines).
2. Gaussian accepts the file without input errors (you test this on your machine; the agent cannot run Gaussian).
3. Atom-number constraints that don't exist (e.g. atom 12 in an 8-atom molecule) give a clear warning instead of a broken file.
