# Molecule viewer with step editor

**One page, like a small Molden, that opens a Gaussian .log, shows the energy plot and the 3D structure side by side, and lets you remove outlier points so the plot is readable again.**

## IN (first version)

- One HTML file, opens with a double-click, works offline (no internet, no installation).
- On first open it already shows a made-up sample scan.
- Open a file with a button or by dragging it onto the page.
- **Gaussian .log / .out:**
  - Relaxed scan: one point per converged scan point (a 10-step scan gives 11 points).
  - Optimisation / TS search: every optimisation step.
  - Other jobs (e.g. single point): the structure and its energy.
- **Multi-frame .xyz** (energy in the comment line) and plain single .xyz.
- Left: energy plot (relative energy in kcal/mol, the lowest point = 0). Right: 3D structure you can rotate and zoom.
- Click a point on the plot to show its structure. Arrow keys step through the points.
- **Remove points:** select one or more points and remove them. The plot rescales right away. "Undo" brings them back.
- Measure: click 2, 3 or 4 atoms to see the bond length (Å), angle or dihedral (°).
- Atom numbers can be shown on the structure (same numbers as in Gaussian).
- **Save:** cleaned multi-frame .xyz with the remaining points and their energies in the comment lines. It opens in Molden and back in this tool.
- Short hints on the page so colleagues and students understand it.

## Next (step 7, one at a time)

- ORCA output (.out) and ORCA `_trj.xyz`.
- Single-structure formats: .mol / .sdf, .pdb.

## OUT (not now)

- Building or editing atoms (moving, adding, deleting atoms).
- Orbitals, densities, vibrations / frequency animation.
- High-quality images for papers.
- Editing the original .log file itself.

## Sample data

- A made-up Gaussian relaxed-scan .log (a small molecule, ~8 atoms, 10 steps = 11 points) with one clear outlier point, stored as `sample-scan.log` in this folder and built into the page.
- A made-up multi-frame .xyz for testing the .xyz reader.

## Done when

1. The sample scan .log opens with 11 points, and the energies on the plot match the energies in the file.
2. Removing the outlier point immediately shows a clean, readable energy plot with the structure beside it; Undo brings it back.
3. The saved .xyz contains only the remaining points with their energies, and opens again in this tool (you check that it also opens in Molden).
4. A Gaussian optimisation .log shows every optimisation step.
