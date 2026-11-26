# Molecule studio (local Molden-like viewer and editor)

**A Python desktop program on my Ubuntu Linux computer that opens Gaussian and ORCA outputs from my synced folder and does most of what Molden does, plus a few extras: converged-only scans, an editor with undo, my own fragments and displacement along the imaginary mode for IRC.**

Separate from `projects/molecule-viewer/` (that one stays as it is). It reuses ideas and code from `molview.py`.

## IN (first version, built in small steps)

**Program**
- Python + tkinter desktop app for my PC: Ubuntu 24.04, Python 3.12, tkinter 8.6 (installed with `sudo apt install python3-tk`).
- Python standard library only: no numpy or other extra packages, so nothing more to install. Anything extra is asked first.
- Starts with `python3 molstudio.py`, everything in this folder (`projects/molecule-studio/`), Unix line endings.
- Opens Gaussian .log/.out, ORCA .out / `_trj.xyz`, and .xyz (single and multi-frame). Open with a button or from the command line (`python3 molstudio.py file.log`).
- On first start it shows made-up sample data.

**View**
- Rotatable, zoomable 3D structure; measure bond, angle, dihedral.
- **Atom numbers on / off** (same numbers as in the program).
- **Display model:** ball-and-stick, sticks, spacefill (van der Waals spheres), wireframe.

**Geometry and energy plots**
- Optimisations and TS searches: every step on the energy plot.
- **Scans: only converged structures and their energies** (a 10-step scan gives 11 points).
- Click a point to see its structure; remove outlier points, Undo brings them back.

**Vibrations / IR**
- List of frequencies; imaginary ones clearly marked.
- Animate a normal mode; IR spectrum plot.
- **Displace along the imaginary mode** (forward and backward, chosen step size) and save the two structures as .xyz, ready to use as IRC or follow-up starting structures.

**Editor (Z-matrix / structure)**
- Change bond lengths, angles, dihedrals; add or delete atoms.
- **Undo / Redo** for every edit.
- **Fragment library** with more organic fragments (e.g. methyl, phenyl, carbonyl, hydroxyl, amine, carboxyl, ester, amide, nitro, common rings).
- **My fragments:** save a selected part of a structure as my own fragment, kept inside the program.
- Save the structure as .xyz.

## OUT (not now)
- Job builder and sending jobs to the cluster: a **separate program** later.
- Orbitals and densities.
- Paper-quality images.
- Job queue monitoring.
- Editing the original output files.

## Sample data (all made up)
- A Gaussian relaxed scan .log with non-converged steps in between and one outlier.
- A Gaussian TS + frequency .log with one imaginary frequency.
- An ORCA optimisation .out and `_trj.xyz`.
- A multi-frame .xyz.

## Done when
1. **Scan & TS:** the sample scan shows only converged points with the energies from the file; the sample TS shows its imaginary mode animated and the IR spectrum.
2. **IRC start:** displacing along the imaginary mode saves correct forward and backward structures (atoms moved by the chosen step along the mode).
3. **Editor + undo:** I add a fragment, change a bond, undo it, redo it, and save my own fragment that is still there after restarting the program.
4. **View options:** atom numbers turn on and off, and all four display models work on the sample structures.

## Build order (one step at a time, check after each)
1. Window + open files + 3D view + display models + atom numbers on/off + energy plot (converged-only scans).
2. Vibrations, IR spectrum, displace along imaginary mode.
3. Editor with undo/redo and fragment library + my fragments.
