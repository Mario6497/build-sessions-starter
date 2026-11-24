# What I learned

The agent adds two short sentences after each step.

- **Step 1:** We wrote down my job, weekly tools and AI level in `me.md`. Remembering: the agent uses this to suggest ideas that fit my real work.
- **Step 4:** We wrote `projects/gaussian-input/PLAN.md` for a Gaussian input builder (constrained opt + TS search). Remembering: writing down what is OUT keeps the first version small.
- **Step 5:** The agent built `gaussian-input.html`, one page that turns an .xyz file into a Gaussian .gjf with constraints, scans and my defaults. Remembering: it shows the current bond/angle value next to each constraint, so wrong atom numbers are easy to spot.
- **Step 6:** The agent ran the page in a hidden browser with 7 test cases and compared the .gjf output with the plan. Remembering: the agent can't run Gaussian, so the final test (Gaussian accepts the file) is mine. I ran one job in Gaussian and it worked, so all "done when" points are met.
- **Step 7:** We added a PBS job script next to the input (now saved as .com), built from my real script's layout with PBS memory = %mem + 1GB. Remembering: the script must use Unix line endings, or the cluster can fail to read a file made on Windows.
- **Step 7 (2):** Optimisation now writes opt+freq and then a single point in one file with `--Link1--`, reading geometry and orbitals from the .chk. Remembering: `Geom=AllCheck` takes title, charge and multiplicity from the .chk, so the second part has no coordinates.
- **Step 7 (3):** TS search now also gets the `--Link1--` single point, with its own on/off switch. Remembering: when a feature changes on purpose, an old check can "fail" — the check has to be updated too.
- **New project, step 4:** We wrote `projects/molecule-viewer/PLAN.md` for a small Molden-like viewer that opens a Gaussian .log and lets me remove outlier points from the energy plot. Remembering: a 10-step scan has 11 converged points, and the cleaned result is saved as a multi-frame .xyz that Molden can open.
- **New project, step 5:** The agent built `molecule-viewer.html`. It reads a Gaussian scan/opt .log or a multi-frame .xyz and shows the energy plot, a rotatable 3D structure and a table where I can remove and bring back points. Remembering: the tool reads only the first job in a .log, so `--Link1--` parts (freq, single point) don't show up as extra points.
- **New project, step 6:** The agent ran 15 checks on the viewer in a hidden browser, and all passed. They covered 11 scan points, matching energies, the outlier removal rescaling the plot from 240 to 33 kcal/mol, the saved .xyz reopening, and measurements. Remembering: the agent can't see the page or run Molden, so how it looks and whether Molden opens the saved file are mine to check.
