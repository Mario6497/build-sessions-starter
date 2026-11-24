# What I learned

The agent adds two short sentences after each step.

- **Step 1:** We wrote down my job, weekly tools and AI level in `me.md`. Remembering: the agent uses this to suggest ideas that fit my real work.
- **Step 4:** We wrote `projects/gaussian-input/PLAN.md` for a Gaussian input builder (constrained opt + TS search). Remembering: writing down what is OUT keeps the first version small.
- **Step 5:** The agent built `gaussian-input.html`, one page that turns an .xyz file into a Gaussian .gjf with constraints, scans and my defaults. Remembering: it shows the current bond/angle value next to each constraint, so wrong atom numbers are easy to spot.
- **Step 6:** The agent ran the page in a hidden browser with 7 test cases and compared the .gjf output with the plan. Remembering: the agent can't run Gaussian, so the final test (Gaussian accepts the file) is mine. I ran one job in Gaussian and it worked, so all "done when" points are met.
- **Step 7:** We added a PBS job script next to the input (now saved as .com), built from my real script's layout with PBS memory = %mem + 1GB. Remembering: the script must use Unix line endings, or the cluster can fail to read a file made on Windows.
