# Tuned_Inhibition_PLOS_Comp_Bio_2021
Code for the Tuned Inhibition simulations of Maniscalco, Odegaard, Grimaldi, Cho, Basso, Lau, &amp; Peters 2021 PLOS Comp Bio
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008779

#### Simulating individual trials
`TI_sim_trial.m` is the function used for simulating individual trials. Please refer to the detailed help comments in the file for details on usage.

`TI_sim_trial_example.m` can be used to plot example output for a single simulated trial.

#### Simulating full experiments
Subfolders correspond to different experiments simulated in Maniscalco et al. 2021, as follows.

`MPL2016` : Code for simulating the data of Maniscalco, Peters, & Lau 2016 (Figure 2 in Maniscalco et al. 2021)

`KML2015_1A` : Code for simulating the data of Koizumi, Maniscalco, & Lau 2015 expt 1A (Figures 3-4 in Maniscalco et al. 2021)

`KML2015_2B` : Code for simulating the data of Koizumi, Maniscalco, & Lau 2015 expt 2B (Figures 5-6 in Maniscalco et al. 2021)

In each folder, there is a corresponding Word document (e.g. `MPL2016_simulation_instructions.docx` in the `MPL2016` folder) that contains detailed instructions on how to conduct every step of the simulations.
