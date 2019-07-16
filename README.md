# Dissertation

This contains the MATLAB code (excluding its dependencies) for generating and analyzing all the simulations and data used in the creation of Chapter 2 of my Dissertation.

### Dependencies: 

- The "coupling_addition" branch of my PERSONAL fork of DynaSim installed, [available here](https://github.com/asoplata/DynaSim/tree/coupling_addition) . See the instructions for installation at that repo, but only install my custom version, not the default version.
- The "dynasim-extended-benita-model" model files for DynaSim installed, [available here](https://github.com/asoplata/dynasim-extended-benita-model) . See the instructions for installation at that repo.
- A "recent" version of Matlab (both 2017a and 2018b were used originally)
- A "modern" desktop/laptop CPU (i.e. built in the 2010's)
- At LEAST 16GB of RAM (preferably more)
- At least 80GB of disk space
- Many hours (possibly up to 48 or more) of time to run the simulations (probably only 20 minutes 
each), run the analysis (long), and save the data (long). Depending on your
CPU and disk write speed, the entire data generation process may take 3 to 
6 hours or more **per simulation**. Of course, the time to generate each individual plot after loading
the proper data should take only a couple seconds.

### How to use:

1. First, run "run_sims_and_calc_PAC_signals.m", which may take up to 48 hours of continuous computation time. If MATLAB crashes in the middle, just rerun the remaining simulations by commenting out the previous ones.
2. Second, run "plotting_script.m", which will generate almost all of the BASE plots used in the original figures. That's it! Note that for the figures in the actual manuscript, they were prettified and organized using Inkscape. After the Dissertation is published I will likely add the SVG files I used in this repo as well.