# Puixeu_Hayward_2025

### This is the code for the simulations from [Puixeu & Hayward 2025](https://www.biorxiv.org/content/10.1101/2024.11.29.626061v1)

We simulate 2-sex populations with separate sexes under a stabilizing selection, mutation and drift regime in order to examine the relationship between intersex correlation and sexual dimorphism at equilibrium but also during sex-specific adaptation after a shift in optima.
We define three types of simulations (exact or individual-based, Wright-Fisher and Wright-Fisher Hardy-Weinberg simulations) and consider two types of highly-polygenic genetic architectures: approximately infinitesimal and multigenic.
For all details regarding the simulation framework, parameter ranges and simulated scenarios please refer to the Simulations sections of the manuscript.

---

Here is a summary of the content of this respository and how to run the code.

There are two python notebooks with code to reproduce the figures in the manuscript:
- `code_Fig1.ipynb` for main __Figure 1__ and __Supplementary Figures S3-S8__
- `code_Fig2-4.ipynb` for main __Figures 2-4__ and __Supplementary Figures S9-S12__

The simulation results used to run the code in the two notebooks are found in the folders `RESULTS_FIG1/` and `RESULTS_FIG2-4/`, respectively. More concrete details on which results are used for which figures is provided in the relevant parts of each notebook. Detailed information on the meaning of each field and the generation of these results are found in the relevant scripts (see below).

To generate the results in the two folders above, scripts are provided in the `SCRIPTS/` folder:
- `script_fig1.py` to produce results used in Figure 1 and Supplementary Figures S3-S8
- `script_fig2u4.py`: to produce results used in Figure 2, Figure 4C-F and Supplementary Figures S9-S12
- `script_fig3u4.py`: to prpoduce results used in Figure 3 and Figure 4A-B
Concrete details on how to run each script, as well as a description of input parameters and output files are provided in each script.

The `MODULEs` folder contains the core code of the simulations, split into three documented modules:
- `IBS_2sexes_mutation_class.py`: contains objects and methods to implement and keep track of mutations.
- `IBS_2sexes_population_class.py`: contains objects and methods to store and evolve a population.
- `IBS_2sexes_mutation_class.py`: contains objects and methods to run simulations of a given population.
