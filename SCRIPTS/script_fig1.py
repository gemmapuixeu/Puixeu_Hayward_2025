'''
(Modified from IBS_2sexes_runAna1_24.py)

This script is used to run various sets of simulations, used to produce results for Figure 1 and Supplementary Figures S3-S5, S7-S8.
To obtain results for Figure S6 one can simply use the run_simulations() function (as a single replicate is enough) and change the generations parameter in line 101 to collect the output variables across the desired amount of generations.

Essentially, we let the population burn in for 'Ngen' generations and then report a few variables for 1 generation across 'reps' replicates (in the output file).

To run the script in the command line, type:
python script_fig1.py N Ngen VA oF oM E2Ns rfm reps simtype path

A usage example would be:
python script_fig1.py 1000 100 9 0 0 1 0.5 200 WFHW full_path_to_current_directory/RESULTS/

### INPUT

The arguments it takes correspond to:
    - N: population size
    - Ngen: the population will be burned in for Ngen*N generations
    - VA: genetic variance
    - oF, oM: sex-specific optima
    - E2Ns: parameter modulating the average squared effect size of incoming mutations. In the manuscript we use E2Ns=1 for approximately infinitesimal and E2Ns=16 for multigenic genetic architectures, but any other value can be used.
    - rfm: intersex correlation
    - reps: number of replicates
    - simtype: the type of simulations we want to run. Options are: exact_fertility, exact_viability, WF, WFHW, WFxd2 (see Supplementary Section 4). WFxd2 simulations correspond to WF simulations.
    - path: the path where the results will be stored. It should be of the type full_path_to_current_directory/results_folder/

# OUTPUT

It produces an output file containing various fields, each row corresponding to a separate replicate. Each field:
    - rfm: intersex correlation
    - SDv: empirical variance-normalized variance in SD
    - meanF: average female trait
    - meanM: average male trait
    - meanA: average phenotype across both sexes
    - varF: female phenotypic variance
    - varM: male phenotypic variance
    - varT: total penotpyic variance (see Supplementary Section 6)
    - varA: average variance across sexes
    - meanF_FIX: contribution to the average phenotype from fixed mutations (female)
    - meanM_FIX: contribution to the average phenotype from fixed mutations (male)
    - meanF_FIXshared: contribution to the aveage phenotype from fixed mutations with equal effects in both sexes (female)
    - meanM_FIXshared: contribution to the aveage phenotype from fixed mutations with equal effects in both sexes (male)
    - meanF_FIXspec: contribution to the aveage phenotype from fixed mutations with sex-specific effects (female)
    - meanM_FIXspec: contribution to the aveage phenotype from fixed mutations with sex-specific effects (male)
    - meanF_SEGshared: contribution to the aveage phenotype from segregating mutations with equal effects in both sexes (female)
    - meanM_SEGshared: contribution to the aveage phenotype from segregating mutations with equal effects in both sexes (male)
    - meanF_SEGspec: contribution to the aveage phenotype from segregating mutations with sex-specific effects (female)
    - meanM_SEGspec: contribution to the aveage phenotype from segregating mutations with sex-specific effects (male)
    - varF_shared: contribution to the phenotypic variance from segregating mutations with equal effects in both sexes (female)
    - varM_shared: contribution to the phenotypic variance from segregating mutations with equal effects in both sexes (male)
    - varF_spec: contribution to the phenotypic variance from segregating mutations with sex-specific effects (female)
    - varM_spec: contribution to the phenotypic variance from segregating mutations with sex-specific effects (male)
All elements above are calculated using allele frequencies and sex-specific allele effects. The three below correspond to variances
computed across individuals in individual-based (exact) simulations (See Supplementary Sections 4 and 5).
    - varF_C: female variance computed across individuals in exact simulations
    - varM_C: male variance computed across individuals in exact simulations
    - varA_C: average variance across sexes computed across individuals in exact simulations (not corresponding to the total variance, see Supplementary Section 6).
'''

# IMPORT

# Import the relevant modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import re
import timeit
import datetime

# We also need to import our modules, which include the functions used in our simulations
import sys
sys.path.append('../MODULES/')
from IBS_2sexes_simul_class import SimulatePopulation, Analytics

# FUNCTIONS

'''Function to empirically calculate the variance-normalized variance in sexual dimorphism'''
def calc_SDbyVA (s):
    return float((np.array(s.meanS_f[-1])-np.array(s.meanS_m[-1]))**2/(0.5*(np.array(s.varS_f[-1])+np.array(s.varS_m[-1])+1e-5)))

'''Function to run our simulations'''
def run_simulations (N, Ngen, oF, oM, VA, rfm, E2Ns, simtype):
    start_time = timeit.default_timer()    
    # Initialize the analytics
    a0 = Analytics(N=N, U=0, optF=oM, optM=oM, q=1-rfm, E2Ns=E2Ns)
    # Get mutation rate given variance
    U = a0.get_U_from_VA(VA)
    # Run the model with the new mutation rate:
    
    # Initiate the population given the simluation type
    s = SimulatePopulation(N=N, U=U, optF=oF, optM=oM, q=1-rfm, E2Ns=E2Ns)
    if simtype=="exact_fertility": s.initiate_basic_population_classExact_fertilitySelection()
    elif simtype=="exact_viability": s.initiate_basic_population_classExact_viabilitySelection()
    elif simtype=="WF": s.initiate_basic_population_classWF()
    elif simtype=="WFHW": s.initiate_basic_population_classWF_HW()
    elif simtype=="WFxd2": s.initiate_basic_population_classWF_xd2()
    else: raise Exception ("Invalid simtype: exact_fertility, exact_viability, WF, WFHW, WFxd2")

    # Let the population burn-in for Ngen*N generations (so it attains balancing selection-mutation-drift balance)
    # This burn-in function does not record any variables
    s.burn_basic_population_to_reach_steady_state(burn_time_N=Ngen)
    # Run the simulations for one generation, as we want to collect variables for a single generation after the population has burned in
    generations=1
    s.run_simulations(generations)

    # Compute standardized SD (following function above)
    SDv = calc_SDbyVA(s)
    # Collect all the variables we are interested in
    df = {'rfm':rfm, 'SDv':SDv, 'meanF':s.meanS_f[-1], 'meanM':s.meanS_m[-1], 'meanA':0.5*(s.meanS_f[-1]+s.meanS_m[-1]),
          'varF':s.varS_f[-1], 'varM':s.varS_m[-1], 'varT': s.varS_t[-1], 'varA':0.5*(s.varS_f[-1]+s.varS_m[-1]),
          'meanF_FIX':s.fixedBckg_f[-1], 'meanM_FIX':s.fixedBckg_m[-1], 'meanF_FIXshared':s.fixedBckg_f_shared[-1], 
          'meanM_FIXshared':s.fixedBckg_m_shared[-1], 'meanF_FIXspec':s.fixedBckg_f_spec[-1], 'meanM_FIXspec':s.fixedBckg_m_spec[-1], 
          'meanF_SEGshared': s.meanSeg_f_shared[-1], 'meanM_SEGshared': s.meanSeg_m_shared[-1],'meanF_SEGspec': s.meanSeg_f_spec[-1], 'meanM_SEGspec': s.meanSeg_m_spec[-1],
          'varF_shared':s.varS_f_shared[-1], 'varM_shared':s.varS_m_shared[-1], 'varF_spec':s.varS_f_spec[-1], 'varM_spec':s.varS_m_spec[-1]}

    if 'exact' in simtype:
        df['varF_C'] = s.varC_f[-1]
        df['varM_C'] = s.varC_m[-1]
        df['varA_C'] = 0.5*(s.varC_f[-1]+s.varC_m[-1])
    else:
        df['varF_C'] = np.nan
        df['varM_C'] = np.nan
        df['varA_C'] = np.nan

    print ('Replicate runtime:', timeit.default_timer() - start_time)
    
    return df

''' Function to save the simulation results into a file
If the given name already exists (meaning, we have already run simulations with the same set of parameters), it creates a new file
(with a new version) to ensure we are not overwriting previous results.'''
def save_version (df, path, name):
    c=1
    tosearch = path+name+"_v"+str(c)+".csv"
    files_present = os.path.isfile(tosearch) 
    while files_present:
        c+=1
        tosearch = path+name+"_v"+str(c)+".csv"
        files_present = os.path.isfile(tosearch) 
    df.to_csv (tosearch, index=None)
    
    return tosearch
    
'''Function to record results for various replicates
and store them in a separate file'''
def run_simulations_replicates (N, Ngen, oF, oM, VA, E2Ns, rfm, reps, simtype, path):

    # Create an empty dataframe
    df = pd.DataFrame(index=range(reps), columns = ['rfm', 'SDv', 'meanF', 'meanM', 'meanA', 'varF', 'varM', 'varT', 'varA',
              'meanF_FIX', 'meanM_FIX', 'meanF_FIXshared', 'meanM_FIXshared', 'meanF_FIXspec', 'meanM_FIXspec', 
              'meanF_SEGshared', 'meanM_SEGshared', 'meanF_SEGspec', 'meanM_SEGspec',
              'varF_shared', 'varM_shared', 'varF_spec', 'varM_spec', 'varF_C', 'varM_C', 'varA_C'])
    
    # Collect the relevant variables across 'reps' replicates
    for i in range(reps):
        res = run_simulations(N, Ngen, oF, oM, VA, rfm, E2Ns, simtype)
        df.iloc[i] = res

    # Define the name of the file and save the results as a separate version
    name = "simtypes_%iN_%iNgen_%iVA_%iE2Ns_%ireps_%.3frfm_%.3foF_%.3foM_%s"%(N,Ngen,VA,E2Ns,reps,rfm,oF,oM, simtype)
    save_version(df, path, name)


# RUN CODE

# Collect the simulation parameters as input arguments
N = int(sys.argv[1])
Ngen = int(sys.argv[2])
VA = float(sys.argv[3])
oF = float(sys.argv[4])
oM = float(sys.argv[5])
E2Ns = float(sys.argv[6])
rfm = float(sys.argv[7])
reps = int(sys.argv[8])
simtype = str(sys.argv[9])
path = str(sys.argv[10])

# Print the time at which simulations start running
print (datetime.date.today(), datetime.datetime.now())
startat = timeit.default_timer()
# Print the parameters used in the simulations
print ('N=%i; Ngen=%i; oF=%.7f; oM=%.7f; VA=%i; E2Ns=%i; rfm=%.2f; reps=%i; simtype=%s; path=%s'%(N,Ngen,oF,oM,VA,E2Ns,rfm,reps,simtype,path))
# Run the simulations
run_simulations_replicates (N, Ngen, oF, oM, VA, E2Ns, rfm, reps, simtype, path)
# Print the time taken to run the simulations
print ("ALL RUN TIME: %f"%(timeit.default_timer() - startat))
