'''
(Modified from ANA8_24_2/IBS_2sexes_runAna8.py)


This script is used to run various sets of simulations, used to produce results for Figure 3 and Figure 4A-B.

Essentially, the population is burned in for 3/2*'generations' generations where sex-specific optima are at 0 (set by oF0, oM0). After that, a shift in sex-specific optima (to oF1, oM1) occurs and the population adapts for 'generations' generations. Then, sex-specific optima are switched back to 0 (or to values oF2, oM2) and the population adapts for another 'generations' generations. The output file reports the average and standadrd error across 'replicates' replicates of a few variables across this process (excluding the first 'generations' generations of the initial burn-in process).
Adaptation follows a 2-sex Wright-Fisher Hardy-Weinberg selection process.

To run the script in the command line, type:
python script_fig3u4.py replicates N VA rfm oi generationsN E2Ns path

A usage example would be:
module load python
python script_fig3u4.py 200 1000 40 0.95 0.15 20 1 full_path_to_ALREADY_CREATED_results_directory/

#### INPUT

The arguments it takes correspond to:
    - replicates: number of replicates (200)
    - N: population size (1000)
    - VA: genetic variance (40)
    - rfm: intersex correlation. It takes values from 0 to 1 (rfm)
    - oFi: coefficient for female optimum after the shift. The optimum itself is then calculated as oF = oFi*np.sqrt(2*N) (0.15)
    - oMi: coefficient for male optimum after the shift. The optimum itself is then calculated as oM = oMi*np.sqrt(2*N) (0.15)
    - generationsN: number of generations scaled by population size. The number of generations corresponds to generations = generationsN*N (20)
    - E2Ns: parameter modulating the average squared effect size of incoming mutations. It can take any value (E2Ns)
    - path: the path where the results will be stored. The directory should be created before running the simulations. It should be of the type full_path_to_results_directory/

We use the default values indicated in brackets for each argument, and run simulations for various combinations of those parameters where no default value is indicated:
    - rfm: typically, for rfm=0.5,0.8,0.95
    - E2Ns: E2Ns=1 for approximately infinitesimal and E2Ns=16 for multigenic genetic architectures


#### OUTPUT

It produces an output file containing various fields, each row corresponding to a different generation. The output file contains two fields for each of the described variables, one for the average (VARIABLENAME) and another for the standard error of the mean (VARIABLENAMEstderr):
    - 'meanF': female trait mean
    - 'meanM': male trait mean
    - 'varF': female genetic variance
    - 'varM': male genetic variance
    - 'covar': intersex covariance
    - 'rFM': intersex correlation
    - 'varAshared': contribution of shared mutations to average variance
    - 'varAspec': contribution of sex-specific mutations to average variance

'''

# IMPORT

# Import the relevant modules
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import timeit
import random
import os

# We also need to import our modules, which include all the functions used in our simulations
import sys
sys.path.append('../MODULES/')
from IBS_2sexes_simul_class import SimulatePopulation, Analytics


# FUNCTIONS

''' Function to run the simulations for a single replicate'''

def run_simulations_divergent_convergent (generations, q, N, VA, E2Ns, oF0, oM0, oF1, oM1, oF2, oM2):
    start_time = timeit.default_timer()
    # initialize analytics with scalemut
    a0 = Analytics(N=N, U=0, optF=oF0, optM=oM0, q=q, E2Ns=E2Ns)
    # get mutation rate given variance
    U = a0.get_U_from_VA(VA)
    # run the model with the mutation rate
    # initiate population as WF HW
    # with optima oF0, oM0
    s = SimulatePopulation(N=N, U=U, optF=oF0, optM=oM0, q=q, E2Ns=E2Ns)
    s.initiate_basic_population_classWF_HW()
    # burn population in (without keeping track of it)
    s.burn_basic_population_to_reach_steady_state(burn_time_N=int(generations/N))
    # run simulations for generations/2 generations before optima shift (recording output)
    s.run_simulations(int(generations/2))
    # shift sex-specific optima to oF1, oM1
    s.change_optima(optF=oF1, optM=oM1)
    # run simulations for generations generations
    s.run_simulations(generations, timeAt0=0)
    # shift sex-specific optima to oF2, oM2 
    # run simulations for generations generations
    s.change_optima(optF=oM2, optM=oM2)
    s.run_simulations(generations, timeAt0=0)
    print ("WF-HW simulations run time: %f"%(timeit.default_timer() - start_time))
    
    # store relevant variables
    df = pd.DataFrame(np.array([s.meanS_f, s.meanS_m, s.varS_f, s.varS_m, s.covarS, s.rFMS, 0.5*(np.array(s.varS_f_shared)+np.array(s.varS_m_shared)), 0.5*(np.array(s.varS_f_spec)+np.array(s.varS_m_spec))]).T,
                      columns = ['meanF', 'meanM', 'varF', 'varM', 'covar', 'rFM', 'varAshared', 'varAspec'])
    return df


'''Function to run simulations across various replicates. 
Since to output all the variables across generations for each replicate we'd need a 3D table, we are instead reporting averages and standard errors of the mean (SEM)
across replicates for each variable, and storing them in a separate file.'''

def run_simulations_divergent_convergent_reps (replicates, generations, q, N, VA, E2Ns, oF0, oM0, oF1, oM1, oF2, oM2, path):
    
    ST = timeit.default_timer()
    
    #run the simulations
    dfs = []
    for i in range(replicates):
        dfs += [run_simulations_divergent_convergent (generations, q, N, VA, E2Ns, oF0, oM0, oF1, oM1, oF2, oM2)]
    
    #get means and standard deviations across replicates
    means = pd.concat(dfs).groupby(level=0).mean()
    s = pd.concat(dfs).groupby(level=0).std(ddof=1) 
    # with ddof=1 we are calculating s; with ddof=0 we are calculating sigma. We want s since we don't know the standard deviation of the population
    # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.std.html
    stderrs = s/np.sqrt(replicates)
    stderrs.columns=['meanFstderr','meanMstderr','varFstderr','varMstderr','covarstderr','rFMstderr', 'varAsharedstderr', 'varAspecstderr']
    df = pd.concat([means, stderrs], axis=1)
    
    #path = "RESULTS/"
    name = "ana8inf_%iN_%iNgen_%iReps_%.3fq_%iVA_%iE2Ns"%(N,int(generations/N),replicates,q,VA,E2Ns)
    
    version = save_version(df, path, name)
    
    print ("%s took this long to run:"%version)
    print (timeit.default_timer() - ST)
    
    return df

''' Function to save the simulation results into a file
If the given name already exists (meaning, we have already run simulations with the same set of parameters), it creates a new file
(with a new version) to ensure we are not removing previous results.'''

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

# RUN CODE

# Collect the simulation parameters as input arguments
replicates = int(sys.argv[1])
N = int(sys.argv[2])
VA = int(sys.argv[3])
rfm = float(sys.argv[4])
oFMi = float(sys.argv[5])
generationsN = float(sys.argv[6])
E2Ns = int(sys.argv[7])
path = str(sys.argv[8])

generations = int(generationsN * N)

# Set optima. oF0,oM0,oF2,oM2 are set to 0 by default but can be changed
oF0=0
oM0=0
oF2=0
oM2=0
oF1 = oFMi * np.sqrt(2*N)
oM1 = -oFMi * np.sqrt(2*N)
# We parametrize rfm in terms of q (for historical reasons in the project which are not conceptually relevant anymore)
q = (1-rfm)

# Print the time at which simulations start running
#print (datetime.date.today(), datetime.datetime.now())
startat = timeit.default_timer()
# Run the simulations using the functions above and store results in file
df = run_simulations_divergent_convergent_reps (replicates, generations, q, N, VA, E2Ns, oF0, oM0, oF1, oM1, oF2, oM2, path)
# Print the time taken to run the simulations
print ("ALL RUN TIME: %f"%(timeit.default_timer() - startat))


