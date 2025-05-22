'''
(Modified from ANA10_24/IBS_2sexes_runAna10.py)

This script is used to run various sets of simulations, used to produce results for Figure 2, Figure 4C-F and Supplementary Figures S9-S12.

Essentially, the population is burned in for 'btN' generations where sex-specific optima are both at 0. After that, a shift in sex-specific optima occurs and the population adapts. The output file reports the average and standadrd error across 'replicate' replicates of a few variables for 'generationsN'*N generations after the shift.

To run the script in the command line, type:
python script_fig2u4.py replicates N VA rfm oFi oMi generationsN E2Ns selmode btN samplehow

A usage example would be:
python script_fig2u4.py 200 1000 40 0.95 0.15 0.15 10 1 WFHW 10 full

#### INPUT

The arguments it takes correspond to:
    - replicates: number of replicates (200)
    - N: population size (1000)
    - VA: genetic variance (40)
    - rfm: intersex correlation. It takes values from 0 to 1 (rfm)
    - oFi: coefficient for female optimum after the shift. The optimum itself is then calculated as oF = oFi*np.sqrt(2*N) (oFi)
    - oMi: coefficient for male optimum after the shift. The optimum itself is then calculated as oM = oMi*np.sqrt(2*N) (oMi)
    - generationsN: after burn in and optima shift, we will keep track of the population for generationsN*N generations (10)
    - E2Ns: parameter modulating the average squared effect size of incoming mutations. It can take any value (E2Ns)
    - selmode: the type of simulations we want to run. Options are 'exact_fertility', 'exact_viability', 'WF', 'WFHW' (WFHW)
    - btN: the population will be burned-in for Ngen*N generations (10)
    - samplehow: choose whether we want to sample the population every generation or only sparsely. Options are 'full' and 'sparse' (full)

We generally use the default values indicated in brackets for each argument, and run simulations for various combinations of those parameters where no default value is indicated:
    - rfm: typically, for rfm=0.5,0.8,0.95
    - oFi, oMi: typically oFi=0.15,0.25,0.5 and (respectively) oMi=0.15,0.25,0.5 (3 scenarios of concordant adaptation with different shift sizes) and oMi=-0.15,-0.25,-0.5 (3 scenarios of dimorphic adaptation with different shift sizes)
- E2Ns: E2Ns=1 for approximately infinitesimal and E2Ns=16 for multigenic genetic architectures


#### OUTPUT

It produces an output file containing various fields, each row corresponding to a different generation. The output file contains two fields for each of the described variables, one for the average (VARIABLENAME) and another for the standard error of the mean (VARIABLENAME_se):
    - 'sampledtimes': the times (in generations) the variables have been sampled at
    - 'meanF': female trait mean
    - 'meanM': male trait mean
    - 'varF': female genetic variance
    - 'varM': male genetic variance
    - 'covar': intersex covariance
    - 'rFM': intersex correlation
    - 'm3F': 3rd central moment (female)
    - 'm3M': 3rd central moment (male)
    - 'fixedBckgF': fixed background (female), corresponding to the contributed of fixed mutations to the average phenotype
    - 'fixedBckgM': fixed background (male), corresponding to the contributed of fixed mutations to the average phenotype
    - 'varF_shared': contribution of shared mutations to female variance
    - 'varM_shared': contribution of shared mutations to male variance
    - 'varF_spec': contribution of sex-specific mutations to female variance
    - 'varM_spec': contribution of sex-specific mutations to female variance
    - 'varA_shared': contribution of shared mutations to average variance
    - 'varA_spec': contribution of sex-specific mutations to average variance
    - 'Da': distance between average phenotype to the average optimum
    - 'Dd': distance between average distance phenotype to the average distance optimum
    - 'Va': average genetic variance
    - 'Vd': average distance genetic variance
    - 'Fa': distance of the average fixed background from the average optimum
    - 'Fd': distance of the average distance fixed background from the average distance optimum
    - 'm3a': average 3rd central moment
    - 'm3d': average distance 3rd central moment
    - 'VfVm': product of sex-specific variances
    - 'epsilon': epsilon term in quasi-static approximation (Eq S.54)
    - 'Dta1': first term in quasi-static approximation for Da (Eq S.54)
    - 'Dta2': second term in quasi-static approximation for Da (Eq S.54)
    - 'Dta': quasi-stataic approximation for Da (Eq S.54)
    - 'Dtd1': first term in quasi-static approximation for Dd (Eq S.54)
    - 'Dtd2': second term in quasi-static approximation for Dd (Eq S.54)
    - 'Dtd': quasi-stataic approximation for Da (Eq S.54)
    - 'Dt2a': simplified quasi-static approximation for Da (Eq S.55)
    - 'Dt2d': simplified quasi-static approximation for Dd (Eq S.55)
    - 'intV': average across 10N generations of the average genetic variance
    - 'intVshared': average across 10N generations of the contribution of shared mutations to average genetic variance
    - 'intVspec': average across 10N generations of the contribution of sex-specific mutations to average genetic variance
    - 'intB': average across 10N generations of the intersex covariance
    - 'intrFM': average across 10N generations of the intersex correlation
    - 'intV0': average across 10N generations of the difference between average genetic variance and analytical genetic variance
    - 'intV0shared': average across 10N generations of the difference between the contribution of shared mutations to average genetic variance and analytical genetic variance
    - 'intV0spec': average across 10N generations of the difference between the contribution of sex-specific mutations to average genetic variance and analytical genetic variance
    - 'intB0': average across 10N generations of the difference between intersex covariance and analytical intersex covariance
    - 'intrFM0': average across 10N generations of the difference between intersex correlation and analytical intersex correlation
    - 'intV_half': average across 5N generations of the average genetic variance
    - 'intVshared_half': average across 5N generations of the contribution of shared mutations to average genetic variance
    - 'intVspec_half': average across 5N generations of the contribution of sex-specific mutations to average genetic variance
    - 'intB_half': average across 5N generations of the intersex covariance
    - 'intrFM_half': average across 5N generations of the intersex correlation
    - 'intV0_half': average across 5N generations of the difference between average genetic variance and analytical genetic variance
    - 'intV0shared_half': average across 5N generations of the difference between the contribution of shared mutations to average genetic variance and analytical genetic variance
    - 'intV0spec_half': average across 5N generations of the difference between the contribution of sex-specific mutations to average genetic variance and analytical genetic variance
    - 'intB0_half': average across 5N generations of the difference between intersex covariance and analytical intersex covariance
    - 'intrFM0_half': average across 5N generations of the difference between intersex correlation and analytical intersex correlation
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

def runSimulations_2sex (N, VA, q, oF, oM, generations, E2Ns, selmode, btN, samplehow):
    start_time = timeit.default_timer()    
    a0 = Analytics(N=N, U=0, optF=oM, optM=oM, q=q, E2Ns=E2Ns)
    # get mutation rate given variance
    U = a0.get_U_from_VA(VA)
    a = Analytics(N=N, U=U, optF=oF, optM=oM, q=q, E2Ns=E2Ns)
    # Run the model with the mutation rate:
    
    # Initiate the population given the chosen selection mode
    s = SimulatePopulation(N=N, U=U, optF=0, optM=0, q=q, E2Ns=E2Ns, samplehow=samplehow)
    if selmode=="exact_fertility": s.initiate_basic_population_classExact_fertilitySelection()
    elif selmode=="exact_viability": s.initiate_basic_population_classExact_viabilitySelection()
    elif selmode=="WF": s.initiate_basic_population_classWF()
    elif selmode=="WFHW": s.initiate_basic_population_classWF_HW()
    # Burn in the population for btN simulations
    s.burn_basic_population_to_reach_steady_state(burn_time_N=btN)
    # Change sex-specific optima
    s.change_optima(optF=oF, optM=oM)
    # Run simlations for 'generations' generations
    s.run_simulations(generations)
    print ("Simulations run time: %f"%(timeit.default_timer() - start_time))

    # Get analytical estimtes
    V = a.calculate_Va()
    B = a.calculate_covar()
    rFM = a.calculate_rFM()

    if samplehow=='full': #if run for all generations:
        gens = range(generations)
        l = len(gens)
        l_5N = int(len(gens)/2)
    elif samplehow=='sparse': #if run for sample times
        gens = s._used_sample_times
        l = len(gens)
        l_5N = int(len([i for i in gens if i < generations/2]))

    # Define the 'empirical integrals' (used in Figure 4C-E and Figure S12)
    # All of the fields below are a single value, which we are repeating in a list of l elements ([]*l) so that it fits in the final table
    # intV: the average across 10N generations of the average genetic variance
    intV = [sum(0.5*(np.array(s.varS_f)+np.array(s.varS_m)))/l]*l
    # intVshared: the average across 10N generations of the average genetic variance contributed by shared mutations
    intVshared = [sum(0.5*(np.array(s.varS_f_shared)+np.array(s.varS_m_shared)))/l]*l
    # intVspec: the average across 10N generations of the average genetic variance contributed by sex-specific mutations
    intVspec = [sum(0.5*(np.array(s.varS_f_spec)+np.array(s.varS_m_spec)))/l]*l
    # intB: the average across 10N generations of the intersex covariance
    intB = [sum(s.covarS)/l]*l
    # intrFM: the average across 10N generations of the intersex correlation
    intrFM = [sum(s.rFMS)/l]*l
    # intV0: the average across 10N generations of the difference between average genetic variance and the analytical genetic variance
    intV0 = [sum(0.5*(np.array(s.varS_f)+np.array(s.varS_m))-V)/l]*l
    # intV0shared: the average across 10N generations of the difference between average genetic variance contributed by shared mutations and the analytical genetic variance
    intV0shared = [sum(0.5*(np.array(s.varS_f_shared)+np.array(s.varS_m_shared))-V)/l]*l
    # intV0spec: the average across 10N generations of the difference between average genetic variance contributed by sex-specific mutations and the analytical genetic variance
    intV0spec = [sum(0.5*(np.array(s.varS_f_spec)+np.array(s.varS_m_spec))-V)/l]*l
    # intB: the average across 10N generations of the difference between intersex covariance and the analytical intersex covariance
    intB0 = [sum(s.covarS-B)/l]*l
    # intB: the average across 10N generations of the difference between intersex correlation and the analytical intersex correlation
    intrFM0 = [sum(s.rFMS-rFM)/l]*l

    # The variables below correspond to the ones defined above, except that the average is across the first 5N (instead of 10N) generations after the shift
    intV_half = [sum(0.5*(np.array(s.varS_f[:l_5N])+np.array(s.varS_m[:l_5N])))/l_5N]*l
    intVshared_half = [sum(0.5*(np.array(s.varS_f_shared[:l_5N])+np.array(s.varS_m_shared[:l_5N])))/l_5N]*l
    intVspec_half = [sum(0.5*(np.array(s.varS_f_spec[:l_5N])+np.array(s.varS_m_spec[:l_5N])))/l_5N]*l
    intB_half = [sum(s.covarS[:l_5N])/l_5N]*l
    intrFM_half = [sum(s.rFMS[:l_5N])/l_5N]*l
    intV0_half = [sum(0.5*(np.array(s.varS_f[:l_5N])+np.array(s.varS_m[:l_5N]))-V)/l_5N]*l
    intV0shared_half = [sum(0.5*(np.array(s.varS_f_shared[:l_5N])+np.array(s.varS_m_shared[:l_5N]))-V)/l_5N]*l
    intV0spec_half = [sum(0.5*(np.array(s.varS_f_spec[:l_5N])+np.array(s.varS_m_spec[:l_5N]))-V)/l_5N]*l
    intB0_half = [sum(s.covarS[:l_5N]-B)/l_5N]*l
    intrFM0_half = [sum(s.rFMS[:l_5N]-rFM)/l_5N]*l
    
    # Return all relevant variables
    df = pd.DataFrame(np.array([gens, s.meanS_f, s.meanS_m, s.varS_f, s.varS_m, s.covarS, s.rFMS, s.m3_f, s.m3_m, s.fixedBckg_f, s.fixedBckg_m,
                                s.varS_f_shared, s.varS_m_shared, s.varS_f_spec, s.varS_m_spec, 0.5*(np.array(s.varS_f_shared)+np.array(s.varS_m_shared)), 0.5*(np.array(s.varS_f_spec)+np.array(s.varS_m_spec)), 
                                s.D_a, s.D_d, s.V_a, s.V_d, s.F_a, s.F_d, s.m3_a, s.m3_d, s.VfVm,
                                s.epsilon, s.Dt_a1, s.Dt_a2, s.Dt_a, s.Dt_d1, s.Dt_d2, s.Dt_d, s.Dt2_a, s.Dt2_d,
                                intV, intVshared, intVspec, intB, intrFM, intV0, intV0shared, intV0spec, intB0, intrFM0,
                                intV_half, intVshared_half, intVspec_half, intB_half, intrFM_half, intV0_half, intV0shared_half, intV0spec_half, intB0_half, intrFM0_half]).T, 
                      columns = ['sampledtimes','meanF', 'meanM', 'varF', 'varM', 'covar', 'rFM', 'm3F', 'm3M', 'fixedBckgF', 'fixedBckgM',
                                 'varF_shared', 'varM_shared', 'varF_spec', 'varM_spec', 'varA_shared', 'varA_spec',
                                 'Da', 'Dd', 'Va', 'Vd', 'Fa', 'Fd', 'm3a', 'm3d', 'VfVm', 
                                 'epsilon', 'Dta1', 'Dta2', 'Dta', 'Dtd1', 'Dtd2', 'Dtd', 'Dt2a', 'Dt2d',
                                 'intV', 'intVshared', 'intVspec', 'intB', 'intrFM', 'intV0', 'intV0shared', 'intV0spec', 'intB0', 'intrFM0',
                                 'intV_half', 'intVshared_half', 'intVspec_half', 'intB_half', 'intrFM_half', 'intV0_half', 'intV0shared_half', 'intV0spec_half', 'intB0_half', 'intrFM0_half'])
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

'''Function to run simulations across various replicates. 
Since to output all the variables across generations for each replicate we'd need a 3D table, we are instead reporting averages and standard errors of the mean (SEM)
across replicates for each variable, and storing them in a separate file.'''
def runSimulations_2sex_reps (replicates, N, VA, q, oFi, oMi, generations, E2Ns, selmode, btN, samplehow, path):
    
    ST = timeit.default_timer()
    
    #get sex-specific optima from the coefficients
    oF = oFi * np.sqrt(2*N)
    oM = oMi * np.sqrt(2*N)
    
    #run the simulations
    dfs = []
    for i in range(replicates):
        print ('Replicate',i)
        dfs += [runSimulations_2sex (N, VA, q, oF, oM, generations, E2Ns, selmode, btN, samplehow)]
    
    #get means and SEM across replicates
    means = pd.concat(dfs).groupby(level=0).mean()
    s = pd.concat(dfs).groupby(level=0).std(ddof=1) 
    stderrs = s/np.sqrt(replicates)
    stderrs.columns=['sampledtimes_se', 'meanF_se','meanM_se','varF_se','varM_se','covar_se','rFM_se', 'm3F_se', 'm3M_se', 'fixedBckgF_se', 'fixedBckgM_se', 
                     'varF_shared_se', 'varM_shared_se', 'varF_spec_se', 'varM_spec_se', 'varA_shared_se', 'varA_spec_se',
                     'Da_se', 'Dd_se', 'Va_se', 'Vd_se', 'Fa_se', 'Fd_se', 'm3a_se', 'm3d_se', 'VfVm_se', 
                     'epsilon_se', 'Dta1_se', 'Dta2_se', 'Dta_se', 'Dtd1_se', 'Dtd2_se', 'Dtd_se', 'Dt2a_se', 'Dt2d_se',
                     'intV_se', 'intVshared_se', 'intVspec_se', 'intB_se', 'intrFM_se', 'intV0_se', 'intV0shared_se', 'intV0spec_se', 'intB0_se', 'intrFM0_se',
                     'intV_half_se', 'intVshared_half_se', 'intVspec_half_se', 'intB_half_se', 'intrFM_half_se', 'intV0_half_se', 'intV0shared_half_se', 'intV0spec_half_se', 'intB0_half_se', 'intrFM0_half_se']
    df = pd.concat([means, stderrs], axis=1)
           
    nameRes = "ana10_"+selmode+"_DF_%iN_%ibtN_%iNgen_%iReps_%.4fq_%iVA_%iE2Ns_%.4foFi_%.4foMi_%s"%(N,btN,int(generations/N),replicates,q,VA, E2Ns,oFi,oMi,samplehow)

    version = save_version(df, path, nameRes)
    
    print ("%s took this long to run:"%version)
    print (timeit.default_timer() - ST)
    
    return df


    
# RUN CODE

# Collect the simulation parameters as input arguments
replicates = int(sys.argv[1])
N = int(sys.argv[2])
VA = float(sys.argv[3])
rfm = float(sys.argv[4])
oFi = float(sys.argv[5])
oMi = float(sys.argv[6])
generationsN = int(sys.argv[7])
E2Ns = int(sys.argv[8])
selmode = str(sys.argv[9])
btN = int(sys.argv[10])
samplehow = str(sys.argv[11])
path = str(sys.argv[12])

generations = int(generationsN * N)
# We parametrize rfm in terms of q (for historical reasons in the project which are not conceptually relevant anymore)
q = (1-rfm)

# Run the simulations using the functions above and store results in file
res = runSimulations_2sex_reps(replicates, N, VA, q, oFi, oMi, generations, E2Ns, selmode, btN, samplehow, path)