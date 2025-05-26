import numpy as np
import math
from scipy.stats import moment
import matplotlib.pyplot as plt
import math
from scipy.stats import gamma
from scipy.stats import expon
from scipy.special import dawsn
from scipy.integrate import quad
from sklearn.linear_model import LinearRegression

from IBS_2sexes_population_class import _PopulationBasic, PopulationExact, PopulationWF


class SimulatePopulation(object):
    '''
    Class to store one set of simulation parameters
    Parameters:
        - N: Population size
        - U: Mutational rate (per gamete per generation)
        - optF, optM: sex-specific optima
        - q: proportion of sex-specific mutations. rfm = 1-q
        - angle_s: angle determining the proportion of selection acting on each sex
                   Vsf/Vs =cos(angle_s) and Vsm/Vs=sin(angle_s)
                   By default angle_s=pi/4 so that it is equal in both sexes
        - E2Ns: Expected steady-state scaled selection coeffecient of incoming mutations.
            We draw them from a gamma distribution with this expected value
            We typically work with two regimes: approximately infinitesimal (E2Ns=1) and multigenic (E2Ns=16)
    '''

    
    def __init__(self, N, U, optF=0, optM=0, q=0, angle_s=math.pi/4, E2Ns=16):
        
        ''' Function to initialize the simulator'''
    
        # Initialize general parameters
        
        self.N = N
        self.U = U
        self.E2Ns = E2Ns
        self.Vs = 2 * self.N # with of fitness fuction
        self._BURNED = False # True if the basic population has already had a burn in to reach steady-state
    
        
        # Initialize sex-specific parameters
    
        self.opt_F = optF
        self.opt_M = optM
        self.angle_s = angle_s
        self.q = q

        # means and variances calculated across individuals (only in Exact --individual-based-- simulations)
        self.meanC_f = [];
        self.varC_f = [];
        self.meanC_m = [];
        self.varC_m = [];
        # number of segregating mutations
        self.nM = []
        # genic quantities, computed from allele frequencies and allelic effects of segregating mutations
        self.meanS_f = [];
        self.meanS_m = [];
        self.varS_f = [];
        self.varS_m = []
        self.VfVm = []
        self.meanS_T = []; # phenotypic average across sexes
        self.varS_t =[] # total genic variance calculated as the sum of the variance within and between sexes
        self.varS_0 =[] # overall additive genetic variance (variance from the total allele's magnitude)
        self.fixedBckg_f = []  # the fixed background
        self.fixedBckg_m = [] 
        self.covarS = []
        self.covarMatrixS = []
        self.D_f = []; # distance between optimum and average phenotype
        self.D_m = []
        # aspects related to third order central moments
        self.m3pos_f = [];
        self.m3neg_f = [];
        self.m3_f = [];
        self.m3pos_m = [];
        self.m3neg_m = [];
        self.m3_m = []
        self.m3C_f = [];
        self.m3C_m = []
        # contributions from shared and sex-specific mutations
        self.meanSeg_f_shared = []
        self.meanSeg_m_shared = []
        self.meanSeg_f_spec = []
        self.meanSeg_m_spec =  []
        self.fixedBckg_f_shared = []
        self.fixedBckg_m_shared = []
        self.fixedBckg_f_spec = []
        self.fixedBckg_m_spec =  []
        self.varS_f_shared =  []
        self.varS_m_shared =  []
        self.varS_f_spec =  []
        self.varS_m_spec =  []
    
        # average and average difference:
        self.opt_a = 0.5*(self.opt_F+self.opt_M)
        self.opt_d = 0.5*(self.opt_F-self.opt_M)
        self.mean_a = [];
        self.mean_d = [];
        self.D_a = [];
        self.D_d = [];
        self.F_a = [];
        self.F_d = [];
        self.V_a = [];
        self.V_d = []
        self.m3_a = []
        self.m3_d = []
        self.epsilon = []
        self.Dt_a = []
        self.Dt_d = []
        self.Dt_a1 = []
        self.Dt_a2 = []
        self.Dt_d1 = []
        self.Dt_d2 = []
        self.Dt2_a = []
        self.Dt2_d = []
        
    def initiate_basic_population_classExact_fertilitySelection(self):
        ''' Fuction to initiate the population for Exact simulations with fertility selection'''
        self._pop = PopulationExact(N=self.N, U=self.U, optF=self.opt_F, optM=self.opt_M, q=self.q, angle_s=self.angle_s, E2Ns=self.E2Ns)
        # record that we have a basic population
        self._BASIC_POPULATION = True
        
    def initiate_basic_population_classExact_viabilitySelection(self):
        ''' Fuction to initiate the population for Exact simulations with viability selection'''
        self._pop = PopulationExact(N=self.N, U=self.U, optF=self.opt_F, optM=self.opt_M, q=self.q, angle_s=self.angle_s,  E2Ns=self.E2Ns, selection_mode='offspring') 
        self._BASIC_POPULATION = True
        
    def initiate_basic_population_classWF(self):
        ''' Fuction to initiate the population for Wright-Fisher simulations'''
        self._pop = PopulationWF(N=self.N, U=self.U, optF=self.opt_F, optM=self.opt_M, q=self.q, angle_s=self.angle_s, E2Ns=self.E2Ns)
        self._BASIC_POPULATION = True

    def initiate_basic_population_classWF_HW(self):
        ''' Fuction to initiate the population for Wright-Fisher Hardy-Weinberg simulations'''
        self._pop = PopulationWF(N=self.N, U=self.U, optF=self.opt_F, optM=self.opt_M, HW=1, q=self.q, angle_s=self.angle_s, E2Ns=self.E2Ns)
        self._BASIC_POPULATION = True

    def initiate_basic_population_classWF_xd2(self):
        ''' Fuction to initiate the population for Wright-Fisher simulations parametrized in terms of the distance between sex-specific allele frequencies'''
        self._pop = PopulationWF(N=self.N, U=self.U, optF=self.opt_F, optM=self.opt_M, HW=1, xd2=1, q=self.q, angle_s=self.angle_s, E2Ns=self.E2Ns)
        self._BASIC_POPULATION = True
    
    def run_simulations(self, generations, timeAt0=1):
        
        '''
        Runs the simulations 'generations' generations and stores various variables.
        if timeAt0==1: resets all the variables
        else: appends to them (eg. after changing means)
        '''
        
        if timeAt0==1:
            self.meanC_f = [];
            self.varC_f = [];
            self.meanC_m = [];
            self.varC_m = [];
            self.meanC_T = []  
            self.varC_T = []   # total empirical variance calculated across individuals
            self.varC_t = []   # total empirical variance calculated as the sum of the variance within and between sexes
            self.varC_w = []   # variance 'within' sexes
            self.varC_b = []   # variance 'between' sexes
            self.rFMC = []
            self.covarC = []

            self.meanS_f = []
            self.varS_f = []
            self.meanS_m = []
            self.fixedBckg_f = [] # the fixed background
            self.fixedBckg_m = []
            self.varS_m = []
            self.varS_t = [] # total genic variance calculated as the sum of the variance within and between sexes
            self.varS_w = [] # variance 'within' sexes
            self.varS_b = [] # variance 'between' sexes
            self.covarS = []
            self.covarMatrixS = []
            self.rFMS = []
            self.VfVm = []

            self.meanSeg_f_shared = []
            self.meanSeg_m_shared = []
            self.meanSeg_f_spec = []
            self.meanSeg_m_spec =  []
            self.varS_f_shared =  []
            self.varS_m_shared =  []
            self.varS_f_spec =  []
            self.varS_m_spec =  []
            self.fixedBckg_f_shared = []
            self.fixedBckg_m_shared = []
            self.fixedBckg_f_spec = []
            self.fixedBckg_m_spec =  []

            self.nM = []
            self.D_f = [];
            self.D_m = []
            self.m3pos_f = [];
            self.m3neg_f = [];
            self.m3_f = [];
            self.m3pos_m = [];
            self.m3neg_m = [];
            self.m3_m = []
            self.m3C_f = [];
            self.m3C_m = []

            # average and difference:
            self.mean_a = [];
            self.mean_d = [];
            self.D_a = [];
            self.D_d = [];
            self.F_a = [];
            self.F_d = [];
            self.V_a = [];
            self.V_d = []
            self.m3_a = []
            self.m3_d = []
            self.epsilon = []
            self.Dt_a = []
            self.Dt_d = []
            self.Dt_a1 = []
            self.Dt_a2 = []
            self.Dt_d1 = []
            self.Dt_d2 = []
            self.Dt2_a = []
            self.Dt2_d = []
  
        for g in range(generations):
            '''
            # THIS PART ONLY APPLIES TO SIMULATIONS WITH PopulationExact
            '''
            if 'Exact' in str(self._pop):
                # phenotypic mean and variance
                phenos_f = [self._pop._phenotype(ind, "F") for ind in self._pop._individuals_f]
                phenos_m = [self._pop._phenotype(ind, "M") for ind in self._pop._individuals_m]
                phenos_both = phenos_f + phenos_m
                # phenotypic mean
                self.meanC_f += [np.mean(phenos_f)]
                self.meanC_m += [np.mean(phenos_m)]
                self.meanC_T += [np.mean(phenos_both)]
                # phenotypic variance
                self.varC_f += [np.var(phenos_f)]
                self.varC_m += [np.var(phenos_m)]
                self.varC_T += [np.var(phenos_both)]
                # within, between and total variance
                SSB = 0.5*self.N * ((self._pop._mean_f-0.5*(self._pop._mean_f+self._pop._mean_m))**2 +\
                                    (self._pop._mean_m-0.5*(self._pop._mean_f+self._pop._mean_m))**2)
                Vb = SSB / self.N
                Vw = 0.5*(self.varC_f[-1]+self.varC_m[-1])
                Vt = Vb + Vw
                self.varC_b += [Vb]
                self.varC_w += [Vw]
                self.varC_t += [Vt]
                # third moment (skewness)
                self.m3C_f += [moment(phenos_f, moment=3)]
                self.m3C_m += [moment(phenos_m, moment=3)]

                # Intersex correlation from parent-offspring regressions

                # calculate phenotypes of parents and offspring
                zD = [self._pop._phenotype(i, "F") for i in self._pop._individuals_f] # phenotypes of daughters
                zMD = [self._pop._phenotype(i, "F") for i in self._pop._mothers_f] # phenotypes of daughters' mothers
                zFD = [self._pop._phenotype(i, "M") for i in self._pop._fathers_f] # phenotypes of daughters' fathers
                zS = [self._pop._phenotype(i, "M") for i in self._pop._individuals_m] # phenotypes of sons
                zMS = [self._pop._phenotype(i, "F") for i in self._pop._mothers_m] # phenotypes of sons' mothers
                zFS = [self._pop._phenotype(i, "M") for i in self._pop._fathers_m] # phenotypes of sons' fathers
                # calculate rfm from parent-offspring regressions
                self.rFMC += [self.calc_rFM (zD, zMD, zFD, zS, zMS, zFS)] 
                # calculate covariance
                self.covarC += [self.rFMC[-1] * np.sqrt(self.varC_f[-1]*self.varC_m[-1])]


            ''' GENIC quantities
            calculated from allelic frequencies and effect sizes and recovered from the simulator'''
            
            self.meanS_f += [self._pop._mean_f]
            self.meanS_m += [self._pop._mean_m]
            self.meanS_T += [0.5*(self._pop._mean_f+self._pop._mean_m)] #for Nf!=Nm (self.meanS_f*Nf+self.meanS_m*Nm)/(Nf+Nm)
            self.fixedBckg_f += [self._pop._mean_fixed_ess_f]
            self.fixedBckg_m += [self._pop._mean_fixed_ess_m]
            self.varS_f += [self._pop._var_ess_f]
            self.varS_m += [self._pop._var_ess_m]
            self.VfVm += [self.varS_f[-1]*self.varS_m[-1]]
            self.varS_0 += [self._pop._var_ess_0] # Overall additive genetic variance

            # Contributions of shared and sex-specific mutations to means and variances
            self.meanSeg_f_shared += [self._pop._mean_f_shared]
            self.meanSeg_m_shared += [self._pop._mean_m_shared]
            self.meanSeg_f_spec += [self._pop._mean_f_spec]
            self.meanSeg_m_spec += [self._pop._mean_m_spec]
            self.varS_f_shared += [self._pop._var_f_shared]
            self.varS_m_shared += [self._pop._var_m_shared]
            self.varS_f_spec += [self._pop._var_f_spec]
            self.varS_m_spec += [self._pop._var_m_spec]
            self.fixedBckg_f_shared += [self._pop._mean_fixed_ess_f_shared]
            self.fixedBckg_m_shared += [self._pop._mean_fixed_ess_m_shared]
            self.fixedBckg_f_spec += [self._pop._mean_fixed_ess_f_spec]
            self.fixedBckg_m_spec += [self._pop._mean_fixed_ess_m_spec]
            
            #Calculate total variance as sum of within and between-sex variance
            SSB = 0.5*self.N * ((self._pop._mean_f-0.5*(self._pop._mean_f+self._pop._mean_m))**2 +\
                                (self._pop._mean_m-0.5*(self._pop._mean_f+self._pop._mean_m))**2)
            Vb = SSB / self.N
            Vw = 0.5*(self._pop._var_ess_f+self._pop._var_ess_m)
            Vt = Vb + Vw
            self.varS_b += [Vb]
            self.varS_w += [Vw]
            self.varS_t += [Vt]
            
            #Covariance
            self.covarS += [self._pop._covar_ess]
            self.covarMatrixS += [np.array([[self._pop._var_ess_f, self._pop._covar_ess], [self._pop._covar_ess, self._pop._var_ess_m]])]
            if self._pop._var_ess_f*self._pop._var_ess_m >0: self.rFMS += [self._pop._covar_ess/(np.sqrt(self._pop._var_ess_f)*np.sqrt(self._pop._var_ess_m))]
            else: self.rFMS += [0]
            self.D_f += [self._pop._dist_ess_f]
            self.D_m += [self._pop._dist_ess_m]
            # third moment
            self.m3pos_f += [self._pop._mu3_pos_ess_f]
            self.m3neg_f += [self._pop._mu3_neg_ess_f]
            self.m3_f += [0.5*(self._pop._mu3_fff+self._pop._mu3_fmm)] # this works too: [0.5*(self._pop._mu3_ess_f+self._pop._mu3_fmm)] because _mu3_ess_f = mu3_fff # this not: [self._pop._mu3_ess_f] updated on 300323
            self.m3pos_m += [self._pop._mu3_pos_ess_m]
            self.m3neg_m += [self._pop._mu3_neg_ess_m]
            self.m3_m += [0.5*(self._pop._mu3_mmm+self._pop._mu3_ffm)] # this works too: [0.5*(self._pop._mu3_ess_m+self._pop._mu3_ffm)] because _mu3_ess_m = mu3_mmm # this not [self._pop._mu3_ess_m] updated on 300323
            # number of segregating mutations
            self.nM += [self._pop._numSeg] #[len(self._pop._segregating)]

            # Average and average difference quantities
            self.mean_a += [0.5*(self.meanS_f[-1]+self.meanS_m[-1])]
            self.mean_d += [0.5*(self.meanS_f[-1]-self.meanS_m[-1])]
            self.D_a += [self.opt_a - self.mean_a[-1]]
            self.D_d += [self.opt_d - self.mean_d[-1]]
            self.F_a += [self.opt_a - (0.5*(self.fixedBckg_f[-1]+self.fixedBckg_m[-1]))]
            self.F_d += [self.opt_d - (0.5*(self.fixedBckg_f[-1]-self.fixedBckg_m[-1]))]
            self.V_a += [0.5*(self.varS_f[-1]+self.varS_m[-1])]
            self.V_d += [0.5*(self.varS_f[-1]-self.varS_m[-1])]
            self.m3_a += [0.25*(self._pop._mu3_fff + self._pop._mu3_fmm + self._pop._mu3_ffm + self._pop._mu3_mmm)] #[0.25*(self._pop._mu3_ess_f + self._pop._mu3_fmm + self._pop._mu3_ffm + self._pop._mu3_ess_m)]
            self.m3_d += [0.25*(self._pop._mu3_fff + self._pop._mu3_fmm - self._pop._mu3_ffm - self._pop._mu3_mmm)] #[0.25*(self._pop._mu3_ess_f + self._pop._mu3_fmm - self._pop._mu3_ffm - self._pop._mu3_ess_m)]

            # Quasi-static approximation
            if (self.V_a[-1]**2-self.covarS[-1]**2-self.V_d[-1]**2) !=0:
                self.epsilon += [self.V_d[-1]**2/(self.V_a[-1]**2-self.covarS[-1]**2-self.V_d[-1]**2)]
                self.Dt_a2 += [self.m3_d[-1]*self.V_d[-1]/(self.V_a[-1]**2-self.covarS[-1]**2-self.V_d[-1]**2)]
                self.Dt_d2 += [self.m3_a[-1]*self.V_d[-1]/(self.V_a[-1]**2-self.covarS[-1]**2-self.V_d[-1]**2)]
            else:
                self.epsilon += [np.NaN]
                self.Dt_a2 += [np.NaN]
                self.Dt_d2 += [np.NaN]
            # Without 1e-7 the division sometimes gives an error (division by zero)
            self.Dt_a1 += [self.m3_a[-1]/(self.V_a[-1]+self.covarS[-1]+1e-7)*(1-self.epsilon[-1])]
            self.Dt_d1 += [self.m3_d[-1]/(self.V_a[-1]-self.covarS[-1]+1e-7)*(1-self.epsilon[-1])]
            self.Dt_a += [self.Dt_a1[-1]-self.Dt_a2[-1]]
            self.Dt_d += [self.Dt_d1[-1]-self.Dt_d2[-1]]
            self.Dt2_a += [self.m3_a[-1]/(self.V_a[-1]+self.covarS[-1]+1e-7)]
            self.Dt2_d += [self.m3_d[-1]/(self.V_a[-1]-self.covarS[-1]+1e-7)]            

            ''' Simulate next generation'''
            self._pop.next_gen()

    def calc_rFM (self,zD, zMD, zFD, zS, zMS, zFS):
        ''' 
        Function to calculate intersex correlation empiricially (in Exact simulations) from parent-offspring regressions
        Following Bonduriansky & Rowe, 2005
        '''
        
        def calc_h2_singleparent (pheno_offspring, pheno_parent):
            # offspring are 50% similar to one parent, so regression coefficient: b=0.5*h2
            # h2=2*b of offspring to single parent
            model = LinearRegression().fit(np.array(pheno_parent).reshape((-1, 1)), np.array(pheno_offspring))
            return 2*float(model.coef_)

        h2MD = calc_h2_singleparent(zD, zMD)
        h2FD = calc_h2_singleparent(zD, zFD)
        h2MS = calc_h2_singleparent(zS, zMS)
        h2FS = calc_h2_singleparent(zS, zFS)

        if h2MD*h2FS == 0 or (h2FD*h2MS)/(h2MD*h2FS) < 0: rFM = 0
        else: rFM = np.sqrt((h2FD*h2MS)/(h2MD*h2FS))
        return rFM
                
    def burn_basic_population_to_reach_steady_state(self, burn_time_N):
        '''
        Evolves the basic pop burn_time_N*N gens to reach steady-state
        '''
        
        # if there was already a burn in period, no need to burn in
        if not self._BURNED:
            for time in range(self.N*burn_time_N):
                # advance one generation
                self._pop.next_gen()
            self._BURNED = True
    
    def change_optima(self, optF=0, optM=0, opt1sex=0):
        '''
        Changes the sex-specific optima
        '''
        # in the simulator
        self.opt_F = optF
        self.opt_M = optM
        self.opt_a = 0.5*(optF+optM)
        self.opt_d = 0.5*(optF-optM)
        # in the population
        self._pop._FITNESS_OPTIMUM_F = optF
        self._pop._FITNESS_OPTIMUM_M = optM
    
    '''
    Calculate variance average and std across the last "avgXgen" generations
    if avgXgen = "all" use all data (eg. after a burn-in)
    '''
    def get_mean_varS_f(self, avgXgen="all"): 
        if avgXgen=="all": return np.mean(self.varS_f)
        else: return np.mean(self.varS_f[avgXgen:])
    def get_mean_varS_m(self, avgXgen="all"): 
        if avgXgen=="all": return np.mean(self.varS_m)
        else: return np.mean(self.varS_m[avgXgen:])
    def get_mean_varS_t(self, avgXgen="all"): 
        if avgXgen=="all": return np.mean(self.varS_t)
        else: return np.mean(self.varS_t[avgXgen:])
    def get_sderror_varS_f(self, avgXgen="all"):
        if avgXgen=="all": np.std(self.varS_f) / np.sqrt(len(self.varS_f))
        else: return np.std(self.varS_f[avgXgen:]) / np.sqrt(avgXgen)
    def get_sderror_varS_m(self, avgXgen="all"):
        if avgXgen=="all": np.std(self.varS_m) / np.sqrt(len(self.varS_m))
        else: return np.std(self.varS_m[avgXgen:]) / np.sqrt(avgXgen)

    

    def _get_ta_td_lande(self):
        '''
        Returns the Lande time (duration of the directional selection-dominated rapid phase of adaptation)
        '''
        a = Analytics(N=self.N, U=self.U, q=self.q)
        
        Va = a.calculate_Va()
        B = a.calculate_covar()
        Vs = a.Vs
        
        k = 1

        ta = 0; td = 0;
        
        if self.opt_a!=0: ta = 2*Vs/(Va+B) * (np.log(abs(self.opt_a)) - np.log(k))
        if self.opt_d!=0 and Va!=B: td = 2*Vs/(Va-B) * (np.log(abs(self.opt_d)) - np.log(k))
        #return int(-np.log(self._DELTA_UNIT / self.shift_s0)*self.Vs / self._VAR_0)
        
        return ta,td





    
    
class Analytics(_PopulationBasic):
    '''
    Class to calculate and store various analytical calculations from the simulator. 
    See _PopulationBasic for a description of the parameters.
    '''

    def __init__(self, N, U, optF=0, optM=0, q=None, E2Ns=16):
        super(self.__class__, self).__init__(N=N, U=U, optF=optF, optM=optM, q=q, E2Ns=E2Ns)
        
    '''
    Phenotypic variance
    '''
    def calculate_VA(self):

        return 2 * self.N * self.U * quad(self._phenovar_tointegrate, 0, np.inf)[0]
    
    def _phenovar_tointegrate(self, S):
        # probability of mutation of size a
        g = gamma(a=self._shape_s, loc=0, scale=self._scale_s)
        gg = g.pdf(S)
        # marginal density of alleles with a given effect size is
        vv = 4 * math.sqrt(S) * dawsn(math.sqrt(S) / 2)
        return gg * vv
                
    def get_U_from_VA(self, VA):
        ''' 
        Calculate U for a given VA
        '''
        return VA / ( 2 * self.N * quad(self._phenovar_tointegrate, 0, np.inf)[0])
    
    def calculate_VAperU (self):
        '''
        Calculate VA corresponding to U
        '''
        return quad(self._phenovar_tointegrate, 0, np.inf)[0]

    '''
    Sex-specific variances and covariance from overall variance
    '''
    def calculate_VA_f(self):
        f = 1/np.cos(self.angle_s)**2*\
             ( (1-self.q)/2*np.cos(np.pi/4)**2 + self.q/4*np.cos(np.pi/2)**2+ self.q/4*np.cos(0)**2+
               (1-self.q)/2*np.cos(5*np.pi/4)**2 + self.q/4*np.cos(3*np.pi/2)**2+ self.q/4*np.cos(np.pi)**2)

        return f*self.calculate_VA()

    def calculate_VA_m(self):
        m = 1/np.sin(self.angle_s)**2*\
             ( (1-self.q)/2*np.sin(np.pi/4)**2 + self.q/4*np.sin(np.pi/2)**2+ self.q/4*np.sin(0)**2+
               (1-self.q)/2*np.sin(5*np.pi/4)**2 + self.q/4*np.sin(3*np.pi/2)**2+ self.q/4*np.sin(np.pi)**2)

        return m*self.calculate_VA()

    def calculate_covar(self):
        c = 1/(np.cos(self.angle_s)*np.sin(self.angle_s))*\
             ( (1-self.q)/2*np.cos(np.pi/4)*np.sin(np.pi/4) + self.q/4*np.cos(np.pi/2)*np.sin(np.pi/2) + self.q/4*np.cos(0)*np.sin(0)+
               (1-self.q)/2*np.cos(5*np.pi/4)*np.sin(5*np.pi/4) + self.q/4*np.cos(3*np.pi/2)*np.sin(3*np.pi/2) + self.q/4*np.cos(np.pi)*np.sin(np.pi))

        return c*self.calculate_VA()

    def calculate_Gmatrix(self):
        covar = self.calculate_covar()
        varF = self.calculate_VA_f()
        varM = self.calculate_VA_m()
        return np.array([[varF,covar],[covar,varM]])

    '''
    Total variance from the sum of within-sex and between-sex variance
    '''
    def calculate_VA_w(self):
        return 0.5 * (self.calculate_VA_f() + self.calculate_VA_m())
    
    def calculate_VA_b(self):
        Vb = 0.5*self.N * ((self.optF-0.5*(self.optF+self.optM))**2 +\
                             (self.optM-0.5*(self.optF+self.optM))**2)
        return Vb / self.N
    
    def calculate_VA_t(self):
        return self.calculate_VA_w() + self.calculate_VA_b()


    def calculate_rFM(self):
        '''
        Intersex correlation
        '''
        return self.calculate_covar()/np.sqrt(self.calculate_VA_f()*self.calculate_VA_m())
    

    '''
    Average and average distance variances 
    '''
    def calculate_Va(self):
        varF = self.calculate_VA_f()
        varM = self.calculate_VA_m()
        return 0.5*(varF+varM)
    
    def calculate_Vd(self):
        varF = self.calculate_VA_f()
        varM = self.calculate_VA_m()
        return 0.5*(varF-varM)
    

    def calculate_dist2(self):
        '''
        Distance from optimum
        '''
        return 1 * self.Vs / (2 * self.N)


    '''
    Calculate the third moment
    '''
    def calculate_thirdmoment(self):
        return 2 * self.N * self.U * quad(self._thirdmoment_tointegrate, 0, np.inf)[0]
        # in Mathematica: 2*1000*0.01* Integrate[2*Sqrt[S]*(1 - Exp[-S/4])*PDF[ExponentialDistribution[1/16]][S], {S, 0, \[Infinity]}]

    def _thirdmoment_tointegrate(self, S):
        # probability of mutation of size a
        g = gamma(a=self._shape_s, loc=0, scale=self._scale_s)
        gg = g.pdf(S)
        # marginal density of alleles with a given effect size is
        vv = 2 * math.sqrt(S) * (1 - np.exp(-S / 4))
        return gg * vv

    def n_segmuts (self):
        '''
        Expected number of segregating mutations
        the formula for (MAF) Soujorn time are from Hayward&Sella2022_SI
        '''
        t_xp_lower = lambda x,S: 2*self.N*x * 2*np.exp(-S*x*(1-x)) / (x*(1-x))
        t_xp_higher = lambda x,S: 2*np.exp(-S*x*(1-x)) / (x*(1-x))

        #soujorn_time = lambda x,S: t_xp_lower(x,S) if x<(1/(2*self.N)) else t_xp_higher(x,S)
        soujorn_time = lambda x,S: t_xp_lower(x,S) if x<(1/(2*self.N)) else t_xp_higher(x,S)
        soujorn_time_integ = lambda S: quad(lambda x: soujorn_time(x,S), 0,0.5)[0]

        gS = lambda S: gamma(a=self._shape_s, loc=0, scale=self._scale_s).pdf(S)

        n_segmuts = quad (lambda S: 2*self.N*self.U*gS(S)*soujorn_time_integ(S), 0, np.inf)[0]
        
        return n_segmuts
