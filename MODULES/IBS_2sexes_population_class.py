from collections import defaultdict, namedtuple
import math
import numpy as np
import random
from scipy.stats import gamma
from scipy.special import dawsn
from scipy.integrate import quad

from IBS_2sexes_mutation_class import Mutation

MutationalProcess = namedtuple('MutationalProcess', 'mu shape scale')

class _PopulationBasic(object):
    '''
    Class to store and evolve the population
    Parameters:
        - N: Population size
        - U: Mutational rate (per gamete per generation)
        - optF, optM: sex-specific optima
        - q: proportion of sex-specific mutations. rfm = 1-q
        - angle_s: angle determining the proportion of selection acting on each sex
        - E2Ns: Expected steady-state scaled selection coeffecient of incoming mutations.
            We draw them from a gamma distribution with this expected value
            We typically work with two regimes: approximately infinitesimal (E2Ns=1) and multigenic (E2Ns=16)
    '''
    
    def __init__(self, N, U, optF, optM, q, angle_s=math.pi/4, E2Ns=16):


        self.N = N
        self.halfN = int(N/2)
        self.Vs = 2 * self.N # width of fitness function

        # mutational process parameters
        self.E2Ns= E2Ns
        self.V2Ns= E2Ns**2 # if var=mean**2 gamma -> exponential
        self.U = U
        # these are scale and shape of the gamma distr so that mean and var are E2Ns and V2Ns:
        self._scale_s = float(self.V2Ns) / float(self.E2Ns) 
        self._shape_s = float(self.E2Ns) / float(self._scale_s)
        self._mu = MutationalProcess(mu=self.U, shape=self._shape_s, scale=self._scale_s)
        self._MU_SCALING_COEFF = float(self.Vs) / (2*float(self.N)) # HETEROZYGOTE mutation scaling parameter

        # set of segregating mutations
        self._segregating = set()
        # number of segregating mutations
        self._numSeg = 0
               
        self.q = q
        self.angle_s = angle_s #by default angle_s=pi/4
        self.gamma_sf = np.cos(self.angle_s) #for angle_s=pi/4, equals 1/sqrt(2), so that Vsf=Vsm
        self.gamma_sm = np.sin(self.angle_s)
        self.Vsf = self.Vs / (2*self.gamma_sm**2) #(1-self.gamma_sf**2))
        self.Vsm = self.Vs / (2*self.gamma_sf**2)

        self._FITNESS_OPTIMUM_F = optF
        self._FITNESS_OPTIMUM_M = optM
        self._FITNESSCOEF_F = 0.5 / float(self.Vsf) # selection model paramter
        self._FITNESSCOEF_M = 0.5 / float(self.Vsm) # selection model paramter
         
        # Moments that are required to update the segregating variants
        # 1st order central moments
        self._mean_f = 0
        self._mean_m = 0
        self._dist_ess_f = 0
        self._dist_ess_m = 0
        self._mean_fixed_ess_f = 0
        self._mean_fixed_ess_m = 0
        # 2nd order central moments
        self._var_ess_f = 0
        self._var_ess_m = 0
        self._var_ess_0 = 0 # overall variance
        self._covar_ess = 0
        # 3rd order central moments
        self._mu3_pos_ess_f = 0
        self._mu3_neg_ess_f = 0
        self._mu3_ess_f = 0
        self._mu3_pos_ess_m = 0
        self._mu3_neg_ess_m = 0
        self._mu3_ess_m = 0
        self._mu3_fff = 0
        self._mu3_ffm = 0
        self._mu3_fmm = 0
        self._mu3_mmm = 0
        # contributions of segregating (shared vs sex-spec mutations) and fixed background to mean and variance
        self._mean_f_shared = 0
        self._mean_m_shared = 0
        self._mean_f_spec = 0
        self._mean_m_spec = 0
        self._var_f_shared = 0
        self._var_m_shared = 0
        self._var_f_spec = 0
        self._var_m_spec = 0
        self._mean_fixed_ess_f_spec = self._mean_fixed_ess_f
        self._mean_fixed_ess_m_spec = self._mean_fixed_ess_m
        self._mean_fixed_ess_f_shared = 0
        self._mean_fixed_ess_m_shared = 0  

   
    def _update_essential_moments(self):
        '''
        Updates the following statistics
        '''

        mean_seg_f = 0.0
        mean_seg_m = 0.0
        var_f = 0.0
        var_m = 0.0
        var_0 = 0.0
        covar = 0.0
        mean_f_shared = 0.0
        mean_m_shared = 0.0
        mean_f_spec = 0.0
        mean_m_spec = 0.0
        var_f_shared = 0.0
        var_m_shared = 0.0
        var_f_spec = 0.0
        var_m_spec = 0.0
        mu3_pos_f = 0.0
        mu3_neg_f = 0.0
        mu3_pos_m = 0.0
        mu3_neg_m = 0.0
        fixed_state_minus_mean_fixed = 0.0
        mu3_fmm = 0.0
        mu3_ffm = 0.0
        mu3_fff = 0.0
        mu3_mmm = 0.0

        # add up contributions across mutations
        for mut in self._segregating:
            # determine sign
            sign_f=mut.derived_sign_f()
            sign_m=mut.derived_sign_m()
            # contributions to mean
            mean_seg_f += mut.segregating_effect(sex="F")
            mean_seg_m += mut.segregating_effect(sex="M")               
            # constributions to variance
            var_f += mut._var(sex="F")
            var_m += mut._var(sex="M")
            var_0 += mut._var()
            # contributions to covariance
            covar += mut._covar()

            # contributions of shared vs sex-specific mutations to mean and variance
            if round(mut.PHENO_SIZE_HET_F,7)-round(mut.PHENO_SIZE_HET_M,7)==0: #mutation is shared
                mean_f_shared += mut.segregating_effect(sex="F")
                mean_m_shared += mut.segregating_effect(sex="M")
                var_f_shared += mut._var(sex="F")
                var_m_shared += mut._var(sex="M")
            elif abs(round(mut.PHENO_SIZE_HET_F,7)-round(mut.PHENO_SIZE_HET_M,7))>0: #mutation is sex-specific
                mean_f_spec += mut.segregating_effect(sex="F")
                mean_m_spec += mut.segregating_effect(sex="M")
                var_f_spec += mut._var(sex="F")
                var_m_spec += mut._var(sex="M")                
            else: raise  Exception("Could not tell whether mutation is shared or sex-specific.")
            
            # contributions to third moment
            con_mu3_f = mut._mu3(sex="F")
            con_mu3_m = mut._mu3(sex="M")
            ## add them to the positive or negative lists
            if con_mu3_f>0: mu3_pos_f += con_mu3_f
            else: mu3_neg_f += con_mu3_f
            if con_mu3_m >0: mu3_pos_m += con_mu3_m
            else: mu3_neg_m += con_mu3_m
            mu3_fff += con_mu3_f
            mu3_mmm += con_mu3_m
            # hybrid combinations
            mu3_fmm += mut._mu3(sex='FMM')
            mu3_ffm += mut._mu3(sex='FFM')

        mean_f = mean_seg_f +self._mean_fixed_ess_f
        mean_m = mean_seg_m +self._mean_fixed_ess_m
        dist_f = self._FITNESS_OPTIMUM_F - mean_f
        dist_m = self._FITNESS_OPTIMUM_M - mean_m

        # update the statistics
        self._mean_f = mean_f
        self._mean_m = mean_m
        self._dist_ess_f = dist_f
        self._dist_ess_m = dist_m
        self._var_ess_f = var_f
        self._var_ess_m = var_m
        self._var_ess_0 = var_0
        self._covar_ess = covar
        self._mu3_pos_ess_f = mu3_pos_f
        self._mu3_neg_ess_f = mu3_neg_f
        self._mu3_ess_f = mu3_pos_f + mu3_neg_f
        self._mu3_pos_ess_m = mu3_pos_m
        self._mu3_neg_ess_m = mu3_neg_m
        self._mu3_ess_m = mu3_pos_m + mu3_neg_m
        self._mu3_fmm = mu3_fmm
        self._mu3_ffm = mu3_ffm
        self._mu3_fff = mu3_fff
        self._mu3_mmm = mu3_mmm
        #NOTE: _mu3_ess_f = _mu3_fff and _mu3_ess_m = _mu3_mmm
        
        # contributions of shared vs sex-specific mutations
        self._mean_f_shared = mean_f_shared
        self._mean_m_shared = mean_m_shared
        self._mean_f_spec = mean_f_spec
        self._mean_m_spec = mean_m_spec
        self._var_f_shared = var_f_shared
        self._var_m_shared = var_m_shared
        self._var_f_spec = var_f_spec
        self._var_m_spec = var_m_spec
        
        #number of segregating mutations
        self._numSeg = len(self._segregating)

    
    def _remove_extinct_fixed_from_seg(self):
        '''
        Calculates which mutations are extinct or fixed and removes them
        from the segregating list. Adds effect to fixed background
        '''
        # iterate over segregating mutation list and search for extinct or fixed mutations
        extinct = [mu for mu in self._segregating if mu.is_extinct()]
        fixed   = [mu for mu in self._segregating if mu.is_fixed()]

        # remove extinct and fixed mutations from the segregating list:
        for mu in (extinct + fixed):
            self._segregating.remove(mu) ## this is a python builtin function

        # for fixed mutations
        for mu in fixed:
            self._mean_fixed_ess_f += 2*mu.pheno_size_het_f()
            self._mean_fixed_ess_m += 2*mu.pheno_size_het_m()
            if round(mu.PHENO_SIZE_HET_F,7)-round(mu.PHENO_SIZE_HET_M,7)==0: #mutation is shared
                self._mean_fixed_ess_f_shared += 2*mu.pheno_size_het_f()
                self._mean_fixed_ess_m_shared += 2*mu.pheno_size_het_m()
            elif abs(round(mu.PHENO_SIZE_HET_F,7)-round(mu.PHENO_SIZE_HET_M,7))>0: #mutation is sex-specific
                self._mean_fixed_ess_f_spec += 2*mu.pheno_size_het_f()
                self._mean_fixed_ess_m_spec += 2*mu.pheno_size_het_m()
            # remove mutation from all offspring (only for Exact simulations):
            self._remove_mut(mu) # HERE: DOES THIS WORK IF _REMOVE_MUT IS DEFINED IN THE DAUGHTER CLASS?


###### MUTATION

    def _get_mutation(self, sex):
        '''
        Returns a new mutation with sex-specific effect sizes
        - sample effective scaled size (2Nea^2/Vs) from a gamma distribution
        - get the heterozygote effect size
        - get sex-specific effect sizes given angle (which also determines sign)
        
        Parameter:
        - q: which type of mutation regime:
           - random: random sex-specific effect sizes 
           - float (between 0-1): partially_sexspec. q (q=1-rfm) represents the proportion of mutations that are sex-specific. q/2 will then be female- and male-specific.
               If q=0, all effect sizes are the same in both sexes and rfm=1
        '''
        
        # the magnitude of the scaled squared effect size has gamma distribution
        scaled_size = self._get_random_selection_coefficient()

        # turn the scaled selection coefficient into the effect of a HETEROZYGOTE
        pheno_size_het  = math.sqrt(self._MU_SCALING_COEFF * scaled_size)
        

        q = self.q
        # the angle between sexes (which in turn determines the sign)
        #random angle
        if q=="random": angle_mut = math.radians(random.uniform(0,365))
        #q/2 proportion of female-/male-specific mutations
        elif 0 <= q <= 1:  
            rnd = random.random()
            if rnd < q/2: angle_mut = 0
            elif rnd < q: angle_mut = math.pi/2
            else: angle_mut = math.pi/4

        else: raise Exception ("Invalid q (proportion of sex-specific mutations): should be between 0 and 1")

        #negative sign with 0.5 probability
        if random.getrandbits(1):
            angle_mut += math.pi

        # sex-specific effect sizes
        '''
        To compute a_f and a_m one needs to take gamma_sf=cos(angle_s), gamma_sm=sin(angle_s) into account. 
        ATTENTION: This works when angle_s=pi/4 (assuming selection strength on females and males is equal), which is our default value.
        '''
        if self.angle_s == math.pi/4: 
            pheno_size_het_f = pheno_size_het * math.cos(angle_mut) / np.cos(self.angle_s)
            pheno_size_het_m = pheno_size_het * math.sin(angle_mut) / np.sin(self.angle_s)
        else: raise Exception("angle_s is not pi/4")
        
        
        # we add the mutation to the segregating list
        mu = Mutation(scaled_size=scaled_size,pheno_size_het_f=pheno_size_het_f,pheno_size_het_m=pheno_size_het_m,N=self.N, sex=sex)

        return mu

    def _get_random_selection_coefficient(self):
        '''
        Returns the squared effect size distribution
        '''
        return np.random.gamma(self._mu.shape, self._mu.scale)



class PopulationExact(_PopulationBasic):
    '''
    Class to store and evolve a population where we keep track of all individuals (individual-based simulations)
    Parameters specific to this class:
        - Selection mode: (string) determines whether selection acts on parents fitness (fertility selection)
                            or on offspring fitness (viability selection)
    
    '''
    
    def __init__(self, N, U, optF, optM, q=None, angle_s=math.pi/4, E2Ns=16, selection_mode='parents'):
        
        super(self.__class__, self).__init__(N,U, optF, optM, q, angle_s, E2Ns)
        
        # list of mutations carried by each individual of the population:
        # self._individuals_f[i], self._individuals_m[i] are dictionaries representing the genotype of individual i in the population, separately for both sexes:
        #   - The keys are mutations present in the individual
        #   - Values are the ploidity of the mutation (either 1 or 2)
        self._individuals_f = [defaultdict(int) for _ in range(self.halfN)]
        self._individuals_m = [defaultdict(int) for _ in range(self.halfN)]

        # same data structure is used as a temporary data-structure when constructing next generation
        self._offspring_f = [defaultdict(int) for _ in range(self.halfN)]
        self._offspring_m = [defaultdict(int) for _ in range(self.halfN)]
        
        # parents (mothers and fathers of females and males), to calculate rFM
        self._mothers_f = [defaultdict(int) for _ in range(self.halfN)]
        self._fathers_f = [defaultdict(int) for _ in range(self.halfN)]
        self._mothers_m = [defaultdict(int) for _ in range(self.halfN)]
        self._fathers_m = [defaultdict(int) for _ in range(self.halfN)]

        # read selection mode
        if selection_mode == 'parents': # selection acts on parents
            self._sample_offspring = self._sample_offspring_by_parents_fitness
        elif selection_mode == 'offspring': # selection acts on offspring
            self._sample_offspring = self._sample_offspring_by_children_fitness

    def next_gen(self):
        '''
        - Generates the next generation of offspring
        - Updates the essential moments
        - Updates the segregating list
        '''

        # create all offspring, using the selection mode specified above:
        self._sample_offspring()
        
        # advance generation:
        # swap individuals and offspring
        self._individuals_f, self._offspring_f = self._offspring_f, self._individuals_f
        self._individuals_m, self._offspring_m = self._offspring_m, self._individuals_m
            
        # update the stats
        #super(self.__class__, self)._update_essential_moments() # HERE: CHANGE BACK IF THIS GIVES AN ERROR
        self._update_essential_moments()
        # update the segregating list
        self._update_segregating_list()
        

    def _update_segregating_list(self):
        """
        Updates the segregating mutations list.
        """
        # Recalculate frequencies of muations
        self._reset_frequencies()
        # Remove extinct and fixed mutations from segregating list
        self._remove_extinct_fixed_from_seg()


    
    def _remove_mut(self, mut):
        """
        Removes the mutation mut from all offspring (used in _remove_extinct_fixed_from_seg()).
        """
        for indv in self._individuals_f:
            del indv[mut]
        for indv in self._individuals_m:
            del indv[mut]

    def _reset_frequencies(self):
        """
        Recalcualtes the frequencies of all the mutations in the population
        """
        # reset counts in storage frequency
        for mu in self._segregating:
            mu._frequency_class.store_freq_f = 0
            mu._frequency_class.store_freq_m = 0

        # update counts of "current freq"
        for indv in self._individuals_f:
            for mu,ploidity in indv.items():
                mu._frequency_class.store_freq_f += ploidity
        for indv in self._individuals_m:
            for mu,ploidity in indv.items():
                mu._frequency_class.store_freq_m += ploidity

        # update counts
        for mut in self._segregating:
            mut._frequency_class.update_freqs(mut._frequency_class.store_freq_f, mut._frequency_class.store_freq_m)
            
    def _phenotype(self,indv,sex=0):
        """
        Computes the phenotype of an individual
        Parameters:
            - indv: (dict) The member of the population whose phenotype is to be computed
            - sex: (string) the sex of the individual (F, M)
        Returns: (float) The phenotype of the indvidual
        """

        a2_total = 0
        
        # sum the effect of segregating mutations
        for mu,ploidity in indv.items():
            if sex == "F": a_total = mu.pheno_size_het_f()
            elif sex == "M": a_total = mu.pheno_size_het_m()
            a2_total += ploidity*a_total
            
        # add the fixed background
        if sex == "F": pheno_total = a2_total + self._mean_fixed_ess_f
        elif sex == "M": pheno_total = a2_total + self._mean_fixed_ess_m

        return pheno_total
    
    
    def _fitness_of_pheno(self, phenotype, sex=0):
        """
        Computes the fitness of a phenotype
        Parameters:
            - phenotype: (float) The phenotype of the member of the population whose fitness is to be computed
            - sex: (string) the sex of the individual (F, M)
        Returns: (float) The fitness of the indvidual
        """
        a2_total_square = 0.0
        
        if sex=="F": 
            a2_total_square += (phenotype - self._FITNESS_OPTIMUM_F) ** 2
            return math.exp(- a2_total_square * self._FITNESSCOEF_F)
        elif sex=="M": 
            a2_total_square += (phenotype - self._FITNESS_OPTIMUM_M) ** 2
            return math.exp(- a2_total_square * self._FITNESSCOEF_M)
    
            
    def _fitness_of_indv(self,indv,sex=0):
        """
        Computes the fitness of an individual
        Parameters:
            - indv: (dict) The member of the population whose fitness is to be computed
            - sex: (string) the sex of the individual (F, M)
        Returns: (float) The fitness of the indvidual
        """
        pheno_total = self._phenotype(indv,sex)
        return self._fitness_of_pheno(pheno_total,sex)


    def _sample_offspring_by_parents_fitness(self):
        """
        Sample self._offspring_f, self.offspring_m, based on parents' fitness.

        We first compute a list of the fitness of each female and male in the population
        Then normalize each to a probability vector and compute the CDF
        For every empty female (male) child slot (N/2 of them) in offspring pick random mothers and fathers based on their fitness
        Then use method _sample_child to create the female (male) child
        """
        
        # Compute the fitness
        fit_vals_f = [self._fitness_of_indv(self._individuals_f[i],"F") for i in range(self.halfN)]
        fit_vals_m = [self._fitness_of_indv(self._individuals_m[i],"M") for i in range(self.halfN)]
        # Normalize to a probability vector and compute the CDF.
        tmp_f = 1.0 / np.sum(fit_vals_f)
        cdf_f = np.add.accumulate([x*tmp_f for x in fit_vals_f])
        tmp_m = 1.0 / np.sum(fit_vals_m)
        cdf_m = np.add.accumulate([x*tmp_m for x in fit_vals_m])
        
        # generate female offspring one by one and store the parents
        self._mothers_f = []
        self._fathers_f = []
        for child in self._offspring_f:

            # sample two non-identical random parents from the cdf
            p1, p2 = np.searchsorted(cdf_f,random.random()), np.searchsorted(cdf_m,random.random())

            # generate daughter child of these two parents
            self._sample_child(child, self._individuals_f[p1], self._individuals_m[p2], 'F')
            
            # store the parents
            self._mothers_f += [self._individuals_f[p1]]
            self._fathers_f += [self._individuals_m[p2]]
        
        # generate male offspring one by one and store the parents
        self._mothers_m = []
        self._fathers_m = []
        for child in self._offspring_m:

            # sample two non-identical random parents from the cdf
            p1, p2 = np.searchsorted(cdf_f,random.random()), np.searchsorted(cdf_m,random.random())

            # generate random son of these two parents
            self._sample_child(child, self._individuals_f[p1], self._individuals_m[p2], 'M')
            
            # store the parents
            self._mothers_m += [self._individuals_f[p1]]
            self._fathers_m += [self._individuals_m[p2]]

    def _sample_offspring_by_children_fitness(self):
        """
        We use viability selection to accept or reject randomly generated offspring to sample self._offspring_f and self._offspring_m.
        - Pick 2 parents randomly
        - Then use method _sampleChild to create the child
        - Reject or accept the resulting child with probability proportional to its fitness
        Keep going till we have N little brats (half female and half male)
        """

        # generate FEMALE OFFSPRING one by one
        index = 0
        while index < self.halfN:

            # sample distinct random parents uniformly
            p1, p2 = random.randint(0,self.N/2 - 1), random.randint(0,self.N/2 - 1)

            # generate random child of these two parents
            self._sample_child(self._offspring_f[index], self._individuals_f[p1], self._individuals_m[p2], 'F')

            # reject or accept the resulting child with probability proportional to its fitness
            if random.random() < self._fitness_of_indv(self._offspring_f[index],'F'):
                index+= 1

        # generate MALE OFFSPRING one by one
        index = 0
        while index < self.halfN:

            # sample distinct random parents uniformly
            p1, p2 = random.randint(0,self.N/2 - 1), random.randint(0,self.N/2 - 1)

            # generate random child of these two parents
            self._sample_child(self._offspring_m[index], self._individuals_f[p1], self._individuals_m[p2], 'M')

            # reject or accept the resulting child with probability proportional to its fitness
            if random.random() < self._fitness_of_indv(self._offspring_m[index],'M'):
                index+= 1

                
        
   # samples one random child of parents p1 and p2
    def _sample_child(self, child, p1, p2, sex):
        """
        Takes a child and 2 parents and generates a random child with new mutations
        - The child's dictionary gets cleard
        - Each mutation of the parents is given to it with probability 0.5 if ploidity 1 or definitely if ploidity 2
        - De novo mutations are added according to a poisson distribution with mean 2*mu
        """

        # clear offspring dictionary
        child.clear()

        # for each parent
        for p in [p1, p2]:

            # for each mutation carried by the parent
            for mu,ploidity in p.items():

                # if the parent has two copies of the mutation he is bound to pass it on
                if ploidity == 2:
                    child[mu] += 1

                # if the parent is heterozygous, s/he passes it on with probabilty 0.5
                elif random.getrandbits(1):
                    child[mu] += 1
                    
        # add random de-novo mutations:
        # number of de-novo mutations is a poisson variable
        for _ in range(np.random.poisson(2*self._mu.mu)):

            # we add the mutation to the segregating list and to the new offspring (in heterozygous state)
            mu = self._get_mutation(sex)

            self._segregating.add(mu)
            child[mu] = 1
        
        


        
class PopulationWF(_PopulationBasic):
    """
    Class to store and evolve a population where we keep track of the frequencies of each mutation across individuals of a population.
    We do not keep track of every population member each generation. Instead we merely keep track of a list of segregating mutations
    in the populations updating using a Wright-Fisher process

    Parameters specific to this class:
        HW: (0 or 1). Determines whether we consider Hardy-Weinberg (HW) equilibrium:
            - If 0 (not): we do not assume HWE. We sample next generations' allele frequencies using the _wright_fisher() method
                          and update sex-specific frequencies (xf, xm)
            - If 1 (yes): we might assume HWE. We sample next generations' allele frequencies using the _wright_fisher_HW() method.
                          wether we consider HWE depends on parameter xd2
        xd2: (0 or 1): Using method _wright_fisher_HW(), determines whether we consider the squared difference between
                       sex-specific allele frequencies (WF) or set it to zero (WFHW)
            - If 0: HW is assumed (xd2=0)
            - If 1: no HW is assumed (xd2 =(xf-xm)**2)
    """

    def __init__(self, N, U, optF, optM, HW=0, xd2=0, q=None, angle_s=math.pi/4, E2Ns=16):
        super(self.__class__, self).__init__(N=N, U=U, optF=optF, optM=optM, q=q, angle_s=angle_s, E2Ns=E2Ns)
        
        self.HW = HW
        self.xd2 = xd2

    def next_gen(self):
        """
        Progresses to the next generation.
        """
        
        #Update mutation frequencies via Wright-Fisher process, potentially assuming Hardy-Weinberg equilibrium
        if self.HW==0: self._wright_fisher()
        if self.HW==1: 
            if self.xd2==0: self._wright_fisher_HW()
            elif self.xd2==1: self._wright_fisher_HW(XD2=1)
            else: raise Exception("xd2 should be 0 or 1")

        # Add de novo mutations
        self._new_mutations()

        # Remove fixed and extinct mutations from list and updates fixed effect (meanFixed):
        self._remove_extinct_fixed_from_seg()

        #update mean phenotype and variance in phenotype
        #super(self.__class__, self)._update_essential_moments() # HERE: CHANGE BACK IF THIS GIVES AN ERROR
        self._update_essential_moments()

    def _power(self, r_diff, varies, p, c, a):
        """
        The power of the fitness equation, used below. 
        a is the phenotypic effect of the mutation. c is 0,0.5 or 1. r_diff is the
        signed difference ( optimum phenotype -mean phenotype)
        """            
        return -(r_diff - a*(c-2*p))**2*float(varies)



    def _wright_fisher(self):
        """
        Updates the list of segregating mutations according to a Wright-Fisher process
        (updates sex-specific frequencies separately)
        """
        varies_f = 1.0 / (2.0*float(self.Vsf))
        varies_m = 1.0 / (2.0*float(self.Vsm))
        r_diff_f = self._dist_ess_f
        r_diff_m = self._dist_ess_m

        # for each mutation:
        for mut in self._segregating:
            # get allele frequencies
            p = mut._frequency_class.x()
            pf = mut._frequency_class.x(sex="F")
            pm = mut._frequency_class.x(sex="M")
            qf = 1-pf
            qm = 1-pm
            q = 1-p

            # get heterozygous SIGNED phenotypic effects 
            af = mut.pheno_size_het_f()
            am = mut.pheno_size_het_m()
            c_00, c_01, c_11 = 0, 1, 2

            # calculate fitness of each genotype
            power_00_f = self._power(r_diff_f, varies_f, pf, c_00, af)
            power_01_f = self._power(r_diff_f, varies_f, pf, c_01, af)
            power_11_f = self._power(r_diff_f, varies_f, pf, c_11, af)
            power_00_m = self._power(r_diff_m, varies_m, pm, c_00, am)
            power_01_m = self._power(r_diff_m, varies_m, pm, c_01, am)
            power_11_m = self._power(r_diff_m, varies_m, pm, c_11, am)

            w_00_f = math.exp(power_00_f)
            w_01_f = math.exp(power_01_f)
            w_11_f = math.exp(power_11_f)
            w_00_m = math.exp(power_00_m)
            w_01_m = math.exp(power_01_m)
            w_11_m = math.exp(power_11_m)

            # calculate average fitness
            meanfit_f = pf*pm*w_11_f + qf*qm*w_00_f + (pf*qm+qf*pm)*w_01_f
            meanfit_m = pf*pm*w_11_m + qf*qm*w_00_m + (pf*qm+qf*pm)*w_01_m

            # HERE: THIS WAS NOT COMMENTED OUT...
            #weight_f = (p**2*w_11_f + p*q*w_01_f) / meanfit_f
            #weight_m = (p**2*w_11_m + p*q*w_01_m) / meanfit_m

            # calculate expected sex-specific allele frequencies
            weight_f = (pf*pm*w_11_f + 0.5*(pf*qm+qf*pm)*w_01_f) / meanfit_f
            weight_m = (pf*pm*w_11_m + 0.5*(pf*qm+qf*pm)*w_01_m) / meanfit_m
            
            if weight_f>1:
                print('weight_f:', weight_f, ' pf:', pf,' r_diff_f:', r_diff_f, ' phenoSize:', af)
                raise Exception("The female weight is bigger than one.")
            if weight_m>1:
                print('weight_m:', weight_m, ' pm:', pm,'r_diff_m:', r_diff_m, ' phenoSize:', am)
                raise Exception("The male weight is bigger than one.")

            # update frequencies using binomial sampling
            mut._frequency_class.update_freqs(np.random.binomial(self.N, weight_f), np.random.binomial(self.N, weight_m))

            
    def _wright_fisher_HW(self, XD2=0):
        """
        Updates the list of segregating mutations according to a Wright-Fisher process, potentially assuming HWE
        Parameters:
            XD2.  If 0 (default): HW is assumed (xd2=0)
                  If 1: no HW is assumed (xd2 =(xf-xm)**2)
        """
        varies_f = 1.0 / (2.0*float(self.Vsf))
        varies_m = 1.0 / (2.0*float(self.Vsm))
        r_diff_f = self._dist_ess_f
        r_diff_m = self._dist_ess_m

        # for each mutation
        for mut in self._segregating:
            # get overall allele frequency
            p = mut._frequency_class.x()
            q = 1-p
            
            # get sex-specific SIGNED phenotypic effects
            af = mut.pheno_size_het_f()
            am = mut.pheno_size_het_m()
            c_00, c_01, c_11 = 0, 1, 2

            # calculate fitness for each genotype
            power_00_f = self._power(r_diff_f, varies_f, p, c_00, af)
            power_01_f = self._power(r_diff_f, varies_f, p, c_01, af)
            power_11_f = self._power(r_diff_f, varies_f, p, c_11, af)
            power_00_m = self._power(r_diff_m, varies_m, p, c_00, am)
            power_01_m = self._power(r_diff_m, varies_m, p, c_01, am)
            power_11_m = self._power(r_diff_m, varies_m, p, c_11, am)

            w_00_f = math.exp(power_00_f)
            w_01_f = math.exp(power_01_f)
            w_11_f = math.exp(power_11_f)
            w_00_m = math.exp(power_00_m)
            w_01_m = math.exp(power_01_m)
            w_11_m = math.exp(power_11_m)

            # calculate expected allele frequency in next generation using an equation that considers overall sex-specific frequencies
            # and the square difference between sex-specific frequencies (xd2)
            if XD2==0: xd2=0 # xd2=0 implies HW
            else: xd2 = mut._frequency_class.xd2() # does not assume HW
            
            meanfit_f = w_00_f + 2*p*(w_11_f-w_00_f)/2 + 2*p*q * (1+xd2/(4*p*q)) * (2*w_01_f-w_11_f-w_00_f)/2
            meanfit_m = w_00_m + 2*p*(w_11_m-w_00_m)/2 + 2*p*q * (1+xd2/(4*p*q)) * (2*w_01_m-w_11_m-w_00_m)/2
            
            A = ((w_11_f-w_00_f)/(2*meanfit_f) + (w_11_m-w_00_m)/(2*meanfit_m)) 
            H = (2*w_01_f-w_11_f-w_00_f)/(2*meanfit_f)+(2*w_01_m-w_11_m-w_00_m)/(2*meanfit_m)
            
            pnext = p + p*q/2 * ((1-xd2/(4*p*q))*A + (1-2*p)*(1+xd2/(4*p*q))*H)
            
            if pnext>1 or pnext<0:
                print('pnext:', pnext, ' p:', p,'r_diff_m:', r_diff_m, ' phenoSize:', am)
                raise Exception("pnext is bigger than one.")

            # obtain the new frequency using binomial sampling
            # update sex-specific allele frquencies considering that xf=xm=x
            newfreq = np.random.binomial(2*self.N, pnext)
            newfreq_f = 0.5*newfreq
            newfreq_m = 0.5*newfreq
            
            mut._frequency_class.update_freqs(newfreq_f, newfreq_m)

    
    def _new_mutations(self):
        """
        Add random de-novo mutations to segregating list.
        The number of de-novo mutations is a poisson variable with mean 2Nmu
        With two sexes: for each sex, the number of the-novo mutations is a poisson variable with mean Nmu
        """
        
        for _ in range(np.random.poisson(self.N*self.U)):
            mu = self._get_mutation(sex="F")
            self._segregating.add(mu)

        for _ in range(np.random.poisson(self.N*self.U)):
            mu = self._get_mutation(sex="M")
            self._segregating.add(mu)
        
    def _remove_mut(self, mut):
        """
        (Called in _remove_extinc_fixed_from_sex())
        Does nothing when there are no individual population members
        """
        return
            

        
        
        
        
    