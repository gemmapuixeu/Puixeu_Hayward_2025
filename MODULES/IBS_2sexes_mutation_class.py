import numpy as np
import math
import random

class Mutation(object):
    """
    Container class for mutations.

    Parameters:
        scaled_size: (float) The scaled size of the mutation
        pheno_size_f, pheno_size_m: (float) The phenotypic effect size of the mutation in both sexes
        sex: (string: f, m) The sex of the individual in which the mutation appears
        N: (int) The size of the population in which the mutation appears
    """
    
    def __init__(self, scaled_size, pheno_size_het_f, pheno_size_het_m, sex, N):
        
        # population size
        self._N = N
        # useful factor to multiply stuff with
        self._DENOMINATOR = 1.0 / (2.0 * float(self._N))
        
        self.SCALED_SIZE = scaled_size
        self.PHENO_SIZE_HET_F = pheno_size_het_f
        self.PHENO_SIZE_HET_M = pheno_size_het_m
        self.DERIVED_SIGN_F = np.sign(pheno_size_het_f) # sign of the mutation (+/-)
        self.DERIVED_SIGN_M = np.sign(pheno_size_het_m)
        self.PHENO_SIZE_HET_ABS_F = abs(pheno_size_het_f) # absolute heterozygote phenotypic effect
        self.PHENO_SIZE_HET_ABS_M = abs(pheno_size_het_m)
        
        self.Vs = 2*self._N # width of the fitness function
        
        # overall (unsexed) parameters
        self._MU_SCALING_COEFF = float(self.Vs) / (2*float(self._N))
        if self.PHENO_SIZE_HET_F == 0: self.DERIVED_SIGN = self.DERIVED_SIGN_M
        else: self.DERIVED_SIGN = self.DERIVED_SIGN_F
        self.PHENO_SIZE_HET  = math.sqrt(self._MU_SCALING_COEFF * self.SCALED_SIZE) * self.DERIVED_SIGN
        
        # the class that stores the mutation's frequency information
        self._frequency_class = _MutationFrequencies(N=N,sex=sex)
    
    def segregating_effect(self,sex=None):
        """
        Returns the contribution of the mutant to the average phenotype.
        Still needs to be multiplied with unit vector
        """
        if sex=="F": 
            x = self._frequency_class.x(sex="F", minor=False)
            effect = x * self.PHENO_SIZE_HET_F * 2
        elif sex=="M": 
            x = self._frequency_class.x(sex="M", minor=False)
            effect = x * self.PHENO_SIZE_HET_M * 2
        elif sex==None: 
            x = self._frequency_class.x()
            effect = x * self.PHENO_SIZE_HET * 2
        return effect

    def _var(self,sex=None):
        """
        Returns the current contribution of the mutant to phenotypic variance in both sexes.
        """
        if sex=="F": 
            x = self._frequency_class.x(sex="F", minor=False)
            a = self.PHENO_SIZE_HET_F
        elif sex=="M": 
            x = self._frequency_class.x(sex="M", minor=False)
            a = self.PHENO_SIZE_HET_M
        elif sex==None: 
            x = self._frequency_class.x()
            a = self.PHENO_SIZE_HET
        return 2 * (a)**2 * x * (1-x)

    def _covar(self):
        """
        Returns the current contribution of the mutant to phenotypic intersex covariance in both sexes.
        This formula assumes HW (considers only overall x). We still need to derive the more general one.
        """
        af = self.PHENO_SIZE_HET_F
        am = self.PHENO_SIZE_HET_M
        x = self._frequency_class.x()
        return 2 * af * am * x*(1-x)
    
    def _mu3(self,sex=None):#: minor=False):
        """
        Returns the current contribution of the mutant to the third moment in the phenotype distribution.
        """   
        xf = self._frequency_class.x(sex="F", minor=False)
        xm = self._frequency_class.x(sex="M", minor=False)
        x = self._frequency_class.x()
        af = self.PHENO_SIZE_HET_F
        am = self.PHENO_SIZE_HET_M
        a = self.PHENO_SIZE_HET
        
        if sex=="F": return 2 * x * (1-x) * (1 - 2*x) * af ** 3
        elif sex=="M": return 2 * x * (1-x) * (1 - 2*x) * am ** 3 
        elif sex=="FMM": return 2 * x * (1-x) * (1 - 2*x) * af*am*am 
        elif sex=="FFM": return 2 * x * (1-x) * (1 - 2*x) * af*af*am 
        elif sex==None: return 2 * x * (1-x) * (1 - 2*x) * a ** 3 

    
    def pheno_size_het_f(self):
        return self.PHENO_SIZE_HET_F

    def pheno_size_het_m(self):
        return self.PHENO_SIZE_HET_M

    def pheno_size_het(self):
        return self.PHENO_SIZE_HET
    
    def derived_sign_f(self):
        return self.DERIVED_SIGN_F
    
    def derived_sign_m(self):
        return self.DERIVED_SIGN_M

    def derived_sign(self):
        return self.DERIVED_SIGN

    def is_segregating(self):
        return self._frequency_class.is_segregating()

    def is_extinct(self, minor=False, frozen=False):
        return self._frequency_class.is_extinct(minor=minor, frozen=frozen)

    def is_fixed(self, minor=False, frozen=False):
        return self._frequency_class.is_fixed(minor=minor, frozen=frozen)

class _MutationFrequencies(object):
    '''
    Class that stores the mutation's frequency information
    '''
    
    def __init__(self, N, sex):

        # useful factor to multiply stuff with
        self._N = N
        self._DENOMINATOR = 1.0 / (2.0 * float(self._N))

        # initialize this mutation only in females or males
        self.frequency = 1
        if sex==None: self.frequency_f = -1; self.frequency_m = -1
        elif sex=="F": self.frequency_f = 1; self.frequency_m = 0
        elif sex=="M": self.frequency_f = 0; self.frequency_m = 1
            
        # to store the frequencies. Used to update them in the exact simulations in `_reset_frequencies`
        # the simulations seemed to work even without specifying this earlier...
        self.store_freq_f = 0
        self.store_freq_m = 0
        
        # How many times muts has had its frequency updated
        self.lifetime = 0
       
    def x(self, sex=None, minor=False):
        """
        Returns the fraction haplotypes with the mutation.(=frequency/2N)

        Parameters:
            sex: can be "F", "M" (for sex-specific frequencies) or False (for the average across sexes; default).
            minor: (boolean) True if we are interested in the fraction of haplos with the minor allele
        Returns:
            (float) If not minor, the fraction of the 2*_N haplotypes with the mutation.
            Otherwise, the fraction of the 2*_N haplotypes with the minor mutation.
        """
        if sex==None: frequency = self.frequency # no factor of two here, it should be divided by 2N
        elif sex=="F": frequency = self.frequency_f * 2 #need factor of two because it will be divided by 2N (and should only be divided by N)
        elif sex=="M": frequency = self.frequency_m * 2
        else: print("Sex not properly specified"); return -1
        
        if not minor: frequency = frequency
        else: frequency = self._minor_frequency(frequency) #not used in my version
        return frequency * self._DENOMINATOR

    
    def xd2(self):
        """
        Returns the squared difference between sex-specific allele frequencies
        """
        xf = self.x(sex="F")
        xm = self.x(sex="M")
        return (xf-xm)**2
        

    def _minor_frequency(self,frequency):
        """
        Returns the frequency of the minor allele
        """
        if frequency <= self._N:
            return frequency
        else:
            return 2*self._N- frequency
        
    def update_freqs(self, freq_f, freq_m, update=None):
        """
        Updates allele frequency (as counts) in females, male sand overall.
        Update: to include counts of the times every mutation has been updated (lifetime)
        """

        if self.is_segregating():

            self.frequency_f = freq_f
            self.frequency_m = freq_m
            self.frequency = freq_f+freq_m

            if update:
                self.lifetime += 1

                
    def is_segregating(self):
        """
        Tells us if the variant is not yet fixed or exctinct.
        """
        if self.frequency == 0 or self.frequency == 2 * self._N:
            return False
        else:
            return True

    def is_extinct(self, minor=False, frozen=False):
        """
        Tells us if the mutant is extinct.

        Parameters:
            minor: (boolean) True if we are interested in the minor mutant
            frozen: (boolean) True if we are interested the frozen mutant.
        Returns:
            (boolean) If not minor or frozen, returns True if the mutant is extinct.
            If minor and frozen, returns True if the mutant that was minor at freeze
            my_time is now extinct.
        """
        if minor and frozen:
            if self.frozen_freq > self._N:
                if self.frequency == 2 * self._N:
                    return True
                else:
                    return False
        if self.frequency == 0:
            return True
        else:
            return False

    def is_fixed(self, minor=False, frozen=False):
        """
        Tells us if the mutant is fixed.

        Parameters:
            minor: (boolean) True if we are interested in the minor mutant
            frozen: (boolean) True if we are interested in the frozen mutant.
        Returns:
            (boolean) If not minor or frozen, returns True if the mutant is fixed.
            If minor and frozen, returns True if the mutant that was minor at freeze
            my_time is now fixed.
        """
        if minor and frozen:
            if self.frozen_freq > self._N:
                if self.frequency == 0:
                    return True
                else:
                    return False
        if self.frequency == 2 * self._N:
            return True
        else:
            return False
