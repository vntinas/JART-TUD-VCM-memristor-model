############################ JART-TUD VCM device model ############################
# JART-TUD VCM Library
# author: Vasileios Ntinas, Ph.D.
# affiliation: Technische Universität Dresden (TUD), Germany
# version: v0.1
# version date: June 23, 2024
# 
# Cite us
# If you use the material of this repository, you should cite Ref. [1]
# 
# Abbreviations
# JART - Jülich Aachen Resistive Switching Tools
# TUD - Technische Universität Dresden
# VCM - Valence Change Mechanism
# 
# References
# [1] V. Ntinas, D. Patel, Y. Wang, I. Messaris, V. Rana, S. Menzel, A. Ascoli, R. Tetzlaff 
#     "A Simplified Variability-Aware VCM Memristor Model for Efficient Circuit Simulation",
#     In Proc. 19th Int. Conf. on Synthesis, Modeling, Analysis and Simulation Methods and Applications to Circuit Design (SMACD'23), 2023,
#     Best Paper Award on Emerging Technologies and Applications,
#     DOI: 10.1109/SMACD58065.2023.10192107 (https://doi.org/10.1109/SMACD58065.2023.10192107)
#     BibTex: https://scholar.googleusercontent.com/scholar.bib?q=info:_LUZqGsBOhIJ:scholar.google.com/&output=citation&scisdr=ClHj7u-eEILHgMG4jJk:AFWwaeYAAAAAZnO-lJkeFrHMqTmQVnEzTEZTh4I&scisig=AFWwaeYAAAAAZnO-lL6dqc9kBy6ATnmgBWOtHx8&scisf=4&ct=citation&cd=-1
# [2] C. Bengel, A. Siemon, F. Cüppers, S. Hoffmann-Eifert, A. Hardtdegen, M. von Witzleben, L. Hellmich, R. Waser, S. Menzel,
#     "Variability-aware modeling of filamentary oxide-based bipolar resistive switching cells using SPICE level compact models." 
#     IEEE Transactions on Circuits and Systems I: Regular Papers 67.12 (2020): 4618-4630.
#     DOI: 10.1109/TCSI.2020.3018502 (https://doi.org/10.1109/TCSI.2020.3018502)
#     BibTex: https://scholar.googleusercontent.com/scholar.bib?q=info:RWKGlE_Uuf0J:scholar.google.com/&output=citation&scisdr=ClHj7u-eEILHgMHEzPM:AFWwaeYAAAAAZnPC1POLdT9gn9D4YySrGGbz7H4&scisig=AFWwaeYAAAAAZnPC1OLLASh_1P8RjeL1UOvNB8U&scisf=4&ct=citation&cd=-1
###################################################################################
import math
import numpy as np
from JART_TUD_params import *
from decimal import Decimal


#################################
# SIMPLIFIED JART-TUD VCM MODEL #
#################################
# This simplified version of the JART VCM v1b var model describes the behavior of the
# device using the following explicit formulation:
#
#           Im = f(Nd,Vm,rd,ld)                        (1)
#           dNd/dt = g(Nd,Vm,rd,ld,Ndmin,Ndmax)        (2)
# 
# where
#   Im denotes the current through the device
#   Nd denotes the state variable of the device model
#   Vm denotes the externally applied voltage across the device
#   f is an explicit formula that expresses directly Im [1]
#   g is the state equation of the device model (expressed by the ionic current of the device [2])
#   rd, ld, Ndmin, Ndmax are the variability-aware parameters of the model (capturing the device-to-device (spatial) and the cycle-to-cycle (temporal) variability of the device)
#
# - This library provides two implementations of the JART-TUD VCM device model for different programming approaches
#   1. Object-Oriented Programming approach
#       Each memristor device should be instantiated as an object of the class "JART_TUD_memristor",
#       which carries the current value of the state variable of each memristor, along with the values of the 4 variability parameters.
#       The functions 'f' and 'g' are members of the class and should be called by each memristor object individually.
#   2. Functional Programming approach
#       Only one memristor has to be instantiated just to carry all the fixed parameters of the model (this object operates as a Dictionary).
#       The functions 'f' and 'g' are NOT members of the class.
#       The non-fixed parameters (Nd, rd, ld, Ndmin, Ndmax) should be provided as arguments to the functions 'f' and 'g'.
# 
# - All the fixed parameters of the class are imported from the JART_TUD_params file

#######################################################################
# 1. Object-Oriented Programming approach
#######################################################################
class JART_TUD_memristor:
    def __init__(self,
                 Ninit=0.010,       #Default Initial Memristor State close to HRS
                 Ndiscmin=0.008,    #Default Memristor State Low Boundary at the nominal value
                 Ndiscmax=20,       #Default Memristor State High Boundary at the nominal value
                 rvar=45e-9,        #Default Conductive Filament radius at the nominal value
                 lvar=0.4           #Default Disc Region Length at the nominal value
                 ):
        self.Ndiscmin = Ndiscmin
        self.Ndiscmax = Ndiscmax
        self.Ndisc = Ninit
        self.rvar = rvar
        self.lvar = lvar
    
    def __setattr__(self, name, value):
        if name == 'Ndiscmin':
            # Check Ndiscmin value -> Hard-coded limits with values from [2]
            if value < 4e-3 or value > 25e-2:
                raise ValueError('Specified Ndiscmin out of range.')
        elif name == 'Ndiscmax':
            # Check Ndiscmax value -> Hard-coded limits with values from [2]
            if value < 18 or value > 22:
                raise ValueError('Specified Ndiscmax out of range.')
        elif name == 'Ndisc':
            # Check Ndisc value -> Changing limits due to the variability parameters from [2]
            if value < self.__dict__[f"_Ndiscmin"]*(1-1e-8) or value > self.__dict__[f"_Ndiscmax"]*(1+1e-8):    #Note: This attribute check is called inside the ODE solver.
                raise ValueError('Specified Ndisc out of range.')                                               ###### There might be a small numerical error near the boundaries.
        elif name == 'rvar':                                                                                    ###### The term 1e-8 is used to handle this kind of numerical errors.
            # Check rvar value -> Hard-coded limits with values from [2]
            if value < 40.5e-9 or value > 49.5e-9:
                raise ValueError('Specified rvar out of range.')
        elif name == 'lvar':
            # Check lvar value -> Hard-coded limits with values from [2]
            if value < 0.36 or value > 0.44:
                raise ValueError('Specified lvar out of range.')
        self.__dict__[f"_{name}"] = value
    
    def __getattr__(self, name):
        return self.__dict__[f"_{name}"]
    
    def Imem(self, V_m):
        if V_m < 0:
            if self.Ndisc <= 0: #This check helps with the convergence of the ODE solver when this function is used in simulations
                I_mem = math.nan
            else:
                d_r = (self.rvar - rdet) / (delta_r * rdet)     #Eq.(16) in [1]
                d_l = (self.lvar - ldet) / (delta_l * ldet)     #Eq.(16) in [1]
                
                # Calculate p_X,Y parameters for V_m<0 (ie. variability dependence)
                p1_0 = p1_0_n + Dp1_0_r_n*d_r + Dp1_0_l_n*d_l
                p1_1 = p1_1_n + Dp1_1_r_n*d_r + Dp1_1_l_n*d_l
                p1_2 = p1_2_n + Dp1_2_r_n*d_r + Dp1_2_l_n*d_l
                p1_3 = p1_3_n + Dp1_3_r_n*d_r + Dp1_3_l_n*d_l
                p1_4 = p1_4_n + Dp1_4_r_n*d_r + Dp1_4_l_n*d_l

                p2_0 = p2_0_n                                       # All p_X,Y terms are assigned for completeness 
                                                                    # *The assignment of non-variability-dependent p_X,Y 
                                                                    # terms can be discarded for code efficiency

                p3_0 = p3_0_n + Dp3_0_r_n*d_r + Dp3_0_l_n*d_l
                p3_1 = p3_1_n + Dp3_1_r_n*d_r + Dp3_1_l_n*d_l


                p4_0 = p4_0_n
                p4_1 = p4_1_n + Dp4_1_r_n*d_r + Dp4_1_l_n*d_l
                p4_2 = p4_2_n

                p5_0 = 0
                p5_1 = p5_1_n + Dp5_1_r_n*d_r + Dp5_1_l_n*d_l
                p5_2 = p5_2_n + Dp5_2_r_n*d_r + Dp5_2_r2_n*(d_r**2)

                p7_0 = p7_0_n + Dp7_0_r_n*d_r + Dp7_0_l_n*d_l

                p9_0 = p9_0_n + Dp9_0_r_n*d_r + Dp9_0_l_n*d_l
                p9_1 = p9_1_n
                p9_2 = p9_2_n
                p9_3 = p9_3_n + Dp9_3_r_n*d_r

                p10_0 = p10_0_n
                p10_1 = p10_1_n + Dp10_1_r_n*d_r
                p10_2 = p10_2_n
                p10_3 = p10_3_n

                p11_0 = p11_0_n + Dp11_0_r_n*d_r
                p11_1 = p11_1_n
                p11_2 = p11_2_n
                p11_3 = p11_3_n
                
                # Calculate p_X parameters for V_m<0 (ie. voltage dependence)
                p1 = p1_0*(p1_1*V_m + p1_2*(V_m**2))/(1 + p1_3*V_m + p1_4*(V_m**2))
                p2 = p2_0
                p3 = p3_0 + p3_1*V_m
                p4 = p4_0 - p4_1*np.exp(-p4_2*V_m)
                p5 = p5_0 + p5_1*V_m + p5_2*(V_m**2)
                p6 = 1
                p7 = p7_0
                p8 = 1
                p9 = p9_0 + (p9_1-p9_0)/(1+np.exp((V_m-p9_2)/p9_3))
                p10 = p10_0 + (p10_1-p10_0)/(1+np.exp((V_m-p10_2)/p10_3))
                p11 = 1/(p11_0 + (p11_1-p11_0)/(1+np.exp((V_m-p11_2)/p11_3)))
                
                # Calculate I_mem for V_m<0 (ie. state dependence)
                I_mem = p1*(p2*(np.exp((np.log(self.Ndisc/Ndmin)-p3)/p4)-1)+(np.log(self.Ndisc/Ndmin)-p3))  #ExpLin part
                I_mem = I_mem + p5/(p6 + p7*(p8*np.exp(np.log(self.Ndisc/Ndmin)-p9))**(-p10))**(1./p11)     #Generalized logistic function part
        else:
            if self.Ndisc <= 0: #This check helps with the convergence of the ODE solver when this function is used in simulations
                I_mem = math.nan
            else:
                d_r = (self.rvar - rdet) / (delta_r * rdet)     #Eq.(16) in [1]
                d_l = (self.lvar - ldet) / (delta_l * ldet)     #Eq.(16) in [1]
                
                # Calculate p_X,Y parameters for V_m>0 (ie. variability dependence)
                p5_0 = p5_1_p + Dp5_1_r_p * d_r + Dp5_1_l_p * d_l
                p5_1 = p5_1_p + Dp5_1_r_p * d_r + Dp5_1_l_p * d_l   # p5_0=p5_1 in [1]
                p5_2 = p5_2_p                                       # All p_X,Y terms are assigned for completeness 
                                                                    # *The assignment of non-variability-dependent p_X,Y 
                                                                    # terms can be discarded for code efficiency
                p6_0 = p6_0_p + Dp6_0_r_p * d_r
                p6_1 = p6_1_p
                
                p7_0 = p7_0_p + Dp7_0_r_p*d_r + Dp7_0_r2_p*(d_r**2) + Dp7_0_l_p*d_l
                p7_1 = p7_1_p
                Dp7_2_l = Dp7_2_l_p + Dp7_2_l_r_p*d_r + Dp7_2_l_r2_p*(d_r**2)     #Special case: Dp7_2_l dependency to d_r
                p7_2 = p7_2_p + Dp7_2_r_p*d_r + Dp7_2_r2_p*(d_r**2) + Dp7_2_l*d_l
                p7_3 = p7_3_p

                p8_0 = p8_0_p + Dp8_0_l_p*d_l
                p8_1 = p8_1_p

                p10_0 = p10_0_p + Dp10_0_l_p*d_l
                p10_1 = p10_1_p
                p10_2 = p10_2_p

                p11_0 = p11_0_p + Dp11_0_l_p*d_l
                p11_1 = p11_1_p
                p11_2 = p11_2_p
                
                # Calculate p_X parameters for V_m>0 (ie. voltage dependence)
                p1 = 0
                p2 = 0
                p3 = 0
                p4 = 1
                p5 = p5_0 - p5_1*np.exp(-p5_2*V_m)                  #Eq.(12) in [1]
                p6 = p6_0 + p6_1*V_m                                #First-order V_m dependency
                p7 = p7_0 + p7_1*V_m + p7_2*np.exp(-p7_3*V_m)       #Eq.(7) in [1]
                p8 = p8_0 + p8_1*V_m                                #First-order V_m dependency
                p9 = 0
                p10 = p10_0 + p10_1*V_m + p10_2*(V_m**2)            #Second-order V_m dependency
                p11 = p11_0 + p11_1*V_m + p11_2*(V_m**2)            #Second-order V_m dependency
                
                # Calculate I_mem for V_m>0 (ie. state dependence)
                I_mem = p1*(p2*(np.exp((np.log(self.Ndisc/Ndmin)-p3)/p4)-1)+(np.log(self.Ndisc/Ndmin)-p3))  #ExpLin part
                I_mem = I_mem + p5/(p6 + p7*(p8*np.exp(np.log(self.Ndisc/Ndmin)-p9))**(-p10))**(1./p11)     #Generalized logistic function part
        return I_mem
    
    
    def dNdisc_dt(self, t, y, my_pulse):
        if y < self.Ndiscmin*(1-1e-8) or y > self.Ndiscmax*(1+1e-8):
            return math.nan
        self.Ndisc = y
        V_m = my_pulse.pulse_gen(np.array([t]))
        I_mem = self.Imem(V_m)
        A = math.pi * (self.rvar ** 2)
        Rseries = RseriesTiOx + R0 * (1 + R0 * alphaline * (I_mem ** 2) * Rthline)
        Vseries = I_mem * Rseries
        Rdisc = self.lvar * 1e-9 / (self.Ndisc * 1e26 * zvo * e * un * A)
        Vdisc = I_mem * Rdisc

        if ((self.Ndisc < self.Ndiscmin) and (V_m > 0)) or ((self.Ndisc > self.Ndiscmax) and (V_m < 0)):
            Ndisc_update = 0
            return Ndisc_update
        else:
            cvo = (Nplug + self.Ndisc) / 2 * 1e26
            if V_m < 0:
                E_ion = Vdisc / (self.lvar * 1e-9)
                Rtheff = Rth0 * (rdet/self.rvar)**2
                Flim = 1 - (self.Ndisc / self.Ndiscmax) ** 10
            else:
                E_ion = (V_m - Vseries) / (lcell * 1e-9)
                Rtheff = Rth0 * Rtheff_scaling * (rdet/self.rvar)**2
                Flim = 1 - (self.Ndiscmin / self.Ndisc) ** 10
            gamma = zvo * E_ion * a / (math.pi * dWa)
            dWamin = dWa * e * (math.sqrt(1 - gamma ** 2) - gamma * math.pi / 2 + gamma * math.asin(gamma))
            dWamax = dWa * e * (math.sqrt(1 - gamma ** 2) + gamma * math.pi / 2 + gamma * math.asin(gamma))
            T = I_mem * (V_m - Vseries) * Rtheff + T0
            I_ion = zvo * e * cvo * a * ny0 * A * (
                        math.exp(-dWamin / (kb * T)) - math.exp(-dWamax / (kb * T))) * Flim
            Ndisc_update = -I_ion / (A * self.lvar * 1e-9 * e * zvo) / 1e26
            return Ndisc_update


#######################################################################
# 2. Functional Programming approach
#######################################################################
def Imem(V_m, Ndisc, rvar=45e-9, lvar=0.4):
    return np.where(V_m < 0, 
                    Imem_neg(V_m, Ndisc, rvar, lvar), 
                    Imem_pos(V_m, Ndisc, rvar, lvar))

def dNdisc_dt(V_m, Ndisc, rvar=45e-9, lvar=0.4, Ndiscmin=8e-3, Ndiscmax=20):
        if Ndisc < Ndiscmin*(1-1e-8) or Ndisc > Ndiscmax*(1+1e-8):
            return math.nan
        I_mem = Imem(V_m)
        A = math.pi * (rvar ** 2)
        Rseries = RseriesTiOx + R0 * (1 + R0 * alphaline * (I_mem ** 2) * Rthline)
        Vseries = I_mem * Rseries
        Rdisc = lvar * 1e-9 / (Ndisc * 1e26 * zvo * e * un * A)
        Vdisc = I_mem * Rdisc

        if ((Ndisc < Ndiscmin) and (V_m > 0)) or ((Ndisc > Ndiscmax) and (V_m < 0)):
            Ndisc_update = 0
            return Ndisc_update
        else:
            cvo = (Nplug + Ndisc) / 2 * 1e26
            if V_m < 0:
                E_ion = Vdisc / (lvar * 1e-9)
                Rtheff = Rth0 * (rdet/rvar)**2
                Flim = 1 - (Ndisc / Ndiscmax) ** 10
            else:
                E_ion = (V_m - Vseries) / (lcell * 1e-9)
                Rtheff = Rth0 * Rtheff_scaling * (rdet/rvar)**2
                Flim = 1 - (Ndiscmin / Ndisc) ** 10
            gamma = zvo * E_ion * a / (math.pi * dWa)
            dWamin = dWa * e * (math.sqrt(1 - gamma ** 2) - gamma * math.pi / 2 + gamma * math.asin(gamma))
            dWamax = dWa * e * (math.sqrt(1 - gamma ** 2) + gamma * math.pi / 2 + gamma * math.asin(gamma))
            T = I_mem * (V_m - Vseries) * Rtheff + T0
            I_ion = zvo * e * cvo * a * ny0 * A * (
                        math.exp(-dWamin / (kb * T)) - math.exp(-dWamax / (kb * T))) * Flim
            Ndisc_update = -I_ion / (A * lvar * 1e-9 * e * zvo) / 1e26
            return Ndisc_update


def dNdisc_dt_test(V_m, Ndisc, rvar=45e-9, lvar=0.4, Ndiscmin=8e-3, Ndiscmax=20):

    cond_nan = (Ndisc < Ndiscmin * (1 - 1e-8)) | (Ndisc > Ndiscmax * (1 + 1e-8))
    I_mem = Imem(V_m, Ndisc, rvar, lvar)
    A = np.pi * (rvar ** 2)
    Rseries = RseriesTiOx + R0 * (1 + R0 * alphaline * (I_mem ** 2) * Rthline)
    Vseries = I_mem * Rseries
    Rdisc = lvar * 1e-9 / (Ndisc * 1e26 * zvo * e * un * A)
    Vdisc = I_mem * Rdisc

    cond_update = ((Ndisc < Ndiscmin) & (V_m > 0)) | ((Ndisc > Ndiscmax) & (V_m < 0))
    Ndisc_update_0 = 0.0

    cvo = (Nplug + Ndisc) / 2 * 1e26
    E_ion = np.where(V_m < 0, Vdisc / (lvar * 1e-9), (V_m - Vseries) / (lcell * 1e-9))
    Rtheff = np.where(V_m < 0, Rth0 * (rdet/rvar)**2, Rth0 * Rtheff_scaling * (rdet/rvar)**2)
    Flim = np.where(V_m < 0, 1 - (Ndisc / Ndiscmax) ** 10, 1 - (Ndiscmin / Ndisc) ** 10)

    gamma = zvo * E_ion * a / (np.pi * dWa)
    dWamin = dWa * e * (np.sqrt(1 - gamma ** 2) - gamma * np.pi / 2 + gamma * np.arcsin(gamma))
    dWamax = dWa * e * (np.sqrt(1 - gamma ** 2) + gamma * np.pi / 2 + gamma * np.arcsin(gamma))
    T = I_mem * (V_m - Vseries) * Rtheff + T0
    I_ion = zvo * e * cvo * a * ny0 * A * (
                np.exp(-dWamin / (kb * T)) - np.exp(-dWamax / (kb * T))) * Flim
    Ndisc_update = -I_ion / (A * lvar * 1e-9 * e * zvo) / 1e26

    return np.where(cond_nan, np.nan, np.where(cond_update, Ndisc_update_0, Ndisc_update))

def Imem_neg(V_m, Ndisc, rvar=45e-9, lvar=0.4):
    d_r = (rvar - rdet) / (delta_r * rdet)     #Eq.(16) in [1]
    d_l = (lvar - ldet) / (delta_l * ldet)     #Eq.(16) in [1]
    
    # Calculate p_X,Y parameters for V_m<0 (ie. variability dependence)
    p1_0 = p1_0_n + Dp1_0_r_n*d_r + Dp1_0_l_n*d_l
    p1_1 = p1_1_n + Dp1_1_r_n*d_r + Dp1_1_l_n*d_l
    p1_2 = p1_2_n + Dp1_2_r_n*d_r + Dp1_2_l_n*d_l
    p1_3 = p1_3_n + Dp1_3_r_n*d_r + Dp1_3_l_n*d_l
    p1_4 = p1_4_n + Dp1_4_r_n*d_r + Dp1_4_l_n*d_l

    p2_0 = p2_0_n                                       # All p_X,Y terms are assigned for completeness 
                                                        # *The assignment of non-variability-dependent p_X,Y 
                                                        # terms can be discarded for code efficiency

    p3_0 = p3_0_n + Dp3_0_r_n*d_r + Dp3_0_l_n*d_l
    p3_1 = p3_1_n + Dp3_1_r_n*d_r + Dp3_1_l_n*d_l


    p4_0 = p4_0_n
    p4_1 = p4_1_n + Dp4_1_r_n*d_r + Dp4_1_l_n*d_l
    p4_2 = p4_2_n

    p5_0 = 0
    p5_1 = p5_1_n + Dp5_1_r_n*d_r + Dp5_1_l_n*d_l
    p5_2 = p5_2_n + Dp5_2_r_n*d_r + Dp5_2_r2_n*(d_r**2)

    p7_0 = p7_0_n + Dp7_0_r_n*d_r + Dp7_0_l_n*d_l

    p9_0 = p9_0_n + Dp9_0_r_n*d_r + Dp9_0_l_n*d_l
    p9_1 = p9_1_n
    p9_2 = p9_2_n
    p9_3 = p9_3_n + Dp9_3_r_n*d_r

    p10_0 = p10_0_n
    p10_1 = p10_1_n + Dp10_1_r_n*d_r
    p10_2 = p10_2_n
    p10_3 = p10_3_n

    p11_0 = p11_0_n + Dp11_0_r_n*d_r
    p11_1 = p11_1_n
    p11_2 = p11_2_n
    p11_3 = p11_3_n
    
    # Calculate p_X parameters for V_m<0 (ie. voltage dependence)
    p1 = p1_0*(p1_1*V_m + p1_2*(V_m**2))/(1 + p1_3*V_m + p1_4*(V_m**2))
    p2 = p2_0
    p3 = p3_0 + p3_1*V_m
    p4 = p4_0 - p4_1*np.exp(-p4_2*V_m)
    p5 = p5_0 + p5_1*V_m + p5_2*(V_m**2)
    p6 = 1
    p7 = p7_0
    p8 = 1
    p9 = p9_0 + (p9_1-p9_0)/(1+np.exp((V_m-p9_2)/p9_3))
    p10 = p10_0 + (p10_1-p10_0)/(1+np.exp((V_m-p10_2)/p10_3))
    p11 = 1/(p11_0 + (p11_1-p11_0)/(1+np.exp((V_m-p11_2)/p11_3)))
    
    I_mem = p1*(p2*(np.exp((np.log(Ndisc/Ndmin)-p3)/p4)-1)+(np.log(Ndisc/Ndmin)-p3))        #ExpLin part
    I_mem = I_mem + p5/(p6 + p7*(p8*np.exp(np.log(Ndisc/Ndmin)-p9))**(-p10))**(1./p11)      #Generalized logistic function part
    return np.where(V_m>0, np.nan, I_mem)


def Imem_pos(V_m, Ndisc, rvar=45e-9, lvar=0.4):
    d_r = (rvar - rdet) / (delta_r * rdet)     #Eq.(16) in [1]
    d_l = (lvar - ldet) / (delta_l * ldet)     #Eq.(16) in [1]
    
    # Calculate p_X,Y parameters for V_m>0 (ie. variability dependence)
    p5_0 = p5_1_p + Dp5_1_r_p * d_r + Dp5_1_l_p * d_l
    p5_1 = p5_1_p + Dp5_1_r_p * d_r + Dp5_1_l_p * d_l   # p5_0=p5_1 in [1]
    p5_2 = p5_2_p                                       # All p_X,Y terms are assigned for completeness 
                                                        # *The assignment of non-variability-dependent p_X,Y 
                                                        # terms can be discarded for code efficiency
    p6_0 = p6_0_p + Dp6_0_r_p * d_r
    p6_1 = p6_1_p
    
    p7_0 = p7_0_p + Dp7_0_r_p*d_r + Dp7_0_r2_p*(d_r**2) + Dp7_0_l_p*d_l
    p7_1 = p7_1_p
    Dp7_2_l = Dp7_2_l_p + Dp7_2_l_r_p*d_r + Dp7_2_l_r2_p*(d_r**2)     #Special case: Dp7_2_l dependency to d_r
    p7_2 = p7_2_p + Dp7_2_r_p*d_r + Dp7_2_r2_p*(d_r**2) + Dp7_2_l*d_l
    p7_3 = p7_3_p

    p8_0 = p8_0_p + Dp8_0_l_p*d_l
    p8_1 = p8_1_p

    p10_0 = p10_0_p + Dp10_0_l_p*d_l
    p10_1 = p10_1_p
    p10_2 = p10_2_p

    p11_0 = p11_0_p + Dp11_0_l_p*d_l
    p11_1 = p11_1_p
    p11_2 = p11_2_p
    
    # Calculate p_X parameters for V_m>0 (ie. voltage dependence)
    p1 = 0
    p2 = 0
    p3 = 0
    p4 = 1
    p5 = p5_0 - p5_1*np.exp(-p5_2*V_m)                  #Eq.(12) in [1]
    p6 = p6_0 + p6_1*V_m                                #First-order V_m dependency
    p7 = p7_0 + p7_1*V_m + p7_2*np.exp(-p7_3*V_m)       #Eq.(7) in [1]
    p8 = p8_0 + p8_1*V_m                                #First-order V_m dependency
    p9 = 0
    p10 = p10_0 + p10_1*V_m + p10_2*(V_m**2)            #Second-order V_m dependency
    p11 = p11_0 + p11_1*V_m + p11_2*(V_m**2)            #Second-order V_m dependency
    
    I_mem = p1*(p2*(np.exp((np.log(Ndisc/Ndmin)-p3)/p4)-1)+(np.log(Ndisc/Ndmin)-p3))        #ExpLin part
    I_mem = I_mem + p5/(p6 + p7*(p8*np.exp(np.log(Ndisc/Ndmin)-p9))**(-p10))**(1./p11)      #Generalized logistic function part
    I_mem = np.where(Ndisc <= 0, np.nan, I_mem)
    return np.where(V_m<0, np.nan, I_mem)


##########################
# PULSE GENERATION CLASS #
##########################
class pulse:
    
    def __init__(self,
        ## Pulse Parameters ##
        p_form,
        # DC Pulse #
        V_dc=math.nan,
        # Triangular Pulse #
        V_tr_p=math.nan,
        t_tr_p=math.nan,
        V_tr_n=math.nan,
        t_tr_n=math.nan,
        # Square Pulse #
        V_sq_low=math.nan,
        V_sq_high=math.nan,
        t_sq_d=math.nan,
        t_sq_r=math.nan,
        t_sq_width=math.nan,
        t_sq_f=math.nan,
        t_sq_wait=math.nan,
        # Sinusoidal Pulse #
        V_sin_A=math.nan,
        t_sin_per=math.nan,
        phi_0=math.nan,
        V_sin_offset=math.nan):
        self.__p_form = p_form
        if p_form=='DC':
            if math.isnan(V_dc):
                raise NameError("\nDC input is selected but the amplitude is missing.\nPlease insert the V_dc argument!")
            else:
                self.__V_dc = V_dc
        elif p_form=="trig":
            if math.isnan(V_tr_p):
                raise NameError("\nTriangular input is selected but the positive amplitude is missing.\nPlease insert the V_tr_p argument!")
            elif math.isnan(t_tr_p):
                raise NameError("\nTriangular input is selected but the positive duration is missing.\nPlease insert the t_tr_p argument!")
            elif math.isnan(V_tr_n):
                raise NameError("\nTriangular input is selected but the negative amplitude is missing.\nPlease insert the V_tr_n argument!")
            elif math.isnan(t_tr_n):
                raise NameError("\nTriangular input is selected but the negative duration is missing.\nPlease insert the t_tr_n argument!")
            else:
                self.__V_tr_p = V_tr_p
                self.__t_tr_p = t_tr_p
                self.__V_tr_n = V_tr_n
                self.__t_tr_n = t_tr_n
        elif p_form=="square":
            if math.isnan(V_sq_low):
                raise NameError("\nSquare input is selected but the lowest amplitude is missing.\nPlease insert the V_sq_low argument!")
            elif math.isnan(V_sq_high):
                raise NameError("\nSquare input is selected but the highest amplitude is missing.\nPlease insert the V_sq_high argument!")
            elif math.isnan(t_sq_d):
                raise NameError("\nSquare input is selected but the delay duration is missing.\nPlease insert the t_sq_d argument!")
            elif math.isnan(t_sq_r):
                raise NameError("\nSquare input is selected but the rising edge duration is missing.\nPlease insert the t_sq_r argument!")
            elif math.isnan(t_sq_width):
                raise NameError("\nSquare input is selected but the pulse duration is missing.\nPlease insert the t_sq_width argument!")
            elif math.isnan(t_sq_f):
                raise NameError("\nSquare input is selected but the falling edge duration is missing.\nPlease insert the t_sq_f argument!")
            elif math.isnan(t_sq_wait):
                raise NameError("\nSquare input is selected but the after-pulse wait duration is missing.\nPlease insert the t_sq_wait argument!")
            else:
                self.__V_sq_low = V_sq_low
                self.__V_sq_high = V_sq_high
                self.__t_sq_d = t_sq_d
                self.__t_sq_r = t_sq_r
                self.__t_sq_width = t_sq_width
                self.__t_sq_f = t_sq_f
                self.__t_sq_wait = t_sq_wait
        elif p_form=="sin":
            if math.isnan(V_sin_A):
                raise NameError("\nSinusoidal input is selected but the amplitude is missing.\nPlease insert the V_sin_A argument!")
            elif math.isnan(t_sin_per):
                raise NameError("\nSinusoidal input is selected but the period is missing.\nPlease insert the t_sin_per argument!")
            elif math.isnan(phi_0):
                raise NameError("\nSinusoidal input is selected but the initial phase is missing.\nPlease insert the phi_0 argument!")
            elif math.isnan(V_sin_offset):
                raise NameError("\nSinusoidal input is selected but the voltage offset is missing.\nPlease insert the V_sin_offset argument!")
            else:
                self.__V_sin_A = V_sin_A
                self.__t_sin_per = t_sin_per
                self.__phi_0 = phi_0
                self.__V_sin_offset = V_sin_offset
        else:
            raise NameError("Invalid pulse form. Please insert a valid p_form argument (DC / trig / square / sin)!")

    def pulse_gen(self, t):
        V_ms = np.zeros(len(t))
        if self.__p_form == 'DC':
            for i in range(len(t)):
                V_ms[i] = self.__V_dc
        elif self.__p_form == 'trig':
            for i in range(len(t)):
                t_mod = float(Decimal(t[i]) % Decimal(2*self.__t_tr_p+2*self.__t_tr_n))
                if t_mod<=self.__t_tr_p:
                    V_m = self.__V_tr_p * (t_mod / self.__t_tr_p)
                elif t_mod<=2*self.__t_tr_p:
                    V_m = 2*self.__V_tr_p - self.__V_tr_p * (t_mod / self.__t_tr_p)
                elif t_mod<=2*self.__t_tr_p+self.__t_tr_n:
                    V_m = self.__V_tr_n * ((t_mod-2*self.__t_tr_p) / self.__t_tr_n)
                else:
                    V_m = 2*self.__V_tr_n - self.__V_tr_n * ((t_mod-2*self.__t_tr_p) / self.__t_tr_n)
                V_ms[i] = V_m
        elif self.__p_form == 'square':
            for i in range(len(t)):
                t_mod = float(Decimal(t[i]) % Decimal(self.__t_sq_d+self.__t_sq_r+self.__t_sq_width+self.__t_sq_f+self.__t_sq_wait))
                if t_mod<=self.__t_sq_d:
                    V_m = self.__V_sq_low
                elif t_mod<=self.__t_sq_d+self.__t_sq_r:
                    V_m = self.__V_sq_low + (self.__V_sq_high-self.__V_sq_low) * (t_mod-self.__t_sq_d) / self.__t_sq_r
                elif t_mod<=self.__t_sq_d+self.__t_sq_r+self.__t_sq_width:
                    V_m = self.__V_sq_high
                elif t_mod<=self.__t_sq_d+self.__t_sq_r+self.__t_sq_width+self.__t_sq_f:
                    V_m = self.__V_sq_high - (self.__V_sq_high-self.__V_sq_low) * (t_mod-self.__t_sq_d-self.__t_sq_r-self.__t_sq_width) / self.__t_sq_f
                else:
                    V_m = self.__V_sq_low
                V_ms[i] = V_m
        elif self.__p_form=='sin':
            for i in range(len(t)):
                V_ms[i] = self.__V_sin_A*math.sin(2*math.pi*t[i]/self.__t_sin_per+self.__phi_0) + self.__V_sin_offset
        return V_ms