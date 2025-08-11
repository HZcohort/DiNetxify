# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 17:20:37 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

import pandas as pd
import numpy as np
from scipy.stats import binomtest
from .utility import write_log


def binomial(d1:float,d2:float,n_com:int,n_d1d2:int,n_d2d1:int,enforce:bool,log_file:str):
    """
    Perform binomial test to determine whether significantly more invidiauls (>50%) have D1/D2 diagnosed before D2/D1 among those with D1-D2 disease pair.
    If 

    Parameters
    ----------
    d1 : float, Disease 1
    d2 : float, Disease 2
    n_com : int, Number of individuals with non-temporal D1-D2
    n_d1d2 : int, Number of individuals with temporal D1->D2
    n_d2d1 : int, Number of individuals with temporal D2->D1
    enforce : bool, If true, exclude those with non-temporal D1-D2
    log_file : str, Path and prefix for the log file

    Returns
    -------
        list, binomail test results

    """
    n_com,n_d1d2,n_d2d1 = int(n_com),int(n_d1d2),int(n_d2d1)
    if enforce:
        n_sum = n_d1d2+n_d2d1
    else:
        n_sum = n_d1d2+n_d2d1+n_com
    
    if n_d1d2 == n_d2d1 == 0:
        return [d1,d2,f'{d1}->{d2}',n_com,n_d1d2,n_d2d1,np.nan,np.nan,'NA']
    elif n_d1d2 > n_d2d1:
        test_result = binomtest(n_d1d2,n_sum,alternative='greater')
        p_value = test_result.pvalue
        proportion = test_result.proportion_estimate
        proportion_low = test_result.proportion_ci().low
        proportion_high = test_result.proportion_ci().high
        write_log(log_file,f'{d1} to {d2}: done\n')
        return [d1,d2,f'{d1}->{d2}',n_com,n_d1d2,n_d2d1,p_value,proportion,f'{proportion_low:.3f}-{proportion_high:.3f}']
    else:
        test_result = binomtest(n_d2d1,n_sum,alternative='greater')
        p_value = test_result.pvalue
        proportion = test_result.proportion_estimate
        proportion_low = test_result.proportion_ci().low
        proportion_high = test_result.proportion_ci().high
        write_log(log_file,f'{d2} to {d1}: done\n')
        return [d2,d1,f'{d2}->{d1}',n_com,n_d2d1,n_d1d2,p_value,proportion,f'{proportion_low:.3f}-{proportion_high:.3f}']
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
