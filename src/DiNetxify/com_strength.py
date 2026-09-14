# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 02:25:20 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

import pandas as pd
import numpy as np
from scipy.stats import t, norm
#from .data_management import DiseaseNetworkData
from .utility import write_log

def com_rr(n:int, c:int, p1:int, p2:int, method:str):
    """
    Calculate relative risk (RR; may also be referred to as the
    observed-to-expected ratio, OER) and its two-sided significance test.

    RR/OER is calculated as:

        RR = (n * c) / (p1 * p2)

    where p1*p2/n is the expected number of individuals with both diseases
    under independence.

    Parameters
    ----------
    n : int
        Total number of individuals.
    c : int
        Number of individuals with temporal/non-temporal d1 and d2 disease pair.
    p1 : int
        Number of individuals with d1 diagnosis.
    p2 : int
        Number of individuals with d2 diagnosis.
    method : {'multinomial', 'original'}
        Method used to calculate the standard error and P-value.

        'multinomial':
            Uses the multinomial delta method to estimate the standard
            error of log(RR/OER).

        'original':
            Uses the original standard-error expression,
            with a Student's t distribution with n degrees of freedom.

    Returns
    -------
    rr : float
        Relative risk / observed-to-expected ratio.
    theta : float
        Standard error of log(RR/OER).
    p : float
        Two-sided P-value.
    """

    rr = (n*c) / (p1*p2)
    if method == 'multinomial':
        # Multinomial delta-method SE of log(RR/OER)
        theta = (1/c - 1/p1 - 1/p2 - 1/n + 2*c/(p1*p2))**0.5
        z_ = abs(np.log(rr) / theta)
        # Two-sided Wald z-test
        p = norm.sf(z_) * 2
    elif method == 'original':
        # Original method
        theta = (1/c + 1/((p1*p2)/n) - 1/n - 1/n)**0.5
        t_ = abs(np.log(rr) / theta)
        p = t.sf(t_, n) * 2
    else:
        raise ValueError("method must be 'multinomial' or 'original'")
    return rr, theta, p

def com_phi(n:int, c:int, p1:int, p2:int, method:str):
    """
    Calculate phi-correlation and its two-sided significance test.

    Parameters
    ----------
    n : int
        Total number of individuals.
    c : int
        Number of individuals with both d1 and d2 diagnosis.
    p1 : int
        Number of individuals with d1 diagnosis.
    p2 : int
        Number of individuals with d2 diagnosis.
    method : {'multinomial', 'original'}, default='multinomial'
        Method used to calculate the standard error and P-value.

        'multinomial':
            Uses the multinomial delta method to estimate the standard
            error of phi and a two-sided Wald z-test.

        'original':
            Uses the original Pearson-correlation-based approach.

    Returns
    -------
    phi : float
        Phi correlation coefficient.
    phi_theta : float
        Standard error of phi.
    p_phi : float
        Two-sided P-value.
    """

    try:
        phi = (c*n - p1*p2) / (((p1*p2)*(n-p1)*(n-p2))**0.5)
    except:
        raise ValueError('phi correlation calculation error, either number of '
                         'individuals with d1/d2 diagnosis is zero or equal to '
                         'total number of individuals.')
    if method == 'multinomial':
        # 2 x 2 cell probabilities
        p11 = c / n
        p10 = (p1 - c) / n
        p01 = (p2 - c) / n
        p00 = (n - p1 - p2 + c) / n
        probs = np.array([p00, p10, p01, p11])
        # Marginal probabilities
        a = p10 + p11
        b = p01 + p11
        # Denominator of phi
        d = np.sqrt(a * (1-a) * b * (1-b))
        # Gradient of phi
        da = np.array([0., 1., 0., 1.])
        db = np.array([0., 0., 1., 1.])
        du = np.array([0., -b, -a, 1-a-b])
        dlogd_da = (1 - 2*a) / (2*a*(1-a))
        dlogd_db = (1 - 2*b) / (2*b*(1-b))
        grad = (du / d - phi * (dlogd_da * da + dlogd_db * db))
        # Multinomial covariance matrix
        sigma = np.diag(probs) - np.outer(probs, probs)
        # Delta-method standard error
        phi_theta = np.sqrt(grad @ sigma @ grad / n)
        # Wald z-test
        if phi_theta == 0:
            phi_t = np.inf
        else:
            phi_t = abs(phi / phi_theta)
        p_phi = norm.sf(phi_t) * 2
    elif method == 'original':
        # Pearson correlation
        if abs(phi) == 1:
            phi_theta = 0.0
            phi_t = np.inf
        else:
            phi_theta = ((1 - phi**2) / (n - 2))**0.5
            phi_t = abs(phi / phi_theta)
        p_phi = t.sf(phi_t, n - 2) * 2
    else:
        raise ValueError("method must be 'multinomial' or 'original'")
    return phi, phi_theta, p_phi
    

def com_phi_rr(args) -> list:
    """
    Estimate comorbidity strength for a disease pair using phi-correlation and RR.

    Parameters:
    ----------
    d1 : float
        Disease 1

    d2 : float
        Disease 2

    message : string
        additional comment

    Global Variables:
    ----------
    trajectory : dictionary
    DiseaseNetworkData.trajectory dictionary

    threshold_config : tuple
        ``(proportion_threshold, n_threshold)`` configuration. The
        proportional threshold is calculated from this disease pair's
        eligible sub-cohort after history and sex restrictions; the absolute
        threshold is used unchanged.
    
    log_file : str
        Path and prefix for the log file
    Returns:
    ----------
    result : list
        list, comorbidity strength estimation results
    """
    # shared global data
    global trajectory_
    global threshold_config_
    global log_file_
    global se_method_

    d1, d2, message = args
    ineligible_d_dict = trajectory_['ineligible_disease']
    eligible_d_dict_withdate = trajectory_['eligible_disease_withdate']
    temporal_pair_dict = trajectory_['d1d2_temporal_pair']
    com_pair_dict = trajectory_['d1d2_com_pair']
    disease_pair_index = trajectory_['disease_pair_index']
    d1d2_index = disease_pair_index[f'{d1}_{d2}']
    d2d1_index = disease_pair_index[f'{d2}_{d1}']
    
    #get number of individuals
    N = len(ineligible_d_dict) #total number of exposed individuals
    sub_individual = [id_ for id_,x in ineligible_d_dict.items() if d1 not in x and d2 not in x]
    n = len(sub_individual) #total number of sub-cohort
    pair_threshold = (int(n * threshold_config_[0])
                      if threshold_config_[0] is not None else threshold_config_[1])
    #filter eligible_d_dict_withdate
    n_p1p2 = sum([d1 in x and d2 in x for x in eligible_d_dict_withdate.values()]) #number of individuals with both d1 and d2 diagnosis.
    p1 = sum([d1 in eligible_d_dict_withdate[id_] for id_ in sub_individual]) #number of individuals with d1 diagnosis.
    p2 = sum([d2 in eligible_d_dict_withdate[id_] for id_ in sub_individual]) #number of individuals with d2 diagnosis.
    n_com = sum([d1d2_index in x or d2d1_index in x for x in com_pair_dict.values()]) #number of individuals with non-temporal d1-d2 disease pair
    n_tra_d1_d2 = sum([d1d2_index in x for x in temporal_pair_dict.values()]) #number of individuals with temporal d1->d2 disease pair
    n_tra_d2_d1 = sum([d2d1_index in x for x in temporal_pair_dict.values()]) #number of individuals with temporal d2->d1 disease pair
    c = sum([n_com,n_tra_d1_d2,n_tra_d2_d1]) #number of individuals with temporal/non-temporal d1 and d2 disease pair
    
    if message:
        write_log(log_file_,f'{d1} and {d2}: {message}\n')
        return [d1,d2,f'{d1}-{d2}',N,n,n_p1p2,p1,p2,n_com,n_tra_d1_d2,n_tra_d2_d1,c]
    elif c<pair_threshold:
        write_log(log_file_,f'{d1} and {d2}: Less than threshold of {pair_threshold}\n')
        return [d1,d2,f'{d1}-{d2}',N,n,n_p1p2,p1,p2,n_com,n_tra_d1_d2,n_tra_d2_d1,c,f'Less than threshold of {pair_threshold}']
    else:
        phi,phi_theta,phi_p = com_phi(n,c,p1,p2,se_method_)
        rr,rr_theta,rr_p = com_rr(n,c,p1,p2,se_method_)
        write_log(log_file_,f'{d1} and {d2}: Done\n')
        return [d1,d2,f'{d1}-{d2}',N,n,n_p1p2,p1,p2,n_com,n_tra_d1_d2,n_tra_d2_d1,c,np.nan,phi,phi_theta,phi_p,rr,rr_theta,rr_p]


def com_phi_rr_wrapper(trajectory:dict,
                       d1:float,
                       d2:float,
                       message:str,
                       threshold_config:tuple,
                       log_file:str,
                       se_method:str) -> list:
    """
    Wrapper for com_phi_rr that assigns default values to global variables if needed.

    Parameters:
    ----------
    trajectory : dictionary
        DiseaseNetworkData.trajectory dictionary
    
    d1 : float
        Disease 1

    d2 : float
        Disease 2

    message : string
        additional comment
    
    threshold_config : tuple
        ``(proportion_threshold, n_threshold)`` configuration. The
        proportional threshold is calculated from this disease pair's
        eligible sub-cohort after history and sex restrictions.
    
    log_file : str
        Path and prefix for the log file
    
    Returns:
    ----------
    result : list
        list, comorbidity strength estimation results

    """
    # shared global data
    global trajectory_
    global threshold_config_
    global log_file_
    global se_method_
    # set global variables if not already defined
    trajectory_ = trajectory
    threshold_config_ = threshold_config
    log_file_ = log_file
    se_method_ = se_method
    # call the original function
    return com_phi_rr((d1,d2,message))

def init_worker(trajectory:dict,
                threshold_config:tuple,
                log_file:str,
                se_method:str):
    """
    This function sets up the necessary global variables for a worker process in a multiprocessing environment.
    It assigns the provided parameters to global variables that can be accessed by com_phi_rr function in the worker process.

    Parameters:
    ----------
    trajectory : dictionary
        DiseaseNetworkData.trajectory dictionary
    
    threshold_config : tuple
        ``(proportion_threshold, n_threshold)`` configuration passed to each
        worker for pair-specific threshold calculation.
    
    log_file : str
        Path and prefix for the log file

    Returns:
    ----------
    None

    """
    # shared global data
    global trajectory_
    global threshold_config_
    global log_file_
    global se_method_
    # set global variables if not already defined
    trajectory_ = trajectory
    threshold_config_ = threshold_config
    log_file_ = log_file
    se_method_ = se_method



















