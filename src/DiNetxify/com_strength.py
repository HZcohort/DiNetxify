# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 02:25:20 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

import pandas as pd
import numpy as np
from scipy.stats import t
#from .data_management import DiseaseNetworkData
from .utility import write_log

def com_rr(n:int,c:int,p1:int,p2:int):
    """
    Parameters
    ----------
    n : int
        total number of individuals.
    c : int
        number of individuals with temporal/non-temporal d1 and d2 disease pair
    p1 : int
        number of individuals with d1 diagnosis.
    p2 : int
        number of individuals with d2 diagnosis.

    Returns
    -------
        RR, RR_se and p-value

    """
    rr = (n*c)/(p1*p2)
    theta = (1/c + 1/((p1*p2)/n) - 1/n - 1/n)**0.5
    t_ = abs(np.log(rr)/theta)
    p = (1-t.cdf(t_,n))*2
    return rr,theta,p

def com_phi(n:int,c:int,p1:int,p2:int):
    """
    Parameters
    ----------
    n : int
        total number of individuals
    c : int
        number of individuals with temporal/non-temporal d1 and d2 disease pair
    p1 : int
        number of individuals with d1 diagnosis
    p2 : int
        number of individuals with d2 diagnosis

    Returns
    -------
        phi and p-value

    """
    phi = (c*n-p1*p2)/(((p1*p2)*(n-p1)*(n-p2))**0.5)
    try:
        z_phi = 0.5*np.log((1+phi)/(1-phi))
    except:
        phi -= 1e-3 #when phi == exactly 1 in some cases
        z_phi = 0.5*np.log((1+phi)/(1-phi))
    z_phi_theta = (1/(n-3))**0.5
    z_phi_t = abs(z_phi/z_phi_theta)
    p_phi = (1-t.cdf(z_phi_t,n))*2
    return phi,z_phi_theta,p_phi
    

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

    n_threshold : int
        Number of individuals threshold
    
    log_file : str
        Path and prefix for the log file
    Returns:
    ----------
    result : list
        list, comorbidity strength estimation results
    """
    # shared global data
    global trajectory_
    global n_threshold_
    global log_file_

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
    elif c<n_threshold_:
        write_log(log_file_,f'{d1} and {d2}: Less than threshold of {n_threshold_}\n')
        return [d1,d2,f'{d1}-{d2}',N,n,n_p1p2,p1,p2,n_com,n_tra_d1_d2,n_tra_d2_d1,c,f'Less than threshold of {n_threshold_}']
    else:
        phi,phi_theta,phi_p = com_phi(n,c,p1,p2)
        rr,rr_theta,rr_p = com_rr(n,c,p1,p2)
        write_log(log_file_,f'{d1} and {d2}: Done\n')
        return [d1,d2,f'{d1}-{d2}',N,n,n_p1p2,p1,p2,n_com,n_tra_d1_d2,n_tra_d2_d1,c,np.nan,phi,phi_theta,phi_p,rr,rr_theta,rr_p]


def com_phi_rr_wrapper(trajectory:dict,
                       d1:float,
                       d2:float,
                       message:str,
                       n_threshold:int,
                       log_file:str) -> list:
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
    
    n_threshold : int
        Number of individuals threshold
    
    log_file : str
        Path and prefix for the log file
    
    Returns:
    ----------
    result : list
        list, comorbidity strength estimation results

    """
    # shared global data
    global trajectory_
    global n_threshold_
    global log_file_
    # set global variables if not already defined
    trajectory_ = trajectory
    n_threshold_ = n_threshold
    log_file_ = log_file
    # call the original function
    return com_phi_rr((d1,d2,message))

def init_worker(trajectory:dict,
                n_threshold:int,
                log_file:str):
    """
    This function sets up the necessary global variables for a worker process in a multiprocessing environment.
    It assigns the provided parameters to global variables that can be accessed by com_phi_rr function in the worker process.

    Parameters:
    ----------
    trajectory : dictionary
        DiseaseNetworkData.trajectory dictionary
    
    n_threshold : int
        Number of individuals threshold
    
    log_file : str
        Path and prefix for the log file

    Returns:
    ----------
    None

    """
    # shared global data
    global trajectory_
    global n_threshold_
    global log_file_
    # set global variables if not already defined
    trajectory_ = trajectory
    n_threshold_ = n_threshold
    log_file_ = log_file





















