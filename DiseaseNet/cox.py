# -*- coding: utf-8 -*-
"""
Created on Thu Dec 5 19:48:09 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

#import sys
import pandas as pd
import numpy as np
import time
from .data_management import DiseaseNetworkData
from .utility import write_log
from statsmodels.duration.hazard_regression import PHReg
import warnings
warnings.filterwarnings('ignore')


def cox_conditional(data:DiseaseNetworkData,n_threshold:int,phecode:float,
                    covariates:list,log_file:str,lifelines_disable:bool):
    """
    Perfoming Cox conditional analysis based on the provided DiseaseNetworkData object.

    Parameters:
    ----------
    data : DiseaseNetworkData
        An DiseaseNetworkData object.
    
    n_threshold : int
        Number of cases threshold. Cox analysis are only conducted when number of cases larger than this 
        threshold among exposed group.

    phecode : float
        The outcome phecode for running the Cox analysis.
    
    covariates : list
        List of covariates to be adjusted in the Cox model.
    
    log_file : str
        Path and prefix for the log file.
    
    lifelines_disable : bool
        Whether to disable the use of lifelines. 
        While lifelines generally require a longer fitting time, they are more resilient to violations of model assumptions.

    Returns:
    ----------
    result : list
        A list of the Cox analysis results.

    """
    
    if lifelines_disable:
        cph = None
    else:
        from lifelines import CoxPHFitter
        cph = CoxPHFitter()
    
    #phecode information
    phecode_dict = data.phecode_info[phecode]
    disease_name = phecode_dict['phenotype']
    system = phecode_dict['category']
    level = phecode_dict['level']
    sex_specific = phecode_dict['sex']
    if sex_specific == 'Female':
        sex_code = 1
    elif sex_specific == 'Male':
        sex_code = 0
    else:
        sex_code = None
    exl_range = phecode_dict['phecode_exclude_range']
    
    #result list
    result = [phecode,disease_name,system,sex_specific]
    
    #information about the dataframe
    info_dict = data.get_attribute('phenotype_info')
    id_col = info_dict['phenotype_col_dict']['Participant ID']
    exp_col = info_dict['phenotype_col_dict']['Exposure']
    sex_col = info_dict['phenotype_col_dict']['Sex']
    index_date_col = info_dict['phenotype_col_dict']['Index date']
    end_date_col = info_dict['phenotype_col_dict']['End date']
    matching_col = info_dict['phenotype_col_dict']['Matching identifier']
    
    #history and diagnosis dict
    history = data.history
    diagnosis = data.diagnosis
    
    #default columns
    exl_flag_col = 'flag_exl'
    outcome_time_col = 'outcome_date'
    outcome_col = 'outcome'
    time_col = 'time_years'
    
    #time start
    time_start = time.time()
    
    #outcome disease list
    d_lst = phecode_dict['leaf_list']
    
    #df processing
    #make sure sex_col is always included but not duplicated
    if sex_col in covariates:
        dataset_analysis = data.phenotype_df[covariates+[id_col,index_date_col,end_date_col,exp_col,matching_col]]
    else:
        dataset_analysis = data.phenotype_df[covariates+[id_col,index_date_col,end_date_col,exp_col,sex_col,matching_col]]
    
    if pd.isna(exl_range):
        dataset_analysis[exl_flag_col] = 0
    else:
        exl_list = phecode_dict['exclude_list']
        #start check each individual eligibility
        exl_flag_lst = []
        for id_ in dataset_analysis[id_col].values:
            history_ = history[id_]
            if len(set(history_).intersection(exl_list))>0:
                exl_flag_lst.append(1)
            else:
                exl_flag_lst.append(0)
        dataset_analysis[exl_flag_col] = exl_flag_lst
    
    #sex specific
    if sex_code is not None:
        dataset_analysis[exl_flag_col] = dataset_analysis.apply(lambda row: 1 if row[sex_col] != sex_code 
                                                                else row[exl_flag_col],axis=1)
    #exclude eligible individuals
    dataset_analysis = dataset_analysis.loc[dataset_analysis[exl_flag_col]==0]
    
    #check number
    if len(dataset_analysis) == 0:
        result += [0,'Potentially sex specific']
        write_log(log_file,f'No individuals remaining after filtering for phecode {phecode}\n')
        return result
    
    #check number
    number_exposed = len(dataset_analysis[dataset_analysis[exp_col]==1])
    number_unexposed = len(dataset_analysis[dataset_analysis[exp_col]==0])
    if number_exposed == 0:
        result += [0,'Disease specific (zero exposed)']
        write_log(log_file,f'No exposed individuals remaining after filtering for phecode {phecode}\n')
        return result
    
    if number_unexposed == 0:
        result += [0,'Disease specific (zero unexposed)']
        write_log(log_file,f'No unexposed individuals remaining after filtering for phecode {phecode}\n')
        return result
    
    #define diagnosis time and outcome
    outcome_time_lst = []
    for id_ in dataset_analysis[id_col].values:
        diagnosis_ = diagnosis[id_]
        try:
            date = min([diagnosis_[x] for x in d_lst if x in diagnosis_])
        except:
            date = pd.NaT
        outcome_time_lst.append(date)
    dataset_analysis[outcome_time_col] = outcome_time_lst
    dataset_analysis[outcome_col] = dataset_analysis[outcome_time_col].apply(lambda x: 0 if pd.isna(x) else 1)
    dataset_analysis[end_date_col] = dataset_analysis[[end_date_col,outcome_time_col]].min(axis=1)
    
    #length
    length = len(dataset_analysis[(dataset_analysis[exp_col]==1) & (dataset_analysis[outcome_col]==1)])
    result += [length]
    
    #calculate time in years
    dataset_analysis[time_col] = (dataset_analysis[end_date_col] - dataset_analysis[index_date_col]).dt.days/365.25
    
    #calculate time at risk
    n_exp = len(dataset_analysis.loc[(dataset_analysis[exp_col]==1) & (dataset_analysis[outcome_col]==1)])
    n_unexp = len(dataset_analysis.loc[(dataset_analysis[exp_col]==0) & (dataset_analysis[outcome_col]==1)])
    time_exp = dataset_analysis.groupby(by=exp_col)[time_col].sum().loc[1]/1000
    time_unexp = dataset_analysis.groupby(by=exp_col)[time_col].sum().loc[0]/1000
    str_exp = '%i/%.2f (%.2f)' % (n_exp,time_exp,n_exp/time_exp)
    str_noexp = '%i/%.2f (%.2f)' % (n_unexp,time_unexp,n_unexp/time_unexp)
    
    #return and save results if less than threshold
    if length < n_threshold:
        result += [f'Less than threshold of {n_threshold}',str_exp,str_noexp]
        write_log(log_file,f'Number of cases {length} less than threshold {n_threshold} for phecode {phecode}\n')
        return result
    
    #exclude those with negative time
    dataset_analysis = dataset_analysis[dataset_analysis[time_col]>0]
    
    #restricted to groups with at least one case
    match_id = dataset_analysis[dataset_analysis[outcome_col]==1][matching_col].to_list()
    dataset_analysis = dataset_analysis[dataset_analysis[matching_col].isin(match_id)]
    
    #check var of covariates, remove these with var()==0
    for var in covariates:
        dataset_analysis[var] = dataset_analysis[var].astype(float)
        if dataset_analysis.groupby(by=matching_col)[var].var().mean() <= 0: #lowest var() allowed
            covariates.remove(var)

    #error message
    e_stats = None
    e_lifelines = None
    error_message = None
    
    try:
        model = PHReg(np.asarray(dataset_analysis[time_col],dtype=float),
                    np.asarray(dataset_analysis[[exp_col]+covariates],dtype=float),
                    status=np.asarray(dataset_analysis[outcome_col],dtype=int), 
                    strata=np.asarray(dataset_analysis[matching_col]))
        model_result = model.fit(method='bfgs',maxiter=300,disp=0)
        if pd.isna(model_result.params[0]) or pd.isna(model_result.bse[0]):
            e_stats = 'No converge for statsmodels Cox'
            model = cph.fit(dataset_analysis[[time_col,outcome_col,exp_col,matching_col]+covariates],
                            fit_options=dict(step_size=0.2), duration_col=time_col, event_col=outcome_col,strata=[matching_col])
            result_temp = model.summary.loc[exp_col]
            result += ['fitted_lifelines',str_exp,str_noexp]
            result += [x for x in result_temp[['coef','se(coef)','p']]]
        else:
            result += ['fitted',str_exp,str_noexp]
            result += [model_result.params[0],model_result.bse[0],model_result.pvalues[0]]
    except Exception as e:
        if e_stats:
            e_lifelines = e
        else:
            e_stats = e
        try:
            model = cph.fit(dataset_analysis[[time_col,outcome_col,exp_col,matching_col]+covariates],
                            fit_options=dict(step_size=0.2), duration_col=time_col, event_col=outcome_col,strata=[matching_col])
            result_temp = model.summary.loc[exp_col]
            result += ['fitted_lifelines',str_exp,str_noexp]
            result += [x for x in result_temp[['coef','se(coef)','p']]]
        except Exception as e:
            if e_lifelines:
                None
            else:
                e_lifelines = e
            if lifelines_disable:
                error_message = e_stats
            else:
                error_message = f'{e_stats} (statsmodels); {e_lifelines} (lifelines)'
            result += [error_message,str_exp,str_noexp]
    #print
    time_end = time.time()
    time_spend = time_end - time_start
    if error_message:
        write_log(log_file,f'An error occurred during the Cox model fitting for phecode {phecode} (elapsed {time_spend:.2f}s)\n{error_message}\n')
    else:
        write_log(log_file,f'Cox model successfully fitted for phecode {phecode} (elapsed {time_spend:.2f}s)\n')
            
    return result

def cox_unconditional(data:DiseaseNetworkData,n_threshold:int,phecode:float,
                      covariates:list,log_file:str,lifelines_disable:bool):
    """
    Perfoming Cox unconditional analysis based on the provided DiseaseNetworkData object.

    Parameters:
    ----------
    data : DiseaseNetworkData
        An DiseaseNetworkData object.
    
    n_threshold : int
        Number of cases threshold. Cox analysis are only conducted when number of cases larger than this 
        threshold among exposed group.

    phecode : float
        The outcome phecode for running the Cox analysis.
    
    covariates : list
        List of covariates to be adjusted in the Cox model.
    
    log_file : str
        Path and prefix for the log file where output will be written.
    
    lifelines_disable : bool
        Whether to disable the use of lifelines. 

    Returns:
    ----------
    result : list
        A list of the Cox analysis results.

    """
    
    if lifelines_disable:
        cph = None
    else:
        from lifelines import CoxPHFitter
        cph = CoxPHFitter()
    
    #phecode information
    phecode_dict = data.phecode_info[phecode]
    disease_name = phecode_dict['phenotype']
    system = phecode_dict['category']
    level = phecode_dict['level']
    sex_specific = phecode_dict['sex']
    if sex_specific == 'Female':
        sex_code = 1
    elif sex_specific == 'Male':
        sex_code = 0
    else:
        sex_code = None
    exl_range = phecode_dict['phecode_exclude_range']
    
    #result list
    result = [phecode,disease_name,system]
    
    #information about the dataframe
    info_dict = data.get_attribute('phenotype_info')
    id_col = info_dict['phenotype_col_dict']['Participant ID']
    exp_col = info_dict['phenotype_col_dict']['Exposure']
    sex_col = info_dict['phenotype_col_dict']['Sex']
    index_date_col = info_dict['phenotype_col_dict']['Index date']
    end_date_col = info_dict['phenotype_col_dict']['End date']
    
    #history and diagnosis dict
    history = data.history
    diagnosis = data.diagnosis
    
    #default columns
    exl_flag_col = 'flag_exl'
    outcome_time_col = 'outcome_date'
    outcome_col = 'outcome'
    time_col = 'time_years'
    
    #time start
    time_start = time.time()
    
    #outcome disease list
    d_lst = phecode_dict['leaf_list']
    
    #df processing
    #make sure sex_col is always included but not duplicated
    if sex_col in covariates:
        dataset_analysis = data.phenotype_df[covariates+[id_col,index_date_col,end_date_col,exp_col]]
    else:
        dataset_analysis = data.phenotype_df[covariates+[id_col,index_date_col,end_date_col,exp_col,sex_col]]

    if pd.isna(exl_range):
        dataset_analysis[exl_flag_col] = 0
    else:
        exl_list = phecode_dict['exclude_list']
        #start check
        exl_flag_lst = []
        for id_ in dataset_analysis[id_col].values:
            history_ = history[id_]
            if len(set(history_).intersection(exl_list))>0:
                exl_flag_lst.append(1)
            else:
                exl_flag_lst.append(0)
        dataset_analysis[exl_flag_col] = exl_flag_lst
    
    #sex specific
    if sex_code:
        dataset_analysis[exl_flag_col] = dataset_analysis.apply(lambda row: 1 if row[sex_col] != sex_code 
                                                                else row[exl_flag_col],axis=1)
    #exclude eligible individuals
    dataset_analysis = dataset_analysis.loc[dataset_analysis[exl_flag_col]==0]
    
    #check number
    if len(dataset_analysis) == 0:
        result += [0,'Sex specific potentially']
        write_log(log_file,f'No individuals remaining after filtering for phecode {phecode}\n')
        return result
    
    #check number
    number_exposed = len(dataset_analysis[dataset_analysis[exp_col]==1])
    number_unexposed = len(dataset_analysis[dataset_analysis[exp_col]==0])
    if number_exposed == 0:
        result += [0,'Disease specific (zero exposed)']
        write_log(log_file,f'No exposed individuals remaining after filtering for phecode {phecode}\n')
        return result
    if number_unexposed == 0:
        result += [0,'Disease specific (zero unexposed)']
        write_log(log_file,f'No unexposed individuals remaining after filtering for phecode {phecode}\n')
        return result
    
    #define diagnosis time and outcome
    outcome_time_lst = []
    for id_ in dataset_analysis[id_col].values:
        diagnosis_ = diagnosis[id_]
        try:
            date = min([diagnosis_[x] for x in d_lst if x in diagnosis_])
        except:
            date = pd.NaT
        outcome_time_lst.append(date)
    dataset_analysis[outcome_time_col] = outcome_time_lst
    dataset_analysis[outcome_col] = dataset_analysis[outcome_time_col].apply(lambda x: 0 if pd.isna(x) else 1)
    dataset_analysis[end_date_col] = dataset_analysis[[end_date_col,outcome_time_col]].min(axis=1)
    
    #length
    length = len(dataset_analysis[(dataset_analysis[exp_col]==1) & (dataset_analysis[outcome_col]==1)])
    result += [length]
    
    #calculate time in years
    dataset_analysis[time_col] = (dataset_analysis[end_date_col] - dataset_analysis[index_date_col]).dt.days/365.25
    
    #calculate time at risk
    n_exp = len(dataset_analysis.loc[(dataset_analysis[exp_col]==1) & (dataset_analysis[outcome_col]==1)])
    n_unexp = len(dataset_analysis.loc[(dataset_analysis[exp_col]==0) & (dataset_analysis[outcome_col]==1)])
    time_exp = dataset_analysis.groupby(by=exp_col)[time_col].sum().loc[1]/1000
    time_unexp = dataset_analysis.groupby(by=exp_col)[time_col].sum().loc[0]/1000
    str_exp = '%i/%.2f (%.2f)' % (n_exp,time_exp,n_exp/time_exp)
    str_noexp = '%i/%.2f (%.2f)' % (n_unexp,time_unexp,n_unexp/time_unexp)
    
    #return and save results if less than threshold
    if length < n_threshold:
        result += [f'Less than threshold of {n_threshold}',str_exp,str_noexp]
        write_log(log_file,f'Number of cases {length} less than threshold {n_threshold} for phecode {phecode}\n')
        return result
    
    #exclude those with negative time
    dataset_analysis = dataset_analysis[dataset_analysis[time_col]>0]
    
    #check var of covariates, remove these with var()==0
    for var in covariates:
        if dataset_analysis[var].var() <= 0: #lowest var() allowed
            covariates.remove(var)

    #error message
    e_stats = None
    e_lifelines = None
    error_message = None
    
    try:
        model = PHReg(np.asarray(dataset_analysis[time_col],dtype=float),
                    np.asarray(dataset_analysis[[exp_col]+covariates],dtype=float),
                    status=np.asarray(dataset_analysis[outcome_col],dtype=int))
        model_result = model.fit(method='bfgs',maxiter=300,disp=0)
        if pd.isna(model_result.params[0]) or pd.isna(model_result.bse[0]):
            e_stats = 'No converge for statsmodels Cox'
            model = cph.fit(dataset_analysis[[time_col,outcome_col,exp_col]+covariates],
                            fit_options=dict(step_size=0.2), duration_col=time_col, event_col=outcome_col)
            result_temp = model.summary.loc[exp_col]
            result += ['fitted_lifelines',str_exp,str_noexp]
            result += [x for x in result_temp[['coef','se(coef)','p']]]
        else:
            result += ['fitted',str_exp,str_noexp]
            result += [model_result.params[0],model_result.bse[0],model_result.pvalues[0]]
    except Exception as e:
        if e_stats:
            e_lifelines = e
        else:
            e_stats = e
        try:
            model = cph.fit(dataset_analysis[[time_col,outcome_col,exp_col]+covariates],
                            fit_options=dict(step_size=0.2), duration_col=time_col, event_col=outcome_col)
            result_temp = model.summary.loc[exp_col]
            result += ['fitted_lifelines',str_exp,str_noexp]
            result += [x for x in result_temp[['coef','se(coef)','p']]]
        except Exception as e:
            if e_lifelines:
                None
            else:
                e_lifelines = e
            if lifelines_disable:
                error_message = e_stats
            else:
                error_message = f'{e_stats} (statsmodels); {e_lifelines} (lifelines)'
            result += [error_message,str_exp,str_noexp]
    #print
    time_end = time.time()
    time_spend = time_end - time_start
    if error_message:
        write_log(log_file,f'An error occurred during the Cox model fitting for phecode {phecode} (elapsed {time_spend:.2f}s)\n{error_message}\n')
    else:
        write_log(log_file,f'Cox model successfully fitted for phecode {phecode} (elapsed {time_spend:.2f}s)\n')
            
    return result







