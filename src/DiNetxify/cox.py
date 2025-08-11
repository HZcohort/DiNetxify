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
from .utility import write_log, check_variance_vif_single, check_history_exclusion, time_first_diagnosis
from statsmodels.duration.hazard_regression import PHReg
import warnings
warnings.filterwarnings('ignore')


def cox_conditional(phecode: float):
    """
    Perfoming Cox conditional analysis based on the provided DiseaseNetworkData object.

    Parameters:
    ----------
    phecode : float
        The outcome phecode for running the Cox analysis.

    Global Variables:
    ----------
    data : DiseaseNetworkData
        An DiseaseNetworkData object.

    n_threshold : int
        Number of cases threshold. Cox analysis are only conducted when number of cases larger than this threshold among exposed group.
    
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
    #shared global variables
    global data_
    global n_threshold_
    global covariates_
    global log_file_
    global lifelines_disable_

    if lifelines_disable_:
        cph = None
    else:
        from lifelines import CoxPHFitter
        cph = CoxPHFitter()
    
    #phecode information
    phecode_dict = data_.phecode_info[phecode]
    disease_name = phecode_dict['phenotype']
    system = phecode_dict['category']
    sex_specific = phecode_dict['sex']
    if sex_specific == 'Female':
        sex_code = 1
    elif sex_specific == 'Male':
        sex_code = 0
    else:
        sex_code = None
    exl_list = phecode_dict['exclude_list']
    min_icd_num = data_.min_required_icd_codes
    
    #result list
    result = [phecode,disease_name,system,sex_specific]
    
    #information about the dataframe
    info_dict = data_.get_attribute('phenotype_info')
    id_col = info_dict['phenotype_col_dict']['Participant ID']
    exp_col = info_dict['phenotype_col_dict']['Exposure']
    sex_col = info_dict['phenotype_col_dict']['Sex']
    index_date_col = info_dict['phenotype_col_dict']['Index date']
    end_date_col = info_dict['phenotype_col_dict']['End date']
    matching_col = info_dict['phenotype_col_dict']['Match ID']
    
    #history and diagnosis dict
    history = data_.history
    diagnosis = data_.diagnosis
    n_diagnosis = data_.n_diagnosis
    
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
    if sex_col in covariates_:
        dataset_analysis = data_.phenotype_df[covariates_+[id_col,index_date_col,end_date_col,exp_col,matching_col]]
    else:
        dataset_analysis = data_.phenotype_df[covariates_+[id_col,index_date_col,end_date_col,exp_col,sex_col,matching_col]]
    
    if len(exl_list)==0: #no list of exlusion phecode
        dataset_analysis[exl_flag_col] = 0
    else:
        #start check each individual eligibility
        exl_flag_lst = []
        for id_ in dataset_analysis[id_col].values:
            exl_flag_lst.append(check_history_exclusion(exl_list,history[id_],n_diagnosis[id_],min_icd_num))
        dataset_analysis[exl_flag_col] = exl_flag_lst
    
    #sex specific
    if sex_code is not None:
        dataset_analysis[exl_flag_col] = dataset_analysis.apply(lambda row: 1 if row[sex_col] != sex_code 
                                                                else row[exl_flag_col],axis=1)

    #for mathced cohort study, exclude eligible exposed along with their matched unexposed, and unexposed
    match_id_exl = dataset_analysis[(dataset_analysis[exl_flag_col]==1) & 
                                    (dataset_analysis[exp_col]==1)][matching_col].to_list()
    dataset_analysis = dataset_analysis[(dataset_analysis[exl_flag_col]==0) & 
                                        ~(dataset_analysis[matching_col].isin(match_id_exl))]
    
    #check number
    if len(dataset_analysis) == 0:
        result += [0, 'Potentially sex specific']
        write_log(log_file_,f'No individuals remaining after filtering for phecode {phecode}\n')
        return result
    
    #check number
    number_exposed = len(dataset_analysis[dataset_analysis[exp_col]==1])
    number_unexposed = len(dataset_analysis[dataset_analysis[exp_col]==0])
    if number_exposed == 0:
        result += [0, 'Disease specific (zero exposed)']
        write_log(log_file_,f'No exposed individuals remaining after filtering for phecode {phecode}\n')
        return result
    
    if number_unexposed == 0:
        result += [0, 'Disease specific (zero unexposed)']
        write_log(log_file_,f'No unexposed individuals remaining after filtering for phecode {phecode}\n')
        return result
    
    #define diagnosis time and outcome
    outcome_time_lst = []
    for id_ in dataset_analysis[id_col].values:
        date = time_first_diagnosis(d_lst, diagnosis[id_], n_diagnosis[id_], min_icd_num)
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
    if length < n_threshold_:
        result += [f'Less than threshold of {n_threshold_}',str_exp,str_noexp]
        write_log(log_file_,f'Number of cases {length} less than threshold {n_threshold_} for phecode {phecode}\n')
        return result
    
    #exclude those with negative time
    dataset_analysis = dataset_analysis[dataset_analysis[time_col]>0]
    
    #restricted to groups with at least one case
    match_id = dataset_analysis[dataset_analysis[outcome_col]==1][matching_col].to_list()
    dataset_analysis = dataset_analysis[dataset_analysis[matching_col].isin(match_id)]
    
    #check the covariates vif
    del_covariates = check_variance_vif_single(
        dataset_analysis,
        [exp_col],
        covariates_,
        vif_cutoff='phenotypic_covar', 
        group_col=matching_col
    )
    final_covariates = [x for x in covariates_ if x not in del_covariates]
    
    #error message
    e_stats = None
    e_lifelines = None
    error_message = None

    try:
        model = PHReg(
            np.asarray(dataset_analysis[time_col],dtype=float),
            np.asarray(dataset_analysis[[exp_col]+final_covariates],dtype=float),
            status=np.asarray(dataset_analysis[outcome_col],dtype=int), 
            strata=np.asarray(dataset_analysis[matching_col])
        )
        model_result = model.fit(method='bfgs',maxiter=300,disp=0)
        if pd.isna(model_result.params[0]) or pd.isna(model_result.bse[0]):
            e_stats = 'No converge for statsmodels Cox'
            model = cph.fit(
                dataset_analysis[[time_col,outcome_col,exp_col,matching_col]+final_covariates],
                fit_options=dict(step_size=0.2), 
                duration_col=time_col, 
                event_col=outcome_col,
                strata=[matching_col]
            )
            result_temp = model.summary.loc[exp_col]
            result += [f'fitted_lifelines and delete the covariate(s): {del_covariates}',str_exp,str_noexp]
            result += [x for x in result_temp[['coef','se(coef)','p']]]
        else:
            result += [f'fitted and delete the covariate(s): {del_covariates}',str_exp,str_noexp]
            result += [model_result.params[0],model_result.bse[0],model_result.pvalues[0]]
    except Exception as e:
        if e_stats:
            e_lifelines = e
        else:
            e_stats = e
        try:
            model = cph.fit(
                dataset_analysis[[time_col,outcome_col,exp_col,matching_col]+final_covariates],
                fit_options=dict(step_size=0.2), 
                duration_col=time_col, 
                event_col=outcome_col,
                strata=[matching_col]
            )
            result_temp = model.summary.loc[exp_col]
            result += [f'fitted_lifelines and delete the covariate(s): {del_covariates}',str_exp,str_noexp]
            result += [x for x in result_temp[['coef','se(coef)','p']]]
        except Exception as e:
            if e_lifelines:
                None
            else:
                e_lifelines = e
            if lifelines_disable_:
                error_message = e_stats
            else:
                error_message = f'{e_stats} (statsmodels); {e_lifelines} (lifelines)'
            result += [error_message,str_exp,str_noexp]
    #print
    time_end = time.time()
    time_spend = time_end - time_start
    if error_message:
        write_log(log_file_,f'An error occurred during the Cox model fitting for phecode {phecode} (elapsed {time_spend:.2f}s)\n{error_message}\n')
    else:
        write_log(log_file_,f'Cox model successfully fitted for phecode {phecode} (elapsed {time_spend:.2f}s)\n')  
    return result

def cox_conditional_wrapper(
    phecode: str, 
    data: DiseaseNetworkData, 
    covariates: list, 
    n_threshold: int, 
    log_file: str, 
    lifelines_disable: bool
) -> pd.DataFrame:
    """
    Wrapper for cox_conditional that assigns default values to global variables if needed.

    Parameters:
    ----------
    phecode : float
        The outcome phecode for running the Cox analysis.
    data : DiseaseNetworkData
        An DiseaseNetworkData object.
    n_threshold : int
        Number of cases threshold. Cox analysis are only conducted when number of cases larger than this threshold among exposed group.
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
    # Assign global variables if not already set
    global data_
    global covariates_
    global n_threshold_
    global log_file_
    global lifelines_disable_
    # Set global variables if not already defined
    data_ = data
    covariates_ = covariates
    n_threshold_ = n_threshold
    log_file_ = log_file
    lifelines_disable_ = lifelines_disable
    # Call the original function
    return cox_conditional(phecode)

def cox_unconditional(phecode:float):
    """
    Perfoming Cox unconditional analysis based on the provided DiseaseNetworkData object.

    Parameters:
    ----------
    phecode : float
            The outcome phecode for running the Cox analysis.

    Global Variables:
    ----------
    data : DiseaseNetworkData
        An DiseaseNetworkData object.
    
    n_threshold : int
        Number of cases threshold. Cox analysis are only conducted when number of cases larger than this threshold among exposed group.
    
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
    #shared global data
    global data_
    global covariates_
    global n_threshold_
    global log_file_
    global lifelines_disable_
    
    if lifelines_disable_:
        cph = None
    else:
        from lifelines import CoxPHFitter
        cph = CoxPHFitter()
    
    #phecode information
    phecode_dict = data_.phecode_info[phecode]
    disease_name = phecode_dict['phenotype']
    system = phecode_dict['category']
    sex_specific = phecode_dict['sex']
    if sex_specific == 'Female':
        sex_code = 1
    elif sex_specific == 'Male':
        sex_code = 0
    else:
        sex_code = None
    exl_list = phecode_dict['exclude_list']
    min_icd_num = data_.min_required_icd_codes
    
    #result list
    result = [phecode,disease_name,system,sex_specific]
    
    #information about the dataframe
    info_dict = data_.get_attribute('phenotype_info')
    id_col = info_dict['phenotype_col_dict']['Participant ID']
    exp_col = info_dict['phenotype_col_dict']['Exposure']
    sex_col = info_dict['phenotype_col_dict']['Sex']
    index_date_col = info_dict['phenotype_col_dict']['Index date']
    end_date_col = info_dict['phenotype_col_dict']['End date']
    
    #history and diagnosis dict
    history = data_.history
    diagnosis = data_.diagnosis
    n_diagnosis = data_.n_diagnosis
    
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
    if sex_col in covariates_:
        dataset_analysis = data_.phenotype_df[covariates_+[id_col,index_date_col,end_date_col,exp_col]]
    else:
        dataset_analysis = data_.phenotype_df[covariates_+[id_col,index_date_col,end_date_col,exp_col,sex_col]]

    if len(exl_list)==0: #no list of exlusion phecode
        dataset_analysis[exl_flag_col] = 0
    else:
        #start check each individual eligibility
        exl_flag_lst = []
        for id_ in dataset_analysis[id_col].values:
            exl_flag_lst.append(check_history_exclusion(exl_list,history[id_],n_diagnosis[id_],min_icd_num))
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
        write_log(log_file_,f'No individuals remaining after filtering for phecode {phecode}\n')
        return result

    if data_.study_design != "exposed-only cohort":
    #check number
        number_exposed = len(dataset_analysis[dataset_analysis[exp_col]==1])
        number_unexposed = len(dataset_analysis[dataset_analysis[exp_col]==0])
        if number_exposed == 0:
            result += [0,'Disease specific (zero exposed)']
            write_log(log_file_,f'No exposed individuals remaining after filtering for phecode {phecode}\n')
            return result
        if number_unexposed == 0:
            result += [0,'Disease specific (zero unexposed)']
            write_log(log_file_,f'No unexposed individuals remaining after filtering for phecode {phecode}\n')
            return result
    
    #define diagnosis time and outcome
    outcome_time_lst = []
    for id_ in dataset_analysis[id_col].values:
        date = time_first_diagnosis(d_lst,diagnosis[id_],n_diagnosis[id_],min_icd_num)
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
    str_exp = '%i/%.2f (%.2f)' % (n_exp,time_exp,n_exp/time_exp)
    if data_.study_design != "exposed-only cohort":
        time_unexp = dataset_analysis.groupby(by=exp_col)[time_col].sum().loc[0]/1000
        str_noexp = '%i/%.2f (%.2f)' % (n_unexp,time_unexp,n_unexp/time_unexp)
    
    #return and save results if less than threshold
    if length < n_threshold_ and data_.study_design != "exposed-only cohort":
        result += [f'Less than the threshold of {n_threshold_}',str_exp,str_noexp]
        write_log(log_file_, f'Number of cases {length} less than the threshold of {n_threshold_} for phecode {phecode}\n')
        return result
    elif length < n_threshold_ and data_.study_design == "exposed-only cohort":
        result += [f'Less than the threshold of {n_threshold_}',str_exp]
        write_log(log_file_, f'Number of cases {length} less than the threshold of {n_threshold_} for phecode {phecode}\n')
        return result
    elif length >= n_threshold_ and data_.study_design == "exposed-only cohort":
        result += [f"Reached the threshold of {n_threshold_}",str_exp]
        write_log(log_file_, f"Number of cases reached the threshold of {n_threshold_} for phecode {phecode}\n")
        return result
    
    #exclude those with negative time
    dataset_analysis = dataset_analysis[dataset_analysis[time_col]>0]
    
    #check the covariates vif
    del_covariates = check_variance_vif_single(dataset_analysis,
                                               [exp_col],covariates_,
                                               vif_cutoff='phenotypic_covar')
    final_covariates = [x for x in covariates_ if x not in del_covariates]

    #error message
    e_stats = None
    e_lifelines = None
    error_message = None
    
    try:
        model = PHReg(
            np.asarray(dataset_analysis[time_col],dtype=float),
            np.asarray(dataset_analysis[[exp_col]+final_covariates],dtype=float),
            status=np.asarray(dataset_analysis[outcome_col],dtype=int)
        )
        model_result = model.fit(method='bfgs',maxiter=300,disp=0)
        if pd.isna(model_result.params[0]) or pd.isna(model_result.bse[0]):
            e_stats = 'No converge for statsmodels Cox'
            model = cph.fit(
                dataset_analysis[[time_col,outcome_col,exp_col]+final_covariates],
                fit_options=dict(step_size=0.2), 
                duration_col=time_col, 
                event_col=outcome_col
            )
            result_temp = model.summary.loc[exp_col]
            result += [f'fitted_lifelines and delete the covariate(s): {del_covariates}',str_exp,str_noexp]
            result += [x for x in result_temp[['coef','se(coef)','p']]]
        else:
            result += [f'fitted and delete the covariate(s): {del_covariates}',str_exp,str_noexp]
            result += [model_result.params[0],model_result.bse[0],model_result.pvalues[0]]
    except Exception as e:
        if e_stats:
            e_lifelines = e
        else:
            e_stats = e
        try:
            model = cph.fit(
                dataset_analysis[[time_col,outcome_col,exp_col]+final_covariates],
                fit_options=dict(step_size=0.2), 
                duration_col=time_col, 
                event_col=outcome_col
            )
            result_temp = model.summary.loc[exp_col]
            result += [f'fitted_lifelines and delete the covariate(s): {del_covariates}',str_exp,str_noexp]
            result += [x for x in result_temp[['coef','se(coef)','p']]]
        except Exception as e:
            if e_lifelines:
                None
            else:
                e_lifelines = e
            if lifelines_disable_:
                error_message = e_stats
            else:
                error_message = f'{e_stats} (statsmodels); {e_lifelines} (lifelines)'
            result += [error_message,str_exp,str_noexp]
    #print
    time_end = time.time()
    time_spend = time_end - time_start
    if error_message:
        write_log(log_file_,f'An error occurred during the Cox model fitting for phecode {phecode} (elapsed {time_spend:.2f}s)\n{error_message}\n')
    else:
        write_log(log_file_,f'Cox model successfully fitted for phecode {phecode} (elapsed {time_spend:.2f}s)\n')
            
    return result

def cox_unconditional_wrapper(
    phecode:str, 
    data:DiseaseNetworkData, 
    covariates:list, 
    n_threshold:int, 
    log_file:str, 
    lifelines_disable:bool
) -> None:
    """
    Wrapper for cox_unconditional that assigns default values to global variables if needed.

    Parameters:
    ----------
    phecode : float
        The outcome phecode for running the Cox analysis.
    data : DiseaseNetworkData
        An DiseaseNetworkData object.
    n_threshold : int
        Number of cases threshold. Cox analysis are only conducted when number of cases larger than this threshold among exposed group.
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
    # Assign global variables if not already set
    global data_
    global covariates_
    global n_threshold_
    global log_file_
    global lifelines_disable_
    # Set global variables if not already defined
    data_ = data
    covariates_ = covariates
    n_threshold_ = n_threshold
    log_file_ = log_file
    lifelines_disable_ = lifelines_disable
    # Call the original function
    return cox_unconditional(phecode)

def init_worker(
    data:str,
    covariates:list,
    n_threshold:int,
    log_file:str,
    lifelines_disable:bool
) -> None:
    """
    This function sets up the necessary global variables for a worker process in a multiprocessing environment.
    It assigns the provided parameters to global variables that can be accessed by cox_unconditional and cox_conditional function in the worker process.

    Parameters:
    ----------
    data : DiseaseNetworkData
        An DiseaseNetworkData object.
    n_threshold : int
        Number of cases threshold. Cox analysis are only conducted when number of cases larger than this threshold among exposed group.
    covariates : list
        List of covariates to be adjusted in the Cox model.
    log_file : str
        Path and prefix for the log file where output will be written.
    lifelines_disable : bool
        Whether to disable the use of lifelines. 

    Returns:
    ----------
    None

    """
    #shared global data
    global data_
    global covariates_
    global n_threshold_
    global log_file_
    global lifelines_disable_
    #assign values
    data_ = data
    covariates_ = covariates
    n_threshold_ = n_threshold
    log_file_ = log_file
    lifelines_disable_ = lifelines_disable


    



