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
from lifelines import CoxPHFitter
from statsmodels.duration.hazard_regression import PHReg
cph = CoxPHFitter()
import warnings
warnings.filterwarnings('ignore')

def cox_conditional(data:DiseaseNetworkData,n_threshold:int,phecode:float,sex_adjustment:bool,log_file:str):
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
    
    sex_adjustment : bool
        Whether sex should be included as an additional covariate in the Cox model.
    
    log_file : str
        Path and prefix for the log file.

    Returns:
    ----------
    result : list
        A list of the Cox analysis results.

    """
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
    covars = info_dict['phenotype_covariates_list']
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
    if level == 1:
        range_upper = round(phecode+0.99, 2)
    elif level == 2:
        range_upper = round(phecode+0.09, 2)
    n_step = int((range_upper - phecode) / 0.01 + 1)
    d_lst = np.linspace(phecode, range_upper, n_step)
    
    #df processing
    dataset_analysis = data.phenotype_df[covars+[id_col,index_date_col,end_date_col,exp_col,sex_col,matching_col]]
    if pd.isna(exl_range):
        dataset_analysis[exl_flag_col] = 0
    else:
        exl_list = []
        for range_ in exl_range.split(','):
            exl_lower,exl_higher = float(range_.split('-')[0]), float(range_.split('-')[1])
            n_step = int((exl_higher - exl_lower) / 0.01 + 1)
            exl_list += list(np.round(np.linspace(exl_lower,exl_higher,n_step),2))
        exl_list = set(exl_list)
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
    if len(dataset_analysis) == 0:
        result += [0,'Sex specific']
        with open(log_file,'ab') as f:
            f.write(f'No individuals remaining after filtering for phecode {phecode}\n'.encode())
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
    
    #return and save results if less than threshold
    if length < n_threshold:
        result += ['less than threshold','%i/%.2f (%.2f)' % (n_exp,time_exp,n_exp/time_exp),
                '%i/%.2f (%.2f)' % (n_unexp,time_unexp,n_unexp/time_unexp)]
        with open(log_file,'ab') as f:
            f.write(f'Number of cases less than threshold for phecode {phecode}\n'.encode())
        return result
    
    #restricted to groups with at least one case
    match_id = dataset_analysis[dataset_analysis[outcome_col]==1][matching_col].to_list()
    dataset_analysis = dataset_analysis[dataset_analysis[matching_col].isin(match_id)]
    
    #additionally include sex
    if sex_adjustment:
        covars += [sex_col]

    #error message none
    e = None
    
    try:
        model = PHReg(np.asarray(dataset_analysis[time_col],dtype=np.float32),
                    np.asarray(dataset_analysis[[exp_col]+covars],dtype=np.float32),
                    status=np.asarray(dataset_analysis[outcome_col],dtype=np.int32), 
                    strata=np.asarray(dataset_analysis[matching_col],dtype=np.int32))
        model_result = model.fit(method='bfgs',maxiter=300,disp=0)
        if pd.isna(model_result.params[0]) or pd.isna(model_result.bse[0]):
            model = cph.fit(dataset_analysis[[time_col,outcome_col,exp_col,matching_col]+covars],
                            fit_options=dict(step_size=0.2), duration_col=time_col, event_col=outcome_col,strata=[matching_col])
            result_temp = model.summary.loc[exp_col]
            result += ['fitted_lifelines','%i/%.2f (%.2f)' % (n_exp,time_exp,n_exp/time_exp),
                        '%i/%.2f (%.2f)' % (n_unexp,time_unexp,n_unexp/time_unexp)]
            result += [x for x in result_temp[['coef','se(coef)','p']]]
        else:
            result += ['fitted','%i/%.2f (%.2f)' % (n_exp,time_exp,n_exp/time_exp),
                                '%i/%.2f (%.2f)' % (n_unexp,time_unexp,n_unexp/time_unexp)]
            result += [model_result.params[0],model_result.bse[0],model_result.pvalues[0]]
    except:
        try:
            model = cph.fit(dataset_analysis[[time_col,outcome_col,exp_col,matching_col]+covars],
                            fit_options=dict(step_size=0.2), duration_col=time_col, event_col=outcome_col,strata=[matching_col])
            result_temp = model.summary.loc[exp_col]
            result += ['fitted_lifelines','%i/%.2f (%.2f)' % (n_exp,time_exp,n_exp/time_exp),
                        '%i/%.2f (%.2f)' % (n_unexp,time_unexp,n_unexp/time_unexp)]
            result += [x for x in result_temp[['coef','se(coef)','p']]]
        except Exception as e:
            result += [e,'%i/%.2f (%.2f)' % (n_exp,time_exp,n_exp/time_exp),
                        '%i/%.2f (%.2f)' % (n_unexp,time_unexp,n_unexp/time_unexp)]
    #print
    time_end = time.time()
    time_spend = time_end - time_start
    if e:
        with open(log_file,'ab') as f:
            f.write(f'An error occurred during the Cox model fitting for phecode {phecode} (elapsed {time_spend:.2f}s)\n'.encode())
    else:
        with open(log_file,'ab') as f:
            f.write(f'Cox model successfully fitted for phecode {phecode} (elapsed {time_spend:.2f}s)\n'.encode())
            
    return result

def cox_unconditional(data:DiseaseNetworkData,n_threshold:int,phecode:float,sex_adjustment:bool,log_file:str):
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
    
    sex_adjustment : bool
        Whether sex should be included as an additional covariate in the Cox model.
    
    log_file : str
        Path and prefix for the log file where output will be written.

    Returns:
    ----------
    result : list
        A list of the Cox analysis results.

    """
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
    covars = info_dict['phenotype_covariates_list']
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
    if level == 1:
        range_upper = round(phecode+0.99, 2)
    elif level == 2:
        range_upper = round(phecode+0.09, 2)
    n_step = int((range_upper - phecode) / 0.01 + 1)
    d_lst = np.linspace(phecode, range_upper, n_step)
    
    #df processing
    dataset_analysis = data.phenotype_df[covars+[id_col,index_date_col,end_date_col,exp_col,sex_col]]
    if pd.isna(exl_range):
        dataset_analysis[exl_flag_col] = 0
    else:
        exl_list = []
        for range_ in exl_range.split(','):
            exl_lower,exl_higher = float(range_.split('-')[0]), float(range_.split('-')[1])
            n_step = int((exl_higher - exl_lower) / 0.01 + 1)
            exl_list += list(np.round(np.linspace(exl_lower,exl_higher,n_step),2))
        exl_list = set(exl_list)
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
    if len(dataset_analysis) == 0:
        result += [0,'Sex specific']
        with open(log_file,'ab') as f:
            f.write(f'No individuals remaining after filtering for phecode {phecode}\n'.encode())
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
    
    #return and save results if less than threshold
    if length < n_threshold:
        result += ['less than threshold','%i/%.2f (%.2f)' % (n_exp,time_exp,n_exp/time_exp),
                '%i/%.2f (%.2f)' % (n_unexp,time_unexp,n_unexp/time_unexp)]
        with open(log_file,'ab') as f:
            f.write(f'Number of cases less than threshold for phecode {phecode}\n'.encode())
        return result
    
    #additionally include sex
    if sex_adjustment:
        covars += [sex_col]

    #error message none
    e = None
    
    try:
        model = PHReg(np.asarray(dataset_analysis[time_col],dtype=np.float32),
                    np.asarray(dataset_analysis[[exp_col]+covars],dtype=np.float32),
                    status=np.asarray(dataset_analysis[outcome_col],dtype=np.int32))
        model_result = model.fit(method='bfgs',maxiter=300,disp=0)
        if pd.isna(model_result.params[0]) or pd.isna(model_result.bse[0]):
            model = cph.fit(dataset_analysis[[time_col,outcome_col,exp_col]+covars],
                            fit_options=dict(step_size=0.2), duration_col=time_col, event_col=outcome_col)
            result_temp = model.summary.loc[exp_col]
            result += ['fitted_lifelines','%i/%.2f (%.2f)' % (n_exp,time_exp,n_exp/time_exp),
                        '%i/%.2f (%.2f)' % (n_unexp,time_unexp,n_unexp/time_unexp)]
            result += [x for x in result_temp[['coef','se(coef)','p']]]
        else:
            result += ['fitted','%i/%.2f (%.2f)' % (n_exp,time_exp,n_exp/time_exp),
                                '%i/%.2f (%.2f)' % (n_unexp,time_unexp,n_unexp/time_unexp)]
            result += [model_result.params[0],model_result.bse[0],model_result.pvalues[0]]
    except:
        try:
            model = cph.fit(dataset_analysis[[time_col,outcome_col,exp_col]+covars],
                            fit_options=dict(step_size=0.2), duration_col=time_col, event_col=outcome_col)
            result_temp = model.summary.loc[exp_col]
            result += ['fitted_lifelines','%i/%.2f (%.2f)' % (n_exp,time_exp,n_exp/time_exp),
                        '%i/%.2f (%.2f)' % (n_unexp,time_unexp,n_unexp/time_unexp)]
            result += [x for x in result_temp[['coef','se(coef)','p']]]
        except Exception as e:
            result += [e,'%i/%.2f (%.2f)' % (n_exp,time_exp,n_exp/time_exp),
                        '%i/%.2f (%.2f)' % (n_unexp,time_unexp,n_unexp/time_unexp)]
    #print
    time_end = time.time()
    time_spend = time_end - time_start
    if e:
        with open(log_file,'ab') as f:
            f.write(f'An error occurred during the Cox model fitting for phecode {phecode} (elapsed {time_spend:.2f}s)\n'.encode())
    else:
        with open(log_file,'ab') as f:
            f.write(f'Cox model successfully fitted for phecode {phecode} (elapsed {time_spend:.2f}s)\n'.encode())
            
    return result


    
    
    
    
    
    
    