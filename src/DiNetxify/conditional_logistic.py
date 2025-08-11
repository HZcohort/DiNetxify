# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 14:56:02 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

import pandas as pd
import numpy as np
from statsmodels.discrete.discrete_model import Logit
from statsmodels.discrete.conditional_models import ConditionalLogit,ConditionalResultsWrapper
import time
import gc
from .utility import write_log,find_best_alpha_and_vars,check_variance_vif_single

import warnings
warnings.filterwarnings('ignore')

def logistic_model(args):
    """
    Fit a conditional LR model to verify the comorbidity association between a temporal disease pair.

    Parameters
    ----------
    d1_lst : list, list of disease 1
    d2 : float, phecode 2.

    Global Variables:
    ----------
    phenotype_df_exposed : pd.DataFrame, phenotypic data for exposed individuals only.
    id_col : str, id column
    end_date_col : str, date of end follow-up
    trajectory_ineligible : dict, trajectory ineligible disease dictionary.
    trajectory_temporal : dict, temporal disease pair dictionary
    trajectory_eligible_withdate : dict, trajectory eligible disease (with date) dictionary.
    all_diagnosis_level : list, list of all diagnosed phecodes, with phecode truncated to corresponding level
    covariates : list, list of covariates to be included in the model.
    all_diseases_lst : list, list of other diseases to be included.
    matching_var_dict : dict, matching variables and the criteria used for incidence density sampling.
    matching_n : int, the maximum number of matched controls for each case.
    log_file : str, Path and prefix for the log file
    parameters : dict, other arguments, including method and the associated parameters.

    Returns
    -------
    Result list

    """
    #shared global data
    global phenotype_df_exposed_
    global id_col_
    global end_date_col_
    global trajectory_ineligible_
    global min_interval_
    global max_interval_
    global trajectory_eligible_withdate_
    global all_diagnosis_level_
    global covariates_
    global all_diseases_lst_
    global matching_var_dict_
    global matching_n_
    global max_n_cases_
    global log_file_
    global parameters_

    d1_lst, d2 = args
    #default columns
    d1_date_col = 'd1_date'
    d2_date_col = 'd2_date'
    d1_col = 'd1'
    d2_col = 'd2'
    outcome_date_col = 'outcome_date'
    mathcing_id_col = 'group_matching_ids'
    constant_col = 'constant'

    #method and parameters
    enforce_time_interval = parameters_['enforce_time_interval']
    method = parameters_['method']
    if method == 'RPCN':
        #alpha_initial = [1, 10, 20, 30, 40, 50] #alpha starting value range if using auto_penalty
        auto_penalty = parameters_['auto_penalty']
        alpha_single = parameters_['alpha']
        alpha_range = parameters_['alpha_range']
        scaling_factor = parameters_['scaling_factor']
    elif method == 'PCN_PCA':
        pca_number = parameters_.get('explained_variance',parameters_.get('n_PC')) #retrive explained_variance first if given, otherwise use n_PC
    
    #time and message
    time_start = time.time()
    message = f'Performing incidence density matching for {d2}: '
    
    #filtering the dataframe first
    N = len(phenotype_df_exposed_)
    vars_for_matching = [id_col_,end_date_col_]+list(matching_var_dict_.keys())
    d2_eligible_lst = [id_ for id_,vals in trajectory_ineligible_.items() if d2 not in vals]
    df_for_matching = phenotype_df_exposed_[phenotype_df_exposed_[id_col_].isin(d2_eligible_lst)][vars_for_matching]

    #check the number of cases
    df_for_matching[d2_date_col] = df_for_matching[id_col_].apply(lambda x: trajectory_eligible_withdate_[x].get(d2,pd.NaT))
    n_cases = len(df_for_matching[~df_for_matching[d2_date_col].isna()])
    if n_cases > max_n_cases_:
        #randomly select max_n_cases from cases and concatenate with the rest controls
        cases = df_for_matching[~df_for_matching[d2_date_col].isna()].sample(max_n_cases_)
        controls = df_for_matching[df_for_matching[d2_date_col].isna()]
        df_for_matching = pd.concat([cases,controls])

    #matching
    df_matched = matching_ids(df_for_matching,matching_var_dict_,matching_n_,id_col_,d2_date_col,end_date_col_,
                              d2_col,outcome_date_col,mathcing_id_col)
    del df_for_matching
    gc.collect()

    # time inforamtion
    time_end = time.time()
    #calculate number of cases and matched controls
    n_cases = len(df_matched[df_matched[d2_col]==1])
    n_controls = len(df_matched[df_matched[d2_col]==0])
    message += f'time spent on matching: {time_end-time_start:.2f}s; number of cases {n_cases:,}, number of matched controls {n_controls:,}\n'

    #list for saving results
    result_all = []

    #------------------------------------loop all the d1
    for d1 in d1_lst:
        #variables to be forced in the model
        forcedin_vars = [d1_col]

        #restart time for model fitting
        time_start = time.time()
        message += f'Fitting model for {d1}-{d2}: '

        #filter the df_matched based on d1
        d1_eligible_lst = [id_ for id_,vals in trajectory_ineligible_.items() if d1 not in vals]
        df_matched_d1 = df_matched[df_matched[id_col_].isin(d1_eligible_lst)]

        #remove groups with no outcomes
        group_sd = df_matched_d1.groupby(mathcing_id_col)[d2_col].std()
        group_id = group_sd[group_sd>0].index
        df_matched_d1 = df_matched_d1[df_matched_d1[mathcing_id_col].isin(group_id)]

        #create other diseases variable (already generated in the main function)
        if method in ['RPCN','PCN_PCA']:
            diseases_lst_ = [x for x in all_diseases_lst_ if x!=d1 and x!=d2]
            all_diseases_var = [str(x) for x in diseases_lst_]
        else:
            all_diseases_var = []
    
        df_analysis = pd.merge(df_matched_d1,phenotype_df_exposed_[[id_col_]+covariates_+all_diseases_var],
                               on=id_col_,how='left')
        del df_matched_d1
        gc.collect()
        df_analysis[d1_date_col] = df_analysis[id_col_].apply(lambda x: trajectory_eligible_withdate_[x].get(d1,pd.NaT))
        if enforce_time_interval==False:
            df_analysis[d1_col] = df_analysis.apply(lambda row: 1 if row[d1_date_col]<row[outcome_date_col] else 0, axis=1)
        else:
            df_analysis[d1_col] = df_analysis.apply(lambda row: 1 if row[d1_date_col]<row[outcome_date_col] and 
                                                    (row[outcome_date_col]-row[d1_date_col]).days<=max_interval_ else 0, axis=1)
            #exclude those with d1 date <= min interval date
            df_analysis['flag']  = df_analysis.apply(lambda row: 1 if (row[outcome_date_col]-row[d1_date_col]).days<=min_interval_ else 0, axis=1)
            #exclude those with flag=1 and d1_col=1
            df_analysis = df_analysis[(df_analysis['flag']==0) | (df_analysis[d1_col]==0)]
        
        #statistics
        n = len(df_analysis) #number of individuals in the matched case-control study
        N_d2 = len(df_analysis[df_analysis[d2_col]==1])
        N_nod2 = len(df_analysis[df_analysis[d2_col]==0])
        N_d1_withd2 = len(df_analysis[(df_analysis[d2_col]==1) & (df_analysis[d1_col]==1)])
        N_d1_nod2 = len(df_analysis[(df_analysis[d2_col]==0) & (df_analysis[d1_col]==1)])
    
        #result list
        result_lst = [d1,d2,f'{d1}-{d2}',N,n,f'{N_d1_withd2}/{N_d2}',f'{N_d1_nod2}/{N_nod2}']
        
        del_covariates = check_variance_vif_single(
            df_analysis,
            forcedin_vars,
            covariates_,
            vif_cutoff=100,
            group_col=mathcing_id_col
        )
        final_covariates = [x for x in covariates_ if x not in del_covariates]

        #simple method
        if method == 'CN':
            try:
                #remove constant variables before final model fitting
                forcedin_vars = [x for x in forcedin_vars if x!=constant_col]
                final_model_vars = forcedin_vars+final_covariates
                model = ConditionalLogit(np.asarray(df_analysis[d2_col],dtype=int),
                                        np.asarray(df_analysis[final_model_vars],dtype=float),
                                        groups=df_analysis[mathcing_id_col].values)
                result = model.fit(disp=False, method='bfgs')
                result_final = MyConditionalResultsWrapper(result)
                beta,se,p,aic = result_final.params[0], result_final.bse[0],result_final.pvalues[0],result_final.aic
                zvalue_dict = {var:z for var,z in zip(final_model_vars,result_final.tvalues)}
                result_lst += [method,f'fitted and delete the covariate(s): {del_covariates}',
                                f'{final_model_vars}',f'{zvalue_dict}',beta,se,p,aic]
                message += f'method={method}; successfully fitted; '
            except Exception as e:
                message += f'method={method}; error encountered: {e}; '
                result_lst += [method,str(e)]
            del df_analysis
            gc.collect()
        
        #partial correlation method
        elif method == 'RPCN':
            del_diseases_var = check_variance_vif_single(df_analysis,
                                                        forcedin_vars,all_diseases_var,
                                                        vif_cutoff='disease_covar',
                                                        group_col=mathcing_id_col)
            all_diseases_var = [x for x in all_diseases_var if x not in del_diseases_var]
            #generate alpha lst after removing these variables
            alpha_lst = np.array([0]*len(forcedin_vars) + [1]*len(all_diseases_var))
            #variables for first model for selecting disease variables
            model_1_vars = forcedin_vars+all_diseases_var
            if auto_penalty:
                alpha_lst = alpha_lst * scaling_factor
                try:
                    # model
                    model = Logit(np.asarray(df_analysis[d2_col], dtype=int),
                                np.asarray(df_analysis[model_1_vars], dtype=float))
                    #search within the defined range
                    final_best_alpha, final_disease_forcedin_vars = find_best_alpha_and_vars(model,alpha_range,alpha_lst,model_1_vars)
                    final_disease_vars = [x for x in final_disease_forcedin_vars if x not in forcedin_vars]
                    #fit the final model
                    #remove constant variables before final model fitting
                    final_disease_forcedin_vars = [x for x in final_disease_forcedin_vars if x!=constant_col]
                    final_model_vars = final_disease_forcedin_vars+final_covariates
                    model_final = ConditionalLogit(np.asarray(df_analysis[d2_col],dtype=int),
                                                np.asarray(df_analysis[final_model_vars],dtype=float),
                                                groups=df_analysis[mathcing_id_col].values)
                    result_final = model_final.fit(disp=False, method='bfgs')
                    result_final = MyConditionalResultsWrapper(result_final) #add aic property
                    beta,se,p,aic = result_final.params[0], result_final.bse[0],result_final.pvalues[0],result_final.aic
                    #get the z-value dictionary
                    z_value_dict = {var:z for var,z in zip(final_model_vars,result_final.tvalues)}
                    result_lst += [f'{method}_auto',f'fitted and delete the diseases variable(s): {del_diseases_var} and covariate(s): {del_covariates}',
                                f'{final_model_vars}',f'{z_value_dict}',final_best_alpha,beta,se,p,aic]
                    message += f'method={method}_auto (alpha={final_best_alpha}, number of other disease included as covariates: {len(final_disease_vars)}); successfully fitted; '
                except Exception as e:
                    result_lst += [f'{method}_auto', str(e)]
                    message += f'method={method}_auto; error encountered: {e}; '

            else:
                try:
                    #fit the initial model to get the non-zero disease list
                    model = Logit(np.asarray(df_analysis[d2_col],dtype=int),
                                np.asarray(df_analysis[model_1_vars],dtype=float)) #use unconditional model for selcting disease variables
                    result = model.fit_regularized(method='l1', alpha=alpha_lst*alpha_single, disp=False)
                    #get list of variables with none zero coef
                    non_zero_indices = np.nonzero(result.params != 0)[0]
                    final_disease_forcedin_vars = [model_1_vars[i] for i in non_zero_indices]
                    final_disease_vars = [x for x in final_disease_forcedin_vars if x not in forcedin_vars]
                    #fit the final conditional model
                    #remove constant variables before final model fitting
                    final_disease_forcedin_vars = [x for x in final_disease_forcedin_vars if x!=constant_col]
                    final_model_vars = final_disease_forcedin_vars+final_covariates
                    model_final = ConditionalLogit(np.asarray(df_analysis[d2_col],dtype=int),
                                                np.asarray(df_analysis[final_model_vars],dtype=float),
                                                groups=df_analysis[mathcing_id_col].values)
                    result_final = model_final.fit(disp=False, method='bfgs')
                    result_final = MyConditionalResultsWrapper(result_final) #add aic property
                    beta,se,p,aic = result_final.params[0],result_final.bse[0],result_final.pvalues[0],result_final.aic
                    #get the z-value dictionary
                    z_value_dict = {var:z for var,z in zip(final_model_vars,result_final.tvalues)}
                    result_lst += [f'{method}_auto',f'fitted and delete the diseases variable(s): {del_diseases_var} and covariate(s): {del_covariates}',
                                f'{final_model_vars}',f'{zvalue_dict}',alpha_single,beta,se,p,aic]
                    message += f'method={method}_fixed_alpha (alpha={alpha_single}, number of other disease included as covariates: {len(final_disease_vars)}); successfully fitted; '
                except Exception as e:
                    result_lst += [f'{method}_auto', str(e)]
                    message += f'method={method}_fixed_alpha (alpha={alpha_single}); error encountered: {e}; '
            del df_analysis
            gc.collect()

        elif method == 'PCN_PCA':
            from sklearn.decomposition import PCA # type: ignore
            try:
                #generate PC from other diseases variables
                pca = PCA(n_components=pca_number)
                #fit PCA model with unduplicated samples
                pca_model = pca.fit(np.asarray(df_analysis.drop_duplicates(subset=[id_col_])[all_diseases_var],dtype=int))
                #transform the original variables
                disease_vars_transformed = pca_model.transform(np.asarray(df_analysis[all_diseases_var],dtype=int))
                all_pca_vars = [f'PCA_{i}' for i in range(disease_vars_transformed.shape[1])]
                disease_vars_transformed = pd.DataFrame(disease_vars_transformed,columns=all_pca_vars)
                disease_vars_transformed.index = df_analysis.index
                df_analysis = pd.concat([df_analysis,disease_vars_transformed],axis=1)
                del disease_vars_transformed
                gc.collect()
                variance_explained = sum(pca.explained_variance_ratio_)
                #fit model with PCA covariates
                del_pca_var = check_variance_vif_single(df_analysis,
                                                        forcedin_vars,all_pca_vars,
                                                        vif_cutoff='pca_covar',
                                                        group_col=mathcing_id_col)
                final_pca_var = [x for x in all_pca_vars if x not in del_pca_var]
                #fit the final model
                #remove constant variables before final model fitting
                forcedin_vars = [x for x in forcedin_vars if x!=constant_col]
                final_model_vars = forcedin_vars+final_pca_var+final_covariates
                model_final = ConditionalLogit(np.asarray(df_analysis[d2_col],dtype=int),
                                            np.asarray(df_analysis[final_model_vars],dtype=float),
                                            groups=df_analysis[mathcing_id_col].values)
                result_final = model_final.fit(disp=False, method='bfgs')
                result_final = MyConditionalResultsWrapper(result_final) #add aic property
                beta,se,p,aic = result_final.params[0], result_final.bse[0], result_final.pvalues[0], result_final.aic
                z_value_dict = {var:z for var,z in zip(final_model_vars,result_final.tvalues)}
                result_lst += [f'{method}_n_components={pca_number}',f'fitted and delete the pca variable(s): {del_pca_var} and covariate(s): {del_covariates}',
                            f'{final_model_vars}',f'{z_value_dict}',variance_explained,beta,se,p,aic]
                message += f'method={method}_n_components={pca_number} (number of PC included as covariates: {len(final_pca_var)}, total variance explained by PC: {variance_explained:.3f}); successfully fitted; '
            except Exception as e:
                result_lst += [f'{method}_n_components={pca_number}',str(e)]
                message += f'method={method}_n_components={pca_number}; error encountered: {e}; '
            del df_analysis
            gc.collect()
        
        #print and return
        time_end = time.time()
        time_spend = time_end - time_start
        message += f'(elapsed {time_spend:.2f}s)\n'
        result_all.append(result_lst)
        #restart time for next loop
        time_start = time.time()

    write_log(log_file_,message)
    return result_all

def logistic_model_wrapper(d1_lst:list,d2:float,phenotype_df_exposed:pd.DataFrame,id_col,end_date_col,trajectory_ineligible:dict,
                            min_interval:int,max_interval:int,trajectory_eligible_withdate:dict,all_diagnosis_level:dict,covariates:list,
                            all_diseases_lst:list,matching_var_dict:dict,matching_n:int,max_n_cases:int,log_file:str,parameters:dict):
    """
    Wrapper for logistic_model that assigns default values to global variables if needed.

    Parameters
    ----------
    phenotype_df_exposed : pd.DataFrame, phenotypic data for exposed individuals only.
    trajectory_ineligible : dict, trajectory ineligible disease dictionary.
    id_col : str, id column
    end_date_col : str, date of end follow-up
    min_interval : int, minimum interval required for d1-d2 disease pair construction.
    max_interval : int, maximum interval allowed for d1-d2 disease pair construction
    trajectory_eligible_withdate : dict, trajectory eligible disease (with date) dictionary.
    all_diagnosis_level : list, list of all diagnosed phecodes, with phecode truncated to corresponding level
    covariates : list, list of covariates to be included in the model.
    all_diseases_lst : list, list of other diseases to be included.
    log_file : str, Path and prefix for the log file
    parameters : dict, other arguments, including method and the associated parameters.

    Returns
    -------
    Result list

    """
    #shared global data
    global phenotype_df_exposed_
    global id_col_
    global end_date_col_
    global trajectory_ineligible_
    global min_interval_
    global max_interval_
    global trajectory_eligible_withdate_
    global all_diagnosis_level_
    global covariates_
    global all_diseases_lst_
    global matching_var_dict_
    global matching_n_
    global max_n_cases_
    global log_file_
    global parameters_

    #assign values
    phenotype_df_exposed_ = phenotype_df_exposed
    id_col_ = id_col
    end_date_col_ = end_date_col
    trajectory_ineligible_ = trajectory_ineligible
    min_interval_ = min_interval
    max_interval_ = max_interval
    trajectory_eligible_withdate_ = trajectory_eligible_withdate
    all_diagnosis_level_ = all_diagnosis_level
    covariates_ = covariates
    all_diseases_lst_ = all_diseases_lst
    matching_var_dict_ = matching_var_dict
    matching_n_ = matching_n
    max_n_cases_ = max_n_cases
    log_file_ = log_file
    parameters_ = parameters
    # Call the original function
    return logistic_model((d1_lst, d2))

def init_worker(phenotype_df_exposed:pd.DataFrame,id_col,end_date_col,trajectory_ineligible:dict,
                min_interval:int,max_interval:int,trajectory_eligible_withdate:dict,all_diagnosis_level:dict,covariates:list,
                all_diseases_lst:list,matching_var_dict:dict,matching_n:int,max_n_cases:int,log_file:str,parameters:dict):
    """
    This function sets up the necessary global variables for a worker process in a multiprocessing environment.
    It assigns the provided parameters to global variables that can be accessed by logistic_model function in the worker process.

    Parameters:
    ----------
    phenotype_df_exposed : pd.DataFrame, phenotypic data for exposed individuals only.
    trajectory_ineligible : dict, trajectory ineligible disease dictionary.
    id_col : str, id column
    end_date_col : str, date of end follow-up
    min_interval : int, minimum interval required for d1-d2 disease pair construction.
    max_interval : int, maximum interval allowed for d1-d2 disease pair construction
    trajectory_eligible_withdate : dict, trajectory eligible disease (with date) dictionary.
    all_diagnosis_level : list, list of all diagnosed phecodes, with phecode truncated to corresponding level
    covariates : list, list of covariates to be included in the model.
    all_diseases_lst : list, list of other diseases to be included.
    log_file : str, Path and prefix for the log file
    parameters : dict, other arguments, including method and the associated parameters.

    Returns:
    ----------
    None

    """
    #shared global data
    global phenotype_df_exposed_
    global id_col_
    global end_date_col_
    global trajectory_ineligible_
    global min_interval_
    global max_interval_
    global trajectory_eligible_withdate_
    global all_diagnosis_level_
    global covariates_
    global all_diseases_lst_
    global matching_var_dict_
    global matching_n_
    global max_n_cases_
    global log_file_
    global parameters_
    #assign values
    phenotype_df_exposed_ = phenotype_df_exposed
    id_col_ = id_col
    end_date_col_ = end_date_col
    trajectory_ineligible_ = trajectory_ineligible
    min_interval_ = min_interval
    max_interval_ = max_interval
    trajectory_eligible_withdate_ = trajectory_eligible_withdate
    all_diagnosis_level_ = all_diagnosis_level
    covariates_ = covariates
    all_diseases_lst_ = all_diseases_lst
    matching_var_dict_ = matching_var_dict
    matching_n_ = matching_n
    max_n_cases_ = max_n_cases
    log_file_ = log_file
    parameters_ = parameters

def determine_best_range(aic_dict):
    """
    Determines the best range of alpha values based on the AIC values.

    Parameters:
    ----------
    aic_dict (dict): Dictionary where keys are alpha values and values are the corresponding AIC values.

    Returns:
    ----------
    tuple: The range (start, end) of alpha values that likely contains the best alpha.
    """
    # Sort the dictionary by alpha values to ensure the order
    sorted_aic = sorted(aic_dict.items())
    alpha_values = [item[0] for item in sorted_aic]
    aic_values = [item[1] for item in sorted_aic]
    # Find the index of the minimum AIC value
    min_aic_index = aic_values.index(min(aic_values))
    # Determine the range based on the surrounding alpha values
    if min_aic_index == 0:  # Best alpha is at the first value
        return (alpha_values[0], alpha_values[1])
    elif min_aic_index == len(aic_values) - 1:  # Best alpha is at the last value
        return (alpha_values[-2], alpha_values[-1])
    else:  # Best alpha is between two values
        return (alpha_values[min_aic_index - 1], alpha_values[min_aic_index + 1])

def matching_ids(df:pd.DataFrame,matching_var_dict:dict,matching_n:int,id_col,outcome_date_col:str,end_date_col:str,
                 save_outcome_col:str,save_outcome_date_col:str, save_matching_col:str):
    """
    Incidence density sampling matching.

    Parameters
    ----------
    df : pd.DataFrame, dataframe for matching
    matching_var_dict : dict, matching variable and matching
    matching_n : int, number of matched controls
    id_col : str, id column
    outcome_date_col : str, date of outcome
    end_date_col : str, date of end follow-up
    
    #the following columns are used for saving the results
    save_outcome_col : str, column name for outcome
    save_outcome_date_col : str, column name for outcome date
    save_matching_col : str, column name for matching

    Returns
    -------
    None.

    """
    result = []
    iter_ = 0
    case = df.loc[~df[outcome_date_col].isna()]
    for index in case.index:
        outcome_time = case.loc[index,outcome_date_col]
        sample = df[(df[end_date_col]>outcome_time) & ~(df[outcome_date_col]<=outcome_time)]
        for var in matching_var_dict:
            var_value = case.loc[index,var]
            if matching_var_dict[var] == 'exact':
                sample = sample[(sample[var]==var_value)]
            else:
                diff_range = matching_var_dict[var]
                sample = sample[(np.abs(sample[var]-var_value) <= diff_range)]
        try:
            sampleed = sample.sample(matching_n)
        except:
            sampleed = sample
        result += [[case.loc[index,id_col],1,outcome_time,iter_]]
        result += [[eid,0,outcome_time,iter_] for eid in sampleed[id_col].to_list()]
        iter_ += 1
    temp = pd.DataFrame(result,columns=[id_col,save_outcome_col,save_outcome_date_col,save_matching_col])
    return temp


#the added aic property can only be used for local comparison (fitted model for a same dataset but different set of variables)
class MyConditionalResultsWrapper(ConditionalResultsWrapper):
    @property
    def aic(self):
        # Example implementation of AIC
        # AIC formula: 2k - 2ln(L)
        # Where k is the number of parameters in the model, and L is the likelihood of the model
        k = len(self.params)  # Number of parameters
        log_likelihood = self.llf  # Log-likelihood of the model
        aic_value = -2*(log_likelihood - k)
        return aic_value
    








