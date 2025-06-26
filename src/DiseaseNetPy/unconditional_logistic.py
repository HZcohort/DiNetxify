# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 02:13:51 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

import pandas as pd
import numpy as np
from statsmodels.discrete.discrete_model import Logit
import time
import gc
from .utility import write_log, find_best_alpha_and_vars, check_variance_vif_single

import warnings
warnings.filterwarnings('ignore')

def logistic_model(args):
    """
    Fit a LR model to verify the comorbidity association between a disease pair.

    Parameters
    ----------
    d1 : str, phecode 1.
    d2 : str, phecode 2.

    Global Variables:
    ----------
    phenotype_df_exposed : pd.DataFrame, phenotypic data for exposed individuals only.
    trajectory_ineligible : dict, trajectory ineligible disease dictionary.
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
    global trajectory_ineligible_
    global trajectory_eligible_withdate_
    global all_diagnosis_level_
    global covariates_
    global all_diseases_lst_
    global log_file_
    global parameters_

    d1, d2 = args
    #default columns
    d1_col = 'd1'
    d2_col = 'd2'
    constant_col = 'constant'
    forcedin_vars = [d1_col,constant_col]

    #method and parameters
    method = parameters_['method']
    if method == 'RPCN':
        #alpha_initial = [1, 10, 20, 30, 40, 50] #alpha starting value range if using auto_penalty
        auto_penalty = parameters_['auto_penalty']
        alpha_single = parameters_['alpha']
        alpha_range = parameters_['alpha_range']
        scaling_factor = parameters_['scaling_factor']
    elif method == 'PCN_PCA':
        pca_number = parameters_.get('explained_variance',parameters_.get('n_PC')) #retrive explained_variance first if given, otherwise use n_PC
    
    #filtering the dataframe first
    N = len(phenotype_df_exposed_)
    d1d2_eligible_lst = [id_ for id_,vals in trajectory_ineligible_.items() if d1 not in vals and d2 not in vals]
    df_analysis = phenotype_df_exposed_[phenotype_df_exposed_[id_col_].isin(d1d2_eligible_lst)]
    
    #create other diseases variable (already generated in the main function)
    if method in ['RPCN','PCN_PCA']:
        diseases_lst_ = [x for x in all_diseases_lst_ if x!=d1 and x!=d2]
        all_diseases_var = [str(x) for x in diseases_lst_]
    
    #d1 and d2 variable
    df_analysis[d1_col] = df_analysis[id_col_].apply(lambda x: 1 if d1 in trajectory_eligible_withdate_[x] else 0)
    df_analysis[d2_col] = df_analysis[id_col_].apply(lambda x: 1 if d2 in trajectory_eligible_withdate_[x] else 0)
    df_analysis[constant_col] = 1
    #statistics
    n = len(df_analysis) #number of individuals in the matched case-control study
    N_d2 = len(df_analysis[df_analysis[d2_col]==1])
    N_nod2 = len(df_analysis[df_analysis[d2_col]==0])
    N_d1_withd2 = len(df_analysis[(df_analysis[d2_col]==1) & (df_analysis[d1_col]==1)])
    N_d1_nod2 = len(df_analysis[(df_analysis[d2_col]==0) & (df_analysis[d1_col]==1)])

    #result list
    result_lst = [d1,d2,f'{d1}-{d2}',N,n,f'{N_d1_withd2}/{N_d2}',f'{N_d1_nod2}/{N_nod2}']
    
    #time and message
    time_start = time.time()
    message = f'{d1} and {d2}: '

    #check phnotypic covariates (applied to all methods)
    del_covariates = check_variance_vif_single(df_analysis,
                                               forcedin_vars,covariates_,
                                               vif_cutoff='phenotypic_covar')
    final_covariates= [x for x in covariates_ if x not in del_covariates]
    
    #simple method
    if method == 'CN':
        try:
            #list of final variables to be included
            final_model_vars = forcedin_vars+final_covariates
            model = Logit(np.asarray(df_analysis[d2_col], dtype=int),
                          np.asarray(df_analysis[final_model_vars], dtype=float))
            result_final = model.fit(disp=False, method='bfgs')
            beta,se,p,aic = result_final.params[0], result_final.bse[0],result_final.pvalues[0],result_final.aic
            zvalue_dict = {var:z for var,z in zip(final_model_vars,result_final.tvalues)}
            result_lst += [method,f'fitted and delete the covariate(s): {del_covariates}',
                            f'{final_model_vars}',f'{zvalue_dict}',beta,se,p,aic]
            message += f'method={method}; successfully fitted; '
        except Exception as e:
            message += f'method={method}; error encountered: {e}; '
            result_lst += [method,str(e)]
    
    #partial correlation method
    elif method == 'RPCN':
        #check disease variables
        del_diseases_var = check_variance_vif_single(df_analysis,
                                                     forcedin_vars,all_diseases_var,
                                                     vif_cutoff='disease_covar')
        all_diseases_var = [x for x in all_diseases_var if x not in del_diseases_var]
        #generate alpha lst after removing these variables
        alpha_lst = np.array([0]*len(forcedin_vars) + [1]*len(all_diseases_var))
        #variables for first model for selecting disease variables
        model_1_vars = forcedin_vars+all_diseases_var
        #method 1
        if auto_penalty:
            alpha_lst = alpha_lst * scaling_factor #scaling_factor only for auto_penalty
            try:
                #model
                model = Logit(np.asarray(df_analysis[d2_col], dtype=int),
                              np.asarray(df_analysis[model_1_vars], dtype=float))
                #search within the defined range, get the diseases variables that will be included plus forced-in variables
                final_best_alpha, final_disease_forcedin_vars = find_best_alpha_and_vars(model,alpha_range,alpha_lst,model_1_vars)
                final_disease_vars = [x for x in final_disease_forcedin_vars if x not in forcedin_vars]
                #fit the final model
                final_model_vars = final_disease_forcedin_vars+final_covariates
                #restrict the dataset to the final model variables
                df_analysis = df_analysis[final_model_vars+[d2_col]]
                model_final = Logit(np.asarray(df_analysis[d2_col],dtype=int),
                                    np.asarray(df_analysis[final_model_vars],dtype=float))
                result_final = model_final.fit(disp=False, method='bfgs')
                beta,se,p,aic = result_final.params[0],result_final.bse[0],result_final.pvalues[0],result_final.aic
                #get the z-value dictionary
                zvalue_dict = {var:z for var,z in zip(final_model_vars,result_final.tvalues)}
                result_lst += [f'{method}_auto',f'fitted and delete the diseases variable(s): {del_diseases_var} and covariate(s): {del_covariates}',
                               f'{final_model_vars}',f'{zvalue_dict}',final_best_alpha,beta,se,p,aic]
                message += f'method={method}_auto (alpha={final_best_alpha}, number of other disease included as covariates: {len(final_disease_vars)}); successfully fitted; '
            except Exception as e:
                result_lst += [f'{method}_auto', str(e)]
                message += f'method={method}_auto; error encountered: {e}; '

        else:
            try:
                #fit the initial model to get the non-zero disease list
                model = Logit(np.asarray(df_analysis[d2_col],dtype=int),
                              np.asarray(df_analysis[model_1_vars],dtype=float))
                result = model.fit_regularized(method='l1', alpha=alpha_lst*alpha_single, disp=False)
                #get list of variables with none zero coef
                non_zero_indices = np.nonzero(result.params != 0)[0]
                final_disease_forcedin_vars = [model_1_vars[i] for i in non_zero_indices]
                final_disease_vars = [x for x in final_disease_forcedin_vars if x not in forcedin_vars]
                #fit the final model
                final_model_vars = final_disease_forcedin_vars+final_covariates
                #restrict the dataset to the final model variables
                df_analysis = df_analysis[final_model_vars+[d2_col]]
                model_final = Logit(np.asarray(df_analysis[d2_col],dtype=int),
                                    np.asarray(df_analysis[final_model_vars]),dtype=float)
                result_final = model_final.fit(disp=False,method='bfgs')
                beta,se,p,aic = result_final.params[0], result_final.bse[0],result_final.pvalues[0],result_final.aic
                #get the z-value dictionary
                zvalue_dict = {var:z for var,z in zip(final_model_vars,result_final.tvalues)}
                result_lst += [f'{method}_fixed_alpha',f'fitted and delete the diseases variable(s): {del_diseases_var} and covariate(s): {del_covariates}',
                               f'{final_model_vars}',f'{zvalue_dict}',alpha_single,beta,se,p,aic]
                message += f'method={method}_fixed_alpha (alpha={alpha_single}, number of other disease included as covariates: {len(final_disease_vars)}); successfully fitted; '
            except Exception as e:
                result_lst += [f'{method}_fixed_alpha', str(e)]
                message += f'method={method}_fixed_alpha (alpha={alpha_single}); error encountered: {e}; '

    elif method == 'PCN_PCA':
        from sklearn.decomposition import PCA # type: ignore
        try:
            #generate PC from other diseases variables
            pca = PCA(n_components=pca_number)
            disease_vars_transformed = pca.fit_transform(np.asarray(df_analysis[all_diseases_var],dtype=int))
            #delete the original disease variables
            df_analysis = df_analysis.drop(all_diseases_var,axis=1)
            all_pca_vars = [f'PCA_{i}' for i in range(disease_vars_transformed.shape[1])]
            disease_vars_transformed = pd.DataFrame(disease_vars_transformed,columns=all_pca_vars)
            disease_vars_transformed.index = df_analysis.index
            df_analysis = pd.concat([df_analysis,disease_vars_transformed],axis=1)
            #delete the disease_vars_transformed
            del disease_vars_transformed
            gc.collect()
            variance_explained = sum(pca.explained_variance_ratio_)
            #check the VIF of PCA variables
            del_pca_var = check_variance_vif_single(df_analysis,
                                                    forcedin_vars,all_pca_vars,
                                                    vif_cutoff='pca_covar')
            final_pca_var = [x for x in all_pca_vars if x not in del_pca_var]
            #fit the final model
            final_model_vars = forcedin_vars+final_pca_var+final_covariates
            model_final = Logit(np.asarray(df_analysis[d2_col],dtype=int),
                                np.asarray(df_analysis[final_model_vars],dtype=float))
            result_final = model_final.fit(disp=False,method='bfgs')
            beta,se,p,aic = result_final.params[0], result_final.bse[0],result_final.pvalues[0],result_final.aic
            zvalue_dict = {var:z for var,z in zip(final_model_vars,result_final.tvalues)}
            result_lst += [f'{method}_n_components={pca_number}',f'fitted and delete the PC variable(s): {del_pca_var} and covariate(s): {del_covariates}',
                           f'{final_model_vars}',f'{zvalue_dict}',variance_explained,beta,se,p,aic]
            message += f'method={method}_n_components={pca_number} (number of PC included as covariates: {len(final_pca_var)}, total variance explained by PC: {variance_explained:.3f}); successfully fitted; '
        except Exception as e:
            result_lst += [f'{method}_n_components={pca_number}',str(e)]
            message += f'method={method}_n_components={pca_number}; error encountered: {e}; '
    
    #print and return
    time_end = time.time()
    time_spend = time_end - time_start
    message += f'(elapsed {time_spend:.2f}s)\n'
    write_log(log_file_,message)
    return result_lst

def logistic_model_wrapper(d1:float,d2:float,phenotype_df_exposed:pd.DataFrame,id_col,trajectory_ineligible:dict,
                            trajectory_eligible_withdate:dict,all_diagnosis_level:dict,covariates:list,
                            all_diseases_lst:list,log_file:str,parameters:dict):
    """
    Wrapper for logistic_model that assigns default values to global variables if needed.

    Parameters
    ----------
    d1 : float, phecode 1.
    d2 : float, phecode 2.
    phenotype_df_exposed : pd.DataFrame, phenotypic data for exposed individuals only.
    trajectory_ineligible : dict, trajectory ineligible disease dictionary.
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
    global trajectory_ineligible_
    global trajectory_eligible_withdate_
    global all_diagnosis_level_
    global covariates_
    global all_diseases_lst_
    global log_file_
    global parameters_
    #assign values
    phenotype_df_exposed_ = phenotype_df_exposed
    id_col_ = id_col
    trajectory_ineligible_ = trajectory_ineligible
    trajectory_eligible_withdate_ = trajectory_eligible_withdate
    all_diagnosis_level_ = all_diagnosis_level
    covariates_ = covariates
    all_diseases_lst_ = all_diseases_lst
    log_file_ = log_file
    parameters_ = parameters
    # Call the original function
    return logistic_model((d1, d2))

def init_worker(phenotype_df_exposed:pd.DataFrame,id_col,trajectory_ineligible:dict,
                trajectory_eligible_withdate:dict,all_diagnosis_level:dict,covariates:list,
                all_diseases_lst:list,log_file:str,parameters:dict):
    """
    This function sets up the necessary global variables for a worker process in a multiprocessing environment.
    It assigns the provided parameters to global variables that can be accessed by logistic_model function in the worker process.

    Parameters:
    ----------
    phenotype_df_exposed : pd.DataFrame, phenotypic data for exposed individuals only.
    trajectory_ineligible : dict, trajectory ineligible disease dictionary.
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
    global trajectory_ineligible_
    global trajectory_eligible_withdate_
    global all_diagnosis_level_
    global covariates_
    global all_diseases_lst_
    global log_file_
    global parameters_
    #assign values
    phenotype_df_exposed_ = phenotype_df_exposed
    id_col_ = id_col
    trajectory_ineligible_ = trajectory_ineligible
    trajectory_eligible_withdate_ = trajectory_eligible_withdate
    all_diagnosis_level_ = all_diagnosis_level
    covariates_ = covariates
    all_diseases_lst_ = all_diseases_lst
    log_file_ = log_file
    parameters_ = parameters
    
    
    