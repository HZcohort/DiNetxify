# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 02:13:51 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

import pandas as pd
import numpy as np
#from datetime import datetime
import statsmodels.api as sm
import time
from .utility import write_log,find_best_alpha_and_vars

import warnings
warnings.filterwarnings('ignore')

def logistic_model(d1:float,d2:float,phenotype_df_exposed:pd.DataFrame,id_col,trajectory_eligible:dict,
                   trajectory_eligible_withdate:dict,history_level:dict,covariates:list,
                   all_diseases_lst:list,log_file:str,parameters:dict):
    """
    Fit a LR model to verify the comorbidity association between a disease pair.

    Parameters
    ----------
    d1 : str, phecode 1.
    d2 : str, phecode 2.
    phenotype_df_exposed : pd.DataFrame, phenotypic data for exposed individuals only.
    trajectory_eligible : dict, trajectory eligible disease dictionary.
    trajectory_eligible_withdate : dict, trajectory eligible disease (with date) dictionary.
    history_level :dict, history dictionary, with phecode truncated to corresponding level
    covariates : list, list of covariates to be included in the model.
    all_diseases_lst : list, list of other diseases to be included.
    log_file : str, Path and prefix for the log file
    parameters : dict, other arguments, including method and the associated parameters.

    Returns
    -------
    Result list

    """
    #method and parameters
    method = parameters['method']
    if method == 'RPCN':
        #alpha_initial = [1, 10, 20, 30, 40, 50] #alpha starting value range if using auto_penalty
        auto_penalty = parameters['auto_penalty']
        alpha_single = parameters['alpha']
        alpha_range = parameters['alpha_range']
        scaling_factor = parameters['scaling_factor']
    elif method == 'PCN_PCA':
        pca_number = parameters.get('explained_variance',parameters.get('n_PC')) #retrive explained_variance first if given, otherwise use n_PC
    
    #filtering the dataframe first
    N = len(phenotype_df_exposed)
    d1d2_eligible_lst = [id_ for id_,vals in trajectory_eligible.items() if d1 in vals and d2 in vals]
    phenotype_df_exposed = phenotype_df_exposed[phenotype_df_exposed[id_col].isin(d1d2_eligible_lst)]
    
    #create other diseases variable
    if method in ['RPCN','PCN_PCA']:        
        all_diseases_lst = [x for x in all_diseases_lst if x!=d1 and x!=d2]
        all_diseases_var = []
        for disease in all_diseases_lst:
            phenotype_df_exposed[str(disease)] = phenotype_df_exposed[id_col].apply(lambda x: 1 if disease in history_level[x] or disease in trajectory_eligible_withdate[x] else 0)
            all_diseases_var.append(str(disease))
        if auto_penalty:
            alpha_lst = np.array([0]*(2) + [1]*len(all_diseases_var)) * scaling_factor #consider the scaling factor when using auto_penalty
        else:
            alpha_lst = np.array([0]*(2) + [1]*len(all_diseases_var))
    
    #d1 and d2 variable
    phenotype_df_exposed['d1'] = phenotype_df_exposed[id_col].apply(lambda x: 1 if d1 in trajectory_eligible_withdate[x] else 0)
    phenotype_df_exposed['d2'] = phenotype_df_exposed[id_col].apply(lambda x: 1 if d2 in trajectory_eligible_withdate[x] else 0)
    phenotype_df_exposed['constant'] = 1
    #statistics
    n = len(phenotype_df_exposed) #number of individuals in the sub-cohort
    N_d1 = len(phenotype_df_exposed[phenotype_df_exposed['d1']==1])
    N_d2_withd1 = len(phenotype_df_exposed[(phenotype_df_exposed['d2']==1) & (phenotype_df_exposed['d1']==1)])
    N_d2_nod1 = len(phenotype_df_exposed[(phenotype_df_exposed['d2']==1) & (phenotype_df_exposed['d1']==0)])
    
    #check var of covariates, remove these with var()==0
    for var in covariates:
        phenotype_df_exposed[var] = phenotype_df_exposed[var].astype(float)
        if phenotype_df_exposed[var].var() <= 0: #lowest var() allowed
            covariates.remove(var)
    
    #result list
    result_lst = [d1,d2,f'{d1}-{d2}',N,n,f'{N_d2_withd1}/{N_d1}',f'{N_d2_nod1}/{n-N_d1}']
    
    #time and message
    time_start = time.time()
    message = f'{d1} and {d2}: '
    
    #simple method
    if method == 'CN':
        try:
            model = sm.Logit(np.asarray(phenotype_df_exposed['d2'],dtype=int),
                             np.asarray(phenotype_df_exposed[['d1','constant']+covariates]),dtype=float)
            result = model.fit(disp=False,method='bfgs')
            beta,se,p,aic = result.params[0], result.bse[0],result.pvalues[0],result.aic
            result_lst += [method,'fitted',beta,se,p,aic]
            message += f'method={method}; successfully fitted; '
        except Exception as e:
            message += f'method={method}; error encountered: {e}; '
            result_lst += [method,e]
    
    #partial correlation method
    elif method == 'RPCN':
        if auto_penalty:
            try:
                #model
                model_1_vars = ['d1','constant']+all_diseases_var #only disease variables
                model = sm.Logit(np.asarray(phenotype_df_exposed['d2'],dtype=int),
                                 np.asarray(phenotype_df_exposed[model_1_vars],dtype=float))
                
                # Initial alphas to check
                """
                aic_dict = {}
                # Check AIC for initial alpha values
                for alpha in alpha_initial:
                    result = model.fit_regularized(method='l1', alpha=alpha_lst*alpha, disp=False)
                    aic_dict[alpha] = result.aic
                # Determine the range with the best AIC
                best_range = determine_best_range(aic_dict)
                """
                #search within the defined range
                final_best_alpha, final_disease_vars = find_best_alpha_and_vars(model,alpha_range,alpha_lst,model_1_vars)
                #fit the final model
                model_final = sm.Logit(np.asarray(phenotype_df_exposed['d2'],dtype=int),
                                       np.asarray(phenotype_df_exposed[final_disease_vars+covariates],dtype=float))
                result_final = model_final.fit(disp=False,method='bfgs')
                beta,se,p,aic = result_final.params[0], result_final.bse[0],result_final.pvalues[0],result_final.aic
                z_value_dict = {var:z for var,z in zip(final_disease_vars+covariates,result_final.tvalues)}
                disease_z_value = {var:z_value_dict[var] for var in final_disease_vars[2::]} #z-value dictionary for other disease variables
                result_lst += [f'{method}_auto','fitted',f'{final_disease_vars[2::]}',f'{disease_z_value}',final_best_alpha,beta,se,p,aic]
                message += f'method={method}_auto (alpha={final_best_alpha}, number of other disease included as covariates: {len(final_disease_vars[2::])}); successfully fitted; '
            except Exception as e:
                result_lst += [f'{method}_auto',e]
                message += f'method={method}_auto; error encountered: {e}; '
                phenotype_df_exposed[['d2']+final_disease_vars+covariates].to_csv(f'{log_file}_{d1}_{d2}.csv',index=False)

        else:
            try:
                #fit the initial model to get the non-zero disease list
                model_1_vars = ['d1','constant']+all_diseases_var #only disease variables
                model = sm.Logit(np.asarray(phenotype_df_exposed['d2'],dtype=int),
                                 np.asarray(phenotype_df_exposed[model_1_vars],dtype=float))
                result = model.fit_regularized(method='l1', alpha=alpha_lst*alpha_single, disp=False)
                non_zero_indices = np.nonzero(result.params != 0)[0]
                final_disease_vars = [model_1_vars[i] for i in non_zero_indices]
                
                #fit the final model
                model_final = sm.Logit(np.asarray(phenotype_df_exposed['d2'],dtype=int),
                                       np.asarray(phenotype_df_exposed[final_disease_vars+covariates]),dtype=float)
                result_final = model_final.fit(disp=False,method='bfgs')
                beta,se,p,aic = result_final.params[0], result_final.bse[0],result_final.pvalues[0],result_final.aic
                z_value_dict = {var:z for var,z in zip(final_disease_vars+covariates,result_final.tvalues)}
                disease_z_value = {var:z_value_dict[var] for var in final_disease_vars[2::]} #z-value dictionary for other disease variables
                result_lst += [f'{method}_fixed_alpha','fitted',f'{final_disease_vars[2::]}',f'{disease_z_value}',alpha_single,beta,se,p,aic]
                message += f'method={method}_fixed_alpha (alpha={alpha_single}, number of other disease included as covariates: {len(final_disease_vars[2::])}); successfully fitted; '
            except Exception as e:
                result_lst += [f'{method}_fixed_alpha',e]
                message += f'method={method}_fixed_alpha (alpha={alpha_single}); error encountered: {e}; '

        
    elif method == 'PCN_PCA':
        from sklearn.decomposition import PCA
        try:
            #generate PC from other diseases variables
            pca = PCA(n_components=pca_number)
            disease_vars_transformed = pca.fit_transform(np.asarray(phenotype_df_exposed[all_diseases_var],dtype=float))
            pca_cols = [f'PCA_{i}' for i in range(disease_vars_transformed.shape[1])]
            disease_vars_transformed = pd.DataFrame(disease_vars_transformed,columns=pca_cols)
            disease_vars_transformed.index = phenotype_df_exposed.index
            phenotype_df_exposed_PCA = pd.concat([phenotype_df_exposed[['d1','d2','constant']+covariates],
                                                  disease_vars_transformed],axis=1)
            variance_explained = sum(pca.explained_variance_ratio_)
            
            #fit model with PCA covariates
            model_final = sm.Logit(np.asarray(phenotype_df_exposed_PCA['d2'],dtype=int),
                                   np.asarray(phenotype_df_exposed_PCA[['d1','constant']+covariates+pca_cols],dtype=float))
            result_final = model_final.fit(disp=False,method='bfgs')
            beta,se,p,aic = result_final.params[0], result_final.bse[0],result_final.pvalues[0],result_final.aic
            z_value_dict = {var:z for var,z in zip(['d1','constant']+covariates+pca_cols,result_final.tvalues)}
            pca_z_value = {var:z_value_dict[var] for var in pca_cols} #z-value dictionary for other disease variables
            result_lst += [f'{method}_n_components={pca_number}','fitted',f'{pca_cols}',f'{pca_z_value}',variance_explained,beta,se,p,aic]
            message += f'method={method}_n_components={pca_number} (number of PC included as covariates: {len(pca_cols)}, total variance explained by PC: {variance_explained:.3f}); successfully fitted; '
        except Exception as e:
            result_lst += [f'{method}_n_components={pca_number}',e]
            message += f'method={method}_n_components={pca_number}; error encountered: {e}; '
    
    #print and return
    time_end = time.time()
    time_spend = time_end - time_start
    message += f'(elapsed {time_spend:.2f}s)\n'
    write_log(log_file,message)
    return result_lst
    

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



    
    
    
    