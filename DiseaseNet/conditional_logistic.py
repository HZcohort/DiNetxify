# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 14:56:02 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

import pandas as pd
import numpy as np
from datetime import datetime
import statsmodels.api as sm
from statsmodels.discrete.conditional_models import ConditionalResultsWrapper
import time
from .utility import write_log

import warnings
warnings.filterwarnings('ignore')

def logistic_model(d1:float,d2:float,phenotype_df_exposed:pd.DataFrame,id_col,end_date_col,trajectory_eligible:dict,
                   trajectory_temporal:dict,trajectory_eligible_withdate:dict,history_level:dict,covariates:list,
                   all_diseases_lst:list,matching_var_dict:dict,matching_n:int,log_file:str,parameters:dict):
    """
    Fit a conditional LR model to verify the comorbidity association between a temporal disease pair.

    Parameters
    ----------
    d1 : str, phecode 1.
    d2 : str, phecode 2.
    phenotype_df_exposed : pd.DataFrame, phenotypic data for exposed individuals only.
    id_col : str, id column
    end_date_col : str, date of end follow-up
    trajectory_eligible : dict, trajectory eligible disease dictionary.
    trajectory_temporal : dict, temporal disease pair dictionary
    trajectory_eligible_withdate : dict, trajectory eligible disease (with date) dictionary.
    history_level :dict, history dictionary, with phecode truncated to corresponding level
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
    #method and parameters
    enforce_time_interval = parameters['enforce_time_interval']
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
    
    #matching
    phenotype_df_exposed['d2_date'] = phenotype_df_exposed[id_col].apply(lambda x: trajectory_eligible_withdate[x].get(d2,pd.NaT))
    df_matched = matching_ids(phenotype_df_exposed,matching_var_dict,matching_n,id_col,'d2_date',end_date_col)
    phenotype_df_exposed = pd.merge(df_matched,phenotype_df_exposed[[id_col,end_date_col,'d2_date']+covariates+list(matching_var_dict.keys())],on=id_col,how='left')
    phenotype_df_exposed['d1_date'] = phenotype_df_exposed[id_col].apply(lambda x: trajectory_eligible_withdate[x].get(d1,pd.NaT))
    phenotype_df_exposed['d1'] = phenotype_df_exposed.apply(lambda row: 1 if row['d1_date']<row['outcome_date'] else 0, axis=1)
    phenotype_df_exposed['constant'] = 1 #not used in condtional model fitting
    
    if enforce_time_interval==True:
        #for those with both d1 exposure and d2 outcome, further verify time interval requirement, as specified in disease pair construction
        d1_d2 = phenotype_df_exposed[(phenotype_df_exposed['d2']==1) & (phenotype_df_exposed['d1']==1)]
        d1_d2_eid = [x for x in d1_d2[id_col].values if (d1,d2) not in trajectory_temporal[x]]
        d1_d2_index = d1_d2[d1_d2[id_col].isin(d1_d2_eid)].index
        phenotype_df_exposed.loc[d1_d2_index,'d2'] = 0 #invalid cases
    
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

    #statistics
    n = len(phenotype_df_exposed) #number of individuals in the matched case-control study
    N_d2 = len(phenotype_df_exposed[phenotype_df_exposed['d2']==1])
    N_nod2 = len(phenotype_df_exposed[phenotype_df_exposed['d2']==0])
    N_d1_withd2 = len(phenotype_df_exposed[(phenotype_df_exposed['d2']==1) & (phenotype_df_exposed['d1']==1)])
    N_d1_nod2 = len(phenotype_df_exposed[(phenotype_df_exposed['d2']==0) & (phenotype_df_exposed['d1']==1)])
    
    #check var of covariates, remove these with var()==0
    for var in covariates:
        phenotype_df_exposed[var] = phenotype_df_exposed[var].astype(float)
        if phenotype_df_exposed.groupby(by='group_matching_ids')[var].var().mean() <= 0:
            covariates.remove(var)
            del phenotype_df_exposed[var]
    
    #result list
    result_lst = [d1,d2,f'{d1}-{d2}',N,n,f'{N_d1_withd2}/{N_d2}',f'{N_d1_nod2}/{N_nod2}']
    
    #time and message
    time_start = time.time()
    message = f'{d1} and {d2}: '
    
    #simple method
    if method == 'CN':
        try:
            model = sm.ConditionalLogit(np.asarray(phenotype_df_exposed['d2'],dtype=int),
                                        phenotype_df_exposed[['d1']+covariates].values,
                                        groups=phenotype_df_exposed['group_matching_ids'].values)
            result = model.fit(disp=False,method='bfgs')
            result = MyConditionalResultsWrapper(result)
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
                                 phenotype_df_exposed[model_1_vars].values) #use unconditional model for selcting disease variables
                
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
                model_final = sm.ConditionalLogit(np.asarray(phenotype_df_exposed['d2'],dtype=int),
                                                  phenotype_df_exposed[final_disease_vars+covariates].values,
                                                  groups=phenotype_df_exposed['group_matching_ids'].values)
                result_final = model_final.fit(disp=False,method='bfgs')
                result_final = MyConditionalResultsWrapper(result_final)
                beta,se,p,aic = result_final.params[0], result_final.bse[0],result_final.pvalues[0],result_final.aic
                z_value_dict = {var:z for var,z in zip(final_disease_vars+covariates,result_final.tvalues)}
                disease_z_value = {var:z_value_dict[var] for var in final_disease_vars[2::]} #z-value dictionary for other disease variables
                result_lst += [f'{method}_auto','fitted',f'{final_disease_vars[2::]}',f'{disease_z_value}',final_best_alpha,beta,se,p,aic]
                message += f'method={method}_auto (alpha={final_best_alpha}, number of other disease included as covariates: {len(final_disease_vars[2::])}); successfully fitted; '
            except Exception as e:
                result_lst += [f'{method}_auto',e]
                message += f'method={method}_auto; error encountered: {e}; '

        else:
            try:
                #fit the initial model to get the non-zero disease list
                model_1_vars = ['d1','constant']+all_diseases_var #only disease variables
                model = sm.Logit(np.asarray(phenotype_df_exposed['d2'],dtype=int),
                                 phenotype_df_exposed[model_1_vars].values)
                result = model.fit_regularized(method='l1', alpha=alpha_lst*alpha_single, disp=False)
                non_zero_indices = np.nonzero(result.params != 0)[0]
                final_disease_vars = [model_1_vars[i] for i in non_zero_indices]
                #fit the final conditional model
                model_final = sm.ConditionalLogit(np.asarray(phenotype_df_exposed['d2'],dtype=int),
                                                  phenotype_df_exposed[final_disease_vars+covariates].values,
                                                  groups=phenotype_df_exposed['group_matching_ids'].values)
                result_final = model_final.fit(disp=False,method='bfgs')
                result_final = MyConditionalResultsWrapper(result_final)
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
            phenotype_df_exposed_PCA = pd.concat([phenotype_df_exposed[['d1','d2']+covariates],
                                                  disease_vars_transformed],axis=1)
            variance_explained = sum(pca.explained_variance_ratio_)
            
            #fit model with PCA covariates
            model_final = sm.ConditionalLogit(np.asarray(phenotype_df_exposed_PCA['d2'],dtype=int),
                                   phenotype_df_exposed_PCA[['d1']+covariates+pca_cols].values,
                                   groups=phenotype_df_exposed['group_matching_ids'].values)
            result_final = model_final.fit(disp=False,method='bfgs')
            result_final = MyConditionalResultsWrapper(result_final)
            beta,se,p,aic = result_final.params[0], result_final.bse[0],result_final.pvalues[0],result_final.aic
            z_value_dict = {var:z for var,z in zip(['d1']+covariates+pca_cols,result_final.tvalues)}
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


def find_best_alpha_and_vars(model, best_range, alpha_lst, co_vars):
    """
    Function to find the best alpha for L1 regularization and the corresponding non-zero variables using an early stopping rule based on consecutive AIC increases.
    
    Parameters:
        model (statsmodels object): The statistical model to be fitted.
        best_range (tuple): A tuple (min_alpha, max_alpha) defining the range to explore.
        alpha_lst (float): The alpha multiplier applied during regularization.
        co_vars (list): List of variable names in the model.
    
    Returns:
        tuple: (final_best_alpha, final_disease_vars) where 'final_best_alpha' is the alpha value that
               results in the lowest AIC before AIC starts to increase consistently, and 'final_disease_vars'
               is a list of variables that are non-zero at this alpha level.
    """
    refined_alphas = np.linspace(best_range[0], best_range[1], num=best_range[1]-best_range[0]+1)
    refined_aic_dict = {}
    refined_vars_dict = {}
    min_aic = float('inf')
    counter = 0  # Counter to track the number of increases after a minimum
    thresold = 3 # early stop threshold

    for alpha in refined_alphas:
        try:
            result = model.fit_regularized(method='l1', alpha=alpha_lst*alpha, disp=False)
            non_zero_indices = np.nonzero(result.params != 0)[0]
            refined_vars_dict[alpha] = [co_vars[i] for i in non_zero_indices if co_vars[i]!='constant'] #constant should not be included for the final model
            refined_aic_dict[alpha] = result.aic
        except:
            # If the model fails to converge, set AIC to infinity
            refined_aic_dict[alpha] = float('inf')
            refined_vars_dict[alpha] = []
        
        # Check for AIC minimum and count increases
        if refined_aic_dict[alpha] < min_aic:
            min_aic = refined_aic_dict[alpha]
            counter = 0  # Reset counter on new minimum
        else:
            counter += 1  # Increment counter on increase
        
        # Break loop if AIC increases 5 times consecutively after a minimum
        if counter >= thresold:
            break

    final_best_alpha = min(refined_aic_dict, key=refined_aic_dict.get) * alpha_lst[-1]
    final_disease_vars = refined_vars_dict[final_best_alpha]
    if len(final_disease_vars) == 0:
        raise ValueError(f"All models failed when trying to find the best alpha for L1 regularization (stoped at {alpha*alpha_lst[-1]}).")
    
    return final_best_alpha, final_disease_vars

def matching_ids(df:pd.DataFrame,matching_var_dict:dict,matching_n:int,id_col,outcome_date_col:str,end_date_col:str):
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

    Returns
    -------
    None.

    """
    result = []
    iter_ = 0
    case = df.loc[~df[outcome_date_col].isna()]
    for index in case.index:
        outcome_time = case.loc[index,outcome_date_col]
        judge = (df[end_date_col]>outcome_time) & ~(df[outcome_date_col]<=outcome_time)
        for var in matching_var_dict:
            var_value = case.loc[index,var]
            if matching_var_dict[var] == 'exact':
                judge_ = (df[var]==var_value)
                judge = judge_ & judge
            else:
                diff_range = matching_var_dict[var]
                judge_ = (np.abs(df[var]-var_value) <= diff_range)
                judge = judge_ & judge
        try:
            sample = df[judge].sample(matching_n)
        except:
            sample = df[judge]
        result += [[case.loc[index,id_col],1,outcome_time,iter_]]
        result += [[eid,0,outcome_time,iter_] for eid in sample[id_col].to_list()]
        iter_ += 1
    temp = pd.DataFrame(result,columns=[id_col,'d2','outcome_date','group_matching_ids'])
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
    








