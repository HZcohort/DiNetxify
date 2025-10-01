# -*- coding: utf-8 -*-
"""
Created on Thu Dec 5 19:48:09 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

import pandas as pd
import numpy as np
import time
import random
from tqdm import tqdm
from .data_management import DiseaseNetworkData
from .utility import (
    log_file_detect,
    filter_phecodes,
    threshold_check,
    n_process_check,
    correction_method_check,
    states_p_adjust,
    check_kwargs_com_tra,
    covariates_check,
    matching_var_check
)

import warnings
warnings.filterwarnings('ignore')

def phewas(
    data: DiseaseNetworkData, 
    covariates: list=None, 
    proportion_threshold: float=None, 
    n_threshold: int=None, 
    n_process: int=1, 
    correction: str='bonferroni', 
    cutoff: float=0.05, 
    system_inc: list=None, 
    system_exl: list=None, 
    phecode_inc: list=None, 
    phecode_exl: list=None, 
    log_file: str=None, 
    lifelines_disable: bool=False
) -> pd.DataFrame:
    """
    Conducts Phenome-wide association studies (PheWAS) using the specified DiseaseNetworkData object.

    Parameters:
    ----------
    data : DiseaseNetworkData
        DiseaseNetworkData object.
    
    covariates : list, default=None
        List of phenotypic covariates to include in the model.
        By default, includes 'sex' and all covariates specified in the 'DiNetxify.DiseaseNetworkData.phenotype_data()' function.
        If you want to include the required variable sex as covariate, always use 'sex' rather than its original column name. 
        For other covariates you specified in the 'DiNetxify.DiseaseNetworkData.phenotype_data()' function, use their original column name.
        For matched cohort study, including a matching variable as covariate could cause issue of Singular Matrix in model fitting.
    
    proportion_threshold : float
        The minimum proportion of cases within the exposed group required for a phecode to be included in the PheWAS analysis.
        If the proportion of cases is below this threshold, the phecode is excluded from the analysis.
        proportion_threshold and n_threshold are mutually exclusive.
    
    n_threshold : int
        The minimum number of cases within the exposed group required for a phecode to be included in the PheWAS analysis.
        If the number of cases is below this threshold, the phecode is excluded from the analysis.
        n_threshold and proportion_threshold are mutually exclusive.      

    n_process : int, default=1
        Specifies the number of parallel processes to use for the analysis.
        Multiprocessing is enabled when `n_process` is set to a value greater than one.

    correction : str, default='bonferroni'
        Method for p-value correction from the statsmodels.stats.multitest.multipletests.
        Available methods are:
        none : no correction
        bonferroni : one-step correction
        sidak : one-step correction
        holm-sidak : step down method using Sidak adjustments
        holm : step-down method using Bonferroni adjustments
        simes-hochberg : step-up method (independent)
        hommel : closed method based on Simes tests (non-negative)
        fdr_bh : Benjamini/Hochberg (non-negative)
        fdr_by : Benjamini/Yekutieli (negative)
        fdr_tsbh : two stage fdr correction (non-negative)
        fdr_tsbky : two stage fdr correction (non-negative)
        See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
    
    cutoff : float, default=0.05
        The significance threshold for adjusted phewas p-values.
    
    system_inc : list, default=None
        List of phecode systems to include in the analysis. 
        system_inc and system_exl are mutually exclusive.
        List of eligible phecode systems: 
            circulatory system; congenital anomalies; dermatologic; digestive; 
            endocrine/metabolic; genitourinary; hematopoietic; infectious diseases; injuries & poisonings; 
            mental disorders; musculoskeletal; neoplasms; neurological; pregnancy complications; 
        respiratory; sense organs; symptoms; others.
    
    system_exl : list, default=None
        List of phecode systems to exclude from the analysis. 
        system_inc and system_exl are mutually exclusive.
        List of eligible phecode systems: 
            circulatory system; congenital anomalies; dermatologic; digestive; 
            endocrine/metabolic; genitourinary; hematopoietic; infectious diseases; injuries & poisonings; 
            mental disorders; musculoskeletal; neoplasms; neurological; pregnancy complications; 
        respiratory; sense organs; symptoms; others.
    
    phecode_inc : list, default=None
        Specific phecodes to include in the analysis. 
        phecode_inc and phecode_exl are mutually exclusive.
        
    phecode_exl : list, default=None
        Specific phecodes to exclude from the analysis. 
        phecode_inc and phecode_exl are mutually exclusive.
    
    log_file : str, default=None
        Path and prefix for the text file where log will be recorded.
        If None, the log will be written to the temporary files directory with file prefix of DiseaseNet_phewas_.
    
    lifelines_disable : bool, default=False
        Whether to disable the use of lifelines. 
        While lifelines generally require a longer fitting time, they are more robust to violations of model assumptions.

    Returns:
    ----------
    pd.DataFrame
        A pandas DataFrame object that contains the results of PheWAS analysis.

    """
    
    #data type check
    if not isinstance(data, DiseaseNetworkData):
        raise TypeError("The input 'data' must be a DiseaseNetworkData object.")
    
    #attribute check
    data_attrs = ['phenotype_df', 'diagnosis', 'n_diagnosis', 'history']
    for attr in data_attrs:
        if getattr(data, attr) is None:
            raise ValueError(f"Attribute '{attr}' is empty.")

    #retrieve phecode information
    phecode_info = data.phecode_info
 
    # check covariates
    covariates = covariates_check(covariates,data.get_attribute('phenotype_info'))
    
    #check lifelines_disable
    if not isinstance(lifelines_disable,bool):
        raise TypeError("The input 'lifelines_disable' must be a bool.")
    
    #check threshold
    n_exposed = data.get_attribute('phenotype_statistics')['n_exposed']
    n_threshold = threshold_check(proportion_threshold,n_threshold,n_exposed)
    
    #check p-value correction method and cutoff
    correction_method_check(correction, cutoff)
    
    #check inclusion and exclusion list
    phecode_lst_all = filter_phecodes(phecode_info, system_inc, system_exl, phecode_inc, phecode_exl)
    print(f'A total of {len(phecode_lst_all)} phecodes included in the PheWAS analysis.')
    
    #check log files
    log_file_final,message = log_file_detect(log_file,'phewas')
    print(message)

    #check number of process
    n_process,start_mehtod = n_process_check(n_process,'PheWAS')
    if n_process>1:
        import multiprocessing
        from .cox import cox_conditional,cox_unconditional,init_worker #use original function as main function and init_worker to initialize global variables
    else:
        from .cox import cox_unconditional_wrapper,cox_conditional_wrapper #use wrapper function as main function

    time_start = time.time()
    #list of phecode to run
    result_all = []
    if n_process == 1:
        for phecode in tqdm(phecode_lst_all, mininterval=15,smoothing=0):
            if data.study_design == 'matched cohort':
                result_all.append(cox_conditional_wrapper(phecode,data,covariates,n_threshold,log_file_final,lifelines_disable))
            else:
                result_all.append(cox_unconditional_wrapper(phecode,data,covariates,n_threshold,log_file_final,lifelines_disable))
    elif n_process > 1:
        with multiprocessing.get_context(start_mehtod).Pool(n_process, initializer=init_worker, initargs=(data,covariates,n_threshold,log_file_final,lifelines_disable)) as p:
            if data.study_design == 'matched cohort':
                result_all = list(
                    tqdm(
                        p.imap(cox_conditional, phecode_lst_all), 
                        total=len(phecode_lst_all),
                        mininterval=15,
                        smoothing=0
                    )
                )
            else:
                result_all = list(
                    tqdm(
                        p.imap(cox_unconditional, phecode_lst_all), 
                        total=len(phecode_lst_all),
                        mininterval=15,
                        smoothing=0
                    )
                )

    time_end = time.time()
    time_spent = (time_end - time_start)/60
    print(f'PheWAS analysis finished (elapsed {time_spent:.1f} mins)')
    
    #generate result dataframe
    max_columns = max([len(x) for x in result_all])
    columns = ['phecode','disease','system','sex','N_cases_exposed','describe','exposed_group','unexposed_group',
               'phewas_coef','phewas_se','phewas_p']
    columns_selected = columns[0:max_columns]
    phewas_df = pd.DataFrame(result_all, columns=columns_selected)
    
    if data.study_design == "exposed-only cohort":
        phewas_df["phewas_p_significance"] = phewas_df["N_cases_exposed"].apply(lambda x:True if x>=n_threshold else False)
        return phewas_df

    #p-value correction
    phewas_df = phewas_multipletests(phewas_df, correction=correction, cutoff=cutoff)
    return phewas_df


def phewas_multipletests(
    df:pd.DataFrame, 
    correction:str='bonferroni', 
    cutoff:float=0.05
) -> pd.DataFrame:
    """
    Adjusts PheWAS p-values for multiple comparisons using specified correction methods.

    Parameters:
    ----------
    df : pd.DataFrame
        DataFrame containing the results from the phewas function.

    correction : str, default='bonferroni'
        Method for p-value correction from the statsmodels.stats.multitest.multipletests.
        Available methods are:
        none : no correction
        bonferroni : one-step correction
        sidak : one-step correction
        holm-sidak : step down method using Sidak adjustments
        holm : step-down method using Bonferroni adjustments
        simes-hochberg : step-up method (independent)
        hommel : closed method based on Simes tests (non-negative)
        fdr_bh : Benjamini/Hochberg (non-negative)
        fdr_by : Benjamini/Yekutieli (negative)
        fdr_tsbh : two stage fdr correction (non-negative)
        fdr_tsbky : two stage fdr correction (non-negative)
        See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
    
    cutoff : float, default=0.05
        The significance threshold for adjusted phewas p-values.

    Returns:
    ----------
    pd.DataFrame
        A DataFrame that contains the PheWAS results with applied p-value corrections.

    """
    #data type check
    if not isinstance(df,pd.DataFrame):
        raise TypeError("The input 'df' must be a pandas DataFrame.")
    
    #check p-value correction method and cutoff
    correction_method_check(correction,cutoff)

    #multiple adjustment
    # loop to adjust
    for p_col,correction_,cutoff_ in [('phewas_p',correction,cutoff)]:
        #if p-value were not presented
        if p_col not in df.columns or len(df[~df[p_col].isna()])==0:
            print(f'No valid {p_col} found in the provided DataFrame, no p-value correction made.')
            continue
        else:
            df[p_col] = pd.to_numeric(df[p_col],errors='coerce')
            
            if correction_ == 'none':
                df[f'{p_col}_significance'] = df[p_col].apply(lambda x: True if x<=cutoff_ else False)
                df[f'{p_col}_adjusted'] = np.nan
            else:
                df = states_p_adjust(df,p_col,correction_,cutoff_,p_col,p_col)
            
    return df


def comorbidity_strength(
    data:DiseaseNetworkData, 
    proportion_threshold:float=None, 
    n_threshold:int=None, 
    n_process:int=1, 
    log_file:str=None, 
    correction_phi:str='bonferroni', 
    cutoff_phi:float=0.05, 
    correction_RR:str='bonferroni', 
    cutoff_RR:float=0.05
) -> pd.DataFrame:
    """
    Conducts comorbidity strength estimation among exposed individuals on all possible disease pairs using the specified DiseaseNetworkData object.
    For each disease pair, we evaluated its relative risk (RR) and phi-correlation as measurement of comorbidity strength.

    Parameters:
    ----------
    data : DiseaseNetworkData
        DiseaseNetworkData object.

    proportion_threshold : float
        The minimum proportion of individuals in the exposed group in which a disease pair must co-occur (temporal or non-temporal) to be included in the comorbidity strength estimation.
        If the proportion of co-occurrence is below this threshold, the disease pair is excluded from the analysis.
        proportion_threshold and n_threshold are mutually exclusive.
    
    n_threshold : int
        The minimum number of individuals in the exposed group in which a disease pair must co-occur (temporal or non-temporal) to be included in the comorbidity strength estimation.
        If the number of co-occurrences is below this threshold, the disease pair is excluded from the analysis.
        n_threshold and proportion_threshold are mutually exclusive.         

    n_process : int, default=1
        Specifies the number of parallel processes to use for the analysis.
        Multiprocessing is enabled when `n_process` is set to a value greater than one.

    correction_phi : str, default='bonferroni'
        Method for phi-correlation p-value correction from the statsmodels.stats.multitest.multipletests.
        Available methods are:
        none : no correction
        bonferroni : one-step correction
        sidak : one-step correction
        holm-sidak : step down method using Sidak adjustments
        holm : step-down method using Bonferroni adjustments
        simes-hochberg : step-up method (independent)
        hommel : closed method based on Simes tests (non-negative)
        fdr_bh : Benjamini/Hochberg (non-negative)
        fdr_by : Benjamini/Yekutieli (negative)
        fdr_tsbh : two stage fdr correction (non-negative)
        fdr_tsbky : two stage fdr correction (non-negative)
        See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
    
    cutoff_phi : float, default=0.05
        The significance threshold for adjusted phi-correlatio p-values.

    correction_RR : str, default='bonferroni'
        Method for RR p-value correction from the statsmodels.stats.multitest.multipletests.
        Available methods are:
        none : no correction
        bonferroni : one-step correction
        sidak : one-step correction
        holm-sidak : step down method using Sidak adjustments
        holm : step-down method using Bonferroni adjustments
        simes-hochberg : step-up method (independent)
        hommel : closed method based on Simes tests (non-negative)
        fdr_bh : Benjamini/Hochberg (non-negative)
        fdr_by : Benjamini/Yekutieli (negative)
        fdr_tsbh : two stage fdr correction (non-negative)
        fdr_tsbky : two stage fdr correction (non-negative)
        See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
    
    cutoff_RR : float, default=0.05
        The significance threshold for adjusted RR p-values.
    
    log_file : str, default=None
        Path and prefix for the text file where log will be recorded.
        If None, the log will be written to the temporary files directory with file prefix of DiseaseNet_com_strength_.

    Returns:
    ----------
    pd.DataFrame
        A pandas DataFrame object that contains the results of PheWAS analysis.

    """
    from itertools import combinations
    
    #data type check
    if not isinstance(data,DiseaseNetworkData):
        raise TypeError("The input 'data' must be a DiseaseNetworkData object.")
    
    #attribute check
    data_attrs = ['trajectory']
    for attr in data_attrs:
        if getattr(data, attr) is None:
            raise ValueError(f"Attribute '{attr}' is empty.")
    
    #retrieve phecode information and others
    phecode_info = data.phecode_info
    trajectory_dict = data.trajectory

    #check threshold
    n_exposed = data.get_attribute('phenotype_statistics')['n_exposed']
    n_threshold = threshold_check(proportion_threshold,n_threshold,n_exposed)
    
    #check p-value correction method and cutoff
    correction_method_check(correction_phi,cutoff_phi)
    correction_method_check(correction_RR,cutoff_RR)
    
    #check log files
    log_file_final,message = log_file_detect(log_file,'com_strength')
    print(message)

    #check number of process
    n_process,start_mehtod = n_process_check(n_process,'comorbidity_strength')
    if n_process>1:
        import multiprocessing
        from .com_strength import com_phi_rr, init_worker #use original function as main function and init_worker to initialize global variables
    else:
        from .com_strength import com_phi_rr_wrapper #use wrapper function as main function
    
    #get all significant phecodes
    phecodes_sig = data.get_attribute('significant_phecodes')
    #all possible disease pairs, highlight those with different sex specificity
    d1d2_pair_lst = []
    for c in combinations(phecodes_sig, 2):
        d1,d2 = c
        sex_d1,sex_d2 = phecode_info[d1]['sex'], phecode_info[d2]['sex']
        if {sex_d1,sex_d2}=={'Female','Male'}:
            d1d2_pair_lst.append((d1,d2,'Disease pair with different sex specificity'))
        else:
            d1d2_pair_lst.append((d1,d2,None))
    random.shuffle(d1d2_pair_lst)
    
    time_start = time.time()
    #list of phecode
    result_all = []
    if n_process == 1:
        for d1,d2,describe in tqdm(d1d2_pair_lst, mininterval=15,smoothing=0):
            result_all.append(com_phi_rr_wrapper(trajectory_dict,d1,d2,describe,n_threshold,log_file_final))
    elif n_process > 1:
        parameters_all = []
        for d1,d2,describe in d1d2_pair_lst:
            # parameters_all.append([trajectory_dict,d1,d2,describe,n_threshold,log_file_final])
            parameters_all.append((d1,d2,describe))
        with multiprocessing.get_context(start_mehtod).Pool(n_process, initializer=init_worker, initargs=(trajectory_dict,n_threshold,log_file_final)) as p:
            result_all = list(
                tqdm(
                    p.imap(com_phi_rr, parameters_all),
                    total=len(d1d2_pair_lst),
                    mininterval=15,
                    smoothing=0
                )
            )

    time_end = time.time()
    time_spent = (time_end - time_start)/60
    print(f'Comorbidity strength estimation finished (elapsed {time_spent:.1f} mins)')

    #generate result dataframe
    max_columns = max([len(x) for x in result_all])
    columns = ['phecode_d1','phecode_d2','name_disease_pair','N_exposed','n_total',
               'n_d1d2_diagnosis','n_d1_diagnosis','n_d2_diagnosis',
               'n_d1d2_nontemporal','n_d1d2_temporal','n_d2d1_temporal','n_d1d2_pair',
               'description','phi','phi_theta','phi_p','RR','RR_theta','RR_p']
    columns_selected = columns[0:max_columns]
    com_df = pd.DataFrame(result_all, columns=columns_selected)
    #annotate disease name and system
    for d in ['d1','d2']:
        com_df[f'disease_{d}'] = com_df[f'phecode_{d}'].apply(lambda x: phecode_info[x]['phenotype'])
        com_df[f'system_{d}'] = com_df[f'phecode_{d}'].apply(lambda x: phecode_info[x]['category'])
        com_df[f'sex_{d}'] = com_df[f'phecode_{d}'].apply(lambda x: phecode_info[x]['sex'])
    
    #p-value correction
    com_df = comorbidity_strength_multipletests(com_df, correction_phi=correction_phi,cutoff_phi=cutoff_phi,
                                                correction_RR=correction_RR,cutoff_RR=cutoff_RR)
    
    return com_df


def comorbidity_strength_multipletests(df:pd.DataFrame, correction_phi:str='bonferroni', cutoff_phi:float=0.05, 
                                       correction_RR:str='bonferroni', cutoff_RR:float=0.05) -> pd.DataFrame:
    """
    Adjusts comorbidity strength p-values (phi-correlation and RR) for multiple comparisons using specified correction methods.

    Parameters:
    ----------
    df : pd.DataFrame
        DataFrame containing the results from the comorbidity_strength function.

    correction_phi : str, default='bonferroni'
        Method for phi-correlation p-value correction from the statsmodels.stats.multitest.multipletests.
        Available methods are:
        none : no correction
        bonferroni : one-step correction
        sidak : one-step correction
        holm-sidak : step down method using Sidak adjustments
        holm : step-down method using Bonferroni adjustments
        simes-hochberg : step-up method (independent)
        hommel : closed method based on Simes tests (non-negative)
        fdr_bh : Benjamini/Hochberg (non-negative)
        fdr_by : Benjamini/Yekutieli (negative)
        fdr_tsbh : two stage fdr correction (non-negative)
        fdr_tsbky : two stage fdr correction (non-negative)
        See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
    
    cutoff_phi : float, default=0.05
        The significance threshold for adjusted phi-correlatio p-values.

    correction_RR : str, default='bonferroni'
        Method for RR p-value correction from the statsmodels.stats.multitest.multipletests.
        Available methods are:
        none : no correction
        bonferroni : one-step correction
        sidak : one-step correction
        holm-sidak : step down method using Sidak adjustments
        holm : step-down method using Bonferroni adjustments
        simes-hochberg : step-up method (independent)
        hommel : closed method based on Simes tests (non-negative)
        fdr_bh : Benjamini/Hochberg (non-negative)
        fdr_by : Benjamini/Yekutieli (negative)
        fdr_tsbh : two stage fdr correction (non-negative)
        fdr_tsbky : two stage fdr correction (non-negative)
        See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
    
    cutoff_RR : float, default=0.05
        The significance threshold for adjusted RR p-values.

    Returns:
    ----------
    pd.DataFrame
        A DataFrame that contains the comorbidity strength estimation results with applied p-value corrections.

    """
    #data type check
    if not isinstance(df,pd.DataFrame):
        raise TypeError("The input 'df' must be a pandas DataFrame.")
    
    #check p-value correction method and cutoff
    correction_method_check(correction_phi,cutoff_phi)
    correction_method_check(correction_RR,cutoff_RR)

    #multiple adjustment
    # loop to adjust
    for p_col,correction_,cutoff_ in [('phi_p',correction_phi,cutoff_phi),('RR_p',correction_RR,cutoff_RR)]:
        #if p-value were not presented
        if p_col not in df.columns or len(df[~df[p_col].isna()])==0:
            print(f'No valid {p_col} found in the provided DataFrame, no p-value correction made.')
            return df
        else:
            df[p_col] = pd.to_numeric(df[p_col],errors='coerce')
            
            if correction_ == 'none':
                df[f'{p_col}_significance'] = df[p_col].apply(lambda x: True if x<=cutoff_ else False)
                df[f'{p_col}_adjusted'] = np.nan
            else:
                df = states_p_adjust(df,p_col,correction_,cutoff_,p_col,p_col)
            
    return df


def binomial_test(
    data:DiseaseNetworkData, 
    comorbidity_strength_result:pd.DataFrame, 
    comorbidity_network_result: pd.DataFrame=None,
    n_process:int=1, 
    log_file:str=None, 
    correction:str='bonferroni', 
    cutoff:float=0.05, 
    enforce_temporal_order:bool=False, 
    **kwargs
) -> pd.DataFrame:
    """
    Conduct binomial test for disease pairs with significant comorbidity stregnth to select those with significant temporal orders (i.e., D1 -> D2).

    Parameters:
    ----------
    data : DiseaseNetworkData
        DiseaseNetworkData object.
    
    comorbidity_strength_result : pd.DataFrame
        DataFrame containing comorbidity strength analysis results produced by the 'DiNetxify.comorbidity_strength' function.
    
    comorbidity_network_result : pd.DataFrame, default=None
        DataFrame containing comorbidity network analysis results produced by the 'DiNetxify.comorbidity_network' function.
        When provided, the binomial test is limited to disease pairs deemed significant in the comorbidity network analysis.

    n_process : int, default=1
        Multiprocessing is disabled for this analysis.

    correction : str, default='bonferroni'
        Method for binomial test p-value correction from the statsmodels.stats.multitest.multipletests.
        Available methods are:
        none : no correction
        bonferroni : one-step correction
        sidak : one-step correction
        holm-sidak : step down method using Sidak adjustments
        holm : step-down method using Bonferroni adjustments
        simes-hochberg : step-up method (independent)
        hommel : closed method based on Simes tests (non-negative)
        fdr_bh : Benjamini/Hochberg (non-negative)
        fdr_by : Benjamini/Yekutieli (negative)
        fdr_tsbh : two stage fdr correction (non-negative)
        fdr_tsbky : two stage fdr correction (non-negative)
        See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
    
    cutoff : float, default=0.05
        The significance threshold for adjusted binomial p-values.
    
    log_file : str, default=None
        Path and prefix for the text file where log will be recorded.
        If None, the log will be written to the temporary files directory with file prefix of DiseaseNet_binomial_test_.
    
    enforce_temporal_order : bool, default=False
        If True, exclude individuals with non-temporal D1-D2 pair when performing the test.
        If False, include all individuals, including those with non-temporal D1-D2 pair.
    
    **kwargs
        Additional keyword argument to define the required columns in 'comorbidity_strength_result':
            phecode_d1_col : str, default='phecode_d1'
                Name of the column in 'comorbidity_strength_result' and 'comorbidity_network_result' that specifies the phecode identifiers for disease 1 of the disease pair.
            phecode_d2_col : str, default='phecode_d2'
                Name of the column in 'comorbidity_strength_result' and 'comorbidity_network_result' that specifies the phecode identifiers for disease 2 of the disease pair.
            n_nontemporal_col : str, default='n_d1d2_nontemporal'
                Name of the column in 'comorbidity_strength_result' that specifies the number of individuals with non-temporal d1-d2 disease pair
            n_temporal_d1d2_col : str, default='n_d1d2_temporal'
                Name of the column in 'comorbidity_strength_result' that specifies the number of individuals with temporal d1->d2 disease pair.
            n_temporal_d2d1_col : str, default='n_d2d1_temporal'
                Name of the column in 'comorbidity_strength_result' that specifies the number of individuals with temporal d2->d1 disease pair.
            significance_phi_col : str, default='phi_p_significance'
                Name of the column in 'comorbidity_strength_result' that indicates the significance of phi-correlation for each disease pair.
            significance_RR_col : str, default='RR_p_significance'
                Name of the column in 'comorbidity_strength_result' that indicates the significance of RR for each disease pair.
            significance_coef_col : str, default='comorbidity_p_significance'
                Name of the column in 'comorbidity_network_result' that indicates the significance of comorbidity network analysis for each disease pair.

    Returns:
    ----------
    pd.DataFrame
        A pandas DataFrame object that contains the results of binomial test.

    """
    from .binomial import binomial
    
    #data type check
    if not isinstance(data,DiseaseNetworkData):
        raise TypeError("The input 'data' must be a DiseaseNetworkData object.")
    
    #attribute check
    data_attrs = ['trajectory']
    for attr in data_attrs:
        if getattr(data, attr) is None:
            raise ValueError(f"Attribute '{attr}' is empty.")
    
    #check comorbidity strength estimation result
    if not isinstance(comorbidity_strength_result,pd.DataFrame):
        raise TypeError("The provided input 'comorbidity_strength_result' must be a pandas DataFrame.")
    
    #check column existence for comorbidity strength result
    phecode_d1_col = kwargs.get('phecode_d1_col', 'phecode_d1')
    phecode_d2_col = kwargs.get('phecode_d2_col', 'phecode_d2')
    n_nontemporal_col = kwargs.get('n_nontemporal_col', 'n_d1d2_nontemporal')
    n_temporal_d1d2_col = kwargs.get('n_temporal_d1d2_col', 'n_d1d2_temporal')
    n_temporal_d2d1_col = kwargs.get('n_temporal_d2d1_col', 'n_d2d1_temporal')
    significance_phi_col = kwargs.get('significance_phi_col', 'phi_p_significance')
    significance_RR_col = kwargs.get('significance_coef_col', 'RR_p_significance')
    for col in [phecode_d1_col, phecode_d2_col, significance_phi_col, significance_RR_col, 
                n_nontemporal_col, n_temporal_d1d2_col, n_temporal_d2d1_col]:
        if col not in comorbidity_strength_result.columns:
            raise ValueError(f"Column {col} not in 'comorbidity_strength_result' DataFrame.")
    
    #check comorbidity network result
    if comorbidity_network_result is not None:
        if not isinstance(comorbidity_network_result,pd.DataFrame):
            raise TypeError("The provided input 'comorbidity_network_result' must be a pandas DataFrame.")
        significance_coef_col = kwargs.get('significance_coef_col', 'comorbidity_p_significance')
        if significance_coef_col not in comorbidity_network_result.columns:
            raise ValueError(f"Column {significance_coef_col} not in 'comorbidity_network_result' DataFrame.")

    #check number of process
    n_process, start_mehtod = n_process_check(n_process,'binomial_test')
    if n_process>1:
        raise ValueError('Multiprocessing has been disabled for this analysis.')

    #check p-value correction method and cutoff
    correction_method_check(correction,cutoff)
    
    #check log files
    log_file_final,message = log_file_detect(log_file,'binomial_test')
    print(message)
    
    #check bool type
    if not isinstance(enforce_temporal_order,bool):
        raise TypeError("The provided input 'enforce_temporal_order' must be a bool.")
    
    #retrieve phecode information
    phecode_info = data.phecode_info
    
    #get all disease pairs with significant comorbidity strength
    comorbidity_sig = comorbidity_strength_result[(comorbidity_strength_result[significance_phi_col]==True) & 
                                                  (comorbidity_strength_result[significance_RR_col]==True)]
    if len(comorbidity_sig) == 0:
        raise ValueError("No disease pair remained after filtering on significance of phi-correlation and RR.")
    
    #filter disease pairs with significant comorbidity network analysis
    if comorbidity_network_result is not None:
        comorbidity_network_sig = comorbidity_network_result[comorbidity_network_result[significance_coef_col]==True]
        if len(comorbidity_network_sig) == 0:
            raise ValueError("No disease pair remained after filtering on results of comorbidity network analysis.")
        else:
            # filter disease pairs with significant comorbidity network analysis
            set_comorbidity_sig = set(comorbidity_sig[[phecode_d1_col,phecode_d2_col]].apply(frozenset, axis=1))
            set_comorbidity_network_sig = set(comorbidity_network_sig[[phecode_d1_col,phecode_d2_col]].apply(frozenset, axis=1))
            # raise error if any disease pairs presented 'comorbidity_network_sig' are not in 'comorbidity_sig'
            # and print out the disease pairs that are not in 'comorbidity_sig'
            if not set_comorbidity_network_sig.issubset(set_comorbidity_sig):
                diff = set_comorbidity_network_sig.difference(set_comorbidity_sig)
                raise ValueError(f"Disease pairs {diff} in 'comorbidity_network_result' are not in 'comorbidity_strength_result'.")
            comorbidity_sig = comorbidity_sig[comorbidity_sig[[phecode_d1_col,phecode_d2_col]].apply(frozenset, axis=1).isin(set_comorbidity_network_sig)]

    time_start = time.time()
    #list of disease pair
    result_all = []

    for d1,d2,n_com,n_d1d2,n_d2d1 in tqdm(comorbidity_sig[[phecode_d1_col,phecode_d2_col,n_nontemporal_col,n_temporal_d1d2_col,n_temporal_d2d1_col]].values,
        mininterval=15,smoothing=0):
        result_all.append(binomial(d1,d2,n_com,n_d1d2,n_d2d1,enforce_temporal_order,log_file_final))
    time_end = time.time()
    time_spent = (time_end - time_start)/60
    print(f'Binomial test finished (elapsed {time_spent:.1f} mins)')
    
    #generate result dataframe
    max_columns = max([len(x) for x in result_all])
    columns = ['phecode_d1','phecode_d2','name_disease_pair','n_d1d2_nontemporal','n_d1d2_temporal','n_d2d1_temporal',
               'binomial_p','binomial_proportion','binomial_proportion_ci']
    columns_selected = columns[0:max_columns]
    bino_df = pd.DataFrame(result_all, columns=columns_selected)
    #annotate disease name and system
    for d in ['d1','d2']:
        bino_df[f'disease_{d}'] = bino_df[f'phecode_{d}'].apply(lambda x: phecode_info[x]['phenotype'])
        bino_df[f'system_{d}'] = bino_df[f'phecode_{d}'].apply(lambda x: phecode_info[x]['category'])
        bino_df[f'sex_{d}'] = bino_df[f'phecode_{d}'].apply(lambda x: phecode_info[x]['sex'])
    
    #p-value correction
    bino_df = binomial_multipletests(bino_df, correction=correction,cutoff=cutoff)
    
    return bino_df


def binomial_multipletests(
    df:pd.DataFrame,
    correction:str='bonferroni', 
    cutoff:float=0.05
) -> pd.DataFrame:
    """
    Adjusts binomial p-values for multiple comparisons using specified correction methods.

    Parameters:
    ----------
    df : pd.DataFrame
        DataFrame containing the results from the comorbidity_strength function.

    correction : str, default='bonferroni'
        Method for binomial p-value correction from the statsmodels.stats.multitest.multipletests.
        Available methods are:
        none : no correction
        bonferroni : one-step correction
        sidak : one-step correction
        holm-sidak : step down method using Sidak adjustments
        holm : step-down method using Bonferroni adjustments
        simes-hochberg : step-up method (independent)
        hommel : closed method based on Simes tests (non-negative)
        fdr_bh : Benjamini/Hochberg (non-negative)
        fdr_by : Benjamini/Yekutieli (negative)
        fdr_tsbh : two stage fdr correction (non-negative)
        fdr_tsbky : two stage fdr correction (non-negative)
        See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
    
    cutoff : float, default=0.05
        The significance threshold for adjusted binomial p-values.

    Returns:
    ----------
    pd.DataFrame
        A DataFrame that contains the binomial test results with applied p-value corrections.

    """
    #data type check
    if not isinstance(df,pd.DataFrame):
        raise TypeError("The input 'df' must be a pandas DataFrame.")
    
    #check p-value correction method and cutoff
    correction_method_check(correction,cutoff)

    #multiple adjustment
    # loop to adjust
    for p_col,correction_,cutoff_ in [('binomial_p',correction,cutoff)]:
        #if p-value were not presented
        if p_col not in df.columns or len(df[~df[p_col].isna()])==0:
            print(f'No valid {p_col} found in the provided DataFrame, no p-value correction made.')
            return df
        else:
            df[p_col] = pd.to_numeric(df[p_col],errors='coerce')
            
            if correction_ == 'none':
                df[f'{p_col}_significance'] = df[p_col].apply(lambda x: True if x<=cutoff_ else False)
                df[f'{p_col}_adjusted'] = np.nan
            else:
                df = states_p_adjust(df,p_col,correction_,cutoff_,p_col,p_col)
    return df


def comorbidity_network(
    data:DiseaseNetworkData,
    comorbidity_strength_result:pd.DataFrame, 
    binomial_test_result:pd.DataFrame=None, 
    method:str='RPCN', 
    covariates:list=None, 
    n_process:int=1, 
    log_file:str=None, 
    correction:str='bonferroni', 
    cutoff:float=0.05, 
    **kwargs
) -> pd.DataFrame:
    """
    Perform comorbidity network analysis on disease pairs with significant comorbidity strength to identify pairs with confirmed comorbidity associations.

    Depending on the selected 'method', the function applies different statistical models to estimate the correlations for each disease pair:
    - **RPCN (Regularized Partial Correlation Network):**
        Utilizes L1-regularized logistic regression to estimate partial correlations for each disease pair.
        Includes both phenotypic variables and other diseases present in the network as covariates.
        The L1 regularization term selects important confounding disease variables.
        After variable selection, a standard logistic regression model is refitted to accurately estimate partial correlations, adjusting for phenotypic variables and the selected confounding diseases.
    - **PCN_PCA (Partial Correlation Network with PCA):**
        Applies a standard logistic regression model for each disease pair.
        Adjusts for phenotypic variables and the top principal components (PCs) of other diseases in the network to estimate partial correlations.
    - **CN (Correlation Network):**
        Uses a standard logistic regression to estimate simple correlations for each disease pair. Adjusts only for phenotypic variables.

    Parameters
    ----------
    data : DiseaseNetworkData
        DiseaseNetworkData object.

    comorbidity_strength_result : pd.DataFrame
        DataFrame containing comorbidity strength analysis results produced by the 'DiNetxify.comorbidity_strength' function.
    
    binomial_test_result : pd.DataFrame, default=None
        DataFrame containing binomial test analysis results produced by the 'DiNetxify.binomial_test' function.

    method : str, default='RPCN'
        Specifies the comorbidity network analysis method to use. Choices are:
        - 'RPCN: Regularized Partial Correlation Network.
        - 'PCN_PCA: Partial Correlation Network with PCA.
        - 'CN': Correlation Network.
        
        **Additional Options for RPCN:**
        - 'alpha' : non-negative scalar
            The weight multiplying the l1 penalty term for other diseases covariates. 
            Ignored if 'auto_penalty' is enabled.
        - 'auto_penalty' : bool, default=True
            If 'True', automatically determine the optimal 'alpha' based on model AIC value.
        - 'alpha_range' : tuple, default=(1,15)
            When 'auto_penalty' is True, search the optimal 'alpha' in this range.
        - 'scaling_factor' : positive scalar, default=1
            The scaling factor for the alpha when 'auto_penalty' is True.
        
        **Additional Options for PCN_PCA:**
        - 'n_PC' : int, default=5
            Fixed number of principal components to include in each model.
        - 'explained_variance' : float
            Determines the number of principal components based on the cumulative explained variance. 
            Overrides 'n_PC' if specified.

    covariates : list, default=None
        List of phenotypic covariates to include in the model.
        By default, includes ['sex'] and all covariates specified in the 'DiNetxify.DiseaseNetworkData.phenotype_data()' function.
        To include the required variable sex as a covariate, always use 'sex' instead of its original column name.
        For other covariates specified in the 'DiNetxify.DiseaseNetworkData.phenotype_data()' function, use their original column names.

    n_process : int, default=1
        Specifies the number of parallel processes to use for the analysis.
        Multiprocessing is enabled when `n_process` is set to a value greater than one.

    correction : str, default='bonferroni'
        Method for comorbidity network analysis p-value correction from the statsmodels.stats.multitest.multipletests.
        Available methods are:
        none : no correction
        bonferroni : one-step correction
        sidak : one-step correction
        holm-sidak : step down method using Sidak adjustments
        holm : step-down method using Bonferroni adjustments
        simes-hochberg : step-up method (independent)
        hommel : closed method based on Simes tests (non-negative)
        fdr_bh : Benjamini/Hochberg (non-negative)
        fdr_by : Benjamini/Yekutieli (negative)
        fdr_tsbh : two stage fdr correction (non-negative)
        fdr_tsbky : two stage fdr correction (non-negative)
        See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
    
    cutoff : float, default=0.05
        The significance threshold for adjusted comorbidity network analysis p-values.
    
    log_file : str, default=None
        Path and prefix for the text file where log will be recorded.
        If None, the log will be written to the temporary files directory with file prefix of DiseaseNet_comorbidity_network_.
    
    **kwargs
        Additional keyword argument to define the required columns in 'comorbidity_strength_result' and 'binomial_test_result':
            phecode_d1_col : str, default='phecode_d1'
                Name of the column in 'comorbidity_strength_result' and 'binomial_test_result' that specifies the phecode identifiers for disease 1 of the disease pair.
            phecode_d2_col : str, default='phecode_d2'
                Name of the column in 'comorbidity_strength_result' and 'binomial_test_result' that specifies the phecode identifiers for disease 2 of the disease pair.
            significance_phi_col : str, default='phi_p_significance'
                Name of the column in 'comorbidity_strength_result' that indicates the significance of phi-correlation for each disease pair.
            significance_RR_col : str, default='RR_p_significance'
                Name of the column in 'comorbidity_strength_result' that indicates the significance of RR for each disease pair.
            significance_binomial_col : str default='binomial_p_significance'
                Name of the column in 'binomial_test_result' that indicates the significance of binomial test for each disease pair.
        
        RPCN Method Parameters:
            alpha : non-negative scalar
                The weight multiplying the l1 penalty term for other diseases covariates. 
                Ignored if 'auto_penalty' is enabled.
            auto_penalty : bool, default=True
                If 'True', automatically determines the best 'alpha' based on model AIC value.
            alpha_range : tuple, default=(1,15)
                When 'auto_penalty' is True, search the optimal 'alpha' in this range.
            scaling_factor : positive scalar, default=1
                The scaling factor for the alpha when 'auto_penalty' is True.

        PCN_PCA Method Parameters:
            n_PC : int, default=5
                Fixed number of principal components to include in each model.
            explained_variance : float
                Cumulative explained variance threshold to determine the number of principal components. 
                Overrides 'n_PC' if specified.
                 
    Returns:
    ----------
    pd.DataFrame
        A pandas DataFrame object that contains the results of comorbidity network analysis.

    """
    
    #data type check
    if not isinstance(data,DiseaseNetworkData):
        raise TypeError("The input 'data' must be a DiseaseNetworkData object.")
    
    #attribute check
    data_attrs = ['history','diagnosis', 'n_diagnosis','trajectory','phenotype_df']
    for attr in data_attrs:
        if getattr(data, attr) is None:
            raise ValueError(f"Attribute '{attr}' is empty.")
    
    #check comorbidity strength estimation result
    if not isinstance(comorbidity_strength_result, pd.DataFrame):
        raise TypeError("The provided input 'comorbidity_strength_result' must be a pandas DataFrame.")
    comorbidity_strength_result_cols = comorbidity_strength_result.columns
    
    #check binomial test result
    if binomial_test_result is not None:
        if not isinstance(binomial_test_result,pd.DataFrame):
            raise TypeError("The provided input 'binomial_test_result' must be a pandas DataFrame.")
    binomial_test_result_cols = binomial_test_result.columns if binomial_test_result is not None else None

    # Allowed methods
    allowed_methods = {'RPCN', 'PCN_PCA', 'CN'}
    if method not in allowed_methods:
        raise ValueError(f"Invalid method '{method}'. Allowed methods are: {allowed_methods}.")

    # check covariates
    covariates = covariates_check(covariates,data.get_attribute('phenotype_info'))
    
    #check number of process
    n_process,start_mehtod = n_process_check(n_process,'comorbidity_network')
    if n_process>1:
        import multiprocessing
        from .unconditional_logistic import logistic_model,init_worker #use original function as main function and init_worker to initialize global variables
    else:
        from .unconditional_logistic import logistic_model_wrapper #use wrapper function as main function

    #check p-value correction method and cutoff
    correction_method_check(correction,cutoff)
    
    #check log files
    log_file_final,message = log_file_detect(log_file,'comorbidity_network_analysis')
    print(message)
    
    #check **kwargs and get parameters
    parameter_dict,phecode_d1_col,phecode_d2_col,significance_phi_col,significance_RR_col,significance_binomial_col = check_kwargs_com_tra(method,comorbidity_strength_result_cols,binomial_test_result_cols,**kwargs)
    
    #get necessary data for model fitting
    phecode_info = data.phecode_info
    trajectory_ineligible = data.trajectory['ineligible_disease']
    trajectory_eligible_withdate = data.trajectory['eligible_disease_withdate']
    all_diagnosis_level = data.trajectory['all_diagnosis_level'] #extract the new history list
    phenotype_df = data.phenotype_df
    exp_col = data.get_attribute('phenotype_info')['phenotype_col_dict']['Exposure']
    id_col = data.get_attribute('phenotype_info')['phenotype_col_dict']['Participant ID']
    exposed_index = phenotype_df[phenotype_df[exp_col]==1].index
    phenotype_df_exposed = data.phenotype_df.loc[exposed_index,[id_col]+covariates]
    
    #get all disease pairs with significant comorbidity strength
    comorbidity_sig = comorbidity_strength_result[(comorbidity_strength_result[significance_phi_col]==True) & 
                                                  (comorbidity_strength_result[significance_RR_col]==True)]
    if len(comorbidity_sig) == 0:
        raise ValueError("No disease pair remained after filtering on significance of phi-correlation and RR.")
    all_diseases_lst = list(set(comorbidity_sig[phecode_d1_col].to_list() + comorbidity_sig[phecode_d2_col].to_list()))
    phecode_sig = data.get_attribute('significant_phecodes')
    invalid_disease = [x for x in all_diseases_lst if x not in phecode_sig]
    if invalid_disease:
        raise ValueError(f"The following phecode from the 'comorbidity_strength_result' are not in the list of PheWAS significant phecode: {invalid_disease}.")
    
    #create other diseases variables
    if parameter_dict['method'] in ['RPCN','PCN_PCA']:        
        for disease in all_diseases_lst:
            phenotype_df_exposed[str(disease)] = phenotype_df_exposed[id_col].apply(lambda x: 1 if disease in all_diagnosis_level[x] else 0)

    time_start = time.time()
    #list of disease pair
    result_all = []
    if n_process == 1:
        for d1,d2 in tqdm(comorbidity_sig[[phecode_d1_col,phecode_d2_col]].values, miniters=20,mininterval=60,smoothing=0):
            result_all.append(logistic_model_wrapper(d1,d2,phenotype_df_exposed,id_col,trajectory_ineligible,trajectory_eligible_withdate,
                                                     all_diagnosis_level,covariates,all_diseases_lst,
                                                     log_file_final,parameter_dict))
    elif n_process > 1:
        parameters_all = []
        for d1,d2 in comorbidity_sig[[phecode_d1_col,phecode_d2_col]].values:
            parameters_all.append([d1,d2])
        with multiprocessing.get_context(start_mehtod).Pool(n_process, initializer=init_worker, initargs=(phenotype_df_exposed,id_col,trajectory_ineligible,trajectory_eligible_withdate,
                                                                                                            all_diagnosis_level,covariates,all_diseases_lst,
                                                                                                            log_file_final,parameter_dict)) as p:
            result_all = list(
                tqdm(
                    p.imap(logistic_model, parameters_all), 
                    total=len(parameters_all),
                    miniters=20,
                    mininterval=60,
                    smoothing=0,
                )
            )

    time_end = time.time()
    time_spent = (time_end - time_start)/60
    print(f'Comorbidity network analysis finished (elapsed {time_spent:.1f} mins)')
    
    #result column names based on different method
    df_columns_common = ['phecode_d1','phecode_d2','name_disease_pair','N_exposed','n_total','n_exposed/n_cases','n_exposed/n_controls',
                         'comorbidity_network_method','describe','co_vars_list','co_vars_zvalues']
    df_columns_dict = {'CN': df_columns_common+['comorbidity_beta','comorbidity_se','comorbidity_p','comorbidity_aic'],
                       'RPCN': df_columns_common+['alpha',
                                                  'comorbidity_beta','comorbidity_se','comorbidity_p','comorbidity_aic'],
                       'PCN_PCA': df_columns_common+['pc_sum_variance_explained',
                                                     'comorbidity_beta','comorbidity_se','comorbidity_p','comorbidity_aic']}
    df_columns = df_columns_dict[parameter_dict['method']]
    
    #generate result dataframe
    max_columns = max([len(x) for x in result_all])
    columns_selected = df_columns[0:max_columns]
    comorbidity_df = pd.DataFrame(result_all, columns=columns_selected)
    #annotate disease name and system
    for d in ['d1','d2']:
        comorbidity_df[f'disease_{d}'] = comorbidity_df[f'phecode_{d}'].apply(lambda x: phecode_info[x]['phenotype'])
        comorbidity_df[f'system_{d}'] = comorbidity_df[f'phecode_{d}'].apply(lambda x: phecode_info[x]['category'])
        comorbidity_df[f'sex_{d}'] = comorbidity_df[f'phecode_{d}'].apply(lambda x: phecode_info[x]['sex'])
    
    #p-value correction
    comorbidity_df = comorbidity_multipletests(comorbidity_df, correction=correction,cutoff=cutoff)
    
    return comorbidity_df


def comorbidity_multipletests(
    df:pd.DataFrame, 
    correction:str='bonferroni', 
    cutoff:float=0.05
) -> pd.DataFrame:
    """
    Adjusts comorbidity network analysis p-values for multiple comparisons using specified correction methods.

    Parameters:
    ----------
    df : pd.DataFrame
        DataFrame containing the results from the 'comorbidity_network' function.

    correction : str, default='bonferroni'
        Method for binomial p-value correction from the statsmodels.stats.multitest.multipletests.
        Available methods are:
        none : no correction
        bonferroni : one-step correction
        sidak : one-step correction
        holm-sidak : step down method using Sidak adjustments
        holm : step-down method using Bonferroni adjustments
        simes-hochberg : step-up method (independent)
        hommel : closed method based on Simes tests (non-negative)
        fdr_bh : Benjamini/Hochberg (non-negative)
        fdr_by : Benjamini/Yekutieli (negative)
        fdr_tsbh : two stage fdr correction (non-negative)
        fdr_tsbky : two stage fdr correction (non-negative)
        See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
    
    cutoff : float, default=0.05
        The significance threshold for adjusted binomial p-values.

    Returns:
    ----------
    pd.DataFrame
        A DataFrame that contains the comorbidity network analysis results with applied p-value corrections.

    """
    #data type check
    if not isinstance(df,pd.DataFrame):
        raise TypeError("The input 'df' must be a pandas DataFrame.")
    
    #check p-value correction method and cutoff
    correction_method_check(correction,cutoff)

    #multiple adjustment
    # loop to adjust
    for p_col,correction_,cutoff_ in [('comorbidity_p',correction,cutoff)]:
        #if p-value were not presented
        if p_col not in df.columns or len(df[~df[p_col].isna()])==0:
            print(f'No valid {p_col} found in the provided DataFrame, no p-value correction made.')
            return df
        else:
            df[p_col] = pd.to_numeric(df[p_col],errors='coerce')
            
            if correction_ == 'none':
                df[f'{p_col}_significance'] = df[p_col].apply(lambda x: True if x<=cutoff_ else False)
                df[f'{p_col}_adjusted'] = np.nan
            else:
                df = states_p_adjust(df,p_col,correction_,cutoff_,p_col,p_col)
    return df


def disease_trajectory(
    data:DiseaseNetworkData, 
    comorbidity_strength_result:pd.DataFrame, 
    binomial_test_result:pd.DataFrame, 
    method:str='RPCN', 
    matching_var_dict:dict={'sex':'exact'}, 
    matching_n:int=2, 
    max_n_cases:int=np.inf,
    global_sampling: bool=False, 
    covariates:list=None, 
    n_process:int=1, 
    log_file:str=None, 
    correction:str='bonferroni', 
    cutoff:float=0.05, 
    **kwargs
) -> pd.DataFrame:
    """
    Perform temporal comorbidity network (disease trajectory) analysis on disease pairs with significant comorbidity strength and temporal order, to identify pairs with confirmed temporal comorbidity associations.
    For each disease pair D1 → D2, a nested case-control dataset is constructed using incidence density sampling, treating D2 as the outcome and D1 as the exposure. 
    A logistic regression model is then applied to estimate the correlations between the diseases.

    Depending on the selected 'method', the function applies different statistical models to estimate the correlations for each disease pair:
    - **RPCN (Regularized Partial Correlation Network):**
        Utilizes L1-regularized conditional logistic regression to estimate partial correlations for each disease pair.
        Includes both phenotypic variables and other diseases present in the network as covariates.
        The L1 regularization term selects important confounding disease variables.
        After variable selection, a standard conditional logistic regression model is refitted to accurately estimate partial correlations, adjusting for phenotypic variables and the selected confounding diseases.
    - **PCN_PCA (Partial Correlation Network with PCA):**
        Applies a standard conditional logistic regression model for each disease pair.
        Adjusts for phenotypic variables and the top principal components (PCs) of other diseases in the network to estimate partial correlations.
    - **CN (Correlation Network):**
        Uses a standard conditional logistic regression to estimate simple correlations for each disease pair. Adjusts only for phenotypic variables.

    Parameters
    ----------
    data : DiseaseNetworkData
        DESCRIPTION.

    comorbidity_strength_result : pd.DataFrame
        DataFrame containing comorbidity strength analysis results produced by the 'DiNetxify.comorbidity_strength' function.
    
    binomial_test_result : pd.DataFrame
        DataFrame containing binomial test analysis results produced by the 'DiNetxify.binomial_test' function.

    method : str, default='RPCN'
        Specifies the comorbidity network analysis method to use. Choices are:
        - 'RPCN: Regularized Partial Correlation Network.
        - 'PCN_PCA: Partial Correlation Network with PCA.
        - 'CN': Correlation Network.
        
        **Additional Options for RPCN:**
        - 'alpha' : non-negative scalar
            The weight multiplying the l1 penalty term for other diseases covariates. 
            Ignored if 'auto_penalty' is enabled.
        - 'auto_penalty' : bool, default=True
            If 'True', automatically determine the optimal 'alpha' based on model AIC value.
        - 'alpha_range' : tuple, default=(1,15)
            When 'auto_penalty' is True, search the optimal 'alpha' in this range.
        - 'scaling_factor' : positive scalar, default=1
            The scaling factor for the alpha when 'auto_penalty' is True.
        
        **Additional Options for PCN_PCA:**
        - 'n_PC' : int, default=5
            Fixed number of principal components to include in each model.
        - 'explained_variance' : float
            Determines the number of principal components based on the cumulative explained variance. 
            Overrides 'n_PC' if specified.
    
    matching_var_dict : dict, default={'sex':'exact'}
        Specifies the matching variables and the criteria used for incidence density sampling.
        For categorical and binary variables, the matching criteria should always be 'exact'.
        For continuous variables, provide a scalar greater than 0 as the matching criterion, indicating the maximum allowed difference when matching.
        To include the required variable sex as a matching variable, always use 'sex' instead of its original column name.
        For other covariates specified in the 'DiNetxify.DiseaseNetworkData.phenotype_data()' function, use their original column names.
    
    matching_n : int, default=2
        Specifies the maximum number of matched controls for each case.
    
    max_n_cases : int, default=np.inf
        Specifies the maximum number of D2 cases to include in the analysis.
        If the number of D2 cases exceeds this value, a random sample of cases will be selected.
    
    global_sampling : bool, default=False
        Indicates whether to perform independent incidence density sampling for each D1→D2 pair (if False),
        or to perform a single incidence density sampling for all Dx→D2 pairs with separate regression models for each D1→D2 pair (if True).
        Global sampling is recommended when processing large datasets, though it might reduce result heterogeneity.
    
    covariates : list, default=None
        List of phenotypic covariates to include in the model.
        By default, includes all covariates specified in the 'DiNetxify.DiseaseNetworkData.phenotype_data()' function.
        Categorical and binary variables used for matching should not be included as covariates.
        Continuous variables used for matching can be included as covariates, but caution is advised.
        To include the required variable sex as a covariate, always use 'sex' instead of its original column name.
        For other covariates specified in the 'DiNetxify.DiseaseNetworkData.phenotype_data()' function, use their original column names.

    n_process : int, default=1
        Specifies the number of parallel processes to use for the analysis.
        Multiprocessing is enabled when `n_process` is set to a value greater than one.

    correction : str, default='bonferroni'
        Method for comorbidity network analysis p-value correction from the statsmodels.stats.multitest.multipletests.
        Available methods are:
        none : no correction
        bonferroni : one-step correction
        sidak : one-step correction
        holm-sidak : step down method using Sidak adjustments
        holm : step-down method using Bonferroni adjustments
        simes-hochberg : step-up method (independent)
        hommel : closed method based on Simes tests (non-negative)
        fdr_bh : Benjamini/Hochberg (non-negative)
        fdr_by : Benjamini/Yekutieli (negative)
        fdr_tsbh : two stage fdr correction (non-negative)
        fdr_tsbky : two stage fdr correction (non-negative)
        See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
    
    cutoff : float, default=0.05
        The significance threshold for adjusted comorbidity network analysis p-values.
    
    log_file : str, default=None
        Path and prefix for the text file where log will be recorded.
        If None, the log will be written to the temporary files directory with file prefix of DiseaseNet_trajectory_.
    
    **kwargs
        Analysis option
            enforce_time_interval : bool, default=True
                If set to True, applies the specified minimum and maximum time intervals when determining the D2 outcome among individuals diagnosed with D1. 
                These time interval requirements have been defined when calling the DiNetxify.DiseaseNetworkData.disease_pair() function.
    
        Additional keyword argument to define the required columns in 'comorbidity_strength_result' and 'binomial_test_result':
            phecode_d1_col : str, default='phecode_d1'
                Name of the column in 'comorbidity_strength_result' and 'binomial_test_result' that specifies the phecode identifiers for disease 1 of the disease pair.
            phecode_d2_col : str, default='phecode_d2'
                Name of the column in 'comorbidity_strength_result' and 'binomial_test_result' that specifies the phecode identifiers for disease 2 of the disease pair.
            significance_phi_col : str, default='phi_p_significance'
                Name of the column in 'comorbidity_strength_result' that indicates the significance of phi-correlation for each disease pair.
            significance_RR_col : str, default='RR_p_significance'
                Name of the column in 'comorbidity_strength_result' that indicates the significance of RR for each disease pair.
            significance_binomial_col : str default='binomial_p_significance'
                Name of the column in 'binomial_test_result' that indicates the significance of binomial test for each disease pair.

        RPCN Method Parameters:
            alpha : non-negative scalar
                The weight multiplying the l1 penalty term for other diseases covariates. 
                Ignored if 'auto_penalty' is enabled.
            auto_penalty : bool, default=True
                If 'True', automatically determines the best 'alpha' based on model AIC value.
            alpha_range : tuple, default=(1,15)
                When 'auto_penalty' is True, search the optimal 'alpha' in this range.
            scaling_factor : positive scalar, default=1
                The scaling factor for the alpha when 'auto_penalty' is True.

        PCN_PCA Method Parameters:
            n_PC : int, default=5
                Fixed number of principal components to include in each model.
            explained_variance : float
                Cumulative explained variance threshold to determine the number of principal components. 
                Overrides 'n_PC' if specified.

    Returns:
    ----------
    pd.DataFrame
        A pandas DataFrame object that contains the results of trajectory analysis.

    """
    
    #data type check
    if not isinstance(data,DiseaseNetworkData):
        raise TypeError("The input 'data' must be a DiseaseNetworkData object.")
    
    #attribute check
    data_attrs = ['history','diagnosis', 'n_diagnosis','trajectory','phenotype_df']
    for attr in data_attrs:
        if getattr(data, attr) is None:
            raise ValueError(f"Attribute '{attr}' is empty.")
    
    #check comorbidity strength estimation result
    if not isinstance(comorbidity_strength_result,pd.DataFrame):
        raise TypeError("The provided input 'comorbidity_strength_result' must be a pandas DataFrame.")
    
    #check binomial test result
    if not isinstance(binomial_test_result,pd.DataFrame):
        raise TypeError("The provided input 'binomial_test_result' must be a pandas DataFrame.")
    
    #allowed methods
    allowed_methods = {'RPCN', 'PCN_PCA', 'CN'}
    if method not in allowed_methods:
        raise ValueError(f"Invalid method '{method}'. Allowed methods are: {allowed_methods}.")
    
    #check the matching variables and matching criteria
    matching_var_check(matching_var_dict,data.get_attribute('phenotype_info'))
    
    #check covariates
    covariates = covariates_check(covariates,data.get_attribute('phenotype_info'),matching_var_dict)

    #check matching number
    if not isinstance(matching_n,int) or matching_n<1:
        raise ValueError("The input 'matching_n' must be a positive integer.")
    
    #check max number of cases
    if max_n_cases == np.inf:
        pass
    elif not isinstance(max_n_cases,int) or max_n_cases<1:
        raise ValueError("The input 'max_n_cases' must be a positive integer.")
    
    #check global sampling
    if not isinstance(global_sampling,bool):
        raise TypeError("The input 'global_sampling' must be a boolean.")
    
    #check number of process
    n_process,start_mehtod = n_process_check(n_process,'trajectory')
    if n_process>1:
        import multiprocessing
        from .conditional_logistic import logistic_model,init_worker #use original function as main function and init_worker to initialize global variables
    else:
        from .conditional_logistic import logistic_model_wrapper #use wrapper function as main function

    #check p-value correction method and cutoff
    correction_method_check(correction,cutoff)
    
    #check log files
    log_file_final,message = log_file_detect(log_file,'disease_trajectory')
    print(message)
    
    #check **kwargs and get parameters
    parameter_dict,phecode_d1_col,phecode_d2_col,significance_phi_col,significance_RR_col,significance_binomial_col = check_kwargs_com_tra(method,comorbidity_strength_result.columns,binomial_test_result.columns,**kwargs)
    
    #get necessary data for model fitting
    phecode_info = data.phecode_info
    
    trajectory_ineligible = data.trajectory['ineligible_disease']
    all_diagnosis_level = data.trajectory['all_diagnosis_level'] #extract the new history list
    trajectory_eligible_withdate = data.trajectory['eligible_disease_withdate']
    phenotype_df = data.phenotype_df
    exp_col = data.get_attribute('phenotype_info')['phenotype_col_dict']['Exposure']
    id_col = data.get_attribute('phenotype_info')['phenotype_col_dict']['Participant ID']
    end_date_col = data.get_attribute('phenotype_info')['phenotype_col_dict']['End date']
    exposed_index = phenotype_df[phenotype_df[exp_col]==1].index
    phenotype_df_exposed = data.phenotype_df.loc[exposed_index,[id_col,end_date_col]+covariates+list(matching_var_dict.keys())]
    min_interval = data.min_interval_days
    max_interval = data.max_interval_days
    
    #get all disease pairs with significant temporal orders
    trajectory_sig = binomial_test_result[binomial_test_result[significance_binomial_col]==True]
    #get all disease pairs with significant comorbidity strength
    comorbidity_sig = comorbidity_strength_result[(comorbidity_strength_result[significance_phi_col]==True) & 
                                                  (comorbidity_strength_result[significance_RR_col]==True)]
    
    if len(trajectory_sig) == 0 or len(comorbidity_sig) == 0:
        raise ValueError("No disease pair remained after filtering on significance of phi-correlation/RR or binomial test.")
    all_diseases_lst = list(set(comorbidity_sig[phecode_d1_col].to_list() + comorbidity_sig[phecode_d2_col].to_list()))
    phecode_sig = data.get_attribute('significant_phecodes')
    invalid_disease = [x for x in all_diseases_lst if x not in phecode_sig]
    if invalid_disease:
        raise ValueError(f"The following phecode from the 'comorbidity_strength_result' are not in the list of PheWAS significant phecode: {invalid_disease}.")
    
    #create other diseases variables
    if parameter_dict['method'] in ['RPCN','PCN_PCA']:        
        for disease in all_diseases_lst:
            phenotype_df_exposed[str(disease)] = phenotype_df_exposed[id_col].apply(lambda x: 1 if disease in all_diagnosis_level[x] else 0)

    #create list of disease pairs for loop
    #if global sampling is True, the list of disease pairs will be based on the unique d2
    #if global sampling is False, ignore the unique d2 and use all disease pairs
    time_start = time.time()
    parameters_all = []
    if global_sampling==True:
        d2_lst = trajectory_sig[phecode_d2_col].unique()
        for d2 in d2_lst:
            d1_lst = trajectory_sig[trajectory_sig[phecode_d2_col]==d2][phecode_d1_col].to_list()
            parameters_all.append([d1_lst,d2])
    else:
        for d1,d2 in trajectory_sig[[phecode_d1_col,phecode_d2_col]].values:
            parameters_all.append([[d1],d2])

    result_all = []
    if n_process == 1:
        for d1_lst, d2 in tqdm(
            parameters_all,
            miniters=10,
            mininterval=60,
            smoothing=0
        ):
            result_all.append(logistic_model_wrapper(d1_lst,d2,phenotype_df_exposed,id_col,end_date_col,trajectory_ineligible,min_interval,max_interval,
                                                    trajectory_eligible_withdate,all_diagnosis_level,covariates,all_diseases_lst,
                                                    matching_var_dict,matching_n,max_n_cases,log_file_final,parameter_dict))
    elif n_process > 1:
        with multiprocessing.get_context(start_mehtod).Pool(n_process, initializer=init_worker, initargs=(phenotype_df_exposed,id_col,end_date_col,trajectory_ineligible,min_interval,max_interval,
                                                                                                          trajectory_eligible_withdate,all_diagnosis_level,covariates,all_diseases_lst,
                                                                                                          matching_var_dict,matching_n,max_n_cases,log_file_final,parameter_dict)) as p:
            result_all = list(
                tqdm(
                    p.imap(logistic_model, parameters_all), 
                    total=len(parameters_all),
                    miniters=10,
                    mininterval=60,
                    smoothing=0,
                )
            )

    time_end = time.time()
    time_spent = (time_end - time_start)/60
    print(f'Disease trajectory analysis finished (elapsed {time_spent:.1f} mins)')

    #unpacked the result
    result_all = [result for sublist in result_all for result in sublist]
    
    #result column names based on different method
    df_columns_common = ['phecode_d1','phecode_d2','name_disease_pair','N_exposed','n_total','n_exposed/n_cases','n_exposed/n_controls',
                         'trajectory_method','describe','co_vars_list','co_vars_zvalues']
    df_columns_dict = {'CN': df_columns_common+['trajectory_beta','trajectory_se','trajectory_p','trajectory_aic'],
                       'RPCN': df_columns_common+['alpha',
                                                  'trajectory_beta','trajectory_se','trajectory_p','trajectory_aic'],
                       'PCN_PCA': df_columns_common+['pc_sum_variance_explained',
                                                     'trajectory_beta','trajectory_se','trajectory_p','trajectory_aic']}
    df_columns = df_columns_dict[parameter_dict['method']]
    
    #generate result dataframe
    max_columns = max([len(x) for x in result_all])
    columns_selected = df_columns[0:max_columns]
    comorbidity_df = pd.DataFrame(result_all, columns=columns_selected)
    #annotate disease name and system
    for d in ['d1','d2']:
        comorbidity_df[f'disease_{d}'] = comorbidity_df[f'phecode_{d}'].apply(lambda x: phecode_info[x]['phenotype'])
        comorbidity_df[f'system_{d}'] = comorbidity_df[f'phecode_{d}'].apply(lambda x: phecode_info[x]['category'])
        comorbidity_df[f'sex_{d}'] = comorbidity_df[f'phecode_{d}'].apply(lambda x: phecode_info[x]['sex'])
    
    #p-value correction
    comorbidity_df = trajectory_multipletests(comorbidity_df, correction=correction, cutoff=cutoff)
    return comorbidity_df


def trajectory_multipletests(
    df:pd.DataFrame, 
    correction:str='bonferroni', 
    cutoff:float=0.05
) -> pd.DataFrame:
    """
    Adjusts trajectory analysis p-values for multiple comparisons using specified correction methods.

    Parameters:
    ----------
    df : pd.DataFrame
        DataFrame containing the results from the 'disease_trajectory' function.

    correction : str, default='bonferroni'
        Method for binomial p-value correction from the statsmodels.stats.multitest.multipletests.
        Available methods are:
        none : no correction
        bonferroni : one-step correction
        sidak : one-step correction
        holm-sidak : step down method using Sidak adjustments
        holm : step-down method using Bonferroni adjustments
        simes-hochberg : step-up method (independent)
        hommel : closed method based on Simes tests (non-negative)
        fdr_bh : Benjamini/Hochberg (non-negative)
        fdr_by : Benjamini/Yekutieli (negative)
        fdr_tsbh : two stage fdr correction (non-negative)
        fdr_tsbky : two stage fdr correction (non-negative)
        See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
    
    cutoff : float, default=0.05
        The significance threshold for adjusted binomial p-values.

    Returns:
    ----------
    pd.DataFrame
        A DataFrame that contains the trajectory test results with applied p-value corrections.

    """
    #data type check
    if not isinstance(df,pd.DataFrame):
        raise TypeError("The input 'df' must be a pandas DataFrame.")
    
    #check p-value correction method and cutoff
    correction_method_check(correction,cutoff)

    #multiple adjustment
    # loop to adjust
    for p_col,correction_,cutoff_ in [('trajectory_p',correction,cutoff)]:
        #if p-value were not presented
        if p_col not in df.columns or len(df[~df[p_col].isna()])==0:
            print(f'No valid {p_col} found in the provided DataFrame, no p-value correction made.')
            return df
        else:
            df[p_col] = pd.to_numeric(df[p_col],errors='coerce')
            
            if correction_ == 'none':
                df[f'{p_col}_significance'] = df[p_col].apply(lambda x: True if x<=cutoff_ else False)
                df[f'{p_col}_adjusted'] = np.nan
            else:
                df = states_p_adjust(df,p_col,correction_,cutoff_,p_col,p_col)
    return df















