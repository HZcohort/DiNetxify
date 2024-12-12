# -*- coding: utf-8 -*-
"""
Created on Thu Dec 5 19:48:09 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

import pandas as pd
import numpy as np
import time
from .data_management import DiseaseNetworkData
from .utility import log_file_detect,filter_phecodes,validate_threshold,validate_n_cpus,validate_correction_method,states_p_adjust

import warnings
warnings.filterwarnings('ignore')

def phewas(data:DiseaseNetworkData, sex_adjustment:bool=True, proportion_threshold:float=None, n_threshold:int=None, 
           n_cpus:int=1, correction:str='bonferroni', cutoff:float=0.05, system_inc:list=None, system_exl:list=None, 
           phecode_inc:list=None, phecode_exl:list=None, log_file:str=None, lifelines_disable:bool=False) -> pd.DataFrame:
    """
    Conducts Phenome-wide association studies (PheWAS) using the specified DiseaseNetworkData object.

    Parameters:
    ----------
    data : DiseaseNetworkData
        DiseaseNetworkData object.
    
    sex_adjustment : bool, default=True
        Whether sex should be included as an additional covariate in the Cox model.
        For example, for matched-cohort study where sex is one of the matching variables, it should not be included as a covariate.
        For unmatched cohort study on the other hand, sex is normally included as a covariate.
    
    proportion_threshold : float
        The minimum proportion of cases within the exposed group required for a phecode to be included in the PheWAS analysis.
        If the proportion of cases is below this threshold, the phecode is excluded from the analysis.
        proportion_threshold and n_threshold are mutually exclusive.
    
    n_threshold : int
        The minimum number of cases within the exposed group required for a phecode to be included in the PheWAS analysis.
        If the number of cases is below this threshold, the phecode is excluded from the analysis.
        n_threshold and proportion_threshold are mutually exclusive.      

    n_cpus : int, default=1
        Number of CPU cores to utilize for the analysis. 
        Multiprocessing is engaged if more than one core is specified.

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
        The significance threshold for adjusted p-values.
    
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
        If None, the log will be written to the temporary files directory with file prefix of DiseaseNet_.
    
    lifelines_disable : bool, default=False
        Whether to disable the use of lifelines. 
        While lifelines generally require a longer fitting time, they are more resilient to violations of model assumptions.

    Returns:
    ----------
    pd.DataFrame
        A pandas DataFrame object that contains the results of PheWAS analysis.

    """
    from .cox import cox_conditional,cox_unconditional
    
    #data type check
    if not isinstance(data,DiseaseNetworkData):
        raise TypeError("The input 'data' must be a DiseaseNetworkData object.")
    
    #attribute check
    data_attrs = ['phenotype_df', 'diagnosis', 'history']
    for attr in data_attrs:
        if getattr(data, attr) is None:
            raise ValueError(f"Attribute '{attr}' is empty.")

    #retrieve phecode information
    phecode_info = data.phecode_info
 
    #check sex adjustment
    if not isinstance(sex_adjustment,bool):
        raise TypeError("The input 'sex_adjustment' must be a bool.")
    
    #check lifelines_disable
    if not isinstance(lifelines_disable,bool):
        raise TypeError("The input 'lifelines_disable' must be a bool.")
    
    #check threshold
    n_exposed = data.get_attribute('phenotype_statistics')['n_exposed']
    n_threshold = validate_threshold(proportion_threshold,n_threshold,n_exposed)
    
    #check number of CPUs
    validate_n_cpus(n_cpus,'PheWAS')
    if n_cpus>1:
        import multiprocessing

    #check p-value correction method and cutoff
    validate_correction_method(correction,cutoff)
    
    #check inclusion and exclusion list
    phecode_lst_all = filter_phecodes(phecode_info,system_inc,system_exl,phecode_inc,phecode_exl)
    print(f'A total of {len(phecode_lst_all)} phecodes included in the PheWAS analysis.')
    
    #check log files
    log_file_final,message = log_file_detect(log_file,'phewas')
    print(message)

    time_start = time.time()
    #list of phecode
    result_all = []
    if data.study_design == 'matched cohort':
        if n_cpus == 1:
            for phecode in phecode_lst_all:
                result_all.append(cox_conditional(data,n_threshold,phecode,sex_adjustment,log_file_final,lifelines_disable))
        elif n_cpus > 1:
            with multiprocessing.get_context('spawn').Pool(n_cpus) as p:
                parameters_all = []
                for phecode in phecode_lst_all:
                    parameters_all.append([data,n_threshold,phecode,sex_adjustment,log_file_final,lifelines_disable])
                result_all = p.starmap(cox_conditional, parameters_all)
    if data.study_design == 'cohort':
        if n_cpus == 1:
            for phecode in phecode_lst_all:
                result_all.append(cox_unconditional(data,n_threshold,phecode,sex_adjustment,log_file_final,lifelines_disable))
        elif n_cpus > 1:
            with multiprocessing.get_context('spawn').Pool(n_cpus) as p:
                parameters_all = []
                for phecode in phecode_lst_all:
                    parameters_all.append([data,n_threshold,phecode,sex_adjustment,log_file_final,lifelines_disable])
                result_all = p.starmap(cox_unconditional, parameters_all)    
    
    time_end = time.time()
    time_spent = (time_end - time_start)/60
    print(f'PheWAS analysis finished (elapsed {time_spent:.1f} mins)')
    
    #generate result dataframe
    max_columns = max([len(x) for x in result_all])
    columns = ['phecode','disease','system','sex','N_cases_exposed','describe','exposed_group','unexposed_group','coef','se','p']
    columns_selected = columns[0:max_columns]
    phewas_df = pd.DataFrame(result_all, columns=columns_selected)
    
    #p-value correction
    phewas_df = phewas_multipletests(phewas_df, correction=correction,cutoff=cutoff)

    return phewas_df

def phewas_multipletests(df:pd.DataFrame, correction:str='bonferroni', cutoff:float=0.05) -> pd.DataFrame:
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
        The significance threshold for adjusted p-values.

    Returns:
    ----------
    pd.DataFrame
        A DataFrame that contains the PheWAS results with applied p-value corrections.

    """
    #data type check
    if not isinstance(df,pd.DataFrame):
        raise TypeError("The input 'df' must be a pandas DataFrame.")
    
    #check p-value correction method and cutoff
    validate_correction_method(correction,cutoff)

    #multiple adjustment
    # loop to adjust
    for p_col,correction_,cutoff_ in [('p',correction,cutoff)]:
        #if p-value were not presented
        if p_col not in df.columns or len(df[~df[p_col].isna()])==0:
            print(f'No valid {p_col} found in the provided DataFrame, no p-value correction made.')
            continue
        else:
            df[p_col] = pd.to_numeric(df[p_col],errors='coerce')
            
            if correction_ == 'none':
                df[f'{p_col}_significance'] = df[p_col].apply(lambda x: True if x<=cutoff_ else False)
                df[f'{p_col}_adjusted'] = np.NaN
            else:
                df = states_p_adjust(df,p_col,correction_,cutoff_,p_col,p_col)
            
    return df



def comorbidity_strength(data:DiseaseNetworkData, proportion_threshold:float=None, n_threshold:int=None, 
                         n_cpus:int=1, log_file:str=None, correction_phi:str='bonferroni', cutoff_phi:float=0.05, 
                         correction_RR:str='bonferroni', cutoff_RR:float=0.05) -> pd.DataFrame:
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

    n_cpus : int, default=1
        Number of CPU cores to utilize for the analysis. 
        Multiprocessing is engaged if more than one core is specified.

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
        The significance threshold for phi-correlatio adjusted p-values.

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
        The significance threshold for RR adjusted p-values.
    
    log_file : str, default=None
        Path and prefix for the text file where log will be recorded.
        If None, the log will be written to the temporary files directory with file prefix of DiseaseNet_.

    Returns:
    ----------
    pd.DataFrame
        A pandas DataFrame object that contains the results of PheWAS analysis.

    """
    from itertools import combinations
    from .comorbidity_strength import com_phi_rr
    
    #data type check
    if not isinstance(data,DiseaseNetworkData):
        raise TypeError("The input 'data' must be a DiseaseNetworkData object.")
    
    #attribute check
    data_attrs = ['trajectory']
    for attr in data_attrs:
        if getattr(data, attr) is None:
            raise ValueError(f"Attribute '{attr}' is empty.")
    
    #retrieve phecode information
    phecode_info = data.phecode_info
    trajectory_dict = data.trajectory
    
    #check threshold
    n_exposed = data.get_attribute('phenotype_statistics')['n_exposed']
    n_threshold = validate_threshold(proportion_threshold,n_threshold,n_exposed)
    
    #check number of CPUs
    validate_n_cpus(n_cpus,'PheWAS')
    if n_cpus>1:
        import multiprocessing

    #check p-value correction method and cutoff
    validate_correction_method(correction_phi,cutoff_phi)
    validate_correction_method(correction_RR,cutoff_RR)
    
    #check log files
    log_file_final,message = log_file_detect(log_file,'com_strength')
    print(message)
    
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
    
    time_start = time.time()
    #list of phecode
    result_all = []
    if n_cpus == 1:
        for d1,d2,describe in d1d2_pair_lst:
            result_all.append(com_phi_rr(trajectory_dict,d1,d2,describe,n_threshold,log_file_final))
    elif n_cpus > 1:
        with multiprocessing.get_context('spawn').Pool(n_cpus) as p:
            parameters_all = []
            for d1,d2,describe in d1d2_pair_lst:
                parameters_all.append([trajectory_dict,d1,d2,describe,n_threshold,log_file_final])
            result_all = p.starmap(com_phi_rr, parameters_all)

    time_end = time.time()
    time_spent = (time_end - time_start)/60
    print(f'Comorbidity strength estimation finished (elapsed {time_spent:.1f} mins)')

    #generate result dataframe
    max_columns = max([len(x) for x in result_all])
    columns = ['phecode_d1','phecode_d2','name_disease_pair','N_exposed','n_total'
               ,'n_d1d2_diagnosis','n_d1_diagnosis','n_d2_diagnosis','n_d1d2_pair',
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
        The significance threshold for phi-correlatio adjusted p-values.

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
        The significance threshold for RR adjusted p-values.

    Returns:
    ----------
    pd.DataFrame
        A DataFrame that contains the comorbidity strength estimation results with applied p-value corrections.

    """
    #data type check
    if not isinstance(df,pd.DataFrame):
        raise TypeError("The input 'df' must be a pandas DataFrame.")
    
    #check p-value correction method and cutoff
    validate_correction_method(correction_phi,cutoff_phi)
    validate_correction_method(correction_RR,cutoff_RR)

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
                df[f'{p_col}_adjusted'] = np.NaN
            else:
                df = states_p_adjust(df,p_col,correction_,cutoff_,p_col,p_col)
            
    return df



















