# -*- coding: utf-8 -*-
"""
Created on Thu Dec 5 19:48:09 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

import pandas as pd
import numpy as np
import time
from .data_management import DiseaseNetworkData
from statsmodels.stats.multitest import multipletests
from .cox import cox_conditional,cox_unconditional
from .utility import log_file_detect

import warnings
warnings.filterwarnings('ignore')

def phewas(data:DiseaseNetworkData, sex_adjustment:bool=True, proportion_threshold:float=None, n_threshold:int=None, 
           n_cpus:int=1, correction:str='bonferroni', cutoff:float=0.05, system_inc:list=None, system_exl:list=None, 
           phecode_inc:list=None, phecode_exl:list=None, log_file:str=None) -> pd.DataFrame:
    """
    Conducts Phenome-wide association studies (PheWAS) using the specified DiseaseNetworkData object.

    Parameters:
    ----------
    data : DiseaseNetworkData
        DiseaseNetworkData object.
    
    sex_adjustment : bool, default=True
        Whether sex should be included as an additional covariate in the Cox model.
        For example, for matched-cohort study where sex is one of the matching variables, it should not be included as a covariate.
        For unmatched cohort study on the hand, sex is normally included as a covariate.
    
    proportion_threshold : float
        Minimum proportion of cases among the exposed group to include a phecode in the analysis. 
        proportion_threshold and n_threshold are mutually exclusive.
    
    n_threshold : int
        Minimum number of cases among the exposed group to include a phecode in the analysis. 
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

    Returns:
    ----------
    pd.DataFrame
        A pandas DataFrame object that contains the results of PheWAS analysis.

    """

    #data type check
    if not isinstance(data,DiseaseNetworkData):
        raise ValueError("Invalid 'data' type: expected a DiseaseNetworkData object.")
    #retrieve phecode information
    phecode_info = data.phecode_info
    
    #check threshold
    n_exposed = data.get_attribute('phenotype_statistics')['n_exposed']
    if proportion_threshold and n_threshold:
        raise ValueError("'n_threshold' and 'proportion_threshold' cannot be specified at the same time.")
    if proportion_threshold:
        if not isinstance(proportion_threshold,float):
            raise ValueError("The 'proportion_threshold' must be a floating-point number.")
        if proportion_threshold<0 or proportion_threshold>1:
            raise ValueError("'proportion_threshold' must be between 0 and 1.")
        else:
            n_threshold = int(n_exposed)*proportion_threshold
    elif n_threshold:
        if not isinstance(n_threshold,int):
            raise ValueError("The 'n_threshold' must be an integer.")
        if n_threshold<0 or n_threshold>n_exposed:
            raise ValueError("'n_threshold' must be a non-negative integer less than or equal to the number of exposed indivduals.")
    else:
        raise ValueError("Either 'n_threshold' nor 'proportion_threshold' is specified.")
    
    #check number of CPUs
    if not isinstance(n_cpus,int):
        raise ValueError("The 'n_cpus' must be an integer.")
    if n_cpus == 1:
        print('Multi-threading is not used.')
    elif n_cpus > 1:
        print(f'Use {n_cpus} CPU cores for PheWAS analysis.')
        import multiprocessing
    else:
        raise ValueError("The specified number of CPUs is not valid. Please enter a positive integer.")
    
    #check p-value correction method and cutoff
    if not isinstance(correction,str):
        raise ValueError("The 'adjustment' must be a string.")
    methods_lst = ['bonferroni','sidak', 'holm-sidak','holm','simes-hochberg',
                   'hommel','fdr_bh','fdr_by','fdr_tsbh','fdr_tsbky','none'] #list of available p-value correction method
    if correction not in methods_lst:
        raise ValueError(f"Choose from the following p-value correction methods: {methods_lst}")
    if not isinstance(cutoff,float):
        raise ValueError("The 'cutoff' must be a floating-point number.")
    if cutoff<=0 or cutoff>=1:
        raise ValueError("'cutoff' must be between 0 and 1, exclusive.")
    
    #check inclusion and exclusion list
    phecode_lst_all = list(phecode_info.keys())
    system_all = set([phecode_info[x]['category'] for x in phecode_lst_all])
    if system_inc and system_exl:
        raise ValueError("'system_inc' and 'system_exl' cannot both be specified.")
    if phecode_inc and phecode_exl:
        raise ValueError("'phecode_inc' and 'phecode_exl' cannot both be specified.")
    if (system_inc or system_exl) and (phecode_inc or phecode_exl):
        print('Warning: both phecode and system level filters applied may result in redundant or ambiguous outcomes.')
    if system_inc:
        if not isinstance(system_inc,list):
            raise ValueError("The 'system_inc' must be a list.")
        if len(system_inc) == 0:
            raise ValueError("The 'system_inc' list is empty.")
        system_unidentified = [x for x in system_inc if x not in system_all]
        if len(system_unidentified)>0:
            raise ValueError(f"The following phecode systems from 'system_inc' are invalid: {system_unidentified}")
        phecode_lst_all = [x for x in phecode_lst_all if phecode_info[x]['category'] in system_inc]
    if system_exl:
        if not isinstance(system_exl,list):
            raise ValueError("The 'system_exl' must be a list.")
        if len(system_exl) == 0:
            raise ValueError("The 'system_exl' list is empty.")
        system_unidentified = [x for x in system_exl if x not in system_all]
        if len(system_unidentified)>0:
            raise ValueError(f"The following phecode systems from 'system_exl' are invalid: {system_unidentified}")
        phecode_lst_all = [x for x in phecode_lst_all if phecode_info[x]['category'] not in system_exl]  
    if phecode_inc:
        if not isinstance(phecode_inc,list):
            raise ValueError("The 'phecode_inc' must be a list.")
        if len(phecode_inc) == 0:
            raise ValueError("The 'phecode_inc' list is empty.")
        phecode_unidentified = [x for x in phecode_inc if x not in phecode_lst_all]
        if len(phecode_unidentified)>0:
            raise ValueError(f"The following phecodes from 'phecode_inc' are invalid: {phecode_unidentified}")
        phecode_lst_all = [x for x in phecode_lst_all if x in phecode_inc]
    if phecode_exl:
        if not isinstance(phecode_exl,list):
            raise ValueError("The 'phecode_exl' must be a list.")
        if len(phecode_exl) == 0:
            raise ValueError("The 'phecode_exl' list is empty.")
        phecode_unidentified = [x for x in phecode_exl if x not in phecode_lst_all]
        if len(phecode_unidentified)>0:
            raise ValueError(f"The following phecodes from 'phecode_exl' are invalid: {phecode_unidentified}")
        phecode_lst_all = [x for x in phecode_lst_all if x not in phecode_exl]
    if len(phecode_lst_all) == 0:
        raise ValueError("No phecodes remain after applying filtering at the phecode and system levels.")
    
    #check log files
    log_file_final,message = log_file_detect(log_file)
    print(message)
    
    time_start = time.time()
    #list of phecode
    result_all = []
    if data.study_design == 'matched cohort':
        if n_cpus == 1:
            for phecode in phecode_lst_all:
                result_all.append(cox_conditional(data,n_threshold,phecode,sex_adjustment,log_file_final))
        elif n_cpus > 1:
            with multiprocessing.get_context('spawn').Pool(n_cpus) as p:
                parameters_all = []
                for phecode in phecode_lst_all:
                    parameters_all.append([data,n_threshold,phecode,sex_adjustment,log_file_final])
                result_all = p.starmap(cox_conditional, parameters_all)
    if data.study_design == 'cohort':
        if n_cpus == 1:
            for phecode in phecode_lst_all:
                result_all.append(cox_unconditional(data,n_threshold,phecode,sex_adjustment,log_file_final))
        elif n_cpus > 1:
            with multiprocessing.get_context('spawn').Pool(n_cpus) as p:
                parameters_all = []
                for phecode in phecode_lst_all:
                    parameters_all.append([data,n_threshold,phecode,sex_adjustment,log_file_final])
                result_all = p.starmap(cox_unconditional, parameters_all)    
    
    time_end = time.time()
    time_spent = (time_end - time_start)/60
    print(f'PheWAS analysis finished (elapsed {time_spent:.1f} mins)')
    
    #generate result dataframe
    max_columns = max([len(x) for x in result_all])
    columns = ['phecode','disease','system','N_cases_exposed','describe','exposed_group','unexposed_group','coef','se','p']
    columns_selected = columns[0:max_columns]
    phewas_df = pd.DataFrame(result_all, columns=columns_selected)
    
    #if p-value were not presented
    if 'p' not in phewas_df.columns or len(phewas_df[~phewas_df['p'].isna()])==0:
        print('Warning: none of the phecodes yielded a valid p-value from the Cox analysis.')
        return phewas_df
    
    #multiple adjustment
    phewas_df['p'] = pd.to_numeric(phewas_df['p'],errors='coerce')
    if correction == 'none':
        phewas_df['significance'] = phewas_df['p'].apply(lambda x: True if x<=cutoff else False)
    else:
        phewas_df_na = phewas_df[phewas_df['p'].isna()]
        phewas_df_nona = phewas_df[~phewas_df['p'].isna()]
        reject_, corrected_p, _, _ = multipletests(phewas_df_nona['p'],method=correction,alpha=cutoff)
        phewas_df_nona['significance'] = reject_
        phewas_df_nona['p_adjusted'] = corrected_p
        phewas_df_na['significance'] = False
        phewas_df_na['p_adjusted'] = np.NaN
        phewas_df = pd.concat([phewas_df_nona,phewas_df_na])
    
    return phewas_df

def phewas_multipletests(phewas_result:pd.DataFrame, correction:str='bonferroni', cutoff:float=0.05) -> pd.DataFrame:
    """
    Adjusts PheWAS p-values for multiple comparisons using specified correction methods.

    Parameters:
    ----------
    phewas_result : pd.DataFrame
        DataFrame containing the results from the phewas function.

    correction : str, default='bonferroni'
        Method for p-value correction from the statsmodels.stats.multitest.multipletests.
        Available methods are:
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
    if not isinstance(phewas_result,pd.DataFrame):
        raise ValueError("The provided input 'df' must be a pandas DataFrame.")
    
    #check p-value correction method and cutoff
    if not isinstance(correction,str):
        raise ValueError("The 'adjustment' parameter must be a string.")
    methods_lst = ['bonferroni','sidak', 'holm-sidak','holm','simes-hochberg',
                   'hommel','fdr_bh','fdr_by','fdr_tsbh','fdr_tsbky'] #list of available p-value correction method
    if correction not in methods_lst:
        raise ValueError(f"Choose from the following p-value correction methods: {methods_lst}")
    if not isinstance(cutoff,float):
        raise ValueError("The 'cutoff' must be a floating-point number.")
    if cutoff<=0 or cutoff>=1:
        raise ValueError("'cutoff' must be between 0 and 1, exclusive.")

    #if p-value were not presented
    if 'p' not in phewas_result.columns or len(phewas_result[~phewas_result['p'].isna()])==0:
        raise ValueError('No valid p-values found in the provided DataFrame.')
    
    #multiple adjustment
    phewas_result['p'] = pd.to_numeric(phewas_result['p'],errors='coerce')
    phewas_df_na = phewas_result[phewas_result['p'].isna()]
    phewas_df_nona = phewas_result[~phewas_result['p'].isna()]
    reject_, corrected_p, _, _ = multipletests(phewas_df_nona['p'],method=correction,alpha=cutoff)
    phewas_df_nona['significance'] = reject_
    phewas_df_nona['p_adjusted'] = corrected_p
    phewas_df_na['significance'] = False
    phewas_df_na['p_adjusted'] = np.NaN
    phewas_result = pd.concat([phewas_df_nona,phewas_df_na])
    
    return phewas_result


"""
def comorbidity_analysis(data:pd.DataFrame, n_cpus:int, adjustment:str) -> pd.DataFrame:
    threshold = 1
    d1d2_lst = []
    result = []
    d1_list = list(data['d1'])
    d2_list = list(data['d2'])
    for idx in range(len(d2_list)):
        d1d2_lst.append([d1_list[idx],d2_list[idx]])

    def comorbidity(d1d2_lst:list) -> list:
        for d1,d2 in d1d2_lst:
            data['flag'] = data['d_eligible'].apply(lambda x: d1 in x and d2 in x)
            df_ = data.loc[data['flag']==True]
            n = len(df_)
            c = sum([d1 in x and d2 in x for x in df_['inpatient_level1'].values]) #d1d2
            if c<=threshold:
                result.append([d1,d2,c,n,np.NaN,np.NaN,np.NaN,np.NaN,np.NaN])
                continue
            p1 = sum([d1 in x for x in df_['inpatient_level1'].values])
            p2 = sum([d2 in x for x in df_['inpatient_level1'].values])
            rr = (n*c)/(p1*p2) #
            theta = (1/c + 1/((p1*p2)/n) - 1/n - 1/n)**0.5 #RR
            t_ = abs(np.log(rr)/theta) #RRt
            p = (1-t.cdf(t_,n))*2 #RRP
            phi = (c*n-p1*p2)/(((p1*p2)*(n-p1)*(n-p2))**0.5) #phi
            if phi == 1:
                phi = 0.99
            z_phi = 0.5*np.log((1+phi)/(1-phi))
            z_phi_theta = (1/(n-3))**0.5
            z_phi_t = abs(z_phi/z_phi_theta)
            p_phi = (1-t.cdf(z_phi_t,n))*2
            result.append([d1,d2,c,n,rr,theta,p,phi,p_phi])
    pool = multiprocessing.Pool(n_cpus)
    pool.map(comorbidity())


def trajectory_analysis(data:pd.DataFrame, n_cpus:int, adjustment:str) -> pd.DataFrame:
    pool = multiprocessing.Pool(n_cpus)

    def exc_lst(disease:str, phecode_cate:pd.DataFrame, phecode_lst:list) -> set:
        lst = []
        exl_range = phecode_cate.loc[phecode_cate['phecode']==disease]['phecode_exclude_range'].values[0]
        if pd.isna(exl_range):
            return set(lst)
        else:
            for range_ in exl_range.split(','):
                exl_lower, exl_higher = float(range_.split('-')[0]), float(range_.split('-')[1])
                exl_list_index = np.where(np.all([phecode_lst>=exl_lower, phecode_lst<=exl_higher], axis=0))[0]
                exl_list_temp = phecode_lst[exl_list_index]
                lst += [x for x in exl_list_temp]
            return set(lst)
        
    def d1_d2(data:pd.DataFrame) -> pd.DataFrame:
        exposed = data.copy()
        id_ = 'eid'
        inpatient_variable = 'inpatient_level1'
        date_start_variable = 'dia_date'
        eligible_variable = 'd_eligible'

        array = exposed[[id_, date_start_variable, inpatient_variable, eligible_variable]].values
        d1d2_result = []
        total = len(array)
        time0 = time.time()
        for i in range(len(array)):
            if i%3000 == 0:
                print("Progress: %.1f%%" % (i/total*100))
                print("Time spent: %.1f mins" % ((time.time()-time0)/60))
            d1d2 = []
            dict_temp = array[i][2]
            eligible = array[i][-1]
            d_list = [x for x in dict_temp.keys() if x in eligible]
            length = len(d_list)
            if length <= 1:
                d1d2_result.append([])
                continue
            for j in range(length-1):
                for k in (range(j+1,length)):
                    d1 = d_list[j]
                    d2 = d_list[k]
                    if dict_temp.get(d2) > dict_temp.get(d1):
                        d1d2.append("%s-%s" % (d1,d2))
            d1d2_result.append(d1d2)
        exposed['d1d2'] = d1d2_result
        return exposed

    def main():
        pass
"""