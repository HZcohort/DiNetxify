# -*- coding: utf-8 -*-
"""
Created on Thu Dec 5 19:48:09 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

import pandas as pd
import numpy as np
import time
import random
from .data_management import DiseaseNetworkData
from .utility import log_file_detect,filter_phecodes,threshold_check,n_process_check,correction_method_check,states_p_adjust
from .utility import check_kwargs_com_tra,covariates_check,matching_var_check

import warnings
warnings.filterwarnings('ignore')


class DiseaseAnalysis:
    """Base class for disease analysis implementations"""
    
    def __init__(self, data: DiseaseNetworkData):
        self.data = data
        self.phecode_info = data.phecode_info

    def _setup_multiprocessing(self):
        """Configure multiprocessing environment"""
        if self.n_process > 1:
            import multiprocessing
            from .cox import init_worker
            return multiprocessing.get_context(self.start_method).Pool(
                self.n_process,
                initializer=init_worker,
                initargs=(self.data, self.covariates, self.n_threshold,
                         self.log_file, self.lifelines_disable)
            )
        return None

    def threshold_check(self, proportion_threshold, n_threshold, n_exposed):
        """Validate and calculate thresholds"""
        return threshold_check(proportion_threshold, n_threshold, n_exposed)


class PheWASAnalysis(DiseaseAnalysis):
    """PheWAS analysis implementation using OOP approach"""
    
    def __init__(self, data: DiseaseNetworkData):
        super().__init__(data)
        self.covariates = None
        self.system_inc = None
        self.system_exl = None  
        self.phecode_inc = None
        self.phecode_exl = None
        self.lifelines_disable = False

    def analyze(self, covariates=None, proportion_threshold=None, n_threshold=None,
                n_process=1, correction='bonferroni', cutoff=0.05, system_inc=None,
                system_exl=None, phecode_inc=None, phecode_exl=None, log_file=None,
                lifelines_disable=False, start_method=None):
        """Main analysis entry point"""
        
        # Data type check
        if not isinstance(self.data, DiseaseNetworkData):
            raise TypeError("The input 'data' must be a DiseaseNetworkData object.")
        
        # Attribute check
        data_attrs = ['phenotype_df', 'diagnosis', 'history']
        for attr in data_attrs:
            if getattr(self.data, attr) is None:
                raise ValueError(f"Attribute '{attr}' is empty.")
        
        # Parameter validation
        self.covariates = covariates_check(covariates, self.data.get_attribute('phenotype_info'))
        
        if not isinstance(lifelines_disable, bool):
            raise TypeError("The input 'lifelines_disable' must be a bool.")
        
        correction_method_check(correction, cutoff)
        
        self.n_process, self.start_method = n_process_check(n_process, 'PheWAS')
        self.correction = correction
        self.cutoff = cutoff
        self.system_inc = system_inc
        self.system_exl = system_exl
        self.phecode_inc = phecode_inc
        self.phecode_exl = phecode_exl
        self.lifelines_disable = lifelines_disable
        
        self.log_file, message = log_file_detect(log_file, 'phewas')
        print(message)
        
        # Threshold calculation
        n_exposed = self.data.get_attribute('phenotype_statistics')['n_exposed']
        self.n_threshold = self.threshold_check(proportion_threshold, n_threshold, n_exposed)
        
        # Filter phecodes
        phecode_lst_all = self._filter_phecodes()
        
        # Run analysis
        raw_results = self._run_analysis(phecode_lst_all)
        
        # Process results
        phewas_df = self._process_results(raw_results)
        
        return phewas_df

    def _filter_phecodes(self):
        """Filter phecodes based on inclusion/exclusion criteria"""
        phecode_list = filter_phecodes(self.phecode_info, self.system_inc, self.system_exl,
                                       self.phecode_inc, self.phecode_exl)
        print(f'A total of {len(phecode_list)} phecodes included in the PheWAS analysis.')
        return phecode_list

    def _run_analysis(self, phecode_list):
        """Core analysis logic"""
        # Setup multiprocessing
        pool = self._setup_multiprocessing()
        random.shuffle(phecode_list)
        from .cox import cox_conditional, cox_unconditional, cox_conditional_wrapper, cox_unconditional_wrapper
        
        try:
            if pool:
                if self.data.study_design == 'matched cohort':
                    results = pool.map(cox_unconditional, phecode_list)
                else:
                    results = pool.map(cox_conditional, phecode_list)
            else:
                results = []
                for phecode in phecode_list:
                    if self.data.study_design == 'matched cohort':
                        results.append(cox_unconditional_wrapper(
                            phecode, self.data, self.covariates, self.n_threshold,
                            self.log_file, self.lifelines_disable))
                    else:
                        results.append(cox_conditional_wrapper(
                            phecode, self.data, self.covariates, self.n_threshold,
                            self.log_file, self.lifelines_disable))
            return results
        finally:
            if pool:
                pool.close()

    def _process_results(self, raw_results):
        """Process results into final dataframe"""
        # Create dataframe
        max_columns = max(len(x) for x in raw_results)
        columns = ['phecode','disease','system','sex','N_cases_exposed','describe',
                  'exposed_group','unexposed_group','phewas_coef','phewas_se','phewas_p']
        columns_selected = columns[:max_columns]
        phewas_df = pd.DataFrame(raw_results, columns=columns_selected)

        # Handle registry study special case
        if self.data.study_design == "registry":
            phewas_df["phewas_p_significance"] = phewas_df["N_cases_exposed"].apply(
                lambda x: x >= self.n_threshold)
            return phewas_df

        # Apply p-value correction
        return phewas_multipletests(
            phewas_df, 
            correction=self.correction, 
            cutoff=self.cutoff)

# Maintain original function for backward compatibility

# Maintain original function for backward compatibility
def phewas(data:DiseaseNetworkData, *args, **kwargs) -> pd.DataFrame:
    """
    Conducts Phenome-wide association studies (PheWAS) using the specified DiseaseNetworkData object.

    Parameters:
    ----------
    data : DiseaseNetworkData
        DiseaseNetworkData object.
    
    covariates : list, default=None
        List of phenotypic covariates to include in the model.
        By default, includes 'sex' and all covariates specified in the 'DiseaseNet.DiseaseNetworkData.phenotype_data()' function.
        If you want to include the required variable sex as covariate, always use 'sex' rather than its original column name. 
        For other covariates you specified in the 'DiseaseNet.DiseaseNetworkData.phenotype_data()' function, use their original column name.
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
        If None, the log will be written to the temporary files directory with file prefix of DiseaseNet_.
    
    lifelines_disable : bool, default=False
        Whether to disable the use of lifelines. 
        While lifelines generally require a longer fitting time, they are more robust to violations of model assumptions.

    Returns:
    ----------
    pd.DataFrame
        A pandas DataFrame object that contains the results of PheWAS analysis.

    """
    analysis = PheWASAnalysis(data)
    return analysis.analyze(*args, **kwargs)

def phewas_multipletests(df:pd.DataFrame, 
                         correction:str='bonferroni', 
                         cutoff:float=0.05) -> pd.DataFrame:
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


class ComorbidityStrengthAnalysis(DiseaseAnalysis):
    """Comorbidity strength analysis implementation using OOP approach"""
    
    def __init__(self, data: DiseaseNetworkData):
        super().__init__(data)
        self.correction_phi = 'bonferroni'
        self.cutoff_phi = 0.05
        self.correction_RR = 'bonferroni'
        self.cutoff_RR = 0.05
        self.log_file = None

    def analyze(self, proportion_threshold=None, n_threshold=None, n_process=1,
                log_file=None, correction_phi='bonferroni', cutoff_phi=0.05,
                correction_RR='bonferroni', cutoff_RR=0.05):
        """Main analysis entry point"""
        # Store parameters
        self.n_process = n_process
        self.correction_phi = correction_phi
        self.cutoff_phi = cutoff_phi
        self.correction_RR = correction_RR 
        self.cutoff_RR = cutoff_RR
        self.log_file, _ = log_file_detect(log_file, 'com_strength')

        # Threshold calculation
        n_exposed = self.data.get_attribute('phenotype_statistics')['n_exposed']
        self.n_threshold = self.threshold_check(proportion_threshold, n_threshold, n_exposed)

        # Get all significant phecodes and generate pairs
        phecode_sig = self.data.get_attribute('significant_phecodes')
        disease_pairs = self._generate_disease_pairs(phecode_sig)
        
        # Run analysis
        raw_results = self._run_analysis(disease_pairs)
        
        # Process and return results
        return self._process_results(raw_results)

    def _generate_disease_pairs(self, phecodes_sig):
        """Generate all possible disease pairs with metadata"""
        from itertools import combinations
        
        pairs = []
        for d1, d2 in combinations(phecodes_sig, 2):
            sex1 = self.phecode_info[d1]['sex']
            sex2 = self.phecode_info[d2]['sex']
            desc = 'Disease pair with different sex specificity' if {sex1, sex2} == {'Female','Male'} else None
            pairs.append((d1, d2, desc))
        random.shuffle(pairs)
        return pairs

    def _run_analysis(self, disease_pairs):
        """Core analysis logic"""
        pool = self._setup_multiprocessing()
        
        try:
            if pool:
                from .comorbidity_strength import com_phi_rr
                return pool.starmap(com_phi_rr, [(d1, d2, desc) for d1, d2, desc in disease_pairs])
            else:
                from .comorbidity_strength import com_phi_rr_wrapper
                return [com_phi_rr_wrapper(d1, d2, desc, self.n_threshold, self.log_file) 
                       for d1, d2, desc in disease_pairs]
        finally:
            if pool:
                pool.close()

    def _process_results(self, raw_results):
        """Process results into final dataframe"""
        # Create dataframe
        columns = ['phecode_d1','phecode_d2','name_disease_pair','N_exposed','n_total',
                  'n_d1d2_diagnosis','n_d1_diagnosis','n_d2_diagnosis',
                  'n_d1d2_nontemporal','n_d1d2_temporal','n_d2d1_temporal','n_d1d2_pair',
                  'description','phi','phi_theta','phi_p','RR','RR_theta','RR_p']
        com_df = pd.DataFrame(raw_results, columns=columns[:len(raw_results[0])])
        
        # Annotate disease info
        for d in ['d1','d2']:
            com_df[f'disease_{d}'] = com_df[f'phecode_{d}'].map(lambda x: self.phecode_info[x]['phenotype'])
            com_df[f'system_{d}'] = com_df[f'phecode_{d}'].map(lambda x: self.phecode_info[x]['category'])
            com_df[f'sex_{d}'] = com_df[f'phecode_{d}'].map(lambda x: self.phecode_info[x]['sex'])
        
        # Apply p-value correction
        return self._apply_multipletests(com_df)

    def _apply_multipletests(self, df):
        """Handle p-value corrections for both phi and RR"""
        for p_col, correction, cutoff in [('phi_p', self.correction_phi, self.cutoff_phi),
                                        ('RR_p', self.correction_RR, self.cutoff_RR)]:
            if p_col in df.columns:
                df = states_p_adjust(df, p_col, correction, cutoff, p_col, p_col)
        return df

# Maintain original function for backward compatibility
def comorbidity_strength(data:DiseaseNetworkData, *args, **kwargs) -> pd.DataFrame:
    """Legacy wrapper for Comorbidity Strength analysis"""
    analysis = ComorbidityStrengthAnalysis(data)
    return analysis.analyze(*args, **kwargs)
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
        If None, the log will be written to the temporary files directory with file prefix of DiseaseNet_.

    Returns:
    ----------
    pd.DataFrame
        A pandas DataFrame object that contains the results of PheWAS analysis.

    """

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

class BinomialTestAnalysis(DiseaseAnalysis):
    """Binomial test analysis implementation using OOP approach"""
    
    def __init__(self, data: DiseaseNetworkData):
        super().__init__(data)
        self.correction = 'bonferroni'
        self.cutoff = 0.05
        self.log_file = None
        self.enforce_temporal_order = False

    def analyze(self, comorbidity_strength_result: pd.DataFrame, n_process=1,
                log_file=None, correction='bonferroni', cutoff=0.05,
                enforce_temporal_order=False, **kwargs):
        """Main analysis entry point"""
        # Store parameters
        self.n_process = n_process
        self.correction = correction
        self.cutoff = cutoff
        self.enforce_temporal_order = enforce_temporal_order
        self.log_file, _ = log_file_detect(log_file, 'binomial_test')
        
        # Validate inputs
        self._validate_inputs(comorbidity_strength_result, kwargs)
        
        # Get significant pairs
        sig_pairs = self._get_significant_pairs(comorbidity_strength_result, kwargs)
        
        # Run analysis
        raw_results = self._run_analysis(sig_pairs, comorbidity_strength_result, kwargs)
        
        # Process and return results
        return self._process_results(raw_results)

    def _validate_inputs(self, comorbidity_result, kwargs):
        """Validate input dataframe and parameters"""
        if not isinstance(comorbidity_result, pd.DataFrame):
            raise TypeError("comorbidity_strength_result must be a pandas DataFrame")
            
        # Check column existence
        cols = [
            kwargs.get('phecode_d1_col', 'phecode_d1'),
            kwargs.get('phecode_d2_col', 'phecode_d2'),
            kwargs.get('n_nontemporal_col', 'n_d1d2_nontemporal'),
            kwargs.get('n_temporal_d1d2_col', 'n_d1d2_temporal'),
            kwargs.get('n_temporal_d2d1_col', 'n_d2d1_temporal'),
            kwargs.get('significance_phi_col', 'phi_p_significance'),
            kwargs.get('significance_RR_col', 'RR_p_significance')
        ]
        for col in cols:
            if col not in comorbidity_result.columns:
                raise ValueError(f"Column {col} not found in comorbidity_strength_result")

        # Validate correction method
        correction_method_check(self.correction, self.cutoff)

    def _get_significant_pairs(self, comorbidity_result, kwargs):
        """Filter significant disease pairs"""
        phi_col = kwargs.get('significance_phi_col', 'phi_p_significance')
        rr_col = kwargs.get('significance_RR_col', 'RR_p_significance')
        
        return comorbidity_result[
            (comorbidity_result[phi_col] == True) & 
            (comorbidity_result[rr_col] == True)
        ]

    def _run_analysis(self, sig_pairs, comorbidity_result, kwargs):
        """Core analysis logic"""
        pool = self._setup_multiprocessing()
        n_temp_col = kwargs.get('n_temporal_d1d2_col', 'n_d1d2_temporal')
        n_temp2_col = kwargs.get('n_temporal_d2d1_col', 'n_d2d1_temporal')
        n_non_col = kwargs.get('n_nontemporal_col', 'n_d1d2_nontemporal')

        try:
            if pool:
                from .binomial import binomial
                params = [(row[0], row[1], row[2], row[3], row[4], self.enforce_temporal_order, self.log_file)
                         for row in sig_pairs[[kwargs.get('phecode_d1_col', 'phecode_d1'),
                                             kwargs.get('phecode_d2_col', 'phecode_d2'),
                                             n_non_col, n_temp_col, n_temp2_col]].values]
                return pool.starmap(binomial, params)
            else:
                from .binomial import binomial_wrapper
                return [binomial_wrapper(d1, d2, n_com, n_d1d2, n_d2d1, 
                                       self.enforce_temporal_order, self.log_file)
                      for d1, d2, n_com, n_d1d2, n_d2d1 in 
                      sig_pairs[[kwargs.get('phecode_d1_col', 'phecode_d1'),
                               kwargs.get('phecode_d2_col', 'phecode_d2'),
                               n_non_col, n_temp_col, n_temp2_col]].values]
        finally:
            if pool:
                pool.close()

    def _process_results(self, raw_results):
        """Process results into final dataframe"""
        columns = ['phecode_d1','phecode_d2','name_disease_pair',
                 'n_d1d2_nontemporal','n_d1d2_temporal','n_d2d1_temporal',
                 'binomial_p','binomial_proportion','binomial_proportion_ci']
        bino_df = pd.DataFrame(raw_results, columns=columns[:len(raw_results[0])])
        
        # Annotate disease info
        for d in ['d1','d2']:
            bino_df[f'disease_{d}'] = bino_df[f'phecode_{d}'].map(lambda x: self.phecode_info[x]['phenotype'])
            bino_df[f'system_{d}'] = bino_df[f'phecode_{d}'].map(lambda x: self.phecode_info[x]['category'])
            bino_df[f'sex_{d}'] = bino_df[f'phecode_{d}'].map(lambda x: self.phecode_info[x]['sex'])
        
        # Apply p-value correction
        return self._apply_multipletests(bino_df)

    def _apply_multipletests(self, df):
        """Handle p-value corrections"""
        if 'binomial_p' in df.columns:
            df = states_p_adjust(df, 'binomial_p', self.correction, self.cutoff, 'binomial_p', 'binomial_p')
        return df

# Maintain original function for backward compatibility
def binomial_test(data:DiseaseNetworkData, *args, **kwargs) -> pd.DataFrame:
    """Legacy wrapper for Binomial Test analysis"""
    analysis = BinomialTestAnalysis(data)
    return analysis.analyze(*args, **kwargs)
    """
    Conduct binomial test for disease pairs with significant comorbidity stregnth to select those with significant temporal orders (i.e., D1 -> D2).

    Parameters:
    ----------
    data : DiseaseNetworkData
        DiseaseNetworkData object.
    
    comorbidity_strength_result : pd.DataFrame
        DataFrame containing comorbidity strength analysis results produced by the 'DiseaseNet.comorbidity_strength' function.

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
        If None, the log will be written to the temporary files directory with file prefix of DiseaseNet_.
    
    enforce_temporal_order : bool, default=False
        If True, exclude individuals with non-temporal D1-D2 pair when performing the test.
        If False, include all individuals, including those with non-temporal D1-D2 pair.
    
    **kwargs
        Additional keyword argument to define the required columns in 'comorbidity_strength_result':
            phecode_d1_col : str, default='phecode_d1'
                Name of the column in 'comorbidity_strength_result' that specifies the phecode identifiers for disease 1 of the disease pair.
            phecode_d2_col : str, default='phecode_d2'
                Name of the column in 'comorbidity_strength_result' that specifies the phecode identifiers for disease 2 of the disease pair.
            n_nontemporal_col : str, default='n_d1d2_nontemporal'
                Name of the column in 'comorbidity_strength_result' that specifies the number of individuals with non-temporal d1-d2 disease pair
            n_temporal_d1d2_col : str, default='n_d1d2_temporal'
                Name of the column in 'comorbidity_strength_result' that specifies the number of individuals with temporal d1->d2 disease pair.
            n_temporal_d2d1_col : str, default='n_d2d1_temporal'
                Name of the column in 'comorbidity_strength_result' that specifies the number of individuals with temporal d2->d1 disease pair.
            significance_phi_col : str, default='phi_p_significance'
                Name of the column in 'comorbidity_strength_result' that indicates the significance of phi-correlation for each disease pair.
            significance_RR_col : str, default='phi_p_significance'
                Name of the column in 'comorbidity_strength_result' that indicates the significance of RR for each disease pair.

    Returns:
    ----------
    pd.DataFrame
        A pandas DataFrame object that contains the results of binomial test.

    """

def binomial_multipletests(df:pd.DataFrame, correction:str='bonferroni', cutoff:float=0.05) -> pd.DataFrame:
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


class ComorbidityNetworkAnalysis(DiseaseAnalysis):
    """Comorbidity network analysis implementation using OOP approach"""
    
    def __init__(self, data: DiseaseNetworkData):
        super().__init__(data)
        self.method = 'RPCN'
        self.covariates = None
        self.correction = 'bonferroni'
        self.cutoff = 0.05
        self.log_file = None
        self.parameter_dict = {}
        # Initialize method-specific parameters like original code
        if self.method == 'RPCN':
            self.parameter_dict.update({
                'alpha': 1.0,
                'auto_penalty': True,
                'alpha_range': (1, 15),
                'scaling_factor': 1.0
            })
        elif self.method == 'PCN_PCA':
            self.parameter_dict.update({
                'n_PC': 5,
                'explained_variance': None
            })

    def analyze(self, comorbidity_strength_result: pd.DataFrame, 
                binomial_test_result: pd.DataFrame, method='RPCN',
                covariates=None, n_process=1, log_file=None,
                correction='bonferroni', cutoff=0.05, **kwargs):
        """Main analysis entry point"""
        # Store parameters
        self.method = method
        self.n_process = n_process
        self.covariates = covariates_check(covariates, self.data.get_attribute('phenotype_info'))
        self.correction = correction
        self.cutoff = cutoff
        self.log_file, _ = log_file_detect(log_file, 'comorbidity_network_analysis')
        self.parameter_dict = check_kwargs_com_tra(method, 
                                                 comorbidity_strength_result.columns,
                                                 binomial_test_result.columns,
                                                 **kwargs)[0]

        # Validate inputs
        self._validate_inputs(comorbidity_strength_result, binomial_test_result)

        # Get significant pairs
        sig_pairs = self._get_significant_pairs(comorbidity_strength_result, binomial_test_result)

        # Run analysis
        raw_results = self._run_analysis(sig_pairs)

        # Process and return results
        return self._process_results(raw_results)

    def _validate_inputs(self, comorbidity_result, binomial_result):
        """Validate input dataframes"""
        if not isinstance(comorbidity_result, pd.DataFrame) or not isinstance(binomial_result, pd.DataFrame):
            raise TypeError("Input results must be pandas DataFrames")
            
        # Check required columns exist
        cols = [
            self.parameter_dict.get('phecode_d1_col', 'phecode_d1'),
            self.parameter_dict.get('phecode_d2_col', 'phecode_d2'),
            self.parameter_dict.get('significance_phi_col', 'phi_p_significance'),
            self.parameter_dict.get('significance_RR_col', 'RR_p_significance'),
            self.parameter_dict.get('significance_binomial_col', 'binomial_p_significance')
        ]
        for col, df in [(cols[0], comorbidity_result), 
                       (cols[1], comorbidity_result),
                       (cols[4], binomial_result)]:
            if col not in df.columns:
                raise ValueError(f"Column {col} not found in input DataFrame")

    def _get_significant_pairs(self, comorbidity_result, binomial_result):
        """Filter significant disease pairs"""
        phi_col = self.parameter_dict.get('significance_phi_col', 'phi_p_significance')
        rr_col = self.parameter_dict.get('significance_RR_col', 'RR_p_significance')
        bin_col = self.parameter_dict.get('significance_binomial_col', 'binomial_p_significance')

        # Get intersection of significant pairs
        com_sig = comorbidity_result[
            (comorbidity_result[phi_col] == True) & 
            (comorbidity_result[rr_col] == True)
        ]
        bin_sig = binomial_result[binomial_result[bin_col] == True]
        
        return pd.merge(
            com_sig,
            bin_sig,
            on=[self.parameter_dict.get('phecode_d1_col', 'phecode_d1'),
                self.parameter_dict.get('phecode_d2_col', 'phecode_d2')]
        )

    def _run_analysis(self, sig_pairs):
        """Core analysis logic"""
        pool = self._setup_multiprocessing()
        phecode_d1_col = self.parameter_dict.get('phecode_d1_col', 'phecode_d1')
        phecode_d2_col = self.parameter_dict.get('phecode_d2_col', 'phecode_d2')

        try:
            if pool:
                from .unconditional_logistic import logistic_model
                params = [(row[0], row[1]) for row in sig_pairs[[phecode_d1_col, phecode_d2_col]].values]
                return pool.starmap(logistic_model, params)
            else:
                from .unconditional_logistic import logistic_model_wrapper
                return [logistic_model_wrapper(d1, d2, self.data, self.covariates,
                                             self.log_file, self.parameter_dict)
                      for d1, d2 in sig_pairs[[phecode_d1_col, phecode_d2_col]].values]
        finally:
            if pool:
                pool.close()

    def _process_results(self, raw_results):
        """Process results into final dataframe"""
        columns = ['phecode_d1','phecode_d2','name_disease_pair','N_exposed','n_total',
                 'n_with_d2/n_with_d1','n_with_d2/n_without_d1','comorbidity_network_method',
                 'describe','comorbidity_beta','comorbidity_se','comorbidity_p','comorbidity_aic']
        net_df = pd.DataFrame(raw_results, columns=columns[:len(raw_results[0])])
        
        # Annotate disease info
        for d in ['d1','d2']:
            net_df[f'disease_{d}'] = net_df[f'phecode_{d}'].map(lambda x: self.phecode_info[x]['phenotype'])
            net_df[f'system_{d}'] = net_df[f'phecode_{d}'].map(lambda x: self.phecode_info[x]['category'])
            net_df[f'sex_{d}'] = net_df[f'phecode_{d}'].map(lambda x: self.phecode_info[x]['sex'])
        
        # Apply p-value correction
        return self._apply_multipletests(net_df)

    def _apply_multipletests(self, df):
        """Handle p-value corrections"""
        if 'comorbidity_p' in df.columns:
            df = states_p_adjust(df, 'comorbidity_p', self.correction, self.cutoff, 'comorbidity_p', 'comorbidity_p')
        return df

# Maintain original function for backward compatibility
def comorbidity_network(data:DiseaseNetworkData, *args, **kwargs) -> pd.DataFrame:
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
        DESCRIPTION.

    comorbidity_strength_result : pd.DataFrame
        DataFrame containing comorbidity strength analysis results produced by the 'DiseaseNet.comorbidity_strength' function.
    
    binomial_test_result : pd.DataFrame
        DataFrame containing binomial test analysis results produced by the 'DiseaseNet.binomial_test' function.

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
        By default, includes ['sex'] and all covariates specified in the 'DiseaseNet.DiseaseNetworkData.phenotype_data()' function.
        To include the required variable sex as a covariate, always use 'sex' instead of its original column name.
        For other covariates specified in the 'DiseaseNet.DiseaseNetworkData.phenotype_data()' function, use their original column names.

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
        If None, the log will be written to the temporary files directory with file prefix of DiseaseNet_.
    
    **kwargs
        Additional keyword argument to define the required columns in 'comorbidity_strength_result' and 'binomial_test_result':
            phecode_d1_col : str, default='phecode_d1'
                Name of the column in 'comorbidity_strength_result' and 'binomial_test_result' that specifies the phecode identifiers for disease 1 of the disease pair.
            phecode_d2_col : str, default='phecode_d2'
                Name of the column in 'comorbidity_strength_result' and 'binomial_test_result' that specifies the phecode identifiers for disease 2 of the disease pair.
            significance_phi_col : str, default='phi_p_significance'
                Name of the column in 'comorbidity_strength_result' that indicates the significance of phi-correlation for each disease pair.
            significance_RR_col : str, default='phi_RR_significance'
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
    analysis = ComorbidityNetworkAnalysis(data)
    return analysis.analyze(*args, **kwargs)


def comorbidity_multipletests(df:pd.DataFrame, correction:str='bonferroni', cutoff:float=0.05) -> pd.DataFrame:
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

class DiseaseTrajectoryAnalysis(DiseaseAnalysis):
    """Disease trajectory analysis implementation using OOP approach"""
    
    def __init__(self, data: DiseaseNetworkData):
        super().__init__(data)
        self.method = 'RPCN'
        self.covariates = None
        self.matching_var_dict = {'sex': 'exact'}  # Default matching
        self.matching_n = 2
        self.correction = 'bonferroni'
        self.cutoff = 0.05
        self.log_file = None
        self.parameter_dict = {}
        
        # Validate matching variables on initialization like original code
        if 'sex' not in self.matching_var_dict:
            raise ValueError("'sex' must be included in matching_var_dict")
        if self.matching_var_dict['sex'] != 'exact':
            raise ValueError("Sex matching must use 'exact' method")

    def analyze(self, comorbidity_strength_result: pd.DataFrame,
                binomial_test_result: pd.DataFrame, method='RPCN',
                matching_var_dict=None, matching_n=2, covariates=None,
                n_process=1, log_file=None, correction='bonferroni',
                cutoff=0.05, enforce_time_interval=True, **kwargs):
        """Main analysis entry point"""
        # Store parameters
        self.method = method
        self.n_process = n_process
        self.matching_var_dict = matching_var_dict or {'sex': 'exact'}
        self.matching_n = matching_n
        self.covariates = covariates_check(covariates, self.data.get_attribute('phenotype_info'), self.matching_var_dict)
        self.correction = correction
        self.cutoff = cutoff
        self.enforce_time_interval = enforce_time_interval
        self.log_file, _ = log_file_detect(log_file, 'disease_trajectory') 
        self.parameter_dict = check_kwargs_com_tra(method,
                                                 comorbidity_strength_result.columns,
                                                 binomial_test_result.columns,
                                                 **kwargs)[0]

        # Validate inputs
        self._validate_inputs(comorbidity_strength_result, binomial_test_result)

        # Get significant pairs
        sig_pairs = self._get_significant_pairs(comorbidity_strength_result, binomial_test_result)

        # Run analysis
        raw_results = self._run_analysis(sig_pairs)

        # Process and return results
        return self._process_results(raw_results)

    def _validate_inputs(self, comorbidity_result, binomial_result):
        """Validate input dataframes and parameters"""
        if not isinstance(comorbidity_result, pd.DataFrame) or not isinstance(binomial_result, pd.DataFrame):
            raise TypeError("Input results must be pandas DataFrames")
            
        # Check required columns exist
        cols = [
            self.parameter_dict.get('phecode_d1_col', 'phecode_d1'),
            self.parameter_dict.get('phecode_d2_col', 'phecode_d2'),
            self.parameter_dict.get('significance_phi_col', 'phi_p_significance'),
            self.parameter_dict.get('significance_RR_col', 'RR_p_significance'),
            self.parameter_dict.get('significance_binomial_col', 'binomial_p_significance')
        ]
        for col, df in [(cols[0], comorbidity_result),
                       (cols[1], comorbidity_result),
                       (cols[4], binomial_result)]:
            if col not in df.columns:
                raise ValueError(f"Column {col} not found in input DataFrame")

        # Validate matching variables
        matching_var_check(self.matching_var_dict, self.data.get_attribute('phenotype_info'))

    def _get_significant_pairs(self, comorbidity_result, binomial_result):
        """Filter significant disease pairs with temporal order"""
        phi_col = self.parameter_dict.get('significance_phi_col', 'phi_p_significance')
        rr_col = self.parameter_dict.get('significance_RR_col', 'RR_p_significance')
        bin_col = self.parameter_dict.get('significance_binomial_col', 'binomial_p_significance')

        # Get intersection of significant pairs
        com_sig = comorbidity_result[
            (comorbidity_result[phi_col] == True) & 
            (comorbidity_result[rr_col] == True)
        ]
        bin_sig = binomial_result[binomial_result[bin_col] == True]
        
        return pd.merge(
            com_sig,
            bin_sig,
            on=[self.parameter_dict.get('phecode_d1_col', 'phecode_d1'),
                self.parameter_dict.get('phecode_d2_col', 'phecode_d2')]
        )

    def _run_analysis(self, sig_pairs):
        """Core analysis logic"""
        pool = self._setup_multiprocessing()
        phecode_d1_col = self.parameter_dict.get('phecode_d1_col', 'phecode_d1')
        phecode_d2_col = self.parameter_dict.get('phecode_d2_col', 'phecode_d2')

        try:
            if pool:
                from .conditional_logistic import logistic_model
                params = [(row[0], row[1]) for row in sig_pairs[[phecode_d1_col, phecode_d2_col]].values]
                return pool.starmap(logistic_model, params)
            else:
                from .conditional_logistic import logistic_model_wrapper
                return [logistic_model_wrapper(d1, d2, self.data, self.covariates,
                                             self.matching_var_dict, self.matching_n,
                                             self.log_file, self.parameter_dict)
                      for d1, d2 in sig_pairs[[phecode_d1_col, phecode_d2_col]].values]
        finally:
            if pool:
                pool.close()

    def _process_results(self, raw_results):
        """Process results into final dataframe"""
        columns = ['phecode_d1','phecode_d2','name_disease_pair','N_exposed','n_total',
                 'n_with_d2/n_with_d1','n_with_d2/n_without_d1','trajectory_method',
                 'describe','trajectory_beta','trajectory_se','trajectory_p','trajectory_aic']
        traj_df = pd.DataFrame(raw_results, columns=columns[:len(raw_results[0])])
        
        # Annotate disease info
        for d in ['d1','d2']:
            traj_df[f'disease_{d}'] = traj_df[f'phecode_{d}'].map(lambda x: self.phecode_info[x]['phenotype'])
            traj_df[f'system_{d}'] = traj_df[f'phecode_{d}'].map(lambda x: self.phecode_info[x]['category'])
            traj_df[f'sex_{d}'] = traj_df[f'phecode_{d}'].map(lambda x: self.phecode_info[x]['sex'])
        
        # Apply p-value correction
        return self._apply_multipletests(traj_df)

    def _apply_multipletests(self, df):
        """Handle p-value corrections"""
        if 'trajectory_p' in df.columns:
            df = states_p_adjust(df, 'trajectory_p', self.correction, self.cutoff, 'trajectory_p', 'trajectory_p')
        return df

# Maintain original function for backward compatibility
def disease_trajectory(data:DiseaseNetworkData, *args, **kwargs) -> pd.DataFrame:
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
        The DiseaseNetworkData object containing processed phenotype and diagnosis data.

    comorbidity_strength_result : pd.DataFrame
        DataFrame containing comorbidity strength analysis results from 'comorbidity_strength'.

    binomial_test_result : pd.DataFrame  
        DataFrame containing binomial test results from 'binomial_test'.

    method : str, default='RPCN'
        Analysis method: 'RPCN', 'PCN_PCA' or 'CN'

    matching_var_dict : dict, default={'sex':'exact'}
        Matching variables and criteria for incidence density sampling

    matching_n : int, default=2
        Maximum number of matched controls per case

    covariates : list, optional
        Phenotypic covariates to include in models

    n_process : int, default=1
        Number of parallel processes

    correction : str, default='bonferroni'
        Multiple testing correction method

    cutoff : float, default=0.05
        Significance threshold

    log_file : str, optional
        Path for log file output

    **kwargs
        Additional method-specific parameters:
        - RPCN: alpha, auto_penalty, alpha_range, scaling_factor
        - PCN_PCA: n_PC, explained_variance
        - enforce_time_interval: bool

    Returns
    -------
    pd.DataFrame
        DataFrame with trajectory analysis results including:
        - phecode pairs
        - effect sizes
        - p-values
        - significance indicators
    """
    analysis = DiseaseTrajectoryAnalysis(data)
    return analysis.analyze(*args, **kwargs)
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
        DataFrame containing comorbidity strength analysis results produced by the 'DiseaseNet.comorbidity_strength' function.
    
    binomial_test_result : pd.DataFrame
        DataFrame containing binomial test analysis results produced by the 'DiseaseNet.binomial_test' function.

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
        For other covariates specified in the 'DiseaseNet.DiseaseNetworkData.phenotype_data()' function, use their original column names.
    
    matching_n : int, default=2
        Specifies the maximum number of matched controls for each case.
    
    covariates : list, default=None
        List of phenotypic covariates to include in the model.
        By default, includes all covariates specified in the 'DiseaseNet.DiseaseNetworkData.phenotype_data()' function.
        Categorical and binary variables used for matching should not be included as covariates.
        Continuous variables used for matching can be included as covariates, but caution is advised.
        To include the required variable sex as a covariate, always use 'sex' instead of its original column name.
        For other covariates specified in the 'DiseaseNet.DiseaseNetworkData.phenotype_data()' function, use their original column names.

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
        If None, the log will be written to the temporary files directory with file prefix of DiseaseNet_.
    
    **kwargs
        Analysis option
            enforce_time_interval : bool, default=True
                If set to True, applies the specified minimum and maximum time intervals when determining the D2 outcome among individuals diagnosed with D1. 
                These time interval requirements should be defined using the DiseaseNet.DiseaseNetworkData.disease_pair() function.
    
        Additional keyword argument to define the required columns in 'comorbidity_strength_result' and 'binomial_test_result':
            phecode_d1_col : str, default='phecode_d1'
                Name of the column in 'comorbidity_strength_result' and 'binomial_test_result' that specifies the phecode identifiers for disease 1 of the disease pair.
            phecode_d2_col : str, default='phecode_d2'
                Name of the column in 'comorbidity_strength_result' and 'binomial_test_result' that specifies the phecode identifiers for disease 2 of the disease pair.
            significance_phi_col : str, default='phi_p_significance'
                Name of the column in 'comorbidity_strength_result' that indicates the significance of phi-correlation for each disease pair.
            significance_RR_col : str, default='phi_RR_significance'
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

def trajectory_multipletests(df:pd.DataFrame, correction:str='bonferroni', cutoff:float=0.05) -> pd.DataFrame:
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
