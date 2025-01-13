# -*- coding: utf-8 -*-
"""
Created on Wed Dec 25 18:08:19 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

import pandas as pd
import numpy as np
import time
from .data_management import DiseaseNetworkData
from .utility import log_file_detect,filter_phecodes,threshold_check,n_cpus_check,correction_method_check,states_p_adjust
from .utility import check_kwargs_com_tra,covariates_check,matching_var_check


class DiseaseNetworkAnalysis:
    def __init__(self, data):
        if not isinstance(data, DiseaseNetworkData):
            raise TypeError("The input 'data' must be a DiseaseNetworkData object.")
        self.data = data
        self.phecode_info = data.phecode_info

    def _check_attributes(self, attrs):
        for attr in attrs:
            if getattr(self.data, attr) is None:
                raise ValueError(f"Attribute '{attr}' is empty.")

    def _check_covariates(self, covariates):
        return covariates_check(covariates, self.data.get_attribute('phenotype_info'))

    def _check_threshold(self, proportion_threshold, n_threshold):
        n_exposed = self.data.get_attribute('phenotype_statistics')['n_exposed']
        return threshold_check(proportion_threshold, n_threshold, n_exposed)

    def _check_n_cpus(self, n_cpus, analysis_type):
        n_cpus_check(n_cpus, analysis_type)
        if n_cpus > 1:
            import multiprocessing

    def _check_correction_method(self, correction, cutoff):
        correction_method_check(correction, cutoff)

    def _check_log_file(self, log_file, analysis_type):
        log_file_final, message = log_file_detect(log_file, analysis_type)
        print(message)
        return log_file_final

    def _generate_result_df(self, result_all, columns):
        max_columns = max([len(x) for x in result_all])
        columns_selected = columns[0:max_columns]
        return pd.DataFrame(result_all, columns=columns_selected)

    def phewas(self, covariates=None, proportion_threshold=None, n_threshold=None, n_cpus=1, correction='bonferroni', cutoff=0.05, system_inc=None, system_exl=None, phecode_inc=None, phecode_exl=None, log_file=None, lifelines_disable=False):
        from .cox import cox_conditional,cox_unconditional
        
        self._check_attributes(['phenotype_df', 'diagnosis', 'history'])
        covariates = self._check_covariates(covariates)
        n_threshold = self._check_threshold(proportion_threshold, n_threshold)
        self._check_n_cpus(n_cpus, 'PheWAS')
        self._check_correction_method(correction, cutoff)
        log_file_final = self._check_log_file(log_file, 'phewas')

        phecode_lst_all = filter_phecodes(self.phecode_info, system_inc, system_exl, phecode_inc, phecode_exl)
        print(f'A total of {len(phecode_lst_all)} phecodes included in the PheWAS analysis.')

        time_start = time.time()
        result_all = []
        if self.data.study_design == 'matched cohort':
            if n_cpus == 1:
                for phecode in phecode_lst_all:
                    result_all.append(cox_conditional(self.data, n_threshold, phecode, covariates, log_file_final, lifelines_disable))
            elif n_cpus > 1:
                with multiprocessing.get_context(start_mehtod).Pool(n_cpus) as p:
                    parameters_all = [[self.data, n_threshold, phecode, covariates, log_file_final, lifelines_disable] for phecode in phecode_lst_all]
                    result_all = p.starmap(cox_conditional, parameters_all)
        elif self.data.study_design == 'cohort':
            if n_cpus == 1:
                for phecode in phecode_lst_all:
                    result_all.append(cox_unconditional(self.data, n_threshold, phecode, covariates, log_file_final, lifelines_disable))
            elif n_cpus > 1:
                with multiprocessing.get_context(start_mehtod).Pool(n_cpus) as p:
                    parameters_all = [[self.data, n_threshold, phecode, covariates, log_file_final, lifelines_disable] for phecode in phecode_lst_all]
                    result_all = p.starmap(cox_unconditional, parameters_all)

        time_end = time.time()
        time_spent = (time_end - time_start) / 60
        print(f'PheWAS analysis finished (elapsed {time_spent:.1f} mins)')

        columns = ['phecode', 'disease', 'system', 'sex', 'N_cases_exposed', 'describe', 'exposed_group', 'unexposed_group', 'phewas_coef', 'phewas_se', 'phewas_p']
        phewas_df = self._generate_result_df(result_all, columns)
        phewas_df = phewas_multipletests(phewas_df, correction=correction, cutoff=cutoff)

        return phewas_df

    def comorbidity_strength(self, proportion_threshold=None, n_threshold=None, n_cpus=1, log_file=None, correction_phi='bonferroni', cutoff_phi=0.05, correction_RR='bonferroni', cutoff_RR=0.05):
        self._check_attributes(['trajectory'])
        n_threshold = self._check_threshold(proportion_threshold, n_threshold)
        self._check_n_cpus(n_cpus, 'comorbidity_strength')
        self._check_correction_method(correction_phi, cutoff_phi)
        self._check_correction_method(correction_RR, cutoff_RR)
        log_file_final = self._check_log_file(log_file, 'com_strength')

        phecodes_sig = self.data.get_attribute('significant_phecodes')
        d1d2_pair_lst = [(d1, d2, 'Disease pair with different sex specificity' if {self.phecode_info[d1]['sex'], self.phecode_info[d2]['sex']} == {'Female', 'Male'} else None) for d1, d2 in combinations(phecodes_sig, 2)]

        time_start = time.time()
        result_all = []
        if n_cpus == 1:
            for d1, d2, describe in d1d2_pair_lst:
                result_all.append(com_phi_rr(self.data.trajectory, d1, d2, describe, n_threshold, log_file_final))
        elif n_cpus > 1:
            with multiprocessing.get_context(start_mehtod).Pool(n_cpus) as p:
                parameters_all = [[self.data.trajectory, d1, d2, describe, n_threshold, log_file_final] for d1, d2, describe in d1d2_pair_lst]
                result_all = p.starmap(com_phi_rr, parameters_all)

        time_end = time.time()
        time_spent = (time_end - time_start) / 60
        print(f'Comorbidity strength estimation finished (elapsed {time_spent:.1f} mins)')

        columns = ['phecode_d1', 'phecode_d2', 'name_disease_pair', 'N_exposed', 'n_total', 'n_d1d2_diagnosis', 'n_d1_diagnosis', 'n_d2_diagnosis', 'n_d1d2_nontemporal', 'n_d1d2_temporal', 'n_d2d1_temporal', 'n_d1d2_pair', 'description', 'phi', 'phi_theta', 'phi_p', 'RR', 'RR_theta', 'RR_p']
        com_df = self._generate_result_df(result_all, columns)
        for d in ['d1', 'd2']:
            com_df[f'disease_{d}'] = com_df[f'phecode_{d}'].apply(lambda x: self.phecode_info[x]['phenotype'])
            com_df[f'system_{d}'] = com_df[f'phecode_{d}'].apply(lambda x: self.phecode_info[x]['category'])
            com_df[f'sex_{d}'] = com_df[f'phecode_{d}'].apply(lambda x: self.phecode_info[x]['sex'])

        com_df = comorbidity_strength_multipletests(com_df, correction_phi=correction_phi, cutoff_phi=cutoff_phi, correction_RR=correction_RR, cutoff_RR=cutoff_RR)

        return com_df

    def binomial_test(self, comorbidity_strength_result, n_cpus=1, log_file=None, correction='bonferroni', cutoff=0.05, enforce_temporal_order=False, **kwargs):
        self._check_attributes(['trajectory'])
        if not isinstance(comorbidity_strength_result, pd.DataFrame):
            raise TypeError("The provided input 'comorbidity_strength_result' must be a pandas DataFrame.")
        self._check_n_cpus(n_cpus, 'binomial_test')
        self._check_correction_method(correction, cutoff)
        log_file_final = self._check_log_file(log_file, 'binomial_test')

        phecode_d1_col = kwargs.get('phecode_d1_col', 'phecode_d1')
        phecode_d2_col = kwargs.get('phecode_d2_col', 'phecode_d2')
        n_nontemporal_col = kwargs.get('n_nontemporal_col', 'n_d1d2_nontemporal')
        n_temporal_d1d2_col = kwargs.get('n_temporal_d1d2_col', 'n_d1d2_temporal')
        n_temporal_d2d1_col = kwargs.get('n_temporal_d2d1_col', 'n_d2d1_temporal')
        significance_phi_col = kwargs.get('phi_p_significance', 'phi_p_significance')
        significance_RR_col = kwargs.get('phi_p_significance', 'phi_p_significance')

        comorbidity_sig = comorbidity_strength_result[(comorbidity_strength_result[significance_phi_col] == True) & (comorbidity_strength_result[significance_RR_col] == True)]
        if len(comorbidity_sig) == 0:
            raise ValueError("No disease pair remained after filtering on significance of phi-correlation and RR.")

        time_start = time.time()
        result_all = []
        if n_cpus == 1:
            for d1, d2, n_com, n_d1d2, n_d2d1 in comorbidity_sig[[phecode_d1_col, phecode_d2_col, n_nontemporal_col, n_temporal_d1d2_col, n_temporal_d2d1_col]].values:
                result_all.append(binomial(d1, d2, n_com, n_d1d2, n_d2d1, enforce_temporal_order, log_file_final))
        elif n_cpus > 1:
            raise ValueError('Multiprocessing has been disabled for this analysis.')

        time_end = time.time()
        time_spent = (time_end - time_start) / 60
        print(f'Binomial test finished (elapsed {time_spent:.1f} mins)')

        columns = ['phecode_d1', 'phecode_d2', 'name_disease_pair', 'n_d1d2_nontemporal', 'n_d1d2_temporal', 'n_d2d1_temporal', 'binomial_p', 'binomial_proportion', 'binomial_proportion_ci']
        bino_df = self._generate_result_df(result_all, columns)
        for d in ['d1', 'd2']:
            bino_df[f'disease_{d}'] = bino_df[f'phecode_{d}'].apply(lambda x: self.phecode_info[x]['phenotype'])
            bino_df[f'system_{d}'] = bino_df[f'phecode_{d}'].apply(lambda x: self.phecode_info[x]['category'])
            bino_df[f'sex_{d}'] = bino_df[f'phecode_{d}'].apply(lambda x: self.phecode_info[x]['sex'])

        bino_df = binomial_multipletests(bino_df, correction=correction, cutoff=cutoff)

        return bino_df

    # Similar methods for comorbidity_network and disease_trajectory can be added here following the same pattern.