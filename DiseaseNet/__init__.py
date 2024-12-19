# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 23:48:05 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

from .data_management import DiseaseNetworkData
from .analysis import phewas,phewas_multipletests, comorbidity_strength, comorbidity_strength_multipletests
from .analysis import binomial_test, binomial_multipletests, comorbidity_network, comorbidity_multipletests
from .analysis import disease_trajectory, trajectory_multipletests
#from .analysis import phewas, comorbidity_analysis, trajectory_analysis
#from .regression import unconditional_logistic, conditional_logistic
#from .visualization import create_3d_network


"""
import DiseaseNet as dnt

#create data objective
col_dict = {'Participant ID': 'new_index',
            'Exposure': 'outcome',
            'Sex':'sex',
            'Index date': 'date_start',
            'End date': 'time_end',
            'Matching identifier': 'match_2'}
vars_lst = ['age','social','BMI','smoking', 'drinking']

data = dnt.DiseaseNetworkData(study_design='matched cohort',phecode_level=1,date_fmt='%Y-%m-%d')
data.phenotype_data(phenotype_data_path='phenotype.csv',
                         column_names=col_dict, covariates=vars_lst)

data.merge_medical_records('inp_1.csv',diagnosis_code='ICD-10-WHO',
                           column_names={'Participant ID': 'eid',
                                         'Diagnosis code': 'diag_icd10',
                                         'Date of diagnosis': 'date'})

data.merge_medical_records('inp_2.csv',diagnosis_code='ICD-10-WHO',
                           column_names={'Participant ID': 'eid',
                                         'Diagnosis code': 'diag_icd10',
                                         'Date of diagnosis': 'date'})

#save
data.save('module_test\dep')

#load
data = dnt.DiseaseNetworkData(study_design='matched cohort',phecode_level=1,date_fmt='%Y-%m-%d')
data.load('dep.npy')
data.modify_phecode_level(2)

#phewas
phewas_result = dnt.phewas(data,n_threshold=200,n_cpus=5,system_inc=['digestive'],sex_adjustment=False,
                           lifelines_disable=True,log_file='phewas.log')
phewas_result = phewas_result[phewas_result['phewas_coef']>=0]
phewas_result = dnt.phewas_multipletests(phewas_result,correction='fdr_bh')

#generate disease pairs for only exposed group
data.disease_pair(phewas_result,min_interval_days=30,max_interval_days=365*5) #min_interval_days and max_interval_days added
data.save('module_test\dep_withtra') #save again new data

#comorbidity strength estimation
com_strength_result = dnt.comorbidity_strength(data,proportion_threshold=0.001,n_cpus=5,
                                               log_file='com.log')
com_strength_result = com_strength_result[(com_strength_result['phi']>=0) & (com_strength_result['RR']>=1)] #further filtering as required
com_strength_result = dnt.comorbidity_strength_multipletests(com_strength_result,correction_phi='fdr_bh',correction_RR='fdr_bh')


#binomial test
binomial_result = dnt.binomial_test(data, com_strength_result,n_cpus=1,enforce_temporal_order=True,
                                    log_file='com.log')
binomial_result = dnt.binomial_multipletests(binomial_result,correction='fdr_bh')

#comorbidity network analysis - traditional method
comorbidity_result = dnt.comorbidity_network(data, com_strength_result, binomial_result, n_cpus=6, method='CN',log_file='com.log')
#comorbidity network analysis - partial correlation network method
comorbidity_result = dnt.comorbidity_network(data, com_strength_result, binomial_result, n_cpus=6, method='RPCN',log_file='com.log')
#comorbidity network analysis - new PCA method
comorbidity_result = dnt.comorbidity_network(data, com_strength_result, binomial_result, n_cpus=6, method='PCN_PCA',n_PC=15,
                                             log_file='com.log')
comorbidity_result = dnt.comorbidity_multipletests(comorbidity_result,correction='fdr_bh')

#trajectory analysis  - traditional method
trajectory_result = dnt.disease_trajectory(data, com_strength_result, binomial_result, method='CN',n_cpus=6,
                                           matching_var_dict={'age':2,'sex':'exact'},matching_n=5,
                                           covariates=['social','BMI','smoking','drinking'],
                                           log_file='com.log')
#trajectory analysis  - partial correlation network method
trajectory_result = dnt.disease_trajectory(data, com_strength_result, binomial_result, method='RPCN',n_cpus=6,
                                           matching_var_dict={'age':2,'sex':'exact'},matching_n=5,
                                           covariates=['social','BMI','smoking','drinking'],
                                           log_file='com.log')
#trajectory analysis  - new PCA method
trajectory_result = dnt.disease_trajectory(data, com_strength_result, binomial_result, method='PCN_PCA',n_PC=15,
                                           n_cpus=6,matching_var_dict={'age':2,'sex':'exact'},matching_n=5,
                                           covariates=['social','BMI','smoking','drinking'],
                                           log_file='com.log')
trajectory_result = dnt.trajectory_multipletests(trajectory_result,correction='fdr_bh')



#3D visulization
dnt.create_3d_network(unconditional_logistic_result=uncond,conditional_logistic_result=cond,
                      save_html='./3d.html')
test again

"""

















