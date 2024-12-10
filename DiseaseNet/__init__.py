# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 23:48:05 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

from .data_management import DiseaseNetworkData
from .analysis import phewas,phewas_multipletests
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

#phewas
phewas_result = dnt.phewas(data,n_threshold=200,n_cpus=5,system_inc=['digestive'],sex_adjustment=False,
                           log_file='phewas.log')
dnt.phewas_multipletests(phewas_result,adjustment='fdr_bh')

#generate trajectory for only exposed group
full_data.trajectory(phewas_result=phewas_level1)

#comorbidity strength estimation
com = dnt.comorbidity_analysis(full_data,phewas_result=phewas_level1,n_cpus=10,adjustment='FDR')

#trajectory analysis
tra = dnt.trajectory_analysis(full_data,comorbidity_result=com,n_cpus=10,adjustment='FDR')

#unconditional logistic regression
uncond = dnt.unconditional_logistic(full_data,comorbidity_result=com,n_cpus=10,adjustment='FDR',coe=0.1)

#conditional logistic regression
cond = dnt.conditional_logistic(full_data,trajectory_result=tra,n_cpus=10,adjustment='FDR',coe=0.1)

#3D visulization
dnt.create_3d_network(unconditional_logistic_result=uncond,conditional_logistic_result=cond,
                      save_html='./3d.html')

"""












