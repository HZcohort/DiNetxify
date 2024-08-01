# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 23:48:05 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

from .data_management import DiseaseNetworkData
#from .analysis import phewas, comorbidity_analysis, trajectory_analysis
#from .regression import unconditional_logistic, conditional_logistic
#from .visualization import create_3d_network

"""
import DiseaseNet as dnt

#create data objective
full_data = dnt.DiseaseNetworkData(study_design='matched cohort',phecode_level=1,date_fmt='%Y-%m-%d')
full_data.phenotype_data('phe.csv', column_names={'Participant ID': 'eid',
                                                  'Exposure': 'status',
                                                  'Index date': 'index_date',
                                                  'End date': 'final_date',
                                                  'Matching identifier': 'group'}, 
                         covariates=['age', 'sex', 'BMI'])
full_data.merge_medical_records('icd10.csv', diagnosis_code='ICD-10-WHO', column_names={'Participant ID': 'eid',
                                                                                         'Diagnosis code': 'ICD10',
                                                                                         'Date of diagnosis': 'date'})
full_data.merge_medical_records(r'icd9.csv', diagnosis_code='ICD-9-WHO', column_names={'Participant ID': 'eid',
                                                                                         'Diagnosis code': 'ICD9',
                                                                                         'Date of diagnosis': 'date'})
full_data.save(r'./prefix') #saved to prefix.npz

#phewas
phewas_level1 = dnt.phewas(full_data,n_cpus=10,adjustment='FDR')

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














