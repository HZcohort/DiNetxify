# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 12:49:35 2024

@author: hc
"""

from DiseaseNetPy import DiseaseNetworkData as dnd

if __name__ == "__main__":
    data = dnd(study_design="registry")
    data.phenotype_data(
        "data/data.csv",
        {"Participant ID": "eid",
        "Sex": "sex",
        "Index date": "date_start",
        "End date": "date_end"},
        ["bmi", "income"]
    )
    
    data.merge_medical_records(
        "data/inp.csv", 
        "ICD-10-WHO",
        {"Participant ID": "eid",
        "Diagnosis code": "icd10",
        "Date of diagnosis": "date_dignosed"}
    )

    phewas_result = dnd.phewas(
        data, 
        proportion_threshold=0.01, 
        n_process=8
    )

    data.disease_pair(phewas_result)

    comorbidity_strength_result = dnd.comorbidity_strength(
        data, 
        proportion_threshold=0.01, 
        n_process=8
    )

    binomial_test_result = dnd.binomial_test(
        data, 
        comorbidity_strength_result
    )

    comorbidity_network_result = dnd.comorbidity_network(
        data, 
        comorbidity_strength_result, 
        binomial_test_result,
        n_process=8
    )

    disease_trajectory_result = dnd.disease_trajectory(
        data, 
        comorbidity_strength_result, 
        binomial_test_result,
        n_process=8
    )

    my_network = dnd.ThreeDimensionalDiseaseNetwork(
        phewas_result,
        comorbidity_network_result,
        disease_trajectory_result
    )