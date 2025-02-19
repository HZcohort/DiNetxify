# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 12:49:35 2024

@author: hc
"""

import DiseaseNetPy as dnp
import pandas as pd

if __name__ == "__main__":
    # data = dnp(study_design="registry")
    # data.phenotype_data(
    #     "data/data.csv",
    #     {"Participant ID": "eid",
    #     "Sex": "sex",
    #     "Index date": "date_start",
    #     "End date": "date_end"},
    #     ["bmi", "income"]
    # )
    
    # data.merge_medical_records(
    #     "data/inp.csv", 
    #     "ICD-10-WHO",
    #     {"Participant ID": "eid",
    #     "Diagnosis code": "icd10",
    #     "Date of diagnosis": "date_dignosed"}
    # )

    # phewas_result = dnp.phewas(
    #     data, 
    #     proportion_threshold=0.01, 
    #     n_process=8
    # )

    # data.disease_pair(phewas_result)

    # comorbidity_strength_result = dnp.comorbidity_strength(
    #     data, 
    #     proportion_threshold=0.01, 
    #     n_process=8
    # )

    # binomial_test_result = dnp.binomial_test(
    #     data, 
    #     comorbidity_strength_result
    # )

    # comorbidity_network_result = dnp.comorbidity_network(
    #     data, 
    #     comorbidity_strength_result, 
    #     binomial_test_result,
    #     n_process=8
    # )

    # disease_trajectory_result = dnp.disease_trajectory(
    #     data, 
    #     comorbidity_strength_result, 
    #     binomial_test_result,
    #     n_process=8
    # )
    phewas = pd.read_csv(
        "C:/Users/bovin/Desktop/data/phewas.csv"
    )

    phewas = phewas.loc[
        phewas["phewas_p_significance"]==True
    ]

    comorbidity = pd.read_csv(
        "C:/Users/bovin/Desktop/data/comorbidity.csv"
    )
    
    comorbidity = comorbidity.loc[
        comorbidity["comorbidity_p_significance"]==True
    ]

    trajectory = pd.read_csv(
        "C:/Users/bovin/Desktop/data/trajectory.csv"
    )
    
    trajectory = trajectory.loc[
        trajectory["trajectory_p_significance"]==True
    ]

    my_network = dnp.ThreeDimensionalNetwork(
        phewas,
        comorbidity,
        trajectory
    )

    my_network.threeDimension_plot(
        "/plot.html",
        45,
        15,
        "full",
        "black",
        1.0,
    )

    print(my_network._nodes_attrs)
    print(len(my_network._nodes_attrs.keys()))