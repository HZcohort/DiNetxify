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
        "C:/Users/bovin/Desktop/phewas_summary_L1L2.csv"
    )

    # phewas = phewas.loc[
    #     phewas["phewas_p_significance"]==True
    # ]

    comorbidity = pd.read_csv(
        "C:/Users/bovin/Desktop/com_summary_0.01.csv"
    )
    
    # comorbidity = comorbidity.loc[
    #     comorbidity["comorbidity_p_significance"]==True
    # ]

    trajectory = pd.read_csv(
        "C:/Users/bovin/Desktop/tra_summary_0.01.csv"
    )
    
    # trajectory = trajectory.loc[
    #     trajectory["trajectory_p_significance"]==True
    # ]

    my_network = dnp.ThreeDimensionalNetwork(
        phewas,
        comorbidity,
        trajectory,
        296.2,
        (0,0,0),
        3,
        0.002,
        source="d1",
        target="d2"
    )

    # my_network.threeDimension_plot(
    #     "/compact_plot.html",
    #     90,
    #     5,
    #     "compact",
    #     "black",
    #     1.0,
    # )

    my_network.threeDimension_plot(
        "/full_plot.html",
        90,
        5,
        "full",
        "black",
        1.0,
        0.3,
        layer_distance=40
    )

    # my_network.threeDimension_plot(
    #     "/half_plot.html",
    #     90,
    #     5,
    #     "half",
    #     "black",
    #     1.0,
    # )

    # my_network.comorbidity_network_plot(
    #     "/comorbidity.png",
    #     90,
    #     5
    # )

    # my_network.significant_trajectory_plot(
    #     "/sig_tra",
    #     "black"
    # )