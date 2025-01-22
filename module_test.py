# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 12:49:35 2024

@author: hc
"""

# import pandas as pd
# from DiseaseNet.visualization import ThreeDimensionalDiseaseNetwork as td

# phewas_result = pd.read_csv("C:/Users/bovin/Desktop/phewas_result.csv")
# comorbidity = pd.read_csv("C:/Users/bovin/Desktop/comorbidity_network_result.csv")
# comorbidity = comorbidity.loc[comorbidity["comorbidity_p_significance"]==True]
# trajetory = pd.read_csv("C:/Users/bovin/Desktop/disease_trajectory_result.csv")
# trajetory = trajetory.loc[trajetory["trajectory_p_significance"]==True]

# figure = td(comorbidity, trajetory, 999, (0, 0, 50), 3, phewas_result)
# figure.plot_plotly(55, 15, "full", "black", 1, "plot.html")

# import random
# import pandas as pd

# eid = [x for x in range(10000)]
# sex = random.choices(range(2), k=10000)
# birth_year = random.choices(range(1956, 1967), k=10000)
# birth_month = random.choices(range(1, 13), k=10000)
# birth_day = random.choices(range(1, 28), k=10000)
# income = random.choices(range(3), k=10000)
# bmi = random.choices(range(4), k=10000)
# date_start = [f"{random.choices(range(1967,1975), k=1)[0]}-{random.choices(range(1,13), k=1)[0]}-{random.choices(range(1,28), k=1)[0]}" for i in range(10000)]
# date_end = [f"{random.choices(range(1976,1987), k=1)[0]}-{random.choices(range(1,13), k=1)[0]}-{random.choices(range(1,28), k=1)[0]}" for i in range(10000)]
# data = pd.DataFrame({"eid":eid, 
#                      "sex":sex, 
#                      "birth_year":birth_year, 
#                      "birth_month":birth_month, 
#                      "birth_day":birth_day, 
#                      "income":income, 
#                      "bmi":bmi, 
#                      "date_start":date_start,
#                      "date_end":date_end}, 
#                     columns=["eid", "sex", "birth_year", "birth_month", "birth_day", "income", "bmi", "date_start", "date_end"])
# data.to_csv("C:/Users/bovin/Desktop/data.csv")
# date_dignosed = [f"{random.choices(range(1976,1978), k=1)[0]}-{random.choices(range(1,13), k=1)[0]}-{random.choices(range(1,28), k=1)[0]}" for i in range(10000)]
# disease_icd10 = random.choices(["C809", "C800", "N390", "C859", "Y280", "C221", "C61", "C260"], k=10000)
# inp_data = pd.DataFrame({"eid":eid,
#                          "icd10":disease_icd10,
#                          "date_dignosed":date_dignosed},
#                          columns=["eid", "icd10", "date_dignosed"]).to_csv("C:/Users/bovin/Desktop/inp.csv")

# from DiseaseNet.data_management import DiseaseNetworkData as dnd
# from DiseaseNet.analysis import phewas, comorbidity_strength, binomial_test, comorbidity_network, disease_trajectory
from DiseaseNet.visualization import ThreeDimensionalDiseaseNetwork as tddn
import pandas as pd

if __name__ == "__main__":
    # # 加载数据
    # data = dnd(study_design="registry")
    # data.phenotype_data("C:/Users/bovin/Desktop/data.csv",
    #                     {"Participant ID": "eid",
    #                     "Sex": "sex",
    #                     "Index date": "date_start",
    #                     "End date": "date_end"},
    #                     ["bmi", "income"])
    
    # data.merge_medical_records("C:/Users/bovin/Desktop/inp.csv", 
    #                         "ICD-10-WHO",
    #                         {"Participant ID": "eid",
    #                         "Diagnosis code": "icd10",
    #                         "Date of diagnosis": "date_dignosed"})
    # # 数据分析
    # phewas_result = phewas(data, proportion_threshold=0.1, n_process=8)
    # data.disease_pair(phewas_result)
    # comorbidity_strength_result = comorbidity_strength(data, proportion_threshold=0.1)
    # binomial_test_result = binomial_test(data, comorbidity_strength_result)
    # comorbidity_network_result = comorbidity_network(data, comorbidity_strength_result, binomial_test_result)
    # disease_trajectory_result = disease_trajectory(data, comorbidity_strength_result, binomial_test_result)

    # phewas_result.to_csv("C:/Users/bovin/Desktop/s.csv")
    # data.disease_pair(phewas_result)
    # comorbidity_strength_result = comorbidity_strength(data, proportion_threshold=0.1)
    # comorbidity_strength_result.to_csv("C:/Users/bovin/Desktop/sa.csv")
    # binomial_test_result = binomial_test(data, comorbidity_strength_result)

    phewas_result = pd.read_csv("C:/Users/bovin/Desktop/phewas_result.csv")
    disease_trajecty_result = pd.read_csv("C:/Users/bovin/Desktop/disease_trajectory_result.csv")
    comorbidity_network_result = pd.read_csv("C:/Users/bovin/Desktop/comorbidity_network_result.csv")
    my_network = tddn(comorbidity_network_result, disease_trajecty_result, phewas_result)
    my_network.plot_3d(50, 15, "full", "black", 0.5, "C:/Users/bovin/Desktop/my_network")