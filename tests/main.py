import diseasenetpy as dnt

if __name__ =="__main__":
    # Define required columns and other covariates columns
    col_dict = {
        'Participant ID': 'ID',
        'Exposure': 'exposure',                  
        'Sex': 'sex',                           
        'Index date': 'date_start',             
        'End date': 'date_end',                 
        'Match ID': 'group_id'                   
    }
    vars_lst = ['age', 'BMI']
    
    # Initialize the data object with study design and phecode level
    data = dnt.DiseaseNetworkData(
        study_design="matched cohort",
        phecode_level=1,
        date_fmt="%Y-%m-%d"
    )

    # Load the phenotype CSV file into the data object
    data.phenotype_data(
        phenotype_data_path="data/dummy_cohort.csv",
        column_names=col_dict,
        covariates=vars_lst
    )

    # Merge with the first medical records file (CSV)
    data.merge_medical_records(
        medical_records_data_path="data/dummy_EHR_ICD9.csv",
        diagnosis_code="ICD-9-WHO",
        column_names={
            'Participant ID':'ID',
            'Diagnosis code':'diag_icd9',
            'Date of diagnosis':'dia_date'
        }
    )

    # Merge with the second medical records file (CSV)
    data.merge_medical_records(
        medical_records_data_path="data/dummy_EHR_ICD10.csv",
        diagnosis_code="ICD-10-WHO",
        column_names={
            'Participant ID':'ID',
            'Diagnosis code':'diag_icd10',
            'Date of diagnosis':'dia_date'
        }
    )

    # Describe the basic information fo phenotype data and medical data
    data.Table1(
        continuous_stat_mode="auto" 
    )

    # PheWAS Analysis
    phewas_result = dnt.phewas(
        data=data,                                             
        proportion_threshold=0.01,                            
        n_process=8,                                          
        system_exl=[                                           
            'symptoms', 
            'others', 
            'injuries & poisonings', 
            'pregnancy complications'
        ],
        covariates=['age', 'BMI'],  
        lifelines_disable=True,                                
    )

    # Generate Disease Pair for Each Individual and Update the Data Object
    data.disease_pair(
        phewas_result=phewas_result,               
        force=True
    )

    # Comorbidity Strength Estimation
    com_strength_result = dnt.comorbidity_strength(
        data=data,                                   
        proportion_threshold=0.01,                  
        n_process=8,                                         
    )

    # Binomial Test
    binomial_result = dnt.binomial_test(
        data=data,                                        
        comorbidity_strength_result=com_strength_result,                                 
        enforce_temporal_order=True,                      
    )

    # Comorbidity Network Analysis
    comorbidity_result = dnt.comorbidity_network(
        data=data,                                       
        comorbidity_strength_result=com_strength_result, 
        binomial_test_result=binomial_result,          
        n_process=8,                                   
        covariates=['age', 'BMI'],  
        method='CN',      
    )

    # Trajectory Analysis
    trajectory_result = dnt.disease_trajectory(
        data=data,                                       
        comorbidity_strength_result=com_strength_result, 
        binomial_test_result=binomial_result,           
        method='CN',                                  
        n_process=8,                                      
        matching_n=5,                                    
        enforce_time_interval=False,                     
        covariates=['age', 'BMI'], 
    )

    # visualization
    network = dnt.visualization.ThreeDimensionalNetwork(
        phewas_result=phewas_result,
        comorbidity_result=comorbidity_result,
        trajectory_result=trajectory_result
    )

    network.threeDimension_plot(
        "figure/threeDimensionPlot"
    )
    network.significant_trajectory_plot(
        "figure/significant_trajectory"
    )
    network.phewas_plot(
        "figure/phewas_plot"
    )
    network.comorbidity_network_plot(
        "figure/comorbidity_network_plot"
    )

    # save the results and data of 'DiseaseNetworkData' object
    data.save("data/data.pkl.gz")
    phewas_result.to_csv("result/phewas_result.csv")
    com_strength_result.to_csv("result/com_strength_result.csv")
    binomial_result.to_csv("result/binomial_result.csv")
    comorbidity_result.to_csv("result/comorbidity_result.csv")
    trajectory_result.to_csv("result/trajectory_result.csv")