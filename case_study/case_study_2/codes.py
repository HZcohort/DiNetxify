import pandas as pd
import time
import DiNetxify as dnt

n_cores = 40
min_interval_days = 30
max_interval_days = 5*365

path_case = r'/your/data/path'
path_result = r'/your/result/path'

if __name__ == "__main__":

    data = dnt.DiseaseNetworkData(study_design='exposed-only cohort',phecode_level=2,
                                phecode_version='1.3a')
    time1 = time.time()
    data.phenotype_data(rf'{path_case}/baseline.csv',
                        column_names={'Participant ID':'eid',
                                    'Sex': 'sex',
                                    'Index date': 'date_attend',
                                    'End date': 'date_end'},
                        covariates=['age','social_cate','BMI','income','diet','smoking','drinking','physical'])
    data.merge_medical_records(rf'{path_case}/icd10.csv',
                            diagnosis_code='ICD-10-WHO',column_names={'Participant ID': 'eid',
                                                                        'Diagnosis code': 'diag_icd10',
                                                                        'Date of diagnosis': 'dia_date'})
    data.merge_medical_records(rf'{path_case}/icd9.csv',
                            diagnosis_code='ICD-9-WHO',column_names={'Participant ID': 'eid',
                                                                        'Diagnosis code': 'diag_icd9',
                                                                        'Date of diagnosis': 'dia_date'})
    time2 = time.time()
    print(f"Data harmonization cost {(time2-time1)/60} min")
    #step 1
    phewas_result = dnt.phewas(
        data=data,                                             # DiseaseNetworkData object
        n_threshold=2500,                                       # Minimum proportion of cases to include
        n_process=n_cores,                                     # Number of parallel processes
        system_exl=[                                           # Phecode systems to exclude
            'pregnancy complications','congenital anomalies', 'symptoms', 'others', 'injuries & poisonings'],
        phecode_exl=[681.3, 681.5],
        covariates=['age', 'sex', 'social_cate','BMI','income','diet','smoking','drinking','physical'],
        lifelines_disable=True,                                # Disable lifelines for faster computation
        log_file=path_result+'/phewas.log',                # Path to log file
        correction='none'
    )
    phewas_result.to_csv(f'{path_result}/phewas.csv',index=False)
    
    #step 2
    print('Phewas: ',len(phewas_result[phewas_result['phewas_p_significance']==True]))

    data.disease_pair(phewas_result=phewas_result,min_interval_days=min_interval_days,max_interval_days=max_interval_days,force=True)

    #step 3
    com_strength_result = dnt.comorbidity_strength(
        data=data,                                     # DiseaseNetworkData object
        n_threshold=250,                   # Minimum proportion for comorbidity
        n_process=n_cores,                                   # Number of parallel processes
        log_file=path_result+'/comorbidity.log'          # Path to log file
    )
    com_strength_result.to_csv(f'{path_result}/comorbidity_strength.csv',index=False)

    com_strength_result = pd.read_csv(rf'{path_result}/comorbidity_strength.csv')
    com_strength_result = com_strength_result[(com_strength_result['phi'] > 0) & (com_strength_result['RR'] > 1)]
    com_strength_result = dnt.comorbidity_strength_multipletests(df=com_strength_result,
                                                                 correction_phi='fdr_bh',correction_RR='fdr_bh',
                                                                 cutoff_phi=0.05,cutoff_RR=0.05)
    print('Comorbidity strength: ',len(com_strength_result[(com_strength_result['RR_p_significance']==True) & 
                                                       (com_strength_result['phi_p_significance']==True)]))

    #step 4
    comorbidity_result = dnt.comorbidity_network(
        data=data,                                       # DiseaseNetworkData object
        comorbidity_strength_result=com_strength_result, # Comorbidity strength results
        n_process=n_cores,                                   # Number of parallel processes
        covariates=['age','sex','social_cate','BMI','income','diet','smoking','drinking','physical'],  # Covariates to adjust for
        method='PCN_PCA',                                   # Analysis method ('CN', 'PCN_PCA', 'RPCN')
        n_PC=20,                                     # Number of principal components to use
        log_file=path_result+'/uncond_logistic.log'         # Path to log file
    )
    comorbidity_result = dnt.comorbidity_multipletests(df=comorbidity_result,correction='fdr_bh',cutoff=0.05)
    print('Comorbidity network: ',len(comorbidity_result[comorbidity_result['comorbidity_p_significance']==True]))
    comorbidity_result.to_csv(f'{path_result}/comorbidity_PCA.csv',index=False)
    
    #step 5
    binomial_result = dnt.binomial_test(
    data=data,                                        # DiseaseNetworkData object
    comorbidity_strength_result=com_strength_result,  # Comorbidity strength results
    comorbidity_network_result=comorbidity_result,    # Comorbidity network results
    n_process=1,                                      # Number of CPU cores (1 to disable multiprocessing)
    enforce_temporal_order=True,                      # Enforce temporal order in testing
    log_file=path_result+'/binomial.log'              # Path to log file
    )
    binomial_result = dnt.binomial_multipletests(df=binomial_result, correction='fdr_bh', cutoff=0.05)
    print('Binomial: ',len(binomial_result[binomial_result['binomial_p_significance']==True]))
    binomial_result.to_csv(f'{path_result}/binomial.csv',index=False)

    #step 6
    trajectory_result = dnt.disease_trajectory(
        data=data,                                       # DiseaseNetworkData object
        comorbidity_strength_result=com_strength_result, # Comorbidity strength results
        binomial_test_result=binomial_result,           # Binomial test results
        method='PCN_PCA',                                   # Trajectory analysis method ('CN', 'PCN_PCA', 'RPCN')
        n_PC=20,
        n_process=n_cores,                                     # Number of parallel processes
        matching_var_dict={'sex': 'exact', 'age': 2, 'social_cate': 'exact'},    # Matching variables and criteria
        matching_n=5,                                    # Number of matched controls per case
        enforce_time_interval=True,                     # Enforce time interval in trajectory analysis
        covariates=['age','BMI','income','diet','smoking','drinking','physical'],  # Covariates to adjust for
        log_file=path_result+'/cond_logistic.log',            # Path to log file
        global_sampling=True)
    trajectory_result = dnt.trajectory_multipletests(df=trajectory_result,correction='fdr_bh',cutoff=0.05)
    print('Trajectory analysis: ',len(trajectory_result[trajectory_result['trajectory_p_significance']==True]))
    trajectory_result.to_csv(f'{path_result}/trajectory_PCA_match5.csv',index=False)

