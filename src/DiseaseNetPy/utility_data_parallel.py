import multiprocessing
import pandas as pd
import gc
from functools import partial
import time
from .utility import check_history_exclusion,phecode_leaf_to_root,time_first_diagnosis,n_process_check

# Define the worker function
def process_chunk(chunk, phecode_lst, disease_pair_index, node_dict, phecode_set, history_dict, diagnosis_dict, 
                    n_diagnosis_dict, phecode_info_dict, min_interval_days, max_interval_days, 
                    min_icd_num, sex_value_dict):
    
    # Initialize result dictionaries for this chunk
    chunk_ineligible_disease_dict = {}
    chunk_eligible_withdate_dict = {}
    chunk_d1d2_temporl_pair_dict = {}
    chunk_d1d2_com_pair_dict = {}
    chunk_all_diagnosis_level = {}
    
    from itertools import combinations
    
    # Process each individual in the chunk
    for id_, sex in chunk:
        temp_dineligible_list = []
        temp_dpair_temporal_lst = []
        temp_dpair_com_lst = []
        temp_deligible_dict_withdate = {}
        
        # Get individual's data
        diagnosis_ = diagnosis_dict[id_]
        n_diagnosis_ = n_diagnosis_dict[id_]
        history_ = history_dict[id_]
        
        # Generate a new all diagnosed phecode list correspond to the specified level
        all_diag_ind = set([node_dict.get(x, x) for x in history_ + list(diagnosis_.keys())]).intersection(phecode_set)
        all_diag_ind = [x for x in all_diag_ind if sum([n_diagnosis_.get(leaf, 0) for leaf in phecode_info_dict[x]['leaf_list']]) >= min_icd_num]
        chunk_all_diagnosis_level[id_] = all_diag_ind
        
        # Generate eligible disease dictionary
        for phecode in phecode_lst:
            leaf_lst = phecode_info_dict[phecode]['leaf_list']
            exl_lst = phecode_info_dict[phecode]['exclude_list']
            sex_specific = phecode_info_dict[phecode]['sex']
            
            if (check_history_exclusion(exl_lst, history_, n_diagnosis_, min_icd_num) == 0) and (sex_specific == 'Both' or sex_value_dict[sex_specific] == sex):
                date = time_first_diagnosis(leaf_lst, diagnosis_, n_diagnosis_, min_icd_num)
                if not pd.isna(date):
                    temp_deligible_dict_withdate[phecode] = date
            else:
                temp_dineligible_list.append(phecode)
        
        # Generate disease pair dictionary
        if len(temp_deligible_dict_withdate) <= 1:
            chunk_ineligible_disease_dict[id_] = temp_dineligible_list
            chunk_d1d2_temporl_pair_dict[id_] = temp_dpair_temporal_lst
            chunk_d1d2_com_pair_dict[id_] = temp_dpair_com_lst
            chunk_eligible_withdate_dict[id_] = temp_deligible_dict_withdate
        else:
            for d1, d2 in combinations(temp_deligible_dict_withdate, 2):
                date1, date2 = temp_deligible_dict_withdate[d1], temp_deligible_dict_withdate[d2]
                if abs((date1 - date2).days) > max_interval_days:
                    continue
                elif abs((date1 - date2).days) <= min_interval_days:
                    temp_dpair_com_lst.append(disease_pair_index[f'{d1}_{d2}'])  # order insensitive
                    continue
                else:
                    if date1 > date2:
                        temp_dpair_temporal_lst.append(disease_pair_index[f'{d2}_{d1}'])
                    else:
                        temp_dpair_temporal_lst.append(disease_pair_index[f'{d1}_{d2}'])
            
            # Save for the individual
            chunk_ineligible_disease_dict[id_] = temp_dineligible_list
            chunk_d1d2_temporl_pair_dict[id_] = temp_dpair_temporal_lst
            chunk_d1d2_com_pair_dict[id_] = temp_dpair_com_lst
            chunk_eligible_withdate_dict[id_] = temp_deligible_dict_withdate
    
    # Return this chunk's results
    return {
        'ineligible_disease': chunk_ineligible_disease_dict,
        'eligible_disease_withdate': chunk_eligible_withdate_dict,
        'd1d2_temporal_pair': chunk_d1d2_temporl_pair_dict,
        'd1d2_com_pair': chunk_d1d2_com_pair_dict,
        'all_diagnosis_level': chunk_all_diagnosis_level
    }

def parallel_d1d2_from_diagnosis_history(df:pd.DataFrame, id_col:str, sex_col:str, sex_value_dict:dict, 
                                phecode_lst:list, disease_pair_index:dict, history_dict:dict, diagnosis_dict:dict, n_diagnosis_dict:dict,
                                phecode_info_dict:dict, min_interval_days:int, max_interval_days:int, min_icd_num:int,
                                n_process:int) -> dict:
    """
    Parallel version of d1d2_from_diagnosis_history that constructs d1->d2 disease pairs
    for each individual from a list of significant phecodes using multiprocessing.
    
    Parameters
    ----------
        df : dataframe of phenotype data, contains at id and sex columns
        id_col : id column in the df
        sex_col : sex column in the df
        sex_value_dict: dictionary for coding 'Female' and 'Male' in the phenotype data
        phecode_lst : list of significant phecodes
        disease_pair_index : dictionary containing disease pair index
        history_dict : dictionary containing medical records history
        diagnosis_dict : dictionary containing diagnosis and date
        n_diagnosis_dict : number of phecode occurence
        phecode_info_dict : phecode information
        min_interval_days : minimum interval required for d1-d2 disease pair construction
        max_interval_days : maximum interval allowed for d1-d2 disease pair construction
        min_icd_num : The minimum number of ICD codes mapping to a specific phecode required for the phecode to be considered valid.
        n_process : Number of processes to use.

    Returns
    -------
        D1->D2 dictionary.
    """
    # Check multiprocessing parameters
    n_process,start_mehtod = n_process_check(n_process,'Disease pairs construction')
    
    # Generate a dictionary mapping leaf phecode to root phecode (shared across all individuals)
    node_dict = phecode_leaf_to_root(phecode_info_dict)
    phecode_set = set(phecode_info_dict.keys())
    
    # Extract IDs and sex values
    id_sex_pairs = df[[id_col, sex_col]].values
    chunk_size = int(len(id_sex_pairs)//n_process/2)
    
    # Divide the work into chunks
    chunks = [id_sex_pairs[i:i + chunk_size] for i in range(0, len(id_sex_pairs), chunk_size)]
    
    # Create a partial function with all the fixed parameters
    worker_func = partial(
        process_chunk,
        phecode_lst=phecode_lst,
        disease_pair_index=disease_pair_index,
        node_dict=node_dict,
        phecode_set=phecode_set,
        history_dict=history_dict,
        diagnosis_dict=diagnosis_dict,
        n_diagnosis_dict=n_diagnosis_dict,
        phecode_info_dict=phecode_info_dict,
        min_interval_days=min_interval_days,
        max_interval_days=max_interval_days,
        min_icd_num=min_icd_num,
        sex_value_dict=sex_value_dict
    )
    time_0 = time.time()
    # Set up multiprocessing pool with context manager
    with multiprocessing.get_context(start_mehtod).Pool(processes=n_process) as pool:
        # Map the worker function to the chunks
        results = pool.map(worker_func, chunks)
    time_1 = time.time()
    print(f'Parallel processing time: {time_1-time_0} seconds')
    # Merge results from all chunks
    ineligible_disease_dict = {}
    eligible_withdate_dict = {}
    d1d2_temporl_pair_dict = {}
    d1d2_com_pair_dict = {}
    all_diagnosis_level = {}
    
    while results:
        result = results.pop()  # removes and returns the last element
        ineligible_disease_dict.update(result['ineligible_disease'])
        eligible_withdate_dict.update(result['eligible_disease_withdate'])
        d1d2_temporl_pair_dict.update(result['d1d2_temporal_pair'])
        d1d2_com_pair_dict.update(result['d1d2_com_pair'])
        all_diagnosis_level.update(result['all_diagnosis_level'])
        gc.collect()
    
    # Final dictionary
    return {
        'disease_pair_index' : disease_pair_index,
        'ineligible_disease': ineligible_disease_dict,
        'eligible_disease_withdate': eligible_withdate_dict,
        'd1d2_temporal_pair': d1d2_temporl_pair_dict,
        'd1d2_com_pair': d1d2_com_pair_dict,
        'all_diagnosis_level': all_diagnosis_level
    }