# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 10:07:40 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

import pandas as pd
import numpy as np
from datetime import datetime
from tqdm import tqdm

def decimal_to_short(code:float) -> str:
    """

    Convert an ICD9 code from decimal format to short format.

    Parameters:
    ----------

    """
    parts = str(code).split(".")
    parts[0] = parts[0].zfill(3)
    return "".join(parts)

def phecode_leaf_to_root(phecode_dict:dict):
    """
    Generate a new dictionary mapping all the leaf phecodes to its node phecode

    Parameters:
        phecode_dict : dict, dictionary of phecode information contained in the diseasenetwork data.

    Returns:
        mapping dictionary
    """
    new_dict = {}
    for phecode in phecode_dict:
        leaf_lst = phecode_dict[phecode]['leaf_list']
        for leaf in leaf_lst:
            new_dict[leaf] = phecode
    return new_dict

def read_check_csv(path_file:str, 
                   cols_check:list, 
                   date_cols:list, 
                   date_fmt:str, 
                   separator_lst:list=[',', '\t'],
                   return_df:bool=True):
    """
    
    Reads (or not) a CSV or TSV file into a pandas DataFrame after performing several checks.
    
    Parameters:
    ----------
    path_file : str
        Path to the file to be read.
            
    cols_check : list
        List of columns that must be present in the file.

    date_cols : list
        List of columns that contain dates.

    date_fmt : str
        Date format string, compatible with datetime.strptime.
        
    separator_lst : list
        List of seperators to try, defaul is csv or tab seperator.
    
    Returns:
    -------
    pd.DataFrame (if return_df=True)
        A DataFrame containing the data from the file if all checks pass.
    string
        Seperator used in the file, either TAB or CSV seperator.
    
    Raises:
    -------
    ValueError: If any of the following conditions are met:
        - The file cannot be read with the specified separators.
        - Date columns cannot be converted to the specified date format.
    
    """
    # Try reading the file with comma and tab separators
    for sep in separator_lst:
        try:
            # Attempt to read the first 50 rows
            df = pd.read_csv(path_file, sep=sep, nrows=50) 
            # Check if all required columns are present
            if not all(col in df.columns for col in cols_check):
                cols_not_in = [col for col in cols_check if col not in df.columns]
                raise ValueError(f'Tried with seperator "{sep}", but the required columns {cols_not_in} were not found')
            # Check for missing values in the required columns
            if df[cols_check].isnull().all().all():
                print(df[cols_check].isnull().all())
                raise ValueError("All values in the first 50 rows are missing; unable to proceed. Please verify your input data.")
            # Check date columns with the specified format
            for date_col in date_cols:
                if not pd.to_datetime(df[date_col], format=date_fmt, errors='raise').notnull().all():
                    raise ValueError(f"Date format error in column {date_col}")
            # Read and return the full dataset if needed
            if return_df:
                full_df = pd.read_csv(path_file, sep=sep, usecols=cols_check)
                return full_df,sep
            else:
                return sep
        except (pd.errors.ParserError, ValueError) as e:
            print(f"Error encountered: {e}")  # Print the exception message
            continue  # Try the next separator
    raise ValueError("File cannot be read with the given specifications or data format is incorrect, check the error information above.")


def diff_date_years(dataframe,date1:str,date2:str,rounding:int=4):
    """
    
    Calculate the difference of two dates (datetime format) in year in a dataframe.

    Parameters
    ----------
    dataframe : pd.DataFrame
        The DataFrame that contains two dates
    
    date1 : str
        Start date column.

    date2 : str
        End date column.

    Returns
    -------
    pd.series
        Difference between the two dates in year (4 digits).

    """
    days_per_year = 365.25
    df = dataframe.copy()
    df['diff'] = df[date2] - df[date1]
    df['diff'] = df['diff'].apply(lambda x: round(x.days/days_per_year,rounding))
    return df['diff']
    

def convert_column(dataframe, column:str):
    """
    
    The convert_column function is designed to analyze a specified column in a 
    given DataFrame, detect its data type (binary, continuous, or categorical),
    and convert it accordingly. 
    
    Parameters
    ----------
    dataframe : pd.DataFrame
        The DataFrame containing the column to be analyzed and converted. 
        This DataFrame is not modified in-place.
        
    column : str
        The name of the column in the DataFrame to analyze and convert.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the converted column. 

    str
        A string indicating the detected and treated data type of the column: 
        'binary', 'continuous', or 'categorical'.

    """
    # Detect the type of the column
    df = dataframe.copy()
    #deal with na values
    df[column].replace('', np.nan, inplace=True)
    unique_vals = df[column].dropna().unique()
    n_unique_vals = len(unique_vals)
    #rename the variable to '_{var}'
    new_column = f'_{column}'
    df.rename(columns={column:new_column},inplace=True)

    if n_unique_vals == 2:
        # Check if binary
        if all(isinstance(x, (int, np.integer, float, np.floating)) and x in {0, 1} for x in unique_vals):
            # Already binary and no missing values
            print(f"'{column}' is already binary variable")
            return df[[new_column]],'binary'
        elif df[new_column].isna().any():
            # Binary with missing values
            print(f"Warning: '{column}' is binary variable but has missing values. Treating as categorical variable.")
            return pd.get_dummies(df[new_column], prefix=column, dummy_na=True, drop_first=True),'binary'
        else:
            # Convert to binary
            mapping = {value: idx for idx, value in enumerate(unique_vals)}
            print(f"'{column}' converted to binary variable with mapping: {mapping}")
            df[new_column] = df[new_column].map(mapping)
            return df[[new_column]],'binary'
    elif (df[new_column].dtype == float or df[new_column].dtype == int or df[new_column].dtype == object) and n_unique_vals>=10:
        # Treat as continuous
        n_missing_before = df[new_column].isna().sum()
        df[new_column] = pd.to_numeric(df[new_column], errors='coerce')
        n_missing_after = df[new_column].isna().sum()
        if n_missing_after > n_missing_before:
            print(f"Warning: '{column}' is a continuous variable and {n_missing_after-n_missing_before} set to missing after conversion (total {n_missing_after} missing).")
        else:
            print(f"'{column}' is a continuous variable (total {n_missing_after} missing).")
        return df[[new_column]],'continuous'
    else:
        # Treat as categorical
        if n_unique_vals <= 5:
            print(f"'{column}' is treated as a categorical variable.")
        else:
            print(f"Warning: '{column}' is treated as a categorical variable, but there are {n_unique_vals} unique values.")
        if df[new_column].isna().any():
            dummies = pd.get_dummies(df[new_column], prefix=column, dummy_na=True).astype('int')
        else:
            dummies = pd.get_dummies(df[new_column], prefix=column).astype('int')
        # Drop the column with the lowest variance
        dummy_variances = dummies.var()
        lowest_var_column = dummy_variances.idxmin()
        dummies = dummies.drop(columns=lowest_var_column)
        return dummies,'categorical'

def phenotype_required_columns(dataframe,col_dict:dict,date_fmt:str,study_desgin:str, single_sex:bool, 
                               sex_value_dict:dict):
    """
    
    This function processes required columns in the given phenotype dataframe. 

    Parameters
    ----------
    dataframe : pd.DataFrame
        A pandas DataFrame recording the phenoptype data.
    col_dict : dict
        A dictionary that maps the descriptions to their corresponding column names in the dataframe for the required columns.
    date_fmt : str
        The format string for parsing dates in the date columns.
    study_desgin : str
        The study design of the analysis.
    single_sex : bool
        Indicating whether the study population is restricted to females or males.
    sex_value_dict : dict
        A dictionary that maps coding to specific sex.

    Returns
    -------
    None
        The function does not return a value but assign new attributes, raises exceptions and prints warnings if issues with the data integrity are detected. 

    Raises
    ------
    ValueError if:
        - there are missing values in the required columns; 
        - there are formatting issues with the date columns, 
        - duplicate entries in the Participant ID, 
        - the exposure variable does not consist of exactly two unique binary values as expected.
        - the sex variable does not consist of exactly two unique binary values as expected.

    """
    for col in col_dict.keys():
        if dataframe[col_dict[col]].isna().any():
            raise ValueError(f'Warning: The {col} column contains missing value')
    
    eid_col = col_dict['Participant ID']
    index_date_col = col_dict['Index date']
    end_date_col = col_dict['End date']
    sex_col = col_dict['Sex']
    #reverse the dictionary
    sex_value_dict_r = {v:k for k,v in sex_value_dict.items()}
    
    try:
        dataframe[index_date_col] = dataframe[index_date_col].apply(lambda x: datetime.strptime(x,date_fmt))
    except:
        raise ValueError("Check the date format for column Index date")
    try:
        dataframe[end_date_col] = dataframe[end_date_col].apply(lambda x: datetime.strptime(x,date_fmt))
    except:
        raise ValueError("Check the date format for column End date")
    if dataframe[eid_col].duplicated().any():
        raise ValueError("Duplicates found in Participant ID column, which is not allowed")

    # process the exposure column
    if study_desgin != "exposed-only cohort":
        exposure_col = col_dict['Exposure']
        unique_vals = dataframe[exposure_col].unique()
        n_unique_vals = len(unique_vals)
        if n_unique_vals != 2:
            raise ValueError('The exposure variable does not have 2 unique values')
        else:
            if all(isinstance(x, (int, np.integer, float, np.floating)) and x in {0, 1} for x in unique_vals):
                None
            else:
                raise TypeError("The exposure variable does not have 2 unique values")
    
    #prcess the sex column
    unique_vals = dataframe[sex_col].unique()
    n_unique_vals = len(unique_vals)
    if n_unique_vals != 2 and single_sex == False:
        raise ValueError('The sex variable does not have 2 unique values')
    elif n_unique_vals != 1 and single_sex == True:
        raise ValueError('single_sex is True but the sex variable does not have 1 unique value')
    else:
        if all(isinstance(x, (int, np.integer, float, np.floating)) and x in {0, 1} for x in unique_vals):
            if single_sex == True:
                sex_contained = [sex_value_dict_r[x] for x in unique_vals]
                print(f"Warning: only {sex_contained} individuals are presented in the data")
            else:
                None
        else:
            raise TypeError("The 'Sex' variable must be coded as 1 (female) and 0 (male).")
    
    
def medical_records_process(
    medical_records:str,
    col_dict:dict,
    code_type:str,
    date_fmt:str,
    chunk_n,
    seperator,
    exclusion_list:list,
    all_phecode_dict:dict,
    phecode_map:dict
):
    """
    Read the medical records dataframe (in chunks), mapped to phecode and update the provided nested dictionary.

    Parameters
    ----------
    medical_records : str
        The file path containing medical records data in CSV or TSV format.
    
    col_dict : dict
        A dictionary that maps the descriptions to their corresponding column names in the dataframe for the required columns.
    
    code_type : str
        Type of the diagnosis code, either ICD-9 or ICD-10 (CM or WHO).
    
    date_fmt : str
        The format string for parsing dates in the date columns.
    
    chunk_n : int
        Chunk size for read the medical records data.
    
    seperator : str
        Seperator used for reading.

    exclusion_list : list
        A list of diagnosis codes to exclude.
    
    all_phecode_dict : dict
        A empty nested dictionary for updating. Keys are participant ID, values are empty dictionary
    
    phecode_map : dict
        Phecode map dictionary. Keys are ICD codes, and values are list of mapped phecodes.

    Returns : tuple
    -------
        Update the "all_phecode_dict" with the provided medical records data.
        Also return tuple of some statistics.

    """
    eid_col = col_dict['Participant ID']
    icd_col = col_dict['Diagnosis code']
    date_col = col_dict['Date of diagnosis']
    
    n_total_missing = 0
    n_total_read = 0
    n_total_records = 0
    n_total_trunc_4 = 0
    n_total_trunc_3 = 0
    n_total_no_trunc = 0
    n_total_no_mapping = 0
    no_mapping_list = {}
    
    chunks = pd.read_csv(medical_records,sep=seperator,iterator=True,chunksize=chunk_n,
                         usecols=[eid_col,icd_col,date_col])
    for chunk in chunks:
        len_before = len(chunk)
        #convert the icd_col to string
        chunk[icd_col] = chunk[icd_col].astype(str)
        #filtering the participant ID
        chunk = chunk[chunk[eid_col].isin(all_phecode_dict)]
        #drop records in the exclusion list
        if exclusion_list:
            chunk = chunk[~chunk[icd_col].str[:5].isin(exclusion_list) & 
                          ~chunk[icd_col].str[:4].isin(exclusion_list) & 
                          ~chunk[icd_col].str[:3].isin(exclusion_list)]
        len_valid = len(chunk)
        #drop na values
        chunk.dropna(how='any', inplace=True)
        n_missing = len_valid - len(chunk)
        #comvert the date column to datetime
        chunk[date_col] = chunk[date_col].apply(lambda x: datetime.strptime(x,date_fmt))
        if 'ICD-9' in code_type: 
            chunk[icd_col] = chunk[icd_col].apply(lambda x: decimal_to_short(x))
        elif 'ICD-10' in code_type:
            chunk[icd_col] = chunk[icd_col].apply(lambda x: x.replace('.',''))
        else:
            raise ValueError(f'unrecognized diagnosis code type {code_type}')
        n_total_read += len_before
        n_total_missing += n_missing
        n_total_records += len_valid
        print(f'{n_total_read:,} records read, {n_total_records:,} left after filltering on participant ID/exclusion list of diagnosis codes, {n_total_missing:,} records with missing values excluded.')
        #drop records not in the list
        #sort and drop duplicates
        chunk = chunk.sort_values(by=[date_col],ascending=True).drop_duplicates()

        #mapping
        new_phecode_lst = []
        for patient_id, icd, date in chunk[[eid_col,icd_col,date_col]].values:
            if icd in phecode_map:
                for phecode in phecode_map[icd]:
                    new_phecode_lst.append([patient_id,phecode,date])
                n_total_no_trunc += 1
            #try trunction 
            elif icd[0:4] in phecode_map:
                for phecode in phecode_map[icd[0:4]]:
                    new_phecode_lst.append([patient_id,phecode,date])
                n_total_trunc_4 += 1
            #try trunction 
            elif icd[0:3] in phecode_map:
                for phecode in phecode_map[icd[0:3]]:
                    new_phecode_lst.append([patient_id,phecode,date])
                n_total_trunc_3 += 1
            else:
                n_total_no_mapping += 1
                try:
                    no_mapping_list[icd] += 1
                except:
                    no_mapping_list[icd] = 1
                continue
        del chunk #save memory
        #update the all_phecode_dict with phecode
        for patient_id, phecode, new_date in new_phecode_lst:
            if phecode in all_phecode_dict[patient_id]:
                old_date,n = all_phecode_dict[patient_id][phecode]
                all_phecode_dict[patient_id][phecode] = [min(old_date,new_date),n+1]
            else:
                all_phecode_dict[patient_id][phecode] = [new_date,1]
    #print final report
    print(f'Total: {n_total_records:,} diagnosis records processed, {n_total_missing:,} records with missing values were excluded.')
    print(f'{n_total_no_trunc:,} diagnosis records mapped to phecode without truncating.')
    print(f'{n_total_trunc_4:,} diagnosis records mapped to phecode after truncating to 4 digits.')
    print(f'{n_total_trunc_3:,} diagnosis records mapped to phecode after truncating to 3 digits.')
    print(f'{n_total_no_mapping:,} diagnosis records not mapped to any phecode.')
    return n_total_records,n_total_missing,n_total_trunc_4,n_total_trunc_3,n_total_no_mapping,no_mapping_list
    

def diagnosis_history_update(diagnosis_dict:dict, n_diagnosis_dict:dict, history_dict:dict, start_date_dict:dict, end_date_dict:dict, phecode_dict:dict):
    """
    Update the diagnosis and history dictionary (nested) using provided phecode dictionary.

    Parameters
    ----------
    diagnosis_dict : dict
        DiseaseNetworkData.diagnosis, a nested dictionary with participant ID being the key and a dictionary being the value, where phecode is the key and date of diagnosis is the value.
    
    n_diagnosis_dict : dict
        DiseaseNetworkData.n_diagnosis, a nested dictionary with participant ID being the key and a dictionary being the value, where phecode is the key and number of occurence is the value.

    history_dict : dict
        DiseaseNetworkData.history, a nested dictionary with participant ID being the key and a list of phecode being the value.
    
    start_date_dict : dict
        A dictionary recording the date of follow-up start for each participant.
        
    end_date_dict : dict
        A dictionary recording the date of follow-up end for each participant.
    
    phecode_dict : dict
        A phecode dictionary from func medical_records_process().

    Returns : int
    -------
        Update the diagnosis_dict and history_dict inplace using phecode_dict.
        Return the numbr of phecodes that are invalid (> date of follow-up end) for each participant.

    """
    n_invalid = {}
    for patient_id in phecode_dict:
        for phecode,[date,n] in phecode_dict[patient_id].items():
            #first update the number of phecode occurence
            n_diagnosis_dict[patient_id][phecode] = n_diagnosis_dict[patient_id].get(phecode,0) + n
            #then update the diangosis and history dictionary
            if date > end_date_dict[patient_id]:
                try:
                    n_invalid[patient_id] += 1
                except:
                    n_invalid[patient_id] = 1
            elif date <= start_date_dict[patient_id]:
                if phecode in diagnosis_dict[patient_id]:
                    del diagnosis_dict[patient_id][phecode]
                    history_dict[patient_id].append(phecode)
                elif phecode not in history_dict[patient_id]:
                    history_dict[patient_id].append(phecode)
            elif date > start_date_dict[patient_id]:
                if phecode in history_dict[patient_id]:
                    continue
                elif phecode not in diagnosis_dict[patient_id] or date < diagnosis_dict[patient_id][phecode]:
                    diagnosis_dict[patient_id][phecode] = date
    return n_invalid
    
def threshold_check(proportion_threshold, n_threshold, n_exposed):
    """
    Checks and determines the threshold for analysis based on either a proportion or an absolute count.

    Parameters:
        proportion_threshold (float or None): Proportion of exposed individuals to use as a threshold.
        n_threshold (int or None): Absolute number of exposed individuals to use as a threshold.
        n_exposed (int): Total number of exposed individuals.

    Returns:
        int: The calculated or checked original threshold.

    Raises:
        ValueError: If both 'proportion_threshold' and 'n_threshold' are specified, or if any input
                    is invalid (e.g., types or ranges).
    """
    if proportion_threshold and n_threshold:
        raise ValueError("'n_threshold' and 'proportion_threshold' cannot be specified at the same time.")
    if proportion_threshold is not None:
        if not isinstance(proportion_threshold, float):
            raise TypeError("The 'proportion_threshold' must be a float.")
        if not (0 < proportion_threshold <= 1):
            raise ValueError("'proportion_threshold' must be between 0 and 1.")
        return int(n_exposed * proportion_threshold)
    elif n_threshold is not None:
        if not isinstance(n_threshold, int):
            raise TypeError("The 'n_threshold' must be an int.")
        if not (0 <= n_threshold <= n_exposed):
            raise ValueError("'n_threshold' must be a non-negative integer less than or equal to the number of exposed individuals.")
        return n_threshold
    else:
        raise ValueError("Either 'n_threshold' or 'proportion_threshold' must be specified.")

def n_process_check(n_process:int,analysis_name:str):
    """
    Check the number of process specified for analysis.
    Also check the operation system information to decide the maximum number of process that can be used and the method of start the process.

    Parameters:
        n_process (int): The number of process to use.

    Returns:
        None

    Side Effects:
        Prints a message about the process usage for the analysis.

    Raises:
        ValueError: If 'n_process' is not a positive integer.
    """
    if not isinstance(n_process, int):
        raise TypeError("The 'n_process' must be an int.")
    if n_process == 1:
        print(f'Multiprocessing is not used for {analysis_name} analysis.')
        return n_process,None
    elif n_process > 1:
        import os
        #close multi-threading
        os.environ["MKL_NUM_THREADS"] = '1'
        os.environ["OPENBLAS_NUM_THREADS"] = '1'
        os.environ["OMP_NUM_THREADS"] = '1'
        os.environ["THREADPOOL_LIMIT"] = '1'
        os.environ["VECLIB_MAXIMUM_THREADS"] = '1'
        #add missing environment variables
        os.environ["BLIS_NUM_THREADS"] = '1'
        os.environ["NUMEXPR_NUM_THREADS"] = '1'
        os.environ["MKL_DYNAMIC"] = 'FALSE'  # Disable dynamic threading for MKL
        #system information
        max_process = os.cpu_count()
        operation_system = os.name
        if n_process > max_process:
            raise ValueError(f"The specified number of process is greater than the number of logical cores ({max_process}).")
        if operation_system == 'nt':
            start_method = 'spawn'
        elif operation_system == 'posix':
            start_method = 'fork'
        else:
            raise ValueError(f"Unsupported operation system: {operation_system}")
        print(f'Use {n_process} process and set start method to {start_method} for {analysis_name} analysis.')
        return n_process,start_method
    else:
        raise ValueError("The specified number of process is not valid. Please enter a positive integer.")

def correction_method_check(correction, cutoff):
    """
    Check the p-value correction method and its cutoff threshold.

    Parameters:
        correction (str): The p-value correction method to use.
        cutoff (float): The cutoff threshold for significance.

    Returns:
        None

    Raises:
        ValueError: If 'correction' is not a recognized method, or if 'cutoff' is invalid.
    """
    methods_lst = ['bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg',
                   'hommel', 'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky', 'none']
    if not isinstance(correction, str):
        raise TypeError("The 'correction' must be a string.")
    if correction not in methods_lst:
        raise ValueError(f"Choose from the following p-value correction methods: {methods_lst}")
    if not isinstance(cutoff, float):
        raise TypeError("The 'cutoff' must be a float.")
    if not (0 < cutoff < 1):
        raise ValueError("'cutoff' must be between 0 and 1, exclusive.")

def log_file_detect(file_path,prefix):
    """
    Try to get the log file path and test it.

    Parameters
    ----------
    file_path : str
        Input log file path or None.
    
    prefix : str
        Prefix to be added except for the 12 randomly generated characters.
    
    Returns : str and a message
    -------
        The final log file path and message.

    """
    import tempfile
    import os
    
    if not file_path:
        characters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
        random_file_name = f'DiseaseNet_{prefix}_'+''.join(np.random.choice(list(characters),12))+'.log'
        temp_folder_path = tempfile.gettempdir()
        temp_file = os.path.join(temp_folder_path,random_file_name)
        return temp_file,f'Logging to {temp_file}'
    else:
        if not file_path.endswith('.log'):
            file_path += '.log'
        try:
            with open(file_path,'wb') as f:
                f.write(''.encode())
        except:
            characters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
            random_file_name = f'DiseaseNet_{prefix}_'+''.join(np.random.choice(list(characters),12))+'.log'
            temp_folder_path = tempfile.gettempdir()
            temp_file = os.path.join(temp_folder_path,random_file_name)
            return (temp_file, f'{file_path} does not exist or is not writable, logging to {temp_file}.')
    return file_path,f'Logging to {file_path}'


def filter_phecodes(phecode_info, system_inc=None, system_exl=None, phecode_inc=None, phecode_exl=None):
    """
    Filters a list of phecodes based on inclusion and exclusion criteria for systems and phecodes.

    Parameters:
        phecode_info (dict): Dictionary where keys are phecodes and values are dictionaries containing system information.
        system_inc (list, optional): List of systems to include.
        system_exl (list, optional): List of systems to exclude.
        phecode_inc (list, optional): List of phecodes to include.
        phecode_exl (list, optional): List of phecodes to exclude.

    Returns:
        list: Filtered list of phecodes.

    Raises:
        ValueError: If inclusion and exclusion parameters conflict or contain invalid values.
    """
    # Initialize the full list of phecodes
    phecode_lst_all = list(phecode_info.keys())
    system_all = set([phecode_info[x]['category'] for x in phecode_lst_all])

    # check input criteria
    if system_inc and system_exl:
        raise ValueError("'system_inc' and 'system_exl' cannot both be specified.")
    if phecode_inc and phecode_exl:
        raise ValueError("'phecode_inc' and 'phecode_exl' cannot both be specified.")
    if (system_inc or system_exl) and (phecode_inc or phecode_exl):
        print('Warning: both phecode and system level filters applied may result in redundant or ambiguous outcomes.')

    # Filter based on system inclusion
    if system_inc:
        if not isinstance(system_inc, list):
            raise TypeError("The 'system_inc' must be a list.")
        if len(system_inc) == 0:
            raise ValueError("The 'system_inc' list is empty.")
        system_unidentified = [x for x in system_inc if x not in system_all]
        if system_unidentified:
            raise ValueError(f"The following phecode systems from 'system_inc' are invalid: {system_unidentified}")
        phecode_lst_all = [x for x in phecode_lst_all if phecode_info[x]['category'] in system_inc]

    # Filter based on system exclusion
    if system_exl:
        if not isinstance(system_exl, list):
            raise TypeError("The 'system_exl' must be a list.")
        if len(system_exl) == 0:
            raise ValueError("The 'system_exl' list is empty.")
        system_unidentified = [x for x in system_exl if x not in system_all]
        if system_unidentified:
            raise ValueError(f"The following phecode systems from 'system_exl' are invalid: {system_unidentified}")
        phecode_lst_all = [x for x in phecode_lst_all if phecode_info[x]['category'] not in system_exl]

    # Filter based on phecode inclusion
    if phecode_inc:
        if not isinstance(phecode_inc, list):
            raise TypeError("The 'phecode_inc' must be a list.")
        if len(phecode_inc) == 0:
            raise ValueError("The 'phecode_inc' list is empty.")
        phecode_unidentified = [x for x in phecode_inc if x not in phecode_lst_all]
        if phecode_unidentified:
            raise ValueError(f"The following phecodes from 'phecode_inc' are invalid: {phecode_unidentified}")
        phecode_lst_all = [x for x in phecode_lst_all if x in phecode_inc]

    # Filter based on phecode exclusion
    if phecode_exl:
        if not isinstance(phecode_exl, list):
            raise TypeError("The 'phecode_exl' must be a list.")
        if len(phecode_exl) == 0:
            raise ValueError("The 'phecode_exl' list is empty.")
        phecode_unidentified = [x for x in phecode_exl if x not in phecode_lst_all]
        if phecode_unidentified:
            raise ValueError(f"The following phecodes from 'phecode_exl' are invalid: {phecode_unidentified}")
        phecode_lst_all = [x for x in phecode_lst_all if x not in phecode_exl]

    # Check if any phecodes remain
    if not phecode_lst_all:
        raise ValueError("No phecodes remain after applying filtering at the phecode and system levels.")

    return phecode_lst_all

def check_history_exclusion(exl_lst,history_list,n_diagnosis_dict, threshold):
    """
    Given the nested exclusion list for a phecode, and the history diagnosis list 
    as well as number of phecode occurence dictionary for a individual, return whether this individual is eligible (0) or not (1)

    Parameters:
        exl_lst: nested list of phecode for exlusion
        history_list: dictionary (in the future version) or list of history phecode diagnosis for a individual
        n_diagnosis_dict: dictionary recording the number of phecode occurence for a individual
        threshold: min_required_icd_codes
    
    Returns:
        eligible (0) or not (1)
    """
    for lst in exl_lst:
        if len(set(history_list).intersection(set(lst)))>0 and sum([n_diagnosis_dict.get(d,0) for d in lst])>=threshold:
            return 1
        else:
            continue
    return 0

def time_first_diagnosis(d_lst,diagnosis_dict,n_diagnosis_dict,threshold):
    """
    Given the list of leaf phecodes for a root phecode, and the diagnosis dictionary and number of phecode occurence dictionary for a individual, 
    return the first date (or nan if no records or not eligible) of diagnosis for that root phecode

    Parameters:
        d_lst: list of leaf phecode for a root phecode
        diagnosis_dict: diagnosis dictionary for a individual
        n_diagnosis_dict: dictionary recording the number of phecode occurence for a individual
        threshold: min_required_icd_codes
    
    Returns:
        first diagnosis date or pd.NaT
    """
    if sum([n_diagnosis_dict.get(d, 0) for d in d_lst]) >= threshold:
        dates = [diagnosis_dict.get(x, pd.NaT) for x in d_lst]
        valid_dates = [d for d in dates if pd.notna(d)]
        return min(valid_dates) if valid_dates else pd.NaT
    else:
        return pd.NaT


def states_p_adjust(df,p_col,correction,cutoff,prefix_sig_col,prefix_padj_col):
    """
    Applies p-value adjustment for multiple comparisons and determines significance based on a cutoff.
    
    Parameters:
        df (pd.DataFrame): The input DataFrame containing p-values.
        p_col (str): Column name in 'df' containing p-values to adjust.
        correction (str): The method to use for p-value adjustment.
        cutoff (float): The significance cutoff value for the adjusted p-values.
        prefix_sig_col (str): Prefix for the new column indicating significance.
        prefix_padj_col (str): Prefix for the new column containing adjusted p-values.
    
    Returns:
        pd.DataFrame: The input DataFrame with added columns for adjusted p-values and significance.
    
    Raises:
        ValueError: If the specified column does not exist or other parameter issues arise.
    """
    from statsmodels.stats.multitest import multipletests
    
    df_na = df[df[p_col].isna()]
    df_nona = df[~df[p_col].isna()]
    reject_, corrected_p, _, _ = multipletests(df_nona[p_col],method=correction,alpha=cutoff)
    df_nona[f'{prefix_sig_col}_significance'] = reject_
    df_nona[f'{prefix_padj_col}_adjusted'] = corrected_p
    df_na[f'{prefix_sig_col}_significance'] = False
    df_na[f'{prefix_padj_col}_adjusted'] = np.nan
    result = pd.concat([df_nona,df_na])
    return result

def write_log(log_file, message, retries=50, delay=0.001):
    import time
    """
    Writes a log message to a file with a retry mechanism for simplicity.
    
    Parameters
    ----------
        log_file (str): Path to the log file.
        message (str): Log message to write to the file.
        retries (int): Number of retries in case of a file access conflict.
        delay (float): Delay in seconds between retries.
    """
    for attempt in range(retries):
        try:
            with open(log_file, 'ab') as f:
                f.write(message.encode())
            return
        except PermissionError:
            if attempt < retries - 1:
                time.sleep(delay)
            else:
                raise PermissionError(f"Failed to write to {log_file} after {retries} retries.")

def get_d_lst(lower,upper):
    """
    Get all the possible phecodes between two provided phecodes
    
    Parameters
    ----------
        lower : float the smaller phecode
        upper : float the larger phecode

    Returns
    -------
        A list of phecodes.

    """
    if lower>=upper:
        raise ValueError("The larger phecode is larger or equal to lower phecode.")
    
    n_step = int((upper - lower) / 0.01 + 1)
    d_lst = np.linspace(lower, upper, n_step)

    return d_lst

def get_exclison_lst(exl_range_str):
    """
    Get a list of phecodes for exclision give the exclusion description string.

    Parameters
    ----------
    exl_range_str : string
        String of phecodes exclusion criteria.

    Returns
    -------
    A set of phecodes

    """
    exl_list = []
    if pd.isna(exl_range_str):
        return exl_list
    else:
        for range_ in exl_range_str.split(','):
            exl_lower,exl_higher = float(range_.split('-')[0]), float(range_.split('-')[1])
            exl_list += list(get_d_lst(exl_lower,exl_higher))
        exl_list = set(exl_list)
    
    return set(exl_list)

def d1d2_from_diagnosis_history(df:pd.DataFrame, id_col:str, sex_col:str, sex_value_dict:dict, 
                                phecode_lst:list, disease_pair_index:dict, history_dict:dict, diagnosis_dict:dict, n_diagnosis_dict:dict,
                                phecode_info_dict:dict, min_interval_days:int, max_interval_days:int, min_icd_num:int) -> dict:
    """
    Construct d1->d2 disease pairs for each individual from a list of significant phecodes.
    
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

    Returns
    -------
        D1->D2 dictionary.
    
    """
    ineligible_disease_dict = {}
    eligible_withdate_dict = {}
    d1d2_temporl_pair_dict = {}
    d1d2_com_pair_dict = {}
    all_diagnosis_level = {} #this is the newly generated list correspond to the specified phecode level (consider the occurence as well)
    
    from itertools import combinations

    #generate a dictionary mapping leaf phecode to root phecode 
    node_dict = phecode_leaf_to_root(phecode_info_dict)
    phecode_set = set(phecode_info_dict.keys())
    for id_,sex in tqdm(df[[id_col,sex_col]].values, mininterval=15,smoothing=0):
        temp_dineligible_list = []
        temp_dpair_temporal_lst = []
        temp_dpair_com_lst = []
        temp_deligible_dict_withdate = {}
        diagnosis_ = diagnosis_dict[id_]
        n_diagnosis_ = n_diagnosis_dict[id_]
        history_ = history_dict[id_]
        #generate a new all diagnosed phecode list correspond to the specified level, consider the n. of occurence as well
        all_diag_ind = set([node_dict.get(x,x) for x in history_+list(diagnosis_.keys())]).intersection(phecode_set)
        all_diag_ind = [x for x in all_diag_ind if sum([n_diagnosis_.get(leaf,0) for leaf in phecode_info_dict[x]['leaf_list']])>=min_icd_num]
        all_diagnosis_level[id_] = all_diag_ind
        #generate eligible disease dictionary
        for phecode in phecode_lst:
            leaf_lst = phecode_info_dict[phecode]['leaf_list']
            exl_lst = phecode_info_dict[phecode]['exclude_list']
            sex_specific = phecode_info_dict[phecode]['sex']
            if (check_history_exclusion(exl_lst,history_,n_diagnosis_,min_icd_num)==0) and (sex_specific=='Both' or sex_value_dict[sex_specific]==sex):
                date = time_first_diagnosis(leaf_lst,diagnosis_,n_diagnosis_,min_icd_num)
                if not pd.isna(date):
                    temp_deligible_dict_withdate[phecode] = date
            else:
                temp_dineligible_list.append(phecode)
        #generate disease pair dictionary
        if len(temp_deligible_dict_withdate) <= 1:
            ineligible_disease_dict[id_] = temp_dineligible_list
            d1d2_temporl_pair_dict[id_] = temp_dpair_temporal_lst
            d1d2_com_pair_dict[id_] = temp_dpair_com_lst
            eligible_withdate_dict[id_] = temp_deligible_dict_withdate
        else:
            for d1,d2 in combinations(temp_deligible_dict_withdate,2):
                date1, date2 = temp_deligible_dict_withdate[d1], temp_deligible_dict_withdate[d2]
                if abs((date1 - date2).days) > max_interval_days:
                    temp_dpair_com_lst.append(disease_pair_index[f'{d1}_{d2}']) #order insensitive
                    continue
                elif abs((date1 - date2).days) <= min_interval_days:
                    temp_dpair_com_lst.append(disease_pair_index[f'{d1}_{d2}']) #order insensitive
                    continue
                else:
                    if date1 > date2:
                        temp_dpair_temporal_lst.append(disease_pair_index[f'{d2}_{d1}'])
                    else:
                        temp_dpair_temporal_lst.append(disease_pair_index[f'{d1}_{d2}'])
            #save for the individual
            ineligible_disease_dict[id_] = temp_dineligible_list
            d1d2_temporl_pair_dict[id_] = temp_dpair_temporal_lst
            d1d2_com_pair_dict[id_] = temp_dpair_com_lst
            eligible_withdate_dict[id_] = temp_deligible_dict_withdate
    
    #final dictionary
    trajectory_dict = {'disease_pair_index' : disease_pair_index,
                       'ineligible_disease':ineligible_disease_dict,
                       'eligible_disease_withdate':eligible_withdate_dict,
                       'd1d2_temporal_pair':d1d2_temporl_pair_dict,
                       'd1d2_com_pair':d1d2_com_pair_dict,
                       'all_diagnosis_level':all_diagnosis_level}
    
    return trajectory_dict

def validate_method_specific_kwargs(method:str, kwargs:dict):
    """
    Validate method-specific kwargs for comorbidity_network/disease_trajectory function
    
    Parameters
    ----------
    method : str comorbidity_network/disease_trajectory analysis method
    kwargs : dict method-specific kwargs
    
    Returns
    -------
    dict : parameter_dict containing validated method-specific parameters
    """
    alpha = None
    auto_penalty = None
    scaling_factor = None
    n_PC = None
    explained_variance = None
    
    # Method-specific kwargs validation
    if method == 'RPCN':
        # RPCN-specific parameters
        alpha = kwargs.pop('alpha', None)
        auto_penalty = kwargs.pop('auto_penalty', True)
        alpha_range = kwargs.pop('alpha_range',(1,15))
        scaling_factor = kwargs.pop('scaling_factor', 1)

        if not isinstance(auto_penalty, bool):
            raise TypeError(f"'auto_penalty' should be a bool, got {type(auto_penalty).__name__}.")

        if auto_penalty:
            # If auto_penalty is True, alpha should not be provided
            if 'alpha' in kwargs:
                raise ValueError("When 'auto_penalty' is True, 'alpha' should not be provided.")
            if not isinstance(alpha_range,tuple):
                raise TypeError(f"'alpha_range' should be a tuple, got {type(alpha_range).__name__}.")
            if len(alpha_range)>2 or alpha_range[1] <= alpha_range[0]:
                raise ValueError("The range defined in 'alpha_range' is invalid.")
            for alpha_value in alpha_range:
                if not isinstance(alpha_value, int) or alpha_value<0:
                    raise TypeError(f"Upper and lower bounds defined in 'alpha_range' should be int>=0, got {alpha_value}.")
            if not isinstance(scaling_factor, (int, float)) or scaling_factor <= 0:
                raise TypeError(f"'scaling_factor' should be a positive scalar, got {type(scaling_factor).__name__}.")
            parameter_dict = {'method':'RPCN','auto_penalty':True,'alpha':alpha, 'alpha_range':alpha_range, 'scaling_factor':scaling_factor}
        else:
            # If auto_penalty is False, alpha must be provided, while alpha_range shoud not be provided
            if 'alpha_range' in kwargs:
                raise ValueError("When 'auto_penalty' is False, 'alpha_range' should not be provided.")
            if 'scaling_factor' in kwargs:
                raise ValueError("When 'auto_penalty' is False, 'scaling_factor' should not be provided.")
            if alpha is None:
                raise ValueError("When 'auto_penalty' is False, 'alpha' must be provided.")
            if not isinstance(alpha, (int, float)):
                raise TypeError(f"'alpha' should be a scalar, got {type(alpha).__name__}.")
            if alpha<0:
                raise ValueError("'alpha' must be a positive scaler.")
            elif alpha<=1:
                print("Warning: The provided 'alpha' value is low. The optimal 'alpha' value is normally within the range of (1, 15].")
            parameter_dict = {'method':'RPCN','auto_penalty':False,'alpha':alpha, 'alpha_range':alpha_range, 'scaling_factor':scaling_factor}
            
    elif method == 'PCN_PCA':
        # PCN_PCA-specific parameters
        n_PC = kwargs.pop('n_PC', 5)
        explained_variance = kwargs.pop('explained_variance', None)

        if explained_variance is not None:
            if not isinstance(explained_variance, float):
                raise TypeError(f"'explained_variance' should be a float, got {type(explained_variance).__name__}.")
            if not (0 < explained_variance <= 1):
                raise ValueError("'explained_variance' should be between 0 and 1.")
            if 'n_PC' in kwargs:
                raise ValueError("Cannot specify both 'n_PC' and 'explained_variance'. Please choose one.")
            parameter_dict = {'method':'PCN_PCA','explained_variance':explained_variance}
        else:
            if not isinstance(n_PC, int):
                raise TypeError(f"'n_PC' should be an integer, got {type(n_PC).__name__}.")
            if n_PC <= 0:
                raise ValueError("'n_PC' should be a positive integer.")
            parameter_dict = {'method':'PCN_PCA','n_PC':n_PC}

    elif method == 'CN':
        # CN does not have method-specific parameters, ensure none are provided
        method_specific_params = {'alpha', 'L1_wt', 'auto_penalty', 'n_PC', 'explained_variance'}
        if any(param in kwargs for param in method_specific_params):
            raise ValueError(f"No additional parameters are required for method '{method}'.")
        parameter_dict = {'method':'CN'}
    
    return parameter_dict

def check_kwargs_com_tra(method:str,comorbidity_strength_cols:list,binomial_test_result_cols:list,**kwargs):
    """
    Check the **kwargs from comorbidity_network/disease_trajectory function
    
    Parameters
    ----------
    method : str comorbidity_network/disease_trajectory analysis method
    comorbidity_strength_cols : list columns of comorbidity_strength_result dataframe
    binomial_test_result_cols : list columns of binomial_test_result dataframe
    **kwargs : **kwargs from comorbidity_network/disease_trajectory function

    Returns
    -------
    List of required parameters

    """
    # Default keyword arguments
    default_kwargs = {
        'phecode_d1_col': 'phecode_d1',
        'phecode_d2_col': 'phecode_d2',
        'significance_phi_col': 'phi_p_significance',
        'significance_RR_col': 'RR_p_significance',
        'significance_binomial_col': 'binomial_p_significance'
    }
    # Update default_kwargs with user-provided kwargs
    column_kwargs = {k: kwargs.pop(k, v) for k, v in default_kwargs.items()}
    # check that no unexpected keyword arguments are present for column definitions
    allowed_column_kwargs = set(default_kwargs.keys())
    extra_column_kwargs = set(kwargs.keys()) - set([
        'alpha', 'auto_penalty','alpha_range', 'n_PC', 'explained_variance' ,'enforce_time_interval', 'scaling_factor'])
    invalid_column_kwargs = extra_column_kwargs - allowed_column_kwargs
    if invalid_column_kwargs:
        raise ValueError(f"Invalid keyword arguments: {invalid_column_kwargs}")
    
    # check that required columns exist in the DataFrames
    required_columns_strength = [
        column_kwargs['phecode_d1_col'],
        column_kwargs['phecode_d2_col'],
        column_kwargs['significance_phi_col'],
        column_kwargs['significance_RR_col']
    ]
    required_columns_binomial = [
        column_kwargs['phecode_d1_col'],
        column_kwargs['phecode_d2_col'],
        column_kwargs['significance_binomial_col']
    ]
    missing_columns_strength = [col for col in required_columns_strength if col not in comorbidity_strength_cols]
    if binomial_test_result_cols is None:
        missing_columns_binomial = []
    else:
        missing_columns_binomial = [col for col in required_columns_binomial if col not in binomial_test_result_cols]
    if missing_columns_strength:
        raise ValueError(f"The following required columns are missing in 'comorbidity_strength_result': {missing_columns_strength}")
    if missing_columns_binomial:
        raise ValueError(f"The following required columns are missing in 'binomial_test_result': {missing_columns_binomial}")
        
    # Validate method-specific kwargs
    parameter_dict = validate_method_specific_kwargs(method, kwargs)
    
    #enforce parameter
    enforce_time_interval = kwargs.pop('enforce_time_interval', True)
    if not isinstance(enforce_time_interval, bool):
        raise TypeError(f"'enforce_time_interval' should be a bool, got {type(enforce_time_interval).__name__}.")
    parameter_dict.update({'enforce_time_interval':enforce_time_interval})
    
    # After method-specific validation, ensure no unexpected kwargs remain
    if kwargs:
        raise ValueError(f"Unexpected keyword arguments provided: {kwargs.keys()}")
    

    return [parameter_dict,column_kwargs['phecode_d1_col'],column_kwargs['phecode_d2_col'],
            column_kwargs['significance_phi_col'],column_kwargs['significance_RR_col'],column_kwargs['significance_binomial_col']]

def matching_var_check(matching_var_dict:dict,phenotype_info:dict):
    """
    Check the specified matching_var dictionary.

    Parameters
    ----------
    matching_var_dict : dict, matching variable and matching criteria.
    phenotype_info : dict, phenotype information from DiseaseNetwork data.

    Returns
    -------
    None.

    """
    if not isinstance(matching_var_dict,dict):
        raise TypeError(f"'matching_var_dict' should be a dictionary, got {type(matching_var_dict).__name__}.")
    
    if len(matching_var_dict) == 0:
        raise ValueError("Empty 'matching_var_dict'. At least one matching variable is needed.")
    
    
    for var,val in matching_var_dict.items():
        if var not in phenotype_info['phenotype_covariates_type']:
            raise ValueError(f'Unknown matching variable: {var}')
        var_type = phenotype_info['phenotype_covariates_type'][var]
        if var_type == 'categorical' and val != 'exact':
            raise ValueError(f"Matching variable {var} is categorical, the matching criteria should always be 'exact', got {val}.")
        if var_type == 'continuous':
            if not isinstance(val, (int,float)) or val < 0:
                raise ValueError(f"Invalid matching criteria {val} for matching variable {var}.")
                

def covariates_check(covariates:list,phenotype_info:dict,matching_var_dict:dict=None,
                     exclude:bool=False):
    """
    Check whether the given list of covariates is valid and return the transformed covariates name.
    If matching dictionary is given, also check whether there is any overlap with the matching variables.
    
    Parameters
    ----------
    covariates : list or None, list of covariates.
    matching_var_dict : dict or None, dictionary of matching variables and criteria used
    phenotype_info : dict, phenotype information from DiseaseNetwork data
    exclude : bool, whether to exclude the covariates from the final list or raise an error.

    Returns
    -------
    List of final covariates.

    """
    sex_col = phenotype_info['phenotype_col_dict']['Sex']
    all_possible_covars = [sex_col] + phenotype_info['phenotype_covariates_original']
    if covariates is not None:
        if not isinstance(covariates, list):
            raise TypeError(f"'covariates' should be a list of strings, got {type(covariates).__name__}.")
        invalid_vars = [x for x in covariates if x not in all_possible_covars]
        if invalid_vars:
            raise ValueError(f"Invalid covariates '{invalid_vars}'. Allowed covariates are: {all_possible_covars}.")
    else:
        if matching_var_dict is not None and sex_col in matching_var_dict:
            all_possible_covars.remove(sex_col)
            covariates = all_possible_covars 
        else:
            covariates = all_possible_covars
    #transform the covariates
    var_rename_dict = phenotype_info['phenotype_covariates_converted']
    covariates_final = [x for var in covariates for x in var_rename_dict.get(var,[var])]
    
    if matching_var_dict is None:
        return covariates_final
    else:
        for var in matching_var_dict:
            if var in covariates and phenotype_info['phenotype_covariates_type'][var]=='categorical':
                if exclude == True:
                    covariates_final.remove(var)
                else:
                    raise ValueError(f'Categorical covariate {var} has already been used for matching.')
        return covariates_final

def find_best_alpha_and_vars(model, best_range, alpha_lst, co_vars):
    """
    Function to find the best alpha for L1 regularization and the corresponding non-zero variables using an early stopping rule based on consecutive AIC increases.
    
    Parameters:
        model (statsmodels object): The statistical model to be fitted.
        best_range (tuple): A tuple (min_alpha, max_alpha) defining the range to explore.
        alpha_lst (float): The alpha multiplier applied during regularization.
        co_vars (list): List of variable names in the model.
    
    Returns:
        tuple: (final_best_alpha, final_disease_vars) where 'final_best_alpha' is the alpha value that
               results in the lowest AIC before AIC starts to increase consistently, and 'final_disease_vars'
               is a list of variables that are non-zero at this alpha level.
    """
    refined_alphas = np.linspace(best_range[0], best_range[1], num=best_range[1]-best_range[0]+1)
    refined_aic_dict = {}
    refined_vars_dict = {}
    min_aic = float('inf')
    counter = 0  # Counter to track the number of increases after a minimum
    counter_failed = 0 #counter for failed models
    thresold = thresold_failed = int(len(refined_alphas)*0.5) # early stop threshold and early stop threshold for failed models

    for alpha in refined_alphas:
        try:
            result = model.fit_regularized(method='l1', alpha=alpha_lst*alpha, disp=False)
            non_zero_indices = np.nonzero(result.params != 0)[0]
            refined_vars_dict[alpha] = [co_vars[i] for i in non_zero_indices] #constant is included
            refined_aic_dict[alpha] = result.aic
        except:
            # If the model fails to converge, set AIC to infinity
            refined_aic_dict[alpha] = float('inf')
            refined_vars_dict[alpha] = []
            counter_failed += 1
            continue

        if counter_failed >= thresold_failed: #stop if failed 2 times
            break

        # Check for AIC minimum and count increases
        if refined_aic_dict[alpha] < min_aic:
            min_aic = refined_aic_dict[alpha]
            counter = 0  # Reset counter on new minimum
        else:
            counter += 1  # Increment counter on increase
        
        # Break loop if AIC increases 5 times consecutively after a minimum
        if counter >= thresold:
            break

    final_best_alpha = min(refined_aic_dict, key=refined_aic_dict.get)
    if refined_aic_dict[final_best_alpha] == float('inf'):
        raise ValueError(f"Models failed to fit when trying to find the best alpha for L1 regularization, consider change the 'alpha_range' or 'scaling_factor'.")
    else:
        final_best_alpha = final_best_alpha * alpha_lst[-1]
        final_disease_vars = refined_vars_dict[final_best_alpha]
    return final_best_alpha, final_disease_vars

#decprecated function
def check_variance_vif(df:pd.DataFrame, 
                       covar_lst:list, 
                       disease_var_lst:list=None, 
                       pca_var_lst:list=None, 
                       group_col:str=None,
                       vif_cutoff:int=200) -> list:
    """
    Check within group variance and Variance inflation factor (VIF) for the all the covariates.
    Covariate with 0 within group variance will be removed.
    Covariate with VIF > 5 will be removed.

    Args:
        df (pd.DataFrame): dataframe for analysis
        covar_lst (list): phenotypical covariates, should always be included
        disease_var_lst (list): disease variables
        pca_var_lst (list): pca variables
        group_col (str): group column for within group variance check
        vif_cutoff (int): VIF cutoff value

    Returns:
        list: a nested list of covariates to be removed.
    """
    from statsmodels.stats.outliers_influence import variance_inflation_factor

    #keep the order, make sure covar_lst is the last for VIF check
    all_vars = disease_var_lst if disease_var_lst is not None else [] + pca_var_lst if pca_var_lst is not None else [] + covar_lst
    var_removed = {}
    #return empty list if no variables are provided
    if len(all_vars) == 0:
        return []

    #if group_col is provided, calculate the within group variance
    if group_col is not None:
        for var in all_vars:
            group_means = df.groupby(group_col)[var].transform('mean')
            ss_within = ((df[var] - group_means) ** 2).sum()
            variance_within = ss_within / (len(df) - df[group_col].nunique())
            if variance_within == 0: var_removed[var] = 'within_group_variance==0'
    #if not provided, check the whole dataset variance
    else:
        for var in all_vars:
            if df[var].var() == 0: var_removed[var] = 'variance==0'

    #always check VIF
    all_vars = [x for x in all_vars if x not in var_removed] #remove the variables with 0 within group variance
    vars_all_set = all_vars.copy() #variables that remained unchanged during the VIF check loop
    for var in vars_all_set:
        index_ = all_vars.index(var)
        vif = variance_inflation_factor(df[all_vars],index_)
        if vif >= vif_cutoff:
            all_vars.remove(var)
            var_removed[var] = f'VIF={vif}'
   
    #check the removed variables
    covar_removed = {x:var_removed[x] for x in var_removed if x in covar_lst}
    if disease_var_lst is not None:
        disease_removed = {x:var_removed[x] for x in var_removed if x in disease_var_lst}
        return covar_removed,disease_removed
    elif pca_var_lst is not None:
        pca_removed = {x:var_removed[x] for x in var_removed if x in pca_var_lst}
        return covar_removed,pca_removed
    elif disease_var_lst is None and pca_var_lst is None:
        return covar_removed
    else:
        raise ValueError("Invalid input.")

def compute_vif_sm_exact(df):
    """
    Compute Variance Inflation Factor (VIF) for each feature in the DataFrame.
    This function handles perfectly collinear variables by setting their VIF to infinity.
    """
    result = {}
    X = df.dropna().astype(float).values
    # detect perfectly collinear columns by checking correlation matrix
    corr_matrix = df.corr().abs()
    
    # find groups of perfectly collinear variables
    collinear_groups = []
    processed_columns = set()

    for i, col in enumerate(df.columns):
        if col in processed_columns:
            continue
        # find columns perfectly correlated with this one
        perfect_matches = [c for c in df.columns if c != col and corr_matrix.loc[col, c] >= 0.999]
        if perfect_matches:
            # Create a group including the current column
            group = [col] + perfect_matches
            collinear_groups.append(group)
            processed_columns.update(group)
        else:
            processed_columns.add(col)
    
    #for each group of collinear variables, mark others as infinite VIF
    perfect_collinear = []
    for group in collinear_groups:
        perfect_collinear.extend(group[1:])
    for col in perfect_collinear:
        result[col] = np.inf

    # Calculate VIF for non-collinear variables
    remaining_cols = [col for col in df.columns if col not in perfect_collinear]
    # Add the representative columns from collinear groups back to remaining columns
    for group in collinear_groups:
        if group[0] not in remaining_cols:
            remaining_cols.append(group[0])
    
    if remaining_cols:
        X_remaining = df[remaining_cols].values
        G = X_remaining.T @ X_remaining
        G_inv = np.linalg.inv(G)
        # Calculate VIF for remaining variables
        remaining_vifs = dict(zip(remaining_cols, np.diag(G) * np.diag(G_inv)))
        result.update(remaining_vifs)
    
    return result

def check_variance_vif_single(df:pd.DataFrame, forcedin_var_lst:list,
                              covar_lst:list, group_col:str=None,
                              vif_cutoff:int=None) -> list:
    """
    Check within group variance and Variance inflation factor (VIF) for the all the covariates.
    Covariate with 0 within group variance will be removed.
    Covariate with VIF > 5 will be removed.

    Args:
        df (pd.DataFrame): dataframe for analysis
        forcedin_var_lst (list): covariates should always be included, such as constant and exposure variable
        covar_lst (list): phenotypic covariates
        vif_cutoff (int): VIF cutoff value (int) or some prespecified settings

    Returns:
        list: a nested list of covariates to be removed.
    """
    from statsmodels.stats.outliers_influence import variance_inflation_factor
    #prespecified vif value for some types
    vif_dict = {'phenotypic_covar':500,'disease_covar':20,
                'pca_covar':500}
    vif_cutoff = vif_dict.get(vif_cutoff,vif_cutoff)

    #keep the order, make sure covar_lst is the last for VIF check
    var_removed = {}
    #return empty list if no variables are provided
    if len(covar_lst) == 0:
        return []

    #if group_col is provided, calculate the within group variance
    if group_col is not None:
        for var in covar_lst:
            group_means = df.groupby(group_col)[var].transform('mean')
            ss_within = ((df[var] - group_means) ** 2).sum()
            variance_within = ss_within / (len(df) - df[group_col].nunique())
            if variance_within == 0: var_removed[var] = 'within_group_variance==0'
    #if not provided, check the whole dataset variance
    else:
        for var in covar_lst:
            if df[var].var() == 0: var_removed[var] = 'variance==0'

    #always check VIF
    covar_lst = [x for x in covar_lst if x not in var_removed] #remove the variables with 0 within group variance
    vars_all = covar_lst + forcedin_var_lst #remaining covariates plus forced-in variables
    #vars_all_set = vars_all.copy() #variables that remained unchanged during the VIF check loop

    while True:
        all_vifs = compute_vif_sm_exact(df[vars_all])
        offenders = {v:val for v, val in all_vifs.items()
                     if v not in forcedin_var_lst and val > vif_cutoff}
        if not offenders:
            break
        # Remove the variable with the highest VIF
        worst = max(offenders, key=offenders.get)
        var_removed[worst] = f"VIF={offenders[worst]:.2f}"
        vars_all.remove(worst)

    #return the dictionary of variables to be excluded
    return var_removed
