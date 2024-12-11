# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 10:07:40 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

import pandas as pd
import numpy as np
from datetime import datetime

def decimal_to_short(code:float) -> str:
    """

    Convert an ICD9 code from decimal format to short format.

    Parameters:
    ----------

    """
    parts = str(code).split(".")
    parts[0] = parts[0].zfill(3)
    return "".join(parts)

def read_check_csv(path_file:str, cols_check:list, date_cols:list, date_fmt:str, separator_lst:list=[',', '\t'],
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
            if df[cols_check].isnull().any().any():
                print("Warning: missing values found in the first 50 rows")
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


def diff_date_years(dataframe,date1:str,date2:str,rounding:int=2):
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
        Difference between the two dates in year (2 digits).

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
            return df[[new_column]].map(mapping),'binary'
    elif (df[new_column].dtype == float or df[new_column].dtype == int or df[new_column].dtype == object) and n_unique_vals>=10:
        # Treat as continuous
        n_missing_before = df[new_column].isna().sum()
        df[new_column] = pd.to_numeric(df[new_column], errors='coerce').round(2)
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
            return pd.get_dummies(df[new_column], prefix=column, dummy_na=True, drop_first=True),'categorical'
        else:
            return pd.get_dummies(df[new_column], prefix=column, drop_first=True),'categorical'

def phenotype_required_columns(dataframe,col_dict:dict,date_fmt:str):
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
    exposure_col = col_dict['Exposure']
    sex_col = col_dict['Sex']
    
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
    if n_unique_vals != 2:
        raise ValueError('The sex variable does not have 2 unique values')
    else:
        if all(isinstance(x, (int, np.integer, float, np.floating)) and x in {0, 1} for x in unique_vals):
            None
        else:
            raise TypeError("The sex variable does not have 2 unique values")
    
    
def medical_records_process(medical_records:str,col_dict:dict,code_type:str,date_fmt:str,chunk_n,seperator,
                            all_phecode_dict:dict,phecode_map:dict):
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
    n_total_no_mapping = 0
    no_mapping_list = {}
    
    chunks = pd.read_csv(medical_records,sep=seperator,iterator=True,chunksize=chunk_n,
                         usecols=[eid_col,icd_col,date_col])
    for chunk in chunks:
        len_before = len(chunk)
        chunk = chunk[chunk[eid_col].isin(all_phecode_dict)]
        len_valid = len(chunk)
        chunk.dropna(how='any',inplace=True)
        n_missing = len_valid - len(chunk)
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
        print(f'{n_total_read:,} records read ({n_total_records:,} included after filltering on participant ID), {n_total_missing:,} records with missing values excluded.')
        #drop records not in the list
        #sort and drop duplicates
        chunk = chunk.sort_values(by=[date_col],ascending=True).drop_duplicates(subset=[eid_col,icd_col],keep='first')
        #mapping
        new_phecode_lst = []
        for patient_id, icd, date in chunk[[eid_col,icd_col,date_col]].values:
            if icd in phecode_map:
                for phecode in phecode_map[icd]:
                    new_phecode_lst.append([patient_id,phecode,date])
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
                if new_date < all_phecode_dict[patient_id][phecode]:
                    all_phecode_dict[patient_id][phecode] = new_date
            else:
                all_phecode_dict[patient_id][phecode] = new_date
    #print final report
    print(f'Total: {n_total_records:,} diagnosis records processed, {n_total_missing:,} records with missing values were excluded.')
    print(f'{n_total_trunc_4:,} diagnosis records mapped to phecode after truncating to 4 digits.')
    print(f'{n_total_trunc_3:,} diagnosis records mapped to phecode after truncating to 3 digits.')
    print(f'{n_total_no_mapping:,} diagnosis records not mapped to any phecode.')
    return n_total_records,n_total_missing,n_total_trunc_4,n_total_trunc_3,n_total_no_mapping,no_mapping_list
    

def diagnosis_history_update(diagnosis_dict:dict, history_dict:dict, start_date_dict:dict, end_date_dict:dict, phecode_dict:dict):
    """
    Update the diagnosis and history dictionary (nested) using provided phecode dictionary.

    Parameters
    ----------
    diagnosis_dict : dict
        DiseaseNetworkData.diagnosis, a nested dictionary with participant ID being the key and a dictionary being the value, where phecode is the key and date of diagnosis is the value.
        
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
        for phecode,date in phecode_dict[patient_id].items():
            if date <= start_date_dict[patient_id]:
                if phecode not in history_dict[patient_id]:
                    history_dict[patient_id].append(phecode)
            elif date > start_date_dict[patient_id] and date <= end_date_dict[patient_id]:
                if phecode not in diagnosis_dict[patient_id] or date < diagnosis_dict[patient_id][phecode]:
                    diagnosis_dict[patient_id][phecode] = date
            else:
                try:
                    n_invalid[patient_id] += 1
                except:
                    n_invalid[patient_id] = 1
    return n_invalid   
    
def validate_threshold(proportion_threshold, n_threshold, n_exposed):
    """
    Validates and determines the threshold for analysis based on either a proportion or an absolute count.

    Parameters:
        proportion_threshold (float or None): Proportion of exposed individuals to use as a threshold.
        n_threshold (int or None): Absolute number of exposed individuals to use as a threshold.
        n_exposed (int): Total number of exposed individuals.

    Returns:
        int: The calculated or validated threshold.

    Raises:
        ValueError: If both `proportion_threshold` and `n_threshold` are specified, or if any input
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

def validate_n_cpus(n_cpus,analysis_name):
    """
    Validates the number of CPUs specified for analysis.

    Parameters:
        n_cpus (int): The number of CPUs to use.

    Returns:
        None

    Side Effects:
        Prints a message about the CPU usage for the analysis.

    Raises:
        ValueError: If `n_cpus` is not a positive integer.
    """
    if not isinstance(n_cpus, int):
        raise TypeError("The 'n_cpus' must be an int.")
    if n_cpus == 1:
        print('Multi-threading is not used.')
    elif n_cpus > 1:
        print(f'Use {n_cpus} CPU cores for {analysis_name} analysis.')
    else:
        raise ValueError("The specified number of CPUs is not valid. Please enter a positive integer.")

def validate_correction_method(correction, cutoff):
    """
    Validates the p-value correction method and its cutoff threshold.

    Parameters:
        correction (str): The p-value correction method to use.
        cutoff (float): The cutoff threshold for significance.

    Returns:
        None

    Raises:
        ValueError: If `correction` is not a recognized method, or if `cutoff` is invalid.
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

def log_file_detect(file_path):
    """
    Try to get the log file path and test it.

    Parameters
    ----------
    file_path : str
        Input log file path or None. 
    
    Returns : str and a message
    -------
        The final log file path and message.

    """
    import tempfile
    import os
    
    if not file_path:
        characters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
        random_file_name = 'DiseaseNet_'+''.join(np.random.choice(list(characters),12))+'.log'
        temp_folder_path = tempfile.gettempdir()
        temp_file = os.path.join(temp_folder_path,random_file_name)
        return temp_file,f'log written to {temp_file}'
    else:
        if not file_path.endswith('.log'):
            file_path += '.log'
        try:
            with open(file_path,'wb') as f:
                f.write(''.encode())
        except:
            characters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
            random_file_name = 'DiseaseNet_'+''.join(np.random.choice(list(characters),12))+'.log'
            temp_folder_path = tempfile.gettempdir()
            temp_file = os.path.join(temp_folder_path,random_file_name)
            return (temp_file, f'{file_path} does not exist or is not writable, log written to {temp_file}.')
    return file_path,f'log written to {file_path}'


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

    # Validate input criteria
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

def states_p_adjust(df,p_col,correction,cutoff,prefix_sig_col,prefix_padj_col):
    """
    Applies p-value adjustment for multiple comparisons and determines significance based on a cutoff.
    
    Parameters:
        df (pd.DataFrame): The input DataFrame containing p-values.
        p_col (str): Column name in `df` containing p-values to adjust.
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
    df_na[f'{prefix_padj_col}_adjusted'] = np.NaN
    result = pd.concat([df_nona,df_na])
    return result

def write_log(log_file, message, retries=50, delay=0.1):
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

def d1d2_from_diagnosis_history(df:pd.DataFrame, id_col:str, sex_col:str, phecode_lst:list, history_dict:dict, diagnosis_dict:dict,
                                phecode_info_dict:dict, time_interval_days:int) -> dict:
    """
    Construct d1->d2 disease pairs for each individual from a list of significant phecodes.
    
    Parameters
    ----------
        df : dataframe of phenotype data, contains at id and sex columns
        id_col : id column in the df
        sex_col : sex column in the df
        phecode_lst : list of significant phecodes
        history_dict : dictionary containing medical records history
        diagnosis_dict : dictionary containing diagnosis and date
        phecode_info_dict : phecode information
        time_interval_days : time interval required for d1-d2 disease pair construction

    Returns
    -------
        D1->D2 dictionary.
    
    """
    sex_value_dict = {'Female':1,'Male':0}
    trajectory_dict = {}
    from itertools import combinations
    
    for id_,sex in df[[id_col,sex_col]].values:
        temp_deligible_dict = {}
        temp_dpair_lst = []
        diagnosis_ = diagnosis_dict[id_]
        history_ = history_dict[id_]
        #generate eligible disease dictionary
        for phecode in phecode_lst:
            leaf_lst = phecode_info_dict[phecode]['leaf_list']
            exl_lst = phecode_info_dict[phecode]['exclude_list']
            sex_specific = phecode_info_dict[phecode]['sex']
            if len(exl_lst.intersection(set(history_)))==0 and (sex_specific=='Both' or sex_value_dict[sex_specific]==sex):
                try:
                    date = min([diagnosis_[x] for x in leaf_lst if x in diagnosis_])
                except:
                    date = pd.NaT
            temp_deligible_dict[phecode] = date
        #generate disease pair dictionary
        temp_deligible_dict_withdate = {i:j for i,j in temp_deligible_dict.items() if not pd.isna(j)}
        if len(temp_deligible_dict_withdate) <= 1:
            trajectory_dict[id_] = {'eligible_disease':temp_deligible_dict,
                                    'd1d2_pair':temp_dpair_lst}
        else:
            for d1,d2 in combinations(temp_deligible_dict_withdate,2):
                date1, date2 = temp_deligible_dict_withdate[d1], temp_deligible_dict_withdate[d2]
                if abs((date1 - date2).days) <= time_interval_days:
                    continue
                if date1 > date2:
                    temp_dpair_lst.append((d1,d2))
                else:
                    temp_dpair_lst.append((d2,d1))
        #save for the individual
        trajectory_dict[id_] = {'eligible_disease':temp_deligible_dict,
                                'd1d2_pair':temp_dpair_lst}
        
    return trajectory_dict






























