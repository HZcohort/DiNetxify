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

    """
    for col in col_dict.keys():
        if dataframe[col_dict[col]].isna().any():
            raise ValueError(f'Warning: The {col} column contains missing value')
    
    eid_col = col_dict['Participant ID']
    index_date_col = col_dict['Index date']
    end_date_col = col_dict['End date']
    exposure_col = col_dict['Exposure']
    
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
        raise ValueError('The exposure variable does not have 2 unique vales')
    else:
        if all(isinstance(x, (int, np.integer, float, np.floating)) and x in {0, 1} for x in unique_vals):
            None
        else:
            raise ValueError("The exposure variable does not have 2 unique vales")
    
    
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
    
    chunk_n : TYPE
        Chunk size for read the medical records data.
    
    seperator : TYPE
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
    
    
    
    
    
    
    
    
    
    
    
    
    








