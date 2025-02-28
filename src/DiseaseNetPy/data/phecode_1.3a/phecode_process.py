# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 16:32:53 2024

@author: Can Hou, Elizabeth Unnur Gisladottir
"""

"""
This is a test version of the phecode mapping and is not an official release. 
It is not recommended for general research use. 
The primary update is the addition of mappings for several ICD-9 and ICD-10 codes that appear to be specific to the Swedish inpatient/outpatient registers. 
These new mappings are appended to the end of the corresponding CSV/EXCEL files.
Two new .npy files were generated using these updated mapping files.
"""

# %%phecode data prepare
import pandas as pd
import numpy as np

def level_number(n):
    if int(n) == n:
        return 1
    elif int(round(n*10,3)) == round(n*10,3):
        return 2
    elif int(round(n*100,3)) == round(n*100,3):
        return 3
    else:
        return 9

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
    
    n_step = int(round(upper - lower,3) / 0.01 + 1)
    d_lst = np.linspace(lower, upper, n_step)
    d_lst = [round(x,3) for x in d_lst]

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

def group_by_integer(codes):
    """
    Groups a list of numbers into sublists based on their integer part.
    
    Parameters:
        codes (list of float): The list of numbers to group.
        
    Returns:
        list of list of float: A nested list where each sublist contains numbers sharing the same integer part.
    """
    from collections import defaultdict

    grouped = defaultdict(list)
    for code in codes:
        key = int(code)  # Get the integer part of the number
        grouped[key].append(code)
    
    # Convert dictionary values to a list of lists and return
    return list(grouped.values())

phecode_definition = pd.read_csv(r'src/DiseaseNetPy/data/phecode_1.2/phecode_info.csv')
phecode_definition['level'] = phecode_definition['phecode'].apply(lambda x: level_number(x))
phecode_definition.fillna({'sex':'Both','category':'others'},inplace=True)
all_phecode = set(phecode_definition['phecode'].to_list())

#dict for level 1 list
phecode_level1_lst = set([int(x) for x in phecode_definition['phecode'].values])
phecode_definition_level1 = phecode_definition[phecode_definition['phecode'].isin(phecode_level1_lst)]

phecode_dict = {}
for i in phecode_definition_level1.index:
    temp_dict = {}
    phecode = phecode_definition_level1.loc[i,'phecode']
    for col in ['phenotype','phecode_exclude_range','sex','category','level']:
        temp_dict[col] = phecode_definition_level1.loc[i,col]
    #list of level code
    leaf_lst = [x for x in all_phecode if int(x) == phecode]
    temp_dict['leaf_list'] = set(leaf_lst)
    #list of exclude phecode
    exl_list = get_exclison_lst(temp_dict['phecode_exclude_range'])
    exl_list = [x for x in exl_list if x in all_phecode]
    #group exl_list
    exl_list = set(list(exl_list))
    exl_list = group_by_integer(exl_list)
    temp_dict['exclude_list'] = exl_list
    phecode_dict[phecode] = temp_dict
    
np.save(r'src/DiseaseNetPy/data/phecode_1.2/level1_info.npy',phecode_dict)


#dict for level 2 list
phecode_level2_lst = phecode_definition[phecode_definition['level']==2]['phecode'].to_list()
phecode_level2_lst_1s = [x for x in phecode_level1_lst if x not in [int(x) for x in phecode_level2_lst]]
phecode_level2_lst = phecode_level2_lst + phecode_level2_lst_1s
phecode_definition_level2 = phecode_definition[phecode_definition['phecode'].isin(phecode_level2_lst)]

phecode_dict = {}
for i in phecode_definition_level2.index:
    temp_dict = {}
    phecode = phecode_definition_level2.loc[i,'phecode']
    for col in ['phenotype','phecode_exclude_range','sex','category','level']:
        temp_dict[col] = phecode_definition_level2.loc[i,col]
    #list of level code
    if temp_dict['level'] == 1:
        leaf_lst = [x for x in all_phecode if int(x) == phecode]
        temp_dict['leaf_list'] = set(leaf_lst)
    elif temp_dict['level'] == 2:
        leaf_lst = get_d_lst(phecode,round(phecode+0.09, 2))
        leaf_lst = [x for x in leaf_lst if x in all_phecode]
    temp_dict['leaf_list'] = set(leaf_lst)
    #list of exclude phecode
    exl_list = get_exclison_lst(temp_dict['phecode_exclude_range'])
    exl_list = [x for x in exl_list if x in all_phecode]
    #group exl_list
    exl_list = set(list(exl_list))
    exl_list = group_by_integer(exl_list)
    temp_dict['exclude_list'] = exl_list
    phecode_dict[phecode] = temp_dict

np.save(r'src/DiseaseNetPy/data/phecode_1.2/level2_info.npy',phecode_dict)

# %% mapping files
def decimal_to_short(code):
    """
    Convert an ICD9 code from decimal format to short format.
    """
    parts = code.split(".")
    parts[0] = parts[0].zfill(3)
    return "".join(parts)

#WHO mapping
who_9 = pd.read_excel(r'src\DiseaseNetPy\data\phecode_1.3a\phecode_map_who_icd9.xlsx')
who_9['ICD'] = who_9['icd9'].apply(lambda x: decimal_to_short(str(x)))
#who_9[who_9['ICD'].apply(lambda x: len(x)<=2)]

who_9_dict = {}
for icd,phe in who_9[['ICD','phecode']].values:
    try:
        who_9_dict[icd].append(phe)
    except:
        who_9_dict[icd] = [phe]
np.save(r'src\DiseaseNetPy\data\phecode_1.3a\ICD-9-WHO.npy',who_9_dict)

who_10 = pd.read_excel(r'src\DiseaseNetPy\data\phecode_1.3a\phecode_map_who_icd10.xlsx')
who_10['ICD'] = who_10['ICD10'].apply(lambda x: x.replace('.',''))
#who_10[who_10['ICD'].apply(lambda x: len(x)<=2)]

who_10_dict = {}
for icd,phe in who_10[['ICD','PHECODE']].values:
    try:
        who_10_dict[icd].append(phe)
    except:
        who_10_dict[icd] = [phe]
np.save(r'src\DiseaseNetPy\data\phecode_1.3a\ICD-10-WHO.npy',who_10_dict)
