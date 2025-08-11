import pandas as pd
import numpy as np
import scipy.stats as stats
import warnings
warnings.filterwarnings('ignore')

def categorical(df,var,group_var,group_var_value,table1):
    """
    Function to generate the table rows for categorical variables
    """
    var_values = df[var].value_counts().index
    all_rows = []
    #a header row for the variable
    all_rows.append([f'{var} (n, %)']+['']*len(group_var_value)+[''])
    #a row for each value of the variable
    for value in var_values:
        row = [f'{var}={value}']
        for group in group_var_value:
            group_data = df[df[group_var]==group]
            numbers = group_data[df[var]==value].shape[0]
            percentage = numbers/group_data.shape[0]*100
            row.append(f'{numbers:,} ({percentage:.2f}%)')
        if value == var_values[-1]:
            if len(group_var_value) >= 2:
                cont_table = pd.crosstab(df[var],df[group_var])
                try:
                    chitest = stats.chi2_contingency(cont_table)
                    p_value = f"Chi-squared test p-value={chitest[1]:.3e}"
                except:
                    chitest = stats.fisher_exact(cont_table)
                    p_value = f"Fisher's exact test p-value={chitest[1]:.3e}"
                row.append(p_value)
            else:
                row.append('NA')
        else:
            row.append('')
        all_rows.append(row)
    table1.extend(all_rows)

def continuous_normal(df,var,group_var,group_var_value,table1):
    """
    Function to generate the table rows for normally distributted continuous variables
    """
    row_continous = [f'{var} (mean, SD)']
    for group in group_var_value:
        group_data = df[df[group_var]==group]
        #using f-string to format the number
        mean_sd = f'{group_data[var].mean():,.2f} ({group_data[var].std():,.2f})'
        row_continous.append(mean_sd)
    if len(group_var_value) == 2:
        ttest = stats.ttest_ind(df[df[group_var]==group_var_value[0]][var],df[df[group_var]==group_var_value[1]][var])
        p_value = f"Student's T-test p-value={ttest[1]:.3e}"
        row_continous.append(p_value)
    else: 
        row_continous.append('NA')
    table1.append(row_continous)

def continuous_unnormal(df,var,group_var,group_var_value,table1):
    """
    Function to generate the table rows for unnormally distributted continuous variables
    """
    row_continous = [f'{var} (median, IQR)']
    for group in group_var_value:
        group_data = df[df[group_var]==group]
        #using f-string to format the number
        median_iqr = f'{group_data[var].median():,.2f} ({group_data[var].quantile(0.25):,.2f}-{group_data[var].quantile(0.75):,.2f})'
        row_continous.append(median_iqr)
    if len(group_var_value) == 2:
        utest = stats.mannwhitneyu(df[df[group_var]==group_var_value[0]][var],df[df[group_var]==group_var_value[1]][var])
        p_value = f'Mann-Whitney U test p-value={utest[1]:.3e}'
        row_continous.append(p_value)
    else:
        row_continous.append('NA')
    table1.append(row_continous)

def desceibe_table(df:pd.DataFrame,var_lst:list,var_type_lst:list,group_var:str,
                   sex_value_dict:dict, continuous_stat_mode:str='auto',group_var_value:list=[1,0]):
    """
    Function to generate the simple table1

    df: pandas DataFrame
    var_lst: list of variable names
    var_type_lst: list of variable types
    group_var: the group variable
    sex_value_dict: mapping dictionary for sex variable
    continuous_stat_mode: either 'auto', 'normal' or 'nonnormal'
    group_var_value: the value of group variable
    """
    table1 = []
    #table 1 header
    header = ['Variable']
    for group in group_var_value:
        n_total = df[df[group_var]==group].shape[0]
        header.append(f"{group_var}={group} (n={n_total:,})")
    header.append('Test and p-value')
    #table1.append(header)
    #table 1 content
    for var,type_ in zip(var_lst,var_type_lst):
        if type_ == 'continuous':
            if continuous_stat_mode=='auto':
                p = stats.shapiro(df[var].values)[1]
                if p>0.05:
                    continuous_normal(df,var,group_var,group_var_value,table1)
                else:
                    continuous_unnormal(df,var,group_var,group_var_value,table1)
            elif continuous_stat_mode=='normal':
                continuous_normal(df,var,group_var,group_var_value,table1)
            elif continuous_stat_mode=='nonnormal':
                continuous_unnormal(df,var,group_var,group_var_value,table1)
        elif type_ == 'categorical' or type_ == 'binary':
            if var == 'sex':
                df[var] = df[var].map({j:i for i,j in sex_value_dict.items()})
            #fill the na with 'NA'
            df[var] = df[var].fillna('NAN')
            categorical(df,var,group_var,group_var_value,table1)
    return pd.DataFrame(table1,columns=header)

