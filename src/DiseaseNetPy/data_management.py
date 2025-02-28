# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 23:48:05 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

import pandas as pd
import numpy as np
import os
from datetime import datetime
from .utility import convert_column, phenotype_required_columns, read_check_csv
from .utility import medical_records_process, diagnosis_history_update, d1d2_from_diagnosis_history

import warnings
warnings.filterwarnings('ignore')

class DiseaseNetworkData:
    """
    A class for handling disease network data creation and operations, for use in DiseaseNet module.

    Parameters:
    ----------
    study_design : str
        Specify the type of study design, either "cohort", "matched cohort", or "registry".
    
    phecode_level : int
        The level of phecode to use for analysis, where level 1 (with a total of 585 medical conditions) corresponds to 3-digit ICD-10 codes and level 2 (a total of 1257 medical conditions) to 4-digit ICD-10 codes. 
        Level 2 phecodes offer a more granular analysis with potentially smaller sample sizes per disease category. 
        For larger studies, level 2 phecodes may enhance result interpretation.
        For smaller studies, level 1 is recommended to maintain statistical power.
    
    min_required_icd_codes : int, default=1
        The minimum number of ICD codes mapping to a specific phecode required for the phecode to be considered valid.
        For example, if set to 2, a single diagnosis record will not be sufficient to count as an occurrence.
        Ensure that your medical records are complete (i.e., not limited to only the first occurrence for each code) when using this parameter.
    
    date_fmt : str, default='%Y-%m-%d'
        The format of the date fields in your phenotype and medical records data.
    
    phecode_version : str, default='1.2'
        The version of the phecode system used for converting diagnosis codes. 
        Version 1.2 is the offical version of the phecode system, with mapping files available for ICD-9-CM, ICD-9-WHO, ICD-10-CM, and ICD-10-WHO codes.
        While option 1.3a is provided, its an unofficial version, and not recommended for general use.
    """
    
    def __init__(
        self, 
        study_design:str='cohort', 
        phecode_level:int=1, 
        min_required_icd_codes:int=1,
        date_fmt:str='%Y-%m-%d', 
        phecode_version:str='1.2',
    ):
        #fixed attributes
        #phenotype data
        self.__module_dir = os.path.dirname(__file__)
        self.__study_design_options = ['matched cohort','cohort','registry']
        self.__id_col = 'eid'
        self.__exposure_col = 'exposure'
        self.__sex_col = 'sex'
        self.__index_date_col = 'index_date'
        self.__end_date_col = 'end_date'
        self.__mathcing_identifier_col = 'group'
        self.__phenotype_col_dict = {
            'matched cohort':{'Participant ID':self.__id_col,
                                'Exposure':self.__exposure_col,
                                'Sex':self.__sex_col,
                                'Index date': self.__index_date_col,
                                'End date': self.__end_date_col,
                                'Matching identifier':self.__mathcing_identifier_col},
            'cohort':{'Participant ID':self.__id_col,
                        'Exposure':self.__exposure_col,
                        'Sex':self.__sex_col,
                        'Index date': self.__index_date_col,
                        'End date': self.__end_date_col},
            'registry':{"Participant ID":self.__id_col,
                        "Sex":self.__sex_col,
                        "Index date":self.__index_date_col,
                        "End date":self.__end_date_col}
        }
        #medical records data
        self.__diagnosis_code_options = ['ICD-9-CM', 'ICD-9-WHO', 'ICD-10-CM', 'ICD-10-WHO']
        self.__phecode_level_options = [1, 2]
        self.__phecode_version_options = ['1.2','1.3a']
        self.__medical_records_cols = ['Participant ID','Diagnosis code','Date of diagnosis']
        self.__medical_recods_info = {}
        #default column value
        self.__sex_value_dict = {'Female':1,'Male':0}
        
        #check study design
        if not study_design in self.__study_design_options:
            raise ValueError(f"Choose from the following study design: {self.__study_design_options}")
        #check date format
        try:
            datetime.strftime(datetime.now(), date_fmt)
        except:
            raise ValueError("The specified date format is invalid. Visit https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior for details.")
        #check phecode_level
        if phecode_level not in self.__phecode_level_options:
            raise ValueError(f"Choose from the following phecode level {self.__phecode_level_options}")
        #check phecode version
        if phecode_version not in self.__phecode_version_options:
            raise ValueError(f"Supported phecode version {self.__phecode_version_options}")
        #check minimal require ICD code number
        if not isinstance(min_required_icd_codes,int) or min_required_icd_codes<=0:
            raise TypeError("The input 'min_required_icd_codes' must be a int number > 0.")
        
        #assign parameter attributes
        self.study_design = study_design
        self.date_fmt = date_fmt
        self.phecode_level = phecode_level
        self.phecode_version = phecode_version
        self.min_required_icd_codes = min_required_icd_codes
        #load necessary phecode files
        self.__phecode_info_npyfile = os.path.join(self.__module_dir,f'data/phecode_{self.phecode_version}/level{self.phecode_level}_info.npy')
        self.phecode_info = np.load(self.__phecode_info_npyfile,allow_pickle=True).item()
        #reset key attributes
        self.phenotype_df = None
        self.diagnosis = None
        self.n_diagnosis = None #new dictionary, recording the number of occurence of each phecode
        self.history = None
        self.trajectory = None
        self.__warning_phenotype = []
        self.__warning_medical_records = []
        #variable names placeholder
        self.__varialbe_name_place_holder = list(self.__phenotype_col_dict[self.study_design].values()) #required variables
        self.__varialbe_name_place_holder += ['flag_exl','outcome_date','outcome','time_years'] #variables reserved for cox analysis
        self.__varialbe_name_place_holder += ['d1','d2','constant','d1_date','d2_date','outcome_date','group_matching_ids'] #variables reserved for conditional/unconditional logistic regression
        self.__varialbe_name_place_holder += [f'PCA_{i}' for i in range(100)] #variables reserved for PCA
        self.__varialbe_name_place_holder += [str(x) for x in self.phecode_info] #variables reserved for other diseases
        self.__varialbe_name_place_holder += ['follow_up',self.__exposure_col] #other reserved variables
    
    def phenotype_data(
        self, 
        phenotype_data_path:str, 
        column_names:dict, 
        covariates:list, 
        force:bool=False
    ) -> None:
        """
        
        Merges phenotype and medical records data into the main data attribute.
        
        Parameters:
        ----------
        phenotype_data_path : str
            The file path containing phenotype data in CSV or TSV format. 
            The file must include a header row with column names.
                
        column_names : dict
            A dictionary mapping required variable names to their corresponding identifiers in the dataset. 
            Expected keys include 'Participant ID', 'Index date', 'End date', 'Exposure' (for cohort and matched cohortstudy), 'Sex', and 'Matching identifier' (for matched cohort study). 
            For example:
            column_names={'Participant ID': 'eid',
                          'Exposure': 'status',
                          'Sex': 'sex',
                          'Index date': 'index_date',
                          'End date': 'final_date',
                          'Matching identifier': 'group'}
            The 'Exposure' variable must be coded as 0 (unexposed) and 1 (exposed). 
            The 'Sex' variable must be coded as 1 (female) and 0 (male).
            Dates must be formatted as '%Y-%m-%d' unless specified otherwise.
            Records with missing values in any required columns are not allowed.
        
        covariates : list
            A list of variable names representing additional covariates, such as ['age', 'BMI']. 
            Provide an empty list if no additional covariates are included. 
            The function will automatically detect and convert variable types. 
            Individuals with missing values in continuous variables will be removed, while those missing in categorical variables will be categorized separately.
        
        force : bool, default=False
            If True, the data will be loaded and existing attributes will be overwritten, even if they contain data. 
            The default is False, which will raise an error if data already exists.

        Returns:
        ----------
        None
            This function modifies the object's main data attribute in-place.
        
        """
        if not force:
            data_attrs = ['phenotype_df', 'diagnosis', 'n_diagnosis', 'history', 'trajectory']
            for attr in data_attrs:
                if getattr(self, attr) is not None:
                    raise ValueError(f"Attribute '{attr}' is not empty. Use force=True to overwrite existing data.")
        
        #reset some of the attributes
        #reset key dataset
        self.phenotype_df = None
        self.diagnosis = None
        self.n_diagnosis = None
        self.history = None
        self.trajectory = None
        print('Phenotype and medical records data reset.')
        #attributes related to phenotype data
        self.__warning_phenotype = []
        self.__phenotype_statistics = {}
        self.__phenotype_info = {}
        #attributes related to medical records data
        self.__warning_medical_records = []
        self.__medical_recods_statistics = {}
        self.__medical_recods_info = {}
        self._medical_recods_info = {}
        #attributes of phewas information
        self.__significant_phecodes = None

        #check input
        #----------
        self.__phenotype_info['phenotype_col_dict'] = self.__phenotype_col_dict[self.study_design]
        #check whether the required columns are in the dictionary
        value_notin = [x for x in self.__phenotype_info['phenotype_col_dict'].keys() if x not in column_names.keys()]
        if len(value_notin) > 0:
            raise ValueError(f"{value_notin} not specified in the column_names dictionary")
        #check whether unknow columns are in the dictionary
        value_unknown = [x for x in column_names.keys() if x not in self.__phenotype_info['phenotype_col_dict'].keys()]
        if len(value_unknown) > 0:
            raise ValueError(f"{value_unknown} in the column_names dictionary is not required for the study design {self.study_design}")        
        #check variable names conflict
        duplicate_cols_0 = set(column_names.values()).intersection(set(covariates)) #duplicate of original column names
        duplicate_cols_1 = set(self.__varialbe_name_place_holder).intersection(set(covariates)) #duplicate of renamed column names
        if len(duplicate_cols_0)>0:
            raise ValueError(f"The variable {duplicate_cols_0} is specified both as a required variable and in the covariates.")
        if len(duplicate_cols_1)>0:
            raise ValueError(f"The covariates {duplicate_cols_1} have duplicate names that conflict with reserved variables. Please rename them.")
        #check other columns
        all_cols = list(column_names.values())+covariates
        date_cols = [column_names[x] for x in ['Index date','End date']]
        phenotype_data_,_ = read_check_csv(phenotype_data_path,all_cols,date_cols,self.date_fmt)
        #define other atributes for phenotype data
        self.__phenotype_info['phenotype_covariates_original'] = covariates
        
        #deal with the phenotype data
        #----------
        #reset index
        phenotype_data_.reset_index(inplace=True)
        #deal with required columns
        #rename
        rename_dict = {column_names[k]:self.__phenotype_info['phenotype_col_dict'][k] for k in column_names.keys()}
        phenotype_data_.rename(columns=rename_dict,inplace=True)
        #check and convert, change inplace
        phenotype_required_columns(phenotype_data_,self.__phenotype_info['phenotype_col_dict'],self.date_fmt,self.study_design)
        #convert covariates
        self.__phenotype_info['phenotype_covariates_converted'] = {}
        self.__phenotype_info['phenotype_covariates_list'] = []
        self.__phenotype_info['phenotype_covariates_type'] = {self.__sex_col:'categorical'}
        for var in self.__phenotype_info['phenotype_covariates_original']:
            converted_df,var_type = convert_column(phenotype_data_[[var]],var)
            self.__phenotype_info['phenotype_covariates_converted'][var] = list(converted_df.columns)
            self.__phenotype_info['phenotype_covariates_list'] += list(converted_df.columns)
            self.__phenotype_info['phenotype_covariates_type'][var] = var_type
            phenotype_data_ = pd.merge(phenotype_data_,converted_df,left_index=True, right_index=True)
        #assign to attributes
        all_cols = list(self.__phenotype_info['phenotype_col_dict'].values()) + self.__phenotype_info['phenotype_covariates_list'] + self.__phenotype_info['phenotype_covariates_original']
        self.phenotype_df = phenotype_data_[all_cols] #all other columns were not remained
        del phenotype_data_#save memory
        #remove
        n_before_any_remove = len(self.phenotype_df)
        n_before_remove = n_before_any_remove
        for var in self.__phenotype_info['phenotype_covariates_converted']:
            self.phenotype_df = self.phenotype_df.dropna(subset=self.__phenotype_info['phenotype_covariates_converted'][var])
            n_after_remove = len(self.phenotype_df)
            if n_after_remove < n_before_remove:
                self.__warning_phenotype.append(f'Warning: {n_before_remove-n_after_remove} individuals removed due to missing values in {var}.')
                print(self.__warning_phenotype[-1])
                n_before_remove = n_after_remove
        if n_before_remove < n_before_any_remove:
            self.__warning_phenotype.append(f'Warning: {n_before_any_remove-n_before_remove} individuals removed due to missing values in covariates.')
            print(self.__warning_phenotype[-1])
        if self.study_design == "registry":
            self.phenotype_df[self.__exposure_col] = 1
            self.__phenotype_col_dict["registry"]["Exposure"] = self.__exposure_col #add exposure column for registry study
        #generate basic statistic for printing
        self.__phenotype_statistics['n_cohort'] = len(self.phenotype_df)
        self.__phenotype_statistics['n_exposed'] = len(self.phenotype_df[self.phenotype_df[self.__exposure_col]==1])
        self.__phenotype_statistics['n_unexposed'] = len(self.phenotype_df[self.phenotype_df[self.__exposure_col]==0])
        if self.__phenotype_statistics['n_exposed'] <= 5000:
            self.__warning_phenotype.append('Warning: There are too few exposed individuals for downstream disease network analysis')
            print(self.__warning_phenotype[-1])
        if self.__phenotype_statistics['n_cohort'] <= 10000:
            self.__warning_phenotype.append('Warning: The total sample size might be insufficient for PheWAS analysis')
            print(self.__warning_phenotype[-1])
        if self.study_design == 'matched cohort':
            #number of matched per
            self.__phenotype_statistics['n_per_group'] = self.phenotype_df.groupby(by=self.__mathcing_identifier_col)[self.__id_col].count().mean()
            if self.__phenotype_statistics['n_per_group'] < 2:
                self.__warning_phenotype.append(f"Warning: The average number of individuals in mathced group is small ({self.__phenotype_statistics['n_per_group']:.2f})")
                print(self.__warning_phenotype[-1])
            elif self.__phenotype_statistics['n_per_group'] > 10:
                self.__warning_phenotype.append(f"Warning: The average number of individuals in mathced group is large ({self.__phenotype_statistics['n_per_group']:.2f})")
                print(self.__warning_phenotype[-1])
            #no mathced group
            try:
                self.__phenotype_statistics['n_single_group'] = self.phenotype_df.groupby(by=self.__mathcing_identifier_col)[self.__id_col].count().value_counts().loc[1]
            except:
                self.__phenotype_statistics['n_single_group'] = 0
            if self.__phenotype_statistics['n_single_group'] > 0:
                self.__warning_phenotype.append(f"Warning: {self.__phenotype_statistics['n_single_group']} groups had only one individual")
                print(self.__warning_phenotype[-1])
        
        self.phenotype_df['follow_up'] = (self.phenotype_df[self.__end_date_col] - self.phenotype_df[self.__index_date_col]).dt.days/365.25

        self.__phenotype_statistics['avg_follow_exposed'] = self.phenotype_df[self.phenotype_df[self.__exposure_col]==1]['follow_up'].mean()
        self.__phenotype_statistics['avg_follow_unexposed'] = self.phenotype_df[self.phenotype_df[self.__exposure_col]==0]['follow_up'].mean()
        self.__phenotype_statistics['n_neg_follow_exposed'] = len(self.phenotype_df[(self.phenotype_df[self.__exposure_col]==1) & 
                                                                                    (self.phenotype_df['follow_up']<=0)])
        self.__phenotype_statistics['n_neg_follow_unexposed'] = len(self.phenotype_df[(self.phenotype_df[self.__exposure_col]==0) & 
                                                                                    (self.phenotype_df['follow_up']<=0)])
        if self.__phenotype_statistics['n_neg_follow_unexposed'] > 0 or self.__phenotype_statistics['n_neg_follow_exposed'] > 0:
            self.__warning_phenotype.append(f"Warning: {self.__phenotype_statistics['n_neg_follow_exposed']} exposed individuals and {self.__phenotype_statistics['n_neg_follow_unexposed']} unexposed individuals have negative or zero follow-up time.\nConsider removing them before merge.")
            print(self.__warning_phenotype[-1])
    
    def Table1(self, continuous_stat_mode:str='auto') -> pd.DataFrame:
        """
        Generate a simple Table 1 from the provided phenotpe data.

        Parameters
        ----------
        continuous_stat_mode : str, default='auto'
            Specifie the method for displaying statistics for continuous variables. It accepts three modes:
                auto: Automatically determines the appropriate summary statistics based on the results of a normality test.
                normal: Treats all continuous variables as normally distributed, displaying the mean and standard deviation.
                nonnormal: Treats all continuous variables as non-normally distributed, displaying the median and interquartile range.
        
        Returns
        -------
        pd.DataFrame : Table 1 for the phenotype data.
        """
        from .utility_summary import desceibe_table
        
        #attribute check
        data_attrs = ['phenotype_df']
        for attr in data_attrs:
            if getattr(self, attr) is None:
                raise ValueError(f"Attribute '{attr}' is empty.")
        
        #check continuous_stat_mode
        if continuous_stat_mode not in ['auto','normal','nonnormal']:
            raise ValueError("The input 'continuous_stat_mode' must be one of 'auto', 'normal', or 'nonnormal'.")
        
        #generate the variable list
        var_type_dict = self.__phenotype_info['phenotype_covariates_type']
        #----------all continuous variables
        continuous_vars = [var for var in var_type_dict if var_type_dict[var]=='continuous']
        #get their converted names 
        continuous_vars = [self.__phenotype_info['phenotype_covariates_converted'][var][0] for var in continuous_vars]
        continuous_vars += ['follow_up'] #add generated follow_up variable
        #----------all categorical variables
        categorical_vars = [var for var in var_type_dict if var_type_dict[var]=='categorical' or var_type_dict[var]=='binary']
        #call the table 1 function, make sure continuous variables are in the front, then categorical/binary variables
        all_tab1_vars = continuous_vars + categorical_vars
        all_tab1_vars_type = ['continuous' for _ in continuous_vars] + ['categorical' for _ in categorical_vars]
        #consider the study desgin
        if self.study_design == 'matched cohort' or self.study_design == 'cohort':
            table1 = desceibe_table(self.phenotype_df,all_tab1_vars,all_tab1_vars_type,self.__exposure_col,
                                    self.__sex_value_dict,continuous_stat_mode)
        elif self.study_design == 'registry':
            table1 = desceibe_table(self.phenotype_df,all_tab1_vars,all_tab1_vars_type,self.__exposure_col,
                                    self.__sex_value_dict,continuous_stat_mode,group_var_value=[1])
        return table1

    def merge_medical_records(
        self, 
        medical_records_data_path:str, 
        diagnosis_code:str, 
        column_names:dict, 
        date_fmt:str=None, 
        chunksize:int=1000000) -> None:
        """
        Merge the loaded phenotype data with one or more medical records data.
        If you have multiple medical records data to merge (e.g., with different diagnosis code types), you can call this function multiple times.

        Parameters
        ----------
        medical_records_data_path : str
            The file path containing medical records data in CSV or TSV format. 
            Required columns include Participant ID, Diagnosis code (with the specified type), and Date of diagnosis (default format '%Y-%m-%d').
        
        diagnosis_code : str
            Diagnosis code type used in the medical records data. Valid options include 'ICD-9-CM', 'ICD-9-WHO', 'ICD-10-CM', and 'ICD-10-WHO'. 
            Specify only one format type per medical records data. 
            Note: Mixing CM/WHO codes or ICD-9/10 codes within the same file is not allowed, you need to split the data first by yourself. 
            If you use codes from other systems like ICD-8 or ICD-11, you must first mapped them to ICD-9 or ICD-10 format first by yourself.
            Format of ICD-10 codes: either with decimal point (F32.1) or without (F321).
            Format of ICD-9 codes: either decimal format (9.9) or non-decimal "short" format (0099).
        
        column_names : dict
            Maps required variable names to their corresponding identifiers in the medical records dataset. 
            Required keys include 'Participant ID', 'Diagnosis code', and 'Date of diagnosis'. 
            Example:
            column_names={'Participant ID': 'eid',
                          'Diagnosis code': 'ICD',
                          'Date of diagnosis': 'date'}
            Ensure that 'Participant ID' in medical records matches that in the phenotype data. 
            The 'Diagnosis code' should correspond to the specified ICD type (CM or WHO ICD-9/10).
            Diagnosis dates must follow the format '%Y-%m-%d' unless otherwise specified. 
            Records with missing required data will be excluded.
        
        date_fmt : str, (default to use the same format as phenotype data)
            The format of the date fields in your medical records data.
        
        chunksize : int, defalut=1,000,000
            Read medical records data in chunks. This is useful for large dataset.
            You can increase this value if you get enough memory.
        
        Returns
        -------
        None
            This function modifies the object's main data attribute in-place.

        """
        #attribute check
        data_attrs = ['phenotype_df']
        for attr in data_attrs:
            if getattr(self, attr) is None:
                raise ValueError(f"Attribute '{attr}' is empty.")
        
        if len(self.__medical_recods_info) > 0:
             print(f'{len(self.__medical_recods_info)} medical records data already merged, merging with a new one.')
             self._medical_recods_info = {}
        else:
            #reset
            self.__warning_medical_records = []
            self.__medical_recods_info = {}
            self.__medical_recods_statistics = {}
            self._medical_recods_info = {}
            self.diagnosis = None
            self.n_diagnosis = None
            self.history = None
            
        #check diagnosis code
        if diagnosis_code not in self.__diagnosis_code_options:
            raise ValueError(f"The diagnosis code {diagnosis_code} is not supported, try to map it to one of them: {self.__diagnosis_code_options}")
        self._medical_recods_info['diagnosis_code'] = diagnosis_code
        #mapping files
        self._medical_recods_info['mapping_npyfile'] = os.path.join(self.__module_dir,'data/phecode_%s/%s.npy' % (self.phecode_version,self._medical_recods_info['diagnosis_code']))
        self._phecode_mapping = np.load(self._medical_recods_info['mapping_npyfile'],allow_pickle=True).item()
        #check date format
        if date_fmt is None:
            self._medical_recods_info['date_fmt'] = self.date_fmt #use defalut date format
        else:
            try:
                datetime.strftime(datetime.now(),date_fmt)
                self._medical_recods_info['date_fmt'] = date_fmt
            except:
                raise ValueError("The specified date format is invalid. Visit https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior for details.")
        
        #for medical records data
        #----------
        #check input
        if medical_records_data_path in self.__medical_recods_info:
            raise ValueError(f"{medical_records_data_path} has already been merged.")
        #check whether the required columns are in the dictionary
        value_notin = [x for x in self.__medical_records_cols if x not in column_names.keys()]
        if len(value_notin) > 0:
            raise ValueError(f"{value_notin} not specified in the column_names dictionary")
        #check columns and get seperator
        all_cols = list(column_names.values())

        date_cols = [column_names[x] for x in ['Date of diagnosis']]
        seperator = read_check_csv(medical_records_data_path,all_cols,date_cols,
                                   self._medical_recods_info['date_fmt'],return_df=False)
        self._medical_recods_info['sep'] = seperator
        self._medical_recods_info['column_names'] = column_names
        self._medical_recods_info['chunksize'] = chunksize
        
        #process the medical records data
        #----------
        self._phecode_dict = {id_:{} for id_ in self.phenotype_df[self.__id_col].values} #empty dict
        result_tuple = medical_records_process(
            medical_records_data_path,
            self._medical_recods_info['column_names'],
            self._medical_recods_info['diagnosis_code'],
            self._medical_recods_info['date_fmt'],
            self._medical_recods_info['chunksize'],
            self._medical_recods_info['sep'],
            self._phecode_dict,
            self._phecode_mapping
        )
        self._medical_recods_info['n_total_records'],self._medical_recods_info['n_total_missing'],self._medical_recods_info['n_total_trunc_4'],self._medical_recods_info['n_total_trunc_3'],self._medical_recods_info['n_total_no_mapping'],self._medical_recods_info['no_mapping_list'] = result_tuple
        #print some warning information for the current medical records dataset
        prop_missing = self._medical_recods_info['n_total_missing']/self._medical_recods_info['n_total_records']
        if prop_missing >= 0.001:
            self.__warning_medical_records.append(f"Warning: {prop_missing*100:.2f}% of diagnosis records have missing values for file {medical_records_data_path}.")
            print(self.__warning_medical_records[-1])
        prop_nomapping = self._medical_recods_info['n_total_no_mapping']/self._medical_recods_info['n_total_records']
        if prop_nomapping >= 0.001:
            self.__warning_medical_records.append(f"Warning: {prop_nomapping*100:.2f}% of {self._medical_recods_info['diagnosis_code']} codes were not mapped to phecodes for file {medical_records_data_path}.")
        
        #if not existing create new one
        if len(self.__medical_recods_info)==0:
            self.diagnosis = {id_:{} for id_ in self.phenotype_df[self.__id_col].values}
            self.n_diagnosis = {id_:{} for id_ in self.phenotype_df[self.__id_col].values}
            self.history = {id_:[] for id_ in self.phenotype_df[self.__id_col].values}
        #update new or exsiting diagnosis/diagnosis_n and history data
        self._medical_recods_info['n_invalid'] = diagnosis_history_update(self.diagnosis,self.n_diagnosis,self.history,
                                                                          dict(self.phenotype_df[[self.__id_col,self.__index_date_col]].values),
                                                                          dict(self.phenotype_df[[self.__id_col,self.__end_date_col]].values),
                                                                          self._phecode_dict)
        del self._phecode_dict #save memory
        self.__medical_recods_info[medical_records_data_path] = self._medical_recods_info
        print(f"Phecode diagnosis records successfully merged ({sum(self._medical_recods_info['n_invalid'].values()):,} invalid records were not merged, typically with diagnosis date later than date of follow-up end)\n")
        
        #generate basic statistics
        self.__medical_recods_statistics['n_merged_files'] = len(self.__medical_recods_info)
        self.__medical_recods_statistics['n_merged_records'] = sum([self.__medical_recods_info[x]['n_total_records'] for x in self.__medical_recods_info])
        self.__medical_recods_statistics['n_missing'] = sum([self.__medical_recods_info[x]['n_total_missing'] for x in self.__medical_recods_info])
        #number of diagnosis and history
        exposed_id = self.phenotype_df[self.phenotype_df[self.__exposure_col]==1][self.__id_col].values
        self.__medical_recods_statistics['n_phecode_diagnosis_per_exposed'] = np.mean([len(self.diagnosis[id_]) for id_ in exposed_id])
        self.__medical_recods_statistics['n_phecode_history_per_exposed'] = np.mean([len(self.history[id_]) for id_ in exposed_id])
        unexposed_id = self.phenotype_df[self.phenotype_df[self.__exposure_col]==0][self.__id_col].values
        self.__medical_recods_statistics['n_phecode_diagnosis_per_unexposed'] = np.mean([len(self.diagnosis[id_]) for id_ in unexposed_id])
        self.__medical_recods_statistics['n_phecode_history_per_unexposed'] = np.mean([len(self.history[id_]) for id_ in unexposed_id])

    def get_attribute(
        self, 
        attr_name:str
    ) -> None:
        """
        Retrieves the value of a specified attribute, providing controlled access to the class's private and protected data.
        
        Parameters
        ----------
        attr_name : str
            The name of the attribute to retrieve. The name should correspond to one of the predefined keys in the private_attrs dictionary.
            If the requested attribute is not found, a ValueError is raised.
        
        Returns
        -------
        value : any
            The value of the requested attribute. If the attribute is a mutable data type (like dictionaries), a copy of the data is returned
            to prevent accidental modification of the internal state.

        """

        # Use a dictionary to map attribute names to their private counterparts
        private_attrs = {
            'warning_phenotype': self.__warning_phenotype,
            'phenotype_statistics': self.__phenotype_statistics,
            'phenotype_info': self.__phenotype_info,
            'warning_medical_records': self.__warning_medical_records,
            'medical_records_statistics': self.__medical_recods_statistics,
            'medical_records_info': self.__medical_recods_info,
            'module_dir':self.__module_dir,
            'significant_phecodes':self.__significant_phecodes
        }
        if attr_name in private_attrs:
            value = private_attrs[attr_name]
            if isinstance(value, dict):
                return value.copy()  # Return a copy if it's a dictionary
            return value
        else:
            raise ValueError(f"Attribute {attr_name} not found")
    
    def modify_phecode_level(
        self, 
        phecode_level:int
    ) -> None:
        """
        Modify the phecode level setting.
        
        Parameters
        ----------
        phecode_level : int
            The level of phecode to use for analysis, where level 1 (with a total of 585 medical conditions) corresponds to 3-digit ICD-10 codes and level 2 (a total of 1257 medical conditions) to 4-digit ICD-10 codes. 
            Level 2 phecodes offer a more granular analysis with potentially smaller sample sizes per disease category. 
            For larger studies, level 2 phecodes may enhance result interpretation. 
            For smaller studies, level 1 is recommended to maintain statistical power.

        Returns
        -------
        None.
            This function modifies the object's main data attribute in-place.
        """
        #attribute check
        data_attrs = ['trajectory']
        for attr in data_attrs:
            if getattr(self, attr) is not None:
                raise ValueError(f"Attribute '{attr}' is not empty, the phecode level cannot be modified.")
        
        if phecode_level not in self.__phecode_level_options:
            raise ValueError(f"Choose from the following phecode level {self.__phecode_level_options}")
        
        if phecode_level == self.phecode_level:
            raise ValueError("No modification needed.")
        else:
            self.phecode_level = phecode_level
            #load necessary phecode files
            self.__phecode_info_npyfile = os.path.join(self.__module_dir,f'data/phecode_{self.phecode_version}/level{self.phecode_level}_info.npy')
            self.phecode_info = np.load(self.__phecode_info_npyfile,allow_pickle=True).item()
            print(f'Phecode level set to {self.phecode_level} now.')
    
    
    def disease_pair(
        self, 
        phewas_result:pd.DataFrame, 
        min_interval_days:int=0, 
        max_interval_days:int=np.inf, 
        force:bool=False, 
        **kwargs
    ) -> None:
        """
        This function reads PheWAS results from a DataFrame generated by the 'DiseaseNet.pheewas' method. 
        It filters for phecodes with significant associations and constructs all possible temporal (D1 → D2, where D2 is diagnosed after D1) and non-temporal disease pairs (D1 - D2) for each exposed individual, based on the identified significant phecodes.
        
        Parameters
        ----------
        phewas_result : pd.DataFrame
            DataFrame containing PheWAS analysis results produced by the 'DiseaseNet.pheewas' function.
        
        min_interval_days : int/float, default=0
            Minimum required time interval (in days) between diagnosis dates when constructing temporal D1 → D2 disease pair for each individual.
            Individuals with D1 and D2 diagnoses interval less than or equal to this value are considered to have non-temporal D1-D2 disease pair (without a clear sequence).
        
        max_time_interval_days : int/float, default=np.inf
            Maximum allowed time interval (in days) between diagnosis dates when constructing temporal and non-temporal D1-D2 disease pair for each individual.
            Individuals with D1 and D2 diagnoses interval greater than this value are considered to have either temporal or non-temporal D1-D2 disease pair, although they were diagnosed with both D1 and D2.
    
        force : bool, default=False
            If True, the data will be loaded and existing attributes will be overwritten, even if they contain data. 
            The default is False, which will raise an error if data already exists.
        
        **kwargs
            Additional keyword argument to define the required columns in 'phewas_result':
                phecode_col : str, default='phecode'
                    Name of the column in 'phewas_result' that specifies the phecode identifiers.
                significance_col : str, default='phewas_p_significance'
                    Name of the column in 'phewas_result' that indicates the significance of each phecode in the PheWAS analysis.
        
        Returns
        -------
        None.
            This function modifies the object's main data attribute in-place.

        """
        #attribute check
        if not force:
            data_attrs = ['trajectory']
            for attr in data_attrs:
                if getattr(self, attr) is not None:
                    raise ValueError(f"Attribute '{attr}' is not empty. Use force=True to overwrite existing data.")
        
        data_attrs = ['phenotype_df', 'diagnosis', 'n_diagnosis', 'history']
        for attr in data_attrs:
            if getattr(self, attr) is None:
                raise ValueError(f"Attribute '{attr}' is empty.")

        #check phewas result
        if not isinstance(phewas_result,pd.DataFrame):
            raise TypeError("The provided input 'phewas_result' must be a pandas DataFrame.")
        
        # Default column names
        phecode_col = kwargs.get('phecode_col', 'phecode')
        significance_col = kwargs.get('significance_col', 'phewas_p_significance')
        # Check column existence
        for col in [phecode_col, significance_col]:
            if col not in phewas_result.columns:
                raise ValueError(f"Column {col} not in 'phewas_result' DataFrame.")
        
        #check min and max interval 
        if not isinstance(min_interval_days, (float,int)):
            raise TypeError("The provided input 'min_interval_days' must be a int")
        if min_interval_days<0:
            raise ValueError("The provided input 'min_interval_days' is not valide.")
        if not isinstance(max_interval_days, (float,int)):
            raise TypeError("The provided input 'max_interval_days' must be a int")
        if max_interval_days<0 or max_interval_days<=min_interval_days:
            raise ValueError("The provided input 'max_interval_days' is not valide.")
        
        #get list of significant phecodes
        significant_phecodes = phewas_result[phewas_result[significance_col]==True][phecode_col].to_list()
        invalid_phecodes = [x for x in significant_phecodes if x not in self.phecode_info]
        if len(invalid_phecodes)>0:
            print(f'Warning: the following phecodes are invalid and ignored (possible phecode level disconcordance): {invalid_phecodes}')
        valid_phecodes = [x for x in significant_phecodes if x in self.phecode_info]
        if len(valid_phecodes) == 0:
            raise ValueError("No significant phecodes from 'phewas_result' are found.")
        if len(valid_phecodes) <= 10:
            print(f'Warning: only {len(valid_phecodes)} significant phecodes are found.')
        self.__significant_phecodes = valid_phecodes
        #for this subgroup of people only
        exp_col = self.__phenotype_info['phenotype_col_dict']['Exposure']
        exposed_index = self.phenotype_df[self.phenotype_df[exp_col]==1].index

        self.trajectory = d1d2_from_diagnosis_history(
            self.phenotype_df.loc[exposed_index],
            self.__phenotype_info['phenotype_col_dict']['Participant ID'],
            self.__phenotype_info['phenotype_col_dict']['Sex'],
            self.__sex_value_dict,
            self.__significant_phecodes,
            self.history,
            self.diagnosis,
            self.n_diagnosis,
            self.phecode_info,
            min_interval_days,
            max_interval_days,
            self.min_required_icd_codes
        )

    def load(
        self, 
        file:str, 
        force:bool=False
    ) -> None:
        """
        Load data from a gzip-compressed pickle file and restore the attributes to this DiseaseNet.DiseaseNetworkData object.
        This method is intended for restoring data to an empty object. 
        If data is already present in any attribute and 'force' is not set to True, an error will be raised to prevent accidental data overwrite.

        Parameters
        ----------
        file : str
            The filename (string) from which the data object will be loaded.
            The '.pkl.gz' extension will be automatically appended if not already included.
        
        force : bool, default=False
            If True, the data will be loaded and existing attributes will be overwritten, even if they contain data. 
            The default is False, which will raise an error if data already exists.

        Returns
        -------
        None.
    
        """
        import pickle
        import gzip
        import gc
        
        # Check for existing data if force is not True
        if not force:
            data_attrs = ['phenotype_df', 'diagnosis', 'n_diagnosis', 'history', 'trajectory']
            for attr in data_attrs:
                if getattr(self, attr) is not None:
                    raise ValueError(f"Attribute '{attr}' is not empty. Use force=True to overwrite existing data.")

        # Add '.pkl.gz' extension if not present
        if not file.endswith('.pkl.gz'):
            file += '.pkl.gz'
        # Load the dictionary from .npy file
        with gzip.open(file, 'rb') as f:
            data_dict = pickle.load(f)
        # Restore the pandas DataFrame attribute
        self.phenotype_df = data_dict.pop('phenotype_df')
        # Restore all simple attributes directly from data_dict
        simple_attrs = ['study_design', 'date_fmt', 'phecode_level', 'phecode_version','min_required_icd_codes',
                        'phecode_info', 'diagnosis', 'n_diagnosis', 'history', 'trajectory']
        for attr in simple_attrs:
            setattr(self, attr, data_dict.pop(attr, None))
        # Restore all private attributes from data_dict
        private_attrs = ['__warning_phenotype', '__phenotype_statistics', '__phenotype_info',
                         '__warning_medical_records', '__medical_recods_statistics', '__medical_recods_info',
                         '__module_dir', '__significant_phecodes']
        for attr in private_attrs:
            setattr(self, f'_DiseaseNetworkData{attr}', data_dict.pop(attr, None))
        # Clear the remaining data_dict to free memory
        del data_dict
        gc.collect()

        print("All attributes restored.")

    def save(
        self,
        file:str
    ) -> None:
        """
        Save the DiseaseNet.DiseaseNetworkData object's attributes to a gzip-compressed pickle file,
        which can be restored using the corresponding load method.
    
        Parameters
        ----------
        file : str
            The filename or path prefix where the data object will be saved. 
            The '.pkl.gz' extension will be automatically appended if not already included.
    
        Returns
        -------
        None.
        """
        import pickle
        import gzip

        #attribute check
        data_attrs = ['phenotype_df']
        for attr in data_attrs:
            if getattr(self, attr) is None:
                raise ValueError(f"Attribute '{attr}' is empty, nothing to save")
        
        data_attrs = ['diagnosis','history','trajectory']
        for attr in data_attrs:
            if getattr(self, attr) is None:
                print(f"Attribute '{attr}' is empty.") 
        
        # Create a dictionary to save the data
        save_dict = {
            'study_design': self.study_design,
            'date_fmt': self.date_fmt,
            'phecode_level': self.phecode_level,
            'phecode_version': self.phecode_version,
            'min_required_icd_codes': self.min_required_icd_codes,
            'phecode_info': self.phecode_info,
            'phenotype_df': self.phenotype_df,
            'diagnosis': self.diagnosis,
            'n_diagnosis': self.n_diagnosis,
            'history': self.history,
            'trajectory': self.trajectory,
            '__warning_phenotype': self.__warning_phenotype,
            '__phenotype_statistics': self.__phenotype_statistics,
            '__phenotype_info': self.__phenotype_info,
            '__warning_medical_records': self.__warning_medical_records,
            '__medical_recods_statistics': self.__medical_recods_statistics,
            '__medical_recods_info': self.__medical_recods_info,
            '__module_dir':self.__module_dir,
            '__significant_phecodes':self.__significant_phecodes
        }
        
        # Add '.pkl.gz' extension if not present
        if not file.endswith('.pkl.gz'):
            file += '.pkl.gz'
        
        #save it
        with gzip.open(file, 'wb') as f:
            pickle.dump(save_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
        
        print(f"DiseaseNetworkData is saved to {file}.")

    def __str__(self):
        self.__print_string = 'DiseseNet.DiseaseNetworkData\n\n'
        #add cohort description
        self.__print_string += f'Study design: {self.study_design}\n'
        if self.phenotype_df is not None:
            self.__print_string += '\nPhentype data\n'
            self.__print_string += f"Total number of individuals: {self.__phenotype_statistics['n_cohort']:,} ({self.__phenotype_statistics['n_exposed']:,} exposed and {self.__phenotype_statistics['n_unexposed']:,} unexposed)\n"
            if self.study_design == 'matched cohort':
                self.__print_string += f"Average number of individuals in mathced group: {self.__phenotype_statistics['n_per_group']:.2f}\n"
            self.__print_string += f"Average follow-up years: {self.__phenotype_statistics['avg_follow_exposed']:.2f} (exposed) and {self.__phenotype_statistics['avg_follow_unexposed']:.2f} (unexposed)\n"
        if len(self.__medical_recods_info) > 0:
            self.__print_string += '\nMerged Medical records\n'
            self.__print_string += f"{self.__medical_recods_statistics['n_merged_files']} medical records data with {self.__medical_recods_statistics['n_merged_records']:,} diagnosis records were merged ({self.__medical_recods_statistics['n_missing']:,} with missing values).\n"
            self.__print_string += f"Average number of disease diagnosis during follow-up: {self.__medical_recods_statistics['n_phecode_diagnosis_per_exposed']:.2f} (exposed) and {self.__medical_recods_statistics['n_phecode_diagnosis_per_unexposed']:.2f} (unexposed)\n"
            self.__print_string += f"Average number of disease diagnosis before follow-up: {self.__medical_recods_statistics['n_phecode_history_per_exposed']:.2f} (exposed) and {self.__medical_recods_statistics['n_phecode_history_per_unexposed']:.2f} (unexposed)\n"
        if len(self.__warning_phenotype) > 0:
            self.__print_string += '\n'
            self.__print_string += '\n'.join(self.__warning_phenotype)
        if len(self.__warning_medical_records) > 0:
            self.__print_string += '\n'
            self.__print_string += '\n'.join(self.__warning_medical_records)
        return self.__print_string
    
    def __repr__(self):
        return self.__str__()

    def __len__(self):
        if self.phenotype_df is not None:
            return self.__phenotype_statistics['n_cohort']
        else:
            return 0
    

    

    





















    
    
    
    
    
    
    
    