# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 23:48:05 2024

@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""

import pandas as pd
import numpy as np
import os
from datetime import datetime
from .utility import convert_column,phenotype_required_columns,diff_date_years,read_check_csv
from .utility import medical_records_process,diagnosis_history_update


class DiseaseNetworkData:
    """
    A class for handling disease network data creation and operations, for use in DiseaseNet module.

    Parameters:
    ----------
    study_design : str
        Specify the type of study design, either "cohort" or "matched cohort".
    
    phecode_level : int
        The level of phecode to use for analysis, where level 1 (with a total of 585 medical conditions) corresponds to 3-digit ICD-10 codes and level 2 (a total of 1257 medical conditions) to 4-digit ICD-10 codes.
        Level 2 phecodes offer a more granular analysis with potentially smaller sample sizes per disease category. 
        For larger studies, level 2 phecodes may enhance result interpretation. 
        For smaller studies, level 1 is recommended to maintain statistical power.
    
    date_fmt : str, default='%Y-%m-%d'
        The format of the date fields in your phenotype and medical records data.
    
    phecode_version : str, default='1.2'
        The version of the phecode system used for converting diagnosis codes. Currently, only version 1.2 is supported.
        
    chunksize_medical_records : int, default=100
        Passed to pd.read_csv chunksize for reading medical records data with iterator mode.
        If this value is 0, the iterator mode in pd.read_csv will not be used.
    
    """
    
    def __init__(self,study_design:str='cohort', phecode_level:int=1, date_fmt:str='%Y-%m-%d',phecode_version:str='1.2'):
        #fixed attributes
        #phenotype data
        self.__module_dir = os.path.dirname(__file__)
        self.__test_date = '2023-1-1'
        self.__study_design_options = ['matched cohort','cohort']
        self.__id_col = 'eid'
        self.__exposure_col = 'exposure'
        self.__index_date_col = 'index_date'
        self.__end_date_col = 'end_date'
        self.__mathcing_identifier_col = 'group'
        self.__phenotype_col_dict = {'matched cohort':{'Participant ID':self.__id_col,
                                                     'Exposure':self.__exposure_col,
                                                     'Index date': self.__index_date_col,
                                                     'End date': self.__end_date_col,
                                                     'Matching identifier':self.__mathcing_identifier_col},
                                     'cohort':{'Participant ID':self.__id_col,
                                               'Exposure':self.__index_date_col,
                                               'Index date': self.__index_date_col,
                                               'End date': self.__end_date_col,}}
        #medical records data
        self.__diagnosis_code_options = ['ICD-9-CM', 'ICD-9-WHO', 'ICD-10-CM', 'ICD-10-WHO']
        self.__phecode_level_options = [1,2]
        self.__phecode_version_options = ['1.2']
        self.__medical_records_cols = ['Participant ID','Diagnosis code','Date of diagnosis']
        self.__medical_recods_info = {}
        
        #check study design
        if not study_design in self.__study_design_options:
            raise ValueError(f"Choose from the following study design: {self.__study_design_options}")
        #check date format
        try:
            datetime.strftime(datetime.now(),date_fmt)
        except:
            raise ValueError("The specified date format is invalid. Visit https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior for details.")
        #check phecode_level
        if phecode_level not in self.__phecode_level_options:
            raise ValueError(f"Choose from the following phecode level {self.__phecode_level_options}")
        #check phecode version
        if phecode_version not in self.__phecode_version_options:
            raise ValueError(f"Supported phecode version {self.__phecode_version_options}")
        
        #assign parameter attributes
        self.study_design = study_design
        self.date_fmt = date_fmt
        self.phecode_level = phecode_level
        self.phecode_version = phecode_version
        #load necessary phecode files
        self.__phecode_info_npyfile = os.path.join(self.__module_dir,f'data/phecode_{self.phecode_version}/level{self.phecode_level}_info.npy')
        self.phecode_info = np.load(self.__phecode_info_npyfile,allow_pickle=True).item()
        #reset key attributes
        self.phenotype_df = None
        self.diagnosis = None
        self.history = None
        self.trajectory = None
        self.__warning_phenotype = []
        self.__warning_medical_records = []

    
    def phenotype_data(self, phenotype_data_path:str, column_names:dict, covariates:list):
        """
        
        Merges phenotype and medical records data into the main data attribute.
        
        Parameters:
        ----------
        phenotype_data_path : str
            The file path containing phenotype data in CSV or TSV format. The file must include a header row with columns such as Participant ID, 
            Index date, End date, Exposure, Matching identifier (if a matched cohort design), and other covariates like age and sex.
                
        column_names : dict
            A dictionary mapping required column names to their corresponding identifiers in the dataset. 
            Expected keys include 'Participant ID', 'Index date', 'End date', 'Exposure', and 'Matching identifier' (if applicable). 
            For example:
            column_names={'Participant ID': 'eid',
                          'Exposure': 'status',
                          'Index date': 'index_date',
                          'End date': 'final_date',
                          'Matching identifier': 'group'}
            The 'Exposure' column should be coded as 0 (unexposed) and 1 (exposed). Dates must be formatted as '%Y-%m-%d' unless specified otherwise. 
            Records with missing values in any required columns are not allowed.
        
        covariates : list
            A list of column names representing additional covariates, such as ['age', 'sex', 'BMI']. 
            Provide an empty list if no additional covariates are included. 
            The system will automatically detect and convert variable types. 
            Individuals with missing values in continuous variables will be removed, while those missing in categorical variables will be categorized separately.

        Returns:
        ----------
        None
            This function modifies the object's main data attribute in-place.
        
        """
        #reset some of the attributes
        #reset key dataset
        self.phenotype_df = None
        self.diagnosis = None
        self.history = None
        self.trajectory = None
        print('Phenotype and medical records data reset.')
        #attributes related to phenotype data
        self.__warning_phenotype = []
        self.__warning_phenotype = []
        self.__phenotype_statistics = {}
        self.__phenotype_info = {}
        #attributes related to medical records data
        self.__warning_medical_records = []
        self.__medical_recods_statistics = {}
        self.__medical_recods_info = {}
        self._medical_recods_info = {}


        #check input
        #----------
        self.__phenotype_info['phenotype_col_dict'] = self.__phenotype_col_dict[self.study_design]
        #check whether the required columns are in the dictionary
        value_notin = [x for x in self.__phenotype_info['phenotype_col_dict'].keys() if x not in column_names.keys()]
        if len(value_notin) > 0:
            raise ValueError(f"{value_notin} not specified in the column_names dictionary")
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
        phenotype_required_columns(phenotype_data_,self.__phenotype_info['phenotype_col_dict'],self.date_fmt)
        #convert covariates
        self.__phenotype_info['phenotype_covariates_converted'] = {}
        self.__phenotype_info['phenotype_covariates_list'] = []
        self.__phenotype_info['phenotype_covariates_type'] = {}
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
        n_before_remove = len(self.phenotype_df)
        self.phenotype_df = self.phenotype_df.dropna(how='any')
        n_after_remove = len(self.phenotype_df)
        if n_after_remove < n_before_remove:
            self.__warning_phenotype.append(f'Warning: {n_before_remove-n_after_remove} individuals removed after processing')
            print(self.__warning_phenotype[-1])
        
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
            
        self.phenotype_df['follow_up'] = diff_date_years(self.phenotype_df[[self.__index_date_col,self.__end_date_col]],
                                                         self.__index_date_col,self.__end_date_col)
        self.__phenotype_statistics['avg_follow_exposed'] = self.phenotype_df[self.phenotype_df[self.__exposure_col]==1]['follow_up'].mean()
        self.__phenotype_statistics['avg_follow_unexposed'] = self.phenotype_df[self.phenotype_df[self.__exposure_col]==0]['follow_up'].mean()
        self.__phenotype_statistics['n_neg_follow_exposed'] = len(self.phenotype_df[(self.phenotype_df[self.__exposure_col]==1) & 
                                                                                    (self.phenotype_df['follow_up']<=0)])
        self.__phenotype_statistics['n_neg_follow_unexposed'] = len(self.phenotype_df[(self.phenotype_df[self.__exposure_col]==0) & 
                                                                                      (self.phenotype_df['follow_up']<=0)])
        if self.__phenotype_statistics['n_neg_follow_unexposed'] > 0 or self.__phenotype_statistics['n_neg_follow_exposed'] > 0:
            self.__warning_phenotype.append(f"Warning: {self.__phenotype_statistics['n_neg_follow_exposed']} exposed individuals and {self.__phenotype_statistics['n_neg_follow_unexposed']} unexposed individuals have negative or zero follow-up time.\nConsider removing them before merge.")
            print(self.__warning_phenotype[-1])
        #----------
       
    def merge_medical_records(self, medical_records_data_path:str, diagnosis_code:str, 
                              column_names:dict, date_fmt:str=None, chunksize:int=1000000):
        """
        Merge the loaded phenotype data with one more medical records data.
        If you have multiple medical records data to merge (in most cases, with difference diagnosis code types), you can call this function multiple times.

        Parameters
        ----------
        medical_records_data_path : str
            The file path containing medical records data in CSV or TSV format. Required columns include Participant ID, Diagnosis code (with the specified type), 
            and Date of diagnosis (default format '%Y-%m-%d').
        
        diagnosis_code : str
            Diagnosis code type used in the medical records data. Valid options include 'ICD-9-CM', 'ICD-9-WHO', 'ICD-10-CM', and 'ICD-10-WHO'. 
            Specify only one format type per medical records data. 
            Note: Mixing CM/WHO codes or ICD-9/10 codes within the same file is not allowed, you need to split the data first by yourself. 
            If you use codes from other systems like ICD-8 or ICD-11, you must first mapped them to ICD-9 or ICD-10 format first by yourself.
            Format of ICD-10 codes: either with decimal point (F32.1) or without (F321).
            Format of ICD-9 codes: either decimal format (9.9) or non-decimal "short" format (0099).
        
        column_names : dict
            Maps required column names to their corresponding identifiers in the medical records dataset. 
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
        result_tuple = medical_records_process(medical_records_data_path,
                                               self._medical_recods_info['column_names'],
                                               self._medical_recods_info['diagnosis_code'],
                                               self._medical_recods_info['date_fmt'],
                                               self._medical_recods_info['chunksize'],
                                               self._medical_recods_info['sep'],
                                               self._phecode_dict,
                                               self._phecode_mapping)
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
            self.history = {id_:[] for id_ in self.phenotype_df[self.__id_col].values}
        #update new or exsiting diagnosis and history data
        self._medical_recods_info['n_invalid'] = diagnosis_history_update(self.diagnosis,self.history,
                                                                          dict(self.phenotype_df[[self.__id_col,self.__index_date_col]].values),
                                                                          dict(self.phenotype_df[[self.__id_col,self.__end_date_col]].values),
                                                                          self._phecode_dict)
        del self._phecode_dict #save memory
        self.__medical_recods_info[medical_records_data_path] = self._medical_recods_info
        print(f"Phecode diagnosis records successfully merged ({sum(self._medical_recods_info['n_invalid'].values()):,} invalid records were not merged, typically with diagnosis date later than date of follow-up end)\n")
        print("To obtain more detailed information about the merged medical records data, please use the '.get_medical_records_info()' method on the instance of DiseaseNetworkData you created.")
        
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

    def load(self, filepath):
        None

    def save(self, filepath):
        None
    
    def get_phenotype_info(self):
        return self.__phenotype_info
    
    def get_medical_records_info(self):
        return self.__medical_recods_info
    
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
    

    

    





















    
    
    
    
    
    
    
    