# -*- coding: utf-8 -*-
"""
Created on Wed Jan 1 19:48:09 2025

@author: Haowen Liu - Biomedical Big data center of West China Hospital, Sichuan University
@author: Can Hou - Biomedical Big data center of West China Hospital, Sichuan University
"""
import matplotlib.ticker as ticker
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import plotly.offline as py
import matplotlib.cm as cm
import matplotlib as mpl
import networkx as nx
import pandas as pd
import numpy as np
import itertools
import random
import math
import os
import warnings
import gc
from networkx.algorithms import community as nx_community

from collections import Counter

from typing import (
    Any,
    List,
    Dict,
    Tuple,
    Optional,
)

Df = pd.DataFrame

class Plot(object):
    """Initialize the Plot class.

    This class integrates and visualizes disease relationships from three
    complementary analyses:
    1. Phenome-Wide Association Study (PHEWAS) results
    2. Comorbidity network analysis
    3. Disease trajectory analysis

    Two visualization styles are supported:
    - After-exposure analysis: specify ``exposure_name`` to show the exposure as
      the origin node of downstream disease trajectories.
    - Before-outcome analysis: specify ``outcome_name`` instead of
      ``exposure_name`` to show a terminal outcome node below the disease
      network in the 3D and trajectory plots.

    Disease modules are shared across related network plots. If both
    comorbidity_result and trajectory_result are provided, the module graph is
    prepared from the comorbidity network plus trajectory-only disease pairs.
    Modules are detected by repeated direct NetworkX Louvain runs and the
    partition with highest NetworkX modularity is retained. Node positions for
    3D and comorbidity plots use bounded random placement within cluster
    sectors, matching the traditional visual style while avoiding long retry
    loops.

    Args:
        comorbidity_result (pd.DataFrame, optional): 
            Result dataframe from comorbidity network analysis containing:
            - Non-temporal disease pairs (D1-D2)
            - Association metrics (e.g., beta coefficients, p-values)
            - Significance identifier (True or False)
            Required for comorbidity network, 3D network, and trajectory plots.
            
        trajectory_result (pd.DataFrame, optional): 
            Result dataframe from temporal disease trajectory analysis containing:
            - Temporal disease pairs (source->target)
            - Temporal association metrics (e.g., beta coefficients, p-values)
            - Significance identifier (True or False)
            Required for 3D network and trajectory plots.
            
        phewas_result (pd.DataFrame): 
            Result dataframe from PHEWAS analysis containing:
            - Phecode disease
            - Effect sizes (e.g., hazard ratios)
            - Case counts
            - Disease system classifications
        
        exposure_name (str, optional):
            Identifier for the primary exposure variable of interest.
            Use for the traditional after-exposure visualization style.
            Do not specify together with outcome_name.
            
        exposure_location (Tuple[float], optional): 
            Custom 3D coordinates (x,y,z) for positioning the exposure node.
            If None, will be automatically positioned at (0,0,0).
            If exposure_name is None, this parameter is ignored.
            
        exposure_size (float, optional): 
            Relative size scaling factor for the exposure node in the plot.
            If exposure_name is None, this parameter is ignored.

        outcome_name (str, optional):
            Identifier for the terminal outcome in before-outcome analysis.
            Use instead of exposure_name for the before-outcome visualization
            style. In 3D and trajectory plots, the outcome is drawn as a
            terminal node below the disease network.

        outcome_location (Tuple[float], optional):
            Custom 3D coordinates (x,y,z) for positioning the terminal outcome node.
            If None, the outcome is automatically positioned below the disease network.

        outcome_size (float, optional):
            Relative size scaling factor for the terminal outcome node.
        
        Additional parameters:
        If your result DataFrames use the default column names, keep these parameters as is.

        phecode_col (str, optional):
            Column in the PheWAS DataFrame that stores phecodes.
            Default 'phecode'.

        disease_col (str, optional):
            Column in the PheWAS DataFrame with disease names.
            Default 'disease'.

        system_col (str):
            Column in the PheWAS DataFrame with disease system labels.
            Default 'system'.

        phewas_number_col (str, optional):
            Column in the PheWAS DataFrame with case counts.
            Default ``N_cases_exposed``. When ``outcome_name`` is supplied and
            that default column is absent, ``N_cases_outcome_group`` is selected
            automatically.

        phewas_coef_col (str, optional):
            Column in the PheWAS DataFrame with effect sizes.
            Default 'phewas_coef'.

        phewas_se_col (str, optional):
            Column in the PheWAS DataFrame with standard errors.
            Default 'phewas_se'.

        source_col (str, optional):
            Column in comorbidity_result and trajectory_result for source or antecedent diseases.
            Default 'phecode_d1'.

        target_col (str, optional):
            Column in comorbidity_result and trajectory_result for target or consequent diseases.
            Default 'phecode_d2'.

        disease_pair_col (str):
            Column in comorbidity_result and trajectory_result with disease pair identifiers.
            Default 'name_disease_pair'.

        comorbidity_beta_col (str, optional):
            Column in comorbidity_result with effect sizes.
            Default 'comorbidity_beta'.

        trajectory_beta_col (str, optional):
            Column in trajectory_result with effect sizes.
            Default 'trajectory_beta'.

        phewas_significance_col (str):
            Column in the PheWAS DataFrame used for significance filtering.
            Default 'phewas_p_significance'.

        comorbidity_significance_col (str):
            Column in comorbidity_result used for significance filtering.
            Default 'comorbidity_p_significance'.

        trajectory_significance_col (str):
            Column in trajectory_result used for significance filtering.
            Default 'trajectory_p_significance'.

        **kwargs:
        SYSTEM (List[str]):
            Use together with COLOR to assign colors by phecode system. If not provided, systems and their order are inferred from the PheWAS results. Default order:
            ['neoplasms',
            'genitourinary',
            'digestive',
            'respiratory',
            'infectious diseases',
            'mental disorders',
            'musculoskeletal',
            'hematopoietic',
            'dermatologic',
            'circulatory system',
            'neurological',
            'endocrine/metabolic',
            'sense organs',
            'injuries & poisonings',
            'congenital anomalies',
            'symptoms',
            'others']

        COLOR (List[str]):
            Use together with SYSTEM to map systems to colors, one to one. The length of COLOR must be at least the length of SYSTEM. Supported formats include 'red', '#ED9A8D', or 'rgb(255, 0, 0)'. Default order:
            ['#F46D5A',
            '#5DA5DA',
            '#5EBCD1',
            '#C1D37F',
            '#CE5A57',
            '#A5C5D9',
            '#F5B36D',
            '#7FCDBB',
            '#ED9A8D',
            '#94B447',
            '#8C564B',
            '#E7CB94',
            '#8C9EB2',
            '#E0E0E0',
            '#F1C40F',
            '#9B59B6',
            '#4ECDC4',
            '#6A5ACD']

        Notes:
        - All provided DataFrames should use consistent phecode identifiers.
        - Comorbidity network and disease trajectory results are filtered based on significance identifier and effect size (keep positive effects).
        - Node sizes default to case counts.
        - By default, colors are assigned by disease system.

    Examples:
        After-exposure visualization using standard result columns:

        >>> plot = Plot(
        ...     phewas_df,
        ...     comorbidity_df,
        ...     trajectory_df,
        ...     exposure_name="Study exposure",
        ... )

        Before-outcome visualization:

        >>> plot = Plot(
        ...     outcome_phewas_df,
        ...     outcome_comorbidity_df,
        ...     outcome_trajectory_df,
        ...     outcome_name="Terminal outcome",
        ... )

        For PheWAS-only or exposed-only results, omit network tables and the
        exposure/outcome node as appropriate:

        >>> plot = Plot(phewas_df)
    """
    def __init__(
        self, 
        phewas_result: Df,
        comorbidity_result: Optional[Df]=None, 
        trajectory_result: Optional[Df]=None,
        exposure_name: Optional[str] | None=None,
        exposure_location: Optional[Tuple[float]] | None=None,
        exposure_size: Optional[float] | None=None,
        outcome_name: Optional[str] | None=None,
        outcome_location: Optional[Tuple[float]] | None=None,
        outcome_size: Optional[float] | None=None,
        phecode_col: Optional[str]='phecode',
        disease_col: Optional[str]='disease',
        system_col: Optional[str]='system',
        phewas_number_col: Optional[str]='N_cases_exposed',
        phewas_coef_col: Optional[str]='phewas_coef',
        phewas_se_col: Optional[str]='phewas_se',
        source_col: Optional[str]='phecode_d1',
        target_col: Optional[str]='phecode_d2',
        disease_pair_col: Optional[str]='name_disease_pair',
        comorbidity_beta_col: Optional[str]='comorbidity_beta',
        trajectory_beta_col: Optional[str]='trajectory_beta',
        phewas_significance_col: Optional[str]='phewas_p_significance',
        comorbidity_significance_col: Optional[str]='comorbidity_p_significance',
        trajectory_significance_col: Optional[str]='trajectory_p_significance',
        **kwargs
    ):
        
        if not isinstance(phewas_result, pd.DataFrame):
            raise TypeError(
                f"phewas_result is NOT a pandas.DataFrame (type: {type(phewas_result)})"
            )
        if exposure_name is not None and outcome_name is not None:
            raise ValueError("Specify either 'exposure_name' or 'outcome_name', not both.")
        if (
            outcome_name is not None
            and phewas_number_col == 'N_cases_exposed'
            and phewas_number_col not in phewas_result.columns
            and 'N_cases_outcome_group' in phewas_result.columns
        ):
            phewas_number_col = 'N_cases_outcome_group'
        if outcome_name is not None and outcome_size is None:
            outcome_size = 15

        optional_variables = {
            'comorbidity_result': comorbidity_result,
            'trajectory_result': trajectory_result
        }
        for var_name, var in optional_variables.items():
            if var is None or isinstance(var, pd.DataFrame):
                continue
            raise TypeError(f"{var_name} is NOT a pandas.DataFrame (type: {type(var)})")

        # Dictionary of variables to check (name: value)
        validate_string_params = {
            'phecode_col': phecode_col,
            'disease_col': disease_col,
            'system_col': system_col,
            'phewas_number_col': phewas_number_col,
            'phewas_coef_col': phewas_coef_col,
            'phewas_se_col': phewas_se_col,
            'source_col': source_col,
            'target_col': target_col,
            'disease_pair_col': disease_pair_col,
            'comorbidity_beta_col': comorbidity_beta_col,
            'trajectory_beta_col': trajectory_beta_col,
            'phewas_significance_col': phewas_significance_col,
            'comorbidity_significance_col': comorbidity_significance_col,
            'trajectory_significance_col': trajectory_significance_col
        }
        # Check each variable's type whether is string
        for name, value in validate_string_params.items():
            if isinstance(value, str):
                #assign attribute name
                setattr(self, name, value)
                continue
            else:
                raise TypeError(f"{name} is NOT a string {type(value)}")

        self.__validate_result_columns(
            phewas_result,
            'phewas_result',
            [
                self.phecode_col,
                self.disease_col,
                self.system_col,
                self.phewas_number_col,
                self.phewas_significance_col
            ]
        )
        if self.phewas_coef_col in phewas_result.columns:
            self.__validate_result_columns(
                phewas_result,
                'phewas_result',
                [self.phewas_coef_col]
            )
        if self.phewas_se_col in phewas_result.columns:
            self.__validate_result_columns(
                phewas_result,
                'phewas_result',
                [self.phewas_se_col]
            )

        if comorbidity_result is not None:
            self.__validate_result_columns(
                comorbidity_result,
                'comorbidity_result',
                [
                    self.source_col,
                    self.target_col,
                    self.disease_pair_col,
                    self.comorbidity_beta_col,
                    self.comorbidity_significance_col
                ]
            )

        if trajectory_result is not None:
            self.__validate_result_columns(
                trajectory_result,
                'trajectory_result',
                [
                    self.source_col,
                    self.target_col,
                    self.disease_pair_col,
                    self.trajectory_beta_col,
                    self.trajectory_significance_col
                ]
            )
        
        # all phecodes to analysis
        diseases_phewas = set(phewas_result[self.phecode_col].to_list())
        if comorbidity_result is not None:
            diseases_com = (
                comorbidity_result[self.source_col].to_list() +
                comorbidity_result[self.target_col].to_list()
            )
            for disease in set(diseases_com):
                if disease not in diseases_phewas:
                    raise ValueError(
                        f"{disease} of comorbidity result is NOT in the phewas_result (pandas.DataFrame)"
                    )

        if trajectory_result is not None:
            diseases_tra = (
                trajectory_result[self.source_col].to_list() +
                trajectory_result[self.target_col].to_list()
            )
            for disease in set(diseases_tra):
                if disease not in diseases_phewas:
                    raise ValueError(
                        f"{disease} of trajectory result is NOT in the phewas_result (pandas.DataFrame)"
                    )

        # filter the results meeting some rules
        phewas_result, comorbidity_result, trajectory_result = self.__filter_significant(
            phewas_result,
            comorbidity_result,
            trajectory_result,
        )

        # if 'exposure' is in kwargs, then throw error
        if 'exposure' in kwargs.keys():
            raise ValueError("Unknown parameter 'exposure'. Did you mean 'exposure_name'?")
        # for other parameter in kwargs, other than SYSTEM and COLOR, throw error
        for key in kwargs.keys():
            if key not in ['SYSTEM', 'COLOR']:
                raise ValueError(f"Unknown parameter '{key}'")

        COLOR = [
            '#F46D5A',
            '#5DA5DA',
            '#5EBCD1',
            '#C1D37F',
            '#CE5A57',
            '#A5C5D9',
            '#F5B36D',
            '#7FCDBB',
            '#ED9A8D',
            '#94B447',
            '#8C564B',
            '#E7CB94',
            '#8C9EB2',
            '#E0E0E0',
            '#9B59B6',
            "#F1C40F",
            '#4ECDC4',
            '#6A5ACD' 
        ]
        system = list(phewas_result[self.system_col].unique())
        self.sys_dict = {
            'neoplasms':'Neoplasms', 
            'genitourinary':'Genitourinary diseases', 
            'digestive':'Digestive diseases', 
            'respiratory':'Respiratory diseases',
            'infectious diseases':'Infectious diseases', 
            'mental disorders':'Mental disorders', 
            'musculoskeletal':'Musculoskeletal diseases',
            'hematopoietic':'Hematopoietic diseases', 
            'dermatologic':'Dermatologic diseases', 
            'circulatory system':'Circulatory system diseases',
            'neurological':'Neurological diseases',
            'endocrine/metabolic':'Endocrine/metabolic diseases', 
            'sense organs':'Diseases of the sense organs',
            'injuries & poisonings': 'Injuries & poisonings',
            'congenital anomalies': 'Congenital anomalies',
            'symptoms':'Symptoms',
            'others':'Other diseases',
            'pregnancy complications':'Pregnancy complications'
        }

        SYSTEM = kwargs.get("SYSTEM", system)
        COLOR = kwargs.get("COLOR", COLOR)
        # ensure all systems in the results are included in SYSTEM
        for x in system:
            if x not in SYSTEM:
                raise ValueError(f"Phecode system {x} of is NOT in SYSTEM")
        # check whether all listed systems are valid
        for x in SYSTEM:
            if x not in self.sys_dict.keys():
                raise ValueError(f"{x} is NOT a valid phecode system")

        if len(SYSTEM) > len(COLOR):
            raise ValueError(
                f"the length of SYSTEM is more than that of COLOR"
            )
        else:
            COLOR = COLOR[0: len(SYSTEM)]
            
        system_color = dict(
            zip(
                SYSTEM,
                COLOR
            )
        )
        system_color = kwargs.get("system_color", system_color)

        self.__init_attrs(
            phewas = phewas_result,
            base_comorbidity = comorbidity_result,
            base_trajectory = trajectory_result,
            comorbidity = pd.DataFrame(),
            trajectory = pd.DataFrame(),
            exposure = None,
            exposure_name = exposure_name,
            exposure_location = exposure_location,
            exposure_size = exposure_size,
            outcome_name = outcome_name,
            outcome_location = outcome_location,
            outcome_size = outcome_size,
            source = self.source_col,
            target = self.target_col,
            describe = phewas_result[[self.phecode_col, self.disease_col, self.system_col, self.phewas_number_col]].copy(),
            commorbidity_nodes = set(),
            trajectory_nodes = set(),
            system_color = system_color,
            nodes_attrs = {},
            network_attrs = {}
        )
        self._cluster_cache = {}
        self._prepared_cluster_key = None

    @staticmethod
    def __validate_result_columns(
        result: Df,
        result_name: str,
        columns: List[str]
    ) -> None:
        for col_name in columns:
            if col_name not in result.columns:
                raise ValueError(
                    f"{col_name} is NOT a column name in the {result_name} (pandas.DataFrame)"
                )

    def __require_plot_results(
        self,
        plot_name: str,
        require_comorbidity: bool=False,
        require_trajectory: bool=False
    ) -> None:
        missing = []
        if require_comorbidity and self._base_comorbidity is None:
            missing.append("comorbidity_result")
        if require_trajectory and self._base_trajectory is None:
            missing.append("trajectory_result")

        if missing:
            missing_text = ", ".join(missing)
            raise ValueError(f"{plot_name} requires {missing_text}.")

    def __prepare_network_data(
        self,
        plot_name: str,
        require_trajectory: bool=False
    ) -> None:
        self.__require_plot_results(
            plot_name,
            require_comorbidity=True,
            require_trajectory=require_trajectory
        )

        comorbidity_result = self._base_comorbidity.copy()
        trajectory_result = (
            self._base_trajectory.copy()
            if self._base_trajectory is not None else
            pd.DataFrame()
        )

        if comorbidity_result.empty:
            raise ValueError(
                f"{plot_name} cannot be generated because comorbidity_result has no significant positive associations to plot."
            )

        exposure = None
        clustering_edges = self.__combined_clustering_edges(
            comorbidity_result,
            trajectory_result
        )
        clustering_edges = self.__clean_clustering_edges(
            clustering_edges,
            self.comorbidity_beta_col,
        )
        cluster_key = self.__network_cache_key(
            clustering_edges,
            self.comorbidity_beta_col,
        )

        if require_trajectory:
            if trajectory_result.empty:
                raise ValueError(
                    f"{plot_name} cannot be generated because trajectory_result has no significant positive associations to plot."
                )

            if self._exposure_name:
                exposure = 1000
                trajectory_result = self.__sequence(
                    trajectory_result,
                    exposure,
                    self.source_col,
                    self.target_col,
                    self.disease_pair_col
                )

        self.__init_attrs(
            comorbidity = clustering_edges,
            trajectory = trajectory_result,
            exposure = exposure,
            commorbidity_nodes = (
                self.__get_nodes(
                    clustering_edges,
                    self.source_col,
                    self.target_col
                ) if not clustering_edges.empty else set()
            ),
            trajectory_nodes = (
                self.__get_nodes(
                    trajectory_result,
                    self.source_col,
                    self.target_col
                ) if not trajectory_result.empty else set()
            ),
            nodes_attrs = {},
            network_attrs = {}
        )

        self.__make_node_basic_attrs(
            self.phecode_col,
            self.phewas_number_col,
            self.disease_col,
            self.system_col,
        )
        self._prepared_cluster_key = cluster_key
        cached_cluster = self._cluster_cache.get(cluster_key)
        if cached_cluster is not None:
            self.__apply_cluster(cached_cluster)

    def __combined_clustering_edges(
        self,
        comorbidity_result: Df,
        trajectory_result: Df
    ) -> Df:
        """Build the graph used for Louvain modules across all network plots."""
        clustering_edges = comorbidity_result.copy()
        if trajectory_result is None or trajectory_result.empty:
            return clustering_edges

        trajectory_edges = trajectory_result[
            [
                self.source_col,
                self.target_col,
                self.disease_pair_col,
                self.trajectory_beta_col
            ]
        ].copy()
        trajectory_edges = trajectory_edges.rename(
            columns={self.trajectory_beta_col: self.comorbidity_beta_col}
        )

        pair_key_col = "__dinetxify_pair_key"
        clustering_edges[pair_key_col] = clustering_edges.apply(
            lambda row: frozenset(row[[self.source_col, self.target_col]]),
            axis=1
        )
        trajectory_edges[pair_key_col] = trajectory_edges.apply(
            lambda row: frozenset(row[[self.source_col, self.target_col]]),
            axis=1
        )
        trajectory_edges = trajectory_edges[
            ~trajectory_edges[pair_key_col].isin(clustering_edges[pair_key_col])
        ]

        clustering_edges = pd.concat(
            [clustering_edges, trajectory_edges],
            axis=0,
            ignore_index=True
        )
        clustering_edges.drop_duplicates(
            subset=[pair_key_col],
            inplace=True,
            ignore_index=True,
            keep="first"
        )
        clustering_edges.drop(columns=[pair_key_col], inplace=True)
        return clustering_edges

    def __clean_clustering_edges(
        self,
        clustering_edges: Df,
        weight_col: str
    ) -> Df:
        clustering_edges = clustering_edges.copy()
        cleaned_weights = []
        for idx, value in clustering_edges[weight_col].items():
            cleaned_weights.append(
                self.__coerce_edge_weight(value, row_index=idx, weight_col=weight_col)
            )
        clustering_edges[weight_col] = cleaned_weights
        return clustering_edges

    @staticmethod
    def __coerce_node_key(value: Any) -> Any:
        if isinstance(value, np.generic):
            return value.item()
        return value

    @staticmethod
    def __coerce_edge_weight(
        value: Any,
        row_index: Any,
        weight_col: str
    ) -> float:
        if isinstance(value, dict):
            if weight_col in value:
                value = value[weight_col]
            elif "weight" in value:
                value = value["weight"]
            else:
                raise ValueError(
                    f"Invalid edge weight at row {row_index}: expected a number, "
                    f"but got a dictionary without '{weight_col}' or 'weight'."
                )

        try:
            weight = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid edge weight at row {row_index}: expected a number in "
                f"column '{weight_col}', but got {type(value).__name__}."
            ) from exc

        if not np.isfinite(weight):
            raise ValueError(
                f"Invalid edge weight at row {row_index}: column '{weight_col}' "
                "contains a non-finite value."
            )
        return weight

    def __network_cache_key(
        self,
        clustering_edges: Df,
        weight_col: str
    ) -> Tuple[Tuple[Any, Any, float], ...]:
        edge_records = []
        for _, row in clustering_edges[[self._source, self._target, weight_col]].iterrows():
            source = row[self._source]
            target = row[self._target]
            ordered_nodes = tuple(sorted([source, target], key=str))
            edge_records.append(
                (ordered_nodes[0], ordered_nodes[1], float(row[weight_col]))
            )
        return tuple(sorted(edge_records, key=lambda x: (str(x[0]), str(x[1]), x[2])))

    def __apply_cluster(
        self,
        cluster_ans: Dict[Any, int]
    ) -> None:
        self.__update_node_attrs(
            cluster=dict(cluster_ans)
        )
        if cluster_ans:
            self._network_attrs.update(
                {"cluster number": max(cluster_ans.values()) + 1}
            )

    @staticmethod
    def __check_disease_pairs(
        tra_df: Df,
        com_df: Df,
        source: str,
        target: str
    ) -> None:
        """Verifies that all disease pairs in the trajectory network exist in the comorbidity network.

        This validation method checks whether every source-target disease pair present in the
        result of disease trajectory analysis dataframe (tra_df) also exists in the result of 
        comorbidity network analysis dataframe (com_df). Raises a warning if any pair is missing.

        Args:
            tra_df (Df): DataFrame containing temporal disease pairs with source and target columns
            com_df (Df): DataFrame containing non-temporal disease pairs with source and target columns
            source (str): Column name for source diseases in both dataframes
            target (str): Column name for target diseases in both dataframes

    Raises:
        Warning: If any disease pair in tra_df is not found in com_df
        """
        tra_pairs = [
            [row[source], row[target]]
            for _, row in tra_df.iterrows()
        ]

        com_pairs = [
            [row[source], row[target]]
            for _, row in com_df.iterrows()
        ]

        for pair in tra_pairs:
            if pair not in com_pairs:
                warnings.warn(
                    "Disease pairs in trajectory_result are not all present "
                    "in comorbidity_result.",
                    UserWarning,
                    stacklevel=2
                )
                break

    @staticmethod
    def __sequence(
        df: Df,
        exposure: float,
        source: str,
        target: str,
        col_disease_pair: str
    ) -> Df:
        """Generates temporal disease pairs starting from an exposure point.

        Processes a disease network dataframe to create trajectory sequences beginning
        with the specified exposure disease. The method:
        1. Filters diseases connected to the exposure
        2. Creates new trajectory pairs from exposure to each connected disease
        3. Combines them with existing trajectories

        Args:
            df (Df): Input dataframe containing disease trajectory relationships
            exposure (float): The starting disease/exposure point for new trajectories
            source (str): Column name for source diseases (D1 of D1->D2)
            target (str): Column name for target diseases (D2 of D1->D2)
            col_disease_pair (str): Column name to store formatted disease pairs

        Returns:
            Df: Expanded dataframe containing both exposure and a first disease in order 
                of disease trajectories

        Example:
            Given exposure=1.0 and connections to diseases [2.0, 3.0], creates:
            [[1.0, 2.0, "1.0-2.0"], [1.0, 3.0, "1.0-3.0"]] plus original data
        """
        first_layer = df.loc[~df[source].isin(df[target].values)][source].unique()
        d1_d2 = [[exposure, d, f'{exposure}-{d}'] for d in first_layer]
        d1_d2_df = pd.DataFrame(
            d1_d2,
            columns=[source, target, col_disease_pair]
        )
        trajectory_df = pd.concat([df, d1_d2_df])
        return trajectory_df
    
    @staticmethod
    def __get_nodes(
        df: Df,
        source: str,
        target: str
    ) -> set:
        """Extracts all unique disease from source and target columns of a dataframe.

        Combines values from both specified columns and returns them as a set to ensure
        uniqueness. This is particularly useful for network/graph analysis where you need
        to identify all distinct disease.

        Args:
            df (Df): Input dataframe containing node relationships
            source (str): Column name containing source diseases
            target (str): Column name containing target diseases

        Returns:
            set: A set containing all unique nodes from both source and target columns

        Example:
            >>> df = pd.DataFrame({'from': ['A', 'B'], 'to': ['B', 'C']})
            >>> __get_nodes(df, 'from', 'to')
            set('A', 'B', 'C')
        """
        return set(df[source].to_list() + df[target].to_list())
    
    @staticmethod
    def __split_name(name: str) -> str:
        """Formats disease names by inserting line breaks for better readability.
    
        Splits long disease names into multiple lines when the accumulated character count
        exceeds the specified maximum line length (12). Maintains word boundaries and removes
        any trailing whitespace.

        Args:
            name (str): Disease name to be formatted

        Returns:
            str: Formatted disease name with line breaks for better display

        Examples:
            >>> __split_name("Chronic obstructive pulmonary disease")
            'Chronic obstructive\npulmonary disease'
        """
        words = name.split(' ')
        total_number, new_word = 0, ''
        for word in words:
            if total_number >= 12:
                total_number = 0
                total_number += len(word)
                new_word += '\n%s' % (word)
            else:
                total_number += len(word)
                new_word += ' %s' % (word)
        return new_word.strip(' ')

    @staticmethod
    def __sphere_cordinate(
        center: Tuple[float, float, float],
        r: float
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Generates 3D coordinate arrays for a sphere surface.

        Calculates the (x, y, z) coordinates of points uniformly distributed on the surface
        of a sphere with given center and radius using spherical coordinates.

        Args:
            center: (x, y, z) coordinates of the sphere's center point
            radius: Radius of the sphere (must be positive)
            resolution: Number of points along each angular dimension (default: 50)
                    Higher values create smoother spheres but require more memory

        Returns:
            A tuple containing three 2D numpy arrays representing:
            (x_coordinates, y_coordinates, z_coordinates) of the sphere surface points

        Example:
            >>> x, y, z = __sphere_coordinate((0, 0, 0), 1.0)
            >>> x.shape
            (50, 50)
        """
        theta1 = np.linspace(0, 2*np.pi, 50)
        phi1 = np.linspace(0, np.pi, 50)
        x = r * np.outer(np.sin(theta1), np.sin(phi1))
        y = r * np.outer(np.cos(theta1), np.sin(phi1))
        z = r * np.outer(np.ones(50), np.cos(phi1))
        x += center[0]
        y += center[1]
        z += center[2]
        return (x, y, z)
    
    @staticmethod
    def __calculate_order(
        source_lst :List[float], 
        target_lst :List[float],
        exposure
    ) -> Dict[float, int]:
        """Calculates topological order distances for nodes in a directed graph.

        Performs a topological sort using Kahn's algorithm to determine the longest path
        distance from any source node (in-degree 0) to each node in the graph. This is
        useful for determining hierarchical relationships in dependency graphs.

        Args:
            source_lst: List of source nodes for each edge in the graph
            target_lst: List of target nodes for each edge in the graph
                (must be same length as source_lst)

        Returns:
            A dictionary mapping each node to its maximum distance from a source node.
            Source nodes have distance 1 by default.

        Example:
            >>> source = [1.0, 1.0, 2.0, 3.0]
            >>> target = [2.0, 3.0, 4.0, 4.0]
            >>> __calculate_order(source, target)
            {1.0: 1, 2.0: 2, 3.0: 2, 4.0: 3}
            
        Note:
            - The graph must be a Directed Acyclic Graph (DAG)
            - If cycles exist, the algorithm will only process nodes reachable
            from true source nodes (in-degree 0)
        """
        tra_df = pd.DataFrame({'d1':source_lst,'d2':target_lst})

        if exposure == 0:
            first_layer = tra_df.loc[~tra_df['d1'].isin(tra_df['d2'].values)]['d1'].unique()
            d1_d2 = [[exposure, d] for d in first_layer]
            d1_d2_df = pd.DataFrame(
                d1_d2,
                columns=['d1', 'd2']
            )
            tra_df = pd.concat([tra_df, d1_d2_df])

        G = nx.from_pandas_edgelist(tra_df, source='d1', target='d2', create_using=nx.DiGraph())
        # Compute strongly connected components (SCCs)
        scc = list(nx.strongly_connected_components(G))
        node_to_scc = {}
        for idx, comp in enumerate(scc):
            for node in comp:
                node_to_scc[node] = idx
        # Condense the graph into a DAG where each node is a SCC
        C = nx.condensation(G, scc)
        # Use BFS on the condensed graph to assign layers
        #get the exposure from df
        exposure = tra_df[tra_df['d1']==exposure]['d1'].values[0]
        start_comp = [node_to_scc[exposure]]
        comp_layer = {node: 1 for node in start_comp}
        queue = start_comp
        while queue:
            current = queue.pop(0)
            current_layer = comp_layer[current]
            for neighbor in C.successors(current):
                # Update the layer if this gives a longer path from the start
                if neighbor not in comp_layer or comp_layer[neighbor] < current_layer + 1:
                    comp_layer[neighbor] = current_layer + 1
                    queue.append(neighbor)
        # Map each original node to its component's layer
        d_lst_layer = {node: comp_layer[node_to_scc[node]] for node in G.nodes() if node_to_scc[node] in comp_layer}

        return d_lst_layer
    
    @staticmethod
    def __most_frequent_element(lst: List[Any]) -> Any:
        """Finds the most frequently occurring element in a list.

        Args:
            lst: Input list of elements to analyze

        Returns:
            The most frequent element, or tuple of (element, count) if return_count=True.
            Returns None if input list is empty.

        Raises:
            ValueError: If input list is empty and no default tie_breaker behavior is specified

        Examples:
            >>> __most_frequent_element([1, 2, 2, 3, 3, 3])
            3
        """
        if not lst:
            return None
        
        counter = Counter(lst)
        most_common = counter.most_common(1)
        return most_common[0][0]

    def __check_node_attrs(self, key: str) -> bool:
        """Checks if a given attribute key exists in any node's attributes.

        This method verifies whether the specified attribute key is present in any
        of the node attributes stored in the graph. The check stops at the first
        occurrence of the key.

        Args:
            key: The attribute key to search for in node attributes

        Returns:
            True if the key exists in any node's attributes, False otherwise

        Example:
            >>> graph._nodes_attrs = {1: {'color': 'red'}, 2: {'size': 10}}
            >>> graph.__check_node_attrs('color')
            True
            >>> graph.__check_node_attrs('weight')
            False

        Note:
            This performs a shallow check for key existence only, not the values.
        """
        for _, attrs in self._nodes_attrs.items():
            if key in attrs.keys():
                continue
            else:
                return False
        return True

    def __calculate_ratio(
        self, 
        cluster_nodes: Dict[int, float]
    ) -> Dict[int, float]:
        """Calculates size ratios for clusters relative to the total network size.

        For each cluster, sums the sizes of all nodes in the cluster, then computes
        each cluster's proportion of the total network size. This is useful for
        visualizing cluster importance or dominance in the network.

        Args:
            cluster_nodes: Dictionary mapping cluster IDs to lists of node IDs.
                        Format: {cluster_id: [node_id1, node_id2, ...]}

        Returns:
            Dictionary mapping cluster IDs to their proportional size ratios.
            Ratios sum to 1.0 across all clusters.

        Example:
            >>> self._nodes_attrs = {
                    1: {"size": 10}, 
                    2: {"size": 20},
                    3: {"size": 30}
                }
            >>> cluster_nodes = {0: [1, 2], 1: [3]}
            >>> self.__calculate_ratio(cluster_nodes)
            {0: 0.5, 1: 0.5}  # (10+20)/60 and 30/60

        Raises:
            KeyError: If any node in cluster_nodes is missing from self._nodes_attrs
            ZeroDivisionError: If total network size is 0
        """
        for cluster, nodes in cluster_nodes.items():
            size = [self._nodes_attrs[x]["size"] for x in nodes]
            sum_size = sum(size)
            cluster_nodes.update({cluster:sum_size})

        total_size = sum(cluster_nodes.values())
        for key, value in cluster_nodes.items():
            cluster_nodes.update({key:value/total_size})
        return cluster_nodes
    
    def __init_attrs(self, **kwargs) -> None:
        """Initializes protected attributes from keyword arguments.
    
        Dynamically sets protected instance attributes (prefixed with '_') from provided
        keyword arguments. This provides a flexible way to initialize multiple protected
        attributes at once.

        Args:
            **kwargs: Arbitrary keyword arguments where:
                    key = attribute name (will be prefixed with '_')
                    value = attribute value

        Example:
            >>> obj.__init_attrs(size=10, color='red', visible=True)
            >>> obj._size
            10
            >>> obj._color
            'red'
            >>> obj._visible
            True

        Note:
            - All created attributes will be protected (start with underscore)
            - Existing attributes with same name will be overwritten
            - For public attributes, consider using regular initialization
        """
        for key, value in kwargs.items():
            setattr(self, f"_{key}", value)
    
    def __update_node_attrs(self, **kwargs) -> None:
        """Update stored node-attribute dictionaries from keyed mappings."""
        for key, value in kwargs.items():
            for node, attr in value.items():
                if node in self._nodes_attrs.keys():
                    self._nodes_attrs[node].update({f"{key}":attr})

    def __filter_significant(
        self,
        phewas_result: Df,
        comorbidity_result: Optional[Df],
        trajectory_result: Optional[Df],
    ) -> Tuple[Df, Optional[Df], Optional[Df]]:
        """Filters input dataframes to only include statistically significant results.

        Applies boolean filters to each input dataframe based on specified significance
        columns, returning only rows marked as significant (True) and with positive effect sizes 
        for comorbidity network and disease trajectory reesults.

        Args:
            phewas_result: DataFrame containing PheWAS analysis results
            comorbidity_result: DataFrame containing comorbidity analysis results
            trajectory_result: DataFrame containing disease trajectory results

        Returns:
            A tuple containing the filtered dataframes in order:
            (filtered_phewas, filtered_comorbidity, filtered_trajectory)

        Example:
            >>> phewas_df = pd.DataFrame({'p_value': [0.01, 0.5], 'significant': [True, False]})
            >>> comorbidity_df = pd.DataFrame({'OR': [1.2, 3.4], 'sig': [True, True]})
            >>> filtered = __filter_significant(
                    phewas_df, comorbidity_df, None,
                    filter_phewas_col='significant',
                    filter_comorbidity_col='sig'
                )
            # Returns:
            # - phewas_df with only significant rows
            # - comorbidity_df unchanged (all rows significant)
            # - None (no trajectory input)

        Note:
            - None values for filter columns skip filtering for that dataframe
            - Original dataframes are not modified (returns filtered copies
        """
        phewas_result = phewas_result.copy()

        if comorbidity_result is not None:
            comorbidity_result = comorbidity_result.loc[
                (comorbidity_result[self.comorbidity_significance_col] == True) &
                (comorbidity_result[self.comorbidity_beta_col] > 0)
            ].copy()

        if trajectory_result is not None:
            trajectory_result = trajectory_result.loc[
                (trajectory_result[self.trajectory_significance_col] == True) &
                (trajectory_result[self.trajectory_beta_col] > 0)
            ].copy()
        
        return phewas_result, comorbidity_result, trajectory_result
    
    def __get_same_nodes(self, key_name:str) -> Dict[Any, Any]:
        """Groups node IDs by their attribute values for a specified attribute key.

        Creates a dictionary mapping each unique attribute value to a list of node IDs
        that share that value. This is useful for analyzing or visualizing nodes with
        common characteristics.

        Args:
            key_name: The name of the attribute to group by. Must exist in all nodes'
                    attribute dictionaries.

        Returns:
            A dictionary where:
            - Keys are the unique attribute values found for the specified key
            - Values are lists of node IDs that have each attribute value

        Example:
            >>> graph._nodes_attrs = {
                    1: {'color': 'red', 'size': 10},
                    2: {'color': 'blue', 'size': 20},
                    3: {'color': 'red', 'size': 15}
                }
            >>> graph.group_nodes_by_attribute('color')
            {'red': [1, 3], 'blue': [2]}

        Raises:
            KeyError: If the specified key_name is not found in any node's attributes
        """
        attrs = []
        for attr in self._nodes_attrs.values():
            attrs.append(attr[key_name])
        
        nodes = {x:[] for x in set(attrs)}
        for attr in set(attrs):
            for node, value in self._nodes_attrs.items():
                if value[key_name] == attr:
                    nodes[attr].append(node)
        return nodes

    def __get_edge_attrs(
        self,
        edge_lst: List[tuple[str, float]],
    ) -> Dict[str, List[float]]:
        """Generates edge attribute data for 3D visualization coordinates.

        Processes a list of edges to create coordinate sequences for visualization,
        handling special cases for exposure nodes. Produces three separate coordinate
        streams (x, y, z) with None values separating different edges.

        Args:
            edge_lst: List of edges as tuples (source_node, target_node)
                    where nodes are referenced by their string IDs

        Returns:
            Dictionary with three coordinate streams:
            {
                0: List of x-coordinates [source_x, target_x, None, ...],
                1: List of y-coordinates [source_y, target_y, None, ...],
                2: List of z-coordinates [source_z, target_z, None, ...]
            }
            None values separate different edges in the visualization path.

        Example:
            >>> self._exposure = "E1"
            >>> self._exposure_location = (1.0, 2.0, 3.0)
            >>> self._nodes_attrs = {
                    "N1": {"location": (4.0, 5.0, 6.0)},
                    "N2": {"location": (7.0, 8.0, 9.0)}
                }
            >>> self.__get_edge_attrs([("E1", "N1"), ("N1", "N2")])
            {
                0: [1.0, 4.0, None, 4.0, 7.0, None],
                1: [2.0, 5.0, None, 5.0, 8.0, None],
                2: [3.0, 6.0, None, 6.0, 9.0, None]
            }

        Note:
            - Exposure nodes use special predefined locations
            - None values are inserted between edges for visualization breaks
            - Assumes all nodes have "location" attributes with 3D coordinates
        """
        edge_attrs = {
            0:[],
            1:[],
            2:[]
        }
        for edge in edge_lst:
            source = edge[0]
            target = edge[1]
            if source == self._exposure:
                source_loc = self._exposure_location
            else:
                source_loc = self._nodes_attrs[source]["location"]
            target_loc = self._nodes_attrs[target]["location"]

            for i in range(3):
                edge_attrs[i] += [
                    source_loc[i], 
                    target_loc[i], 
                    None
                ]
        return edge_attrs

    def __terminal_trajectory_nodes(self) -> List[Any]:
        """Return trajectory sink nodes for terminal-outcome visualization."""
        if self._trajectory.empty:
            return []
        source_nodes = set(self._trajectory[self._source].to_list())
        target_nodes = set(self._trajectory[self._target].to_list())
        terminal_nodes = target_nodes - source_nodes
        if len(terminal_nodes) == 0:
            terminal_nodes = target_nodes or source_nodes
        if self._exposure is not None:
            terminal_nodes.discard(self._exposure)
        return [
            node for node in sorted(terminal_nodes, key=str)
            if node in self._nodes_attrs
        ]

    def __outcome_3d_location(self, layer_distance: float) -> Tuple[float, float, float]:
        """Return explicit or automatic 3D location for the terminal outcome node."""
        if self._outcome_location is not None:
            return self._outcome_location
        if not self._nodes_attrs:
            return (0, 0, 0)
        min_z = min(attrs["location"][2] for attrs in self._nodes_attrs.values())
        return (
            self._exposure_location[0],
            self._exposure_location[1],
            min_z - layer_distance,
        )

    def __sig_nodes(self) -> List[List[float]]:
        """Identifies significant disease progression paths in the trajectory network.

        Analyzes the disease trajectory network to find:
        1. All longest paths starting from exposure node
        2. Significant cycles (if no exposure specified)
        3. Paths where all nodes belong to the same cluster

        Returns:
            List of significant paths, where each path is represented as:
            List[node_ids] with nodes in progression order

        Algorithm:
            1. Builds directed graph from trajectory data
            2. Finds all longest paths from exposure node (or node 0 if no exposure)
            3. Identifies significant cycles (when no exposure specified)
            4. Filters paths where all nodes share the same cluster

        Example:
            >>> self._exposure = "D1"
            >>> self._trajectory = pd.DataFrame({
                    'source': ['D1', 'D1', 'D2'],
                    'target': ['D2', 'D3', 'D4']
                })
            >>> self.__sig_nodes()
            [['D2'], ['D3'], ['D2', 'D4']]  # Sample output

        Note:
            - Uses DFS to find longest paths
            - Only considers cycles when no exposure is specified
            - Path significance determined by cluster homogeneity
        """
        def get_all_longest_paths(G, start):
            all_paths = []
            def dfs(current_node, path):
                if current_node in path:
                    return
                path.append(current_node)
                if not G.out_edges(current_node):
                    all_paths.append(path.copy())
                else:
                    for neighbor in G.successors(current_node):
                        dfs(neighbor, path.copy())
            
            dfs(start, [])
            return all_paths

        directed_graph = nx.DiGraph()
        tra_df = self._trajectory
        
        if self._exposure is None:
            exposure = 0
            tra_df = self.__sequence(
                self._trajectory,
                exposure,
                self._source,
                self._target,
                "name_disease_pair"
            )
        else:
            exposure = self._exposure

        pairs = tra_df[
            [
                self._source,
                self._target
            ]
        ].values

        for source, target in pairs:
            directed_graph.add_edge(
                source,
                target
            )

        all_paths = get_all_longest_paths(
            directed_graph,
            exposure
        )

        all_paths = [path[1::] for path in all_paths]

        if exposure==0:
            cycles = nx.simple_cycles(directed_graph)
            for cycle in cycles:
                if self._nodes_attrs[cycle[0]]["order"]==1:
                    all_paths.append(cycle)
        
        sig_paths = []
        for path in all_paths:
            cluster = [
                self._nodes_attrs[node]["cluster"] for node in path
            ]
            if len(set(cluster))==1:
                sig_paths.append(path)
        return sig_paths

    def __sphere_attrs(
        self, 
        center: Tuple[float], 
        r: Tuple[float], 
        color: str, 
        light_dict: Optional[Dict[str, float]]=dict(
            ambient=0.2,
            diffuse=0.8,
            specular=0.4,
            roughness=0.2,
            fresnel=2.0
        ),
        light_position_dict: Optional[Dict[str, float]]=dict(
            x=1.5,
            y=1.5,
            z=1.5
        )
    ):
        """Generates visualization attributes for 3D sphere representation.

        Creates all necessary attributes for plotting a 3D sphere in Plotly, including
        coordinates, colors, and lighting parameters. Used for visualizing disease nodes
        in 3D network visualizations.

        Args:
            center: (x, y, z) coordinates of sphere center
            radius: Radius of the sphere (positive float)
            color: Hex or named color for sphere visualization
            name: Name identifier for the disease/sphere
            label: Display label for the sphere
            light_dict: Lighting properties for Plotly surface rendering. Defaults to:
                    {
                        'ambient': 0.2,
                        'diffuse': 0.8,
                        'specular': 0.4,
                        'roughness': 0.2,
                        'fresnel': 2.0
                    }
            light_position_dict: 3D position of light source. Defaults to:
                                {'x': 1.5, 'y': 1.5, 'z': 1.5}

        Returns:
            A tuple containing:
            - x, y, z: 2D numpy arrays of sphere surface coordinates
            - colorscale: Color gradient definition for Plotly
            - light_dict: Final lighting properties used
            - light_position_dict: Final light position used

        Example:
            >>> x, y, z, colors, lights, light_pos = __sphere_attrs(
                    center=(0, 0, 0),
                    radius=1.5,
                    color='#FF0000',
                    name='Diabetes',
                    label='DM2'
                )
            >>> x.shape  # Returns coordinate arrays
            (50, 50)
        """
        x, y, z = self.__sphere_cordinate(center, r)
        colorscale = [[0.0, color], [0.5, color], [1.0, color]]
        return x, y, z, colorscale, light_dict, light_position_dict

    def __calculate_location_random(
            self, 
            max_radius :float,
            min_radius :float,
            cluster_reduction_ratio :float,
            z_axis: float,
            cluster_ratio :Dict[int, float],
            max_attempts :Optional[int]=200
        ) -> None:
        """Randomly calculates and assigns 3D locations for clustered nodes.

        Distributes nodes in 3D space with the following characteristics:
        - Nodes are positioned in circular sectors based on their cluster
        - Radial distance is randomized between min/max radius
        - Z-coordinate decreases with node order (hierarchy)
        - A bounded collision check avoids obvious overlap without the long
          retry loops used by earlier versions

        Args:
            max_radius: Maximum radial distance from center
            min_radius: Minimum radial distance from center
            cluster_reduction_ratio: Ratio to reduce cluster sector angles (prevents edge crowding)
            z_axis: Z-axis spacing between different orders
            cluster_ratio: Dictionary mapping clusters to their angular proportion
                        (e.g., {0: 0.3, 1: 0.7} for 30%/70% split)
            max_attempts: Maximum random placements tried per node. Defaults
                to 200 to keep layout generation bounded.

        Returns:
            None: Updates node locations directly in self._nodes_attrs

        Algorithm:
            1. Assigns exposure location as (0, 0, 0) if not set
            2. Calculates angular range based on cluster
            3. Randomly samples radius and angle within the cluster sector
            4. Uses a coarse spatial index for fast overlap checks
            5. Updates node locations

        Example:
            >>> self._calculate_location_random(
                    max_radius=10,
                    min_radius=5,
                    cluster_reduction_ratio=0.1,
                    z_axis=0.5,
                    cluster_ratio={0: 0.4, 1: 0.6}
                )
            # Updates self._nodes_attrs with "location" entries like:
            # {node1: {"location": (3.2, 4.1, -1.0)}, ...}
        """
        if self._exposure_location is None:
            self._exposure_location = (0, 0, 0)

        cluster_number = self._network_attrs["cluster number"]
        interval_ratio = min(cluster_number * 0.05 * cluster_reduction_ratio, 0.8)
        gap_angle = interval_ratio / cluster_number * 2 * math.pi
        cluster_location = {x:{} for x in range(cluster_number)}
        cluster_bins = {x:{} for x in range(cluster_number)}
        attempts = max(1, int(max_attempts))
        max_node_size = max(
            [attrs["size"] for attrs in self._nodes_attrs.values()] or [1.0]
        )
        bin_size = max(max_node_size * 2, 1.0)

        def _bin_key(loc: Tuple[float, float, float]) -> Tuple[int, int]:
            return (
                math.floor(loc[0] / bin_size),
                math.floor(loc[1] / bin_size)
            )

        def _nearby_nodes(cluster: int, loc: Tuple[float, float, float]):
            bin_x, bin_y = _bin_key(loc)
            for x_offset in (-1, 0, 1):
                for y_offset in (-1, 0, 1):
                    for oth_node in cluster_bins[cluster].get(
                        (bin_x + x_offset, bin_y + y_offset),
                        []
                    ):
                        yield oth_node, cluster_location[cluster][oth_node]

        for node, attrs in self._nodes_attrs.items():
            order = attrs["order"]
            cluster = attrs["cluster"]
            previous_ratio = sum(
                cluster_ratio.get(i, 0)
                for i in range(cluster)
            )
            min_ang = 2 * math.pi * previous_ratio * (1 - interval_ratio) + cluster * gap_angle
            max_ang = min_ang + 2 * math.pi * cluster_ratio.get(cluster, 0) * (1 - interval_ratio)
            if max_ang <= min_ang:
                max_ang = min_ang + 1e-8

            node_loc = None
            for _ in range(attempts):
                radius = random.uniform(min_radius, max_radius)
                node_ang = random.uniform(min_ang, max_ang)
                candidate_loc = (
                    radius * math.cos(node_ang),
                    radius * math.sin(node_ang),
                    self._exposure_location[2] - order * z_axis
                )
                node_loc = candidate_loc
                is_separate = True
                for oth_node, loc in _nearby_nodes(cluster, candidate_loc):
                    distance = math.sqrt(
                        (loc[0] - candidate_loc[0])**2 +
                        (loc[1] - candidate_loc[1])**2
                    )
                    min_distance = (
                        self._nodes_attrs[node]["size"] +
                        self._nodes_attrs[oth_node]["size"]
                    )
                    if distance < min_distance:
                        is_separate = False
                        break

                if is_separate:
                    break

            cluster_location[cluster].update({node: node_loc})
            cluster_bins[cluster].setdefault(_bin_key(node_loc), []).append(node)
            self._nodes_attrs[node].update({"location": node_loc})
    
    def __make_node_basic_attrs(
        self,
        phewas_phecode: str,
        phewas_number: str,
        disease_col: str,
        system_col: str
    ) -> None:
        """Initializes and updates basic node attributes for network visualization.

        Creates and assigns three core node attributes:
        1. Node names (formatted with phenotype and phecode)
        2. Disease system/category
        3. Node sizes (calculated from phewas results)
        4. System-based colors

        Args:
            phewas_phecode: Column name in self._phewas containing phecode identifiers
            phewas_number: Column name in self._phewas containing numerical values 
                        used for size calculation

        Returns:
            None: Updates self._nodes_attrs and self._network_attrs in place

        Processing Steps:
            1. Creates formatted node names: "Phenotype (Phecode.X)"
            2. Maps each node to its disease system/category
            3. Calculates node sizes from phewas values (cube root scaling)
            4. Assigns system-specific colors to nodes
            5. Updates network-wide system color mapping

        Example:
            >>> self._make_node_basic_attrs(
                    phewas_phecode="phecode",
                    phewas_number="p_value"
                )
            # Results in nodes_attrs containing:
            # {
            #   123.1: {
            #     "name": "Diabetes (123.1)",
            #     "system": "Endocrine",
            #     "size": 4.2,
            #     "color": "#FF0000"
            #   },
            #   ...
            # }
        """
        # disease name attrs
        node_name = {
            node:self.__split_name('%s (%.1f)' % (name, node)) 
            for node, name in self._describe[
                [phewas_phecode, disease_col]
            ].values
        }

        # disease system attrs
        node_system = dict(
            zip(
                self._describe[phewas_phecode], 
                self._describe[system_col]
            )
        )

        # disease size attrs
        node_size = dict(
            zip(
                self._phewas[phewas_phecode], 
                np.cbrt(
                    3*self._phewas[phewas_number]/(4*np.pi)
                )
            )
        )

        for node in self._commorbidity_nodes:
            if node_system[node] in self._system_color.keys():
                self._nodes_attrs.update({node:{}})

        self.__update_node_attrs(
            name = node_name,
            system = node_system,
            size = node_size
        )

        self._network_attrs.update({"system color":self._system_color})
        for _, attrs in self._nodes_attrs.items():
            attrs.update({"color":self._system_color[attrs["system"]]})

    def __cluster(
        self,
        weight: str,
        max_attempts: Optional[int]=100,
    ) -> None:
        """Performs Louvain community detection on the prepared module graph.

        Identifies disease modules from the prepared clustering graph using a
        simple direct NetworkX workflow. Louvain is run repeatedly with seeds
        ``0`` through ``max_attempts - 1`` and the partition with the highest
        NetworkX modularity is retained. No external community package,
        multiprocessing, or subprocess worker is used.

        Args:
            weight: Edge weight column name from the prepared clustering data.
            max_attempts: Number of random seeds evaluated by NetworkX Louvain.
                The partition with highest modularity is retained.
                Defaults to 100.

        Returns:
            None: Updates self._nodes_attrs and self._network_attrs in place with:
                - Individual node cluster assignments
                - Total number of clusters found

        Algorithm Steps:
            1. Constructs an undirected NetworkX graph from prepared disease pairs
            2. Runs NetworkX Louvain with different random seeds
            3. Selects the partition with highest NetworkX modularity
            4. Stores cluster assignments in node attributes
            5. Records total cluster count and modularity in network attributes

        Example:
            >>> self.__cluster(weight='comorbidity_beta')
            # Updates:
            # - self._nodes_attrs with 'cluster' assignments
            # - self._network_attrs['cluster number'] with cluster count
        """
        if max_attempts < 1:
            raise ValueError("max_attempts must be at least 1")

        cluster_key = getattr(self, "_prepared_cluster_key", None)
        cached_cluster = getattr(self, "_cluster_cache", {}).get(cluster_key)
        if cached_cluster is not None:
            self.__apply_cluster(cached_cluster)
            return

        Graph_position = nx.Graph()
        for idx, row in self._comorbidity.iterrows():
            source = self.__coerce_node_key(row[self._source])
            target = self.__coerce_node_key(row[self._target])
            edge_weight = self.__coerce_edge_weight(
                row[weight],
                row_index=idx,
                weight_col=weight,
            )
            Graph_position.add_edge(
                source,
                target,
                weight=edge_weight
            )

        best_modularity = -np.inf
        best_communities = None
        for i in range(max_attempts):
            communities = nx_community.louvain_communities(
                Graph_position,
                weight="weight",
                seed=i
            )
            modularity = nx_community.modularity(
                Graph_position,
                communities,
                weight="weight"
            )
            if modularity >= best_modularity:
                best_modularity = modularity
                best_communities = communities

        cluster_ans = {}
        for cluster_id, community_nodes in enumerate(best_communities):
            for node in community_nodes:
                cluster_ans[node] = cluster_id
        self.__apply_cluster(cluster_ans)
        self._network_attrs.update({"modularity": best_modularity})
        if cluster_key is not None:
            self._cluster_cache[cluster_key] = cluster_ans

    def __ensure_trajectory_nodes_clustered(self, plot_name: str) -> None:
        if self._trajectory.empty:
            return

        trajectory_nodes = self.__get_nodes(
            self._trajectory,
            self._source,
            self._target
        )
        if self._exposure is not None:
            trajectory_nodes.discard(self._exposure)

        missing_nodes = [
            node for node in trajectory_nodes
            if node not in self._nodes_attrs or "cluster" not in self._nodes_attrs[node]
        ]
        if missing_nodes:
            missing_text = ", ".join(str(node) for node in sorted(missing_nodes, key=str))
            raise ValueError(
                f"{plot_name} cannot assign modules to trajectory node(s): {missing_text}. "
                "Some nodes in trajectory_result do not appear in the comorbidity network "
                "used for clustering. Create the Plot object with both comorbidity_result "
                "and trajectory_result, then rerun comorbidity_network_plot() and the "
                "trajectory/3D plot."
            )

    def __make_location_random(
        self, 
        max_radius: float, 
        min_radius: float,
        distance: float,
        cluster_reduction_ratio: float
    ) -> None:
        """Distributes nodes in 3D space with random cluster-sector placement.

        Assigns 3D coordinates to nodes with the following spatial organization:
        - Nodes from the same cluster are grouped in circular sectors in x-y plane
        - Nodes are randomly distributed between min_radius and max_radius
        - Z-coordinate represents hierarchy/layer (lower values for higher orders)
        - Cluster sectors are proportionally sized based on total node size
        - Includes buffer spacing between clusters via reduction ratio
        - Uses a bounded random placement loop with a coarse spatial index

        Args:
            max_radius: Maximum distance from origin in x-y plane (must be > min_radius)
            min_radius: Minimum distance from origin in x-y plane (must be >= 0)
            distance: Vertical spacing between layers along z-axis
            cluster_reduction_ratio: Ratio to reduce cluster sector angles (0-1)
                                Higher values create more space between clusters

        Returns:
            None: Updates node locations directly in self._nodes_attrs

        Spatial Organization:
            - x-y plane: Nodes distributed in radial sectors by cluster
            - z-axis: Represents hierarchical order (exposure on top)
            - Cluster sectors sized proportionally to total node size

        Example:
            >>> self.__make_location_random(
                    max_radius=10.0,
                    min_radius=2.0,
                    distance=0.5,
                    cluster_reduction_ratio=0.1
                )
            # Updates node locations in self._nodes_attrs like:
            # {
            #   "549.4": (3.2, 4.1, -1.5),
            #   "250.1": (5.7, 2.3, -2.0),
            #   ...
            # }
        """
        same_cluster_nodes = self.__get_same_nodes("cluster")
        cluster_ratio = self.__calculate_ratio(same_cluster_nodes)
        self.__calculate_location_random(
            max_radius,
            min_radius,
            cluster_reduction_ratio,
            distance,
            cluster_ratio
        )
    
    def __trajectory_order(self) -> None:
        """Determines hierarchical ordering of nodes in the disease trajectory network.

        Calculates and assigns order/level numbers to nodes based on their position in the
        disease progression paths (D1 → D2 relationships). The ordering follows these rules:
        - Exposure node is considered order 0 (if present)
        - Direct targets of exposure are order 1
        - Each subsequent step in progression increments the order
        - Maximum order is stored as a network attribute

        Returns:
            None: Updates node attributes with 'order' values and network attributes
                with maximum order number

        Algorithm:
            1. Filters trajectory relationships (excluding exposure as source)
            2. Calculates node orders using topological sorting
            3. Updates node attributes with order numbers
            4. Records maximum order in network attributes

        Example:
            Given trajectory D1→D2→D3 and exposure D1:
            - D1: order 0 (exposure)
            - D2: order 1
            - D3: order 2
            - Network attribute "order number" = 2
        """
        tra_df = self._trajectory
        
        if self._exposure:
            node_order = self.__calculate_order(
                tra_df[self._source].to_list(),
                tra_df[self._target].to_list(),
                self._exposure
            )
        else:
            node_order = self.__calculate_order(
                tra_df[self._source].to_list(),
                tra_df[self._target].to_list(),
                0
            )

        self.__update_node_attrs(
            order = node_order
        )

        self._network_attrs.update(
            {"order number":max(node_order.values())}
        )

    def __comorbidity_order(self) -> None:
        """Determines hierarchical ordering for nodes in the comorbidity network.

        Assigns order numbers to nodes based on their comorbidity relationships using:
        1. Direct comorbidity partners' orders (most frequent)
        2. Same-cluster nodes' orders (if no direct partners)
        3. Network maximum order (as fallback)

        The ordering follows these rules:
        - For nodes with comorbidity partners: use most frequent partner order
        - For isolated nodes: use most frequent order from same cluster
        - For completely unconnected nodes: use network maximum order

        Returns:
            None: Updates node attributes with 'order' values

        Processing Steps:
            1. First pass: Assign orders based on direct comorbidity relationships
            2. Second pass: Assign orders based on cluster membership
            3. Final pass: Assign network maximum order to remaining nodes

        Example:
            Given comorbidity pairs:
            - D1 (order 1) - D2
            - D1 - D3 (order 2)
            - D4 (no partners) in cluster with D1
            Then:
            - D2 gets order 1 (from D1)
            - D3 keeps order 2
            - D4 gets order 1 (from cluster)
        """
        comorbidity_nodes = self._comorbidity[[self._source, self._target]].values
        for node, attr in self._nodes_attrs.items():
            if "order" not in attr.keys():
                pairs = [pair.tolist() for pair in comorbidity_nodes if node in pair]
                nodes = [x for pair in pairs for x in pair if x!=node]
                orders = [self._nodes_attrs[x].get("order", None) for x in nodes]
                orders = list(filter(None, orders))
                if orders:
                    order = self.__most_frequent_element(orders)
                    self._nodes_attrs[node].update({"order":order})

        same_cluster_nodes = self.__get_same_nodes("cluster")
        for node, attr in self._nodes_attrs.items():
            if "order" not in attr.keys():
                nodes = same_cluster_nodes[attr["cluster"]]
                orders = [self._nodes_attrs[x].get("order", None) for x in nodes]
                orders = list(filter(None, orders))
                if orders:
                    order = self.__most_frequent_element(orders)
                    self._nodes_attrs[node].update({"order":order})

        for node, attr in self._nodes_attrs.items():
            if "order" not in attr.keys():
                order = self._network_attrs["order number"]
                self._nodes_attrs[node].update({"order":order})

    def __plot(
        self, 
        line_width: float,
        line_color: str, 
        size_reduction: float,
    ) -> List[Any]:
        """Generates 3D plot data for disease trajectories and comorbidities network visualization.

        Creates Plotly graphical objects representing:
        - Disease nodes as colored spheres (grouped by disease system)
        - Trajectory edges as connecting lines between nodes

        Args:
            line_width: Width of trajectory lines (default: 1.0)
            line_color: Color of trajectory lines (default: 'gray')
            size_reduction: Scaling factor for node sizes (default: 1.0)

        Returns:
            List of Plotly graphical objects (Surface and Scatter3d) containing:
            - Spheres for each disease node (phecode), colored by system
            - Lines for all disease trajectories between nodes

        Example:
            >>> plot_data = network.__plot(line_width=2.0, line_color='blue')
            >>> fig = go.Figure(data=plot_data)
            >>> fig.show()

        Note:
            - Nodes are grouped by disease system in the legend
            - First node in each system shows in legend (others hidden)
            - Hover shows disease name
            - Uses spherical coordinates for 3D layout
        """
        plot_data = [] 
        sys_nodes = self.__get_same_nodes("system")
        for sys, nodes in sys_nodes.items():
            is_showlegend = True
            for node in nodes:
                plot_attrs = self.__sphere_attrs(
                    self._nodes_attrs[node]["location"],
                    self._nodes_attrs[node]["size"]*size_reduction,
                    self._nodes_attrs[node]["color"]
                )

                if node != nodes[0]:
                    is_showlegend = False
                
                data = go.Surface(
                    x=plot_attrs[0],
                    y=plot_attrs[1],
                    z=plot_attrs[2],
                    colorscale=plot_attrs[3],
                    showlegend=is_showlegend,
                    lighting=plot_attrs[4],
                    hovertemplate=self._nodes_attrs[node]["name"],
                    name="%s" % (self.sys_dict[sys]),
                    showscale=False,
                    legendgroup="sphere",
                    legendgrouptitle_text="Disease Systems",
                    lightposition=plot_attrs[-1]
                )
                plot_data.append(data)

        tra_edges = zip(
            self._trajectory[self._source],
            self._trajectory[self._target]
        )


        all_edges = list(tra_edges)
        edges_attrs = self.__get_edge_attrs(all_edges)

        trace_data = go.Scatter3d(
            x=edges_attrs[0],
            y=edges_attrs[1],
            z=edges_attrs[2],
            line=dict(
                color=line_color,
                width=line_width
            ),
            mode='lines',
            legendgrouptitle_text='Types of connections',
            name='Trajectories',
            showlegend=True,
            hoverinfo=None
        )
        plot_data.append(trace_data)
        return plot_data
    
    def get_louvain_clusters(
        self,
        max_attempts: Optional[int]=100,
    ) -> Dict[Any, int]:
        """Runs Louvain clustering and returns disease-module labels.

        This public method computes community labels using the internal ``__cluster``
        workflow and the configured comorbidity edge-weight column. If both
        comorbidity and trajectory results are present, the prepared module graph
        includes comorbidity pairs plus trajectory-only pairs. If cluster labels are
        already available in node attributes, existing labels are reused.

        Args:
            max_attempts: Number of random seeds evaluated by NetworkX Louvain.
                The partition with highest modularity is retained.
                Defaults to 100.

        Returns:
            Dictionary mapping each node ID to its assigned cluster ID.

        Example:
            >>> cluster_dict = result_plot.get_louvain_clusters()
        """
        self.__prepare_network_data("get_louvain_clusters")
        if not self.__check_node_attrs("cluster"):
            self.__cluster(
                self.comorbidity_beta_col,
                max_attempts=max_attempts,
            )
        self.__ensure_trajectory_nodes_clustered("get_louvain_clusters")
        return {node: attrs["cluster"] for node, attrs in self._nodes_attrs.items()}

    def interactive_website(self, path: str) -> None:
        """Generate an offline interactive website for available analysis results.

        The generated directory contains a landing page and full PheWAS scatter
        and forest explorers with system filters, disease selection, tooltips,
        zoom, and image export. ``phewas_result`` is always visualized. A
        comorbidity-network section is added when
        ``comorbidity_result`` was supplied to :class:`Plot`; disease-trajectory
        and 3D sections are added only when both comorbidity and trajectory
        results were supplied.

        The comorbidity page uses the same Louvain modules as the existing Plot
        methods. Its two-dimensional layout is generated automatically with
        deterministic module-aware spring layouts and collision resolution.
        Disease-node radius is derived from the PheWAS case/exposed count, and
        its default display scale is normalized to the network canvas. Edge
        width is derived from the adjusted association strength capped at an
        odds ratio of 10. The resulting graph is also saved as
        ``data/comorbidity.gexf``. The trajectory explorer uses the same module
        assignments, count-based node sizing, and directed temporal paths as the
        standalone trajectory view. The responsive 3D page embeds the output of
        :meth:`three_dimension_plot` with disease-sphere sizes normalized to the
        displayed network, preserving the configured ``exposure_name`` or
        ``outcome_name`` behavior.

        The website is self-contained: Plotly, cleaned result payloads, CSS, and
        JavaScript are written under the output directory and can be opened from
        ``index.html`` without a web server or internet connection. On repeated
        calls, the generated ``assets``, ``data``, and ``pages`` directories are
        replaced; unrelated files directly under ``path`` are preserved.

        Args:
            path: Output directory for the generated website.

        Returns:
            None.

        Examples:
            >>> plot = Plot(
            ...     phewas_result=phewas_result,
            ...     comorbidity_result=comorbidity_result,
            ...     trajectory_result=trajectory_result,
            ...     exposure_name="Exposure",
            ... )
            >>> plot.interactive_website("results/interactive_site")

            For a before-outcome analysis, construct ``Plot`` with
            ``outcome_name`` instead of ``exposure_name`` and call the same
            method.
        """
        from .interactive_visualization import build_interactive_website

        build_interactive_website(self, path)

    def three_dimension_plot(
        self,
        path: str,
        max_radius: Optional[float]=180.0,
        min_radius: Optional[float]=35.0,
        line_color: Optional[str]="black",
        line_width: Optional[float]=1.0,
        size_reduction: Optional[float]=0.5,
        cluster_reduction_ratio: Optional[float]=1,
        layer_distance: Optional[float]=40.0,
        layout_width: Optional[float]=900.0,
        layout_height: Optional[float]=900.0,
        font_style: Optional[str]='Times New Roman',
        font_size: Optional[float]=15.0,
        cluster_max_attempts: Optional[int]=100,
    ) -> None:
        """Generates an interactive 3D disease trajectory network.

        The plot uses shared disease-module labels from direct NetworkX Louvain
        clustering. If cluster labels are not already cached, Louvain is run
        ``cluster_max_attempts`` times with different seeds and the highest
        modularity partition is retained. Disease nodes are then positioned with
        bounded random sampling inside cluster sectors.

        For after-exposure analysis, specifying ``exposure_name`` draws the
        exposure as the origin of trajectories. For before-outcome analysis,
        specifying ``outcome_name`` draws a terminal outcome node below the
        disease network and connects terminal disease nodes to it.

        Args:
            path: Output path for the interactive HTML file.
            max_radius: Maximum radial distance for disease nodes.
            min_radius: Minimum radial distance for disease nodes.
            line_color: Color used for trajectory and outcome-link lines.
            line_width: Width used for trajectory and outcome-link lines.
            size_reduction: Scaling factor applied to disease-node sizes.
            cluster_reduction_ratio: Amount of angular spacing between modules.
            layer_distance: Vertical distance between trajectory-order layers.
            layout_width: Plotly layout width.
            layout_height: Plotly layout height.
            font_style: Font family used in the plot.
            font_size: Font size for labels.
            cluster_max_attempts: Number of seeded NetworkX Louvain runs used
                when clusters are not already cached. Defaults to 100.
        """
        import plotly.graph_objects as go

        self.__prepare_network_data("three_dimension_plot", require_trajectory=True)
        cluster_weight = self.comorbidity_beta_col
        if not self.__check_node_attrs("cluster"):
            self.__cluster(cluster_weight, max_attempts=cluster_max_attempts)
        self.__ensure_trajectory_nodes_clustered("three_dimension_plot")
        if not self.__check_node_attrs("order"):
            self.__trajectory_order()
            self.__comorbidity_order()

        self.__make_location_random(
            max_radius,
            min_radius,
            layer_distance,
            cluster_reduction_ratio
        )

        plot_data = []

        # Exposure marker
        if self._exposure:
            plot_data.append(
                go.Scatter3d(
                    x=[self._exposure_location[0]],
                    y=[self._exposure_location[1]],
                    z=[self._exposure_location[2]],
                    mode='markers',
                    marker=dict(symbol='circle', size=self._exposure_size, color='black'),
                    name=self._exposure_name,
                    legendgrouptitle_text='Origin of Trajectories',
                    showlegend=True,
                    customdata=[str(self._exposure)],
                    hovertemplate=f"{self._exposure_name}<extra></extra>",
                )
            )

        #spheres (Surface) + trajectory edges (Scatter3d lines)
        plot_data += self.__plot(line_width, line_color, size_reduction)

        if self._outcome_name:
            outcome_location = self.__outcome_3d_location(layer_distance)
            terminal_nodes = self.__terminal_trajectory_nodes()
            ox, oy, oz = [], [], []
            for node in terminal_nodes:
                node_location = self._nodes_attrs[node]["location"]
                ox += [node_location[0], outcome_location[0], None]
                oy += [node_location[1], outcome_location[1], None]
                oz += [node_location[2], outcome_location[2], None]
            if terminal_nodes:
                plot_data.append(
                    go.Scatter3d(
                        x=ox,
                        y=oy,
                        z=oz,
                        mode="lines",
                        line=dict(color=line_color, width=line_width),
                        name="Outcome links",
                        showlegend=True,
                        hoverinfo="skip",
                    )
                )
            plot_data.append(
                go.Scatter3d(
                    x=[outcome_location[0]],
                    y=[outcome_location[1]],
                    z=[outcome_location[2]],
                    mode="markers",
                    marker=dict(symbol="circle", size=self._outcome_size, color="black"),
                    name=self._outcome_name,
                    legendgrouptitle_text="Terminal Outcome",
                    showlegend=True,
                    customdata=["outcome"],
                    hovertemplate=f"{self._outcome_name}<extra></extra>",
                )
            )

        #locate the base trajectory edge trace created by __plot
        edge_trace_idx = None
        for i, tr in enumerate(plot_data):
            if getattr(tr, "type", None) == "scatter3d" and ("lines" in (tr.mode or "")) and (tr.name == "Trajectories"):
                edge_trace_idx = i
                break
        if edge_trace_idx is None:
            for i, tr in enumerate(plot_data):
                if getattr(tr, "type", None) == "scatter3d" and ("lines" in (tr.mode or "")):
                    edge_trace_idx = i
                    break
        if edge_trace_idx is None:
            raise RuntimeError("Could not find the trajectory edge trace (Scatter3d lines).")

        ex = list(plot_data[edge_trace_idx].x or [])
        ey = list(plot_data[edge_trace_idx].y or [])
        ez = list(plot_data[edge_trace_idx].z or [])

        #build adjacency from your trajectory dataframe
        edges = list(zip(self._trajectory[self._source], self._trajectory[self._target]))
        adj_starts = {}  # node_id(str) -> list of start indices into ex/ey/ez
        for ei, (s, t) in enumerate(edges):
            start = 3 * ei  # because __get_edge_attrs appends [src, tgt, None] per edge :contentReference[oaicite:4]{index=4}
            ss, tt = str(s), str(t)
            adj_starts.setdefault(ss, []).append(start)
            adj_starts.setdefault(tt, []).append(start)

        #prebuild one highlight line trace per node, initially hidden
        hl_trace_index = {}
        hl_line = dict(width=max(2.0, line_width * 4.0), color="rgba(197,76,130,1.0)")

        for nid, starts in adj_starts.items():
            hx, hy, hz = [], [], []
            for s in starts:
                if s + 1 >= len(ex):
                    continue
                hx += [ex[s], ex[s + 1], None]
                hy += [ey[s], ey[s + 1], None]
                hz += [ez[s], ez[s + 1], None]

            plot_data.append(
                go.Scatter3d(
                    x=hx, y=hy, z=hz,
                    mode="lines",
                    line=hl_line,
                    name=f"__hl__{nid}",
                    showlegend=False,
                    hoverinfo="skip",
                    visible=False
                )
            )
            hl_trace_index[nid] = len(plot_data) - 1

        #node centers used to map Surface clicks (sphere surface point) -> nearest node center
        node_ids, node_x, node_y, node_z = [], [], [], []
        if self._exposure:
            node_ids.append(str(self._exposure))
            node_x.append(float(self._exposure_location[0]))
            node_y.append(float(self._exposure_location[1]))
            node_z.append(float(self._exposure_location[2]))

        for n, attrs in self._nodes_attrs.items():
            loc = attrs["location"]
            node_ids.append(str(n))
            node_x.append(float(loc[0]))
            node_y.append(float(loc[1]))
            node_z.append(float(loc[2]))

        axis = dict(
            showbackground=False,
            showline=False,
            zeroline=False,
            showgrid=False,
            showticklabels=False,
            title=''
        )

        layout = go.Layout(
            title=dict(
                text="Three-Dimensional Disease Network",
                font=dict(size=30, family=font_style),
                x=0.45
            ),
            width=layout_width,
            height=layout_height,
            showlegend=True,
            scene=dict(
                xaxis=dict(axis),
                yaxis=dict(axis),
                zaxis=dict(axis)
            ),
            margin=dict(t=100),
            hovermode='closest',
            clickmode="event",
            legend=dict(
                title=dict(text='Figure legend'),
                font=dict(family=font_style, size=font_size),
                itemclick=False
            ),
            font=dict(family=font_style),
            meta=dict(
                __hl_trace_index=hl_trace_index,
                __hl_node_ids=node_ids,
                __hl_node_x=node_x,
                __hl_node_y=node_y,
                __hl_node_z=node_z,
                __hl_edge_trace_idx=edge_trace_idx
            )
        )

        fig = go.Figure(data=plot_data, layout=layout)

    # inject JavaScript for interactive highlighting
        post_script = r"""
    (function() {
    var gd = document.getElementById('{plot_id}');
    if (!gd) return;

    var meta = (gd.layout && gd.layout.meta) ? gd.layout.meta : {};
    var hlMap = meta.__hl_trace_index || {};
    var nodeIds = meta.__hl_node_ids || [];
    var nx = meta.__hl_node_x || [];
    var ny = meta.__hl_node_y || [];
    var nz = meta.__hl_node_z || [];

    function nearestNodeId(x, y, z) {
        var best = -1, bestD = Infinity;
        for (var i = 0; i < nodeIds.length; i++) {
        var dx = nx[i] - x, dy = ny[i] - y, dz = nz[i] - z;
        var d = dx*dx + dy*dy + dz*dz;
        if (d < bestD) { bestD = d; best = i; }
        }
        return (best >= 0) ? String(nodeIds[best]) : null;
    }

    function setVisible(traceIndex, vis) {
        if (traceIndex === null || traceIndex === undefined) return Promise.resolve();
        return Plotly.restyle(gd, {visible: vis}, [traceIndex]);
    }

    function clearHighlight() {
        var prev = gd.__hlPrevIndex;
        gd.__hlPrevIndex = null;
        gd.__hlPrevNode = null;
        return setVisible(prev, false);
    }

    (function addButton(){
        var btn = document.createElement('button');
        btn.textContent = 'Clear highlight';
        btn.style.cssText = 'margin:6px 0; padding:6px 10px; cursor:pointer;';
        gd.parentNode.insertBefore(btn, gd);
        btn.addEventListener('click', function(){ clearHighlight(); });
    })();

    gd.on('plotly_click', function(ev) {
        if (!ev || !ev.points || !ev.points.length) return;

        var p = ev.points[0];
        var tr = p.data || {};

        // Ignore clicking on any line trace (base edges or highlight edges)
        if (tr.type === 'scatter3d' && String(tr.mode || '').indexOf('lines') !== -1) return;

        var nid = null;

        // Exposure marker provides customdata
        if (p.customdata !== undefined && p.customdata !== null) {
        nid = String(p.customdata);
        } else {
        // Surface spheres do not have node ids per point, map by nearest center
        nid = nearestNodeId(p.x, p.y, p.z);
        }

        if (!nid) return;

        var idx = hlMap[nid];
        if (idx === undefined || idx === null) {
        // no incident edges stored (isolated node)
        clearHighlight();
        return;
        }

        // toggle off
        if (gd.__hlPrevNode === nid) {
        clearHighlight();
        return;
        }

        var prevIdx = gd.__hlPrevIndex;
        gd.__hlPrevIndex = idx;
        gd.__hlPrevNode = nid;

        // Only toggle visibility, avoid rewriting x/y/z arrays (much faster in 3D)
        requestAnimationFrame(function() {
        setVisible(prevIdx, false).then(function() {
            return setVisible(idx, true);
        });
        });
    });

    })();
    """
        fig.write_html(
            path,
            include_plotlyjs=True,
            full_html=True,
            post_script=post_script)
        del fig, plot_data
        gc.collect()

    def comorbidity_network_plot(
        self, 
        path :str,
        max_radius: Optional[float]=180.0,
        min_radius: Optional[float]=35.0,
        size_reduction: Optional[float]=0.5,
        cluster_reduction_ratio: Optional[float]=1,
        line_width: Optional[float]=1.0,
        line_color: Optional[str]="black",
        layer_distance: Optional[float]=40.0,
        font_style: Optional[str]="Times New Roman",
        cluster_max_attempts: Optional[int]=100
    ) -> None:
        """Generates a 2D visualization of the comorbidity network.

        Creates a plot showing disease comorbidities as:
        - Disease nodes (phecodes) as colored circles (grouped by disease system)
        - Comorbidity relationships as connecting lines between nodes
        - Node sizes proportional to disease significance
        - Color coding by disease system/category
        - Random cluster-sector positions shared with the 3D plot

        Args:
            path: Output file path for saving HTML visualization
            max_radius: Maximum radial position for nodes (default: 180.0)
            min_radius: Minimum radial position for nodes (default: 35.0)
            size_reduction: Scaling factor for node sizes (default: 0.5)
            cluster_reduction_ratio: Amount of angular spacing between modules
                (default: 1)
            line_width: Width of comorbidity lines (default: 1.0)
            line_color: Color of comorbidity lines (default: "black")
            layer_distance: Vertical spacing used when trajectory order is
                available (default: 40.0)
            font_style: Font family for text elements (default: "Times New Roman")
            cluster_max_attempts: Number of seeded NetworkX Louvain runs used
                when clusters are not already cached. Defaults to 100.

        Workflow:
            1. Checks/calculates shared module assignments if missing
            2. Computes node orders if missing
            3. Generates bounded random node positions
            4. Creates visualization with:
            - Comorbidity edges as connecting lines
            - Disease nodes as colored circles
            5. Saves interactive plot to HTML file

        Example:
            >>> network.comorbidity_network_plot(
                    "comorbidity.html",
                    max_radius=50,
                    line_color="blue",
                    size_reduction=0.7
                )

        Note:
            - Output is an interactive HTML file using Plotly
            - All distance parameters are in arbitrary units
            - First node in each system shows in legend (others hidden)
            - Hover shows disease name
            - ``outcome_name`` does not add an outcome node to this
              comorbidity-only network view

        """
        self.__prepare_network_data("comorbidity_network_plot")
        if not self.__check_node_attrs("cluster"):
            self.__cluster(
                self.comorbidity_beta_col,
                max_attempts=cluster_max_attempts
            )
        self.__ensure_trajectory_nodes_clustered("comorbidity_network_plot")
        if not self.__check_node_attrs("order"):
            if self._trajectory.empty:
                self.__update_node_attrs(
                    order={node: 1 for node in self._nodes_attrs}
                )
                self._network_attrs.update({"order number": 1})
            else:
                self.__trajectory_order()
                self.__comorbidity_order()
        self.__make_location_random(
            max_radius,
            min_radius,
            layer_distance,
            cluster_reduction_ratio
        )

        fig = go.Figure()
        pairs = self._comorbidity[[self._source, self._target]].values

        for source, target in pairs:
            x_values = [
                self._nodes_attrs[source]["location"][0], 
                self._nodes_attrs[target]["location"][0]
            ]
            y_values = [
                self._nodes_attrs[source]["location"][1], 
                self._nodes_attrs[target]["location"][1]
            ]
            fig.add_trace(
                go.Scatter(
                    x=x_values,
                    y=y_values,
                    mode='lines',
                    line=dict(
                        color=line_color,
                        width=line_width
                    ),
                    showlegend=False
                )
            )

        for sys in self._system_color.keys():
            nodes = [
                x for x in self._commorbidity_nodes 
                if self._nodes_attrs[x]["system"]==sys
            ]
            is_showlegend = True
            for node in nodes:
                x_axis = self._nodes_attrs[node]["location"][0]
                y_axis = self._nodes_attrs[node]["location"][1]
                theta = np.linspace(0, 2 * np.pi, 100)
                r = self._nodes_attrs[node]["size"]*size_reduction
                x_axis_values = x_axis + r * np.cos(theta)
                y_axis_values = y_axis + r * np.sin(theta)

                if node != nodes[0]:
                    is_showlegend = False

                fig.add_trace(go.Scatter(
                    x=x_axis_values, 
                    y=y_axis_values,
                    fill='toself',
                    mode="lines",
                    fillcolor=self._nodes_attrs[node]["color"],
                    showlegend=is_showlegend,
                    hovertemplate='{text}<extra></extra>',
                    text=self._nodes_attrs[node]["name"],
                    name='%s' % (self.sys_dict[sys]),
                    legendgroup='sphere', 
                    legendgrouptitle_text='Disease Systems',
                    line=dict(color=self._nodes_attrs[node]["color"], width=1)
                ))

        fig.update_layout(
            title=dict(
                text="Comorbidity Network", 
                font=dict(size=30, family=font_style),
                x=0.45
            ),
            showlegend=True,
            xaxis=dict(visible=False, scaleanchor="y"), 
            yaxis=dict(visible=False),
            plot_bgcolor='white',
            font=dict(family=font_style),
            hovermode='closest',
            margin=dict(t=100),
            width=900.0,
            height=900.0,
        )

        py.plot(fig, filename=path, auto_open=False)
        del fig
        gc.collect()

    def trajectory_plot(
        self, 
        path: str,
        dpi: Optional[float]=500,
        cluster_max_attempts: Optional[int]=100
    ) -> None:
        """Generates and saves trajectory visualizations for each disease cluster.

        Creates 2D network plots showing disease trajectories within each cluster,
        with nodes positioned hierarchically based on trajectory relationships.
        Each cluster is saved as a separate image file. If ``outcome_name`` is
        specified, a terminal outcome node is added below each cluster trajectory
        and terminal disease nodes are linked to it.

        Args:
            path: Directory path to save output images
            dpi: Image resolution in dots per inch for output files (default: 500)
            cluster_max_attempts: Number of seeded NetworkX Louvain runs used
                when clusters are not already cached. Defaults to 100.

        Workflow:
            1. Performs shared NetworkX Louvain module detection if not already done
            2. Identifies significant trajectories
            3. For each cluster:
            - Creates hierarchical layout
            - Generates visualization with:
                * Nodes colored by disease type
                * Edges weighted by trajectory strength
                * Exposure disease marked specially (if exists)
                * Terminal outcome node marked specially (if outcome_name exists)
            - Saves as PNG image

        Example:
            >>> network.trajectory_plot(
                    "output/plots"
                )

        Note:
            - Outputs one PNG file per cluster (named 'cluster_<number>.png')
            - Uses matplotlib for static visualization
            - Exposure disease (if exists) appears as grey node
            - Outcome disease (if exists) appears as a grey terminal node
            - Node sizes proportional to disease significance
            - Edge widths proportional to trajectory strength
        """
        source = self.source_col
        target = self.target_col
        cluster_weight = self.comorbidity_beta_col

        self.__prepare_network_data("trajectory_plot", require_trajectory=True)
        if not self.__check_node_attrs("cluster"):
            self.__cluster(cluster_weight, max_attempts=cluster_max_attempts)
        self.__ensure_trajectory_nodes_clustered("trajectory_plot")
        if self._exposure:
            exposure = self._exposure
        else:
            exposure = 0

        #check the existence of path
        if not os.path.exists(path):
            raise ValueError(f'No such directory: {path}')
        
        def rotate(angle_,valuex,valuey,pointx,pointy):
            valuex = np.array(valuex)
            valuey = np.array(valuey)
            if angle_ > 0:
                angle = math.radians(angle_)
                sRotatex = (valuex-pointx)*math.cos(angle) + (valuey-pointy)*math.sin(angle) + pointx
                sRotatey = (valuey-pointy)*math.cos(angle) - (valuex-pointx)*math.sin(angle) + pointy
                return (sRotatex,sRotatey)     
            elif angle_ < 0:
                angle = math.radians(abs(angle_))
                nRotatex = (valuex-pointx)*math.cos(angle) - (valuey-pointy)*math.sin(angle) + pointx
                nRotatey = (valuex-pointx)*math.sin(angle) + (valuey-pointy)*math.cos(angle) + pointy
                return (nRotatex,nRotatey) 
            else:    
                return (float(valuex),float(valuey))
        
        def angle_lst(n, angle_step=12):
            if n <= 1:
                if n%2 == 1:
                    rl_n = int((n-1)/2)
                    return [i*10 for i in range(rl_n,0,-1)] + [0] + [i*-10 for i in range(1,rl_n+1)]
                else:
                    rl_n = int(n/2)
                    return [5+(i-1)*10 for i in range(rl_n,0,-1)] + [(i-1)*-10-5 for i in range(1,rl_n+1)]
            else:
                if n % 2 == 1:  
                    half = (n - 1) // 2
                    positive = [angle_step * i for i in range(half, 0, -1)]
                    negative = [-angle_step * i for i in range(1, half+1)]
                    return positive + [0] + negative
                else:  
                    half = n // 2
                    offset = angle_step / 2 
                    positive = [offset + angle_step * (i-1) for i in range(half, 0, -1)]
                    negative = [-offset - angle_step * (i-1) for i in range(1, half+1)]
                    return positive + negative

        def sort_arc(lst, y, r=500):
            n_dots = len(lst)
            x_pos = {}
            for dot in lst:
                angle = angle_lst(n_dots)[lst.index(dot)]
                x_pos[dot] = rotate(angle, 500, y, 500, y-r)
            return x_pos
            
        def hierarchy_layout(df, start_node, method='prox'):
            # Create a directed graph from the dataframe
            G = nx.from_pandas_edgelist(df, source=source, target=target, create_using=nx.DiGraph())

            if G.number_of_nodes() == 0:
                return {}, 150
            
            if start_node not in G:
                G.add_node(start_node)
                first_layer = [n for n, deg in G.in_degree() if deg == 0 and n != start_node]
                for n in first_layer:
                    G.add_edge(start_node, n)   
            
            # Compute strongly connected components (SCCs)
            scc = list(nx.strongly_connected_components(G))
            node_to_scc = {}
            for idx, comp in enumerate(scc):
                for node in comp:
                    node_to_scc[node] = idx

            # Condense the graph into a DAG where each node is a SCC
            C = nx.condensation(G, scc)

            # Use BFS on the condensed graph to assign layers
            start_comp = node_to_scc[start_node]
            comp_layer = {start_comp: 1}
            queue = [start_comp]
            while queue:
                current = queue.pop(0)
                current_layer = comp_layer[current]
                for neighbor in C.successors(current):
                    # Update the layer if this gives a longer path from the start
                    if neighbor not in comp_layer or comp_layer[neighbor] < current_layer + 1:
                        comp_layer[neighbor] = current_layer + 1
                        queue.append(neighbor)
            
            # Map each original node to its component's layer {phecode: number_layer}
            d_lst_layer = {node: comp_layer[node_to_scc[node]] for node in G.nodes() if node_to_scc[node] in comp_layer}

            n_layer = max([x for x in d_lst_layer.values()])
            
            if n_layer <= 5:
                height = 150
            else:
                height = 1500/(n_layer-1)

            layer_y = {i:1500-height*(i-1) for i in range(1, n_layer+1)}
            
            new_pos = {}
            for l in range(1, n_layer+1):
                if l == 1:
                    new_pos[start_node] = (500, 1500)
                else:
                    # find all of nodes in this layer
                    d_lst = [x for x in d_lst_layer.keys() if d_lst_layer[x]==l]
                    temp_pos = sort_arc(d_lst, layer_y[l])
                    for d in d_lst:
                        new_pos[d] = temp_pos[d]
            
            for l in range(3,n_layer+1):
                d_lst = [x for x in d_lst_layer.keys() if d_lst_layer[x]==l]
                l_last = np.arange(2,l)
                d_lst_last = [x for x in d_lst_layer.keys() if d_lst_layer[x] in l_last]
                pos_dict = {i:j for i,j in zip(np.arange(len(d_lst)),[new_pos[d] for d in d_lst])}
                dis_dict = {}
                for d in d_lst:
                    d_connected = [x for x in df.loc[df[target]==d][source].values if x in d_lst_last]
                    dis_d_dict = {}
                    for pos_index in pos_dict.keys():
                        current_pos = pos_dict[pos_index]
                        dis_d_dict[pos_index] = distance(current_pos,d_connected,new_pos)
                    dis_d_dict_order = sorted(dis_d_dict.items(), key = lambda kv:(kv[1], kv[0]))
                    dis_d_dict_ = {i:j for i,j in dis_d_dict_order}
                    dis_dict[d] = dis_d_dict_
                    
                if method == 'exact':
                    dis_dict_iter = {}
                    for p in itertools.permutations(pos_dict.keys()):
                        dis_dict_iter[sum([dis_dict[d][i] for d,i in zip(d_lst,p)])] = p
                    pos_min = dis_dict_iter[min([x for x in dis_dict_iter.keys()])]
                    for d,i in zip(d_lst,pos_min):
                        new_pos[d] = pos_dict[i]
                else:
                    d_max_min = {}
                    for d in d_lst:
                        d_max_min[d] = max(dis_dict[d].values())-min(dis_dict[d].values())
                    d_lst_sorted = sorted(d_max_min.items(), key = lambda kv:(kv[1], kv[0]))
                    d_lst_sorted = [x[0] for x in d_lst_sorted][::-1] 
                    pos_index_occupy = []
                    for d in d_lst_sorted:
                        for pos_index in list(dis_dict[d].keys())[::-1]:
                            if pos_index not in pos_index_occupy:
                                pos_index_occupy.append(pos_index)
                                new_pos[d] = pos_dict[pos_index]
                                break
                            else:
                                continue
            return new_pos,height

        def distance(pos_0,dot_lst,pos_dict_):
            total = 0
            for dot in dot_lst:
                pos_1 = pos_dict_[dot]
                x_ = (pos_0[0] - pos_1[0])**2
                y_ = (pos_0[1] - pos_1[1])**2
                total += (x_ + y_)**0.5
            return total
        
        def ratio(number: float):
            return np.min([number/200*10, 5])+0.2

        tra = self._trajectory.copy()
        tra.index = np.arange(len(tra))
        for idx in tra.index:
            if tra.loc[idx, self._source] == exposure:
                tra.drop(idx, inplace=True)
            else:
                source_cluster = self._nodes_attrs[tra.loc[idx, self._source]]["cluster"]
                target_cluster = self._nodes_attrs[tra.loc[idx, self._target]]["cluster"]
                if source_cluster != target_cluster:
                    tra.drop(idx, inplace=True)

        tra["source_cluster"] = tra[self._source].apply(lambda x: self._nodes_attrs[x]["cluster"])

        coef_dict = {
            (tra.loc[i, self._source], tra.loc[i, self._target]): ratio(tra.loc[i,'n_total']) 
            for i in tra.index
        }

        all_nodes = self.__get_nodes(
            tra,
            self._source,
            self._target
        )

        clusters = set(
            [
                self._nodes_attrs[node]["cluster"] 
                for node in all_nodes
            ]
        )

        for cluster in clusters:
            df = tra[tra["source_cluster"]==cluster]
            source_nodes = df[~df[self._source].isin(df[self._target].values)][self._source].values
            for node in source_nodes:
                temp_df = pd.DataFrame(
                    [['%i-%.1f' % (exposure, node), exposure, node]],
                    columns=['name', self._source, self._target]
                )
                df = pd.concat([df, temp_df])
            df.index = np.arange(len(df))

            position, height = hierarchy_layout(df, exposure)
            graph = nx.DiGraph()
            for idx in df.index:
                graph.add_edge(
                    df.loc[idx, self._source],df.loc[idx, self._target])
            
            # Ensure every graph node has a position
            missing = [n for n in graph.nodes if n not in position]
            if missing:
                # place missing/disconnected nodes on an extra lower arc
                if position:
                    y_bottom = min(y for _, y in position.values()) - height
                else:
                    y_bottom = 1500
                position.update(sort_arc(missing, y_bottom))

            outcome_node = "__dinetxify_outcome__"
            outcome_edges = []
            if self._outcome_name:
                terminal_nodes = [
                    node for node in graph.nodes
                    if node != exposure and graph.out_degree(node) == 0 and node in self._nodes_attrs
                ]
                if len(terminal_nodes) == 0:
                    terminal_nodes = [
                        node for node in graph.nodes
                        if node != exposure and node in self._nodes_attrs
                    ]
                if terminal_nodes:
                    y_bottom = min(y for _, y in position.values()) - height
                    position[outcome_node] = (500, y_bottom)
                    outcome_edges = [(node, outcome_node) for node in terminal_nodes]

            fig, ax_nx = plt.subplots(dpi=600, figsize=(6,4))
            plt.axis("off")

            if self._exposure:
                edges = [x for x in list(graph.edges) if x[0]==exposure]
                nodes = list(
                    set([x for edge in edges for x in edge])
                )

                labels = {}
                for node in nodes:
                    if node==exposure:
                        labels.update({node:self._exposure_name})
                    else:
                        labels.update(
                            {node:self._nodes_attrs[node]["name"]}
                        )
                nx.draw_networkx(
                    graph,
                    position,
                    node_color=[
                        "grey" if x==exposure 
                        else self._nodes_attrs[x]["color"] 
                        for x in nodes
                    ],
                    arrowsize=6,
                    width=[1]*len(edges),
                    node_size=[
                        100 if x==exposure 
                        else np.pi*self._nodes_attrs[x]["size"]**(2)
                        for x in nodes
                    ],
                    font_size=4,
                    connectionstyle="arc3,rad=0",
                    labels=labels,
                    edge_color=["grey"]*len(edges),
                    ax=ax_nx,
                    style=":",
                    edgelist=edges,
                    nodelist=nodes
                )

            edges = [x for x in list(graph.edges) if x[0]!=exposure]
            nodes = list(
                set([x for edge in edges for x in edge])
            )

            nx.draw_networkx(
                graph,
                position,
                node_color=[self._nodes_attrs[x]["color"] for x in nodes],
                arrowsize=6,
                width=[coef_dict[x]**(1/3) for x in edges],
                node_size=[np.pi*self._nodes_attrs[x]["size"]**(2) for x in nodes],
                font_size=4,
                connectionstyle="arc3,rad=0",
                labels={x:self._nodes_attrs[x]["name"] for x in nodes},
                edge_color=["lightgray"]*len(edges),
                ax=ax_nx,
                edgelist=edges,
                nodelist=nodes
            )

            if outcome_edges:
                outcome_graph = nx.DiGraph()
                outcome_graph.add_edges_from(outcome_edges)
                nx.draw_networkx_edges(
                    outcome_graph,
                    position,
                    edgelist=outcome_edges,
                    arrowsize=6,
                    width=[1] * len(outcome_edges),
                    edge_color=["grey"] * len(outcome_edges),
                    ax=ax_nx,
                    style=":",
                )
                nx.draw_networkx_nodes(
                    outcome_graph,
                    position,
                    nodelist=[outcome_node],
                    node_color=["grey"],
                    node_size=[np.pi * self._outcome_size ** 2],
                    ax=ax_nx,
                )
                nx.draw_networkx_labels(
                    outcome_graph,
                    position,
                    labels={outcome_node: self._outcome_name},
                    font_size=4,
                    ax=ax_nx,
                )

            fig.savefig(
                path+'/cluster_%i.png' % (cluster),
                bbox_inches='tight',
                dpi=dpi
                )
            plt.close(fig)
        gc.collect()

    def phewas_plot(
        self,
        path: str,
        system_font_size: Optional[float]=17,
        disease_font_size: Optional[float]=10,
        HR_max: Optional[float]=2,
        incident_number_max: Optional[int]=None,
        exposed_only_cohort: Optional[bool]=False,
        dpi: Optional[float]=200
    ) -> None:
        """Generates a circular PheWAS (Phenome-Wide Association Study) plot.

        Creates a polar bar plot visualizing disease associations across different
        disease categories (systems). Before-outcome cohort results display
        hazard ratios from time-varying Cox models, while nested case-control
        results display odds ratios from conditional logistic models. Node
        counts come from ``N_cases_outcome_group`` when automatically selected
        by ``Plot``.

        The plot includes:
        - Outer ring showing individual disease associations
        - Inner segments grouping by disease system
        - Color gradient indicating effect size (hazard ratio or odds ratio)
        - Automatic text rotation for readability

        Args:
            path: Output file path for saving the plot
            system_font_size: Font size for disease system/category labels (default: 17)
            disease_font_size: Font size for disease labels (default: 10)
            HR_max: Upper bound for the exponentiated-effect heatmap. Values greater than or equal to this render as the same red. Affects color only. (default: 2)
            incident_number_max: Upper bound for the incident count heatmap for exposed-only cohort. Values greater than or equal to this render as the same red. None means auto scale to the max observed count. (default: None)
            exposed_only_cohort: Whether the PheWAS is performed in an exposed-only cohort. If True, the incident count heatmap will be used instead of HR heatmap. (default: False)
            dpi: Image resolution in dots per inch for output files (default: 200)

        Example:
            >>> network.phewas_plot(
                    "phewas_plot.png",
                )

        Note:
            - Uses random effects model for system-level estimates
            - Positive associations shown in red, negative in green
            - Output resolution is controlled by ``dpi``
            - Plot includes:
            * Color bar legend for effect sizes
            * System category labels
            * Individual disease labels
        """
        col_significant = self.phewas_significance_col
        col_coef = self.phewas_coef_col
        col_se = self.phewas_se_col
        col_disease = self.disease_col
        col_system = self.system_col
        col_exposure = self.phewas_number_col
        effect_label = 'Hazard ratio'
        if 'model_type' in self._phewas.columns:
            model_types = self._phewas['model_type'].astype(str).str.lower()
            if model_types.str.contains('conditional logistic').any():
                effect_label = 'Odds ratio'
        #for exposed only cohort, HR_max is not used
        if exposed_only_cohort:
            HR_max = None
            if incident_number_max is not None:
                if incident_number_max <=0 or not isinstance(incident_number_max, int)  :
                    raise ValueError("incident_number_max should be int and larger than 0.")
        elif exposed_only_cohort == False:
            incident_number_max = None

        def random_effect(coef_lst, se_lst):
            if len(coef_lst)==1:
                return coef_lst[0]
            
            #calculate fixed effect result
            w = 1/np.square(se_lst)
            u = (np.sum(w*coef_lst))/(np.sum(w))
            
            #random effect result
            w = np.reshape(w,[len(w)])
            q = np.sum(w*np.square((coef_lst-u)))
            df = len(w)-1
            c = np.sum(w)-(np.sum(w*w)/np.sum(w))
            tau2 = (q-df)/c
            w2 = 1/(np.square(se_lst)+tau2)
            u2 = (np.sum(w2*coef_lst))/(np.sum(w2))
            # seu2 = np.sqrt(1/np.sum(w2))
            return u2

        def sys_mean(df):
            sys_dict = {}
            sys_lst = set(df[col_system].values)
            for sys in sys_lst:
                temp = df.loc[df[col_system]==sys]
                if exposed_only_cohort:
                    mean = random_effect(temp[col_exposure].values, temp[col_exposure].values)
                else:
                    mean = random_effect(temp[col_coef].values, temp[col_se].values)
                sys_dict[mean] = sys
            sys_dict_ = [sys_dict[i] for i in sorted([x for x in sys_dict.keys() if not pd.isna(x)])]
            sys_dict_ += [x for x in sys_dict.values() if x not in sys_dict_]
            return sys_dict_
        
        fig, ax = plt.subplots(
            subplot_kw=dict(polar=True),
            figsize=(20, 20),
            nrows=1,
            ncols=1
        )
        cmap = plt.get_cmap("tab20c")
        max_pos = np.log(HR_max)
        cmap_pos = plt.get_cmap("Reds")
        if exposed_only_cohort:
            phe_df = self._phewas[self._phewas[col_significant]==True]
        else:
            phe_df = self._phewas.loc[(self._phewas[col_coef]>0) & (self._phewas[col_significant]==True)]

        phe_df = phe_df.sort_values(by=col_disease)

        if exposed_only_cohort:
            if incident_number_max is None:
                incident_number_max = phe_df[col_exposure].max()
            incident_number_min = phe_df[col_exposure].min()
            q0,q1,q2,q3 = np.linspace(incident_number_min, incident_number_max, 4)
            phe_df["col_color"] = phe_df[col_exposure].apply(
                lambda x: x/incident_number_max
            )

        else:
            q0,q1,q2,q3 = np.linspace(np.log(1), max_pos, 4)
            phe_df["col_color"] = phe_df[col_coef].apply(
                lambda x: x/max_pos
            )

        phe_df['color'] = phe_df["col_color"].apply(
            lambda x: cmap_pos(x)
        )

        size = 0.1
        edge_width_n = 0.4
        start = 0
        n_system = len(set(phe_df[col_system].values))
        n_total = len(phe_df) + n_system*edge_width_n
        sys_order = sys_mean(phe_df)
        for system in sys_order:
            temp_df = phe_df.loc[phe_df[col_system]==system]
            number = len(temp_df)
            width = np.array([2*np.pi/n_total]*number)
            left = np.cumsum(np.append(start, width[:-1]))
            colors_outter = temp_df['color'].values
            ax.bar(
                x=left,
                width=width, 
                bottom=1-size, 
                height=size,
                color=colors_outter, 
                edgecolor='w', 
                linewidth=1, 
                align="edge"
            )
            ax.bar(
                x=start,
                width=width.sum(),
                bottom=0,
                height=1-1.1*size,
                color=cmap(next(iter([19]*16))),
                edgecolor='w', 
                linewidth=1, 
                align="edge",
                alpha=0.5
            )
            x_system = (left[0] + left[-1] + 2*np.pi/n_total)/2
            if x_system<=0.5*np.pi or x_system>=1.5*np.pi:
                ax.text(
                    x_system,
                    0.15,
                    self.sys_dict[system],
                    ha='left',
                    va='center',
                    rotation=np.rad2deg(x_system),
                    rotation_mode="anchor",
                    fontsize=system_font_size
                )
            else:
                ax.text(
                    x_system,
                    0.15,
                    self.sys_dict[system],
                    ha='right',
                    va='center',
                    rotation=np.rad2deg(x_system)+180,
                    rotation_mode="anchor",
                    fontsize=system_font_size
                )       
            left_text = left + width/2
            rotations = np.rad2deg(left_text)

            for x, rotation, label in zip(left_text, rotations, temp_df[col_disease].values):
                if x<=0.5*np.pi or x>=1.5*np.pi:
                    ax.text(
                        x,
                        1.02, 
                        label,
                        ha='left', 
                        va='center', 
                        rotation=rotation, 
                        rotation_mode="anchor",
                        fontsize=disease_font_size
                    )
                else:
                    ax.text(
                        x,
                        1.02, 
                        label, 
                        ha='right', 
                        va='center', 
                        rotation=rotation+180, 
                        rotation_mode="anchor",
                        fontsize=disease_font_size
                    )
            start = left[-1] + width[0]*(1+edge_width_n)

        if exposed_only_cohort:
            norm = mpl.colors.Normalize(
                vmin=incident_number_min, 
                vmax=incident_number_max
            )

            sm = cm.ScalarMappable(norm=norm ,cmap=cmap_pos)
            bar = plt.colorbar(
                sm, 
                ax=ax,
                location='bottom', 
                label='Incident number', 
                shrink=0.4
            )
        else:
            norm = mpl.colors.Normalize(vmin=0, vmax=max_pos)
            sm = cm.ScalarMappable(norm=norm ,cmap=cmap_pos)
            bar = plt.colorbar(
                sm, 
                ax=ax,
                location='bottom', 
                label=effect_label,
                shrink=0.4
            )

        tick_locator = ticker.MaxNLocator(nbins=3)
        bar.locator = tick_locator
        bar.update_ticks()
        if exposed_only_cohort:
            bar.set_ticks([q0,q1,q2,q3])
            bar.set_ticklabels([int(x) for x in [q0,q1,q2,q3]])
        else:
            bar.set_ticks([q0,q1,q2,q3])
            bar.set_ticklabels([f'{np.exp(x):.1f}' for x in [q0,q1,q2,q3]])

        ax.set_axis_off()
        plt.savefig(
            path, 
            dpi=dpi, 
            bbox_inches='tight'
        )
        plt.close(fig)
        gc.collect()
