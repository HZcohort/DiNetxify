# -*- coding: utf-8 -*-
"""
Created on Wed Jan 1 19:48:09 2025

@author: Haowen Liu - Biomedical Big data center of West China Hospital, Sichuan University
"""
import community as community_louvain
import matplotlib.ticker as ticker
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import plotly.offline as py
import matplotlib.cm as cm
import networkx as nx
import pandas as pd
import numpy as np
import itertools
import random
import math
import os

from collections import (
    defaultdict,
    deque,
    Counter
)

from typing import (
    Any,
    List,
    Dict,
    Tuple,
    Optional,
)

Df = pd.DataFrame

class ThreeDimensionalNetwork(object):
    def __init__(
        self, 
        phewas_result: Df,
        comorbidity_result: Df, 
        trajectory_result: Df,
        disease_system: Optional[List[str]] | None=None,
        exposure: Optional[float] | None=None,
        exposure_location: Optional[Tuple[float]] | None=None,
        exposure_size: Optional[float] | None=None,
        source: Optional[str]='phecode_d1',
        target: Optional[str]='phecode_d2',
        phewas_phecode: Optional[str]='phecode',
        phewas_number: Optional[str]='N_cases_exposed',
        system_col: Optional[str]='system',
        col_disease_pair: Optional[str]='name_disease_pair',
        filter_phewas_col: Optional[str]='phewas_p_significance',
        filter_comorbidity_col: Optional[str]='comorbidity_p_significance',
        filter_trajectory_col: Optional[str]='trajectory_p_significance',
        **kwargs
    ):
        """initialize the ThreeDimensionalDiseaseNetwork class.

        Args:
            comorbidity_result (Df): Result of comorbidity network analysis.
            trajectory_result (Df): Result of trajectory network analysis.
            phewas_result (Df): Result of PHEWAS analysis.
            exposure (float, optional): Phecode of exposure in the cohort studies.
                Defaults to None.

            exposure_location (Tuple[float], optional): Three dimension location 
                of exposure to plot. Defaults to None.

            exposure_size (float, optional): Size of exposure. Defaults to None.
            source (str, optional): Column name of D1. Defaults to 'phecode_d1'.
            target (str, optional): Column name of D2. Defaults to 'phecode_d2'.
        """
        # filter the results
        phewas_result, comorbidity_result, trajectory_result = self.__filter_significant(
            phewas_result,
            comorbidity_result,
            trajectory_result,
            filter_phewas_col,
            filter_comorbidity_col,
            filter_trajectory_col
        )

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
            "#F1C40F",
            '#9B59B6',
            '#4ECDC4',
            '#6A5ACD' 
        ]

        if system_col:
            SYSTEM = list(phewas_result[system_col].unique())

        SYSTEM = kwargs.get("SYSTEM", SYSTEM)
        COLOR = kwargs.get("COLOR", COLOR)
        if len(SYSTEM) > len(COLOR):
            raise ValueError(
                f"the length of SYSTEM is more than that of COLOR"
            )
        else:
            COLOR = COLOR[0: len(SYSTEM)]
        
        if disease_system:
            outside_system = [
                sys for sys in disease_system
                if sys not in SYSTEM
            ]
            SYSTEM = disease_system
            if outside_system:
                raise ValueError(
                    f"The system of {outside_system} is not support"
                )
        system_color = dict(
            zip(
                SYSTEM,
                COLOR
            )
        )
        system_color = kwargs.get("system_color", system_color)

        # check the inclusion relation between trajectory and comorbidity
        self.__check_disease_pairs(
            trajectory_result,
            comorbidity_result,
            source,
            target
        )

        # concat the trajectory and comorbidity in vertical level
        df = trajectory_result.copy()
        df.columns = comorbidity_result.columns

        comorbidity_result = pd.concat(
            [comorbidity_result, df],
            axis=0,
            ignore_index=True
        )

        comorbidity_result.drop_duplicates(
            subset=[source, target],
            inplace=True,
            ignore_index=True,
            keep="first"
        )

        if exposure:
            trajectory_result = self.__sequence(
                trajectory_result,
                exposure,
                source,
                target,
                col_disease_pair
            )

        self.__init_attrs(
            comorbidity = comorbidity_result,
            trajectory = trajectory_result,
            phewas = phewas_result,
            exposure = exposure,
            exposure_location = exposure_location,
            exposure_size = exposure_size,
            source = source,
            target = target,
            describe = pd.read_csv(
                os.path.join(
                    os.path.dirname(__file__),
                    "data/phecode_1.2/phecode_info.csv"
                )
            ),
            commorbidity_nodes = self.__get_nodes(
                comorbidity_result,
                source,
                target
            ),
            trajectory_nodes = self.__get_nodes(
                trajectory_result,
                source,
                target
            ),
            system_color = system_color,
            nodes_attrs = {},
            network_attrs = {}
        )
        
        self.__make_node_basic_attrs(
            phewas_phecode,
            phewas_number
        )

    @staticmethod
    def __check_disease_pairs(
        tra_df: Df,
        com_df: Df,
        source: str,
        target: str
    ) -> None:
        """_summary_

        Args:
            tra_df (Df): _description_
            com_df (Df): _description_
            d1_str (str): _description_
            d2_str (str): _description_
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
                Warning(
                    "Disease pairs of trajectory network has \
                    not been included comorbidity network"
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
        """_summary_

        Args:
            df (Df): _description_
            exposure (float): _description_
            d1_str (str): _description_
            d2_str (str): _description_

        Returns:
            Df: _description_
        """
        df = df.loc[df[source].isin(df[target].values)]
        first_layer = set(df[source].values)
        d1_d2 = [[exposure, d, '%f-%f' % (exposure,d)] for d in first_layer]
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
        """_summary_

        Args:
            df (Df): _description_
            d1_str (str): _description_
            d2_str (str): _description_

        Returns:
            set: _description_
        """
        return set(df[source].to_list() + df[target].to_list())
    
    @staticmethod
    def __split_name(name: str) -> str:
        """line feed the name of disease.

        Args:
            name (str): Disease name.

        Returns:
            str: Disease name lined feed.
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
        center: Tuple[float],
        r: float
    ) -> Tuple[float]:
        """get the cordinate of sphere.

        Args:
            center (tuple): the location of centre point.
            r (float): the radius of sphere.

        Returns:
            tuple: the cordinate of sphere.
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
        target_lst :List[float]
    ) -> Dict[float, int]:
        """_summary_

        Args:
            source_lst (List[float]): _description_
            target_lst (List[float]): _description_

        Returns:
            Dict[int]: _description_
        """
        adj = defaultdict(list)
        in_degree = defaultdict(int)
        nodes = set(source_lst) | set(target_lst)
        
        for u, v in zip(source_lst, target_lst):
            adj[u].append(v)
            in_degree[v] += 1
        
        distance = {node: 1 for node in nodes}
        
        queue = deque()
        for node in nodes:
            if in_degree[node] == 0:
                queue.append(node)
        
        while queue:
            u = queue.popleft()
            for v in adj[u]:
                if distance[u] + 1 > distance[v]:
                    distance[v] = distance[u] + 1
                in_degree[v] -= 1
                if in_degree[v] == 0:
                    queue.append(v)
        return distance
    
    @staticmethod
    def __most_frequent_element(lst: List[Any]) -> Any:
        counter = Counter(lst)
        most_common = counter.most_common(1)
        return most_common[0][0]

    def __check_node_attrs(self, key: str) -> bool:
        """_summary_

        Args:
            key (str): _description_

        Returns:
            bool: _description_
        """
        for _, attrs in self._nodes_attrs.items():
            if key in attrs.keys():
                is_exist = True
            else:
                is_exist = False
        return is_exist

    def __calculate_ratio(
        self, 
        cluster_nodes: Dict[int, float]
    ) -> Dict[int, float]:
        for cluster, nodes in cluster_nodes.items():
            size = [self._nodes_attrs[x]["size"] for x in nodes]
            sum_size = sum(size)
            cluster_nodes.update({cluster:sum_size})

        total_size = sum(cluster_nodes.values())
        for key, value in cluster_nodes.items():
            cluster_nodes.update({key:value/total_size})
        return cluster_nodes
    
    def __init_attrs(self, **kwargs) -> None:
        """_summary_
        """
        for key, value in kwargs.items():
            setattr(self, f"_{key}", value)
    
    def __update_node_attrs(self, **kwargs) -> None:
        """_summary_
        """
        for key, value in kwargs.items():
            for node, attr in value.items():
                if node in self._nodes_attrs.keys():
                    self._nodes_attrs[node].update({f"{key}":attr})
    def __filter_significant(
        self,
        phewas_result,
        comorbidity_result,
        trajectory_result,
        filter_phewas_col,
        filter_comorbidity_col,
        filter_trajectory_col
    ) -> Df:
        if filter_phewas_col:
            phewas_result = phewas_result.loc[
                phewas_result[filter_phewas_col] == True
            ]
        
        if filter_comorbidity_col:
            comorbidity_result = comorbidity_result.loc[
                comorbidity_result[filter_comorbidity_col] == True
            ]

        if filter_trajectory_col:
            trajectory_result = trajectory_result.loc[
                trajectory_result[filter_trajectory_col] == True
            ]
        return phewas_result, comorbidity_result, trajectory_result
    
    def __get_same_nodes(self, key_name:str) -> Dict[Any, Any]:
        """_summary_

        Args:
            key_name (str): _description_

        Returns:
            Dict[Any]: _description_
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

    def __sig_nodes(self) -> List[List[float]]:
        """_summary_

        Returns:
            List[List[float]]: _description_
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
        """get the atrribution of sphere to plot.

        Args:
            center (Tuple[float]): Three dimension location of centre point.
            r (float): Radius of sphere.
            name (str): Name of disease.
            label (str): the label of disease.
            color (str): the color of sphere to plot.
            light_dict (dict, optional): the atrributions of light in the ployly module, more details in plotly. 
                                         Defaults to dict(ambient=0.2, diffuse=0.8, specular=0.4, roughness=0.2, fresnel=2.0).
            light_position_dict (dict, optional): the location of the light. Defaults to dict(x=1.5, y=1.5, z=1.5).

        Returns:
            x (iterable object): the cordinates in the x-axis.
            y (iterable object): the cordinates in the y-axis.
            z (iterable object): the cordinates in the z-axis.
            colorscale_ (list[list]): the color of sphere
            light_dict (dict{str:float}): the atrributions of light in the ployly module. 
                                          Defaults to dict(ambient=0.2, diffuse=0.8, specular=0.4, roughness=0.2, fresnel=2.0).
            label (str): the label of disease.
            name (str): the name of disease.
            light_position_dict (dict{str:float}): the location of the light. Defaults to dict(x=1.5, y=1.5, z=1.5).
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
            max_attempts :Optional[int]=10000
        ) -> None:
        """calculate the location of center point of sphere.

        Args:
            node (float): phecode.
            _hash_dict (dict): {cluster:[phecodes in the cluster]}.
            begin_angle (float, optional): the minimum of angle in the arc. Defaults to 0.25.
            end_angele (float, optional): the maximum of angle in the arc. Defaults to 0.25.
            _iter_time (int, optional): the time of iteration. Defaults to 10000.

        Returns:
            _location_dict (dict{float:{tuple}}): the dict of {phecodes:location(x,y,z)}.
            _hash_dict (dict{float:list}): the dict of {cluster:[phecodes in the cluster]}.
        """
        if self._exposure_location is None:
            self._exposure_location = (0, 0, 0)

        cluster_location = {x:{} for x in range(self._network_attrs["cluster number"])}
        for node, attrs in self._nodes_attrs.items():
            is_sep = True
            order = attrs["order"]
            cluster = attrs["cluster"]
            max_ang = 2*math.pi*sum([cluster_ratio[i] for i in range(cluster+1)])
            min_ang = max_ang - 2*math.pi*cluster_ratio[cluster]

            for _ in range(max_attempts):
                radius = random.uniform(min_radius, max_radius)
                node_ang = random.uniform(
                    min_ang + cluster_reduction_ratio/2 * 2*math.pi*cluster_ratio[cluster],
                    max_ang - cluster_reduction_ratio/2 * 2*math.pi*cluster_ratio[cluster]
                )
                node_loc = (
                    radius * math.cos(node_ang),
                    radius * math.sin(node_ang),
                    self._exposure_location[2] - order*z_axis
                )

                for oth_node, loc in cluster_location[cluster].items():
                    distance = math.sqrt(
                        (loc[0]-node_loc[0])**2+
                        (loc[1]-node_loc[1])**2
                    )
                    sum_length = sum(
                        [
                            self._nodes_attrs[node]["size"],
                            self._nodes_attrs[oth_node]["size"],
                        ]
                    )
                    if sum_length < distance:
                        is_sep = False
                        break

                if is_sep:
                    cluster_location[cluster].update({node:node_loc})
                    self._nodes_attrs[node].update({"location":node_loc})
                    break

            cluster_location[cluster].update({node:node_loc})
            self._nodes_attrs[node].update({"location":node_loc})
    
    def __make_node_basic_attrs(
        self,
        phewas_phecode: str,
        phewas_number: str
    ) -> None:

        # disease name attrs
        node_name = {
            node:self.__split_name('%s (%.1f)' % (name, node)) 
            for node, name in self._describe[
                ["phecode", "phenotype"]
            ].values
        }

        # disease system attrs
        node_system = dict(
            zip(
                self._describe["phecode"], 
                self._describe["category"]
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
        max_attempts: Optional[int]=5000,
    ) -> None:
        """Use the louvain algorithm to cluster 
            the disease acording to the relation 
            of comorbidity network.

        Args:
            max_attempts (int, optional): the time of iteration. Defaults to 5000.
            weight (str, optional): the weight of the disease pair (D1 with D2). Defaults to 'comorbidity_beta'.
        """
        # create network class and add the edges
        Graph_position = nx.Graph()
        [
            Graph_position.add_edge(
                row[self._source],
                row[self._target],
                weight=row[weight]
            ) 
            for _, row in self._comorbidity.iterrows()
        ]

        # [
        #     Graph_position.add_edge(
        #         row[self._source],
        #         row[self._target],
        #         weight=row["trajectory_beta"]
        #     ) 
        #     for _, row in self._trajectory.iterrows()
        # ]
        # random and repeated clustering the nodes
        result = []
        for i in range(max_attempts):
            partition = community_louvain.best_partition(
                Graph_position, 
                random_state=i
            )

            result.append([
                i, community_louvain.modularity(
                    partition, 
                    Graph_position
                )
            ])

        result_df = pd.DataFrame(
            result, 
            columns=['rs', 'modularity']
        ).sort_values(by='modularity')
        best_rs = result_df.iloc[-1, 0]

        # final result with the best score
        cluster_ans = community_louvain.best_partition(
            Graph_position, 
            random_state=best_rs
        )

        self.__update_node_attrs(
            cluster = dict(cluster_ans)
        )

        self._network_attrs.update(
            {"cluster number":max(cluster_ans.values()) + 1}
        )

    def __make_location_random(
        self, 
        max_radius: float, 
        min_radius: float,
        distance: float,
        cluster_reduction_ratio: float
    ) -> None:
        """get the three dimension location of nodes(phecodes), using the metod that nodes(phecodes) of 
           one cluster be gathered in one sector in the x-y plane and the latter nodes(phecodes) locates the under layer.

        Args:
            max_radius (float): the maximum of radius in the sector.
            min_radius (float): the minimum of radius in the sector.

        Returns:
            dict: {str:tuple} for example: {"549.4":(1, 2, 3)}
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
        """get the layer number of nodes(phecodes) in the trajectory network (D1->D2).
        """
        tra_df = self._trajectory.loc[
            self._trajectory[self._source]!=self._exposure
        ]
        
        node_order = self.__calculate_order(
            tra_df[self._source].to_list(),
            tra_df[self._target].to_list()
        )

        self.__update_node_attrs(
            order = node_order
        )

        self._network_attrs.update(
            {"order number":max(node_order.values())}
        )

    def __comorbidity_order(self) -> None:
        """get the layer number of nodes(phecodes) in the comorbidity network (D1-D2).
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
        """get the attribution of plot. This method plots the all trajectory(D1->D2), comorbidity(D1-D2), and nodes(phecodes).

        Args:
            line_color (str): the color of line between each node(phecode)
            line_width (float): the width of line between each node(phecode)

        Returns:
            list: the attribution of plot
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
                    name="%s Disease" % (sys.title()),
                    showscale=False,
                    legendgroup="sphere",
                    legendgrouptitle_text="Disease",
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
            legendgrouptitle_text='All Trajectories',
            name='Trajectories',
            showlegend=True,
            hoverinfo=None
        )
        plot_data.append(trace_data)
        return plot_data
    
    def threeDimension_plot(
        self, 
        path: str,
        max_radius: Optional[float]=35.0, 
        min_radius: Optional[float]=180.0,
        line_color: Optional[str]="black", 
        line_width: Optional[float]=1.0,
        size_reduction: Optional[float]=0.5,
        cluster_reduction_ratio: Optional[float]=0.4,
        cluster_weight: str="comorbidity_beta",
        layer_distance: Optional[float]=40.0,
        layout_width: Optional[float]=900.0,
        layout_height: Optional[float]=900.0,
        font_style: Optional[str]='Times New Roman',
        font_size: Optional[float]=15.0,
    ) -> None:
        """plot the three-dimension comorbidity and trajectory network.

        Args:
            max_radius (float): the maximum of radius in the sector.
            min_radius (float): the minimum of radius in the sector.
            plot_method (str): the method of plot, which is one of "full", "half", "compact".
            line_color (str): the color of line in the nodes(phecodes).
            line_width (float): the width of line in the nodes(phecodes).
            layer_distance (float): the distance of two adjoining layers.
            file_name (str): the name of file to save.
            layout_width (float, optional): the width of layout in the figure. Defaults to 900.
            layout_height (float, optional): the height of layout in the figure. Defaults to 900.
            font_style (str, optional): the font style of layout in the figure. Defaults to 'Times New Roman'.
            font_size (float, optional): the font size of layout in the figure. Defaults to 15.

        Raises:
            KeyError: if the augrment plot_method does not be included in "full", "compact", and "half", it will raises KeyError.
        """
        if not self.__check_node_attrs("cluster"):
            self.__cluster(cluster_weight)
        if not self.__check_node_attrs("order"):
            self.__trajectory_order()
            self.__comorbidity_order()
        if not self.__check_node_attrs("location"):
            self.__make_location_random(
                max_radius,
                min_radius,
                layer_distance,
                cluster_reduction_ratio
            )

        plot_data = []

        # plot exposure
        if self._exposure:
            exposure_data = go.Scatter3d(
                x=[self._exposure_location[0]],
                y=[self._exposure_location[1]],
                z=[self._exposure_location[2]],
                mode='markers',
                marker=dict(
                    symbol='circle',
                    size=self._exposure_size,
                    color='black'
                ),
                name="Exposure disease",
                legendgrouptitle_text='Origin of Trajectories',
                showlegend=True
            )
            plot_data += [exposure_data]

        # plot the nodes and edges
        plot_data += self.__plot(
            line_width, 
            line_color, 
            size_reduction
        )
  
        # axis
        axis = dict(
            showbackground=False, 
            showline=False,
            zeroline = False,
            showgrid = False,
            showticklabels = False,
            title=''
        )

        # layout
        layout = go.Layout(
            title=dict(
                text="Three Dimensional Network", 
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
            legend=dict(
                title=dict(text='Trace of clusters'),
                font=dict(family=font_style,size=font_size),
                itemclick=False
            ), 
            font=dict(family=font_style)
        )

        # plot the figure
        fig = go.Figure(data=plot_data, layout=layout)

        # create the file of the figure
        py.plot(fig, filename=path)

    def comorbidity_network_plot(
        self, 
        path :str,
        max_radius: Optional[float]=35.0,
        min_radius: Optional[float]=90.0,
        size_reduction: Optional[float]=0.5,
        cluster_reduction_ratio: Optional[float]=0.4,
        cluster_weight: Optional[str]="comorbidity_beta",
        line_width: Optional[float]=1.0,
        line_color: Optional[str]="black",
        layer_distance: Optional[float]=40.0,
        font_style: Optional[str]="Times New Roman"
    ) -> None:
        """plot the result of commorbidity. The method same to plot_3d just in the plane of x-y.

        Args:
            max_radius (float): the maximum of radius in the sector.
            min_radius (float): the minimum of radius in the sector.
            line_width (float, optional): the width of line in the nodes(phecodes). Defaults to 1.
            line_color (str, optional): the color of line in the nodes(phecodes). Defaults to "black".
        """
        if not self.__check_node_attrs("cluster"):
            self.__cluster(cluster_weight)
        if not self.__check_node_attrs("order"):
            self.__trajectory_order()
            self.__comorbidity_order()
        if not self.__check_node_attrs("location"):
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
                    hovertemplate=self._nodes_attrs[node]["name"],
                    name='%s Disease' % (sys.title()), 
                    legendgroup='sphere', 
                    legendgrouptitle_text='Diseases',
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

        py.plot(fig, filename=path)

    def significant_trajectory_plot(
        self, 
        path: str,
        cluster_weight: Optional[str]="comorbidity_beta",
    ) -> None:
        """plot the incluster trajectory of each cluster.

        Args:
            distance (float): the distance of each nodes(phecodes) in x-y plane.
            layer_distance (float): the distance of two adjoining layers.
            line_width (float): the width of line in the nodes(phecodes).
            line_color (str): the color of line in the nodes(phecodes).
        """
        if not self.__check_node_attrs("cluster"):
            self.__cluster(cluster_weight)

        sig_trajectory = self.__sig_nodes()

        if self._exposure:
            exposure = self._exposure
        else:
            exposure = 0

        if not sig_trajectory:
            raise TypeError("There is no significant trajectory\
                            of network, please try other method of plot")
        
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
        
        def angle_lst(n):
            if n%2 == 1:
                rl_n = int((n-1)/2)
                return [i*10 for i in range(rl_n,0,-1)] + [0] + [i*-10 for i in range(1,rl_n+1)]
            else:
                rl_n = int(n/2)
                return [5+(i-1)*10 for i in range(rl_n,0,-1)] + [(i-1)*-10-5 for i in range(1,rl_n+1)]

        def sort_arc(lst,y,r=300):
            n_dots = len(lst)
            x_pos = {}
            for dot in lst:
                angle = angle_lst(n_dots)[lst.index(dot)]
                x_pos[dot] = rotate(angle,500,y,500,y-r)
            
            return x_pos
            
        def hierarchy_layout(df,method='prox'):
            d_lst_layer = {exposure:1}
            i = 1
            while True:
                d_lst = [x for x in d_lst_layer.keys() if d_lst_layer[x]==i]
                if len(d_lst) == 0:
                    break
                i += 1
                d_next = df.loc[df[self._source].isin(d_lst)][self._target].values
                for d in d_next:
                    d_lst_layer[d] = i
            
            n_layer = max([x for x in d_lst_layer.values()])
            
            if n_layer <= 5:
                height = 150
            else:
                height = 1000/(n_layer-1)
            layer_y = {i:1000-height*(i-1) for i in range(1,n_layer+1)}
            
            new_pos = {}
            for l in range(1,n_layer+1):
                if l == 1:
                    new_pos[exposure] = (500,1000)
                else:
                    d_lst = [x for x in d_lst_layer.keys() if d_lst_layer[x]==l]
                    temp_pos = sort_arc(d_lst,layer_y[l])
                    for d in d_lst:
                        new_pos[d] = temp_pos[d]
            
            for l in range(3,n_layer+1):
                d_lst = [x for x in d_lst_layer.keys() if d_lst_layer[x]==l]
                l_last = np.arange(2,l)
                d_lst_last = [x for x in d_lst_layer.keys() if d_lst_layer[x] in l_last]
                pos_dict = {i:j for i,j in zip(np.arange(len(d_lst)),[new_pos[d] for d in d_lst])}
                dis_dict = {}
                for d in d_lst:
                    d_connected = [x for x in df.loc[df[self._target]==d][self._target].values if x in d_lst_last]
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
            return new_pos

        def distance(pos_0,dot_lst,pos_dict_):
            total = 0
            for dot in dot_lst:
                pos_1 = pos_dict_[dot]
                x_ = (pos_0[0] - pos_1[0])**2
                y_ = (pos_0[1] - pos_1[1])**2
                total += (x_ + y_)**0.5
            return total
        
        def ratio(number: float):
            return np.min([number/200* 5,5])+0.2

        tra = self._trajectory
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
            (tra.loc[i, self._source],tra.loc[i, self._target]):ratio(tra.loc[i,'n_total']) 
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
                    [['%i-%.1f' % (exposure,node), exposure, node]],
                    columns=['name', self._source, self._target]
                )
                df = pd.concat([df, temp_df])
            df.index = np.arange(len(df))

            position = hierarchy_layout(df)
            graph = nx.DiGraph()
            for idx in df.index:
                graph.add_edge(
                    df.loc[idx, self._source],
                    df.loc[idx, self._target]
                )

            fig, ax_nx = plt.subplots(dpi=500,figsize=(7,7))
            plt.axis("off")

            if self._exposure:
                edges = [x for x in list(graph.edges) if x[0]==exposure]
                nodes = list(
                    set([x for edge in edges for x in edge])
                )

                labels = {}
                for node in nodes:
                    if node==exposure:
                        labels.update({node:"Exposure disease"})
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
                    font_size=6,
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
                font_size=6,
                connectionstyle="arc3,rad=0",
                labels={x:self._nodes_attrs[x]["name"] for x in nodes},
                edge_color=["grey"]*len(edges),
                ax=ax_nx,
                edgelist=edges,
                nodelist=nodes
            )

            fig.savefig(path+'/cluster_%i.png' % (cluster))

    def phewas_plot(
        self,
        path: str,
        col_coef: Optional[str]="phewas_coef",
        col_system: Optional[str]="system",
        col_se: Optional[str]="phewas_se",
        col_disease: Optional[str]="disease"
    ) -> None:
        
        def random_effect(coef_lst, se_lst):
            if len(coef_lst)==1:
                return [coef_lst[0], se_lst[0]]
            
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
            seu2 = np.sqrt(1/np.sum(w2))
            return [u2, seu2]

        def sys_mean(df):
            sys_dict = {}
            sys_lst = set(df[col_system].values)
            for sys in sys_lst:
                temp = df.loc[df[col_system]==sys].dropna(subset=[col_coef], how='any')
                mean = random_effect(temp[col_coef].values,temp[col_se].values)[0]
                sys_dict[mean] = sys
            sys_dict_ = [sys_dict[i] for i in sorted([x for x in sys_dict.keys() if not pd.isna(x)])]
            sys_dict_ += [x for x in sys_dict.values() if x not in sys_dict_]
            return sys_dict_
        
        _, ax = plt.subplots(
            subplot_kw=dict(polar=True),
            dpi=500,
            figsize=(20, 20),
            nrows=1,
            ncols=1
        )
        
        cmap = plt.get_cmap("tab20c")
        max_neg = np.log(0.5)
        max_pos = np.log(5.0)
        cmap_neg = plt.get_cmap("Greens")
        cmap_pos = plt.get_cmap("Reds")
        sys_dict = {
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
            'congenital anomalies': 'Congenital anomalies diseases',
            'symptoms':'Symptoms diseases',
            'others':'Others diseases'
        }
        phe_df = self._phewas.loc[self._phewas[col_coef]>0]
        phe_df = phe_df.sort_values(by=col_disease)
        phe_df['color'] = phe_df[col_coef].apply(
            lambda x: cmap_neg(x/max_neg) if x<0 else cmap_pos(x/max_pos)
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
                    sys_dict[system].capitalize(),
                    ha='left',
                    va='center',
                    rotation=np.rad2deg(x_system),
                    rotation_mode="anchor",
                    fontsize=17
                )
            else:
                ax.text(
                    x_system,
                    0.15,
                    sys_dict[system].capitalize(),
                    ha='right',
                    va='center',
                    rotation=np.rad2deg(x_system)+180,
                    rotation_mode="anchor",
                    fontsize=17
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
                        fontsize=10
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
                        fontsize=10
                    )            
            start = left[-1] + width[0]*(1+edge_width_n)
        sm = cm.ScalarMappable(cmap=cmap_pos)
        bar = plt.colorbar(
            sm, 
            ax=ax,
            location='bottom', 
            label='Hazard ratio', 
            shrink=0.4
        )

        tick_locator = ticker.MaxNLocator(nbins=3)
        bar.locator = tick_locator
        bar.update_ticks()
        bar.set_ticklabels(['NA', '1.2', '>2.0'])
        ax.set_axis_off()
        plt.savefig(
            path, 
            dpi=1200, 
            bbox_inches='tight'
        )