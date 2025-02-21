# -*- coding: utf-8 -*-
"""
Created on Wed Jan 1 19:48:09 2025

@author: Haowen Liu - Biomedical Big data center of West China Hospital, Sichuan University
"""
import community as community_louvain
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import plotly.offline as py
import matplotlib.cm as cm
import networkx as nx
import pandas as pd
import numpy as np
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

SYSTEM = [
    'circulatory system', 
    'sense organs', 
    'injuries & poisonings', 
    'neurological',
    'dermatologic', 
    'digestive', 
    'hematopoietic', 
    'musculoskeletal', 
    'endocrine/metabolic', 
    'mental disorders', 
    'infectious diseases',
    'genitourinary',
    'neoplasms',
    'respiratory',
    "symptoms",
    "congenital anomalies",
    "mental disorders",
    "others"
]

COLOR = [
    "#FF5733",
    "#33FF57",
    "#3357FF",
    "#FFFF33",
    "#FF33FF",
    "#33FFFF",
    "#C70039",
    "#900C3F",
    "#581845",
    "#1ABC9C",
    "#2ECC71",
    "#3498DB",
    "#9B59B6",
    "#E74C3C",
    "#F1C40F",
    "#FF7F50",
    "#FFD700", 
    "#A52A2A"
]

class ThreeDimensionalNetwork(object):
    """
    
    """
    def __init__(
        self, 
        phewas_result: Df,
        comorbidity_result: Df, 
        trajectory_result: Df,
        exposure: Optional[float]=None,
        exposure_location: Optional[Tuple[float]]=None,
        exposure_size: Optional[float]=None,
        size_reduction: float=0.1,
        source: Optional[str]='phecode_d1',
        target: Optional[str]='phecode_d2'
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
        self.__check_disease_pairs(
            trajectory_result,
            comorbidity_result,
            source,
            target
        )

        if exposure:
            trajectory_result = self.__sequence(
                trajectory_result,
                exposure,
                source,
                target
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
            trajectory_pairs = [
                [
                    row[source], 
                    row[target]
                ] for _, row in trajectory_result.iterrows()
            ],
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
            nodes_attrs = {},
            network_attrs = {}
        )
        
        for node in self._commorbidity_nodes:
            self._nodes_attrs.update({node:{}})

        for node in self._trajectory_nodes:
            if node not in self._nodes_attrs and node!=exposure:
                self._nodes_attrs.update({node:{}})

        self.__make_node_basic_attrs(size_reduction)

    @staticmethod
    def __check_disease_pairs(
        tra_df: Df,
        com_df: Df,
        d1_str: str,
        d2_str: str
    ) -> None:
        """_summary_

        Args:
            tra_df (Df): _description_
            com_df (Df): _description_
            d1_str (str): _description_
            d2_str (str): _description_
        """
        tra_pairs_lst = [
            [row[d1_str], row[d2_str]]
            for _, row in tra_df.iterrows()
        ]
        com_pairs_lst = [
            [row[d1_str], row[d2_str]]
            for _, row in com_df.iterrows()
        ]
        for pair in tra_pairs_lst:
            if pair not in com_pairs_lst:
                Warning("Disease pair of trajectory network has \
                    not been included comorbidity network")
                break

    @staticmethod
    def __sequence(
        df: Df,
        exposure: float,
        d1_str: str,
        d2_str: str
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
        df = df.loc[df[d1_str].isin(df[d2_str].values)]
        first_layer = set(df[d1_str].values)
        d1_d2 = [[exposure, d, '%f-%f' % (exposure,d)] for d in first_layer]
        d1_d2_df = pd.DataFrame(
            d1_d2,
            columns=[d1_str, d2_str, 'name_disease_pair']
        )
        trajectory_df = pd.concat([df, d1_d2_df])
        return trajectory_df
    
    @staticmethod
    def __get_nodes(
        df: Df,
        d1_str: str,
        d2_str: str
    ) -> set:
        """_summary_

        Args:
            df (Df): _description_
            d1_str (str): _description_
            d2_str (str): _description_

        Returns:
            set: _description_
        """
        return set(df[d1_str].to_list() + df[d2_str].to_list())
    
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

    # @staticmethod
    # def __generate_colors(
    #     n :int, 
    #     colormap :Optional[str]='inferno'
    # ) -> List[str]:
    #     """_summary_

    #     Args:
    #         n (int): _description_
    #         colormap (Optional[str], optional): _description_. Defaults to 'viridis'.

    #     Returns:
    #         List[str]: _description_
    #     """
    #     cmap = cm.get_cmap(colormap, n)
    #     colors = [cmap(i) for i in range(n)]
    #     hex_colors = [
    #         f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}' 
    #         for r, g, b, _ in colors
    #     ]
    #     return hex_colors

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

    # @staticmethod
    # def __calculate_2d_positions(radii: List[float], spacing:float) -> Dict[int, float]:
    #     """calculate the positions of nodes(phecodes).

    #     Args:
    #         radii (list of float): A list of radii of the circles, in the given order.
    #         spacing (float): The additional distance between the circles (i.e., the outer tangent distance).

    #     Returns:
    #         dict: {number: x-coordinate}, where the numbering starts from 1 and corresponds to the input order.
    #     """
    #     circle_positions = {}
    #     position = 0
    #     for i, radius in enumerate(radii):
    #         circle_id = i + 1
    #         if i == 0:
    #             circle_positions[circle_id] = position
    #             position += 2 * radius + spacing  
    #         else:
    #             if i % 2 == 1:
    #                 circle_positions[circle_id] = position  
    #             else:
    #                 circle_positions[circle_id] = -position  
    #             position += 2 * radius + spacing 
    #     return circle_positions

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

    def __calculate_ratio(self, cluster_nodes: Dict[int, float]) -> Dict[int, float]:
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
                self._target
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
            redu_ratio :float,
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
                    min_ang * (1+redu_ratio),
                    max_ang * (1-redu_ratio)
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
                            math.pow(self._nodes_attrs[node]["size"] / 4 * 3 / math.pi, 1/3),
                            math.pow(self._nodes_attrs[oth_node]["size"] / 4 * 3 / math.pi, 1/3)
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
        scale_reduction: float
    ) -> None:
        """_summary_

        Args:
            scale_reduction (Optional[float], optional): _description_. Defaults to 0.1.
        """
        # disease name attrs
        node_name = {
            node:self.__split_name('%s (%.1f)' % (name, node)) 
            for node, name in self._describe[["phecode", "phenotype"]].values
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
                self._phewas["phecode"], 
                self._phewas["N_cases_exposed"]*scale_reduction
            )
        )

        self.__update_node_attrs(
            name = node_name,
            system = node_system,
            size = node_size
        )

        # disease color attrs
        sys_color = dict(
            zip(
                SYSTEM,
                COLOR
            )
        )
        self._network_attrs.update({"system color":sys_color})
        for _, attrs in self._nodes_attrs.items():
            attrs.update({"color":sys_color[attrs["system"]]})

    def __cluster(
        self,
        max_attempts: Optional[int]=5000,
        weight: Optional[str]='comorbidity_beta'
    ) -> None:
        """Use the louvain algorithm to cluster the disease acording to the relation of comorbidity network.

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
        angle: Optional[float]=10.0
    ) -> None:
        """get the three dimension location of nodes(phecodes), using the metod that nodes(phecodes) of 
           one cluster be gathered in one sector in the x-y plane and the latter nodes(phecodes) locates the under layer.

        Args:
            max_radius (float): the maximum of radius in the sector.
            min_radius (float): the minimum of radius in the sector.

        Returns:
            dict: {str:tuple} for example: {549.4:(1, 2, 3)}
        """
        same_cluster_nodes = self.__get_same_nodes("cluster")
        cluster_ratio = self.__calculate_ratio(same_cluster_nodes)
        self.__calculate_location_random(
            max_radius,
            min_radius,
            angle,
            distance,
            cluster_ratio
        )
    
    def __trajectory_order(self) -> None:
        """get the layer number of nodes(phecodes) in the trajectory network (D1->D2).
        """
        node_order = self.__calculate_order(
            self._trajectory[self._source].to_list(),
            self._trajectory[self._target].to_list()
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
        same_cluster_nodes = self.__get_same_nodes("cluster")

        for node, attr in self._nodes_attrs.items():
            if "order" in attr.keys():
                nodes = same_cluster_nodes[attr["cluster"]]
                orders = [self._nodes_attrs[x].get("order", None) for x in nodes]
                orders = list(filter(None, orders))
                if orders:
                    order = self.__most_frequent_element(orders)
                else:
                    order = self._network_attrs["order number"]
                self._nodes_attrs[node].update({"order":order})

        for node, attr in self._nodes_attrs.items():
            if "order" not in attr.keys():
                order = self._network_attrs["order number"]
                self._nodes_attrs[node].update({"order":order})

    def __full_plot(
        self, 
        line_width: float,
        line_color: str, 
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
                    self._nodes_attrs[node]["size"],
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

        com_edges = zip(
            self._comorbidity[self._source],
            self._comorbidity[self._target]
        )

        all_edges = list(tra_edges)
        all_edges.extend(list(com_edges))

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
    
    def __half_plot(
        self, 
        sig_line_width: float, 
        sig_line_color: str,
        no_sig_line_color: Optional[str]="silver", 
        no_sig_line_width: Optional[float]=1.0,
    ) -> List[Any]:
        """get the attribution of plot. This method plots the all trajectory(D1->D2), comorbidity(D1-D2), nodes(phecodes), and highlight their difference.

        Args:
            main_line_width (float): the width of line in the incluster nodes(phecodes).
            nonMain_line_color (str, optional): the color of line in the outcluster nodes(phecodes). Defaults to 'silver'.
            nonMain_line_width (int, optional): the width of line in the outcluster nodes(phecodes). Defaults to 1.

        Returns:
            list: the attribution of plot
        """
        plot_data = []
        sig_trajectory = self.__sig_nodes()

        if sig_trajectory:
            sig_nodes = [
                node for nodes in sig_trajectory for node in nodes 
            ]
        else:
            raise TypeError("There is no significant trajectory \
                            of network, please try other method of plot")

        if self._exposure:
            sig_trajectory = [
                [self._exposure]+path for path in sig_trajectory
            ]

        for sys in SYSTEM:
            nodes = [
                x for x in set(sig_nodes) 
                if self._nodes_attrs[x]["system"]==sys
            ]
            is_showlegend = True
            for node in nodes:
                plot_attrs = self.__sphere_attrs(
                    self._nodes_attrs[node]["location"],
                    self._nodes_attrs[node]["size"],
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

        no_sig_nodes = [x for x in self._nodes_attrs.keys() if x not in sig_nodes]
        is_showlegend = True
        for node in no_sig_nodes:
            plot_attrs = self.__sphere_attrs(
                self._nodes_attrs[node]["location"],
                self._nodes_attrs[node]["size"],
                "grey"
            )

            if node != no_sig_nodes[0]:
                is_showlegend = False

            data = go.Surface(
                x=plot_attrs[0],
                y=plot_attrs[1],
                z=plot_attrs[2],
                colorscale=plot_attrs[3],
                showlegend=is_showlegend,
                lighting=plot_attrs[4],
                hovertemplate=self._nodes_attrs[node]["name"],
                name="Diseases not in significant trajectories",
                showscale=False,
                legendgroup="sphere",
                legendgrouptitle_text="Disease",
                lightposition=plot_attrs[-1]
            )
            plot_data.append(data)

        sig_edges = []
        for nodes in sig_trajectory:
            for i in range(len(nodes)-1):
                sig_edges.append((nodes[i],nodes[i+1]))
        edges_attrs = self.__get_edge_attrs(sig_edges)

        trace_data = go.Scatter3d(
            x=edges_attrs[0],
            y=edges_attrs[1],
            z=edges_attrs[2],
            line=dict(
                color=sig_line_color,
                width=sig_line_width
            ),
            mode='lines',
            legendgroup='trajectories',
            legendgrouptitle_text='All Trajectories',
            name='Significant Trajectories',
            hoverinfo=None
        )

        plot_data.append(trace_data)

        no_sig_edges = []
        pairs = self._comorbidity[[self._source, self._target]].values
        for source, target in pairs:
            if (source, target) not in sig_edges:
                no_sig_edges.append((source, target))

        tra_df = self._trajectory

        pairs = tra_df[[self._source, self._target]].values
        for source, target in pairs:
            if (source, target) not in sig_edges and (source, target) not in no_sig_edges:
                no_sig_edges.append((source, target))
            
        edges_attrs = self.__get_edge_attrs(no_sig_edges)

        trace_data = go.Scatter3d(
            x=edges_attrs[0],
            y=edges_attrs[1],
            z=edges_attrs[2],
            line=dict(
                color=no_sig_line_color,
                width=no_sig_line_width
            ),
            mode='lines',
            legendgroup='trajectories',
            legendgrouptitle_text='All Trajectories',
            name='Insignificant Trajectories',
            hoverinfo=None
        )
        plot_data.append(trace_data)
        return plot_data

    def __compact_plot(
        self,
        sig_line_width: float, 
        sig_line_color: str,
    ) -> List[Any]:
        """get the attribution of plot. This method plots the trajectory(D1->D2) of incluster nodes(phecodes), and incluster nodes(phecodes).

        Args:
            sig_line_width (float): Width of line in the incluster nodes(phecodes).

        Returns:
            list: the attribution of plot
        """
        plot_data = []
        sig_trajectory = self.__sig_nodes()

        if sig_trajectory:
            sig_nodes = [
                node for nodes in sig_trajectory for node in nodes 
            ]
        else:
            raise TypeError("There is no significant trajectory \
                            of network, please try other method of plot")
        
        if self._exposure:
            sig_trajectory = [
                [self._exposure]+path for path in sig_trajectory
            ]

        for sys in SYSTEM:
            nodes = [
                x for x in set(sig_nodes) 
                if self._nodes_attrs[x]["system"]==sys
            ]

            is_showlegend = True
            for node in nodes:
                plot_attrs = self.__sphere_attrs(
                    self._nodes_attrs[node]["location"],
                    self._nodes_attrs[node]["size"],
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

        sig_edges = []
        for nodes in sig_trajectory:
            for i in range(len(nodes)-1):
                sig_edges.append((nodes[i],nodes[i+1]))

        edges_attrs = self.__get_edge_attrs(sig_edges)

        trace_data = go.Scatter3d(
            x=edges_attrs[0],
            y=edges_attrs[1],
            z=edges_attrs[2],
            line=dict(
                color=sig_line_color,
                width=sig_line_width
            ),
            mode='lines',
            legendgroup='trajectories',
            legendgrouptitle_text='All Trajectories',
            name='Significant Trajectories',
            hoverinfo=None
        )

        plot_data.append(trace_data)
        return plot_data

    def threeDimension_plot(
        self, 
        path: str,
        max_radius: float, 
        min_radius: float,
        plot_method: str,
        line_color: str, 
        line_width: float,
        layer_distance: Optional[float]=20.0,
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
            self.__cluster()
        if not self.__check_node_attrs("order"):
            self.__trajectory_order()
            self.__comorbidity_order()
        if not self.__check_node_attrs("location"):
            self.__make_location_random(
                max_radius,
                min_radius,
                layer_distance,
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
        if plot_method == 'full':
            plot_data += self.__full_plot(line_width, line_color)
        elif plot_method == 'compact':
            plot_data += self.__compact_plot(line_width, line_color)
        elif plot_method == 'half':
            plot_data += self.__half_plot(line_width, line_color)
        else:
            raise KeyError(f"This {plot_method} is not exist in the method of plot")

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
                text="", 
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
        max_radius: float,
        min_radius: float,
        line_width: Optional[float]=1.0,
        line_color: Optional[str]="black",
        layer_distance: Optional[float]=20.0,
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
            self.__cluster()
        if not self.__check_node_attrs("order"):
            self.__trajectory_order()
            self.__comorbidity_order()
        if not self.__check_node_attrs("location"):
            self.__make_location_random(
                max_radius,
                min_radius,
                layer_distance,
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

        for sys in SYSTEM:
            nodes = [
                x for x in self._commorbidity_nodes 
                if self._nodes_attrs[x]["system"]==sys
            ]
            is_showlegend = True
            for node in nodes:
                x_axis = self._nodes_attrs[node]["location"][0]
                y_axis = self._nodes_attrs[node]["location"][1]
                theta = np.linspace(0, 2 * np.pi, 100)
                r = self._nodes_attrs[node]["size"]
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
            showlegend=True,
            xaxis=dict(visible=False, scaleanchor="y"), 
            yaxis=dict(visible=False),
            plot_bgcolor='white',
            margin=dict(l=10, r=10, t=10, b=10),
            font=dict(family=font_style),
            hovermode='closest'
        )

        py.plot(fig, filename=path)

    def significant_trajectory_plot(
        self, 
        path: str,
        # distance: float, 
        # layer_distance: float,
        line_color: str,
    ) -> None:
        """plot the incluster trajectory of each cluster.

        Args:
            distance (float): the distance of each nodes(phecodes) in x-y plane.
            layer_distance (float): the distance of two adjoining layers.
            line_width (float): the width of line in the nodes(phecodes).
            line_color (str): the color of line in the nodes(phecodes).
        """
        if not self.__check_node_attrs("cluster"):
            self.__cluster()
        if not self.__check_node_attrs("order"):
            self.__trajectory_order()
            self.__comorbidity_order()

        sig_trajectory = self.__sig_nodes()

        if self._exposure:
            sig_trajectory = [
                [self._exposure]+path for path in sig_trajectory
            ]

        if not sig_trajectory:
            raise TypeError("There is no significant trajectory\
                            of network, please try other method of plot")
        
        for cluster in range(self._network_attrs["cluster number"]):
            nodes = []
            graph = nx.DiGraph()
            # postion = {}
            
            for sig_tra in sig_trajectory:
                if self._nodes_attrs[sig_tra[1]]["cluster"] == cluster:
                    edges = [(sig_tra[i], sig_tra[i+1]) for i in range(len(sig_tra)-1)]
                    graph.add_edges_from(edges)
                    nodes.extend(sig_tra)
            
            # for order in range(self._network_attrs["order number"]):
            #     same_order_nodes = [
            #         node for node in nodes if self._nodes_attrs[node]["order"]==order+1
            #     ]
            #     same_order_nodes_size = [
            #         self._nodes_attrs[node] for node in same_order_nodes
            #     ]
            #     x_axis = self.__calculate_2d_positions(same_order_nodes_size, distance)
            #     for idx, x_axis_value in x_axis.items():
            #         postion.update({
            #             same_order_nodes[idx-1]:(x_axis_value, (order+1)*layer_distance)
            #         })
            
            # color = [self._nodes_attrs[node]["color"] for node in graph.nodes()]
            # size = [self._nodes_attrs[node]["size"] for node in graph.nodes()]
            # label = {node: self._nodes_attrs[node]["name"] for node in graph.nodes()}
            if nodes:
                plt.figure(figsize=(16, 16))
                nx.draw_networkx(
                    graph,
                    # postion,
                    with_labels=True,
                    # labels=label,
                    # node_color=color,
                    node_size=1000,
                    font_size=12, 
                    font_weight='bold',
                    edge_color=line_color,
                    arrows=True,
                    arrowsize=10,
                )

                plt.title("")
                plt.axis("off")
                plt.savefig(f"{path}_cluster{cluster}.png")