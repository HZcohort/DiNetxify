import plotly.offline as py
import plotly.graph_objects as go
import community as community_louvain
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import copy
import random
import pandas as pd
import math
import networkx as nx
import os

class ThreeDimensionalDiseaseNetwork():
    def __init__(self, commorbidity_network_result:pd.DataFrame, 
                 disease_trajectory_result:pd.DataFrame,
                 phewas_result:pd.DataFrame,
                 exposure_disease:float=9999,
                 exposure_disease_location:tuple=(0,0,0),
                 exposure_disease_size:float=1,
                 source:str='phecode_d1',
                 target:str='phecode_d2'):
        # primary attributions
        self.__module_dir = os.path.dirname(__file__)
        self.__commorbidity_network_result = commorbidity_network_result
        self.__phecode_number = dict(zip(phewas_result["phecode"], phewas_result["N_cases_exposed"]))
        self.__describe = pd.read_csv(os.path.join(self.__module_dir, "data/phecode_1.2/phecode_info.csv"))
        self.__exposure_disease = exposure_disease
        self.__exposure_disease_size = exposure_disease_size
        self.__exposure_disease_location = exposure_disease_location
        self.__system = ['circulatory system', 'sense organs', 'injuries & poisonings', 'neurological',
                         'dermatologic', 'digestive', 'hematopoietic', 'musculoskeletal', 
                         'endocrine/metabolic', 'mental disorders', 'infectious diseases', 
                         'genitourinary', 'neoplasms', 'respiratory', "symptoms", "congenital anomalies", 
                         "mental disorders", "others"]
        self.__system_color = ["#FF5733","#33FF57","#3357FF","#FFFF33","#FF33FF","#33FFFF","#C70039",
                               "#900C3F","#581845","#1ABC9C","#2ECC71","#3498DB","#9B59B6","#E74C3C",
                               "#F1C40F","#FF7F50","#FFD700", "#A52A2A"]
        
        # make sure the first layer of phecodes (exposure_disease->phecodes)
        df = disease_trajectory_result.loc[~disease_trajectory_result['phecode_d1'].isin(disease_trajectory_result[target].values)]
        df_lst = set(df[source].values)
        d0_d1 = []
        for d in df_lst: d0_d1.append([exposure_disease, d, '%f-%f' % (exposure_disease,d)])
        d0_d1 = pd.DataFrame(d0_d1,columns=['phecode_d1', 'phecode_d2', 'name_disease_pair'])

        # primary attributions
        self.__disease_trajectory = pd.concat([disease_trajectory_result, d0_d1])
        self.__trajectory_disease_pairs = [[row[source], row[target]] for _, row in self.__disease_trajectory.iterrows()]
        self.__commorbidity_nodes = set(commorbidity_network_result[source].to_list() + commorbidity_network_result[target].to_list())
        self.__trajectory_nodes = set(disease_trajectory_result[source].to_list() + disease_trajectory_result[target].to_list())
        
    @staticmethod
    def __split_name(name:str) -> str:
        """line feed the name of disease.

        Args:
            name (str): disease name.

        Returns:
            str: the disease name lined feed.
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
    def __get_dimension(dimension_dict:dict, node:float) -> float:
        """get the dimension of the phecode.

        Args:
            dimension_dict (dict): {dimension:phecode}.
            node (float): phecode.

        Returns:
            float: the number of dimension.
        """
        for dimension_number, item in dimension_dict.items():
            if node in item: return dimension_number
            
    @staticmethod
    def __except_dimension(nodes:list, dimension_dict:dict) -> list:
        """get the dimension of except phecodes.

        Args:
            nodes (list): phecodes.
            dimension_dict (dict): {dimension:phecode}.

        Returns:
            list: the dimension of nodes.
        """
        nodes_dimension = []
        for item in nodes: nodes_dimension.append(dimension_dict.get(item))
        nodes_dimension = list(filter(None, nodes_dimension))
        return nodes_dimension

    @staticmethod
    def __ball_cordinate(center:tuple, r:float) -> tuple:
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

    def plotly_ball(self, center:tuple, 
                    r:float, 
                    name:str, 
                    label:str, 
                    color:str, 
                    light_dict=dict(ambient=0.2,diffuse=0.8,specular=0.4,roughness=0.2,fresnel=2.0),
                    light_position_dict=dict(x=1.5,y=1.5,z=1.5)):
        """get the atrribution of sphere to plot.

        Args:
            center (tuple): the location of centre point.
            r (float): the radius of sphere.
            name (str): the name of disease.
            label (str): the label of disease.
            color (str): the color of sphere to plot.
            light_dict (dict, optional): the atrributions of light in the ployly module, more details in plotly. Defaults to dict(ambient=0.2, diffuse=0.8, specular=0.4, roughness=0.2, fresnel=2.0).
            light_position_dict (dict, optional): the location of the light. Defaults to dict(x=1.5, y=1.5, z=1.5).

        Returns:
            x_ (iterable object): the cordinates in the x-axis.
            y_ (iterable object): the cordinates in the y-axis.
            z_ (iterable object): the cordinates in the z-axis.
            colorscale_ (list[list]): the color of sphere
            light_dict (dict{str:float}): the atrributions of light in the ployly module. Defaults to dict(ambient=0.2, diffuse=0.8, specular=0.4, roughness=0.2, fresnel=2.0).
            label (str): the label of disease.
            name (str): the name of disease.
            light_position_dict (dict{str:float}): the location of the light. Defaults to dict(x=1.5, y=1.5, z=1.5).
        """
        x, y, z = self.__ball_cordinate(center, r)
        colorscale = [[0.0, color], [0.5, color], [1.0, color]]
        return x, y, z, colorscale, light_dict, label, name, light_position_dict


    def __calculate_location_random(self, node:float, 
                                    _hash_dict:dict,
                                    begin_angle=0.25,
                                    end_angele=0.25,
                                    _iter_time=10000):
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
        # get the cluster number
        cluster_number = self.final_cluster[node]
        # get the dimension number
        for dimension_number, nodes in self.dimension.items():
            if node in nodes: break
        # calculate the corbinate of sphere
        if _hash_dict[cluster_number]:
            iter_time, flag = 0, True
            while flag:
                radius = random.uniform(self.__min_radius, self.__max_radius)
                angle = 2 * math.pi * (self.cumulative_nodes_percentage[self.final_cluster[node]+1]\
                                       -self.cumulative_nodes_percentage[self.final_cluster[node]])\
                                       * random.uniform(begin_angle, 1-end_angele)
                loction = (radius * math.cos(2 * math.pi * self.cumulative_nodes_percentage[self.final_cluster[node]]+angle),
                           radius * math.sin(2 * math.pi * self.cumulative_nodes_percentage[self.final_cluster[node]]+angle),
                           self.height[dimension_number])
                i, hash_list = 0, [x for x in _hash_dict[cluster_number].keys()]
                while i < len(_hash_dict[cluster_number]):
                    distance = math.sqrt((_hash_dict[cluster_number][hash_list[i]][0]-loction[0]) ** 2 +
                                         (_hash_dict[cluster_number][hash_list[i]][1]-loction[1]) ** 2)
                    length = math.pow(self.get_node_size(node) / 4 * 3 / math.pi, 1/3)
                    length += math.pow(self.get_node_size(hash_list[i]) / 4 * 3 / math.pi, 1/3)
                    i += 1
                    if distance < length: break
                iter_time += 1
                if i == len(_hash_dict[cluster_number]) or iter_time == _iter_time: flag = False
        else:
            radius = random.uniform(self.__min_radius, self.__max_radius)
            angle = 2*math.pi*(self.cumulative_nodes_percentage[self.final_cluster[node]+1]\
                               -self.cumulative_nodes_percentage[self.final_cluster[node]])\
                                *random.uniform(begin_angle, 1-end_angele)
            loction = (radius*math.cos(2*math.pi*self.cumulative_nodes_percentage[self.final_cluster[node]]+angle),
                       radius*math.sin(2*math.pi*self.cumulative_nodes_percentage[self.final_cluster[node]]+angle),
                       self.height[dimension_number])
        _location_dict = {node:loction}
        return _location_dict, _hash_dict

    def get_node_size(self, node:float) -> float:
        """get the radius of sphere.

        Args:
            node (float): phecode.

        Returns:
            float: the radius of sphere.
        """
        size = np.sqrt(self.__phecode_number[node])
        return size

    def cluster(self, iter_time:int=5000, 
                source:str='phecode_d1',
                target:str='phecode_d2',
                weight:str='comorbidity_beta') -> dict:
        """Use the louvain algorithm to cluster the disease acording to the relation of comorbidity network.

        Args:
            iter_time (int, optional): the time of iteration. Defaults to 5000.
            source (str, optional): the D1 (from D1 to D2). Defaults to 'phecode_d1'.
            target (str, optional): the D2 (from D1 to D2). Defaults to 'phecode_d2'.
            weight (str, optional): the weight of the disease pair (D1 with D2). Defaults to 'comorbidity_beta'.

        Returns:
            dict: {cluster:list[phecodes]}
        """
        if hasattr(self, "final_cluster"):
            return self.final_cluster
        # create network class and add the edges
        Graph_position = nx.Graph()
        [Graph_position.add_edge(row[source],row[target],weight=row[weight]) 
         for _, row in self.__commorbidity_network_result.iterrows()]
        # random and repeated clustering the nodes
        result = []
        for i in range(iter_time):
            partition = community_louvain.best_partition(Graph_position, random_state=i)
            result.append([i, community_louvain.modularity(partition, Graph_position)])
        random_df = pd.DataFrame(result, columns=['rs', 'modularity']).sort_values(by='modularity')
        best_rs = random_df.iloc[-1, 0]
        # final result with the best score
        self.final_cluster = community_louvain.best_partition(Graph_position, random_state=best_rs)
        self.maximum_class_number = max(self.final_cluster.values()) + 1
        return self.final_cluster

    def location_random(self, max_radius:float, 
                        min_radius:float):
        """_summary_

        Args:
            max_radius (float): _description_
            min_radius (float): _description_

        Returns:
            _type_: _description_
        """
        if not hasattr(self, "final_cluster"):
            self.cluster()
        if not hasattr(self, "dimension"):
            self._dimension()
        self.__max_radius = max_radius
        self.__min_radius = min_radius
        # save the nodes size of each cluster and calculate the sum
        location_dict, hash_dict = {}, {number:{} for number in range(self.maximum_class_number)}
        inclass_nodes_size, sum_size = {number:[] for number in range(self.maximum_class_number)}, 0
        # record nodes size of each cluster
        for dimension_number in range(len(self.dimension)):
            for node in self.dimension[dimension_number]:
                node_size = self.get_node_size(node)
                inclass_nodes_size[self.final_cluster[node]] += [node_size]
                sum_size += node_size
        # calculate the percentage of the sum of nodes size
        nodes_percentage_size_incluster = []
        for cluster_number in range(len(inclass_nodes_size)):
            nodes_percentage_size_incluster.append(sum(inclass_nodes_size[cluster_number]) / sum_size)
        # transform the percentage of the nodes size in one cluster
        inintal_percentage = 0
        the_last_one_cumulative_percentage = -1
        cumulative_nodes_percentage = [inintal_percentage]
        for percentage in nodes_percentage_size_incluster:
            cumulative_percentage = cumulative_nodes_percentage[the_last_one_cumulative_percentage] + percentage
            cumulative_nodes_percentage += [cumulative_percentage]
        self.cumulative_nodes_percentage = cumulative_nodes_percentage
        # calculate the 3D coordinate
        for dimension_number in range(len(self.dimension)):
            for node in self.dimension[dimension_number]:
                one_dimension_location, hash_dict = self.__calculate_location_random(node, hash_dict)
                location_dict.update(one_dimension_location)
                hash_dict[self.final_cluster[node]].update(one_dimension_location)
        location_dict.update({self.__exposure_disease:self.__exposure_disease_location})
        self.location_dict = location_dict
        return self.location_dict
    
    def trajectory(self, source='phecode_d1', 
                   target='phecode_d2'):
        """_summary_

        Args:
            source (str, optional): _description_. Defaults to 'phecode_d1'.
            target (str, optional): _description_. Defaults to 'phecode_d2'.

        Returns:
            _type_: _description_
        """
        if hasattr(self, "trajectory_dimension"):
            return self.trajectory_dimension
        # make sure the dimension of trajectory nodes
        temp_disease, trajectory_dimension, dimension_number = [], {}, -1
        trajectory_df = self.__disease_trajectory.copy()
        while True:
            dimension_number += 1
            # first layers of nodes connecting self.__original_disease
            if len(temp_disease) == 0:
                # nodes of this layer
                temp_disease = [x for x in trajectory_df.loc[trajectory_df[source]==self.__exposure_disease][target]]
                # index of the edges in DataFrame
                index = [x for x in trajectory_df.loc[trajectory_df[source]==self.__exposure_disease].index]
                trajectory_df.drop(index)
                # save in the dict
                trajectory_dimension.update({dimension_number:temp_disease})
            else:
            # further layers
                # nodes of this layer
                trajectory_disease = [x for item in temp_disease for x in trajectory_df.loc[trajectory_df[source]==item][target]]
                # break label
                temp_disease = trajectory_disease.copy()
                if len(temp_disease) == 0: break
                # index of the edges in DataFrame
                index = [x for item in temp_disease for x in trajectory_df.loc[trajectory_df[source]==item].index]
                trajectory_df.drop(index)
                # non-repeated nodes in this layer
                trajectory_disease = list(set(trajectory_disease))
                # save in the dict
                trajectory_dimension.update({dimension_number:trajectory_disease})
        # private attribute
        self.trajectory_dimension = trajectory_dimension
        return self.trajectory_dimension
    
    def commorbidity(self, source='phecode_d1', 
                     target='phecode_d2'):
        """_summary_

        Args:
            source (str, optional): _description_. Defaults to 'phecode_d1'.
            target (str, optional): _description_. Defaults to 'phecode_d2'.

        Returns:
            _type_: _description_
        """
        if not hasattr(self, "trajectory_dimension"):
            self.trajectory()
        # nodes only exist the commorbidity network
        nodes_of_commorbidity_except_trajectory = [x for x in self.__commorbidity_nodes
                                                   if x not in self.__trajectory_nodes]
        # hash_dict records the times of commorbidity nodes connecting nodes of trajectory network
        hash_dict = {x:[0]*len(self.__trajectory_nodes) for x in 
                     nodes_of_commorbidity_except_trajectory}
        # recording the times
        for i, value in enumerate(self.__trajectory_nodes):
            for item in self.__commorbidity_network_result.loc[self.__commorbidity_network_result[source]==value][target]:
                if item in nodes_of_commorbidity_except_trajectory:
                    hash_dict[item][i] += 1
            for item in self.__commorbidity_network_result.loc[self.__commorbidity_network_result[target]==value][source]:
                if item in nodes_of_commorbidity_except_trajectory:
                    hash_dict[item][i] += 1
        # the dimension number of nodes only in commorbidity network and connecting the trajectory nodes
        commorbidity_dimension_dict = {}
        for node in nodes_of_commorbidity_except_trajectory:
            indexs = np.flatnonzero(hash_dict[node])
            try:
                nodes_dimension = [self.__get_dimension(self.trajectory_dimension, 
                                                      self.__trajectory_nodes[x]) for x in indexs]
                if nodes_dimension:
                    dimension_number = stats.mode(nodes_dimension, keepdims=False)[0]
                    commorbidity_dimension_dict.update({node:dimension_number})
                else: continue
            except: continue
        # the dimension number of nodes remain
        while len(hash_dict) != len(commorbidity_dimension_dict):
            other_nodes = [x for x in hash_dict.keys() if x not in commorbidity_dimension_dict]
            for node in other_nodes:
                try:
                    nodes_dimension = self.__except_dimension(list(self.__commorbidity_network_result.loc[self.__commorbidity_network_result[source]==node][target]), 
                                                                 commorbidity_dimension_dict)
                    nodes_dimension += self.__except_dimension(list(self.__commorbidity_network_result.loc[self.__commorbidity_network_result[target]==node][source]), 
                                                                 commorbidity_dimension_dict)
                    dimension_number = stats.mode(nodes_dimension, keepdims=False)[0]
                    commorbidity_dimension_dict.update({node:dimension_number})
                except: continue
        # make cluster:dimension dict
        hash_dict = {n:[] for n in range(self.maximum_class_number)}
        for key, value in commorbidity_dimension_dict.items():
            if np.isnan(value): continue
            hash_dict[self.final_cluster[key]].append(int(value))
        for key, value in hash_dict.items():
            dimension = stats.mode(value, keepdims=False)[0]
            if np.isnan(dimension):
                dimension = max(self.trajectory().keys())
            hash_dict[key] = dimension
        for phecode in commorbidity_dimension_dict.keys():
            if np.isnan(commorbidity_dimension_dict[phecode]):
                cluster_number = self.final_cluster[phecode]
                commorbidity_dimension_dict.update({phecode:hash_dict[cluster_number]})
        self.commorbidity_dimension = commorbidity_dimension_dict
        return self.commorbidity_dimension
    
    def _dimension(self):
        """_summary_

        Returns:
            _type_: _description_
        """
        if not hasattr(self, "commorbidity_dimension"):
            self.commorbidity()
        # all nodes dimension number dict
        all_nodes_dimension = copy.deepcopy(self.trajectory_dimension)
        for node, dimension_number in self.commorbidity_dimension.items():
            if node not in all_nodes_dimension[dimension_number]:
                all_nodes_dimension[dimension_number] += [node]
        # remove the nodes of before layers
        for dimension_number in range(len(all_nodes_dimension)):
            nodes = all_nodes_dimension[dimension_number]
            temp_nodes = []
            for afterward_dimension in range(dimension_number+1, len(all_nodes_dimension)):
                temp_nodes += all_nodes_dimension[afterward_dimension]
            remove_nodes = []
            for node in nodes:
                if node in temp_nodes: remove_nodes.append(node)
            for node in remove_nodes:
                all_nodes_dimension[dimension_number].remove(node)
        self.dimension = all_nodes_dimension
        print("There are total %i layers" % (len(self.dimension)))
        self.end_height = int(input("input the height of the lowest layer:"))
        self.height = np.linspace(self.__exposure_disease_location[-1]-20, self.end_height, len(self.dimension))
        return self.dimension

    def color(self, code='phecode',
              cluster_name='category', 
              describe='phenotype'):
        """_summary_

        Args:
            code (str, optional): _description_. Defaults to 'phecode'.
            cluster_name (str, optional): _description_. Defaults to 'category'.
            describe (str, optional): _description_. Defaults to 'phenotype'.

        Returns:
            _type_: _description_
        """
        code_system = dict(zip(self.__describe[code], self.__describe[cluster_name]))
        system_color = dict(zip(self.__system, self.__system_color))
        code_color = {node:system_color.get(code_system.get(node)) 
                      for node in code_system.keys()}
        self.color_map = code_color
        self.__code_system = code_system
        self.__system_color = system_color
        self.__code_name = {cd:self.__split_name('%s (%.1f)' % (words, cd))
                            for cd, words in self.__describe[[code, describe]].values}
        self.code_name = {cd:self.__split_name('%s (%.1f)' % (words, cd))
                            for cd, words in self.__describe[[code, describe]].values}
        self.__clusterNumber_system = {}
        self.__clusterNumber_color = {}
        from collections import Counter
        def most_frequent_element(lst):
            counter = Counter(lst)
            most_common = counter.most_common(1)
            return most_common[0][0]
        for cluster_number in set(self.final_cluster.values()):
            disease_list = []
            for key, value in self.final_cluster.items():
                if value == cluster_number:
                    disease_list.append(self.__code_system[key])
            self.__clusterNumber_system[cluster_number] = most_frequent_element(disease_list)
            self.__clusterNumber_color[cluster_number] = self.__system_color[most_frequent_element(disease_list)]
        return self.color_map

    def incluster_commorbidity(self, incluster:list[str]) -> list[str]:
        for _, row in self.__commorbidity_network_result.iterrows():
            if row["phecode_d1"] in incluster and \
            self.final_cluster[row["phecode_d2"]]==self.final_cluster[row["phecode_d1"]]:
                incluster.append(row["phecode_d2"])
            if row["phecode_d2"] in incluster and \
            self.final_cluster[row["phecode_d2"]]==self.final_cluster[row["phecode_d1"]]:
                incluster.append(row["phecode_d1"])
        return list(set(incluster))

    def __full_plot(self, line_color:str, 
                    line_width:float,
                    source='phecode_d1', 
                    target='phecode_d2'):
        plot_data = []
        for system in self.__system:
            system_nodes = [x for x in self.__commorbidity_nodes if self.__code_system[x]==system]
            if system_nodes:
                for node in system_nodes:
                    x_, y_, z_, colorscale_, light_dict, _, name, light_position = self.plotly_ball(list(self.location_dict[node]),
                                                                                                        self.get_node_size(node)/13,
                                                                                                        str(self.__code_name[node]),
                                                                                                        'Disease',
                                                                                                        self.__system_color[self.__code_system[node]])
                    if node == system_nodes[0]:
                        data = go.Surface(x=x_, y=y_, z=z_, colorscale=colorscale_, showlegend=True,
                                          lighting=light_dict, hovertemplate=name, name='%s Disease' % (system.title()), showscale=False,
                                          legendgroup='ball', legendgrouptitle_text='Diseases',lightposition=light_position)
                    else:
                        data = go.Surface(x=x_, y=y_, z=z_, colorscale=colorscale_, showlegend=False,
                                          lighting=light_dict, hovertemplate=name, name='%s Disease' % (system.title()), showscale=False,
                                          legendgroup='ball', legendgrouptitle_text='Diseases',lightposition=light_position)
                    plot_data += [data]
        # plot edges
        graph = nx.DiGraph()
        if self.__exposure_disease == 9999:
            [graph.add_edge(row[source], row[target]) for _, row in self.__disease_trajectory.iterrows() if row[source]!=9999]
        else:
            [graph.add_edge(row[source], row[target]) for _, row in self.__disease_trajectory.iterrows()]
        edges_x, edges_y, edges_z = [], [], []
        for edge in graph.edges():
            edges_x += [self.location_dict[edge[0]][0], self.location_dict[edge[1]][0], None]
            edges_y += [self.location_dict[edge[0]][1], self.location_dict[edge[1]][1], None]
            edges_z += [self.location_dict[edge[0]][2], self.location_dict[edge[1]][2], None]

        trace_edges = go.Scatter3d(x=edges_x,
                                   y=edges_y,
                                   z=edges_z,
                                   line=dict(color=line_color,width=line_width),
                                   mode='lines',
                                   hoverinfo='none',
                                   legendgroup='none',
                                   legendgrouptitle_text='All Trajectories',
                                   name='Trajectories',
                                   showlegend=True)
        plot_data += [trace_edges]
        return plot_data
    
    def __half_plot(self, main_line_width:float, 
                    nonMain_line_color='silver', 
                    nonMain_line_width=1,
                    source='phecode_d1', 
                    target='phecode_d2'):
        plot_data = []
        # plot nodes (incluster)
        self.incluster_nodes = []
        for cluster_number in set(self.final_cluster.values()):
            cluster_nodes = []
            source_list = [self.__exposure_disease]
            while True:
                begin_number_nodes = len(cluster_nodes)
                target_list = []
                for node in source_list:
                    target_list += [x for x in self.__disease_trajectory[self.__disease_trajectory[source]\
                                    ==node][target].values if self.final_cluster[x]==cluster_number]
                cluster_nodes += target_list
                cluster_nodes = list(set(cluster_nodes))
                end_number_nodes = len(cluster_nodes)
                if end_number_nodes == begin_number_nodes:
                    break
                source_list = target_list
            self.incluster_nodes += cluster_nodes

        for system in self.__system:
            system_nodes = [x for x in self.incluster_nodes if self.__code_system[x]==system]
            if system_nodes:
                for node in system_nodes:
                    x_, y_, z_, colorscale_, light_dict, _, name, light_position = self.plotly_ball(list(self.location_dict[node]),
                                                                                                        self.get_node_size(node)/13,
                                                                                                        str(self.__code_name[node]),
                                                                                                        'Disease',
                                                                                                        self.__system_color[self.__code_system[node]])
                    if node == system_nodes[0]:
                        data = go.Surface(x=x_, y=y_, z=z_, colorscale=colorscale_, showlegend=True,
                                          lighting=light_dict, hovertemplate=name, name='%s Disease' % (system.title()), showscale=False,
                                          legendgroup='ball', legendgrouptitle_text='Diseases', lightposition=light_position)
                    else:
                        data = go.Surface(x=x_, y=y_, z=z_, colorscale=colorscale_, showlegend=False,
                                          lighting=light_dict, hovertemplate=name, name='%s Disease' % (system.title()), showscale=False,
                                          legendgroup='ball', legendgrouptitle_text='Diseases', lightposition=light_position)
                    plot_data += [data]
        # plot nodes (outcluster)
        self.outcluster_nodes = [x for x in self.__commorbidity_nodes if x not in self.incluster_nodes]
        for node in self.outcluster_nodes:
            x_, y_, z_, colorscale_, light_dict, _, name, light_position = self.plotly_ball(list(self.location_dict[node]),
                                                                                                self.get_node_size(node)/13,
                                                                                                str(self.__code_name[node]),
                                                                                                'Disease',
                                                                                                'grey')
            if node == self.outcluster_nodes[0]:
                data = go.Surface(x=x_, y=y_, z=z_, colorscale=colorscale_, showlegend=True,
                                  lighting=light_dict, hovertemplate=name, name='Outer disease', showscale=False,
                                  legendgroup='ball', legendgrouptitle_text='Diseases', lightposition=light_position)
            else:
                data = go.Surface(x=x_, y=y_, z=z_, colorscale=colorscale_, showlegend=False,
                                  lighting=light_dict, hovertemplate=name, name='Outer disease', showscale=False,
                                  legendgroup='ball', legendgrouptitle_text='Diseases', lightposition=light_position)
            plot_data += [data]
        # plot edges (incluster)
        incluster_trajectory = []
        for system in set(self.final_cluster.values()):
            graph = nx.DiGraph()
            system_nodes = []
            source_list = [self.__exposure_disease]
            while True:
                begin_number_nodes = len(system_nodes)
                target_list = []
                for node in source_list:
                    template_target_list = [x for x in self.__disease_trajectory[self.__disease_trajectory[source]\
                                            ==node][target].values if self.final_cluster[x]==system]
                    target_list += template_target_list
                    if template_target_list:
                        [graph.add_edge(node, target) for target in template_target_list if node!=9999]
                        [incluster_trajectory.append([node, target]) for target in template_target_list]
                system_nodes += target_list
                system_nodes = list(set(system_nodes))
                end_number_nodes = len(system_nodes)
                if end_number_nodes == begin_number_nodes: break
                source_list = target_list
            edges_color, edges_x, edges_y, edges_z = [], [], [], []
            for edge in graph.edges():
                edges_color.append(self.__clusterNumber_color[system])
                edges_x += [self.location_dict[edge[0]][0], self.location_dict[edge[1]][0], None]
                edges_y += [self.location_dict[edge[0]][1], self.location_dict[edge[1]][1], None]
                edges_z += [self.location_dict[edge[0]][2], self.location_dict[edge[1]][2], None]
            trace_edges = go.Scatter3d(x=edges_x,
                                       y=edges_y, 
                                       z=edges_z, 
                                       line=dict(color=edges_color*3,width=main_line_width),
                                       mode='lines',
                                       hoverinfo='none',
                                       legendgroup='1',
                                       legendgrouptitle_text='All Trajectories',
                                       name='Cluster %i Trajectories' % (system+1))
            plot_data += [trace_edges]
        # plot edges (outcluster)
        graph = nx.DiGraph()
        self.incluster_trajectory = []
        for disease_pairs in incluster_trajectory:
            if disease_pairs not in self.incluster_trajectory:
                self.incluster_trajectory.append(disease_pairs)
        [graph.add_edge(disease_pairs[0], disease_pairs[1]) for disease_pairs 
         in self.__trajectory_disease_pairs if (disease_pairs not in incluster_trajectory) & (disease_pairs[0]!=9999)]
        edges_x, edges_y, edges_z = [], [], []
        for edge in graph.edges():
            edges_x += [self.location_dict[edge[0]][0], self.location_dict[edge[1]][0], None]
            edges_y += [self.location_dict[edge[0]][1], self.location_dict[edge[1]][1], None]
            edges_z += [self.location_dict[edge[0]][2], self.location_dict[edge[1]][2], None]
        trace_edges = go.Scatter3d(x=edges_x,
                                   y=edges_y, 
                                   z=edges_z, 
                                   line=dict(color=nonMain_line_color,width=nonMain_line_width),
                                   mode='lines',
                                   hoverinfo='none',
                                   legendgroup='1',
                                   legendgrouptitle_text='All Trajectories',
                                   name="Insignificant Trajectories")
        plot_data += [trace_edges]
        return plot_data

    def __compact_plot(self,main_line_width:float,
                       source='phecode_d1',
                       target='phecode_d2'):
        plot_data = []
        # plot nodes (incluster)
        self.incluster_nodes = []
        for cluster_number in set(self.final_cluster.values()):
            cluster_nodes = []
            source_list = [self.__exposure_disease]
            while True:
                begin_number_nodes = len(cluster_nodes)
                target_list = []
                for node in source_list:
                    target_list += [x for x in self.__disease_trajectory[self.__disease_trajectory[source]\
                                    ==node][target].values if self.final_cluster[x]==cluster_number]
                cluster_nodes += target_list
                cluster_nodes = list(set(cluster_nodes))
                end_number_nodes = len(cluster_nodes)
                if end_number_nodes == begin_number_nodes: break
                source_list = target_list
            self.incluster_nodes += cluster_nodes
        compact_node = []
        for system in self.__system:
            system_nodes = [x for x in self.incluster_nodes if self.__code_system[x]==system]
            if system_nodes:
                for node in system_nodes:
                    compact_node.append(node)
                    x_, y_, z_, colorscale_, light_dict, _, name, light_position = self.plotly_ball(list(self.location_dict[node]),
                                                                                                        self.get_node_size(node)/13,
                                                                                                        str(self.__code_name[node]),
                                                                                                        'Disease',
                                                                                                        self.__system_color[self.__code_system[node]])
                    if node == system_nodes[0]:
                        data = go.Surface(x=x_, y=y_, z=z_, colorscale=colorscale_, showlegend=True,
                                    lighting=light_dict, hovertemplate=name, name='%s Disease' % (system.title()), showscale=False,
                                    legendgroup='ball', legendgrouptitle_text='Diseases', lightposition=light_position)
                    else:
                        data = go.Surface(x=x_, y=y_, z=z_, colorscale=colorscale_, showlegend=False,
                                          lighting=light_dict, hovertemplate=name, name='%s Disease' % (system.title()), showscale=False,
                                          legendgroup='ball', legendgrouptitle_text='Diseases', lightposition=light_position)
                    plot_data += [data]
        self.compact_node = compact_node
        # plot edges (incluster)
        for system in set(self.final_cluster.values()):
            graph = nx.DiGraph()
            system_nodes = []
            source_list = [self.__exposure_disease]
            while True:
                begin_number_nodes = len(system_nodes)
                target_list = []
                for node in source_list:
                    template_target_list = [x for x in self.__disease_trajectory[self.__disease_trajectory[source]\
                                            ==node][target].values if self.final_cluster[x]==system]
                    target_list += template_target_list
                    if template_target_list:
                        [graph.add_edge(node, target) for target in template_target_list if node != 9999]
                system_nodes += target_list
                system_nodes = list(set(system_nodes))
                end_number_nodes = len(system_nodes)
                if end_number_nodes == begin_number_nodes:break
                source_list = target_list
            edges_color, edges_x, edges_y, edges_z = [], [], [], []
            for edge in graph.edges():
                edges_color.append(self.__clusterNumber_color[system])
                edges_x += [self.location_dict[edge[0]][0], self.location_dict[edge[1]][0], None]
                edges_y += [self.location_dict[edge[0]][1], self.location_dict[edge[1]][1], None]
                edges_z += [self.location_dict[edge[0]][2], self.location_dict[edge[1]][2], None]
            trace_edges = go.Scatter3d(x=edges_x,
                                       y=edges_y, 
                                       z=edges_z, 
                                       line=dict(color=edges_color*3, width=main_line_width),
                                       mode='lines',
                                       hoverinfo='none',
                                       legendgroup='1',
                                       legendgrouptitle_text='All Trajectories',
                                       name='Cluster %i Trajectory' % (system+1))
            plot_data += [trace_edges]
        return plot_data

    def plot_3d(self, max_radius:float, 
                min_radius:float,
                plot_method:str,
                line_color:str, 
                line_width:float, 
                file_name:str,
                layout_width:float=900,
                layout_height:float=900,
                font_style:str='Times New Roman',
                font_size:float=15,
                location_method='random'):
        
        if not hasattr(self, "final_cluster"):
            self.cluster()
        if not hasattr(self, "dimension"):
            self._dimension()
        if not hasattr(self, 'location_dict'):
            if location_method == 'random':
                self.location_random(max_radius, min_radius)
        if not hasattr(self, "color_map"):
            self.color()
        plot_data = []
        # plot the origin disease
        if self.__exposure_disease != 9999:
            origin_data = go.Scatter3d(x=[self.__exposure_disease_location[0]],
                                    y=[self.__exposure_disease_location[1]],
                                    z=[self.__exposure_disease_location[2]],
                                    mode='markers',
                                    marker=dict(symbol='circle',size=self.__exposure_disease_size,color='black'),
                                    text=['Depression'],
                                    hoverinfo='text',
                                    legendgroup='origin',
                                    legendgrouptitle_text='Origin of Trajectories',
                                    name='Depression',
                                    showlegend=True)
            plot_data += [origin_data]
        # plot the nodes and edges
        if plot_method == 'full':
            plot_data += self.__full_plot(line_color, line_width)
        elif plot_method == 'compact':
            plot_data += self.__compact_plot(line_width)
        elif plot_method == 'half':
            plot_data += self.__half_plot(line_width)
        else:
            raise KeyError("the plot_method is not exist")
        # axis
        axis = dict(showbackground=False, 
                    showline=False,
                    zeroline = False,
                    showgrid = False,
                    showticklabels = False,
                    title='')
        # layout
        layout = go.Layout(title=dict(text="", 
                                      font=dict(size=30,family=font_style),
                                      x=0.45),
                            width=layout_width,
                            height=layout_height,
                            showlegend=True,
                            scene=dict(xaxis=dict(axis),
                                       yaxis=dict(axis),
                                       zaxis=dict(axis)),
                            margin=dict(t=100),
                            hovermode='closest', 
                            legend=dict(title=dict(text='Trace of clusters'),
                            font=dict(family=font_style,size=font_size),
                            itemclick=False), 
                            font=dict(family=font_style))
        # plot the figure
        fig = go.Figure(data=plot_data, layout=layout)
        # create the file of the figure
        py.plot(fig, filename=file_name)

    def comorbidity_network_plot(self, max_radius:float,
                                 min_radius:float,
                                 line_width:float=1,
                                 source:str="phecode_d1",
                                 target:str="phecode_d2",
                                 location_method:str="random",
                                 line_color:str="black"
                                 ) -> None:

        if not hasattr(self, "final_cluster"):
            self.cluster()
        if not hasattr(self, "dimension"):
            self._dimension()
        if not hasattr(self, 'location_dict'):
            if location_method == 'random':
                self.location_random(max_radius, min_radius)
        if not hasattr(self, "color_map"):
            self.color()

        fig = go.Figure()

        for _, row in self.__commorbidity_network_result.iterrows():
            x_values = [self.location_dict[row[source]][0], self.location_dict[row[target]][0]]
            y_values = [self.location_dict[row[source]][1], self.location_dict[row[target]][1]]
            fig.add_trace(go.Scatter(x=x_values, 
                                     y=y_values, 
                                     mode='lines', 
                                     line=dict(color=line_color, width=line_width)))

        for system in self.__system:
            system_nodes = [x for x in self.__commorbidity_nodes if self.__code_system[x]==system]
            if system_nodes:
                for node in system_nodes:
                    theta = np.linspace(0, 2*np.pi, 100)
                    x_circle = self.location_dict[node][0] + self.get_node_size(node)/13 * np.cos(theta)
                    y_circle = self.location_dict[node][1] + self.get_node_size(node)/13 * np.cos(theta)

                    if node == system_nodes[0]:
                        fig.add_trace(go.Scatter(x=x_circle, y=y_circle,
                                                fill='toself',
                                                fillcolor=self.__system_color[self.__code_system[node]],
                                                showlegend=True,
                                                hovertemplate=self.__code_name[node],
                                                name='%s Disease' % (system.title()), 
                                                showscale=False,
                                                legendgroup='ball', 
                                                legendgrouptitle_text='Diseases',
                                                line=dict(color=self.__system_color[self.__code_system[node]], width=1)))
                    else:
                        fig.add_trace(go.Scatter(x=x_circle, y=y_circle,
                                                fill='toself',
                                                fillcolor=self.__system_color[self.__code_system[node]],
                                                showlegend=False,
                                                hovertemplate=self.__code_name[node],
                                                name='%s Disease' % (system.title()), 
                                                showscale=False,
                                                legendgroup='ball', 
                                                legendgrouptitle_text='Diseases',
                                                line=dict(color=self.__system_color[self.__code_system[node]], width=1)))
        fig.update_layout(
            showlegend=False,
            xaxis=dict(visible=False, scaleanchor="y"),  # 确保x、y比例一致
            yaxis=dict(visible=False),
            plot_bgcolor='white',
            margin=dict(l=10, r=10, t=10, b=10)  # 去除边距
        )

        fig.show()

    def incluster_trajectory_plot(self, distance:float, 
                                  layer_distance:float,
                                  line_width:float,
                                  line_color:float,
                                  max_radius:float,
                                  min_radius:float,
                                  location_method:str="random",
                                  source:str="phecode_d1",
                                  target:str="phecode_d2") -> None:
        
        if not hasattr(self, "final_cluster"):
            self.cluster()
        if not hasattr(self, "dimension"):
            self._dimension()
        if not hasattr(self, 'location_dict'):
            if location_method == 'random':
                self.location_random(max_radius, min_radius)
        if not hasattr(self, "color_map"):
            self.color()
        # nodes incluster
        cluster_nodes_dict = {}
        for cluster_number in set(self.final_cluster.values()):
            cluster_nodes = []
            source_list = [self.__exposure_disease]
            while True:
                begin_number_nodes = len(cluster_nodes)
                target_list = []
                for node in source_list:
                    target_list += [x for x in self.__disease_trajectory[self.__disease_trajectory[source]\
                                    ==node][target].values if self.final_cluster[x]==cluster_number]
                cluster_nodes += target_list
                cluster_nodes = list(set(cluster_nodes))
                end_number_nodes = len(cluster_nodes)
                if end_number_nodes == begin_number_nodes: break
                source_list = target_list
            cluster_nodes_dict[cluster_number] = cluster_nodes

        def calculate_circle_positions_with_ids(radii, spacing):
            """
            计算圆心的 x 坐标，使得圆均匀分布在 y 轴两侧，并且相互之间相离，保持固定的外切距离。

            参数:
                radii (list of float): 圆的半径列表，按给定顺序排列。
                spacing (float): 圆与圆之间的额外距离（即外切距离）。

            返回:
                dict: {编号: x轴坐标}，编号从 1 开始，与输入顺序对应。
            """
            circle_positions = {}
            position = 0
            for i, radius in enumerate(radii):
                circle_id = i + 1
                if i == 0:
                    circle_positions[circle_id] = position
                    position += 2 * radius + spacing  
                else:
                    if i % 2 == 1:
                        circle_positions[circle_id] = position  
                    else:
                        circle_positions[circle_id] = -position  
                    position += 2 * radius + spacing 
            return circle_positions

        incluster_location = {}
        for cluster_number, nodes_lst in cluster_nodes_dict.items():
            temp_dimension = {x:[] for x in self.dimension.keys()}
            for node in nodes_lst:
                temp_dimension[self.commorbidity_dimension[node]].append(node)
            for key, value in temp_dimension.items():
                if value: x_value = calculate_circle_positions_with_ids([math.pow(self.get_node_size(x) / 4 * 3 / math.pi, 1/3) for x in value], distance)
                for i, j in x_value.items():
                    incluster_location[value[i-1]] = (j, self.__exposure_disease_location[1]-key*layer_distance)
        
        for cluster_number, node_lst in cluster_nodes_dict.items():
            edges = [(x, y) for x,y in dict(self.__disease_trajectory[source], self.__disease_trajectory[target]).items() if x in node_lst and y in node_lst]
            G = nx.DiGraph()
            G.add_edges_from(edges)
            pos = {x:incluster_location[x] for x in node_lst}
            color_values = [self.__system_color[self.__code_system[n]] for n in G.nodes()]
            size_values = [self.get_node_size[n] for n in G.nodes()]
            label_values = {n: self.__code_name[n] for n in G.nodes()}
            plt.figure(figsize=(8, 8))
            nx.draw_networkx(
                G,
                pos,
                with_labels=True,
                labels=label_values,
                node_color=color_values,
                node_size=size_values,
                font_size=12, 
                font_weight='bold',
                edge_color=line_color,
                arrows=True,
                arrowsize=20,
            )
            # 关闭坐标轴
            plt.title("dad")
            plt.axis("off")
            plt.show()