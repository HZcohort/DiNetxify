# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 23:48:05 2024

@author: Can Hou/Haowen Liu - Biomedical Big data center of West China Hospital, Sichuan University
"""

from .data_management import DiseaseNetworkData
from .analysis import (
    phewas,
    phewas_multipletests, 
    comorbidity_strength, 
    comorbidity_strength_multipletests,
    binomial_test, 
    binomial_multipletests, 
    comorbidity_network, 
    comorbidity_multipletests,
    disease_trajectory, 
    trajectory_multipletests
)
from .analysis_pipeline import disease_network_pipeline

__version__ = '0.1.7'
















