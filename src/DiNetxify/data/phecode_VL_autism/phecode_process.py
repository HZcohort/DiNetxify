# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 16:32:53 2024

@author: Can Hou
"""

"""
This is a variation from the 1.3a version for a specific use case in the DiNetxify project.
It is not intended to be used as a general-purpose phecode version.
The main change is that 313.3 is not considered as leaf node for 313 in level 1 of phecode hierarchy.
"""

# %%phecode data prepare
import pandas as pd
import numpy as np

level1_info = np.load(r'src\DiNetxify\data\phecode_1.3a\level1_info.npy', allow_pickle=True).item()
level1_info[313]['leaf_list'].remove(313.3)
level1_info[313]['exclude_list'] = [[312.0, 312.3],
                                    [313.0, 313.1, 313.2],
                                    [315.2, 315.0, 315.1, 315.3]]

np.save(r'src\DiNetxify\data\phecode_VL_autism\level1_info.npy', level1_info)