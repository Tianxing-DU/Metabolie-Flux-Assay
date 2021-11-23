import pandas as pd
from numpy import NaN
import pandas as pd, matplotlib.pyplot as plt, seaborn as sns, numpy as np
from sklearn import preprocessing
from scipy.stats import sem
from matplotlib.ticker import PercentFormatter

list_custom = ['C18:2', 'C18:2-OH', 'C18:2-DC', 'C18:1', 'C18:1-OH', 'C18:1-DC', 'C18:0', 'C18:0-OH', 'C18:0-DC',
               'C16:2', 'C16:2-DC', 'C16:1', 'C16:1-OH', 'C16:1-DC', 'C16:0', 'C16:0-OH', 'C16:0-DC', 'C14:2', 'C14:1',
               'C14:0', 'C14:0-OH', 'C12:2', 'C12:1', 'C12:0', 'C12:0-DC', 'C10:2', 'C10:1', 'C10:0', 'C8:1', 'C8:0',
               'C6:0', 'C6:0-OH', 'C6:0-DC', 'C5:1', 'C5:0', 'C5:0-OH', 'C5:0-DC', 'C4:0', 'C4:0-OH', 'C4:0-DC', 'C3:0',
               'C3:0-DC', 'C2:0']

long_chain = ['C18:2', 'C18:2-OH', 'C18:2-DC', 'C18:1', 'C18:1-OH', 'C18:1-DC', 'C18:0', 'C18:0-OH', 'C18:0-DC',
              'C16:2', 'C16:2-DC', 'C16:1', 'C16:1-OH', 'C16:1-DC', 'C16:0', 'C16:0-OH', 'C16:0-DC', 'C14:2', 'C14:1',
              'C14:0', 'C14:0-OH']
middle_chain = ['C12:2', 'C12:1', 'C12:0', 'C12:0-DC', 'C10:2', 'C10:1', 'C10:0', 'C8:1', 'C8:0']
short_chain = ['C6:0', 'C6:0-OH', 'C6:0-DC', 'C5:1', 'C5:0', 'C5:0-OH', 'C5:0-DC', 'C4:0', 'C4:0-OH', 'C4:0-DC', 'C3:0',
               'C3:0-DC', 'C2:0']

LOD_name = ["C0:0", "C2:0", "C3:0", "C3:0-DC", "C4:0", "C4:0-DC", "C4:0-OH", "C5:0", "C5:0-DC", "C5:0-OH", "C5:1",
            "C6:0", "C6:0-DC", "C6:0-OH",
            "C8:0", "C8:1", "C10:0", "C10:1", "C10:2", "C12:0", "C12:0-DC", "C12:1", "C12:2",
            "C14:0", "C14:0-OH", "C14:1", "C14:2", "C16:0", "C16:0-DC", "C16:0-OH", "C16:1", "C16:1-DC", "C16:1-OH",
            "C16:2", "C16:2-DC", "C18:0", "C18:0-DC", "C18:0-OH", "C18:1", "C18:1-DC", "C18:1-OH", "C18:2",
            "C18:2-DC", "C18:2-OH"]

LOD_value = [NaN, 27.24, 15.31, 12.48, 16.68, 2.90, 3.48, 20.11, 1.58, 4.85, 2.33,
             10.61, 17.64, 7.99,
             6.50, 1.63, 2.90, 4.66, 0.79, 3.48, 0.22, 0.12, 0.19,
             0.41, 0.24, 0.19, 0.29, 0.86, 0.12, 0.60, 0.12, 0.02, 0.17,
             0.05, 0.02, 0.36, 0.07, 0.05, 0.07, 0.05, 0.05, 0.07,
             0.05, 0.31]  # comes from R. Ensenauer et al. 2012

LOD_value = [element * 0.001 for element in LOD_value]  # unify the unit

df_LOD = pd.DataFrame(list(zip(LOD_name, LOD_value)), columns=['C-Chains', 'LOD'])                   
# alternative way: df_LOD = pd.DataFrame({'C-Chains': LOD_name, 'LOD': LOD_value})

dict_LOD = dict(zip(LOD_name, LOD_value))  # creat dict for later selection in main codes
