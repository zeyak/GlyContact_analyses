import pandas as pd
from glycowork.motif.analysis import get_heatmap
import matplotlib.pyplot as plt


"""
Plots all glycan binding data in a heatmap
From gylcowork 
"""

# Load glycan binding data
glycan_binding = pd.read_csv('/data/glycan_binding.csv')
print(glycan_binding.head())

get_heatmap(glycan_binding.iloc[:,:-2], motifs = True, feature_set = ['exhaustive'], datatype = 'response', yticklabels = 1,
            xticklabels = False)