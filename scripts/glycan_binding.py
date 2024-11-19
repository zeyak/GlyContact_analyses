import pandas as pd
from glycowork.motif.analysis import get_heatmap
import matplotlib
matplotlib.use('TkAgg')  # Use the TkAgg backend for rendering plots


# Load glycan binding data
glycan_binding = pd.read_csv('/Users/xakdze/PycharmProjects/GlycoShape/data/glycan_binding.csv')
print(glycan_binding.head())

get_heatmap(glycan_binding.iloc[:,:-2], motifs = True, feature_set = ['exhaustive'], datatype = 'response', yticklabels = 1,
            xticklabels = False)