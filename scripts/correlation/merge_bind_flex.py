import pickle
import pandas as pd
from glycowork.motif.analysis import get_heatmap
import matplotlib
matplotlib.use('TkAgg')  # Use the TkAgg backend for rendering plots


#input
binding_data_path = '/Users/xakdze/PycharmProjects/GlycoShape/data/glycan_binding.csv'
flex_data_path = '/Users/xakdze/PycharmProjects/GlycoShape/data/glycan_graphs.pkl'

# Load glycan flexibility data from the pickle file
with open(flex_data_path, 'rb') as file:
    flex_data = pickle.load(file)

# Read files into dataframes
binding_df = pd.read_csv(binding_data_path).T

# Convert the loaded data into a DataFrame
flex_df = pd.DataFrame.from_dict(flex_data, orient='index')

# Merge the two dataframes
merged_df = pd.merge(binding_df, flex_df, left_index=True, right_index=True)
merged_heatmap= merged_df.iloc[:,:-2].T.reset_index(drop=True)

get_heatmap(merged_heatmap, motifs = True, feature_set = ['exhaustive'], datatype = 'response', yticklabels = 1,
            xticklabels = False)