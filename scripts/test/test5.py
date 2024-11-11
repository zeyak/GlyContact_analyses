import pickle
import pandas as pd
import networkx as nx

# Load the data from the pickle file
graph_file_path = '/data/glycan_graphs.pkl'
glycan_binding_file_path = '/data/glycan_binding.csv'

with open(graph_file_path, 'rb') as f:
    data = pickle.load(f)

# convret data to a DataFrame
glycan_df = pd.DataFrame.from_dict(data, orient='index')

# Load glycan binding data
glycan_binding = pd.read_csv(glycan_binding_file_path).T


#merge two dataframes from row index how inner
merged_df = pd.merge(glycan_df, glycan_binding, left_index=True, right_index=True)

