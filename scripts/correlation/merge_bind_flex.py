import pickle
import pandas as pd
from glycowork.motif.analysis import get_heatmap
import matplotlib.pyplot as plt

"""
Merged binding and flexibility data in a heatmap
Plot only a slice of the heatmap
because the full heatmap is too large
"""

# Define file paths
binding_data_path = '/Users/xakdze/PycharmProjects/GlycoShape/data/glycan_binding.csv'
flex_data_path = '/Users/xakdze/PycharmProjects/GlycoShape/data/glycan_graphs.pkl'

# Load glycan flexibility data from the pickle file
with open(flex_data_path, 'rb') as file:
    flex_data = pickle.load(file)
# Load binding data into a DataFrame and transpose it
binding_df = pd.read_csv(binding_data_path)
binding_df_T = binding_df.T

# Convert the flexibility data (dictionary) into a DataFrame
flex_df = pd.DataFrame.from_dict(flex_data, orient='index')
# Reset the index and keep only the glycan names
flex_glycan_df = flex_df.reset_index()[['index']]

# Merge binding and flexibility data on the index
merged_df_glycans = pd.merge(binding_df_T, flex_df, left_index=True, right_index=True).reset_index()[['index']]
# Prepare data for the heatmap by merging with binding data (excluding last two columns)
merged_heatmap = pd.merge(
    merged_df_glycans,
    binding_df_T.iloc[:, :-2],
    left_on='index',
    right_index=True
)

merged_slice = merged_heatmap.iloc[:10, :10]

#Orientation doesnt matter for the cluster heatmap
"""# transpose merged_heatmap and make the first row the header
merged_heatmap = merged_heatmap.T
merged_heatmap.columns = merged_heatmap.iloc[0]
merged_heatmap = merged_heatmap.drop('index')
"""

# Save the heatmap to a file
output_path = '/Users/xakdze/PycharmProjects/GlycoShape/out/merged_heatmap.png'
plt.savefig(output_path)

# Generate the heatmap
get_heatmap(
    merged_heatmap.iloc[:10, :10],
    motifs=True,
    feature_set=['exhaustive'],
    datatype='response',
    yticklabels=1,
    xticklabels=False
)


