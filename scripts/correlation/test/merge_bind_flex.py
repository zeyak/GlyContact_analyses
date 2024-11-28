import pickle
import pandas as pd
from glycowork.motif.analysis import get_heatmap
from glycowork.glycan_data.loader import glycan_binding

from scripts.correlation.merged.merge_bind_flex_SNA import binding_df

"""
Merged binding and flexibility data in a heatmap:

binding_df = 2745 glycans 1465 proteins
flex_df = 582 glycans
merged_df_glycans = 198 glycans 1465 proteins

the plot is still massive!
"""



# Define file paths
binding_data_path = '/data/glycan_binding.csv'
flex_data_path = '/data/glycan_graphs.pkl'

# Load glycan flexibility data from the pickle file
with open(flex_data_path, 'rb') as file:
    flex_data = pickle.load(file)
# Load binding data into a DataFrame and transpose it
binding_df = pd.read_csv(binding_data_path)





binding_df_T = binding_df.T.reset_index()


# Convert the flexibility data (dictionary) into a DataFrame
flex_df = pd.DataFrame.from_dict(flex_data, orient='index')
# Reset the index and keep only the glycan names
flex_glycan_df = flex_df.reset_index()[['index']]

# Merge binding and flexibility data on the index
merged_df_glycans = pd.merge(binding_df_T, flex_glycan_df, on="index")

merged_df_glycans_ = merged_df_glycans.T
# Set the first row as the column headers
merged_df_glycans_.columns = merged_df_glycans_.iloc[0]
# Drop the first row since it's now the header
merged_df_glycans_ = merged_df_glycans_.iloc[1:]

# Generate the heatmap
get_heatmap(
    merged_df_glycans_,
    #binding_df.iloc[:,:-2],
    #glycan_binding.iloc[:, :-2],
    motifs=True,
    feature_set=['exhaustive'],
    datatype='response',
    yticklabels=1,
    xticklabels=False
)


