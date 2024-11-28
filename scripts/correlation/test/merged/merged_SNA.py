import pandas as pd
import pickle
from glycowork.motif.analysis import get_heatmap

from scripts.correlation.binding.binding_SNA import binding_df, binding_df_T

"""
This is not the correct way of filtering SNA
use merged instead of string matching
"""

# Specify file paths for the binding and flexibility data
binding_data_path = 'data/glycan_binding.csv'
flex_data_path = 'data/glycan_graphs.pkl'

# Load glycan flexibility data from the pickle file
with open(flex_data_path, 'rb') as file:
    flex_data = pickle.load(file)

# Load binding data into a DataFrame and transpose it
# Transpose the DataFrame and reset index
#binding_df_T = binding_df.T.reset_index()

# Dynamically create a DataFrame using the index and the last row
#last_row = binding_df_T.iloc[-1, :]  # Extract the last row
#gly_lectin = pd.concat([binding_df_T['index'], last_row], axis=1)


# Convert the flexibility data (a dictionary) into a DataFrame
flex_df = pd.DataFrame.from_dict(flex_data, orient='index')

# Reset the index of the flexibility DataFrame and retain only the glycan names
flex_glycan_df = flex_df.reset_index()[['index']]

# Merge the transposed binding data with the flexibility glycan DataFrame on the index
merged_id = pd.merge(
    binding_df_T,
    flex_glycan_df,
    left_index=True,
    right_on="index"
).reset_index()[['index']]

merged = pd.merge(
    merged_id, binding_df_T,
    left_on='index',
    right_index=True
)

merged.set_index('index', inplace=True)

# Convert the last row into a DataFrame and concatenate
new_row = pd.DataFrame(binding_df_T.iloc[-1, :]).T
merged_ = pd.concat([merged, new_row], axis=0).T


lectin ="SNA"


# look for string lectin on the last column of the merged dataframe(lectin flex & binding data)
lectin= merged_[merged_.iloc[:, -1].str.contains(lectin, na=False)]
lectin_binding_glycans= lectin.dropna(axis=1, how= "all") #drop the columns(glycans) that have NaN values

# setindex to protein names and reasssign it to a vsriable
#lectin = lectin_binding_glycans.set_index('protein')

#melt the dataframe
lectin_melt = lectin_binding_glycans.melt(id_vars = 'protein', var_name = 'glycan', value_name = 'binding')

#filter the values taht are more than 5
lectin_filt = lectin_melt[lectin_melt['binding'] > 1]


# Generate the heatmap
get_heatmap(
    lectin_binding_glycans.iloc[:,:-2],
    motifs=True,
    feature_set=['exhaustive'],
    datatype='response',
    yticklabels=1,
    xticklabels=False
)


