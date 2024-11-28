import pandas as pd
from glycowork.motif.analysis import get_heatmap

# Define file paths
binding_data_path = 'data/glycan_binding.csv'

# Load binding data into a DataFrame and transpose it
binding_df = pd.read_csv(binding_data_path)
binding_df_T = binding_df.T


lectin ="ConA"


# look for string lectin on the last column of the binding_df
lectin= binding_df[binding_df.iloc[:, -1].str.contains(lectin, na=False)]
lectin_binding_glycans= lectin.dropna(axis=1, how="all")



# setindex to protein names and reasssign it to a vsriable
lectin = lectin_binding_glycans.drop(columns= 'target')

#melt the dataframe
lectin_melt = lectin.melt(id_vars = 'protein', var_name = 'glycan', value_name = 'binding')

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


