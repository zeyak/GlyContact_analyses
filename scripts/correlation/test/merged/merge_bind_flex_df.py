import pickle
import pandas as pd
from glycowork.motif.graph import compare_glycans

from scripts.correlation.binding.binding_ConA import lectin

"""
Merged binding and flexibility data in a heatmap:

binding_df = 2745 glycans 1465 proteins
flex_df = 582 glycans

merged_df_glycans = 265 glycans 1465 proteins


This is not the best way to filter lectins,
because merged265 has the glycans as index and the lectins as columns.
therefore lectin names are gone on the columns.

"""


# Load glycan flexibility data from the pickle file
flex_data_path = '/data/glycan_graphs.pkl'
with open(flex_data_path, 'rb') as file:
    flex_data = pickle.load(file)

# Load binding data into a DataFrame and transpose it
binding_data_path = '/data/glycan_binding.csv'
binding_df = pd.read_csv(binding_data_path)

# Convert the flexibility data (a dictionary) into a DataFrame
flex_df = pd.DataFrame.from_dict(flex_data, orient='index')
flex_df = flex_df.reset_index().rename(columns={"index": "glycan"})

# Ensure binding data has a comparable format (index as glycan names)
binding_df_T = binding_df.T.reset_index()
binding_df_T.rename(columns={"index": "glycan"}, inplace=True)

# Initialize an empty list to store rows where compare_glycans returns True
filtered_rows = []

# Loop through the data frames
for _, binding_row in binding_df_T.iterrows():
    for _, flex_row in flex_df.iterrows():
        if compare_glycans(binding_row['glycan'], flex_row['glycan']):
            # Merge the rows if compare_glycans returns True
            #merged_row = pd.concat([binding_row, flex_row], axis=0)
            filtered_rows.append(binding_row)

# Combine all filtered rows into a single DataFrame
merged265_df = pd.DataFrame(filtered_rows)

# Display the resulting DataFrame
print(merged265_df)
