import pickle
import pandas as pd
from glycowork.motif.graph import compare_glycans

# Load glycan flexibility data from the pickle file
flex_data_path = '/data/glycan_graphs.pkl'
with open(flex_data_path, 'rb') as file:
    flex_data = pickle.load(file)

# Load binding data into a DataFrame and transpose it
binding_data_path = '/data/glycan_binding.csv'
binding_df = pd.read_csv(binding_data_path)


# Ensure binding data has a comparable format (index as glycan names)
binding_df_T = binding_df.T.reset_index()
binding_df_T.rename(columns={"index": "glycan"}, inplace=True)

# Initialize an empty list to store rows where compare_glycans returns True
filtered_rows = []

# Loop through the data frames
for _, binding_row in binding_df_T.iterrows():
    #for _, flex_row in flex_df.iterrows():
    for glycan_sequence, graph_nodes in flex_data.items():
        if compare_glycans(binding_row['glycan'], glycan_sequence):
            # Merge the rows if compare_glycans returns True
            #merged_row = pd.concat([binding_row, flex_row], axis=0)
            #filtered_rows.append(binding_row)
            matched_row = {
                "glycan": binding_row['glycan'],
                "protein": binding_row['protein'],
                "binding_score": binding_row['binding_score'],

            }

# Combine all filtered rows into a single DataFrame
result_df = pd.DataFrame(filtered_rows)

# Display the resulting DataFrame
print(result_df)
