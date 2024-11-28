import pickle
import pandas as pd
from glycowork.motif.graph import compare_glycans

# Load glycan flexibility data
flex_data_path = '/data/glycan_graphs.pkl'
with open(flex_data_path, 'rb') as file:
    flex_data = pickle.load(file)

# Load binding data into a DataFrame
binding_data_path = '/data/glycan_binding.csv'
binding_df = pd.read_csv(binding_data_path)

# Convert the flexibility data into a DataFrame
flex_df = pd.DataFrame.from_dict(flex_data, orient='index').reset_index().rename(columns={"index": "glycan"})

# Reset the index to ensure proteins are treated as a column
binding_df = binding_df.reset_index()

# Reshape binding_df to have one row per glycan-protein pair (long format)
binding_long = binding_df.melt(
    id_vars=["index"],  # Treat the index (protein names) as the identifier
    var_name="glycan",  # Glycans were columns, now treated as variables
    value_name="binding_score"  # Binding scores are the values
)

# Rename "index" column to "protein" for clarity
binding_long.rename(columns={"index": "protein"}, inplace=True)


# Initialize an empty list to store rows where compare_glycans returns True
filtered_rows = []

# Loop through the reshaped binding data and flexibility data
for _, binding_row in binding_long.iterrows():
    for _, flex_row in flex_df.iterrows():
        if compare_glycans(binding_row['glycan'], flex_row['glycan']):
            # If glycans match, store the glycan-protein pair with binding and flex data
            matched_row = {
                "glycan": binding_row['glycan'],
                "protein": binding_row['protein'],
                "binding_score": binding_row['binding_score'],
            }
            filtered_rows.append(matched_row)

# Combine filtered rows into a DataFrame
result_df = pd.DataFrame(filtered_rows)

# Display the resulting DataFrame
print(result_df)
