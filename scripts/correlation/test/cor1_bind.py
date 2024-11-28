import pandas as pd

# Load the files to examine their contents
binding_data_path = '/data/glycan_binding.csv'

# Read files into dataframes
binding_df = pd.read_csv(binding_data_path).T


