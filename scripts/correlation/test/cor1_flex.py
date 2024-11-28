import pickle
import pandas as pd

# Define the path to the flexibility data file
data_path = '/data/glycan_graphs.pkl'

# Load glycan flexibility data from the pickle file
with open(data_path, 'rb') as file:
    glycan_data = pickle.load(file)

# Convert the loaded data into a DataFrame
flex_glycan_df = pd.DataFrame.from_dict(glycan_data, orient='index')




