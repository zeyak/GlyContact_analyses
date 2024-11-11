import pandas as pd

# Load glycan binding data
glycan_binding_df = pd.read_csv('/data/glycan_binding.csv')

#remove the columns with only nan values
glycan_binding_df = glycan_binding_df.dropna(axis=1, how='all')




