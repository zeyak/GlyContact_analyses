import pickle
import pandas as pd
import matplotlib.pyplot as plt

# Load the glycan graph data
with open('/data/glycan_graphs.pkl', 'rb') as file:
    glycan_flexibility_data = pickle.load(file)

# Ensure consistency by creating a DataFrame from the flexibility data dictionary
# Assuming glycan_flexibility_data is a dictionary of dictionaries
glycan_flexibility_df = pd.DataFrame.from_dict(glycan_flexibility_data, orient='index')

# Load glycan binding data
glycan_binding_df = pd.read_csv('/data/glycan_binding.csv')

# Handle NaN values in binding data (choose one option)
glycan_binding_df = glycan_binding_df.fillna(0)  # or use dropna()

# Calculate the average binding z-score for each glycan across all lectins
glycan_binding_avg = glycan_binding_df.mean(axis=0)

# Ensure both DataFrames have glycans as indices for proper merging
glycan_flexibility_df.index.name = 'Glycan'
glycan_binding_avg.index.name = 'Glycan'

# Merge the binding data with the flexibility data on glycans
merged_df = glycan_flexibility_df.merge(glycan_binding_avg.rename('Average Binding Score'), left_index=True, right_index=True)

# Check correlation between mean flexibility and binding strength
correlation = merged_df['Mean Flexibility Score'].corr(merged_df['Average Binding Score'])
print("Correlation between Mean Flexibility Score and Average Binding Score:", correlation)

# Plotting to visualize the relationship
plt.scatter(merged_df['Mean Flexibility Score'], merged_df['Average Binding Score'])
plt.xlabel('Mean Flexibility Score')
plt.ylabel('Average Binding Score')
plt.title('Mean Flexibility vs Binding Strength')
plt.show()
