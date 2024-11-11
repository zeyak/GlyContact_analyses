import pickle
import pandas as pd
import networkx as nx

# Load the glycan graph data
with open('/data/glycan_graphs.pkl', 'rb') as file:
    data = pickle.load(file)

# Load glycan binding data
glycan_binding = pd.read_csv('/data/glycan_binding.csv')

# Initialize dictionary to store flexibility summaries
glycan_summary = {}

# Iterate over each glycan and its corresponding graph
for glycan, graph in data.items():
    # Check if graph is a NetworkX Graph object
    if isinstance(graph, nx.Graph):
        scores = []
        # Iterate over nodes to collect flexibility scores
        for _, attributes in graph.nodes(data=True):
            # Ensure 'Mean Score' exists in the node attributes
            if 'Mean Score' in attributes:
                scores.append(attributes['Mean Score'])

        # Only proceed if scores were found for the current glycan
        if scores:
            glycan_summary[glycan] = {
                'Mean Flexibility Score': sum(scores) / len(scores),
                'Max Flexibility Score': max(scores),
                'Min Flexibility Score': min(scores)
            }
        else:
            print(f"No flexibility scores found for glycan: {glycan}")
    else:
        print(f"Data for glycan {glycan} is not a graph.")

# Convert the glycan summary dictionary to a DataFrame for easy analysis
glycan_summary_df = pd.DataFrame(glycan_summary).T  # Transpose to have glycans as rows
print(glycan_summary_df.head())

# Assuming `glycan_summary_df` is your flexibility DataFrame and `glycan_binding` is your binding DataFrame
# Merge the flexibility and binding data on glycan names
combined_df = glycan_summary_df.join(glycan_binding.set_index('Glycan'))

# Correlation between mean flexibility score and each lectin’s z-score
correlations = combined_df.corr()
print("Correlation Matrix:\n", correlations)
