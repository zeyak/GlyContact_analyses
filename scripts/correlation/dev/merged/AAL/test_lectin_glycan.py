import pandas as pd
from glycowork.motif.graph import compare_glycans
import pickle


lectin = "PNA"
binding_motif = ["Gal(b1-3)", "GalNAc"]


# Load glycan flexibility and SASA data from the pickle file
flex_data_path = 'data/glycan_graphs.pkl'
with open(flex_data_path, 'rb') as file:
    flex_data = pickle.load(file)

# Load binding data into a DataFrame and transpose it
binding_data_path = 'data/glycan_binding.csv'
binding_df = pd.read_csv(binding_data_path)


# Filter binding_df for rows where the last column exactly matches "AAL"
binding_df_filtered = binding_df[binding_df.iloc[:, -1].eq(lectin)]
binding_df_filtered = binding_df_filtered.dropna(axis=1, how='all')  # Drop columns with all NaN values

# Initialize a list to store glycan-specific rows
glycan_data = []


"""
there are 4 lectin called AAL in the binding_df
Should I take the average?
or take them all?
"""
# Iterate over glycans in binding_df_filtered
#AAL_row = binding_df_filtered.iloc[0, :-2]  # Exclude the "protein" and "target" column
#glycan_scores = AAL_row.to_dict()  # Convert to dictionary of glycan scores

AAL_df = binding_df_filtered.iloc[:4, :-2]  # First 4 rows, excluding the last 2 columns#
glycan_scores = AAL_df.mean(axis=0).to_dict()  # Calculate the mean score for each glycan
#glycan_scores = AAL_df.to_dict(orient="list")



for glycan, binding_score in glycan_scores.items():
    # Match the glycan using compare_glycan
    matched_flex_glycan = None
    for flex_glycan in flex_data.keys():
        if compare_glycans(glycan, flex_glycan):  # Use the compare_glycan function
            matched_flex_glycan = flex_glycan
            break

    if not matched_flex_glycan:  # Skip if no match found in flex_data
        continue


    graph = flex_data[matched_flex_glycan]
    if not hasattr(graph, "nodes"):  # Skip invalid graphs
        print(f"Skipping invalid graph for glycan: {matched_flex_glycan}")
        continue

    # Initialize lists for filtered metrics
    relevant_nodes = []
    weighted_flex_scores = []

    # Process nodes in the graph
    for node, attributes in graph.nodes(data=True):
        # Check if the node's Monosaccharide matches any binding motif
        if any(motif in attributes.get("Monosaccharide", "") for motif in binding_motif): #plot behavior changes if this is ignored
            relevant_nodes.append(attributes)
            weighted_flex_scores.append(attributes.get("weighted_mean_flexibility", 0))
            # Print relevant node information
            print(f"Node ID: {node}, Attributes: {attributes}")


    # Compute the average SASA metrics for filtered nodes
    def compute_average_metrics(nodes):
        if not nodes:
            return {"mean": None, "median": None, "weighted_score": None}
        total_mean = sum(node["Mean Score"] for node in nodes) / len(nodes)
        total_median = sum(node["Median Score"] for node in nodes) / len(nodes)
        total_weighted = sum(node["Weighted Score"] for node in nodes) / len(nodes)
        return {"mean": total_mean, "median": total_median, "weighted_score": total_weighted}

        """    def compute sum():
        return sum(weighted_flex_scores) / len(weighted_flex_scores) if weighted_flex_scores else None

    def compute max():
        return max(weighted_flex_scores) if weighted_flex_scores else None
        """


    average_metrics = compute_average_metrics(relevant_nodes)

    # Calculate overall weighted mean flexibility score for filtered nodes
    overall_flex_mean = (
        sum(weighted_flex_scores) / len(weighted_flex_scores)
        if weighted_flex_scores
        else None
    )

    # Create a row for the glycan
    glycan_row = {
        "glycan": glycan,
        "binding_score": binding_score,
        "SASA_mean": average_metrics["mean"],
        "SASA_median": average_metrics["median"],
        "SASA_weighted": average_metrics["weighted_score"],
        "weighted_mean_flexibility": overall_flex_mean,
    }
    glycan_data.append(glycan_row)

# Convert the glycan data into a DataFrame
final_df = pd.DataFrame(glycan_data)

# Display the resulting DataFrame
print(final_df)

