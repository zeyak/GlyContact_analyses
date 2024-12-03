import pickle
import pandas as pd
from glycowork.motif.graph import compare_glycans

flex_data_path = 'data/glycan_graphs.pkl'
binding_data_path = 'data/glycan_binding.csv'

"""
updated version of metric_df:

take the mean of the binding scores in case of multiple rows for the same lectin:
glycan_scores = lectin_df.mean(axis=0).to_dict()

glycan_data=> metric_data
glycan_row=> metric_row



"""



def metric_df(lectin, binding_motif):
    # Initialize a list to store glycan-specific rows in the metric_df
    metric_data = []

    # Load glycan flexibility and SASA data from the pickle file
    with open(flex_data_path, 'rb') as file:
        flex_data = pickle.load(file)

    # Load binding data into a DataFrame
    binding_df = pd.read_csv(binding_data_path)
    # Filter binding_df for rows where the last column exactly matches "lectin"
    binding_df_filtered = binding_df[binding_df.iloc[:, -1].eq(lectin)]
    binding_df_filtered = binding_df_filtered.dropna(axis=1, how='all')  # Drop columns with all NaN values
    # Iterate over glycans in binding_df_filtered
    lectin_df = binding_df_filtered.iloc[:, :-2]  # Exclude the "protein" and "target" column
    glycan_scores = lectin_df.mean(axis=0).to_dict()  # Calculate the mean score for each glycan in case of multiple rows

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
            print(graph)
            continue

        # Initialize lists for filtered metrics
        matching_monosaccharides = []  # Nodes matching the binding motif
        flexibility_values = []  # Weighted flexibility values for matching nodes

        # Process nodes in the graph
        for node, attributes in graph.nodes(data=True):
            # Extract the monosaccharide attribute
            monosaccharide = attributes.get("Monosaccharide", "")

            if lectin == "SNA":
                # SNA binds if any of the motifs are present (OR logic)
                match_found = any(motif in monosaccharide for motif in binding_motif)
            elif lectin == "PNA":
                # PNA binds if all motifs are present (AND logic)
                match_found = all(motif in monosaccharide for motif in binding_motif)
            else:
                # Generic case: binds if any motif matches
                match_found = any(motif in monosaccharide for motif in binding_motif)

            # If a match is found, collect the relevant data
            if match_found:
                matching_monosaccharides.append(attributes)
                flexibility_values.append(attributes.get("weighted_mean_flexibility", 0))

        def compute_sasa_metrics(nodes):
            """Compute average SASA metrics (mean, median, and weighted score)."""
            if not nodes:
                return {"SASA_mean": None, "SASA_median": None, "SASA_weighted": None}

            SASA_mean = sum(node["Mean Score"] for node in nodes) / len(nodes)
            SASA_median = sum(node["Median Score"] for node in nodes) / len(nodes)
            SASA_weighted = sum(node["Weighted Score"] for node in nodes) / len(nodes)

            return {"SASA_mean": SASA_mean, "SASA_median": SASA_median, "SASA_weighted": SASA_weighted}

        # Compute SASA metrics and overall flexibility
        sasa_metrics = compute_sasa_metrics(matching_monosaccharides)
        overall_flexibility = (
            sum(flexibility_values) / len(flexibility_values) if flexibility_values else None
        )

        # Create a row summarizing the data for the glycan
        metric_row = {
            "glycan": glycan,
            "binding_score": binding_score,
            "SASA_mean": sasa_metrics["SASA_mean"],
            "SASA_median": sasa_metrics["SASA_median"],
            "SASA_weighted": sasa_metrics["SASA_weighted"],
            "weighted_mean_flexibility": overall_flexibility,
        }

        # Append the metric row to the results
        metric_data.append(metric_row)

    # Convert the glycan data into a DataFrame
    metric_df = pd.DataFrame(metric_data)

    return metric_df
