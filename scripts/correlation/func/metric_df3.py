import pickle
import pandas as pd
from glycowork.motif.graph import compare_glycans, subgraph_isomorphism

flex_data_path = 'data/glycan_graphs.pkl'
binding_data_path = 'data/glycan_binding.csv'


def load_data():
    """Load glycan flexibility data and binding data."""
    with open(flex_data_path, 'rb') as file:
        flex_data = pickle.load(file)
    binding_df = pd.read_csv(binding_data_path)
    return flex_data, binding_df


def filter_binding_data(binding_df, lectin):
    """Filter the binding DataFrame for the given lectin."""
    filtered_df = binding_df[binding_df.iloc[:, -1].eq(lectin)]
    filtered_df = filtered_df.dropna(axis=1, how='all')  # Drop columns with all NaN values
    return filtered_df


def get_glycan_scores(filtered_df):
    """Calculate mean binding scores for glycans."""
    lectin_df = filtered_df.iloc[:, :-2]  # Exclude "protein" and "target" columns
    glycan_scores = lectin_df.mean(axis=0).to_dict()
    return glycan_scores


def find_matching_glycan(flex_data, glycan):
    """Find the matching glycan in flex_data."""
    for flex_glycan in flex_data.keys():
        if compare_glycans(glycan, flex_glycan):
            return flex_glycan
    return None


def compute_sasa_metrics(nodes):
    """Compute average SASA metrics (mean, median, and weighted score)."""
    if not nodes:
        return {"SASA_mean": None, "SASA_median": None, "SASA_weighted": None}

    SASA_mean = sum(node.get("Mean Score", 0) for node in nodes) / len(nodes)
    SASA_median = sum(node.get("Median Score", 0) for node in nodes) / len(nodes)
    SASA_weighted = sum(node.get("Weighted Score", 0) for node in nodes) / len(nodes)

    return {"SASA_mean": SASA_mean, "SASA_median": SASA_median, "SASA_weighted": SASA_weighted}

def process_glycan(matched_flex_glycan, binding_motifs, graph):
    """
    Process a glycan graph to identify nodes matching multiple motifs and calculate metrics.

    Args:
        matched_flex_glycan: The glycan string being processed.
        binding_motifs: A list of motifs to search for.
        graph: The graph representation of the glycan.

    Returns:
        Tuple of node attributes (matching monosaccharides), flexibility values, and found motifs.
    """
    # Initialize lists for filtered metrics and found motifs
    matching_monosaccharides = []  # Nodes matching the binding motifs
    flexibility_values = []  # Weighted flexibility values for matching nodes
    found_motifs = []  # List of motifs found in the glycan

    # Iterate over each binding motif
    for binding_motif in binding_motifs:
        # Use subgraph isomorphism to identify matching nodes for the current motif
        is_present, matched_nodes = subgraph_isomorphism(matched_flex_glycan, binding_motif, return_matches=True)

        if not is_present:
            # If the current motif is not present, skip to the next
            continue

        # Add the current motif to the list of found motifs
        found_motifs.append(binding_motif)

        # Flatten matched nodes (in case of nested lists)
        if isinstance(matched_nodes[0], list):
            matched_nodes = [node for sublist in matched_nodes for node in sublist]

        # Process nodes in the graph to extract attributes for matched nodes
        for node in matched_nodes:
            if node in graph.nodes:
                attributes = graph.nodes[node]
                matching_monosaccharides.append(attributes)
                flexibility_values.append(attributes.get("weighted_mean_flexibility", 0))

    return matching_monosaccharides, flexibility_values, found_motifs


def metric_df(lectin, binding_motif):
    """Main function to compute the metric DataFrame."""
    flex_data, binding_df = load_data()
    filtered_df = filter_binding_data(binding_df, lectin)
    glycan_scores = get_glycan_scores(filtered_df)

    metric_data = []

    for glycan, binding_score in glycan_scores.items():
        matched_flex_glycan = find_matching_glycan(flex_data, glycan)

        if not matched_flex_glycan:
            continue

        graph = flex_data[matched_flex_glycan]
        if not hasattr(graph, "nodes"):  # Skip invalid graphs
            print(f"Skipping invalid graph for glycan: {matched_flex_glycan}")
            continue

        matching_monosaccharides, flexibility_values, found_motifs = process_glycan(
            matched_flex_glycan, binding_motif, graph)

        # Print the motifs that were found
        print(f"Motif(s) found in glycan {matched_flex_glycan}: {', '.join(found_motifs)}")

        # Compute SASA metrics and overall flexibility
        sasa_metrics = compute_sasa_metrics(matching_monosaccharides)
        overall_flexibility = (
            sum(flexibility_values) / len(flexibility_values) if flexibility_values else None
        )

        metric_data.append({
            "glycan": glycan,
            "binding_score": binding_score,
            "SASA_mean": sasa_metrics["SASA_mean"],
            "SASA_median": sasa_metrics["SASA_median"],
            "SASA_weighted": sasa_metrics["SASA_weighted"],
            "weighted_mean_flexibility": overall_flexibility,
        })

    metric_df = pd.DataFrame(metric_data)
    return metric_df #binding_df, filtered_df


# Example usage
#lectin = "PNA"
#binding_motif = "Gal(b1-3)GalNAc"


#lectin = "CMA"
#binding_motif = ["Fuc(a1-2)Gal", "GalNAc"]

#metric_df, binding_df, filtered_df = metric_df(lectin, binding_motif)
