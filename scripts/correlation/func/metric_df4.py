import pickle
import pandas as pd
from glycowork.motif.graph import compare_glycans, subgraph_isomorphism
from glycowork.motif.processing import get_class



# File paths
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


def extract_terminal_and_flexible(motif):
    """
    Extract terminal and flexible monosaccharides from a single glycan motif.

    Args:
        motif (str): Glycan motif.

    Returns:
        tuple: Terminal and flexible monosaccharides.
    """
    # Split the motif into monosaccharides
    terminal = motif.split("(")[0]  # First monosaccharide (before parentheses)
    flexible = []

    # Identify potential flexible monosaccharides
    if ")" in motif:
        parts = motif.split(")")  # Split at closing parentheses
        if len(parts) > 1:
            flexible_part = parts[-1].strip()
            if flexible_part:
                flexible.append(flexible_part)

    return terminal, flexible


def process_glycan(matched_glycan, properties, flex_data):
    """
    Process a glycan string to find nodes matching binding motifs and calculate metrics.
    """
    binding_motif = properties["motif"]
    matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs = [], [], [], []

    for motif in binding_motif:
        try:
            # Extract terminal and flexible nodes for the motif
            #terminal, flexible = extract_terminal_and_flexible(motif)
            #termini_list = [terminal] + flexible
            #print(f"Processing motif: {motif}, Glycan: {matched_glycan}, Termini: {termini_list}")

            is_present, matched_nodes = subgraph_isomorphism(
                matched_glycan, motif,
                return_matches=True,
                termini_list=["terminal", "flexible", "flexible"]
            )

            if is_present:
                found_motifs.append(motif)

                # Flatten matched_nodes if needed
                matched_nodes = [node for sublist in matched_nodes for node in sublist] \
                    if isinstance(matched_nodes[0],list) else matched_nodes
                print(f"Matched nodes: {matched_nodes}")

                # Retrieve graph for the matched glycan
                graph = flex_data.get(matched_glycan)
                print(f"Graph for matched glycan {matched_glycan}: {graph.__dict__ if graph else 'Graph not found'}")

                if graph:
                    graph_nodes = list(graph.nodes)  # Retrieve graph nodes
                    print(f"Graph nodes: {graph_nodes}")

                    # Map matched nodes to graph nodes using index // 2 logic
                    aligned_nodes = [graph_nodes[index // 2] for index in matched_nodes if
                                     (index // 2) < len(graph_nodes)]
                    print(f"Aligned nodes: {aligned_nodes}")

                    # Process every second aligned node in reverse order
                    for graph_node in reversed(aligned_nodes[::2]):  # Reverse and take every second node
                        print(f"Processing graph node: {graph_node}")

                        if graph_node in graph.nodes:  # Ensure the node exists in the graph
                            attributes = graph.nodes[graph_node]  # Access node attributes
                            print(f"Attributes for graph node {graph_node}: {attributes}")
                            print(f" ")

                            matching_monosaccharides.append(attributes.get("Monosaccharide", 0))
                            flexibility_weighted.append(attributes.get("weighted_mean_flexibility", 0))
                            sasa_weighted.append(attributes.get("Weighted Score", 0))
                        else:
                            print(f"Graph node {graph_node} not found for glycan {matched_glycan}")

        except Exception as e:
            print(f"Error processing motif {motif}: {str(e)}")

    return matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs


def compute_sasa_metrics(nodes):
    """
    Compute sum, mean, and max SASA metrics for a list of numeric scores.
    """

    if not nodes:
        return {
            "SASA_weighted": None,
            "SASA_weighted_max": None,
            "SASA_weighted_sum": None,
        }

    SASA_weighted_sum = sum(nodes)
    SASA_weighted = SASA_weighted_sum / len(nodes)
    SASA_weighted_max = max(nodes)

    return {
        "SASA_weighted": SASA_weighted,
        "SASA_weighted_max": SASA_weighted_max,
        "SASA_weighted_sum": SASA_weighted_sum,
    }


def compute_overall_flexibility(flexibility_values):
    """Compute the overall flexibility for a glycan."""
    return (sum(flexibility_values) / len(flexibility_values) if flexibility_values else None)


def metric_df(lectin, properties):
    """
    Generate a metrics DataFrame for a given lectin and its binding properties.

    Args:
        lectin (str): The lectin name.
        properties (dict): The binding properties containing "motif", "terminal", and "flexible".

    Returns:
        pd.DataFrame: A DataFrame containing glycan metrics.
    """
    # Load data
    flex_data, binding_df = load_data()
    filtered_df = filter_binding_data(binding_df, lectin)

    if filtered_df.empty:
        print(f"No binding data found for {lectin}.")
        return pd.DataFrame()

    glycan_scores = get_glycan_scores(filtered_df)
    metric_data = []

    # Filter glycans without branches
    glycans_no_branches = [glycan for glycan in glycan_scores if "]" not in glycan]

    for glycan in glycans_no_branches:
        binding_score = glycan_scores[glycan]

        # Find matching glycan
        matched_glycan = find_matching_glycan(flex_data, glycan)
        if not matched_glycan:
            #print(f"No matching glycan found for {glycan}, skipping.")
            continue

        # Process glycan
        try:
            matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs = process_glycan(
                matched_glycan, properties, flex_data
            )

            # Compute SASA metrics
            if matching_monosaccharides:
                sasa_metrics = compute_sasa_metrics(sasa_weighted)
                overall_flexibility = compute_overall_flexibility(flexibility_weighted)
                glycan_class = get_class(matched_glycan)

                # Debugging logs
                print(f"Processing glycan: {glycan}")
                print(f"Matching monosaccharides: {matching_monosaccharides}")
                print(f"Found motifs: {found_motifs}")
                print(f"SASA weighted: {sasa_weighted}")
                print(f"Flexibility weighted: {flexibility_weighted}")

                # Append computed metrics
                metric_data.append({
                    "glycan": glycan,
                    "binding_score": binding_score,
                    "SASA_weighted": sasa_metrics["SASA_weighted"],
                    "SASA_weighted_max": sasa_metrics.get("SASA_weighted_max", None),
                    "SASA_weighted_sum": sasa_metrics.get("SASA_weighted_sum", None),
                    "weighted_mean_flexibility": overall_flexibility,
                    "class": glycan_class,  # Add class column
                })
            #else:
             #   print(f"No matching monosaccharides found for glycan {glycan}, skipping metrics computation.")
        except Exception as e:
            print(f"Error processing glycan {glycan}: {e}")

    # Convert to DataFrame
    metric_df = pd.DataFrame(metric_data)

    # Debugging log for final DataFrame
    print("Final metrics DataFrame:")
    print(metric_df)

    return metric_df #, glycan_scores, flex_data, binding_df, filtered_df, properties




"""# Example usage
lectin = "CMA"
binding_motif = {
    "CMA": {
        "motif": ["Fuc(a1-2)Gal", "GalNAc"],
        "terminal": ["Fuc", "GalNAc"],
        "flexible": ["Gal"]
    }
}

metric_df, glycan_scores, flex_data, binding_df, filtered_df, properties= metric_df(lectin,binding_motif[lectin])
"""