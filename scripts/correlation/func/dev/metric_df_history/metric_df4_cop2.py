import pickle
import pandas as pd
from glycowork.motif.graph import compare_glycans, subgraph_isomorphism
from glycowork.motif.graph import  glycan_to_nxGraph
from glycowork.motif.processing import get_class
import sys


# File paths
flex_data_path = 'data/glycan_graphs.pkl'
binding_data_path = 'data/20241206_glycan_binding.csv'


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
            #return flex_glycan.keys(glycan)
            return glycan
    return None

def map_daniel_to_luc_graph(matched_glycan, flex_data, daniel_selected_nodes):
    """
    Map selected nodes from a Daniel graph to a corresponding Luc graph based on a predefined mapping.
    """
    daniel_graph = glycan_to_nxGraph(matched_glycan)
    luc_graph = flex_data.get(matched_glycan)

    # Define the mapping logic between Daniel and Luc graphs
    def create_mapping():
        node_count_daniel = len(daniel_graph.nodes)
        node_count_luc = len(luc_graph.nodes)
        if node_count_daniel % 2 == 1:
            max_index = node_count_daniel
        else:
            max_index = node_count_daniel - 1
        return dict(zip(range(0, max_index, 2), range(node_count_luc, 0, -1)))

    # Generate the mapping dictionary
    map_dict = create_mapping()

    # Map the selected nodes from the Daniel graph to the Luc graph
    luc_selected_nodes = []
    for node in daniel_selected_nodes:
        if node in map_dict:
            luc_selected_nodes.append(map_dict[node])

    return luc_selected_nodes

def compute_sasa_metrics(nodes):
    """Compute sum, mean, and max SASA metrics for a list of numeric scores."""
    if not nodes:
        return {"SASA_weighted": None, "SASA_weighted_max": None, "SASA_weighted_sum": None}
    SASA_weighted_sum = sum(nodes)
    SASA_weighted = SASA_weighted_sum / len(nodes)
    SASA_weighted_max = max(nodes)
    return {"SASA_weighted": SASA_weighted, "SASA_weighted_max": SASA_weighted_max, "SASA_weighted_sum": SASA_weighted_sum}

def compute_overall_flexibility(flexibility_values):
    """Compute the overall flexibility for a glycan."""
    return sum(flexibility_values) / len(flexibility_values) if flexibility_values else None

def process_glycan_with_motifs(matched_glycan, binding_motif, flex_data):
    """
    Process a glycan string to find nodes matching binding motifs and calculate metrics.
    Handles both single and multiple binding motifs.
    """
    matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs = [], [], [], []

    # Ensure binding_motif is iterable
    motifs = binding_motif if isinstance(binding_motif, (list, tuple)) else [binding_motif]

    for motif in motifs:
        try:
            # Perform subgraph isomorphism
            try:
                is_present, matched_nodes = subgraph_isomorphism(
                    matched_glycan, motif,
                    return_matches=True,
                    termini_list=["terminal", "flexible", "flexible"]
                )
                if not is_present:
                    #print(f"No isomorphism found for glycan {matched_glycan} with motif {motif}")
                    continue
            except Exception as e:
                print(f"Subgraph isomorphism error for glycan {matched_glycan} with motif {motif}: {e}")
                continue

            found_motifs.append(motif)
            print("")
            print(f"matched_glycan: {matched_glycan}")
            print(f"Processing motif: {motif}")

            # Flatten matched nodes
            try:
                matched_nodes = [node for sublist in matched_nodes for node in sublist] \
                    if isinstance(matched_nodes[0], list) else matched_nodes
                print(f"matched_nodes: {matched_nodes}")
            except Exception as e:
                print(f"Error flattening matched nodes for glycan {matched_glycan}: {e}")
                continue

            # Select monosaccharides
            try:
                matched_mono = [matched_nodes[index] for index in range(0, len(matched_nodes), 2)]
                print(f"matched_mono: {matched_mono}")
                selected_mono = map_daniel_to_luc_graph(matched_glycan, flex_data, matched_mono)
                print(f"selected_mono: {selected_mono}")
            except Exception as e:
                print(f"Mapping error in map_daniel_to_luc_graph for glycan {matched_glycan}: {e}")
                continue

            # Extract attributes from graph nodes
            graph = flex_data.get(matched_glycan)
            if graph and hasattr(graph, "nodes"):
                print (f"graph.nodes: {graph.nodes}")
                print(f"graph.attributes: {graph.nodes(data=True)}")
                for mono in selected_mono:
                    try:
                        attributes = graph.nodes[mono]
                        matching_monosaccharides.append(attributes.get("Monosaccharide", 0))
                        sasa_weighted.append(attributes.get("Weighted Score", 0))
                        flexibility_weighted.append(attributes.get("weighted_mean_flexibility", 0))
                        print(f"matching_monosaccharides: {matching_monosaccharides}")
                        print(f"mono attributes: {attributes}")
                    except Exception as e:
                        print(f"Error extracting attributes for node {mono} in glycan {matched_glycan}: {e}")
            else:
                print(f"Skipping invalid graph or graph with no nodes for glycan: {matched_glycan}")

        except Exception as e:
            print(f"General error processing glycan {matched_glycan} with motif {motif}: {e}")

    return matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs

def generate_metrics_for_glycan(binding_motif, glycan_scores, flex_data):
    """
    Generate metrics for each glycan by processing them one by one and creating metrics at the end.
    """
    metric_data = []
    missing_glycans = []

    # Filter glycans without branches
    glycans_no_branches = [glycan for glycan in glycan_scores if "]" not in glycan]
    print(f"Processing {len(glycans_no_branches)} glycans without branches...")

    for glycan in glycans_no_branches:
        binding_score = glycan_scores[glycan]
        matched_glycan = find_matching_glycan(flex_data, glycan)


        if not matched_glycan:
            missing_glycans.append(glycan)
            continue

        # Process the matched glycan
        matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs = process_glycan_with_motifs(
            matched_glycan, binding_motif, flex_data
        )


        if matching_monosaccharides:
            sasa_metrics = compute_sasa_metrics(sasa_weighted)
            overall_flexibility = compute_overall_flexibility(flexibility_weighted)
            glycan_class = get_class(matched_glycan)

            metric_data.append({
                "glycan": glycan,
                "binding_score": binding_score,
                "SASA_weighted": sasa_metrics["SASA_weighted"],
                "SASA_weighted_max": sasa_metrics.get("SASA_weighted_max"),
                "SASA_weighted_sum": sasa_metrics.get("SASA_weighted_sum"),
                "weighted_mean_flexibility": overall_flexibility,
                "class": glycan_class,
            })
    # sum glycans with no brunches
    print(f"Processed {len(metric_data)} glycans without branches.")

    # Report missing glycans
    if missing_glycans:
     print("")
     print(f"Missing glycans: {missing_glycans}")

    return metric_data

def metric_df(lectin, binding_motif):
    """Generate a metrics DataFrame for a given lectin and its binding binding_motif."""
    flex_data, binding_df = load_data()
    filtered_df = filter_binding_data(binding_df, lectin)
    if filtered_df.empty:
        print(f"No binding data found for {lectin}.")
        return pd.DataFrame()

    glycan_scores = get_glycan_scores(filtered_df)
    metric_data = generate_metrics_for_glycan( binding_motif, glycan_scores, flex_data)

    metric_df = pd.DataFrame(metric_data)
    metric_df.set_index('glycan', inplace=True)
    metric_df.to_excel(f'scripts/correlation/metric_df/{lectin}_metrics.xlsx', index=False, header=True)
    return metric_df


# Example usage
sys.stdout = open('scripts/correlation/metric_df/metric_df.log', 'w')








#deprecated
def map_matched_nodes_to_graph_nodes(matched_nodes, graph_nodes):
    """Align matched nodes with graph nodes."""
    if len(matched_nodes) % 2 == 1:
        new_max = len(matched_nodes)
        old_max = len(graph_nodes)
        map_dict = dict(zip(range(0, new_max, 2), range(old_max, 0, -1)))
    else:
        new_max = len(matched_nodes) - 1
        old_max = len(graph_nodes)
        map_dict = dict(zip(range(0, new_max, 2), range(old_max, 0, -1)))
    return [graph_nodes[index] for index in map_dict.keys() if index < len(graph_nodes)]
def map_daniel_to_luc_graph_(mathced_glycan, flex_data):
    daniel_graph=glycan_to_nxGraph(matched_glycan)
    luc_graph = flex_data.get(matched_glycan)
    if len(daniel_graph.nodes)%2 ==1:
        new_max = len(daniel_graph.nodes)
        old_max = len(luc_graph.nodes)
        map_dict = dict(zip(range(0, new_max, 2),range(old_max,0,-1)))
    else:
        new_max = len(daniel_graph.nodes)-1
        old_max = len(luc_graph.nodes)
        map_dict = dict(zip(range(0, new_max, 2),range(old_max,0,-1)))
    return map_dict
def map_daniel_to_luc_graph__(mathced_glycan, flex_data, daniel_selected_nodes):
    daniel_graph=glycan_to_nxGraph(matched_glycan)
    luc_graph = flex_data.get(matched_glycan)

    def map():
        if len(daniel_graph.nodes)%2 ==1:
            new_max = len(daniel_graph.nodes)
            old_max = len(luc_graph.nodes)
            map_dict = dict(zip(range(0, new_max, 2),range(old_max,0,-1)))
        else:
            new_max = len(daniel_graph.nodes)-1
            old_max = len(luc_graph.nodes)
            map_dict = dict(zip(range(0, new_max, 2),range(old_max,0,-1)))

    luc_selected_nodes = []
    for node in daniel_selected_nodes:
        if node in map_dict:
            luc_selected_nodes.remove(node)
            luc_selected_nodes.append(map_dict[node])

    return luc_selected_nodes
def process_binding_motif(matched_glycan, motif, flex_data):
    """Process a single binding motif and extract metrics."""
    matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs = [], [], [], []
    try:
        is_present, matched_nodes = subgraph_isomorphism(
            matched_glycan, motif,
            return_matches=True,
            termini_list=["terminal", "flexible", "flexible"]
        )
        if is_present:
            found_motifs.append(motif)
            matched_nodes = [node for sublist in matched_nodes for node in sublist] \
                if isinstance(matched_nodes[0], list) else matched_nodes
            graph = flex_data.get(matched_glycan)
            if graph:
                graph_nodes = list(graph.nodes)
                graph_attributes = list(graph.nodes(data=True))
                #mapped_nodes = map_matched_nodes_to_graph_nodes(matched_nodes, graph_nodes)

                mapped_nodes = map_daniel_to_luc_graph(matched_nodes, flex_data)
                # convert the values in aligned_nodes to a list
                aligned_nodes = list(mapped_nodes.values())
                for graph_node in aligned_nodes:
                    attributes = graph.nodes[graph_node]
                    matching_monosaccharides.append(attributes.get("Monosaccharide", 0))
                    sasa_weighted.append(attributes.get("Weighted Score", 0))
                    flexibility_weighted.append(attributes.get("weighted_mean_flexibility", 0))

    except Exception as e:
        print(f"Error processing motif {motif} for glycan {matched_glycan}: {e}")
    return matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs
def process_binding_motif_(matched_glycan, motif, flex_data):
    """Process a single binding motif and extract metrics."""
    matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs = [], [], [], []
    try:
        is_present, matched_nodes = subgraph_isomorphism(
            matched_glycan, motif,
            return_matches=True,
            termini_list=["terminal", "flexible", "flexible"]
        )
        if is_present:
            #flatten the list of matched nodes
            matched_nodes = [node for sublist in matched_nodes for node in sublist] \
                if isinstance(matched_nodes[0], list) else matched_nodes
            # take the monosaccharides
            matched_mono = [matched_nodes[index] for index in range(0, len(matched_nodes), 2)]
            #select teh monosaccharides
            selected_mono = map_daniel_to_luc_graph__(matched_glycan, flex_data, matched_mono)

            for mono in selected_mono:
                graph = flex_data.get(matched_glycan)
                if not hasattr(graph, "nodes"):  # Skip invalid graphs
                    print(f"Skipping invalid graph for glycan: {matched_glycan}")
                    continue

                attributes = graph.nodes[mono]
                matching_monosaccharides.append(attributes.get("Monosaccharide", 0))
                sasa_weighted.append(attributes.get("Weighted Score", 0))
                flexibility_weighted.append(attributes.get("weighted_mean_flexibility", 0))

    except Exception as e:
        print(f"Error processing motif {motif} for glycan {matched_glycan}: {e}")
    return matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs
def process_glycan( binding_motif, flex_data, matched_glycan):
    """Process a glycan string to find nodes matching binding motifs and calculate metrics."""
    matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs = [], [], [], []
    if binding_motif >1:
        for motif in binding_motif:
            mono, sasa, flex, motifs = process_binding_motif_(matched_glycan, motif, flex_data)
            matching_monosaccharides.extend(mono)
            sasa_weighted.extend(sasa)
            flexibility_weighted.extend(flex)
            found_motifs.extend(motifs)
    else:
        mono, sasa, flex, motifs = process_binding_motif_(matched_glycan, binding_motif, flex_data)
    return matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs
def process_glycan_with_motifs_(matched_glycan, binding_motif, flex_data):
    """
    Process a glycan string to find nodes matching binding motifs and calculate metrics.
    Handles both single and multiple binding motifs.
    """
    matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs = [], [], [], []

    # Ensure binding_motif is iterable
    motifs = binding_motif if isinstance(binding_motif, (list, tuple)) else [binding_motif]

    try:
        for motif in motifs:

            # Perform subgraph isomorphism
            is_present, matched_nodes = subgraph_isomorphism(
                matched_glycan, motif,
                return_matches=True,
                termini_list=["terminal", "flexible", "flexible"]
            )
            if is_present:
                found_motifs.append(motif)
                print(f"matched_glycan: {matched_glycan}")
                print(f"Processing motif: {motif}")

                # Flatten matched nodes
                matched_nodes = [node for sublist in matched_nodes for node in sublist] \
                    if isinstance(matched_nodes[0], list) else matched_nodes
                print(f"matched_nodes: {matched_nodes}")

                # Select monosaccharides
                matched_mono = [matched_nodes[index] for index in range(0, len(matched_nodes), 2)]
                print(f"matched_mono: {matched_mono}")
                selected_mono = map_daniel_to_luc_graph__(matched_glycan, flex_data, matched_mono)
                print(f"selected_mono: {selected_mono}")

                # Extract attributes from graph nodes
                graph = flex_data.get(matched_glycan)
                if graph and hasattr(graph, "nodes"):
                    for mono in selected_mono:
                        attributes = graph.nodes[mono]
                        matching_monosaccharides.append(attributes.get("Monosaccharide", 0))
                        sasa_weighted.append(attributes.get("Weighted Score", 0))
                        flexibility_weighted.append(attributes.get("weighted_mean_flexibility", 0))
                        print (f"attributes: {attributes}")
                        print (f"matching_monosaccharides: {matching_monosaccharides}")

                else:
                    print(f"Skipping invalid graph for glycan: {matched_glycan}")

    except Exception as e:
        print(f"Error processing glycan {matched_glycan} with motif {motif}: {e}")

    return matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs
def generate_metrics_for_glycan_( binding_motif, glycan_scores, flex_data):
    """Generate metrics for each glycan."""
    metric_data = []
    missing_glycans = []
    glycans_no_branches = [glycan for glycan in glycan_scores if "]" not in glycan]
    for glycan in glycans_no_branches:
        binding_score = glycan_scores[glycan]
        matched_glycan = find_matching_glycan(flex_data, glycan)
        if not matched_glycan:
            missing_glycans.append(glycan)
            continue

        (matching_monosaccharides,
         sasa_weighted,
         flexibility_weighted,
         found_motifs) = process_glycan(matched_glycan,
             binding_motif, flex_data)

        if matching_monosaccharides:
            sasa_metrics = compute_sasa_metrics(sasa_weighted)
            overall_flexibility = compute_overall_flexibility(flexibility_weighted)
            glycan_class = get_class(matched_glycan)
            metric_data.append({
                "glycan": glycan,
                "binding_score": binding_score,
                "SASA_weighted": sasa_metrics["SASA_weighted"],
                "SASA_weighted_max": sasa_metrics.get("SASA_weighted_max", None),
                "SASA_weighted_sum": sasa_metrics.get("SASA_weighted_sum", None),
                "weighted_mean_flexibility": overall_flexibility,
                "class": glycan_class,
            })
    return metric_data
def generate_metrics_for_glycan__(binding_motif, glycan_scores, flex_data):
    """Generate metrics for each glycan."""
    metric_data = []
    missing_glycans = []

    # Filter out glycans with branches
    glycans_no_branches = [glycan for glycan in glycan_scores if "]" not in glycan]

    for glycan in glycans_no_branches:
        binding_score = glycan_scores[glycan]
        matched_glycan = find_matching_glycan(flex_data, glycan)

        if not matched_glycan:
            missing_glycans.append(glycan)
            continue

        try:
            (matching_monosaccharides,
             sasa_weighted,
             flexibility_weighted,
             found_motifs) = process_glycan(matched_glycan, binding_motif, flex_data)
        except Exception as e:
            print(f"Error processing glycan {glycan}: {e}")
            continue

        if matching_monosaccharides:
            sasa_metrics = compute_sasa_metrics(sasa_weighted) if sasa_weighted else {
                "SASA_weighted": None, "SASA_weighted_max": None, "SASA_weighted_sum": None
            }
            overall_flexibility = compute_overall_flexibility(flexibility_weighted) or None
            glycan_class = get_class(matched_glycan)

            metric_data.append({
                "glycan": glycan,
                "binding_score": binding_score,
                "SASA_weighted": sasa_metrics["SASA_weighted"],
                "SASA_weighted_max": sasa_metrics.get("SASA_weighted_max"),
                "SASA_weighted_sum": sasa_metrics.get("SASA_weighted_sum"),
                "weighted_mean_flexibility": overall_flexibility,
                "class": glycan_class,
            })

    # Report missing glycans if any
    if missing_glycans:
        print(f"Missing glycans: {missing_glycans}")

    return metric_data
def generate_metrics_for_glycan___(binding_motif, glycan_scores, flex_data):
    """
    Generate metrics for each glycan by processing them one by one and creating metrics at the end.
    """
    metric_data = []
    missing_glycans = []
    glycans_no_branches = [glycan for glycan in glycan_scores if "]" not in glycan]

    for glycan in glycans_no_branches:
        binding_score = glycan_scores[glycan]
        matched_glycan_dict = find_matching_glycan(flex_data, glycan)

        if not matched_glycan_dict:
            missing_glycans.append(glycan)
            continue

        # Process each matched glycan in the dictionary
        for matched_glycan in matched_glycan_dict.keys():
            metrics = process_glycan(
                matched_glycan, binding_motif, flex_data
            )
            if metrics:
                metric_data.append(metrics)

    # Report missing glycans if any
    if missing_glycans:
        print(f"Missing glycans: {missing_glycans}")

    return metric_data

"""
matched_glycan = 'Fuc(a1-2)Gal(b1-3)GlcNAc'
lectin = "CMA"
binding_motif = {
    "CMA": {
        "motif": ["Fuc(a1-2)Gal", "GalNAc"],
        "terminal": ["Fuc", "GalNAc"],
        "flexible": ["Gal"]
    }
}

result_df = metric_df(lectin, binding_motif[lectin])
print(result_df)
"""

