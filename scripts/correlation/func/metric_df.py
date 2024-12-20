import pickle
from csv import DictWriter
from typing import Dict
import numpy as np

import pandas as pd
from glycowork.motif.graph import compare_glycans, subgraph_isomorphism
from glycowork.motif.graph import  glycan_to_nxGraph
from glycowork.motif.processing import get_class
import sys

from networkx.classes import Graph
import networkx as nx

# File paths
flex_data_path = 'data/glycan_graphs.pkl'
binding_data_path = 'data/20241206_glycan_binding.csv'


def load_data():
    """Load glycan flexibility data and binding data."""
    with open(flex_data_path, 'rb') as file:
        flex_data = pickle.load(file)
    binding_df = pd.read_csv(binding_data_path)
    return flex_data, binding_df

def filter_binding_data(binding_df: pd.DataFrame, lectin: str) -> pd.DataFrame:
    """Filter the binding DataFrame for the given lectin."""
    filtered_df = binding_df[binding_df.iloc[:, -1].eq(lectin)]
    filtered_df = filtered_df.dropna(axis=1, how='all')  # Drop columns with all NaN values
    return filtered_df

def get_glycan_scores(filtered_df: dict[str, float]) -> Dict[str, float]:
    """Calculate mean binding scores for glycans."""
    lectin_df = filtered_df.iloc[:, :-2]  # Exclude "protein" and "target" columns
    glycan_scores = lectin_df.mean(axis=0).to_dict()
    return glycan_scores

def find_matching_glycan(flex_data , glycan):
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
        if len(daniel_graph.nodes) % 2 == 1:
            dan_max = len(daniel_graph.nodes)
            luc_max = len(luc_graph.nodes)
            map_dict = dict(zip(range(0, dan_max, 2), range(luc_max, 0, -1)))
        else:
            dan_max = len(daniel_graph.nodes) - 1
            luc_max = len(luc_graph.nodes)
            map_dict = dict(zip(range(0, dan_max, 2), range(luc_max, 0, -1)))
        return map_dict
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

def compute_overall_flexibility(flexibility_values ):
    """Compute the overall flexibility for a glycan."""
    return sum(flexibility_values) / len(flexibility_values) if flexibility_values else None

def process_glycan_with_motifs_(matched_glycan: str,
                               properties: str,
                               flex_data: dict[str, Graph]):

    """
    Process a glycan string to find nodes matching binding motifs and calculate metrics.
    Handles both single and multiple binding motifs.
    """
    matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs = [], [], [], []

    motifs = properties["motif"]
    termini_list = properties["termini_list"]

    for motif, termini in zip(motifs, termini_list):
        try:
            # Perform subgraph isomorphism
            try:
                is_present, matched_nodes = subgraph_isomorphism(
                    matched_glycan, motif,
                    return_matches=True,
                    termini_list=termini
                )
                if not is_present:
                    continue
            except Exception as e:
                print(f"Subgraph isomorphism error for glycan {matched_glycan} with motif {motif}: {e}")
                continue

            found_motifs.append(motif)
            print("")
            print(f"matched_glycan: {matched_glycan}")
            print(f"Processing motif: {motif}")

            matched_nodes = [node for sublist in matched_nodes for node in sublist] \
                if isinstance(matched_nodes[0], list) else matched_nodes

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
                print(f"graph.nodes: {graph.nodes}")
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

def convert_pdb_graph_to_glycowork(g):
    """
    Converts a PDB-format graph to glycowork format.
    Args:
        g (nx.Graph): Input graph in PDB format with Monosaccharide and score attributes
    Returns:
        nx.Graph: Converted graph in glycowork format with monosaccharide/linkage split into separate nodes
    """
    new_g = nx.Graph()

    # Find reducing end (node connected to -R node 1)
    reducing_end = None
    for n in g.neighbors(1):
        reducing_end = n
        break

    if reducing_end is None:
        raise ValueError("Could not find reducing end (node connected to -R)")

    # Get node list in DFS order starting from reducing end (excluding -R node)
    nodes = [n for n in nx.dfs_preorder_nodes(g, reducing_end) if n != 1]

    # Create mapping dict
    mapping = {}
    mapping[1] = []  # -R maps to nothing

    # Calculate highest node number
    max_num = 2 * (len(g.nodes) - 2)

    # Map reducing end to highest number
    mapping[reducing_end] = [max_num]

    # Map rest of nodes in reverse order
    curr_num = 0
    for node in reversed(nodes[1:]):  # Skip reducing end
        mapping[node] = [curr_num, curr_num + 1]
        curr_num += 2

    # Create nodes in new graph
    for old_node, new_nodes in mapping.items():
        if not new_nodes:  # Skip -R
            continue

        # Get numerical attributes
        attrs = {k: v for k, v in g.nodes[old_node].items()
                if k not in ['Monosaccharide', 'string_labels']}
        mono_link = g.nodes[old_node]['Monosaccharide']

        if len(new_nodes) == 1:  # Reducing end
            mono = mono_link.split('(')[0]
            new_g.add_node(new_nodes[0], string_labels=mono, **attrs)
        else:  # Normal nodes
            mono = mono_link.split('(')[0]
            link = mono_link.split('(')[1][:-1]

            new_g.add_node(new_nodes[0], string_labels=mono, **attrs)
            new_g.add_node(new_nodes[1], string_labels=link, **attrs)
            new_g.add_edge(new_nodes[0], new_nodes[1])

    # Add edges between components
    for old_u, old_v in g.edges():
        new_u = mapping.get(old_u, [])
        new_v = mapping.get(old_v, [])

        if not new_u or not new_v:  # Skip if either is -R
            continue

        if len(new_u) == 1:  # u is reducing end
            u_node = new_u[0]
        else:
            u_node = new_u[0]  # Use monosaccharide node

        if len(new_v) == 1:  # v is reducing end
            v_node = new_v[0]
        else:
            v_node = new_v[1]  # Use linkage node

        new_g.add_edge(u_node, v_node)

    return new_g

def process_glycan_with_motifs(matched_glycan: str,
                               properties: dict,
                               flex_data: dict[str, nx.Graph]):
    """
    Process a glycan string to find nodes matching binding motifs and calculate metrics.
    Handles both single and multiple binding motifs.
    Includes the integration of convert_pdb_graph_to_glycowork.

    Args:
        matched_glycan (str): Identifier of the glycan to process.
        properties (dict): Properties including motifs and termini lists.
        flex_data (dict): Dictionary mapping glycan identifiers to graphs.

    Returns:
        tuple: Matching monosaccharides, SASA-weighted scores, flexibility-weighted scores, and found motifs.
    """
    matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs = [], [], [], []

    motifs = properties.get("motif", [])
    termini_list = properties.get("termini_list", [])

    for motif, termini in zip(motifs, termini_list):
        try:
            # Perform subgraph isomorphism
            try:
                is_present, matched_nodes = subgraph_isomorphism(
                    matched_glycan, motif,
                    return_matches=True,
                    termini_list=termini
                )
                if not is_present:
                    continue
            except Exception as e:
                print(f"Subgraph isomorphism error for glycan {matched_glycan} with motif {motif}: {e}")
                continue

            found_motifs.append(motif)
            print(f"Processing motif: {motif} for glycan: {matched_glycan}")

            matched_nodes = [node for sublist in matched_nodes for node in sublist] \
                if isinstance(matched_nodes[0], list) else matched_nodes
            print(f"Matched nodes: {matched_nodes}")

            # Convert PDB graph to glycowork format and select nodes directly
            try:
                pdb_graph = flex_data.get(matched_glycan)
                if pdb_graph:
                    glycowork_graph = convert_pdb_graph_to_glycowork(pdb_graph)

                selected_mono = [node for node in matched_nodes if node in glycowork_graph.nodes]
                print(f"Selected monosaccharides: {selected_mono}")

            except Exception as e:
                print(f"Error processing pdb graph of {matched_glycan}: {e}")
                continue

            # Extract attributes from graph nodes
            if glycowork_graph and hasattr(glycowork_graph, "nodes"):
                print(f"Graph nodes: {glycowork_graph.nodes(data=True)}")
                for mono in selected_mono:
                    try:
                        attributes = glycowork_graph.nodes[mono]
                        matching_monosaccharides.append(attributes.get("string_labels", ""))
                        sasa_weighted.append(attributes.get("Weighted Score", 0))
                        flexibility_weighted.append(attributes.get("weighted_mean_flexibility", 0))

                        print(f"Matching monosaccharides: {matching_monosaccharides}")
                        print(f"SASA-weighted scores: {sasa_weighted}")
                        print(f"Flexibility-weighted scores: {flexibility_weighted}")
                        print("")
                    except Exception as e:
                        print(f"Error extracting attributes for node {mono} in glycan {matched_glycan}: {e}")
            else:
                print(f"Skipping invalid graph or graph with no nodes for glycan: {matched_glycan}")

        except Exception as e:
            print(f"General error processing glycan {matched_glycan} with motif {motif}: {e}")

    return matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs






def generate_metrics_for_glycan(properties:str,
                                glycan_scores: dict,
                                flex_data: dict[str, Graph]) -> list[dict]:
    """
    Generate metrics for each glycan by processing them one by one and creating metrics at the end.
    """
    metric_data = []
    missing_glycans = []

    # Filter glycans without branches
    #glycans_no_branches = [glycan for glycan in glycan_scores if "]" not in glycan]
    #print(f"Processing {len(glycans_no_branches)} glycans without branches...")

    #for glycan in glycans_no_branches:
    for glycan in glycan_scores:
        binding_score = glycan_scores[glycan]
        matched_glycan = find_matching_glycan(flex_data, glycan)

        if not matched_glycan:
            missing_glycans.append(glycan)
            continue

        # Process the matched glycan
        matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs = process_glycan_with_motifs(
            matched_glycan, properties, flex_data
        )

        if matching_monosaccharides:
            sasa_metrics = compute_sasa_metrics(sasa_weighted)
            overall_flexibility = compute_overall_flexibility(flexibility_weighted)
            # write nan if returns an empty string
            glycan_class = get_class(matched_glycan)
            if glycan_class == "":
                glycan_class = np.nan

            metric_data.append({
                "glycan": glycan,
                "binding_score": binding_score,
                "SASA_weighted": sasa_metrics["SASA_weighted"],
                "SASA_weighted_max": sasa_metrics.get("SASA_weighted_max"),
                "SASA_weighted_sum": sasa_metrics.get("SASA_weighted_sum"),
                "weighted_mean_flexibility": overall_flexibility,
                "class": glycan_class,
            })

    print(f"Processed {len(metric_data)} glycans without branches.")

    # Report missing glycans
    if missing_glycans:
        print(f"Missing glycans: {missing_glycans}")

    return metric_data

def metric_df(lectin, properties):
    """
    Generate a metrics DataFrame for a given lectin and its properties.
    """
    flex_data, binding_df = load_data()
    filtered_df = filter_binding_data(binding_df, lectin)
    if filtered_df.empty:
        print(f"No binding data found for {lectin}.")
        return pd.DataFrame()

    glycan_scores = get_glycan_scores(filtered_df)
    metric_data = generate_metrics_for_glycan(properties, glycan_scores, flex_data)

    metric_df = pd.DataFrame(metric_data)
    metric_df.set_index('glycan', inplace=True)
    metric_df.to_excel(f'scripts/correlation/metric_df/{lectin}_metrics.xlsx', index=True, header=True)
    return metric_df

# Example usage
sys.stdout = open('scripts/correlation/metric_df/metric_df.log', 'w')

lectin= 'AOL'

properties= {'motif': ['Fuc'], 'termini_list': [['t']]}

matched_glycan= 'Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)[Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-6)]Gal(b1-4)Glc'

df= metric_df(lectin, properties)
