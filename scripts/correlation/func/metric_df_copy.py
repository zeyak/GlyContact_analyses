import pickle
from typing import Dict, Tuple, List
import numpy as np
import pandas as pd
import pingouin as pg
from glycowork.motif.graph import compare_glycans, subgraph_isomorphism
from glycowork.motif.graph import  glycan_to_nxGraph
from glycowork.motif.processing import get_class
import sys
from networkx.classes import Graph
import networkx as nx
from sklearn.preprocessing import LabelEncoder

FLEX_DATA_PATH = 'data/glycan_graphs.pkl'
BINDING_DATA_PATH = 'data/20241206_glycan_binding.csv'


def load_data_pdb():
    """Load glycan flexibility data from PDB source."""
    with open(FLEX_DATA_PATH, 'rb') as file:
        return pickle.load(file)


def remove_and_concatenate_labels(graph):
    """Modify the graph by removing and concatenating labels for odd-indexed nodes."""
    nodes_to_remove = []
    for node in sorted(graph.nodes):
        if node % 2 == 1:
            neighbors = list(graph.neighbors(node))
            if len(neighbors) > 1:
                for i in range(len(neighbors)):
                    for j in range(i + 1, len(neighbors)):
                        graph.add_edge(neighbors[i], neighbors[j])
            predecessor = node - 1
            if predecessor in graph.nodes:
                predecessor_label = graph.nodes[predecessor].get("string_labels", "")
                current_label = graph.nodes[node].get("string_labels", "")
                graph.nodes[predecessor]["string_labels"] = f"{predecessor_label}({current_label})"
            nodes_to_remove.append(node)
    graph.remove_nodes_from(nodes_to_remove)


def trim_gcontact(graph):
    """Trim the G_contact graph by removing node 1 and reconnecting its neighbors."""
    if not isinstance(graph, nx.Graph):
        raise TypeError(f"Expected G_contact to be a networkx.Graph, but got {type(graph).__name__}")
    if 1 in graph:
        neighbors = list(graph.neighbors(1))
        if len(neighbors) > 1:
            for i in range(len(neighbors)):
                for j in range(i + 1, len(neighbors)):
                    graph.add_edge(neighbors[i], neighbors[j])
        graph.remove_node(1)


def compare_graphs_with_attributes(G_contact, G_work):
    """Compare two graphs using node attributes and return a mapping dictionary."""

    def node_match(attrs1, attrs2):
        return 'string_labels' in attrs1 and 'Monosaccharide' in attrs2 and attrs1['string_labels'] in attrs2[
            'Monosaccharide']

    matcher = nx.isomorphism.GraphMatcher(G_work, G_contact, node_match=node_match)
    mapping_dict = {node_g2: node_g for node_g, node_g2 in matcher.mapping.items()} if matcher.is_isomorphic() else {}
    return mapping_dict


def create_glycontact_annotated_graph(glycan: str, mapping_dict, flex_data_pdb_g) -> nx.Graph:
    """Create an annotated glyco-contact graph with flexibility attributes."""
    if glycan not in flex_data_pdb_g or not isinstance(flex_data_pdb_g[glycan], nx.Graph):
        raise ValueError(f"Invalid glycan input or glycan not found in flex_data_pdb_g: {glycan}")

    glycowork_graph = glycan_to_nxGraph(glycan)
    try:
        node_attributes = {node: flex_data_pdb_g[glycan].nodes[node] for node in flex_data_pdb_g[glycan].nodes}
    except KeyError:
        raise KeyError(f'The glycan {glycan} is not present in the flex database')

    mapped_attributes = {mapping_dict[node]: attr for node, attr in node_attributes.items() if node in mapping_dict}
    nx.set_node_attributes(glycowork_graph, mapped_attributes)
    return glycowork_graph


def load_data():
    """Load glycan flexibility and binding data, process graphs, and return results."""
    flex_data = load_data_pdb()
    binding_df = pd.read_csv(BINDING_DATA_PATH)
    invalid_graphs = [glycan for glycan in flex_data if not isinstance(flex_data[glycan], nx.Graph)]
    G_mapped = {}

    for glycan, G_contact in flex_data.items():
        if not hasattr(G_contact, 'neighbors') or not G_contact:
            invalid_graphs.append(glycan)
            continue

        try:
            G_work = glycan_to_nxGraph(glycan)
        except Exception as e:
            print(f"Error converting glycan {glycan} to networkx graph: {e}")
            continue

        remove_and_concatenate_labels(G_work)
        #G_contact = G_work.copy()
        trim_gcontact(G_contact)
        m_dict = compare_graphs_with_attributes(G_contact, G_work)
        G_mapped[glycan] = create_glycontact_annotated_graph(glycan, m_dict, flex_data)

    return G_mapped, binding_df, invalid_graphs


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
            return glycan
    return None

def compute_sasa_metrics(nodes):
    """Compute sum, mean, and max SASA metrics for a list of numeric scores."""
    if not nodes:
        return {"SASA_weighted": None, "SASA_weighted_max": None, "SASA_weighted_sum": None}
    SASA_weighted_sum = sum(nodes)
    SASA_weighted = SASA_weighted_sum / len(nodes)
    SASA_weighted_max = max(nodes)
    return {"SASA_weighted": SASA_weighted, "SASA_weighted_max": SASA_weighted_max, "SASA_weighted_sum": SASA_weighted_sum}

def compute_overall_flexibility(flexibility_values ):

    return sum(flexibility_values) / len(flexibility_values) if flexibility_values else None

def process_glycan_with_motifs(matched_glycan: str,
                               properties: dict,
                               flex_data: dict[str, nx.Graph]):
    """
    Process a glycan string to find nodes matching binding motifs and calculate metrics.
    Handles both single and multiple binding motifs.

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

        matched_nodes = [node for sublist in matched_nodes for node in sublist] if matched_nodes and isinstance(matched_nodes[0], list) else matched_nodes
        print(f"Matched nodes: {matched_nodes}")

        pdb_graph = flex_data.get(matched_glycan)
        if not isinstance(pdb_graph, nx.Graph):
            print(f"No valid graph found for glycan: {matched_glycan}")
            continue

        # Select only monosaccharides (even-indexed nodes)
        selected_mono = [node for node in matched_nodes if node in pdb_graph.nodes and node % 2 == 0]
        print(f"Selected monosaccharides: {selected_mono}")

        if hasattr(pdb_graph, "nodes"):
            print(f"Graph nodes: {pdb_graph.nodes(data=True)}")
            for mono in selected_mono:
                try:
                    attributes = pdb_graph.nodes[mono]
                    monosaccharide = attributes.get('Monosaccharide', "")
                    if monosaccharide:  # Ensure non-empty monosaccharides
                        matching_monosaccharides.append(monosaccharide)
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

    return matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs



def generate_metrics_for_glycan(properties: str,
                                glycan_scores: dict,
                                flex_data: dict[str, Graph]) -> list[dict]:
    """
    Generate metrics for each glycan by processing them one by one and creating metrics at the end.
    """
    metric_data = []
    missing_glycans = []

    for glycan in glycan_scores:
        binding_score = glycan_scores[glycan]
        matched_glycan = find_matching_glycan(flex_data, glycan)

        if not matched_glycan:
            missing_glycans.append(glycan)
            continue

        # Process the matched glycan
        matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs = process_glycan_with_motifs(
            matched_glycan, properties, flex_data)

        # Skip empty monosaccharides
        matching_monosaccharides = [m for m in matching_monosaccharides if m.strip()]

        if matching_monosaccharides:
            sasa_metrics = compute_sasa_metrics(sasa_weighted)
            overall_flexibility = compute_overall_flexibility(flexibility_weighted)
            glycan_class = get_class(matched_glycan) or np.nan

            metric_data.append({
                "glycan": glycan,
                "binding_score": binding_score,
                "SASA_weighted": sasa_metrics["SASA_weighted"],
                "SASA_weighted_max": sasa_metrics.get("SASA_weighted_max"),
                "SASA_weighted_sum": sasa_metrics.get("SASA_weighted_sum"),
                "weighted_mean_flexibility": overall_flexibility,
                "class": glycan_class,
            })

    print(f"Processed {len(metric_data)} glycans with metrics.")

    # Report missing glycans
    if missing_glycans:
        print(f"Not-matched glycan in flex data: {missing_glycans}")

    return metric_data



def metric_df(lectin, properties):
    """
    Generate a metrics DataFrame for a given lectin and its properties.
    """
    flex_data, binding_df, invalid_graphs = load_data()
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

sys.stdout = open('scripts/correlation/metric_df/metric_df.log', 'w')

def perform_mediation_analysis_with_class(metric_df, independent_var, class_var, dependent_var):
    """
    Perform mediation analysis to test if glycan class mediates the relationship between
    the independent variable and the dependent variable. Handles NaN values and single-class cases.

    Parameters:
        metric_df (pd.DataFrame): The DataFrame containing metrics for mediation analysis.
        independent_var (str): Name of the independent variable (e.g., 'weighted_mean_flexibility').
        class_var (str): Name of the mediator variable (categorical, e.g., 'class').
        dependent_var (str): Name of the dependent variable (e.g., 'binding_score').

    Returns:
        dict: A dictionary containing mediation results and total, direct, and indirect effects,
              or a message indicating mediation cannot be performed.
    """

    # Step 1: Drop rows with NaN in the class column
    metric_df = metric_df.dropna(subset=[class_var]).copy()

    # Step 2: Check the number of unique classes in the mediator column
    unique_classes = metric_df[class_var].unique()
    if len(unique_classes) < 2:
        return {
            'message': f"Mediation analysis cannot be performed. Found only one unique class: {unique_classes}."
        }

    # Step 3: Encode class as binary (0 and 1) using LabelEncoder
    label_encoder = LabelEncoder()
    metric_df.loc[:, 'encoded_class'] = label_encoder.fit_transform(metric_df[class_var])

    # Step 4: Perform mediation analysis
    mediation_results = pg.mediation_analysis(
        data=metric_df,
        x=independent_var,  # Independent variable
        m='encoded_class',  # Encoded mediator (binary glycan class)
        y=dependent_var,  # Dependent variable
        alpha=0.05
    )

    # Step 5: Extract mediation effects
    total_effect = mediation_results.loc[mediation_results['path'] == 'Total', 'coef'].values[0]
    direct_effect = mediation_results.loc[mediation_results['path'] == 'Direct', 'coef'].values[0]
    indirect_effect = mediation_results.loc[mediation_results['path'] == 'Indirect', 'coef'].values[0]

    return {
        'results': mediation_results,
        'total_effect': total_effect,
        'direct_effect': direct_effect,
        'indirect_effect': indirect_effect
    }
