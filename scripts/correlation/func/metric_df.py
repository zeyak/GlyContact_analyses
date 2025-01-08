import pickle
from typing import Dict
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


# File paths
flex_data_path = 'data/glycan_graphs.pkl'
binding_data_path = 'data/20241206_glycan_binding.csv'

def load_data_pdb():
    """Load glycan flexibility data and binding data."""
    with open(flex_data_path, 'rb') as file:
        flex_data = pickle.load(file)
    return flex_data

def load_data():
    """Load glycan flexibility data and binding data."""
    with open(flex_data_path, 'rb') as file:
        flex_data_pdb_g = pickle.load(file)

    flex_data = {}
    invalid_graphs = []  # List to store invalid graphs

    for glycan, graph in zip(flex_data_pdb_g.keys(), flex_data_pdb_g.values()):  # Ensure correct glycan-graph pairing
        try:
            # Skip invalid graphs before converting
            if not hasattr(graph, 'neighbors') or not graph:
                invalid_graphs.append(glycan)
                continue

            converted_graph = convert_pdb_graph_to_glycowork(graph)
            flex_data[glycan] = converted_graph  # Store the converted graph for the glycan

        except Exception as e:
            print(f"Error converting graph for glycan {glycan}: {e}")

    # Load binding data
    binding_df = pd.read_csv(binding_data_path)

    return flex_data, binding_df, invalid_graphs

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
    """Compute the overall flexibility for a glycan."""
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

            # Access the graph directly from flex_data
            pdb_graph = flex_data.get(matched_glycan)
            if not pdb_graph:
                print(f"No graph found for glycan: {matched_glycan}")
                continue

            selected_mono = [node for node in matched_nodes if node in pdb_graph.nodes]
            print(f"Selected monosaccharides: {selected_mono}")

            # Extract attributes from graph nodes
            if hasattr(pdb_graph, "nodes"):
                print(f"Graph nodes: {pdb_graph.nodes(data=True)}")
                for mono in selected_mono:
                    try:
                        attributes = pdb_graph.nodes[mono]
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

    for glycan in glycan_scores:
        binding_score = glycan_scores[glycan]
        matched_glycan = find_matching_glycan(flex_data, glycan)

        if not matched_glycan:
            missing_glycans.append(glycan)
            continue

        # Process the matched glycan
        matching_monosaccharides, sasa_weighted, flexibility_weighted, found_motifs = process_glycan_with_motifs(
            matched_glycan, properties, flex_data)

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
