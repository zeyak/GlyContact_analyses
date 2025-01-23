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


# File paths
FLEX_DATA_PATH = 'data/glycan_graphs.pkl'
BINDING_DATA_PATH = 'data/20241206_glycan_binding.csv'


def load_data() -> Tuple[Dict[str, nx.Graph], pd.DataFrame, List[str]]:
    """Load glycan flexibility data, convert graphs, and load binding data."""
    flex_data, invalid_graphs = load_and_convert_flex_data(FLEX_DATA_PATH)
    binding_df = pd.read_csv(BINDING_DATA_PATH)
    return flex_data, binding_df, invalid_graphs


def load_and_convert_flex_data(file_path: str) -> Tuple[Dict[str, nx.Graph], List[str]]:
    """Load glycan flexibility data and convert valid graphs."""
    with open(file_path, 'rb') as file:
        flex_data_pdb_g = pickle.load(file)

    flex_data = {}
    invalid_graphs = []

    for glycan, graph in flex_data_pdb_g.items():
        try:
            # Skip invalid graphs
            if not hasattr(graph, 'neighbors') or not graph:
                invalid_graphs.append(glycan)
                continue
            converted_graph = create_glycontact_annotated_graph(glycan, flex_data_pdb_g)
            flex_data[glycan] = converted_graph
        except Exception as e:
            print(f"Error converting graph for glycan {glycan}: {e}")
            invalid_graphs.append(glycan)

    return flex_data, invalid_graphs


def create_glycontact_annotated_graph(glycan: str, flex_data_pdb_g: Dict[str, nx.Graph]) -> nx.Graph:
    """Create a glyco-contact annotated graph with flexibility attributes."""
    glycowork_graph = glycan_to_nxGraph(glycan)
    num_nodes = len(glycowork_graph) - 1
    mapper = {gcontact: gwork for gwork, gcontact in zip(range(num_nodes, -1, -2), range(2, num_nodes))}

    try:
        node_attributes = {node: flex_data_pdb_g[glycan].nodes[node]
                           for node in flex_data_pdb_g[glycan].nodes}
    except KeyError:
        raise KeyError(f'The glycan {glycan} is not present in the flex database')

    # Map attributes to the glycowork graph nodes
    flex_attribute_mapping = {
        mapper[gcontact_node]: attributes
        for gcontact_node, attributes in node_attributes.items()
        if gcontact_node in mapper
    }

    # Assign the mapped attributes to the glycowork graph
    nx.set_node_attributes(glycowork_graph, flex_attribute_mapping)
    return glycowork_graph

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
                        #matching_monosaccharides.append(attributes.get("string_labels", ""))
                        matching_monosaccharides.append(attributes.get('Monosaccharide', ""))
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
            #if glycan_class == "":
             #   glycan_class = np.nan

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
