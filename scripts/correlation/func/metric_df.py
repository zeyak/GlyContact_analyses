import pickle
from typing import Dict, Tuple, List
import numpy as np
import pandas as pd
import pingouin as pg
from glycowork.motif.graph import compare_glycans, subgraph_isomorphism
from glycowork.motif.graph import  glycan_to_nxGraph
from glycowork.motif.processing import get_class
import sys

from joblib.externals.loky.backend.queues import Queue
from networkx.classes import Graph
import networkx as nx
from sklearn.preprocessing import LabelEncoder

import importlib.resources as pkg_resources
import glycontact

BINDING_DATA_PATH_v3 = 'data/20250216_glycan_binding.csv'
BINDING_DATA_PATH_v2 = 'data/20241206_glycan_binding.csv' #seq,protein
BINDING_DATA_PATH_v1 = 'data/glycan_binding.csv' #seq,protein

def load_data_pdb():
    """Load glycan flexibility data from PDB source."""
    with pkg_resources.files(glycontact).joinpath("glycan_graphs.pkl").open("rb") as f:
        return pickle.load(f)


def load_data_():
    """Load glycan flexibility and binding data, process graphs, and return results."""
    flex_data = load_data_pdb()
    #binding_df = pd.read_csv(BINDING_DATA_PATH_v2)
    invalid_graphs = [glycan for glycan in flex_data if not isinstance(flex_data[glycan], nx.Graph)]

    binding_df_v2 = pd.read_csv(BINDING_DATA_PATH_v2)
    binding_df_v3 = pd.read_csv(BINDING_DATA_PATH_v3)

    # Combine df1 and df2, ensuring no duplicate target rows if applicable
    df_combined = pd.concat([binding_df_v2, binding_df_v3]).drop_duplicates(subset=['target'])

    # Merge df3 with the combined df on target
    binding_df = pd.merge(binding_df_v3, df_combined[['target', 'protein']], how="left", on="target")

    return flex_data, binding_df, invalid_graphs



def load_data():
    """Load glycan flexibility and binding data, process graphs, and return results."""
    flex_data = load_data_pdb()
    invalid_graphs = [glycan for glycan in flex_data if not isinstance(flex_data[glycan], nx.Graph)]

    def map_protein_to_target(df_target, df_map):
        """
        Maps protein names to their corresponding targets (sequences) in df_target
        using mapping from df_map.

        Args:
            df_target (pd.DataFrame): DataFrame that needs the target column updated.
            df_map (pd.DataFrame): DataFrame containing the protein-to-target mapping.

        Returns:
            pd.DataFrame: Updated df_target with mapped target values.
        """
        # Create a mapping dictionary {target -> protein}
        target_to_protein = dict(zip(df_map["target"], df_map["protein"]))

        # Apply mapping to create the "protein" column in df_target
        df_target["protein"] = df_target["target"].map(target_to_protein)

        return df_target


    binding_df_v2 = pd.read_csv(BINDING_DATA_PATH_v2)
    binding_df_v3 = pd.read_csv(BINDING_DATA_PATH_v3)

    binding_df = map_protein_to_target(binding_df_v3, binding_df_v2)

    return flex_data, binding_df, invalid_graphs


def filter_binding_data(binding_df: pd.DataFrame, lectin: str) -> pd.DataFrame:
    """Filter the binding DataFrame for the given lectin."""
    filtered_df = binding_df[binding_df.iloc[:, -1].eq(lectin)]
    filtered_df = filtered_df.dropna(axis=1, how='all')  # Drop columns with all NaN values
    return filtered_df

def get_glycan_scores(filtered_df: dict[str, float]) -> Dict[str, float]:
    """Calculate mean binding scores for glycans."""
    lectin_df = filtered_df.drop(columns=["target", "protein"]) # Exclude "protein" and "target" columns
    glycan_scores = lectin_df.mean(axis=0).to_dict()

    return glycan_scores

def find_matching_glycan(flex_data , glycan):
    """Find the matching glycan in flex_data."""
    for flex_glycan in flex_data.keys():
        if compare_glycans(glycan, flex_glycan):
            return glycan
    return None


def compute_overall_SASA(SASA_values):
    return sum(SASA_values) / len(SASA_values) if SASA_values else None

def compute_overall_flexibility(flexibility_values ):
    return sum(flexibility_values) / len(flexibility_values) if flexibility_values else None

def compute_overall_Q(Q_values):
    return sum(Q_values) / len(Q_values) if Q_values else None

def compute_overall_theta(theta_values):
    return sum(theta_values) / len(theta_values) if theta_values else None


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
    matching_monosaccharides, SASA, flexibility, Q, theta, conformation, found_motifs = [], [], [], [], [], [], []

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
                        SASA.append(attributes.get("SASA", 0))
                        flexibility.append(attributes.get("flexibility", 0))
                        Q.append(attributes.get("Q", 0))
                        theta.append(attributes.get("theta", 0))
                        conformation.append(attributes.get("conformation", 0))
                        

                    print(f"Matching monosaccharides: {matching_monosaccharides}")
                    print(f"SASA scores: {SASA}")
                    print(f"Flexibility: {flexibility}")
                    print(f"Q: {Q}")
                    print(f"theta: {theta}")
                    print(f"conformation: {conformation}")

                    print("")
                except Exception as e:
                    print(f"Error extracting attributes for node {mono} in glycan {matched_glycan}: {e}")
        else:
            print(f"Skipping invalid graph or graph with no nodes for glycan: {matched_glycan}")

    return matching_monosaccharides, SASA, flexibility, Q, theta, conformation, found_motifs



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
        matching_monosaccharides, SASA, flexibility, Q, theta, conformation, found_motifs = process_glycan_with_motifs(
            matched_glycan, properties, flex_data)

        # Skip empty monosaccharides
        matching_monosaccharides = [m for m in matching_monosaccharides if m.strip()]

        if matching_monosaccharides:
            overall_SASA = compute_overall_SASA(SASA)
            overall_flexibility = compute_overall_flexibility(flexibility)
            overall_Q = compute_overall_Q(Q)
            overall_theta = compute_overall_theta(theta)
            glycan_class = get_class(matched_glycan) or np.nan

            metric_data.append({
                "glycan": glycan,
                "binding_score": binding_score,
                "SASA": overall_SASA,
                "flexibility": overall_flexibility,
                "Q": overall_Q,
                "theta": overall_theta,
                "conformation": conformation,
                "monosaccharides": matching_monosaccharides,
                "motifs": found_motifs,
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
    if 'glycan' in metric_df.columns:
        metric_df.set_index('glycan', inplace=True)
    else:
        print(f"⚠️ Warning: 'glycan' column missing for {lectin}. Skipping index setting.")

    metric_df.to_excel(f'scripts/correlation/metric_df/{lectin}_metrics.xlsx', index=True, header=True)
    return metric_df

sys.stdout = open('scripts/correlation/metric_df/metric_df.log', 'w')

def perform_mediation_analysis_with_class(metric_df, independent_var, class_var, dependent_var):
    """
    Perform mediation analysis to test if glycan class mediates the relationship between
    the independent variable and the dependent variable. Handles NaN values and single-class cases.

    Parameters:
        metric_df (pd.DataFrame): The DataFrame containing metrics for mediation analysis.
        independent_var (str): Name of the independent variable (e.g., 'flexibility').
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
