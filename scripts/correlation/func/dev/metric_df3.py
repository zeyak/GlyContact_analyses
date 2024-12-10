import pickle
import pandas as pd
from glycowork.motif.graph import compare_glycans, subgraph_isomorphism
from glycowork.motif.processing import get_class
from glycowork.motif.graph import *
from glycowork.motif.draw import GlycoDraw

GlycoDraw("Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-3)[Neu5Gc(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)][GlcNAc(b1-4)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc", highlight_motif = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc")


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
            return flex_glycan
    return None


def compute_sasa_metrics(nodes):
    """Compute sum, mean, and max SASA metrics (Mean Score, Median Score, and Weighted Score)."""
    if not nodes:
        return {
            "SASA_mean": None,
            "SASA_median": None,
            "SASA_weighted": None,
            "SASA_mean_max": None,
            "SASA_median_max": None,
            "SASA_weighted_max": None,
            "SASA_mean_sum": None,
            "SASA_median_sum": None,
            "SASA_weighted_sum": None,
        }

    SASA_mean_sum = sum(node.get("Mean Score", 0) for node in nodes)
    SASA_median_sum = sum(node.get("Median Score", 0) for node in nodes)
    SASA_weighted_sum = sum(node.get("Weighted Score", 0) for node in nodes)

    SASA_mean = SASA_mean_sum / len(nodes)
    SASA_median = SASA_median_sum / len(nodes)
    SASA_weighted = SASA_weighted_sum / len(nodes)

    SASA_mean_max = max(node.get("Mean Score", 0) for node in nodes)
    SASA_median_max = max(node.get("Median Score", 0) for node in nodes)
    SASA_weighted_max = max(node.get("Weighted Score", 0) for node in nodes)

    return {
        "SASA_mean": SASA_mean,
        "SASA_median": SASA_median,
        "SASA_weighted": SASA_weighted,
        "SASA_mean_max": SASA_mean_max,
        "SASA_median_max": SASA_median_max,
        "SASA_weighted_max": SASA_weighted_max,
        "SASA_mean_sum": SASA_mean_sum,
        "SASA_median_sum": SASA_median_sum,
        "SASA_weighted_sum": SASA_weighted_sum,
    }


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

    def extract_terminal_residue(binding_motif):
        """
        Extract the terminal residue from a binding motif.

        Args:
            binding_motif: A string representing a glycan binding motif (e.g., "Fuc(a1-2)Gal").

        Returns:
            A string representing the terminal residue (e.g., "Fuc(a1-2)").
        """
        # Find the first `(` and include everything up to the matching `)`
        if '(' in binding_motif and ')' in binding_motif:
            end_idx = binding_motif.find(')') + 1
            return binding_motif[:end_idx]  # Extract up to the first `)`
        else:
            return binding_motif  # Return as is if no linkage format exists

    # Preprocess the binding motifs list to extract terminal residues
    terminal_binding_motifs = [extract_terminal_residue(motif) for motif in binding_motifs]

    matching_monosaccharides = []  # Nodes matching the binding motifs
    flexibility_values = []  # Weighted flexibility values for matching nodes
    found_motifs = []  # List of motifs found in the glycan

    # Iterate over the terminal binding motifs
    for binding_motif in terminal_binding_motifs:
        is_present, matched_nodes = subgraph_isomorphism(
            matched_flex_glycan,
            binding_motif,
            return_matches=True,
            termini_list=['terminal', 'flexible']  # Restrict to terminal matches
        )

        if not is_present:
            # If the current motif is not present, skip to the next
            continue

        # Add the current motif to the list of found motifs
        found_motifs.append(binding_motif)

        # Flatten matched nodes (in case of nested lists)
        if matched_nodes:  # Check if there are matches
            if isinstance(matched_nodes[0], list):
                matched_nodes = [node for sublist in matched_nodes for node in sublist]

        # Process nodes in the graph to extract attributes for matched nodes
        for node in matched_nodes:
            if node in graph.nodes:
                attributes = graph.nodes[node]
                #matching_monosaccharides.append(attributes.get("Monosaccharide", 0))
                matching_monosaccharides.append(attributes)
                flexibility_values.append(attributes.get("weighted_mean_flexibility", 0))

                # Debugging: print node data and motif match
                print(f"Node {node} attributes: {attributes}")
                if "name" in attributes:
                    print(f"Matched motif: {binding_motif} in graph node {node} with name {attributes['name']}")
                    print(f"Matched motif: {binding_motif} in graph node {node} with Monosaccharide {attributes['Monosaccharide']}")

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
        overall_flexibility = (sum(flexibility_values) / len(flexibility_values) if flexibility_values else None)

        # Get the glycan's class
        glycan_class = get_class(matched_flex_glycan)

        metric_data.append({
            "glycan": glycan,
            "binding_score": binding_score,
            #"SASA_mean": sasa_metrics["SASA_mean"],
            #"SASA_median": sasa_metrics["SASA_median"],
            "SASA_weighted": sasa_metrics["SASA_weighted"],
            #"SASA_mean_max": sasa_metrics.get("SASA_mean_max", None),
            #"SASA_median_max": sasa_metrics.get("SASA_median_max", None),
            "SASA_weighted_max": sasa_metrics.get("SASA_weighted_max", None),
            #"SASA_mean_sum": sasa_metrics.get("SASA_mean_sum", None),
            #"SASA_median_sum": sasa_metrics.get("SASA_median_sum", None),
            "SASA_weighted_sum": sasa_metrics.get("SASA_weighted_sum", None),
            "weighted_mean_flexibility": overall_flexibility,
            "class": glycan_class,  # Add class column
        })

    metric_df = pd.DataFrame(metric_data)
    return metric_df , binding_df, filtered_df, flex_data


def subgraph_isomorphism_II(glycan, motif, termini_list = [], count = False, return_matches = False):
  if isinstance(glycan, str):
    glycan = glycan_to_nxGraph(glycan, termini='calc' if termini_list else None)
  else:
    g1 = deepcopy(glycan)
  if isinstance(motif, str):
    motif = glycan_to_nxGraph(motif, termini='provided' if termini_list else None, termini_list=termini_list)
  else:
    g2 = deepcopy(motif)
  if isinstance(glycan, str) and isinstance(motif, str):
    if motif.count('(') > glycan.count('('):
      return (0, []) if return_matches else 0 if count else False
    if not count and not return_matches and motif in glycan:
      return True
    motif_comp = min_process_glycans([motif, glycan])
    if 'O' in glycan + motif:
      glycan, motif = [re.sub(r"(?<=[a-zA-Z])\d+(?=[a-zA-Z])", 'O', g).replace('NeuOAc', 'Neu5Ac').replace('NeuOGc', 'Neu5Gc') for g in [glycan, motif]]
  else:
    print(motif.__dict__,glycan.__dict__)
    # if len(glycan.nodes) < len(motif.nodes):
    #   return (0, []) if return_matches else 0 if count else False
    motif_comp = [nx.get_node_attributes(motif, "string_labels").values(), nx.get_node_attributes(glycan, "string_labels").values()]
    if 'O' in ''.join(unwrap(motif_comp)):
      g1, g2 = ptm_wildcard_for_graph(deepcopy(glycan)), ptm_wildcard_for_graph(deepcopy(motif))
    else:
      g1, g2 = deepcopy(glycan), motif
  narrow_wildcard_list = {k: get_possible_linkages(k) if '?' in k else get_possible_monosaccharides(k) for k in set(unwrap(motif_comp))
                          if '?' in k or k in {'Hex', 'HexOS', 'HexNAc', 'HexNAcOS', 'dHex', 'Sia', 'HexA', 'Pen', 'Monosaccharide'} or '!' in k}
  if termini_list or narrow_wildcard_list:
    graph_pair = nx.algorithms.isomorphism.GraphMatcher(g1, g2, node_match = categorical_node_match_wildcard('string_labels', 'unknown', narrow_wildcard_list,
                                                                                                             'termini', 'flexible'))
  else:
    g1_node_attr = set(nx.get_node_attributes(g1, "string_labels").values())
    if all(k in g1_node_attr for k in motif_comp[0]):
      graph_pair = nx.algorithms.isomorphism.GraphMatcher(g1, g2, node_match = nx.algorithms.isomorphism.categorical_node_match('string_labels', 'unknown'))
    else:
      return (0, []) if return_matches else 0 if count else False

  if return_matches:
    mappings = [list(isos.keys()) for isos in graph_pair.subgraph_isomorphisms_iter()]

  # Count motif occurrence
  if count:
    counts = 0
    while graph_pair.subgraph_is_isomorphic():
      mapping = graph_pair.mapping
      inverse_mapping  = {v: k for k, v in mapping.items()}
      if all(inverse_mapping[node] < inverse_mapping[neighbor] for node, neighbor in g2.edges()):
        counts += 1
      g1.remove_nodes_from(mapping.keys())
      if termini_list or narrow_wildcard_list:
        graph_pair = nx.algorithms.isomorphism.GraphMatcher(g1, g2, node_match = categorical_node_match_wildcard('string_labels', 'unknown', narrow_wildcard_list,
                                                                                                             'termini', 'flexible'))
      else:
        g1_node_attr = set(nx.get_node_attributes(g1, "string_labels").values())
        if all(k in g1_node_attr for k in motif_comp[0]):
          graph_pair = nx.algorithms.isomorphism.GraphMatcher(g1, g2, node_match = nx.algorithms.isomorphism.categorical_node_match('string_labels', 'unknown'))
        else:
          return counts if not return_matches else (counts, mappings)
    return counts if not return_matches else (counts, mappings)
  else:
    if graph_pair.subgraph_is_isomorphic():
      for mapping in graph_pair.subgraph_isomorphisms_iter():
        mapping = {v: k for k, v in mapping.items()}
        for node, neighbor in g2.edges():
          if mapping[node] >= mapping[neighbor]:
             return False if not return_matches else (0, [])
        return True if not return_matches else (1, mappings)
  return False if not return_matches else (0, [])


# Example usage
lectin = "CMA"
binding_motif = ["Fuc(a1-2)Gal", "GalNAc"]

metric_df, binding_df, filtered_df, flex_data = metric_df(lectin, binding_motif)

subgraph_isomorphism_II(
            flex_data['Fuc(a1-2)Gal'],
            'Fuc(a1-2)',
            return_matches=True,
            termini_list=['terminal'])
