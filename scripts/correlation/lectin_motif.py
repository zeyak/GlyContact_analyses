from scripts.correlation.func.metric_df2 import metric_df
from scripts.correlation.func.plot_corr_regg import plot_combined

# Set2
lectin_binding_motif = {
    "AAL": ["Fuc"],
    "SNA": ["Neu5Ac(a2-6)", "Neu5Gc(a2-6)"],
    # "ConA": ["Man(a1-2)"],
    "ConA": ["Man"],
    # "MAL-II": ["Neu5Ac"],
    "MAL-II": ["Neu5Ac(a2-3)"],
    # "PNA": ["Gal(b1-3)", "GalNAc"] #this returns empty
    # check how many neigh and if it has more than one.
    "PNA": ["Gal", "GalNAc"],
    "CMA": ["Fuc(a1-2)Gal", "GalNAc"],
    "AOL": ["Fuc"],
    # "CF": ["GalNAc(a1-?)", "GlcNAc(b1-?)"], #subgraph isomorphism
    "HPA": ["GalNAc(a1-?)", "GlcNAc(b1-?)"],
    "LAA": ["Fuc(a1-2)Gal(b1-4)GlcNAc"]

}

for lectin, binding_motif in lectin_binding_motif.items():
    metric_df_instance = metric_df(lectin, binding_motif)
    plot_combined(metric_df_instance, lectin, binding_motif)
