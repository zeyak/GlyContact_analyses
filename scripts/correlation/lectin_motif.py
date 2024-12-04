from scripts.correlation.func.metric_df3 import metric_df
from scripts.correlation.func.plot_corr_regg import plot_combined, plot_corr_binding_SASA_subplots
from scripts.correlation.func.plot_corr_regg import plot_corr_binding_SASA_subplots

set1 = {
    "AAL": ["Fuc"],
    "SNA": ["Neu5Ac(a2-6)", "Neu5Gc(a2-6)"],
     "ConA": ["Man(a1-2)"],
     "MAL-II": ["Neu5Ac"],
     "PNA": ["Gal(b1-3)", "GalNAc"] #this returns empty
}

set2= {
    "AAL": ["Fuc"],
    "SNA": ["Neu5Ac(a2-6)", "Neu5Gc(a2-6)"],
    "ConA": ["Man"],
    "MAL-II": ["Neu5Ac(a2-3)"],
    # check how many neigh and if it has more than one.
    "PNA": ["Gal", "GalNAc"]}

#Set3
lectin_binding_motif = {
    "AAL": ["Fuc"],
    "SNA": ["Neu5Ac(a2-6)", "Neu5Gc(a2-6)"],
    "ConA": ["Man"],
    "MAL-II": ["Neu5Ac(a2-3)"],
"PNA": ["Gal(b1-3)GalNAc"],
"CMA": ["Fuc(a1-2)Gal", "GalNAc"],
# "CF": ["GalNAc(a1-?)", "GlcNAc(b1-?)"], #subgraph isomorphism
# "HPA": ["GalNAc(a1-?)", "GlcNAc(b1-?)"],
#"LAA": ["Fuc(a1-2)Gal(b1-4)GlcNAc"]
                        }

for lectin, binding_motif in lectin_binding_motif.items():
    metric_df_instance = metric_df(lectin, binding_motif)
    #plot_combined(metric_df_instance, lectin, binding_motif)
    plot_corr_binding_SASA_subplots(metric_df_instance, lectin, binding_motif)
