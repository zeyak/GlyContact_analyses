from scripts.correlation.func.metric_df import metric_df
from scripts.correlation.func.plot_corr_regg import plot_combined


lectin_binding_motif = {
    "AAL": ["Fuc"],
     "SNA": ["Neu5Ac(a2-6)", "Neu5Gc(a2-6)"],
     "MAL-II": ["Neu5Ac"],
     "ConA": ["Man(a1-2)"],
     "PNA": ["Gal(b1-3)", "GalNAc"]
}


for lectin, binding_motif in lectin_binding_motif.items():
    metric_df_instance = metric_df(lectin, binding_motif)
    plot_combined(metric_df_instance, lectin, binding_motif)


