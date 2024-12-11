from scripts.correlation.func.metric_df4 import metric_df
from scripts.correlation.func.plot_corr_regg import plot_combined_colors


lectin_binding_motif = {
    "AOL": ["Fuc"],
    "AAL": ["Fuc"],
    "SNA": ["Sia(a2-6)"],
    "ConA": ["Man"],
    "MAL-II": ["Sia(a2-3)"],
    "PNA": ["Gal(b1-3)GalNAc"],
    "CMA": ["Fuc(a1-2)Gal", "GalNAc"],
    "HPA": ["GalNAc(a1-?)", "GlcNAc(b1-?)"]
}


metric_df_ = {}
for lectin, binding_motif in lectin_binding_motif.items():
    print("")
    print(f"Processing lectin: {lectin}")

    metric_df_[lectin] = metric_df(lectin,binding_motif)
    plot_combined_colors(metric_df_[lectin], lectin, binding_motif)

