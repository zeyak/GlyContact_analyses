from scripts.correlation.func.metric_df import metric_df
from scripts.correlation.func.plot_corr_regg import plot_combined_colors


lectin_binding_motif_ = {
    "AAL": {
        "motif": ["Fuc"],
        "terminal": ["Fuc"],
        "flexible": []
    },
    "SNA": {
        "motif": ["Sia(a2-6)"],
        "terminal": ["Sia"],
        "flexible": []
    },
    "ConA": {
        "motif": ["Man"],
        "terminal": ["Man"],
        "flexible": []
    },
    "MAL-II": {
        "motif": ["Sia(a2-3)"],
        "terminal": ["Sia"],
        "flexible": []
    },
    "PNA": {
        "motif": ["Gal(b1-3)GalNAc"],
        "terminal": ["Gal"],
        "flexible": ["GalNAc"]
    },
    "CMA": {
        "motif": ["Fuc(a1-2)Gal", "GalNAc"],
        "terminal": ["Fuc", "GalNAc"],
        "flexible": ["Gal" ]
    },

    "HPA": {
        "motif": ["GalNAc(a1-?)", "GlcNAc(b1-?)"],
        "terminal": ["GalNAc", "GlcNAc"],
        "flexible": []
    }
}

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
#for lectin, properties in lectin_binding_motif.items():
for lectin, binding_motif in lectin_binding_motif.items():
    print("")
    print(f"Processing lectin: {lectin}")

    metric_df_[lectin] = metric_df(lectin,binding_motif)
    plot_combined_colors(metric_df_[lectin], lectin, binding_motif)

