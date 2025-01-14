from scripts.correlation.func.metric_df import metric_df, perform_mediation_analysis_with_class
from scripts.correlation.func.plot_corr_regg import plot_combined_colors, plot_separate_class,visualize_mediation_results_with_class


lectin_binding_motif = {
    "AOL": { "motif": ["Fuc"],
             "termini_list": [["t"]] },

    "AAL": {
        "motif": ["Fuc"],
        "termini_list": [["t"]]
    },
    "SNA": {
        "motif": ["Sia(a2-6)"],
        "termini_list": [["t"]]
    },
    "ConA": {
        "motif": ["Man"],
        "termini_list": [["t"]]
    },
    "MAL-II": {
        "motif": ["Sia(a2-3)"],
        "termini_list": [["t"]]
    },
    "PNA": {
        "motif": ["Gal(b1-3)GalNAc"],
        "termini_list": [["t", "f"]]
    },
    "CMA": {
        "motif": ["Fuc(a1-2)Gal", "GalNAc"],
        "termini_list": [["t", "f"], ["t"]]
    },
    "HPA": {
        "motif": ["GalNAc(a1-?)", "GlcNAc(b1-?)"],
        "termini_list": [["t"], ["t"]]
    }
}


metric_df_ = {}
for lectin, properties in lectin_binding_motif.items():
    print("")
    print(f"Processing lectin: {lectin}")
    print(f"Motif: {properties['motif']}")
    print(f"Termini: {properties['termini_list']}")

    metric_df_[lectin] = metric_df(lectin,properties)
    plot_combined_colors(metric_df_[lectin], lectin, properties["motif"])
    plot_separate_class(metric_df_[lectin], lectin, properties["motif"])

"""    effects = perform_mediation_analysis_with_class(
    metric_df=metric_df_[lectin],
    independent_var='weighted_mean_flexibility',
    class_var='class',
    dependent_var='binding_score')

    visualize_mediation_results_with_class(
        metric_df=metric_df_[lectin],
        lectin=lectin,
        binding_motif=properties["motif"],
        independent_var='weighted_mean_flexibility',
        #independent_var='SASA_weighted_sum',
        class_var='class',
        dependent_var='binding_score',
        effects=effects
    )

"""





