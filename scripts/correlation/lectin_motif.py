from scripts.correlation.func.metric_df import metric_df, perform_mediation_analysis_with_class
from scripts.correlation.func.plot_corr_regg import plot_combined_colors, plot_separate_class,visualize_mediation_results_with_class
from scripts.correlation.func.plot_corr_regg import plot_Binding_vs_Flexibility_and_SASA_with_stats

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
        "motif": ["Man(a1-?)"],
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
    },
    "CMA": {
        "motif": ["Fuc(a1-2)Gal", "GalNAc"],
        "termini_list": [["t", "f"], ["t"]]
    },
    "HPA": {
        "motif": ["GalNAc(a1-?)", "GlcNAc(b1-?)"],
        "termini_list": [["t"], ["t"]]
    },




    "AMA": {
        "motif": ["Man(a1-3)-Man(a1-8)"],
        "termini_list": [["t"], ["t"]] #chitobiose core?
    },
    "GNA": {
        "motif": ["Man(a1-6)", "Man(a1-3)"], #or Man(a1-?)
        "termini_list": [["t"], ["t"]]
    },
    "HHL": {
        "motif": ["Man"],
        "termini_list": [["t"]]
    },
    "MNA-M": {
        "motif": ["Man(a1-3)", "Man(a1-6)"],
        "termini_list": [["t"], ["t"]]
    },
    "NPA": {
        "motif": ["Man(a1-6)", "Man(a1-3)"],
        "termini_list": [["t"], ["t"]]
    },
    "UDA": {
        "motif": ["Man(a1-6)"],
        "termini_list": [["t"]]
    },
    "ABA": {
        "motif": ["Gal(b1-4)GlcNAc"],
        "termini_list": [["t", "f"]]
    },
    "CA": {
        "motif": ["Gal(b1-4)GlcNAc"],
        "termini_list": [["t", "f"]]
    },
    "CAA": {
        "motif": ["Gal(b1-4)GlcNAc"],
        "termini_list": [["t", "f"]]
    },
    "RPA": {
        "motif": ["Fuc"],
        "termini_list": [["t"]]
    },
    "TL": {
        "motif": ["Fuc"],
        "termini_list": [["t"]]
    },
    "ACA": {
        "motif": ["Sia(a2-?)"],
        "termini_list": [["t"]]
    },
    "AIA": {
        "motif": ["Gal", "GlcNAc"],
        "termini_list": [["t"], ["t"]]
    },
    "CF": {
        "motif": ["GalNAc"],
        "termini_list": [["t"]]
    },
    "HAA": {
        "motif": ["GalNAc"],
        "termini_list": [["t"]]
    },
    "MPA": {
        "motif": ["Gal"],
        "termini_list": [["t"]]
    },
    "LAA": {
        "motif": ["Fuc"],
        "termini_list": [["t"]]
    },
    "LcH": {
        "motif": ["Fuc", "Man(a1-2)"],
        "termini_list": [["f"], ["t"]]
    },
    "LTL": {
        "motif": ["Fuc"],
        "termini_list": [["t"]]
    },
    "LTA": {
        "motif": ["Fuc"],
        "termini_list": [["t"]]
    },
    "PSA": {
        "motif": ["Fuc(a1-6)"],
        "termini_list": [["t"]]
    },
    "PTL-I": {
        "motif": ["Fuc", "GalNAc", "Gal"],
        "termini_list": [["t"], ["t"], ["t"]]
    },
    "PTA-II": {
        "motif": ["Fuc", "GalNAc", "Gal"],
        "termini_list": [["t"], ["t"], ["t"]]
    },
    "PTL-II": {
        "motif": ["Fuc"],
        "termini_list": [["t"]]
    },
    "TJA-II": {
        "motif": ["Fuc", "GalNAc"],
        "termini_list": [["t"], ["t"]]
    },
    "UEA-I": {
        "motif": ["Fuc"],
        "termini_list": [["t"]]
    },
    "CTB": {
        "motif": ["Gal", "Fuc(a1-2)"],
        "termini_list": [["t"], ["t"]]
    },
    "MAL-I": {
        "motif": ["Gal"],
        "termini_list": [["t"]]
    },
    "MAA": {
        "motif": ["Gal"],
        "termini_list": [["t"]]
    },
    "MAL": {
        "motif": ["Gal"],
        "termini_list": [["t"]]
    },
    "PSL": {
        "motif": ["Sia(a2-6)"],
        "termini_list": [["f"]]
    },
    "TJA-I": {
        "motif": ["Sia(a2-6)"],
        "termini_list": [["t"]]
    },
    "GS-II": {
        "motif": ["GlcNAc"],
        "termini_list": [["t"]]
    },
    "PWA": {
        "motif": ["GlcNAc"],
        "termini_list": [["t"]]
    },
    "UEA-II": {
        "motif": ["GlcNAc(b1-3)Gal", "GlcNAc(b1-3)GalNAc"],
        "termini_list": [["t", "f"], ["t", "f"]]
    },
    "WGA": {
        "motif": ["GlcNAc"],
        "termini_list": [["t"]]
    },
    "BPA": {
        "motif": ["Gal", "GalNAc"],
        "termini_list": [["t"], ["t"]]
    },
    "BPL": {
        "motif": ["Gal", "GalNAc"],
        "termini_list": [["t"], ["t"]]
    },
    "ECA": {
        "motif": ["LacNAc"],
        "termini_list": [["t"]]
    },
    "GS-I": {
        "motif": ["Gal", "GalNAc"],
        "termini_list": [["t"], ["t"]]
    },
    "LEA": {
        "motif": ["LacNAc"],
        "termini_list": [["t"]]
    },
    "LEL": {
        "motif": ["LacNAc"],
        "termini_list": [["t"]]
    },
    "MOA": {
        "motif": ["Gal(a1-3)Gal", "Gal(a1-3)GalNAc"],
        "termini_list": [["t", "f"], ["t", "f"]]
    },
    "PA-IL": {
        "motif": ["Gal"],
        "termini_list": [["t"]]
    },
    "LecA": {
        "motif": ["Gal"],
        "termini_list": [["t"]]
    },
    "RCA-I": {
        "motif": ["LacNAc"],
        "termini_list": [["t"]]
    },
    "RCA120": {
        "motif": ["LacNAc"],
        "termini_list": [["t"]]
    },
    "SJA": {
        "motif": ["LacNAc"],
        "termini_list": [["t"]]
    },
    "STA": {
        "motif": ["LacNAc"],
        "termini_list": [["t"]]
    },
    "STL": {
        "motif": ["LacNAc"],
        "termini_list": [["t"]]
    },
    "CSA": {
        "motif": ["GalNAc"],
        "termini_list": [["t"]]
    },
    "DBA": {
        "motif": ["GalNAc"],
        "termini_list": [["t"]]
    },
    "SBA": {
        "motif": ["GalNAc"],
        "termini_list": [["t"]]
    },
    "VVL": {
        "motif": ["LacdiNAc", "GalNAc"],
        "termini_list": [["t"], ["t"]]
    },
    "VVA": {
        "motif": ["LacdiNAc", "GalNAc"],
        "termini_list": [["t"], ["t"]]
    },
    "WFA": {
        "motif": ["GalNAc"],
        "termini_list": [["t"]]
    }
}

metric_df_ = {}

for lectin, properties in lectin_binding_motif.items():
    print(f"\nProcessing lectin: {lectin}")
    print(f"Motif: {properties['motif']}")
    print(f"Termini: {properties['termini_list']}")

    # Generate metric DataFrame
    df = metric_df(lectin, properties)

    if df.empty:
        print(f"⚠️ Skipping {lectin} due to missing binding data.")
        continue

    metric_df_[lectin] = df
    print(f"Skipping unknown lectin: {lectin}")
    #plot_combined_colors(metric_df_[lectin], lectin, properties["motif"]) #regression without stats
    #plot_separate_class(metric_df_[lectin], lectin, properties["motif"]) #regression per glycan class without stats

    plot_Binding_vs_Flexibility_and_SASA_with_stats(metric_df_[lectin], lectin, properties["motif"])

    """effects = perform_mediation_analysis_with_class(
    metric_df=metric_df_[lectin],
    independent_var='weighted_mean_flexibility',
    class_var='class',
    dependent_var='binding_score')

    visualize_mediation_results_with_class(
        metric_df=metric_df_[lectin],
        lectin=lectin,
        binding_motif=properties["motif"],
        #independent_var='weighted_mean_flexibility',
        independent_var='SASA_weighted_sum',
        class_var='class',
        dependent_var='binding_score',
        effects=effects
    )
"""

