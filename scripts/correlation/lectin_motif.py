from scripts.correlation.func.metric_df4 import metric_df
from scripts.correlation.func.plot_corr_regg import plot_combined_colors

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

set3= {
    "AAL": ["Fuc"],
    "SNA": ["Sia(a2-6)"],
    "ConA": ["Man"],
    "MAL-II": ["Neu5Ac(a2-3)"],
    "PNA": ["Gal(b1-3)GalNAc"],
    "CMA": ["Fuc(a1-2)Gal", "GalNAc"], # specifiy th terminal node within is teh first monosac if there is more then one
 #"CF": ["GalNAc(a1-?)", "GlcNAc(b1-?)"], # No lectin found
    "AOL": ["Fuc"],
    "HPA": ["GalNAc(a1-?)", "GlcNAc(b1-?)"]
#"LAA": ["Fuc(a1-2)Gal(b1-4)GlcNAc"]#, No lectin found
                        }

#Set4 #Add AOL later on. new glycan data needs df adjustment
lectin_binding_motif = {
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


metric_df_ = {}
for lectin, properties in lectin_binding_motif.items():

    print(f"Processing lectin: {lectin}")
    metric_df_[lectin] = metric_df(lectin, properties)
    #plot_combined_colors(metric_df_[lectin], lectin, properties["motif"])

