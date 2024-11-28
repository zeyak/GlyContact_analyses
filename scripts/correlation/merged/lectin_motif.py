import pandas as pd
from glycowork.motif.graph import compare_glycans
import pickle


lectin = "AAL"
binding_motif = ["Fuc"]

lectin_binding_motif ={
    "AAL": ["Fuc"],
    "SNA": ["Neu5Ac", "Neu5Gc"]
    "MAL-II": ["Neu5Ac"],


}

SNA; Neu5Ac(a2-6)
MAL-II; Neu5Ac(a2-3)
ConA; Man(a1-2)
AAL; Fuc
PNA; Gal(b1-3)GalNAc