#################
# Cell lineages #
#################

# msk_colon
def classify_vectra_colon(activity):
    """ Returns the cell type based on the activity of the markers.
    activity: dict
        Dictionary of marker activity.
    """
    if activity["Foxp3_gt"]==1 or activity["ICOS_gt"]==1 or activity["CD8_gt"]==1 or activity["CD3_gt"]==1:
        out = "Lymphocytes"
    elif activity["panCK+CK7+CAM5.2_gt"]==1:
        out = "Epithelial"
    else:
        out = "Pan-Negative"
    return out

# msk_pancreas
def classify_vectra_pancreas(activity):
    """ Returns the cell type based on the activity of the markers.
    activity: dict
        Dictionary of marker activity.
    """
    if activity["CD40-L"]==1 or activity["PD-1"]==1 or activity["CD8"]==1:
        out = "Lymphocytes"
    elif activity["CD40"]==1 or activity["PD-L1"]==1:
        out = "Other Immune"
    elif activity["panCK"]==1:
        out = "Epithelial"
    else:
        out = "Pan-Negative"
    return out

# Cell Type
codex_colon_assignment = {
    "Lymphocytes": ["CD8+ T", "CD4+ T cell", "B", "NK", "CD7+ Immune", "Plasma"],
    "Myeloids": ["Neutrophil", "M1 Macrophage", "M2 Macrophage"],
    "Other Immune": ["DC", "Paneth"],
    "Epithelial": ["Enterocyte", "MUC1+ Enterocyte", "CD66+ Enterocyte", "CD57+ Enterocyte"],
    "Stroma": ["Stroma"],
    "Cancer": [],
    "Vasculature": ["Endothelial", "Lymphatic"],
    "Muscle": ["Smooth muscle", "ICC",],
    "Precursors": ["TA", "Cycling TA"],
    "Other": ["Goblet", "Neuroendocrine", "Nerve", "Other"]
}
     
# lineage
mibi_decidua_assignment = {
    "Lymphocytes": ["NK1", "NKT",  "CD8T", "NK2", "NK3", "CD4T", "NK4", "Treg"],
    "Myeloids": ["Mac2a", "Mac1a", "Mac1b", "Mac2c", "Mac2b", "Placental_Mac",],
    "Other Immune": ["DC", "Mast"],
    "Epithelial": ["Glandular"],
    "Stroma": ["Fibroblasts", "Myofibroblasts"],
    "Cancer": [],
    "Vasculature": ["Endothelial"],
    "Muscle": ["muscle"],
    "Precursors": ["EVT1a", "EVT1b", "EVT2", "EVT1c"],
    "Other": ["other"],
}


# cell_meta_cluster
mibi_breast_assignment = {
    "Lymphocytes": ["Treg", "CD8T", "CD56", "CD4T", "CD20", "CD4T_HLADR", "CD4T_CD8T_dp", "CD3_DN"],
    "Myeloids": ["CD68", "CD68_CD163_DP", "CD163", "CD4_mono", "CD14"],
    "Other Immune": ["immune_other", "CD11c_HLADR", "ChyTry", "calprotectin"],
    "Epithelial": [],
    "Stroma": ["other_stroma_coll", "VIM", "FAP", "FAP_SMA", "other_stroma_fibronectin"],
    "Cancer": ["tumor_ecad", "tumor_other", "tumor_other_mono", "tumor_vim", "tumor_sma",
               "tumor_ck17", "tumor_CD56"],
    "Vasculature": ["CD31", "CD31_VIM"],
    "Muscle": ["SMA"],
    "Precursors": [],
    "Other": ["other"],
}

#######################
# Marker localization #
#######################
# nuclear, cytoplasm or ECM

marker_localization = {
    "nuclear": [
         "Ki67", "IDO", "SOX9", "Foxp3", "FOXP3", "FoxP3"
    ],
    "membrane": [
        "CD3", "CD8", "ICOS", "panCK+CK7+CAM5.2", "PD-L1","CD40", "CD40-L",
        "panCK", "PD-1", "PD-L1", "HLADR", "ECAD", "CK17", "CD163", "CD20", "CD3", 
        "CD4", "CD45", "CD56", "CD68", "CK7", "HLAG", "PDL1",  "CD16",
        "CD206", "Podoplanin", "CD117", "CD36", "aDefensin5", "Synaptophysin", "CD19", "CD161", 
        "CHGA", "Cytokeratin", "CD15", "CD21", "CD38", "CD34", "BCL2", "DCSIGN", "CD11c",
        "CD7", "CD90", "CD138", "CD66", "ChyTr", "Calprotectin", "CD14", "CD57",
    ],
    "ECM": [
        "VIM", 'Vimentin', "SMA", "aSMA", "Fibronectin", "Collagen1", "MUC2", "MUC1", "CD31",
        "FAP"
    ],
}

def reverse_dict(dictionary):
    rev_dict = {}
    for key in dictionary.keys():
        for ct in dictionary[key]:
            rev_dict[ct] = key
    return rev_dict