import pandas as pd
from pathlib import Path

amr_csv = Path("../../data/card_metadata_aro.csv")
output1 = Path("../data/card_metadata_aro_informativecolumns_minimal.csv")
output2 = Path("../data/card_metadata_aro_extended.csv")


# Define metagenome categories
categories = {
    "human": [
        "human gut metagenome", "human metagenome", "human oral metagenome", "human skin metagenome",
        "human feces metagenome", "human vaginal metagenome", "human nasopharyngeal metagenome", 
        "human lung metagenome", "human saliva metagenome", "human reproductive system metagenome", 
        "human urinary tract metagenome", "human eye metagenome", "human blood metagenome", 
        "human bile metagenome", "human tracheal metagenome", "human brain metagenome", 
        "human milk metagenome", "human semen metagenome", "human skeleton metagenome"
    ],
    "livestock": [
        "bovine gut metagenome", "bovine metagenome", "pig gut metagenome", "pig metagenome", 
        "chicken gut metagenome", "chicken metagenome", "sheep gut metagenome", "sheep metagenome"
    ],
    "marine": ["marine metagenome", "seawater metagenome"],
    "freshwater": ["freshwater metagenome", "lake water metagenome", "groundwater metagenome"],
    "soil": ["soil metagenome"],
    "wastewater": ["wastewater metagenome"]
}

categories_sp = {
    "livestock": [
        "Sus scrofa", "Sus scrofa domesticus", "Sus scrofa affinis", "Bos taurus", "Gallus gallus", "Equus caballus", 
        "Equs caballus", "Ovis aries", "Ovis", "Bos indicus", "Bos mutus", "Bos primigenius", "Bos frontalis",
        "Bos gaurus", "Gallus", "Capra hircus", "Capra aegagrus", "Capra ibex"
    ],
    "human": [
        "Homo sapiens"
    ]  
}

# Function to map organism to metagenome category
def get_category(organism: str) -> str | None:
    if pd.isna(organism):
        return None
    for category, values in categories.items():
        if organism in values:
            return category
    return None

# Function to map organism to metagenome category if organism is sp
def get_category_sp(organism: str) -> str | None:
    if pd.isna(organism):
        return None
    for category, values in categories_sp.items():
        if organism in values:
            return category
    return None


# Assign metagenome_category
# - Logic hierarchy:
#       1. when organism string already contains 'metagenome' → use get_category
#       2. else if librarysource matches the meta(genomic|transcriptomic) → use get_category_sp
#       3. Other metagenome categories → "other"
#    isolates get <NA>

# metagenome_category
def decide_category(row):
    if row["organism_type"] == "Isolate":
        return pd.NA

    # rule 1: organism string already includes 'metagenome'
    cat = get_category(row["organism"])
    if cat is not None:
        return cat

    # rule 2: librarysource categorisation (species names)
    if pd.notna(row["librarysource"]) and pd.Series(row["librarysource"]).str.contains(
        "METAGENOMIC|METATRANSCRIPTOMIC", case=False, regex=True
    ).bool():
        cat = get_category_sp(row["organism"])
        if cat is not None:
            return cat
        
    # rule 3: default bucket
    return "other"


chunk_size = 1_000_000           # adapt to your machine
needed_cols = ["acc", "organism", "librarysource"]  # read minimal set; we create the rest
acc_seen = set()
first_write = True
'''
for chunk in pd.read_csv(amr_csv, dtype=str, chunksize=chunk_size):

    # Organism_type  (Metagenome | Isolate)
    is_metagenome = (
        chunk["organism"].str.contains("metagenome", case=False, na=False) |
        chunk["librarysource"].str.contains("METAGENOMIC|METATRANSCRIPTOMIC", case=False, na=False)
    )
    chunk["organism_type"] = is_metagenome.map({True: "Metagenome", False: "Isolate"})

    # Metagenome_category
    chunk["metagenome_category"] = chunk.apply(decide_category, axis=1)

    # Append to output
    chunk.to_csv(
        output2,
        mode="w" if first_write else "a",
        index=False,
        header=first_write
    )
    first_write = False
'''
# Minimal output file
for chunk in pd.read_csv(amr_csv, dtype=str, usecols=needed_cols, chunksize=chunk_size):

    chunk = chunk.dropna(subset=["acc"])                 # discard rows w/o acc
    chunk = chunk[~chunk["acc"].isin(acc_seen)]          # keep only NEW accs
    acc_seen.update(chunk["acc"])

    # Organism_type  (Metagenome | Isolate)
    is_metagenome = (
        chunk["organism"].str.contains("metagenome", case=False, na=False) |
        chunk["librarysource"].str.contains("METAGENOMIC|METATRANSCRIPTOMIC", case=False, na=False)
    )
    chunk["organism_type"] = is_metagenome.map({True: "Metagenome", False: "Isolate"})

    # Metagenome_category
    chunk["metagenome_category"] = chunk.apply(decide_category, axis=1)

    out_chunk = chunk[["acc", "organism_type", "metagenome_category"]]

    # Append to output
    out_chunk.to_csv(
        output1,
        mode="w" if first_write else "a",
        index=False,
        header=first_write
    )
    first_write = False