import pandas as pd

# Define the file paths
df_total = pd.read_csv("../data/plasmids_sra_metadata.csv", dtype=str, engine='python')
output_file = "../data/plasmids_sra_metadata_extended.csv"

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

# Function to map organism to metagenome category
def get_category(organism: str) -> str | None:
    if pd.isna(organism):
        return None
    for category, values in categories.items():
        if organism in values:
            return category
    return None

# organism_type  (Metagenome | Isolate)
is_meta = df_total["organism"].str.contains("metagenome", case=False, na=False)
df_total["organism_type"] = is_meta.map({True: "Metagenome", False: "Isolate"})

# metagenome_category  (human, soil, …, other, or <NA>)
def decide_category(row):
    if row["organism_type"] == "Isolate":
        return pd.NA                # Rule 2: isolates → NA
    cat = get_category(row["organism"])
    return cat if cat is not None else "other"   # Rules 1 & 3

df_total["metagenome_category"] = df_total.apply(decide_category, axis=1)

# Output file
df_total.to_csv(output_file, index=False)