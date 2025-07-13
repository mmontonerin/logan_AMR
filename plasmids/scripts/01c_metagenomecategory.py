import pandas as pd

# Define the file paths
input_file = "../data/amr_metadata_full.csv"
output_file = "../data/amr_metadata_metagenomes.csv"

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
def get_category(organism):
    for category, values in categories.items():
        if organism in values:
            return category
    return "Other"

df = pd.read_csv(input_file, low_memory=False)

mask     = df["organism"].str.contains("metagenome", case=False, na=False)
df_meta  = df.loc[mask].copy()

df_meta["metagenome_category"] = df_meta["organism"].apply(get_category)
df_meta.to_csv(output_file, index=False)