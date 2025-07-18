import pandas as pd

# Define the file paths
input_file = "../data/full_card_metadata_aro_allfilters.csv"
output_file = "../data/full_card_metadata_aro_allfilters_metagenomes.csv"

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
    "fish": ["fish metagenome", "fish gut metagenome"],
    "marine": ["marine metagenome", "seawater metagenome"],
    "freshwater": ["freshwater metagenome", "lake water metagenome", "groundwater metagenome"],
    "soil": ["soil metagenome"],
    "wastewater": ["wastewater metagenome"],
    "plant": ["plant metagenome", "root metagenome", "leaf metagenome"]
}

# Function to map organism to metagenome category
def get_category(organism):
    for category, values in categories.items():
        if organism in values:
            return category
    return None

df = pd.read_csv(input_file, low_memory=False)

df["metagenome_category"] = df["organism"].apply(get_category)
df_filtered = df.dropna(subset=["metagenome_category"])
df_filtered.to_csv(output_file, index=False)