import pandas as pd

# Define the file paths
input_file = "./data/entries_with_dates_continent.csv"
output_file = "./data/card_aro_entries_filtered_metagenomes.csv"

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

# Process CSV in chunks
chunksize = 1000000
chunks = []
chunk_num = 1

for chunk in pd.read_csv(input_file, chunksize=chunksize, low_memory=False):
    print(f"Processing chunk {chunk_num}...")
    chunk["metagenome_category"] = chunk["organism"].apply(get_category)
    chunk_filtered = chunk.dropna(subset=["metagenome_category"])
    chunks.append(chunk_filtered)
    chunk_num += 1

# Concatenate all filtered chunks
df_filtered = pd.concat(chunks, ignore_index=True)

# Save the filtered data
df_filtered.to_csv(output_file, index=False)

print(f"Filtered data saved to {output_file}")