import pandas as pd
import matplotlib.pyplot as plt

# Load the full organism counts data (update with the correct file path)
full_counts_file = './data/organism_counts_meta_full_with_coldate_and_loc.csv'
organism_counts_df = pd.read_csv(full_counts_file, header=0)

# Define the sets of categories
human_meta = {"human gut metagenome", "human metagenome", "human oral metagenome", "human skin metagenome", "human feces metagenome", 
              "human vaginal metagenome", "human nasopharyngeal metagenome", "human lung metagenome", "human saliva metagenome", 
              "human reproductive system metagenome", "human urinary tract metagenome", "human eye metagenome", "human blood metagenome", 
              "human bile metagenome", "human tracheal metagenome", "human brain metagenome", "human milk metagenome", "human semen metagenome", 
              "human skeleton metagenome"}
farm = {"bovine gut metagenome", "bovine metagenome", "pig gut metagenome", "pig metagenome", "chicken gut metagenome", "chicken metagenome", 
        "sheep gut metagenome", "sheep metagenome"}
fish = {"fish metagenome", "fish gut metagenome"}
sea = {"marine metagenome", "seawater metagenome"}
freshwater = {"freshwater metagenome", "lake water metagenome", "groundwater metagenome"}
soil = {"soil metagenome"}
waste = {"wastewater metagenome"}
plants = {"plant metagenome", "root metagenome", "leaf metagenome"}

# Initialize categories as a dictionary to store sums
categories = {
    "human": human_meta,
    "livestock": farm,
    "fish": fish,
    "marine": sea,
    "freshwater": freshwater,
    "soil": soil,
    "wastewater": waste,
    "plant": plants
}

# Initialize a dictionary to store counts for each category
category_counts = {category: 0 for category in categories}
other_count = 0

# Iterate through the organisms and categorize them
for index, row in organism_counts_df.iterrows():
    organism = row["organism"]
    count = row["count"]
    assigned = False

    for category, category_set in categories.items():
        if organism in category_set:
            category_counts[category] += count
            assigned = True
            break
    
    # If the organism does not belong to any category, count it as "other"
    if not assigned:
        other_count += count

# Add the "other" count
category_counts["other"] = other_count

# Prepare data for plotting
category_labels = list(category_counts.keys())
category_values = list(category_counts.values())

# Plot the bar chart
plt.figure(figsize=(10, 6))
plt.bar(category_labels, category_values, color='skyblue')
plt.xlabel('Metagenome category', fontsize=12)
plt.ylabel('Count', fontsize=12)
plt.title('Metagenomic samples with AMR detection', fontsize=14)
plt.xticks(rotation=45, ha="right")
plt.tight_layout()

# Save the plot
plt.savefig('./data/organism_meta_category_counts.png', dpi=800)

# Second plot: Without "other" category
category_counts_no_other = {k: v for k, v in category_counts.items() if k != "other"}

plt.figure(figsize=(12, 6))
plt.bar(category_counts_no_other.keys(), category_counts_no_other.values(), color="skyblue")
plt.xlabel("Metagenome category", fontsize=14)
plt.ylabel("Count", fontsize=14)
plt.title("Metagenomic samples with AMR detection", fontsize=16)
plt.xticks(rotation=45, ha="right")
plt.grid(axis="y", linestyle="--", alpha=0.7)

plt.tight_layout()
plt.savefig('./data/organism_meta_category_counts_noother.png', dpi=800)