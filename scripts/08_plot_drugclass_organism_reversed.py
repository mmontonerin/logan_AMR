import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

# Read the CSV file
df = pd.read_csv("./data/card_metadata_aro_geolocation_metagenomes_onlydateloc_noamplicon.csv", low_memory=False)

# ---- STEP 1: Process ARO_DrugClass occurrences ----
drug_classes_to_keep = {
    "macrolide antibiotic", "glycopeptide antibiotic", "lincosamide antibiotic", 
    "tetracycline antibiotic", "fluoroquinolone antibiotic", "aminoglycoside antibiotic", 
    "penicillin beta-lactam", "glycylcycline", "cephalosporin"
}

# Expand the drug classes into separate rows
df_exploded = df.assign(ARO_DrugClass=df['ARO_DrugClass'].dropna().str.split(';')).explode('ARO_DrugClass')

# Ensure "colistin" detection in AMR_GeneFamily using regex
def classify_drug(row):
    if pd.notna(row['AMR_GeneFamily']) and re.search(r'\bcolistin\b', row['AMR_GeneFamily'], re.IGNORECASE):
        return "colistin"
    return row['ARO_DrugClass'] if row['ARO_DrugClass'] in drug_classes_to_keep else "Others"

df_exploded['drug_class_grouped'] = df_exploded.apply(classify_drug, axis=1)

# ---- STEP 2: Count occurrences per metagenome_category ----
drug_counts = df_exploded.groupby(['drug_class_grouped', 'metagenome_category']).size().unstack(fill_value=0)

# Define plotting function
def plot_drug_class_distribution(drug_counts, filename, exclude_others=False):
    fig, ax = plt.subplots(figsize=(12, 6))

    # Filter out "Others" if requested
    if exclude_others:
        drug_counts = drug_counts.drop(index="Others", errors='ignore')

    drug_classes = drug_counts.index  # X positions (drug classes)
    x = np.arange(len(drug_classes))  # X positions for bars
    bar_width = 0.15  # Adjust bar width for better visibility

    metagenome_categories = drug_counts.columns  # Different metagenome categories
    colors = plt.colormaps['tab10'](np.linspace(0, 1, len(metagenome_categories)))  # Fixed color mapping

    # Plot bars for each metagenome category
    for i, metagenome in enumerate(metagenome_categories):
        ax.bar(x + i * bar_width, drug_counts[metagenome], width=bar_width, label=metagenome, color=colors[i])

    # Labels and legend
    ax.set_xlabel('ARO Drug Class')
    ax.set_ylabel('Count of Instances')
    title_suffix = " (Excluding 'Others')" if exclude_others else ""
    ax.set_title(f'Distribution of ARO Drug Classes Across Metagenome Categories{title_suffix}')
    ax.set_xticks(x + (len(metagenome_categories) - 1) * bar_width / 2)
    ax.set_xticklabels(drug_classes, rotation=45, ha='right')
    ax.legend(title='Metagenome Category', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()

    # Save the plot in PNG and SVG formats
    plt.savefig(f'./data/{filename}.png', dpi=1200)
    plt.savefig(f'./data/{filename}.svg', format='svg')

# ---- PLOT 1: All drug classes (including "Others") ----
plot_drug_class_distribution(drug_counts, "drugclass_and_metacategory_counts")

# ---- PLOT 2: Excluding "Others" ----
plot_drug_class_distribution(drug_counts, "drugclass_and_metacategory_counts_filtered", exclude_others=True)
