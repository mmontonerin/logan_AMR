import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colormaps

# Read CSV
df = pd.read_csv("./data/card_metadata_aro_geolocation_metagenomes_onlydateloc_noamplicon.csv", low_memory=False)

# Define relevant drug classes
drug_classes_to_keep = {
    "glycopeptide antibiotic", "tetracycline antibiotic", "fluoroquinolone antibiotic",
    "aminoglycoside antibiotic", "penicillin beta-lactam", "glycylcycline", "cephalosporin"
}

# Define classify_drug() function
def classify_drug(row):
    drug = row['ARO_DrugClass']
    if drug in drug_classes_to_keep:
        return drug
    return "Others"

# Expand and filter drug classes
df_exploded = df.assign(ARO_DrugClass=df['ARO_DrugClass'].dropna().str.split(';')).explode('ARO_DrugClass')
df_exploded = df_exploded[df_exploded['ARO_DrugClass'].isin(drug_classes_to_keep)]
df_exploded['drug_class_grouped'] = df_exploded.apply(classify_drug, axis=1)

# Use Set2 colormap for consistent coloring
metagenome_categories = sorted(df_exploded['metagenome_category'].dropna().unique())
cmap = colormaps.get_cmap('Set2')
color_dict = {cat: cmap(i / len(metagenome_categories)) for i, cat in enumerate(metagenome_categories)}

# Count unique accessions
drug_counts = df_exploded.groupby(['drug_class_grouped', 'metagenome_category'])['acc'].nunique().unstack(fill_value=0)

# Plot function
def plot_drug_class_distribution(drug_counts, filename, exclude_others=False):
    fig, ax = plt.subplots(figsize=(12, 6))

    if exclude_others:
        drug_counts = drug_counts.drop(index="Others", errors='ignore')

    drug_classes = drug_counts.index
    x = np.arange(len(drug_classes))
    bar_width = 0.15

    metagenome_categories = drug_counts.columns
    colors = [color_dict[cat] for cat in metagenome_categories]

    for i, metagenome in enumerate(metagenome_categories):
        ax.bar(x + i * bar_width, drug_counts[metagenome], width=bar_width, label=metagenome, color=colors[i])

    ax.set_xlabel('ARO Drug Class')
    ax.set_ylabel('Unique Sample Count')
    title_suffix = " (Excluding 'Others')" if exclude_others else ""
    ax.set_title(f'Distribution of ARO Drug Classes Across Metagenome Categories{title_suffix}')
    ax.set_xticks(x + (len(metagenome_categories) - 1) * bar_width / 2)
    ax.set_xticklabels(drug_classes, rotation=45, ha='right')
    ax.legend(title='Metagenome Category', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig(f'./data/{filename}.png', dpi=1200)
    plt.savefig(f'./data/{filename}.svg', format='svg')

# PLOT 1: All
plot_drug_class_distribution(drug_counts, "side_bar_plot_drugclass_and_metacategory_counts")

# PLOT 2: Exclude "Others"
plot_drug_class_distribution(drug_counts, "side_bar_plot_drugclass_and_metacategory_counts_filtered", exclude_others=True)
