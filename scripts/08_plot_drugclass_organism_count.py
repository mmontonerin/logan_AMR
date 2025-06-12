import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the CSV file
df = pd.read_csv("./data/card_aro_entries_filtered_metagenomes.csv", low_memory=False)

# Filter out rows where assay_type is "AMPLICON"
df = df[df['assay_type'] != 'AMPLICON']

# Extract only the first ARO_DrugClass from the semicolon-separated list
df['first_aro_drug_class'] = df['ARO_DrugClass'].apply(lambda x: x.split(';')[0] if isinstance(x, str) else "")

# ---- STEP 1: Count unique 'acc' per 'metagenome_category' ----
df_unique = df[['acc', 'metagenome_category']].drop_duplicates()
metagenome_counts = df_unique['metagenome_category'].value_counts()

# ---- STEP 2: Get the top 5 most common ARO_DrugClass ----
top_15_drug_classes = df['first_aro_drug_class'].value_counts().nlargest(15).index.tolist()

# Replace all other drug classes with "Other"
df['drug_class_grouped'] = df['first_aro_drug_class'].apply(lambda x: x if x in top_15_drug_classes else "Other")

# ---- STEP 3: Calculate the percentage distribution of ARO_DrugClass ----
drug_distribution = df.groupby(['metagenome_category', 'drug_class_grouped']).size().unstack(fill_value=0)

# Normalize to percentages within each metagenome_category
drug_distribution_percentage = drug_distribution.div(drug_distribution.sum(axis=1), axis=0) * 100

# ---- PLOTTING ----
fig, ax = plt.subplots(figsize=(12, 6))

categories = metagenome_counts.index
x = np.arange(len(categories))  # X positions for bars

bar_width = 1  # Width of each bar

# Bar plot for total unique acc per metagenome_category
#ax.bar(x - bar_width/2, metagenome_counts.values, width=bar_width, label='Total Unique ACCs', color='lightblue')

# Stacked bars for ARO_DrugClass percentages
colors = plt.cm.tab20b(np.linspace(0, 1, len(top_15_drug_classes) + 1))  # +1 for "Other"
bottoms = np.zeros(len(categories))

for i, drug_class in enumerate(top_15_drug_classes + ['Other']):  # Include "Other"
    ax.bar(
        x,
        drug_distribution_percentage[drug_class].reindex(categories, fill_value=0).values * (metagenome_counts.values / 100), 
        width=bar_width, 
        bottom=bottoms, 
        label=drug_class, 
        color=colors[i] if drug_class != "Other" else "grey"  # Grey color for "Other"
    )
    bottoms += drug_distribution_percentage[drug_class].reindex(categories, fill_value=0).values * (metagenome_counts.values / 100)

# Labels and legend
ax.set_xlabel('Metagenome')
ax.set_ylabel('Count / Percentage')
ax.set_title('ARO Drug Class Distribution')
ax.set_xticks(x)
ax.set_xticklabels(categories, rotation=90)
ax.legend(title='ARO Drug Class', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('./data/metacategory_and_drugclass_counts15top.png', dpi=800)
