import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the CSV file
df = pd.read_csv("./data/card_aro_entries_metagenomes_and_drugclass_noamplicon.csv", low_memory=False)

# ---- STEP 1: Count unique 'acc' per 'metagenome_category' ----
df_unique = df[['acc', 'metagenome_category']].drop_duplicates()
metagenome_counts = df_unique['metagenome_category'].value_counts()

# Replace all other drug classes with "Other"
drug_classes = df['drug_action'].dropna().unique()

# ---- STEP 3: Calculate the percentage distribution of ARO_DrugClass ----
drug_distribution = df.groupby(['metagenome_category', 'drug_action']).size().unstack(fill_value=0)

# Normalize to percentages within each metagenome_category
drug_distribution_percentage = drug_distribution.div(drug_distribution.sum(axis=1), axis=0) * 100

# ---- PLOTTING ----
fig, ax = plt.subplots(figsize=(12, 6))

categories = metagenome_counts.index
x = np.arange(len(categories))  # X positions for bars

bar_width = 0.8  # Width of each bar

# Bar plot for total unique acc per metagenome_category
#ax.bar(x - bar_width/2, metagenome_counts.values, width=bar_width, label='Total Unique ACCs', color='lightblue')

# Stacked bars for ARO_DrugClass percentages
colors = plt.cm.CMRmap(np.linspace(0.3, 0.9, len(drug_classes)))
bottoms = np.zeros(len(categories))

for i, drug_class in enumerate(drug_classes): 
    # Check if the drug_class exists in the columns of drug_distribution_percentage
    if drug_class in drug_distribution_percentage.columns:
        ax.bar(
            x,
            drug_distribution_percentage[drug_class].reindex(categories, fill_value=0).values * (metagenome_counts.values / 100), 
            width=bar_width, 
            bottom=bottoms, 
            label=drug_class, 
            color=colors[i]
        )
        bottoms += drug_distribution_percentage[drug_class].reindex(categories, fill_value=0).values * (metagenome_counts.values / 100)

# Labels and legend
ax.set_xlabel('Metagenome')
ax.set_ylabel('Count / Percentage')
ax.set_title('ARO Drug Class Distribution')
ax.set_xticks(x)
ax.set_xticklabels(categories, rotation=90)
ax.legend(title='ARO Drug Class Mechanism of Action', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('./data/metacategory_and_drugclass_action.png', dpi=800)
