import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

# Read the CSV file
df = pd.read_csv("./data/card_metadata_aro_geolocation_metagenomes_onlydateloc_noamplicon.csv", low_memory=False)

# ---- STEP 1: Count unique 'acc' per 'metagenome_category' ----
df_unique = df[['acc', 'metagenome_category']].drop_duplicates()
metagenome_counts = df_unique['metagenome_category'].value_counts()

# ---- STEP 2: Process ARO_DrugClass occurrences ----
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

# ---- STEP 3: Calculate the percentage distribution of ARO_DrugClass ----
drug_distribution = df_exploded.groupby(['metagenome_category', 'drug_class_grouped']).size().unstack(fill_value=0)

# Normalize to percentages within each metagenome_category
drug_distribution_percentage = drug_distribution.div(drug_distribution.sum(axis=1), axis=0) * 100

# ---- PLOTTING ----
fig, ax = plt.subplots(figsize=(12, 6))

categories = metagenome_counts.index
x = np.arange(len(categories))  # X positions for bars

bar_width = 0.9  # Width of each bar

target_drug_classes = list(drug_classes_to_keep) + ["colistin"]
colors = plt.cm.tab20b(np.linspace(0, 1, len(target_drug_classes)))

bottoms = np.zeros(len(categories))
for i, drug_class in enumerate(target_drug_classes):
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
ax.legend(title='ARO Drug Class', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('./data/metacategoryanddrugclass_relevantones_noothers.png', dpi=1200)
plt.savefig('./data/metacategoryanddrugclass_relevantones_noothers.svg', format='svg')
