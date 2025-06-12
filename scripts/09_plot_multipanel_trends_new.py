import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from matplotlib.ticker import MaxNLocator
from matplotlib import colormaps

# Load the dataset
df = pd.read_csv("./data/card_metadata_aro_geolocation_metagenomes_onlydateloc_noamplicon_nouncalculated.csv", low_memory=False)

# Define drug classes to keep
drug_classes_to_keep = {
    "macrolide antibiotic", "glycopeptide antibiotic", "lincosamide antibiotic",
    "tetracycline antibiotic", "fluoroquinolone antibiotic", "aminoglycoside antibiotic",
    "penicillin beta-lactam", "glycylcycline", "cephalosporin"
}

# Expand drug classes into separate rows
df_exploded = df.assign(ARO_DrugClass=df['ARO_DrugClass'].dropna().str.split(';')).explode('ARO_DrugClass')
df_exploded = df_exploded[df_exploded['ARO_DrugClass'].isin(drug_classes_to_keep)]

# Parse collection date
df_exploded['collection_date_sam'] = pd.to_datetime(
    df_exploded['collection_date_sam'].str.strip('[]'),
    format='%Y-%m-%d',
    errors='coerce'
)
df_exploded.dropna(subset=['collection_date_sam', 'metagenome_category'], inplace=True)

# Deduplicate acc per ARO_DrugClass to ensure each accession is only counted once
df_unique = df_exploded.drop_duplicates(subset=['acc', 'ARO_DrugClass'])

# Total unique accessions per drug class
total_unique_samples = df_unique.groupby('ARO_DrugClass')['acc'].nunique()

# Prepare 5-year bin and normalize
df_unique['year'] = df_unique['collection_date_sam'].dt.year
df_unique['year_bin'] = (df_unique['year'] // 5) * 5

aggregated_data = []

continents = sorted(df_unique['geo_loc_name_country_continent_calc'].dropna().unique())
continents.insert(0, 'Worldwide')  # Worldwide first

for drug in sorted(drug_classes_to_keep):
    for continent in continents:
        if continent == 'Worldwide':
            subset = df_unique[df_unique['ARO_DrugClass'] == drug]
        else:
            subset = df_unique[
                (df_unique['ARO_DrugClass'] == drug) &
                (df_unique['geo_loc_name_country_continent_calc'] == continent)
            ]

        for meta_cat in subset['metagenome_category'].unique():
            cat_data = subset[subset['metagenome_category'] == meta_cat]
            grouped = cat_data.groupby('year_bin')['acc'].nunique().reset_index()
            grouped['total'] = total_unique_samples[drug]
            grouped = grouped[grouped['acc'] >= 3]  # Require min 3 samples per bin
            if not grouped.empty:
                grouped['normalized_prevalence'] = grouped['acc'] / grouped['total']
                grouped['log_normalized_prevalence'] = np.log1p(grouped['normalized_prevalence'])
                grouped['ARO_DrugClass'] = drug
                grouped['continent'] = continent
                grouped['metagenome_category'] = meta_cat
                grouped.rename(columns={'year_bin': 'collection_year'}, inplace=True)
                aggregated_data.append(grouped)

# Combine all data
final_df = pd.concat(aggregated_data)

# LOESS smoothing function
def loess_smoothing(x, y, frac=0.5):
    if len(x) < 3:
        return x, y
    loess = sm.nonparametric.lowess(y, x, frac=frac, return_sorted=True)
    return loess[:, 0], loess[:, 1]

# Color map for metagenome categories
categories = sorted(final_df['metagenome_category'].unique())
cmap = colormaps.get_cmap('Set3')
color_dict = {cat: cmap(i / len(categories)) for i, cat in enumerate(categories)}

# Plotting setup
sns.set_theme(style="whitegrid")
fig, axes = plt.subplots(len(drug_classes_to_keep), len(continents),
                         figsize=(4 * len(continents), 3 * len(drug_classes_to_keep)),
                         sharex=True, sharey=True)

if len(drug_classes_to_keep) == 1:
    axes = [axes]

for i, drug in enumerate(sorted(drug_classes_to_keep)):
    for j, continent in enumerate(continents):
        ax = axes[i][j] if len(drug_classes_to_keep) > 1 else axes[j]
        subset = final_df[(final_df['ARO_DrugClass'] == drug) & (final_df['continent'] == continent)]
        if subset.empty:
            ax.set_visible(False)
            continue

        for meta_cat in subset['metagenome_category'].unique():
            cat_data = subset[subset['metagenome_category'] == meta_cat]
            x_vals, y_vals = loess_smoothing(cat_data['collection_year'].values, cat_data['log_normalized_prevalence'].values, frac=0.6)
            ax.plot(x_vals, y_vals, label=meta_cat, color=color_dict[meta_cat], linewidth=2)
            ax.scatter(cat_data['collection_year'], cat_data['log_normalized_prevalence'], color=color_dict[meta_cat], alpha=0.5)

        if i == 0:
            ax.set_title(continent)
        if j == 0:
            ax.set_ylabel(drug)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

# Add a single shared legend
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, title="Metagenome category", bbox_to_anchor=(1.02, 0.5), loc="center left")

# Adjust layout to accommodate legend
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.savefig("./data/multipanel_prevalence_trends_5years.png", dpi=1200)
plt.savefig("./data/multipanel_prevalence_trends_5years.svg")

