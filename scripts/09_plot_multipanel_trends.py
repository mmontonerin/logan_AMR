import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from matplotlib.ticker import MaxNLocator
from matplotlib import colormaps

# Load and prepare data
df = pd.read_csv("./data/card_metadata_aro_geolocation_metagenomes_onlydateloc_noamplicon_nouncalculated.csv", low_memory=False)

# Define drug classes to keep
drug_classes_to_keep = {
    "macrolide antibiotic", "glycopeptide antibiotic", "lincosamide antibiotic", 
    "tetracycline antibiotic", "fluoroquinolone antibiotic", "aminoglycoside antibiotic", 
    "penicillin beta-lactam", "glycylcycline", "cephalosporin"
}

# Expand and filter
df_exploded = df.assign(ARO_DrugClass=df['ARO_DrugClass'].dropna().str.split(';')).explode('ARO_DrugClass')
df_exploded = df_exploded[df_exploded['ARO_DrugClass'].isin(drug_classes_to_keep)]

# Parse dates
df_exploded['collection_date_sam'] = pd.to_datetime(df_exploded['collection_date_sam'].str.strip('[]'), errors='coerce')
df_exploded.dropna(subset=['collection_date_sam', 'metagenome_category'], inplace=True)

# Compute raw prevalence per group
prevalence_df = df_exploded.groupby(['ARO_DrugClass', 'collection_date_sam', 'geo_loc_name_country_continent_calc', 'metagenome_category'])['acc'].nunique().reset_index(name='raw_prevalence')

# Normalize to 0–1 scale within each drug class
prevalence_df['normalized_prevalence'] = prevalence_df.groupby('ARO_DrugClass')['raw_prevalence'].transform(lambda x: (x - x.min()) / (x.max() - x.min()))

# Aggregate data
aggregated_data = []
for (drug, continent, meta_cat), group in prevalence_df.groupby(['ARO_DrugClass', 'geo_loc_name_country_continent_calc', 'metagenome_category']):
    yearly = group.groupby(group['collection_date_sam'].dt.year)['normalized_prevalence'].agg(['mean', 'count']).reset_index()
    yearly['interval'] = 'yearly'
    if (yearly['count'] >= 3).sum() < 1:
        group['year'] = group['collection_date_sam'].dt.year
        group['year_bin'] = (group['year'] // 5) * 5
        binned = group.groupby('year_bin')['normalized_prevalence'].agg(['mean', 'count']).reset_index()
        binned = binned[binned['count'] >= 3]
        if not binned.empty:
            binned = binned.rename(columns={'year_bin': 'collection_year'})
            binned['interval'] = '5year'
            binned['ARO_DrugClass'] = drug
            binned['continent'] = continent
            binned['metagenome_category'] = meta_cat
            aggregated_data.append(binned)
    else:
        yearly = yearly[yearly['count'] >= 3]
        yearly = yearly.rename(columns={'collection_date_sam': 'collection_year'})
        yearly['collection_year'] = yearly['collection_year'].astype(int)
        yearly['ARO_DrugClass'] = drug
        yearly['continent'] = continent
        yearly['metagenome_category'] = meta_cat
        aggregated_data.append(yearly)

# Combine and smooth
final_df = pd.concat(aggregated_data)

def loess_smoothing(x, y, frac=0.6):  # slightly smoother
    if len(x) < 3:
        return x, y
    loess = sm.nonparametric.lowess(y, x, frac=frac, return_sorted=True)
    return loess[:, 0], loess[:, 1]

# Color settings
categories = sorted(final_df['metagenome_category'].unique())
cmap = colormaps.get_cmap('Set3')
color_dict = {cat: cmap(i / len(categories)) for i, cat in enumerate(categories)}

# Plotting setup
sns.set_theme(style="whitegrid")
continents = ['Worldwide'] + sorted(df_exploded['geo_loc_name_country_continent_calc'].dropna().unique())

fig, axes = plt.subplots(len(drug_classes_to_keep), len(continents), figsize=(4*len(continents), 3*len(drug_classes_to_keep)), sharex=True, sharey=True)

if len(drug_classes_to_keep) == 1:
    axes = [axes]

# Plot panels
for i, drug in enumerate(sorted(drug_classes_to_keep)):
    for j, continent in enumerate(continents):
        ax = axes[i][j] if len(drug_classes_to_keep) > 1 else axes[j]
        if continent == 'Worldwide':
            subset = final_df[final_df['ARO_DrugClass'] == drug]
        else:
            subset = final_df[(final_df['ARO_DrugClass'] == drug) & (final_df['continent'] == continent)]

        if subset.empty:
            ax.set_visible(False)
            continue

        for meta_cat in subset['metagenome_category'].unique():
            cat_data = subset[subset['metagenome_category'] == meta_cat]
            x_vals, y_vals = loess_smoothing(cat_data['collection_year'].values, cat_data['mean'].values)
            ax.plot(x_vals, y_vals, color=color_dict[meta_cat], linewidth=2)
            ax.scatter(cat_data['collection_year'], cat_data['mean'], color=color_dict[meta_cat], alpha=0.5)

        if i == 0:
            ax.set_title(continent)
        if j == 0:
            ax.set_ylabel(drug)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

# One shared legend
handles = [plt.Line2D([0], [0], color=color_dict[cat], lw=3) for cat in categories]
fig.legend(handles, categories, title="Metagenome category", bbox_to_anchor=(1.02, 0.5), loc="center left")

plt.tight_layout(rect=[0, 0, 0.85, 1])

plt.savefig("./data/multipanel_prevalence_trends_3.png", dpi=1200)
plt.savefig("./data/multipanel_prevalence_trends_3.svg")