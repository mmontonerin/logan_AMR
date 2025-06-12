import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import make_interp_spline
from matplotlib.ticker import MaxNLocator
from matplotlib import colormaps

# -----------------------------
# Load dataset
# -----------------------------
df = pd.read_csv(
    "./data/card_metadata_aro_geolocation_metagenomes_onlydateloc_noamplicon_nouncalculated.csv",
    low_memory=False
)

# -----------------------------
# Define relevant drug classes
# -----------------------------
drug_classes_to_keep = {
    "macrolide antibiotic", "glycopeptide antibiotic", "lincosamide antibiotic",
    "tetracycline antibiotic", "fluoroquinolone antibiotic", "aminoglycoside antibiotic",
    "penicillin beta-lactam", "glycylcycline", "cephalosporin"
}

# -----------------------------
# Preprocessing
# -----------------------------
df_exploded = df.assign(ARO_DrugClass=df['ARO_DrugClass'].dropna().str.split(';')).explode('ARO_DrugClass')
df_exploded = df_exploded[df_exploded['ARO_DrugClass'].isin(drug_classes_to_keep)]

df_exploded['collection_date_sam'] = pd.to_datetime(df_exploded['collection_date_sam'].str.strip('[]'), errors='coerce')
df_exploded.dropna(subset=['collection_date_sam', 'metagenome_category'], inplace=True)

df_exploded['year'] = df_exploded['collection_date_sam'].dt.year
df_exploded['year_bin'] = (df_exploded['year'] // 5) * 5

# -----------------------------
# Aggregate data with filtering
# -----------------------------
aggregated_data = []

for (drug, continent, meta_cat), group in df_exploded.groupby(['ARO_DrugClass', 'geo_loc_name_country_continent_calc', 'metagenome_category']):
    binned = group.groupby('year_bin')['acc'].nunique().reset_index()
    binned = binned.rename(columns={'acc': 'positive_samples'})
    binned['total_samples'] = group.groupby('year_bin')['acc'].count().values
    binned['prevalence'] = binned['positive_samples'] / binned['total_samples']
    binned = binned[binned['positive_samples'] >= 3]

    if len(binned) >= 3:
        binned['ARO_DrugClass'] = drug
        binned['continent'] = continent
        binned['metagenome_category'] = meta_cat
        aggregated_data.append(binned)

final_df = pd.concat(aggregated_data)

# -----------------------------
# Spline smoothing
# -----------------------------
def spline_smoothing(x, y, num=200):
    x_sorted, y_sorted = zip(*sorted(zip(x, y)))
    n_points = len(x_sorted)

    if n_points < 3:
        return x_sorted, y_sorted  # not enough for spline
    elif n_points == 3:
        k = 2  # quadratic spline
    else:
        k = 3  # cubic spline

    try:
        spline = make_interp_spline(x_sorted, y_sorted, k=k)
        x_new = np.linspace(min(x_sorted), max(x_sorted), num=num)
        y_new = spline(x_new)
        return x_new, y_new
    except Exception:
        return x_sorted, y_sorted  # fallback if spline still fails


# -----------------------------
# Set up plot
# -----------------------------
categories = sorted(final_df['metagenome_category'].unique())
cmap = colormaps.get_cmap('Set3')
color_dict = {cat: cmap(i / len(categories)) for i, cat in enumerate(categories)}

sns.set_theme(style="whitegrid")
continents = ['Worldwide'] + sorted(df_exploded['geo_loc_name_country_continent_calc'].dropna().unique())
fig, axes = plt.subplots(
    len(drug_classes_to_keep), len(continents),
    figsize=(4 * len(continents), 3 * len(drug_classes_to_keep)),
    sharex=True, sharey=True
)

# Handle axis shape in case of single row/column
if len(drug_classes_to_keep) == 1:
    axes = [axes]
if len(continents) == 1:
    axes = [[ax] for ax in axes]

# -----------------------------
# Plot data
# -----------------------------
for i, drug in enumerate(sorted(drug_classes_to_keep)):
    for j, continent in enumerate(continents):
        ax = axes[i][j]
        if continent == 'Worldwide':
            subset = final_df[final_df['ARO_DrugClass'] == drug]
        else:
            subset = final_df[(final_df['ARO_DrugClass'] == drug) & (final_df['continent'] == continent)]

        if subset.empty:
            ax.set_visible(False)
            continue

        for meta_cat in subset['metagenome_category'].unique():
            cat_data = subset[subset['metagenome_category'] == meta_cat]
            binned = cat_data.groupby('year_bin')['prevalence'].mean().reset_index()
            x_vals, y_vals = spline_smoothing(binned['year_bin'].values, binned['prevalence'].values)   
            ax.plot(x_vals, y_vals, label=meta_cat, color=color_dict[meta_cat], linewidth=2)
            ax.scatter(cat_data['year_bin'], cat_data['prevalence'], color=color_dict[meta_cat], alpha=0.5)

        if i == 0:
            ax.set_title(continent)
        if j == 0:
            ax.set_ylabel(drug)

        ax.set_ylim(0, 1)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

# -----------------------------
# Save main plot (without legend)
# -----------------------------
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.savefig("./data/multipanel_prevalence_trends_5years_spline.png", dpi=1200)
plt.savefig("./data/multipanel_prevalence_trends_5years_spline.svg")
plt.close(fig)

# -----------------------------
# Save legend separately
# -----------------------------
# Get handles and labels from one axis
for row in axes:
    for ax in row:
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            break
    if handles:
        break

legend_fig, legend_ax = plt.subplots(figsize=(4, len(categories) * 0.3))
legend_ax.axis('off')
legend = legend_ax.legend(
    handles,
    labels,
    title="Metagenome category",
    loc='center',
    frameon=False,
    ncol=1,
    handlelength=2,
    labelspacing=1.2
)

legend_fig.tight_layout()
legend_fig.savefig("./data/metagenome_category_legend.png", dpi=300, bbox_inches='tight')
legend_fig.savefig("./data/metagenome_category_legend.svg", bbox_inches='tight')
plt.close(legend_fig)
