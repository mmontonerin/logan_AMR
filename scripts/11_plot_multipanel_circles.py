import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess

# --------------------------------------------------------------------
# 0 PRE-PROCESSING DATA
# --------------------------------------------------------------------

# Input table
df = pd.read_csv(
    "./data/card_metadata_aro_geolocation_metagenomes_onlydateloc_noamplicon_nouncalculated.csv",
    low_memory=False
)

# Only keep most relevant drug classes (7)
drug_classes_to_keep = {
    "glycopeptide antibiotic", 
    "tetracycline antibiotic", "fluoroquinolone antibiotic",
    "aminoglycoside antibiotic", "penicillin beta-lactam",
    "glycylcycline", "cephalosporin"
}

# Get individual instances for each drug class 
# (given that most AMR notations will specify different drug resistance for one same marker)
df_exploded = (
    df.assign(ARO_DrugClass=df['ARO_DrugClass'].dropna().str.split(';'))
      .explode('ARO_DrugClass')
      .query("ARO_DrugClass in @drug_classes_to_keep")
)
df_exploded['collection_date_sam'] = pd.to_datetime(
    df_exploded['collection_date_sam'].str.strip('[]'),
    errors='coerce'
)
df_exploded = df_exploded.dropna(subset=['collection_date_sam',
                                         'metagenome_category'])

# Generate yearly bins 
df_exploded['year_bin'] = df_exploded['collection_date_sam'].dt.year


# --------------------------------------------------------------------
# 1 CREATE COUNTS AND SET NORMALIZATION
# --------------------------------------------------------------------

SAMPLE_ID = 'acc'

## Count total number of samples per continent / year / metagenome
# I keep one row per sample
total_samples_df = (
    df_exploded
      .drop_duplicates(subset=[SAMPLE_ID])
      .groupby(['geo_loc_name_country_continent_calc',
                'metagenome_category',
                'year_bin'])
      .agg(total_samples=(SAMPLE_ID, 'count'))
      .reset_index()
)

## Count positive samples for each drug class
# If a specific accession has one hit or more -> included
# We count a specific time in 
positive_samples_df = (
    df_exploded
      .drop_duplicates(subset=[SAMPLE_ID, 'ARO_DrugClass']) # ≥1 hit -> 1 flag
      .groupby(['ARO_DrugClass',
                'geo_loc_name_country_continent_calc',
                'metagenome_category',
                'year_bin'])
      .agg(positive_samples=(SAMPLE_ID, 'count'))
      .reset_index()
)
# Before I had .query("positive_samples >= 3")
# It might be too strict on continents with very few data, and too lax for the big sampled ones

## Merge counts and calculate prevalence
dot_data = (
    positive_samples_df
      .merge(total_samples_df,
             on=['geo_loc_name_country_continent_calc',
                 'metagenome_category',
                 'year_bin'],
             how='left', validate='many_to_one')
      .assign(prevalence=lambda d: d.positive_samples / d.total_samples)
      .dropna(subset=['total_samples'])
)

# Posible things I can try to add only dots with higher statistical significance:
# .query("total_samples >= 20 and positive_samples >= 3") -> in dot_data
# contains at least 20 samples, and at least 3 positives

# Added to section 2: 
# I chose to not filter bins and increase size range of bubbles so that the small data ones will look very very small
# Moreover, the LOWESS line has reduced opacity and one can visualize if there are big areas with no data points.

# --------------------------------------------------------------------
# 2 FACETTED BUBBLE PLOT
# --------------------------------------------------------------------

sns.set_theme(style="ticks")
palette = sns.color_palette("colorblind",
                            n_colors=dot_data['metagenome_category']
                                     .nunique())
color_dict = dict(zip(sorted(dot_data['metagenome_category'].unique()),
                      palette))

g = sns.relplot(
        data     = dot_data,
        x        = 'year_bin',
        y        = 'prevalence',
        hue      = 'metagenome_category',
        size     = 'total_samples',
        sizes    = (10, 800),
        alpha    = .5,
        kind     = 'scatter',
        row      = 'ARO_DrugClass',
        col      = 'geo_loc_name_country_continent_calc',
        facet_kws= dict(sharey=False),
        height   = 3.2,
        aspect   = 1.1,
        palette  = color_dict,
        legend   = 'brief'
)

# --------------------------------------------------------------------
# 3 OVERLAY LOWESS SMOOTH ON EACH drug/continent PAIR
# --------------------------------------------------------------------
# iterate over the facets the way seaborn stores them
for (row_label, col_label), ax in g.axes_dict.items():
    drug  = row_label
    cont  = col_label
    sub   = dot_data[(dot_data['ARO_DrugClass']==drug) &
                     (dot_data['geo_loc_name_country_continent_calc']==cont)]
    if sub.empty:
        ax.set_axis_off() # keep grid tidy if facet is empty
        continue

    # fit + plot one LOWESS curve per metagenome category
    for cat in sub['metagenome_category'].unique():
        cat_data = sub[sub['metagenome_category'] == cat]\
                   .sort_values('year_bin')

        if len(cat_data) < 3: # LOWESS needs ≥3 points
            continue

        z = lowess(
            endog = cat_data['prevalence'],
            exog = cat_data['year_bin'],
            frac = 0.55, # smoothing window (0–1)
            return_sorted = True
        )
        ax.plot(z[:, 0], z[:, 1],
                color     = color_dict[cat],
                linewidth = 1.0,
                zorder    = 1) # behind scatter, above grid

# --------------------------------------------------------------------
# 4 GENERATE PLOT
# --------------------------------------------------------------------
g.set_axis_labels("Year", "Prevalence")
g.set_titles(row_template="{row_name}",
             col_template="{col_name}")
for lh in g.legend.legend_handles:
    lh.set_alpha(1)
    lh.set_color("black")
g.legend.set_bbox_to_anchor((1.02, 0.5))
plt.tight_layout()

g.savefig("./data/plot_dot_size_lowess_faceted.png", dpi=1200)
g.savefig("./data/plot_dot_size_lowess_faceted.svg")