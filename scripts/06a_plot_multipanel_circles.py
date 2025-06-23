import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess

# --------------------------------------------------------------------
# 0 PRE-PROCESSING DATA
# --------------------------------------------------------------------

# Input table
amr_df = pd.read_csv(
    "../data/full_card_metadata_aro_allfilters_metagenomes2.csv",
    low_memory=False)

sra_df = pd.read_csv("../data/SRA_metadata_allfilters.csv", 
    low_memory=False)

# Only keep most relevant drug classes (9)
drug_classes_to_keep = {
    "macrolide antibiotic", "glycopeptide antibiotic",
    "lincosamide antibiotic",
    "tetracycline antibiotic", "fluoroquinolone antibiotic",
    "aminoglycoside antibiotic", "penicillin beta-lactam",
    "glycylcycline", "cephalosporin"
}

# Function to prepare tables to have date / location / metagenome category
def date_collection_metagenome(df):
    df = df.copy()
    df['collection_date_sam'] = pd.to_datetime(
        df['collection_date_sam'].str.strip('[]'),
        errors='coerce')
    df = df.dropna(subset=['collection_date_sam',
                           'metagenome_category'])
    df['year_bin'] = df['collection_date_sam'].dt.year # Generate yearly bins
    return df

amr_df = date_collection_metagenome(amr_df)
sra_df = date_collection_metagenome(sra_df)

# Get individual instances for each drug class 
# (given that most AMR notations will specify different drug resistance for one same marker)
amr_df = (
    amr_df.assign(ARO_DrugClass=amr_df['ARO_DrugClass'].dropna().str.split(';'))
      .explode('ARO_DrugClass')
      .query("ARO_DrugClass in @drug_classes_to_keep")
)

# --------------------------------------------------------------------
# 1 CREATE COUNTS AND SET NORMALIZATION
# --------------------------------------------------------------------

SAMPLE_ID = 'acc'

## Count total number of samples per continent / year / metagenome
total_samples_df = (
    sra_df
      .drop_duplicates(subset=[SAMPLE_ID])
      .groupby(['geo_loc_name_country_continent_calc',
                'metagenome_category',
                'year_bin'], as_index=False)
      .agg(total_samples=(SAMPLE_ID, 'count'))
)

## Count positive samples for each drug class
# If a specific accession has one hit or more -> included
# We count a specific time in 
positive_samples_df = (
    amr_df
      .drop_duplicates(subset=[SAMPLE_ID, 'ARO_DrugClass']) # ≥1 hit -> 1 flag
      .groupby(['ARO_DrugClass',
                'geo_loc_name_country_continent_calc',
                'metagenome_category',
                'year_bin'], as_index=False)
      .agg(positive_samples=(SAMPLE_ID, 'count'))
)
# Before I had .query("positive_samples >= 3")
# It might be too strict on continents with very few data, and too lax for the big sampled ones

## Merge counts and calculate prevalence
#Continent data
continent_data = (
    positive_samples_df
      .merge(total_samples_df,
             on=['geo_loc_name_country_continent_calc',
                 'metagenome_category',
                 'year_bin'],
             how='left', validate='many_to_one')
      .assign(prevalence=lambda d: d.positive_samples / d.total_samples)
)

# Worldwide data
world_total = (
    total_samples_df
      .groupby(['metagenome_category', 'year_bin'], as_index=False)
      .agg(total_samples=('total_samples', 'sum'))
      .assign(geo_loc_name_country_continent_calc='Worldwide')
)

world_pos = (
    positive_samples_df
      .groupby(['ARO_DrugClass', 'metagenome_category', 'year_bin'],
               as_index=False)
      .agg(positive_samples=('positive_samples', 'sum'))
      .assign(geo_loc_name_country_continent_calc='Worldwide')
)

world_data = (world_pos
    .merge(world_total,
           on=['geo_loc_name_country_continent_calc',
               'metagenome_category', 'year_bin'],
           how='left')
    .assign(prevalence=lambda d: d.positive_samples / d.total_samples)
)

# Combine both continent and worldwide data
dot_data = pd.concat([continent_data, world_data], ignore_index=True)

# Posible things I can try to add only dots with higher statistical significance:
# .query("total_samples >= 20 and positive_samples >= 3") -> in dot_data
# contains at least 20 samples, and at least 3 positives

# Added to section 2: 
# I chose to not filter bins and increase size range of bubbles so that the small data ones will look very very small
# Moreover, the LOWESS line has reduced opacity and one can visualize if there are big areas with no data points.

# --------------------------------------------------------------------
# 2 FACETTED BUBBLE PLOT
# --------------------------------------------------------------------

# Specify order of columns
col_order = ["Worldwide",
             "Africa", "Asia", "Europe",
             "North America", "Oceania", "South America"]

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
        sizes    = (20, 600),
        alpha    = .45, # translucent / overlaps visible
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
    # one LOWESS curve per metagenome category
    for cat, grp in sub.groupby('metagenome_category'):
        if len(grp) < 5: # LOWESS needs ≥ 3 points (changed to 5 for better smoothing)
            continue
        pts = grp.sort_values('year_bin')
        z   = lowess(
                endog = pts['prevalence'],
                exog  = pts['year_bin'],
                frac  = 0.55, # smoothing window [0-1]
                return_sorted=True
        )
        ax.plot(z[:, 0], z[:, 1],
                color     = color_dict[cat],
                linewidth = 5,
                alpha     = .85, # translucent / overlaps visible
                zorder    = 1) # behind the scatter bubbles

# --------------------------------------------------------------------
# 4 GENERATE PLOT
# --------------------------------------------------------------------
for ax in g.axes.flatten():
    ax.set_ylim(0, 1.05) # prevalence always in [0,1]

# Generic axis labels
g.set_axis_labels("Year", "Prevalence")

# remove seaborn’s per-panel titles
g.set_titles("")

# continent titles only on the top row
for ax, cont in zip(g.axes[0], g.col_names):
    ax.set_title(cont, fontweight='bold', pad=20)

# drug-class labels only on the first column
for ax, drug in zip(g.axes[:, 0], g.row_names or
                    sorted(dot_data['ARO_DrugClass'].unique())):
    ax.set_ylabel(drug, rotation=90, labelpad=20,
                  fontweight='bold', va='center')

# remove redundant y-labels from the inner columns
for ax in g.axes.flatten():
    if ax not in g.axes[:, 0]:
        ax.set_ylabel("")
    ax.tick_params(axis='y', labelrotation=90)

# x-axis label only on the bottom row
for ax in g.axes.flatten():
    ax.set_xlabel("")
for ax in g.axes[-1]:
    ax.set_xlabel("Year")

# legend: extract once, save separately, then hide from main figure
handles, labels = [], []

# find the first facet that actually carries legend entries
for ax in g.axes.flat:
    handles, labels = ax.get_legend_handles_labels()
    if handles:                                 # stop at the first non-empty
        break

if getattr(g, "_legend", None) is not None:
    g._legend.remove()

# remove every in-figure legend so the grid is clean
for ax in g.axes.flat:
    ax.legend_.remove() if ax.get_legend() else None

# build a tiny canvas just for the legend
leg_fig, leg_ax = plt.subplots(figsize=(2.4, 0.42 * len(labels)))
leg_ax.axis('off')
leg_ax.legend(handles, labels,
              title="Metagenome category",
              frameon=False, ncol=1,
              handlelength=2, labelspacing=1.1)
leg_fig.tight_layout()
leg_fig.savefig("../data/bubbles_legend2_smaller.png", dpi=600, bbox_inches='tight')
leg_fig.savefig("../data/bubbles_legend2.svg", bbox_inches='tight')
plt.close(leg_fig)

# export the main multi-panel plot
plt.tight_layout()
g.savefig("../data/bubbles_lowess_faceted2_smaller.png", dpi=600)
g.savefig("../data/bubbles_lowess_faceted2.svg")


# --------------------------------------------------------------------
# 5  FACETTED GRIDS – ONE FILE PER METAGENOME CATEGORY
# --------------------------------------------------------------------
from pathlib import Path
import re

out_dir = Path("../data")          # change if you prefer another folder
out_dir.mkdir(exist_ok=True)

# Helper to turn a category name into a safe file-name fragment
def slugify(text: str) -> str:
    """simple slugifier -> lower-case alpha-numeric + '_'."""
    text = text.lower()
    text = re.sub(r"[^a-z0-9]+", "_", text).strip("_")
    return text or "unknown"

for cat in sorted(dot_data["metagenome_category"].unique()):
    sub_data = dot_data.query("metagenome_category == @cat").copy()
    
    # --- a. build the basic bubble grid (same params as before, no hue) ---
    g_cat = sns.relplot(
        data     = sub_data,
        x        = "year_bin",
        y        = "prevalence",
        size     = "total_samples",
        sizes    = (20, 600),
        alpha    = .45,
        color    = color_dict[cat],     # single colour for clarity
        kind     = "scatter",
        row      = "ARO_DrugClass",
        col      = "geo_loc_name_country_continent_calc",
        facet_kws= dict(sharey=False),
        height   = 3.2,
        aspect   = 1.1,
        legend   = False               # legend not useful here (only 1 cat)
    )
    
    # --- b. LOWESS overlay (same logic as Section 3) ---
    for (row_label, col_label), ax in g_cat.axes_dict.items():
        drug = row_label
        cont = col_label
        pts  = sub_data[(sub_data["ARO_DrugClass"] == drug) &
                        (sub_data["geo_loc_name_country_continent_calc"] == cont)]
        if pts.empty:
            ax.set_axis_off()
            continue
        
        if len(pts) >= 5:              # need ≥5 for smoother curve
            smoothed = lowess(pts["prevalence"],
                              pts["year_bin"],
                              frac=0.55,
                              return_sorted=True)
            ax.plot(smoothed[:, 0], smoothed[:, 1],
                    color     = color_dict[cat],
                    linewidth = 5,
                    alpha     = .85,
                    zorder    = 1)
    
    # --- c. tidy axes & labels (mirrors Section 4) ---
    for ax in g_cat.axes.flatten():
        ax.set_ylim(0, 1.05)
    
    g_cat.set_axis_labels("Year", "Prevalence")
    g_cat.set_titles("")
    
    for ax, cont in zip(g_cat.axes[0], g_cat.col_names):
        ax.set_title(cont, fontweight="bold", pad=20)
    
    for ax, drug in zip(g_cat.axes[:, 0], g_cat.row_names):
        ax.set_ylabel(drug, rotation=90, labelpad=20,
                      fontweight="bold", va="center")
    
    for ax in g_cat.axes.flatten():
        if ax not in g_cat.axes[:, 0]:
            ax.set_ylabel("")
        ax.tick_params(axis="y", labelrotation=90)
    for ax in g_cat.axes.flatten():
        ax.set_xlabel("")
    for ax in g_cat.axes[-1]:
        ax.set_xlabel("Year")
    
    # --- d. export ----------------------------------------------------
    base = slugify(cat)
    g_cat.savefig(out_dir / f"bubbles_lowess_faceted_{base}_smaller.png", dpi=600)
    g_cat.savefig(out_dir / f"bubbles_lowess_faceted_{base}.svg")
    plt.close(g_cat.figure)          # important to free memory in long loops

print("Per-category figures saved to:", out_dir.resolve())

# --------------------------------------------------------------------
# 6  FACETTED GRIDS – ONLY WORLDWIDE DATA
# --------------------------------------------------------------------

world_dot_data = dot_data.query(
    "geo_loc_name_country_continent_calc == 'Worldwide'"
).copy()

col_order = sorted(world_dot_data["metagenome_category"].unique())

g_world = sns.relplot(
    data      = world_dot_data,
    x         = "year_bin",
    y         = "prevalence",
    hue       = "metagenome_category",
    col       = "metagenome_category",
    col_order = col_order,
    row       = "ARO_DrugClass",
    kind      = "scatter",
    size      = "total_samples",
    sizes     = (20, 600),
    alpha     = .45,
    facet_kws = dict(sharey=False),
    height    = 3.2,
    aspect    = 1.1,
    palette   = color_dict,
    legend    = "brief"          # ← identical to your earlier grid
)

# LOWESS overlay – same logic as Section 3
for (row_lab, col_lab), ax in g_world.axes_dict.items():
    drug = row_lab
    cat  = col_lab
    pts  = world_dot_data[
        (world_dot_data["ARO_DrugClass"] == drug) &
        (world_dot_data["metagenome_category"] == cat)
    ]
    if pts.empty:
        ax.set_axis_off()
        continue
    if len(pts) >= 5:
        smoothed = lowess(
            pts["prevalence"], pts["year_bin"],
            frac=0.55, return_sorted=True
        )
        ax.plot(
            smoothed[:, 0], smoothed[:, 1],
            color     = color_dict[cat],
            linewidth = 5,
            alpha     = .85,
            zorder    = 1
        )

# Cosmetic tweaks – mirrors Section 4
for ax in g_world.axes.flatten():
    ax.set_ylim(0, 1.05)

g_world.set_axis_labels("Year", "Prevalence")
g_world.set_titles("")

for ax, cat in zip(g_world.axes[0], g_world.col_names):
    ax.set_title(cat, fontweight="bold", pad=20)

for ax, drug in zip(g_world.axes[:, 0], g_world.row_names):
    ax.set_ylabel(drug, rotation=90, labelpad=20,
                  fontweight="bold", va="center")

for ax in g_world.axes.flatten():
    if ax not in g_world.axes[:, 0]:
        ax.set_ylabel("")
    ax.tick_params(axis="y", labelrotation=90)
    ax.set_xlabel("")
for ax in g_world.axes[-1]:
    ax.set_xlabel("Year")

# ------------------------ standalone legend -------------------------
handles, labels = [], []
for ax in g_world.axes.flat:
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        break

if getattr(g_world, "_legend", None) is not None:
    g_world._legend.remove()
for ax in g_world.axes.flat:
    ax.legend_.remove() if ax.get_legend() else None

# never let the height go to zero (avoids singular-matrix error)
fig_h = max(0.42 * len(labels), 0.42)

leg_fig, leg_ax = plt.subplots(figsize=(2.4, fig_h))
leg_ax.axis("off")
leg_ax.legend(
    handles, labels,
    title="Metagenome category",
    frameon=False, ncol=1,
    handlelength=2, labelspacing=1.1
)

leg_fig.tight_layout()
leg_fig.savefig("../data/bubbles_worldwide_legend.png", dpi=600, bbox_inches="tight")
leg_fig.savefig("../data/bubbles_worldwide_legend.svg", bbox_inches="tight")
plt.close(leg_fig)

# --- F. export the main figure ----------------------------------------------
plt.tight_layout()
g_world.savefig("../data/bubbles_lowess_worldwide.png", dpi=600)
g_world.savefig("../data/bubbles_lowess_worldwide.svg")