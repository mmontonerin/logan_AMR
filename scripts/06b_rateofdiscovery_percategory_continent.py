import pandas as pd
import seaborn as sns # only for the palette used in circles plot
import matplotlib.pyplot as plt
import numpy as np

# --------------------------------------------------------------------
# 0 PRE-PROCESSING DATA
# --------------------------------------------------------------------

# Input table
amr_df = pd.read_csv(
    "../data/full_card_metadata_aro_allfilters_metagenomes2.csv",
    low_memory=False)

sra_df = pd.read_csv("../data/SRA_metadata_allfilters.csv", low_memory=False)

# Function to prepare tables to have date as year bins
def date_collection(df):
    df = df.copy()
    df['collection_date_sam'] = pd.to_datetime(
        df['collection_date_sam'].str.strip('[]'),
        errors='coerce')
    df = df.dropna(subset=['collection_date_sam'])
    df['quarter_bin'] = (
        df['collection_date_sam']
            .dt.to_period('Q') # Generate 4 bins per year
    )
    return df

amr_df = date_collection(amr_df)
sra_df = date_collection(sra_df)

# Positive samples vs Total
tot_cat_df = (
    sra_df
      .groupby(["quarter_bin", "metagenome_category"])["acc"]
      .nunique()
      .unstack(fill_value=0)          # quarters × category
      .sort_index()
)

pos_cat_df = (
    amr_df
      .groupby(["quarter_bin", "metagenome_category"])["acc"]
      .nunique()
      .unstack(fill_value=0)
      .reindex(index=tot_cat_df.index, columns=tot_cat_df.columns, fill_value=0)
)

# Continent information tables
continent_order = [
    "Africa", "Asia", "Europe",
    "North America", "Oceania", "South America"
]
# Rainbow palette
continent_colors = dict(zip(
    continent_order,
    ["#ff595e", "#ff924c", "#ffca3a", "#8ac926", "#1982c4", "#6a4c93"]
))

def agg_quarter_continent(df, category):
    out = (
        df[df["metagenome_category"] == category]
        .groupby(["quarter_bin", "geo_loc_name_country_continent_calc"])["acc"]
        .nunique()
        .unstack(fill_value=0)
        .reindex(columns=continent_order, fill_value=0)
        .sort_index()
    )
    return out


# --------------------------------------------------------------------
# 1  PLOTS  – loop once per metagenome category
# --------------------------------------------------------------------
categories = sorted(amr_df["metagenome_category"].dropna().unique())

# One fixed colour per metagenome category (same colour-blind palette as bubbles plot)
palette    = sns.color_palette("colorblind", n_colors=len(categories))
cat2color  = dict(zip(categories, palette))

for cat in categories:
    # positives only
    pos_df  = agg_quarter_continent(amr_df, cat)

    fig, ax = plt.subplots(figsize=(11, 3.5))
    pos_df.plot(kind="bar", stacked=True, ax=ax, width=0.9,
                color=[continent_colors[c] for c in pos_df.columns])

    ax.set_title(f"{cat} AMR-positive discovery timeline (stacked by continent)")
    ax.set_xlabel("")
    ax.set_ylabel("Unique AMR-positive accessions")

    period_idx = pos_df.index
    year_pos   = [i for i, p in enumerate(period_idx) if p.quarter == 1]
    year_labs  = [str(p.year) for p in period_idx if p.quarter == 1]

    ax.set_xticks(year_pos)
    ax.set_xticklabels(year_labs, rotation=0)
    ax.legend(title="Continent", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    fig.savefig(f"../data/discovery_timeline_{cat}_positive_continents.png", dpi=600)
    fig.savefig(f"../data/discovery_timeline_{cat}_positive_continents.svg")
    plt.close()

    # Total and positives
    # LEFT BAR  (totals vs positives)
    # draw grey for every quarter; draw coloured layer only if >0
    # always draw the grey layer
    tot_cont_df = agg_quarter_continent(sra_df, cat)
    fig, ax = plt.subplots(figsize=(11, 3.5))

    aligned_idx = tot_cont_df.index          # master index for this figure
    x = np.arange(len(aligned_idx))
    w = 0.42

    # 1) align the Series
    pos_tot = pos_cat_df[cat].reindex(aligned_idx, fill_value=0)
    non_pos = (tot_cat_df[cat].reindex(aligned_idx, fill_value=0) - pos_tot)

    # 2) always draw the grey layer
    ax.bar(x - w/2, non_pos.values,
       width=w, color="#CCCCCC",
       label="No AMR detection samples")

    # 3) coloured overlay only where positive
    mask = pos_tot.values > 0
    ax.bar(x[mask] - w/2, pos_tot.values[mask],
       width=w, bottom=non_pos.values[mask],
       color=cat2color[cat], label="AMR-positive")

    # ----- RIGHT BAR  (continents, totals) ----------------------
    bottoms = np.zeros(len(x))
    for cont in continent_order:
        ax.bar(x + w/2, tot_cont_df[cont], width=w,
               bottom=bottoms, color=continent_colors[cont],
               label=cont if x[0] == 0 else None)
        bottoms += tot_cont_df[cont].to_numpy()

    # ---- Aesthetics -------------------------------------------
    ax.set_title(f"{cat} SRA total and AMR-positive timeline\nLeft: total vs AMR-positive  |  Right: total samples by continent")
    ax.set_xlabel("")
    ax.set_ylabel("Unique accessions")

    period_idx = aligned_idx
    year_pos   = [i for i, p in enumerate(period_idx) if p.quarter == 1]
    year_labs  = [str(p.year) for p in period_idx if p.quarter == 1]
    ax.set_xticks(year_pos)
    ax.set_xticklabels(year_labs, rotation=0)

    ax.set_xlim(-1, len(x))      # a little padding
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    fig.savefig(f"../data/discovery_timeline_{cat}_total_vs_positive_continents.png", dpi=600)
    fig.savefig(f"../data/discovery_timeline_{cat}_total_vs_positive_continents.svg")
    plt.close()
