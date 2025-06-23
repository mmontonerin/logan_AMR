import pandas as pd
import seaborn as sns          # only for the palette used in circles plot
import matplotlib.pyplot as plt
import numpy as np

# --------------------------------------------------------------------
# 0 PRE-PROCESSING DATA
# --------------------------------------------------------------------

# Input table
amr_df = pd.read_csv(
    "../data/full_card_metadata_aro_allfilters_metagenomes2.csv",
    low_memory=False)

#sra_df = pd.read_csv("../data/SRA_metadata_allfilters.csv", low_memory=False) # Old one

sra_df = pd.read_csv("../data/SRA_metadata_allfilters_logan.csv", low_memory=False) # New one with only SRA entries that were assembled in Logan

def date_collection(df):
    df = df.copy()
    df["collection_date_sam"] = pd.to_datetime(
        df["collection_date_sam"].str.strip("[]"), errors="coerce"
    )
    df = df.dropna(subset=["collection_date_sam"])
    df["quarter_bin"] = df["collection_date_sam"].dt.to_period("Q")
    return df

amr_df = date_collection(amr_df)
sra_df = date_collection(sra_df)

# --------------------------------------------------------------------
# 1  USE BioProject AS THE COUNTING ELEMENT
# --------------------------------------------------------------------
BIOPROJECT_COL = "bioproject"

tot_cat_df = (
    sra_df
      .groupby(["quarter_bin", "metagenome_category"])[BIOPROJECT_COL]
      .nunique() # unique BioProjects per quarter & cat
      .unstack(fill_value=0)
      .sort_index()
)

pos_cat_df = (
    amr_df
      .groupby(["quarter_bin", "metagenome_category"])[BIOPROJECT_COL]
      .nunique()
      .unstack(fill_value=0)
      .reindex(index=tot_cat_df.index, columns=tot_cat_df.columns, fill_value=0)
)

# --------------------------------------------------------------------
# 2  CONTINENT TABLE BUILDER
# --------------------------------------------------------------------
continent_order = [
    "Africa", "Asia", "Europe",
    "North America", "Oceania", "South America"
]
continent_colors = dict(zip(
    continent_order,
    ["#ff595e", "#ff924c", "#ffca3a", "#8ac926", "#1982c4", "#6a4c93"]
))

def agg_quarter_continent(df, category):
    """
    Return a quarters per continent table of BioProject *weights*.

    Each BioProject contributes a total weight of 1 per quarter,
    divided equally among the continents it appears in that quarter.
    """
    subset = (
        df[df["metagenome_category"] == category]
          .loc[:, ["quarter_bin",
                   "bioproject",
                   "geo_loc_name_country_continent_calc"]]
    )

    # For every (quarter, BioProject) pair, count its continents
    continent_counts = (
        subset.drop_duplicates()                               # unique rows
              .groupby(["quarter_bin", "bioproject"])
              ["geo_loc_name_country_continent_calc"]
              .nunique()
              .rename("n_cont")
    )

    # Merge back to compute the fractional weight = 1 / n_cont
    subset = subset.drop_duplicates().merge(
        continent_counts,
        on=["quarter_bin", "bioproject"],
        how="left"
    )
    subset["weight"] = 1 / subset["n_cont"]

    out = (
        subset.groupby(
                ["quarter_bin", "geo_loc_name_country_continent_calc"]
              )["weight"]
              .sum()
              .unstack(fill_value=0)
              .reindex(columns=continent_order, fill_value=0)
              .sort_index()
    )
    return out

# --------------------------------------------------------------------
# 3  PLOTS  – loop once per metagenome category
# --------------------------------------------------------------------
categories = sorted(amr_df["metagenome_category"].dropna().unique())
palette    = sns.color_palette("colorblind", n_colors=len(categories))
cat2color  = dict(zip(categories, palette))

for cat in categories:
    # ---------- POSITIVE DISCOVERY   (stacked by continent) ----------
    pos_df = agg_quarter_continent(amr_df, cat)
    fig, ax = plt.subplots(figsize=(11, 3.5))
    pos_df.plot(kind="bar", stacked=True, width=0.9, ax=ax,
                color=[continent_colors[c] for c in pos_df.columns])

    ax.set_title(f"{cat} AMR-positive BioProjects")
    ax.set_ylabel("Unique BioProjects")
    ax.set_xlabel("")
    year_pos  = [i for i, p in enumerate(pos_df.index) if p.quarter == 1]
    year_labs = [str(p.year) for p in pos_df.index if p.quarter == 1]
    ax.set_xticks(year_pos); ax.set_xticklabels(year_labs, rotation=0)
    ax.legend(title="Continent", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    fig.savefig(f"../data/bioproject_positive_{cat}_continents.png", dpi=600)
    fig.savefig(f"../data/bioproject_positive_{cat}_continents.svg")
    plt.close()

    # ---------- TOTAL VS POSITIVE  +  CONTINENT STACK ---------------
    tot_cont_df = agg_quarter_continent(sra_df, cat)      # totals by continent
    aligned_idx = tot_cont_df.index
    x = np.arange(len(aligned_idx)); w = 0.42

    pos_tot = pos_cat_df[cat].reindex(aligned_idx, fill_value=0)
    non_pos = tot_cat_df[cat].reindex(aligned_idx, fill_value=0) - pos_tot

    fig, ax = plt.subplots(figsize=(11, 3.5))

    # LEFT bar  (grey total + coloured positives)
    ax.bar(x - w/2, non_pos.values,
           width=w, color="#CCCCCC", label="AMR-negative")
    mask = pos_tot.values > 0
    ax.bar(x[mask] - w/2, pos_tot.values[mask],
           width=w, bottom=non_pos.values[mask],
           color=cat2color[cat], label="AMR-positive")

    # RIGHT bar  (continent stack of total BioProjects)
    bottoms = np.zeros(len(x))
    for cont in continent_order:
        ax.bar(x + w/2, tot_cont_df[cont], width=w,
               bottom=bottoms, color=continent_colors[cont],
               label=cont if x[0] == 0 else None)
        bottoms += tot_cont_df[cont].to_numpy()

    # cosmetics ------------------------------------------------------
    year_pos  = [i for i, p in enumerate(aligned_idx) if p.quarter == 1]
    year_labs = [str(p.year) for p in aligned_idx if p.quarter == 1]
    ax.set_xticks(year_pos); ax.set_xticklabels(year_labs, rotation=0)
    ax.set_title(f"{cat} BioProjects\nLeft: total vs AMR-positive  |  Right: total by continent")
    ax.set_ylabel("Unique BioProjects")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    fig.savefig(f"../data/bioproject_discovery_timeline_total_vs_positive_{cat}.png", dpi=600)
    fig.savefig(f"../data/bioproject_discovery_timeline_total_vs_positive_{cat}.svg")
    plt.close()

############################################################################
## ADD plot per category to the same script
tot_counts = (
    sra_df
      .groupby(["quarter_bin", "metagenome_category"])[BIOPROJECT_COL]
      .nunique()                        # unique BioProjects
)

pos_counts = (
    amr_df
      .groupby(["quarter_bin", "metagenome_category"])[BIOPROJECT_COL]
      .nunique()
)

tot_df = (
    tot_counts
      .unstack(fill_value=0)
      .sort_index()
)

pos_df = (
    pos_counts
      .unstack(fill_value=0)
      .reindex(index=tot_df.index, columns=tot_df.columns, fill_value=0)
)

# Colour-blind palette, one colour per category
categories     = tot_df.columns.tolist()
palette         = sns.color_palette("colorblind", n_colors=len(categories))
cat2color       = dict(zip(categories, palette))

# --------------------------------------------------------------------
# PLOTS  – loop once per metagenome category
# --------------------------------------------------------------------
for cat in categories:
    ## Discovery timeline (positives only)
    fig, ax = plt.subplots(figsize=(10, 3))
    pos_df[cat].plot.bar(ax=ax, width=0.9, color=cat2color[cat])

    ax.set_title(f"{cat} AMR-positive BioProjects discovery timeline")
    ax.set_xlabel("")
    ax.set_ylabel("Unique AMR-positive BioProjects")

    period_idx = pos_df.index
    year_pos   = [i for i, p in enumerate(period_idx) if p.quarter == 1]
    year_labs  = [str(p.year) for p in period_idx if p.quarter == 1]

    ax.set_xticks(year_pos)
    ax.set_xticklabels(year_labs, rotation=0)
    plt.tight_layout()
    fig.savefig(f"../data/bioproject_discovery_timeline_{cat}_positive_logan.png", dpi=600)
    fig.savefig(f"../data/bioproject_discovery_timeline_{cat}_positive_logan.svg")
    plt.close()

    ## Stacked: total (grey) and positive (colour)
    fig, ax = plt.subplots(figsize=(10, 3))
    pos_tot = pos_df[cat]                 
    non_pos = tot_df[cat] - pos_tot # grey layer = total – positive

    x = np.arange(len(tot_df.index))
    w = 0.42

    # draw grey for every quarter; draw coloured layer only if >0
    # always draw the grey layer
    ax.bar(x,
       non_pos.values,
       color="#CCCCCC",
       label="AMR-negative")

    # draw the coloured overlay only on quarters that have ≥1 positive
    mask = pos_tot.values > 0          # boolean array same length as x
    ax.bar(x[mask],
       pos_tot.values[mask],       # positive heights for those quarters
       bottom=non_pos.values[mask],
       color=cat2color[cat],
       label="AMR-positive")

    ax.set_title(f"{cat} BioProjects total and AMR-positive timeline")
    ax.set_xlabel("")
    ax.set_ylabel("Unique BioProjects")

    ax.set_xticks(year_pos)
    ax.set_xticklabels(year_labs, rotation=0)
    ax.legend(frameon=False, loc="upper left")
    plt.tight_layout()
    fig.savefig(f"../data/bioproject_discovery_timeline_{cat}_stacked_logan.png", dpi=600)
    fig.savefig(f"../data/bioproject_discovery_timeline_{cat}_stacked_logan.svg")
    plt.close()