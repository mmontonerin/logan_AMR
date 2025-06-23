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

# Continent information tables
continent_order = [
    "Africa", "Asia", "Europe",
    "North America", "Oceania", "South America"
]
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

for cat in categories:
    # ---- positives only (AMR) ----
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

    # ---- total vs positive + continent stack (side-by-side) ----
    tot_df = agg_quarter_continent(sra_df, cat)      # all samples
    pos_tot = pos_df.sum(axis=1)                     # total positives (no continent split)
    non_pos = tot_df.sum(axis=1) - pos_tot           # total sampling

    # align indices
    pos_tot = pos_tot.reindex(index=tot_df.index, fill_value=0)
    non_pos = non_pos.reindex(index=tot_df.index, fill_value=0)

    x = np.arange(len(tot_df.index))
    w = 0.42

    fig, ax = plt.subplots(figsize=(11, 3.5))

    # ----- LEFT BAR  (totals vs positives) -----------------------
    ax.bar(x - w/2, non_pos,   width=w, color="#CCCCCC",
           label="No AMR detection samples")
    ax.bar(x - w/2, pos_tot,   width=w, bottom=non_pos,
           color=cat2color[cat], label="AMR-positive")

    # ----- RIGHT BAR  (continents, totals) ----------------------
    bottoms = np.zeros(len(x))
    for cont in continent_order:
        ax.bar(x + w/2, tot_df[cont], width=w,
               bottom=bottoms, color=continent_colors[cont],
               label=cont if x[0] == 0 else None)
        bottoms += tot_df[cont].to_numpy()

    # ---- Aesthetics -------------------------------------------
    ax.set_title(f"{cat} SRA total and AMR-positive timeline\nLeft: total vs AMR-positive  |  Right: total samples by continent")
    ax.set_xlabel("")
    ax.set_ylabel("Unique accessions")

    period_idx = tot_df.index
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
