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

# --------------------------------------------------------------------
# 1  COUNT UNIQUE ACCESSIONS  (total vs. positives)
# --------------------------------------------------------------------
# Total samples per (quarter, category)
tot_counts = (
    sra_df
      .groupby(["quarter_bin", "metagenome_category"])["acc"] 
      .nunique()
)

# Positives per (quarter, category)
pos_counts = (
    amr_df
      .groupby(["quarter_bin", "metagenome_category"])["acc"]
      .nunique()
)

# Total db for stacking
tot_df = (
    tot_counts
      .unstack(fill_value=0)
      .sort_index()
)

# Pos db for plotting
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
# 2  PLOTS  – loop once per metagenome category
# --------------------------------------------------------------------
for cat in categories:
    ## Discovery timeline (positives only)
    fig, ax = plt.subplots(figsize=(10, 3))
    pos_df[cat].plot.bar(ax=ax, width=0.9, color=cat2color[cat])

    ax.set_title(f"{cat} AMR-positive discovery timeline")
    ax.set_xlabel("")
    ax.set_ylabel("Unique AMR-positive accessions")

    period_idx = pos_df.index
    year_pos   = [i for i, p in enumerate(period_idx) if p.quarter == 1]
    year_labs  = [str(p.year) for p in period_idx if p.quarter == 1]

    ax.set_xticks(year_pos)
    ax.set_xticklabels(year_labs, rotation=0)
    plt.tight_layout()
    fig.savefig(f"../data/discovery_timeline_{cat}_positive.png", dpi=600)
    fig.savefig(f"../data/discovery_timeline_{cat}_positive.svg")
    plt.close()

    ## Stacked: total (grey) and positive (colour)
    fig, ax = plt.subplots(figsize=(10, 3))
    non_pos = tot_df[cat] - pos_df[cat]

    ax.bar(tot_df.index.astype(str), non_pos, color="#CCCCCC", width=0.9, label="Total sampling")
    ax.bar(tot_df.index.astype(str), pos_df[cat], bottom=non_pos,
           color=cat2color[cat], width=0.9, label="AMR-positive")

    ax.set_title(f"{cat} SRA total and AMR-positive timeline")
    ax.set_xlabel("")
    ax.set_ylabel("Unique accessions")

    ax.set_xticks(year_pos)
    ax.set_xticklabels(year_labs, rotation=0)
    ax.legend(frameon=False, loc="upper left")
    plt.tight_layout()
    fig.savefig(f"../data/discovery_timeline_{cat}_stacked.png", dpi=600)
    fig.savefig(f"../data/discovery_timeline_{cat}_stacked.svg")
    plt.close()