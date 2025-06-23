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

'''
## Debugging
# after add_year_bin() but BEFORE any merging / summarising
print('sra_df years     →')
print(sra_df['year_bin'].value_counts(dropna=False).sort_index())

print('\namr_df years     →')
print(amr_df['year_bin'].value_counts(dropna=False).sort_index())

# convert to numeric, but keep the original for inspection
sra_df['mbases_num'] = pd.to_numeric(sra_df['mbases'], errors='coerce')

bad = sra_df['mbases_num'].isna() | (sra_df['mbases_num'] <= 0)

print('Rows with bad / zero mbases (per year) →')
print(sra_df.loc[bad, 'year_bin'].value_counts().sort_index())
'''

# --------------------------------------------------------------------
# 1 PLOT DISCOVERY TIMELINE OF AMRs ON SAMPLES WITHIN METAGENOME_CATEGORY
# --------------------------------------------------------------------

per_quarter = (
    amr_df
      .groupby("quarter_bin")["acc"]
      .nunique()                 # unique AMR-positive samples
      .sort_index()              # ensures chronological order
)

fig, ax = plt.subplots(figsize=(12, 4))
per_quarter.plot.bar(ax=ax, width=0.9)

ax.set_title("Discovery timeline of AMR-positive metagenome samples")
ax.set_xlabel("")
ax.set_ylabel("# SRA accessions AMR-positive")

# Years need to be changed from 2008Q1 to 2008, 2009, 2010, etc
year_start_pos  = [i for i, p in enumerate(per_quarter.index) if p.quarter == 1]
year_start_labs = [str(p.year) for p in per_quarter.index if p.quarter == 1]

ax.set_xticks(year_start_pos)
ax.set_xticklabels(year_start_labs, rotation=0)

plt.tight_layout()
fig.savefig("../data/discovery_timeline_amr.png", dpi=600)
fig.savefig("../data/discovery_timeline_amr.svg")
plt.close()

# --------------------------------------------------------------------
# 1 SAME, BUT COLOR STACK BARS BY METAGENOME CATEGORY
# --------------------------------------------------------------------

quarter_cat = (
    amr_df
      .groupby(["quarter_bin", "metagenome_category"])["acc"]
      .nunique()               # count distinct accessions
      .unstack(fill_value=0)   # columns = categories
      .sort_index()            # chronological order
)

n_cols  = quarter_cat.shape[1]
cblind  = sns.color_palette("colorblind", n_colors=n_cols)

fig, ax = plt.subplots(figsize=(12, 4))
quarter_cat.plot(kind="bar", stacked=True, ax=ax, width=0.9, color=cblind)

ax.set_title("Discovery timeline of AMR-positive metagenome samples\n(coloured by metagenome category)")
ax.set_xlabel("")                     # custom x-ticks below
ax.set_ylabel("# SRA accessions AMR-positive")

# --- show only the first quarter of each year on the x-axis -----------
period_idx = quarter_cat.index        # PeriodIndex like 2008Q1 …
year_start_pos  = [i for i, p in enumerate(period_idx) if p.quarter == 1]
year_start_labs = [str(p.year) for p in period_idx if p.quarter == 1]

ax.set_xticks(year_start_pos)
ax.set_xticklabels(year_start_labs, rotation=0)

ax.legend(title="Metagenome category", bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
fig.savefig("../data/discovery_timeline_amr_categories.png", dpi=600)
fig.savefig("../data/discovery_timeline_amr_categories.svg")


