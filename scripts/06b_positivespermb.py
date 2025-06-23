import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# --------------------------------------------------------------------
# 0 PRE-PROCESSING DATA
# --------------------------------------------------------------------

# Input table
amr_df = pd.read_csv(
    "../data/full_card_metadata_aro_allfilters_metagenomes2.csv",
    low_memory=False)

sra_df = pd.read_csv("../data/SRA_metadata_allfilters.csv", 
    low_memory=False)

# Function to prepare tables to have date as year bins
def date_collection(df):
    df = df.copy()
    df['collection_date_sam'] = pd.to_datetime(
        df['collection_date_sam'].str.strip('[]'),
        errors='coerce')
    df = df.dropna(subset=['collection_date_sam'])
    df['year_bin'] = df['collection_date_sam'].dt.year.astype('Int64') # Generate yearly bins
    return df

def clean_mb(df): # Keep only rows with valid mbases > 0
    df = df.copy()
    df['mbases'] = pd.to_numeric(df['mbases'], errors='coerce')
    df = df[df['mbases'].notna() & (df['mbases'] > 0)]
    return df

amr_df = date_collection(amr_df)
sra_df = date_collection(sra_df)
sra_df = clean_mb(sra_df)

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
# 1 COUNT TOTAL MB PER YEAR BIN AND POSITIVE SAMPLES PER MB
# --------------------------------------------------------------------

SAMPLE_ID = 'acc'
SAMPLE_SIZE = 'mbases'

# count hits per sample (amr_df may have ≥1 row per hit)
hits_per_sample = (amr_df
                   .groupby(SAMPLE_ID)
                   .size()
                   .rename('hit_count')
                   .reset_index())

# bring in mbases for every sample; samples with no hits get hit_count = 0
samples = (sra_df[[SAMPLE_ID, SAMPLE_SIZE, 'year_bin']]
           .merge(hits_per_sample, on=SAMPLE_ID, how='left')
           .fillna({'hit_count': 0}))

# rate per sample
samples['hits_per_mb'] = samples['hit_count'] / samples[SAMPLE_SIZE]

# --------------------------------------------------------------------
# 2 STATS PER YEAR BIN
# --------------------------------------------------------------------
def summarise(group):
    n  = len(group)
    mu = group['hits_per_mb'].mean()
    sd = group['hits_per_mb'].std(ddof=1)
    se = sd / np.sqrt(n)
    return pd.Series({'mean': mu, 'sd': sd, 'se': se, 'n': n})

all_stats = (
    samples
    .groupby('year_bin')
    .apply(summarise, include_groups=False)
)

pos_stats = (
    samples[samples['hit_count'] > 0]
    .groupby('year_bin')
    .apply(summarise, include_groups=False)
)

# --------------------------------------------------------------------
# 3  PLOT
# --------------------------------------------------------------------
years      = all_stats.index.astype(int)
x          = np.arange(len(years))
width      = 0.35  # bar width

fig, ax = plt.subplots(figsize=(10, 5))

# bars for ALL samples
ax.bar(x - width/2,
       all_stats['mean'],
       width,
       yerr=all_stats['se'],
       capsize=3,
       label='All samples',
       color='#747474')

# bars for POSITIVE-ONLY samples
ax.bar(x + width/2,
       pos_stats['mean'],
       width,
       yerr=pos_stats['se'],
       capsize=3,
       label='Positive samples',
       color='#eb7c7c')

# aesthetics
ax.set_xlabel('Collection year')
ax.set_ylabel('Positive hits per Mb')
ax.set_xticks(x)
ax.set_xticklabels(years, rotation=45)
ax.set_title('Average positive-hit rate by year')
ax.legend()
fig.tight_layout()

fig.savefig("../data/positives_per_mb.png", dpi=600)
fig.savefig("../data/positives_per_mb.svg")


# --------------------------------------------------------------------
# 4  Inspect SD / sample counts
# --------------------------------------------------------------------
summary = (pd.concat({'all': all_stats, 'positive': pos_stats},
                     axis=1, names=['group', None])
           .swaplevel(axis=1))

summary.to_csv("../data/positives_per_mb_summary.csv")
