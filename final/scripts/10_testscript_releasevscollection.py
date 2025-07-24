import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind

# File paths
sra_path = "../data/SRA_metadata_before20231211_logan_extended.csv"
card_path = "../data/card_metadata_aro_extended.csv"

# Parameters
chunksize = 30_000_000
date_cols = ['collection_date_sam', 'releasedate']
organism_filter = ['Isolate', 'Metagenome']
colors = {"Metagenome": "#FDBF6F", "Isolate": "#C2B2FF"}

# Function to process a chunk
def process_chunk(df, source):
    df = df.copy()
    df = df.dropna(subset=date_cols)
    df['collection_date_sam'] = df['collection_date_sam'].astype(str).str.strip("[]").str.strip() # collection date has brackets
    df['collection_date_sam'] = pd.to_datetime(df['collection_date_sam'], errors='coerce', format='%Y-%m-%d')
    df['releasedate'] = pd.to_datetime(df['releasedate'], errors='coerce', format='%Y-%m-%d')
    df = df.dropna(subset=date_cols)
    df['days_to_release'] = (df['releasedate'] - df['collection_date_sam']).dt.days
    df = df[df['days_to_release'] >= 0]  # Filter out negative days to release
    df = df[df['organism_type'].isin(organism_filter)]
    df['source'] = source
    return df[['acc', 'organism_type', 'days_to_release', 'source']]

# Read SRA in chunks
sra_chunks = pd.read_csv(sra_path, chunksize=chunksize, usecols=lambda col: col in ['acc'] + date_cols + ['organism_type'], low_memory=False)
sra_processed = pd.concat([process_chunk(chunk, 'SRA') for chunk in sra_chunks], ignore_index=True)

# Read CARD in chunks, deduplicating by 'acc'
card_seen = set()
card_rows = []
card_chunks = pd.read_csv(card_path, chunksize=chunksize, usecols=lambda col: col in ['acc'] + date_cols + ['organism_type'], low_memory=False)

for chunk in card_chunks:
    chunk = chunk.dropna(subset=['acc'])
    chunk = chunk[~chunk['acc'].isin(card_seen)]
    card_seen.update(chunk['acc'])
    processed = process_chunk(chunk, 'CARD')
    card_rows.append(processed)

card_processed = pd.concat(card_rows, ignore_index=True)

# Plotting function
def plot_by_source(df, title):
    isolate_days = df[df['organism_type'] == 'Isolate']['days_to_release']
    metagenome_days = df[df['organism_type'] == 'Metagenome']['days_to_release']
    tstat, pval = ttest_ind(isolate_days, metagenome_days, equal_var=False)

    # Add horizontal black lines for the mean
    for i, org_type in enumerate(['Isolate', 'Metagenome']):
        vals = df[df['organism_type'] == org_type]['days_to_release']
        mean = vals.mean()
        plt.hlines(mean, i - 0.2, i + 0.2, color='black', linewidth=2)

    print (title)
    print("Mean Isolate days to release:", isolate_days.mean())
    print("Mean Metagenome days to release:", metagenome_days.mean())
    print("Mean difference (Metagenome - Isolate):", metagenome_days.mean() - isolate_days.mean())

    significance = "*" if pval < 0.05 else "ns"
    label = f"{significance}  (t = {tstat:.2f}, p = {pval:.2e})"

    plt.figure(figsize=(8, 6))

    sns.boxplot(
        data=df,
        x='organism_type',
        y='days_to_release',
        palette=colors,
        showfliers=True,
    )

    # Annotation
    ymax = df['days_to_release'].max()
    plt.text(0.5, ymax + 100, label, ha='center', va='bottom', fontsize=12)

    plt.ylabel("Days from collection to release")
    plt.xlabel("Organism type")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(f"../data/{title.replace(' ', '_').lower()}2.png", dpi=600)
    plt.savefig(f"../data/{title.replace(' ', '_').lower()}2.svg")
    plt.close()

# Plot each separately
plot_by_source(sra_processed, "SRA Samples: Time to Release by Organism Type")
plot_by_source(card_processed, "CARD Samples: Time to Release by Organism Type")


