import pandas as pd
import matplotlib.pyplot as plt

# Helper: bin into quarters
def bin_by_quarter(df):
    df = df.copy()
    df['releasedate'] = pd.to_datetime(
        df['releasedate'].str.strip('[]'), errors='coerce'
    )
    df = df.dropna(subset=['releasedate'])
    df['quarter_bin'] = df['releasedate'].dt.to_period("Q")
    return df

def bin_by_third(df):
    df = df.copy()
    df['releasedate'] = pd.to_datetime(
        df['releasedate'].str.strip('[]'), errors='coerce'
    )
    df = df.dropna(subset=['releasedate'])

    def assign_third(month):
        if 1 <= month <= 4:
            return 'T1'
        elif 5 <= month <= 8:
            return 'T2'
        else:
            return 'T3'

    df['year'] = df['releasedate'].dt.year
    df['third'] = df['releasedate'].dt.month.apply(assign_third)
    df['third_bin'] = df['year'].astype(str) + '-' + df['third']
    return df

# 1. Read the data
df_total = pd.read_csv("../data/plasmids_sra_metadata_extended.csv", dtype=str)
df_amr   = pd.read_csv("../data/amr_metadata_extended.csv", dtype=str)

# 2. Filter and bin
#df_total_q = bin_by_quarter(df_total)
#df_amr_q   = bin_by_quarter(df_amr)

df_total_t = bin_by_third(df_total)
df_amr_t   = bin_by_third(df_amr)

# 3. Build unique accession sets per quarter
#total_counts = df_total_q.groupby("quarter_bin")["seq_name"].nunique()
#amr_counts   = df_amr_q.groupby("quarter_bin")["seq_name"].nunique()

total_counts = df_total_t.groupby("third_bin")["seq_name"].nunique()
amr_counts   = df_amr_t.groupby("third_bin")["seq_name"].nunique()


# Align both series to the same index
all_quarters = sorted(set(total_counts.index) | set(amr_counts.index))
total_counts = total_counts.reindex(all_quarters, fill_value=0)
amr_counts   = amr_counts.reindex(all_quarters, fill_value=0)

# 4. Plot stacked bars
fig, ax = plt.subplots(figsize=(12, 4))
bars_total = ax.bar(
    total_counts.index.astype(str), total_counts, label="Total plasmids", color="#e9c46a"
)
bars_amr = ax.bar(
    amr_counts.index.astype(str), amr_counts, bottom=total_counts - amr_counts,
    label="AMR-positive", color="#F06838"
)

# 2. Identify even-year T1 bins (for regular ticks)
even_year_pos = [
    i for i, b in enumerate(total_counts.index)
    if b.endswith('T1') and int(b.split('-')[0]) % 2 == 0
]
even_year_labs = [
    b.split('-')[0] for b in total_counts.index
    if b.endswith('T1') and int(b.split('-')[0]) % 2 == 0
]

# 3. Force add first and last bin labels
first_pos = 0
first_lab = total_counts.index[0].split('-')[0]

last_pos = len(total_counts.index) - 1
last_lab = total_counts.index[-1].split('-')[0]

# Add them if they’re not already included
tick_positions = sorted(set(even_year_pos + [first_pos, last_pos]))
tick_labels = []
for i in tick_positions:
    year = total_counts.index[i].split('-')[0]
    tick_labels.append(year)

# 4. Apply ticks
ax.set_xticks(tick_positions)
ax.set_xticklabels(tick_labels, rotation=0)


# 6. Labeling
ax.set_title("Discovery timeline of plasmids and AMR-positive subset")
ax.set_ylabel("# unique sequences")
ax.legend()

plt.tight_layout()
fig.savefig("../data/discovery_timeline_plasmids_amr_stack.png", dpi=600)
fig.savefig("../data/discovery_timeline_plasmids_amr_stack.svg")
plt.close()
