import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

# Function to prepare tables to have date as year bins
def date_collection(df: pd.DataFrame) -> pd.DataFrame:
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

# Function to classify sample into Metagenome or Isolate
def classify_sample(df: pd.DataFrame) -> pd.DataFrame:
    # Add 'sample_type' = 'Metagenome' or 'Isolate'
    is_meta = df["organism"].str.contains("metagenome", case=False, na=False)
    df = df.copy()
    df["sample_type"] = is_meta.map({True: "Metagenome", False: "Isolate"})
    return df

# --------------------------------------------------------------------
# 0 PROCESSING DATA
# --------------------------------------------------------------------

# collect accessions, one to collapse only on quarter, and another also on sample type
seen = defaultdict(set)
totals = defaultdict(set)

# Input table
for chunk in pd.read_csv(
    "../data/new_full_card_metadata_aro_dateloc_organism.csv",
    chunksize=80000000,
    low_memory=False
    ):

    chunk = date_collection(chunk)
    chunk = classify_sample(chunk)

    grouped = (
        chunk.groupby(["quarter_bin", "sample_type"])["acc"]
              .apply(set)
    )

    for key, acc_set in grouped.items():
        seen[key].update(acc_set)

# --------------------------------------------------------------------
# 1 PLOT DISCOVERY TIMELINE OF AMRs ON ALL SAMPLES
# --------------------------------------------------------------------

# seen keys are (quarter_bin, sample_type).  Collapse on quarter only
totals = defaultdict(set)

for (q, _type), accs in seen.items():
    totals[q].update(accs)          # merge sets so accessions stay unique

per_quarter = (
    pd.Series({q: len(accs) for q, accs in totals.items()}, name="unique_accessions")
      .sort_index()
)

fig, ax = plt.subplots(figsize=(12, 4))
per_quarter.plot.bar(ax=ax, width=0.9)

ax.set_title("Discovery timeline of AMR-positive total samples")
ax.set_xlabel("")
ax.set_ylabel("# SRA accessions AMR-positive")

# Years need to be changed from 2008Q1 to 2008, 2009, 2010, etc
year_start_pos  = [i for i, p in enumerate(per_quarter.index) if p.quarter == 1]
year_start_labs = [str(p.year) for p in per_quarter.index if p.quarter == 1]

ax.set_xticks(year_start_pos)
ax.set_xticklabels(year_start_labs, rotation=0)

plt.tight_layout()
fig.savefig("../data/discovery_timeline_amr_total.png", dpi=600)
fig.savefig("../data/discovery_timeline_amr_total.svg")
plt.close()

# --------------------------------------------------------------------
# 1 SAME, BUT COLOR STACK BARS BY ORGANISM CATEGORY (METAGENOME OR ISOLATE)
# --------------------------------------------------------------------

records = [
    {"quarter_bin": q, "sample_type": t, "unique_accessions": len(accs)}
    for (q, t), accs in seen.items()
]

quarter_type = (
    pd.DataFrame(records)
      .pivot(index="quarter_bin", columns="sample_type", values="unique_accessions")
      .fillna(0)
      .astype(int)
      .sort_index()
)

for col in ["Metagenome", "Isolate"]:
    quarter_type[col] = quarter_type.get(col, 0)



colors = {"Metagenome": "#FDBF6F",   # soft orange
          "Isolate"   : "#C2B2FF"}   # soft purple

fig, ax = plt.subplots(figsize=(12, 4))
quarter_type[["Metagenome", "Isolate"]].plot(
    kind="bar",
    stacked=True,
    ax=ax,
    width=0.9,
    color=[colors["Metagenome"], colors["Isolate"]],
)

ax.set_title("Discovery timeline of AMR-positive total samples\n(metagenome vs. isolate)")
ax.set_xlabel("")
ax.set_ylabel("# SRA accessions AMR-positive")

# Only show year labels (first quarter of each year)
period_idx = quarter_type.index
year_pos  = [i for i, p in enumerate(period_idx) if p.quarter == 1]
year_labs = [str(p.year) for p in period_idx if p.quarter == 1]

ax.set_xticks(year_pos)
ax.set_xticklabels(year_labs, rotation=0)
ax.legend(title="Sample type", bbox_to_anchor=(1.02, 1), loc="upper left")

plt.tight_layout()
fig.savefig("../data/discovery_timeline_amr_total_org.png", dpi=600)
fig.savefig("../data/discovery_timeline_amr_total_org.svg")
plt.close()


