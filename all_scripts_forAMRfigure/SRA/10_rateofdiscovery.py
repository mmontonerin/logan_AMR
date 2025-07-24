import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
'''
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

# Function to prepare tables to have date as year bins (release date instead of collection date)
def date_collection(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df['releasedate'] = pd.to_datetime(
        df['releasedate'],
        errors='coerce')
    df = df.dropna(subset=['releasedate'])
    df['quarter_bin'] = (
        df['releasedate']
            .dt.to_period('Q') # Generate 4 bins per year
    )
    return df
'''
# Thirds of year instead
def date_collection(df):
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

# --------------------------------------------------------------------
# 0 PROCESSING DATA
# --------------------------------------------------------------------

# collect accessions, one to collapse only on quarter, and another also on sample type
seen = defaultdict(set)
totals = defaultdict(set)

# Input table
for chunk in pd.read_csv(
    "../data/card_metadata_aro_extended.csv",
    chunksize=3000000,
    low_memory=False
    ):

#    chunk = chunk[chunk["collection_date_sam"].notna()]
#   chunk = chunk[chunk["metagenome_category"].notna()]
    chunk = chunk[chunk["releasedate"].notna()]
    chunk = date_collection(chunk)

#    grouped = (
#        chunk.groupby(["quarter_bin", "metagenome_category"])["acc"]
#              .apply(set)
#   )   

 #   grouped = (
 #       chunk.groupby(["quarter_bin", "organism_type"])["acc"]
 #             .apply(set)
 #   )

    grouped = (
    chunk.groupby(["third_bin"])["acc"]
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
per_quarter.plot.bar(ax=ax, width=0.9, color="#e9c46a")

ax.set_title("Discovery timeline of AMR-positive total samples")
ax.set_xlabel("")
ax.set_ylabel("# SRA accessions AMR-positive")

# 2. Identify even-year T1 bins (for regular ticks)
even_year_pos = [
    i for i, b in enumerate(totals.index)
    if b.endswith('T1') and int(b.split('-')[0]) % 2 == 0
]
even_year_labs = [
    b.split('-')[0] for b in totals.index
    if b.endswith('T1') and int(b.split('-')[0]) % 2 == 0
]

# 3. Force add first and last bin labels
first_pos = 0
first_lab = totals.index[0].split('-')[0]

last_pos = len(totals.index) - 1
last_lab = totals.index[-1].split('-')[0]

# Add them if they’re not already included
tick_positions = sorted(set(even_year_pos + [first_pos, last_pos]))
tick_labels = []
for i in tick_positions:
    year = totals.index[i].split('-')[0]
    tick_labels.append(year)

# 4. Apply ticks
ax.set_xticks(tick_positions)
ax.set_xticklabels(tick_labels, rotation=0)

plt.tight_layout()
fig.savefig("../data/discovery_timeline_amr_total_yellow_releasedate.png", dpi=600)
fig.savefig("../data/discovery_timeline_amr_total_yellow_releasedate.svg")
plt.close()
'''
# --------------------------------------------------------------------
# 1 SAME, BUT COLOR STACK BARS BY ORGANISM CATEGORY (METAGENOME OR ISOLATE)
# --------------------------------------------------------------------

records = [
    {"quarter_bin": q, "organism_type": t, "unique_accessions": len(accs)}
    for (q, t), accs in seen.items()
]

quarter_type = (
    pd.DataFrame(records)
      .pivot(index="quarter_bin", columns="organism_type", values="unique_accessions")
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

# --------------------------------------------------------------------
# 1 ONLY AMR-POSITIVE METAGENOME CATEGORIES - COLOR STACK BARS BY METAGENOME CATEGORY
# --------------------------------------------------------------------

# Debugging print
print("Number of entries in seen:", len(seen))
print("First 5 seen keys:", list(seen.keys())[:5])

records_cat = [
    {"quarter_bin": q, "metagenome_category": t, "unique_accessions": len(accs)}
    for (q, t), accs in seen.items()
]

quarter_cat = (
    pd.DataFrame(records_cat)
      .pivot(index="quarter_bin", columns="metagenome_category", values="unique_accessions")
      .fillna(0)
      .astype(int)
      .sort_index()
)

base_col = {
    "livestock":   "#fac723",
    "human":       "#936fac",
    "wastewater":  "#e95e50",
    "marine":      "#0cb2af",
    "freshwater":  "#a1c65d",
    "soil":        "#f29222",
    "other":       "#dedede",
}

category_colors = [base_col.get(cat, "#cccccc") for cat in quarter_cat.columns]

# Debugging print
print("Categories in quarter_cat:", quarter_cat.columns.tolist())
print("Assigned colors:", category_colors)

fig, ax = plt.subplots(figsize=(12, 4))
quarter_cat.plot(kind="bar", stacked=True, ax=ax, width=0.9, color=category_colors)

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
plt.close()

####################################
### No "other" category in the plot
####################################
quarter_cat_no_other = quarter_cat.drop(columns="other", errors="ignore")
colors_no_other = [base_col.get(cat, "#cccccc") for cat in quarter_cat_no_other.columns]

fig, ax = plt.subplots(figsize=(12, 4))
quarter_cat_no_other.plot(
    kind="bar",
    stacked=True,
    ax=ax,
    width=0.9,
    color=colors_no_other,
)

ax.set_title("Discovery timeline of AMR-positive metagenome samples\n(without 'other' category)")
ax.set_xlabel("")
ax.set_ylabel("# SRA accessions AMR-positive")

# Year ticks
period_idx = quarter_cat_no_other.index
year_pos = [i for i, p in enumerate(period_idx) if p.quarter == 1]
year_labs = [str(p.year) for p in period_idx if p.quarter == 1]
ax.set_xticks(year_pos)
ax.set_xticklabels(year_labs, rotation=0)

ax.legend(title="Metagenome category", bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
fig.savefig("../data/discovery_timeline_amr_categories_no_other.png", dpi=600)
fig.savefig("../data/discovery_timeline_amr_categories_no_other.svg")
plt.close()

####################################
### Other at the bottom of the plot
####################################

# Reorder columns: 'other' first, then the rest
cols = quarter_cat.columns.tolist()
if "other" in cols:
    reordered_cols = ["other"] + [c for c in cols if c != "other"]
else:
    reordered_cols = cols

quarter_cat_reordered = quarter_cat[reordered_cols]
colors_reordered = [base_col.get(cat, "#cccccc") for cat in reordered_cols]

fig, ax = plt.subplots(figsize=(12, 4))
quarter_cat_reordered.plot(
    kind="bar",
    stacked=True,
    ax=ax,
    width=0.9,
    color=colors_reordered,
)

ax.set_title("Discovery timeline of AMR-positive metagenome samples\n('other' stacked at bottom)")
ax.set_xlabel("")
ax.set_ylabel("# SRA accessions AMR-positive")

# Year ticks
period_idx = quarter_cat_reordered.index
year_pos = [i for i, p in enumerate(period_idx) if p.quarter == 1]
year_labs = [str(p.year) for p in period_idx if p.quarter == 1]
ax.set_xticks(year_pos)
ax.set_xticklabels(year_labs, rotation=0)

ax.legend(title="Metagenome category", bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
fig.savefig("../data/discovery_timeline_amr_categories_other_bottom.png", dpi=600)
fig.savefig("../data/discovery_timeline_amr_categories_other_bottom.svg")
plt.close()
'''