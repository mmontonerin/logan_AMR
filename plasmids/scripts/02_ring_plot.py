# Doughnut chart (ring) showing the AMR-positive proportion
from matplotlib import pyplot as plt

# 1. Load the tables
import pandas as pd
csv_df = pd.read_csv("../data/plasmids_sra_metadata.csv", dtype=str)       # has column "seq_name"
tsv_df = pd.read_csv("../data/argnorm_output.tsv", sep="\t", dtype=str)  # has column "Contig id"

# 2. Build unique-ID sets
seq_set    = set(csv_df["seq_name"].dropna().unique())
contig_set = set(tsv_df["Contig id"].dropna().unique())

# 3. Quick sanity check / counts
total_seq     = len(seq_set)
total_contig  = len(contig_set)
both_overlap  = len(seq_set & contig_set)
print(f"Plasmids total: {total_seq}")
print(f"AMR-positive: {total_contig}")
print(f"In both: {both_overlap}")

# --- data for the two slices -------------------------------------------------
amr_pos      = both_overlap          # or total_contig (same value here)
amr_negative = total_seq - amr_pos

sizes  = [amr_pos, amr_negative]
labels = [f"AMR-positive ({amr_pos})", f"Other plasmids ({amr_negative})"]

# choose your own palette; the second colour should be neutral / lighter
colors = ["#d95f02",  "#cfd8dc"]     # orange highlight, light grey remainder

fig, ax = plt.subplots(figsize=(4,4), subplot_kw=dict(aspect="equal"))

ax.pie(
    sizes,
    labels     = labels,
    colors     = colors,
    startangle = 90,                 # puts the coloured slice at the top
    counterclock = False,
    wedgeprops = dict(width=0.3, edgecolor="w")  # width<1 turns it into a ring
)

plt.tight_layout()
plt.savefig("../data/ring_plasmids_amr.png", dpi=600)
plt.savefig("../data/ring_plasmids_amr.svg")
plt.close()

###################################################################################

# Function to classify sample into Metagenome or Isolate
def classify_sample(df: pd.DataFrame) -> pd.DataFrame:
    # Add 'sample_type' = 'Metagenome' or 'Isolate'
    is_meta = df["organism"].str.contains("metagenome", case=False, na=False)
    df = df.copy()
    df["sample_type"] = is_meta.map({True: "Metagenome", False: "Isolate"})
    return df

amr_positive_csv = pd.read_csv("../data/amr_metadata_full.csv", dtype=str)

classified_df = classify_sample(amr_positive_csv)

contig_class = set(classified_df["seq_name"].dropna().unique())
contig_class = set(classified_df["sample_type"].dropna())
isolate_class = len(classified_df[classified_df["sample_type"] == "Isolate"])
metagenome_class = len(classified_df[classified_df["sample_type"] == "Metagenome"])

# --- data for the two slices -------------------------------------------------
isolate = isolate_class
metagenome = metagenome_class

sizes  = [isolate, metagenome]
labels = [f"Isolate ({isolate})", f"Metagenome ({metagenome})"]

colors = ["#C2B2FF", "#FDBF6F"]   # soft orange metagenome, soft purple isolate

fig, ax = plt.subplots(figsize=(4,4), subplot_kw=dict(aspect="equal"))

ax.pie(
    sizes,
    labels     = labels,
    colors     = colors,
    startangle = 90,                 # puts the coloured slice at the top
    counterclock = False,
    wedgeprops = dict(width=0.3, edgecolor="w")  # width<1 turns it into a ring
)

plt.tight_layout()
plt.savefig("../data/ring_plasmids_amr_isometa.png", dpi=600)
plt.savefig("../data/ring_plasmids_amr_isometa.svg")
plt.close()

###################################################################
# Metagenome category ring plot
amr_positive_meta = pd.read_csv("../data/amr_metadata_metagenomes.csv", dtype=str)

base_col = {
    "livestock":   "#fac723",
    "human":       "#936fac",
    "wastewater":  "#e95e50",
    "marine":      "#0cb2af",
    "freshwater":  "#a1c65d",
    "soil":        "#f29222",
    "Other":       "#dedede",
}
# Data lists
sizes   = []
labels  = []
colors  = []

for cat, col in base_col.items():
    n = (amr_positive_meta["metagenome_category"] == cat).sum()
    sizes.append(n)
    labels.append(f"{cat} ({n})")
    colors.append(col)

# 2. Create the ring ----------------------------------------------------------
fig, ax = plt.subplots(figsize=(5,5), subplot_kw=dict(aspect="equal"))

ax.pie(
    sizes,
    labels     = labels,
    colors     = colors,
    startangle = 90,
    counterclock = False,
    wedgeprops = dict(width=0.3, edgecolor="w")   # width<1 makes it a ring
)

total = sum(sizes)

plt.tight_layout()
plt.savefig("../data/ring_metagenome_categories.png", dpi=600)
plt.savefig("../data/ring_metagenome_categories.svg")
plt.close()