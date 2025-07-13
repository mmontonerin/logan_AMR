# Doughnut chart (ring) showing the AMR-positive proportion
from matplotlib import pyplot as plt

# 1. Load the tables
import pandas as pd
df_total = pd.read_csv("../data/plasmids_sra_metadata_extended.csv", dtype=str)
df_amr = pd.read_csv("../data/amr_metadata_extended.csv", dtype=str) 

# 2. Build unique-ID sets
df_total_set    = set(df_total["seq_name"].dropna().unique())
df_amr_set = set(df_amr["seq_name"].dropna().unique())

# 3. Quick sanity check / counts
total_seq     = len(df_total_set)
total_amr  = len(df_amr_set)
both_overlap  = len(df_total_set  & df_amr_set)
print(f"Plasmids total: {total_seq}")
print(f"AMR-positive: {total_amr}")
print(f"In both: {both_overlap}")

# --- data for the two slices -------------------------------------------------
amr_pos      = both_overlap
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
    startangle= 0, # coloured slice at 3 o'clock
    counterclock = False,
    wedgeprops = dict(width=0.3, edgecolor="w")  # width<1 turns it into a ring
)

plt.tight_layout()
plt.savefig("../data/ring_plasmids_amr.png", dpi=600)
plt.savefig("../data/ring_plasmids_amr.svg")
plt.close()

###################################################################################

isolate = df_amr.loc[df_amr["organism_type"] == "Isolate", "seq_name"].nunique()
metagenome  = df_amr.loc[df_amr["organism_type"] == "Metagenome", "seq_name"].nunique()

# --- data for the two slices -------------------------------------------------
sizes  = [metagenome, isolate]
labels = [f"Metagenome ({metagenome})", f"Isolate ({isolate})"]

colors = ["#FDBF6F", "#C2B2FF"]   # soft orange metagenome, soft purple isolate

fig, ax = plt.subplots(figsize=(4,4), subplot_kw=dict(aspect="equal"))

ax.pie(
    sizes,
    labels     = labels,
    colors     = colors,
    startangle = 75, # Metagenome towards the right
    counterclock = False,
    wedgeprops = dict(width=0.3, edgecolor="w")  # width<1 turns it into a ring
)

plt.tight_layout()
plt.savefig("../data/ring_plasmids_amr_isometa.png", dpi=600)
plt.savefig("../data/ring_plasmids_amr_isometa.svg")
plt.close()

###################################################################
meta_only = df_amr[df_amr["organism_type"] == "Metagenome"]

base_col = {
    "livestock":   "#fac723",
    "human":       "#936fac",
    "wastewater":  "#e95e50",
    "marine":      "#0cb2af",
    "freshwater":  "#a1c65d",
    "soil":        "#f29222",
    "other":       "#dedede",
}
# Data lists
sizes   = []
labels  = []
colors  = []

for cat, col in base_col.items():
    n = meta_only.loc[meta_only["metagenome_category"] == cat, "seq_name"].nunique()
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