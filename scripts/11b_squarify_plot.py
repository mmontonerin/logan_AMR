from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import squarify

# ────────────────────────────────────────────────────────────────
# 1. file locations
# ────────────────────────────────────────────────────────────────
total_csv = Path("../data/SRA_metadata_before20231211_logan_extended.csv")
amr_csv   = Path("../data/card_metadata_aro_informativecolumns.csv")

total_df = pd.read_csv(total_csv, dtype=str, low_memory=False)
amr_df   = pd.read_csv(amr_csv, dtype=str, low_memory=False)

# ────────────────────────────────────────────────────────────────
# 4. final counts
# ────────────────────────────────────────────────────────────────

total_set = set(total_df["acc"].dropna().unique())
amr_set   = set(amr_df["acc"].dropna().unique())

amr_positive = len(total_set & amr_set)       # overlap only
amr_negative = len(total_set) - amr_positive  # all other

# counts within the AMR table (unique accs again)
isolate_n    = amr_df.loc[amr_df["organism_type"] == "Isolate", "acc"].nunique()
metagenome_n = amr_df.loc[amr_df["organism_type"] == "Metagenome", "acc"].nunique()

meta_only = amr_df[amr_df["organism_type"] == "Metagenome"]

# ────────────────────────────────────────────────────────────────
# 5. helper: draw and save one treemap
# ────────────────────────────────────────────────────────────────
def save_treemap(sizes, labels, colors, outstem, figsize=(4, 4)):
    plt.figure(figsize=figsize)

    # squarify prefers descending areas for a nicer layout
    order = sorted(range(len(sizes)), key=lambda i: sizes[i], reverse=True)
    ordered_sizes  = [sizes[i]  for i in order]
    ordered_labels = [labels[i] for i in order]
    ordered_colors = [colors[i] for i in order]

    squarify.plot(
        sizes=ordered_sizes,
        label=ordered_labels,
        color=ordered_colors,
        text_kwargs=dict(size=10, weight="bold"),
        pad=True
    )
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(f"{outstem}.png", dpi=600)
    plt.savefig(f"{outstem}.svg")
    plt.close()

# ────────────────────────────────────────────────────────────────
# 6. treemap 1 – amr-positive proportion
# ────────────────────────────────────────────────────────────────
sizes  = [amr_positive, amr_negative]
labels = [f"AMR-positive ({amr_positive})",
          f"Other plasmids ({amr_negative})"]
colors = ["#d95f02", "#cfd8dc"]

save_treemap(sizes, labels, colors,
             "../data/treemap_sra_amr")

# ────────────────────────────────────────────────────────────────
# 7. treemap 2 – metagenome vs isolate
# ────────────────────────────────────────────────────────────────
sizes  = [metagenome_n, isolate_n]
labels = [f"Metagenome ({metagenome_n})",
          f"Isolate ({isolate_n})"]
colors = ["#FDBF6F", "#C2B2FF"]

save_treemap(sizes, labels, colors,
             "../data/treemap_sra_amr_isometa")

# ────────────────────────────────────────────────────────────────
# 8. treemap 3 – environment categories (metagenome only)
# ────────────────────────────────────────────────────────────────
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
    n = meta_only.loc[meta_only["metagenome_category"] == cat, "acc"].nunique()
    sizes.append(n)
    labels.append(f"{cat} ({n})")
    colors.append(col)

save_treemap(sizes, labels, colors,
             "../data/treemap_metagenome_categories",
             figsize=(5, 5))