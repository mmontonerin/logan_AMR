import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------
# file paths
# -------------------------------------------------
total_file = "../data/plasmids_sra_metadata_extended.csv" # all plasmids
amr_file   = "../data/amr_metadata_extended.csv" # AMR-positive subset
plot_png1   = "../data/amr_metagenome_enrichment.png"
plot_svg1   = "../data/amr_metagenome_enrichment.svg"
plot_png2   = "../data/amr_organism_enrichment.png"
plot_svg2   = "../data/amr_organism_enrichment.svg"

# -------------------------------------------------
# 1. load data
# -------------------------------------------------
df_total = pd.read_csv(total_file, dtype=str, low_memory=False)
df_amr   = pd.read_csv(amr_file,   dtype=str, low_memory=False)

# -------------------------------------------------
# 3. unique-plasmid counts per organism type
# -------------------------------------------------
totals = (df_total
          .groupby("organism_type")["seq_name"]
          .nunique())

amrs = (df_amr
        .groupby("organism_type")["seq_name"]
        .nunique())

# align indices
cats = sorted(set(totals.index) | set(amrs.index))
totals = totals.reindex(cats, fill_value=0)
amrs   = amrs .reindex(cats, fill_value=0)

# -------------------------------------------------
# 4. log₂ enrichment
# -------------------------------------------------
density_overall = amrs.sum() / totals.sum()
enrichment = np.log2((amrs / totals.replace(0, np.nan)) / density_overall).dropna()

# -------------------------------------------------
# 5. bar plot (colours reused from earlier ring plot)
# -------------------------------------------------
palette = {
    "Isolate":  "#C2B2FF",
    "Metagenome": "#FDBF6F",
} 

cats_plot = enrichment.sort_values().index
vals      = enrichment.loc[cats_plot]
colors    = [palette.get(cat, "#bbbbbb") for cat in cats_plot]

fig, ax = plt.subplots(figsize=(7, 4))
ax.barh(cats_plot, vals, color=colors, edgecolor="black", alpha=0.9)
ax.axvline(0, color="black", linewidth=1)

ax.set_xlabel("log₂ enrichment of AMR-positive plasmids\nvs. mean across organism categories",
              fontsize=11)
ax.set_ylabel("")
ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)

# directional labels
xmax, xmin = vals.max(), vals.min()
pad = (xmax - xmin) * 0.02
ax.text(xmin - pad, len(cats_plot) + 0.3,
        "enriched in AMR-negative plasmids ←", ha="left", va="center")
ax.text(xmax + pad, len(cats_plot) + 0.3,
        "→ enriched in AMR-positive plasmids", ha="right", va="center")

plt.tight_layout()
plt.savefig(plot_png2, dpi=600)
plt.savefig(plot_svg2)
plt.close()

################################
# Metagenome enrichment analysis

# -------------------------------------------------
# 2. drop the uninformative 'other' category
# -------------------------------------------------
mask_amr   = df_amr["metagenome_category"].notna() & \
             (df_amr["metagenome_category"].str.lower() != "other")
mask_total = df_total["metagenome_category"].notna() & \
             (df_total["metagenome_category"].str.lower() != "other")

df_amr   = df_amr[mask_amr]
df_total = df_total[mask_total]

# -------------------------------------------------
# 3. unique-plasmid counts per biome
# -------------------------------------------------
totals = (df_total
          .groupby("metagenome_category")["seq_name"]
          .nunique())

amrs = (df_amr
        .groupby("metagenome_category")["seq_name"]
        .nunique())

# align indices
cats = sorted(set(totals.index) | set(amrs.index))
totals = totals.reindex(cats, fill_value=0)
amrs   = amrs .reindex(cats, fill_value=0)

# -------------------------------------------------
# 4. log₂ enrichment
# -------------------------------------------------
density_overall = amrs.sum() / totals.sum()
enrichment = np.log2((amrs / totals.replace(0, np.nan)) / density_overall).dropna()

# -------------------------------------------------
# 5. bar plot (colours reused from earlier ring plot)
# -------------------------------------------------
palette = {
    "livestock":  "#fac723",
    "human":      "#936fac",
    "wastewater": "#e95e50",
    "marine":     "#0cb2af",
    "freshwater": "#a1c65d",
    "soil":       "#f29222",
}

cats_plot = enrichment.sort_values().index
vals      = enrichment.loc[cats_plot]
colors    = [palette.get(cat, "#bbbbbb") for cat in cats_plot]

fig, ax = plt.subplots(figsize=(7, 4))
ax.barh(cats_plot, vals, color=colors, edgecolor="black", alpha=0.9)
ax.axvline(0, color="black", linewidth=1)

ax.set_xlabel("log₂ enrichment of AMR-positive plasmids\nvs. mean across metagenome categories",
              fontsize=11)
ax.set_ylabel("")
ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)

# directional labels
xmax, xmin = vals.max(), vals.min()
pad = (xmax - xmin) * 0.02
ax.text(xmin - pad, len(cats_plot) + 0.3,
        "enriched in AMR-negative plasmids ←", ha="left", va="center")
ax.text(xmax + pad, len(cats_plot) + 0.3,
        "→ enriched in AMR-positive plasmids", ha="right", va="center")

plt.tight_layout()
plt.savefig(plot_png1, dpi=600)
plt.savefig(plot_svg1)
plt.close()