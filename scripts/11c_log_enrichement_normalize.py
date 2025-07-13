import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------
# file paths
# -------------------------------------------------
total_file = "../data/SRA_metadata_before20231211_logan_extended.csv" # all SRA
amr_file   = "../data/card_metadata_aro_informativecolumns.csv" # AMR-positive subset
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
# -------------------------------------------------
# ISOLATE / METAGENOME  –  balanced sample
# -------------------------------------------------

# 1. unique sets per organism_type
cats = ["Isolate", "Metagenome"]

total_set = {
    cat: set(df_total.loc[df_total["organism_type"] == cat, "acc"])
    for cat in cats
}
amr_set = {
    cat: set(df_amr.loc[df_amr["organism_type"] == cat, "acc"])
    for cat in cats
}

# 2. smallest bucket size  → sample size
n_per_cat = min(len(s) for s in total_set.values())
rng = np.random.default_rng(42)          # reproducible shuffling

totals_sub, amrs_sub = {}, {}
for cat in cats:
    sampled = rng.choice(list(total_set[cat]), n_per_cat, replace=False)
    sampled_set = set(sampled)
    totals_sub[cat] = n_per_cat
    amrs_sub[cat]   = len(amr_set[cat] & sampled_set)

totals = pd.Series(totals_sub)
amrs   = pd.Series(amrs_sub)

# 3. log₂ enrichment
density_overall = amrs.sum() / totals.sum()
enrichment = np.log2((amrs / totals) / density_overall).dropna()

# 4. bar plot
palette = {
    "Isolate":    "#C2B2FF",
    "Metagenome": "#FDBF6F",
}

cats_plot = enrichment.sort_values().index
vals      = enrichment.loc[cats_plot]
colors    = [palette.get(cat, "#bbbbbb") for cat in cats_plot]

fig, ax = plt.subplots(figsize=(7, 4))
ax.barh(cats_plot, vals, color=colors, edgecolor="black", alpha=0.9)
ax.axvline(0, color="black", linewidth=1)

ax.set_xlabel("log₂ enrichment of AMR-positive SRA accessions\n"
              "vs. mean across balanced organism categories",
              fontsize=11)
ax.set_ylabel("")
ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)

xmax, xmin = vals.max(), vals.min()
pad = (xmax - xmin) * 0.02
ax.text(xmin - pad, len(cats_plot) + 0.3,
        "enriched in AMR-negative SRA accessions ←", ha="left", va="center")
ax.text(xmax + pad, len(cats_plot) + 0.3,
        "→ enriched in AMR-positive SRA accessions", ha="right", va="center")

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
# 3. build unique-acc sets per biome
# -------------------------------------------------
cats = sorted(df_total["metagenome_category"].unique())

total_set = {
    cat: set(df_total.loc[df_total["metagenome_category"] == cat, "acc"])
    for cat in cats
}
amr_set = {
    cat: set(df_amr.loc[df_amr["metagenome_category"] == cat, "acc"])
    for cat in cats
}

# -------------------------------------------------
# 4. balance sample sizes
# -------------------------------------------------
n_per_cat = min(len(s) for s in total_set.values())     # smallest bucket

rng = np.random.default_rng(42)                         # reproducible
totals_sub, amrs_sub = {}, {}

for cat in cats:
    sampled = rng.choice(list(total_set[cat]), n_per_cat, replace=False)
    sampled_set = set(sampled)
    totals_sub[cat] = n_per_cat
    amrs_sub[cat]   = len(amr_set[cat] & sampled_set)

totals = pd.Series(totals_sub)
amrs   = pd.Series(amrs_sub)

# -------------------------------------------------
# 5. log₂ enrichment
# -------------------------------------------------
density_overall = amrs.sum() / totals.sum()
enrichment = np.log2((amrs / totals) / density_overall).dropna()

# -------------------------------------------------
# 6. bar plot (same palette and styling)
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

ax.set_xlabel("log₂ enrichment of AMR-positive SRA accessions\n"
              "vs. mean across balanced metagenome categories",
              fontsize=11)
ax.set_ylabel("")
ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)

xmax, xmin = vals.max(), vals.min()
pad = (xmax - xmin) * 0.02
ax.text(xmin - pad, len(cats_plot) + 0.3,
        "enriched in AMR-negative SRA accessions ←", ha="left", va="center")
ax.text(xmax + pad, len(cats_plot) + 0.3,
        "→ enriched in AMR-positive SRA accessions", ha="right", va="center")

plt.tight_layout()
plt.savefig(plot_png1, dpi=600)
plt.savefig(plot_svg1)
plt.close()