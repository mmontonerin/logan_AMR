import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------
# file paths
# -------------------------------------------------
total_file = "../data/plasmids_sra_metadata_extended.csv" # all plasmids
amr_file   = "../data/amr_metadata_extended.csv" # AMR-positive subset
plot_png1   = "../data/amr_metagenome_enrichment_norm.png"
plot_svg1   = "../data/amr_metagenome_enrichment_norm.svg"
plot_png2   = "../data/amr_organism_enrichment_norm.png"
plot_svg2   = "../data/amr_organism_enrichment_norm.svg"

rng = np.random.default_rng(42)          # reproducible random seed

# -------------------------------------------------
# 1. load data
# -------------------------------------------------
df_total = pd.read_csv(total_file, dtype=str, low_memory=False)
df_amr   = pd.read_csv(amr_file,   dtype=str, low_memory=False)

# ────────────────────────────────────────────────────────────────
# PART A  ·  ISOLATE  vs  METAGENOME  (balanced)
# ────────────────────────────────────────────────────────────────
cats_org = ["Isolate", "Metagenome"]

# sets of unique plasmid ids (seq_name) per category
total_sets_org = {c: set(df_total.loc[df_total["organism_type"] == c, "seq_name"])
                  for c in cats_org}
amr_sets_org   = {c: set(df_amr  .loc[df_amr  ["organism_type"] == c, "seq_name"])
                  for c in cats_org}

# drop categories with zero records in df_total
total_sets_org = {c: s for c, s in total_sets_org.items() if len(s) > 0}
amr_sets_org   = {c: amr_sets_org[c] for c in total_sets_org.keys()}
cats_org       = list(total_sets_org.keys())

n_per_cat = min(len(s) for s in total_sets_org.values())

totals_bal, amrs_bal = {}, {}
for c in cats_org:
    sampled = rng.choice(list(total_sets_org[c]), n_per_cat, replace=False)
    sampled_set = set(sampled)
    totals_bal[c] = n_per_cat
    amrs_bal[c]   = len(amr_sets_org[c] & sampled_set)

series_tot = pd.Series(totals_bal)
series_amr = pd.Series(amrs_bal)

density_overall = series_amr.sum() / series_tot.sum()
enrichment_org  = np.log2((series_amr / series_tot) / density_overall)

# ---- plot isolate vs metagenome --------------------------------
palette_org = {
    "Isolate":    "#C2B2FF",
    "Metagenome": "#FDBF6F",
}

cats_plot = enrichment_org.sort_values().index
vals      = enrichment_org.loc[cats_plot]
colors    = [palette_org.get(c, "#bbbbbb") for c in cats_plot]

fig, ax = plt.subplots(figsize=(7, 4))
ax.barh(cats_plot, vals, color=colors, edgecolor="black", alpha=0.9)
ax.axvline(0, color="black", linewidth=1)

ax.set_xlabel("log₂ enrichment of AMR-positive plasmids\n"
              "(balanced isolate vs metagenome)", fontsize=11)
ax.set_ylabel("")
ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)

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

# ────────────────────────────────────────────────────────────────
# PART B  ·  METAGENOME CATEGORIES  (balanced)
# ────────────────────────────────────────────────────────────────
mask_amr   = df_amr["metagenome_category"].notna() & \
             (df_amr["metagenome_category"].str.lower() != "other")
mask_total = df_total["metagenome_category"].notna() & \
             (df_total["metagenome_category"].str.lower() != "other")

df_amr_meta   = df_amr[mask_amr]
df_total_meta = df_total[mask_total]

cats_meta = sorted(df_total_meta["metagenome_category"].unique())

total_sets_meta = {c: set(df_total_meta.loc[df_total_meta["metagenome_category"] == c, "seq_name"])
                   for c in cats_meta}
amr_sets_meta   = {c: set(df_amr_meta  .loc[df_amr_meta  ["metagenome_category"] == c, "seq_name"])
                   for c in cats_meta}

# drop empty categories
total_sets_meta = {c: s for c, s in total_sets_meta.items() if len(s) > 0}
amr_sets_meta   = {c: amr_sets_meta[c] for c in total_sets_meta.keys()}
cats_meta       = list(total_sets_meta.keys())

n_per_cat = min(len(s) for s in total_sets_meta.values())

totals_bal_meta, amrs_bal_meta = {}, {}
for c in cats_meta:
    sampled = rng.choice(list(total_sets_meta[c]), n_per_cat, replace=False)
    sampled_set = set(sampled)
    totals_bal_meta[c] = n_per_cat
    amrs_bal_meta[c]   = len(amr_sets_meta[c] & sampled_set)

series_tot_meta = pd.Series(totals_bal_meta)
series_amr_meta = pd.Series(amrs_bal_meta)

density_overall_meta = series_amr_meta.sum() / series_tot_meta.sum()
enrichment_meta      = np.log2((series_amr_meta / series_tot_meta) /
                               density_overall_meta)

# ---- plot metagenome categories --------------------------------
palette_meta = {
    "livestock":  "#fac723",
    "human":      "#936fac",
    "wastewater": "#e95e50",
    "marine":     "#0cb2af",
    "freshwater": "#a1c65d",
    "soil":       "#f29222",
}

cats_plot = enrichment_meta.sort_values().index
vals      = enrichment_meta.loc[cats_plot]
colors    = [palette_meta.get(c, "#bbbbbb") for c in cats_plot]

fig, ax = plt.subplots(figsize=(7, 4))
ax.barh(cats_plot, vals, color=colors, edgecolor="black", alpha=0.9)
ax.axvline(0, color="black", linewidth=1)

ax.set_xlabel("log₂ enrichment of AMR-positive plasmids\n"
              "(balanced metagenome categories)", fontsize=11)
ax.set_ylabel("")
ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)

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