import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

csv_path = Path("../data/full_card_metadata_aro_allfilters_metagenomes_WHOcategories_latlon.csv")
df = pd.read_csv(csv_path, low_memory=False)

# ------------------------------------------------------------------
# 1. keep only the categories we care about
# ------------------------------------------------------------------
base_col = {
    "livestock":   "#fac723",
    "human":       "#936fac",
    "wastewater":  "#e95e50",
    "marine":      "#0cb2af",
    "freshwater":  "#a1c65d",
    "soil":        "#f29222",
}

wanted = set(base_col)                      # {'livestock', ...}

# make a lower-case helper column so we can match case-insensitively
df["metagenome_category_lc"] = df["metagenome_category"].str.lower()
df = df[df["metagenome_category_lc"].isin(wanted)]

# ------------------------------------------------------------------
# 2. one row-count per accession, *inside* each kept category
# ------------------------------------------------------------------
per_acc = (
    df.groupby(["metagenome_category_lc", "acc"])
      .size()
      .rename("count")
      .reset_index()
)

# you can keep a fixed order or sort by mean; here we sort by mean
cat_order = (
    per_acc.groupby("metagenome_category_lc")["count"]
           .mean()
           .sort_values()
           .index
)

# build the raw (linear) data series
data_linear = [
    per_acc.loc[per_acc["metagenome_category_lc"] == cat, "count"]
    for cat in cat_order
]

# ------------------------------------------------------------------
# 3A. STANDARD LINEAR-SCALE BOX PLOT
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 3 + 0.35 * len(cat_order)))

box = ax.boxplot(
    data_linear,
    vert=False,
    patch_artist=True,
    showmeans=True,
    meanprops=dict(marker="|",
                   markeredgecolor="black",
                   markersize=14,
                   markeredgewidth=2),
    flierprops=dict(marker="o",
                    markersize=3,
                    markerfacecolor="black",
                    alpha=0.25),
    tick_labels=[c.title() for c in cat_order],
)

for patch, cat in zip(box["boxes"], cat_order):
    patch.set_facecolor(base_col[cat])
    patch.set_alpha(0.6)

for median in box["medians"]:
    median.set_color("black")
    median.set_linewidth(1.4)

ax.set_xlabel("Number of AMRs per accession")
ax.set_title("Per-accession count distribution")
ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.savefig("../data/boxplot_amr_per_accession_kept_linear.png", dpi=600)

# ------------------------------------------------------------------
# 3B.  LOG-SCALE VARIANT  (comment-out if you don't need it)
# ------------------------------------------------------------------
# Uncomment this block if you prefer everything on log10(count+1)

data_log = [np.log10(d + 1) for d in data_linear]

fig, ax = plt.subplots(figsize=(8, 3 + 0.35 * len(cat_order)))
box = ax.boxplot(
    data_log,
    vert=False,
    patch_artist=True,
    showmeans=True,
    meanprops=dict(marker="|", markeredgecolor="black", markersize=14, markeredgewidth=2),
    flierprops=dict(marker="o", markersize=3, markerfacecolor="black", alpha=0.25),
    tick_labels=[c.title() for c in cat_order],
)
for patch, cat in zip(box["boxes"], cat_order):
    patch.set_facecolor(base_col[cat])
    patch.set_alpha(0.6)
for median in box["medians"]:
    median.set_color("black")
    median.set_linewidth(1.4)
ax.set_xlabel("log\u2081\u2080(AMRs per accession + 1)")
ax.set_title("Per-accession counts (log scale)")
ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.savefig("../data/boxplot_amr_per_accession_kept_log.png", dpi=600)


# ------------------------------------------------------------------
# 3C.  BROKEN-AXIS VARIANT  (comment-out if you don't need it)
# ------------------------------------------------------------------
# Shows 0–600 AMRs in the top panel and the long tail in the bottom.

lims_main  = (0, 600)    # bulk of the data
lims_tail  = (600, max(map(max, data_linear))*1.05)  # outliers

fig, (ax1, ax2) = plt.subplots(
    nrows=2, sharey=True, figsize=(8, 4 + 0.4 * len(cat_order)),
    gridspec_kw=dict(height_ratios=[3, 1])
)
for ax, (xmin, xmax) in zip((ax1, ax2), (lims_main, lims_tail)):
    box = ax.boxplot(
        data_linear,
        vert=False,
        patch_artist=True,
        showmeans=True,
        meanprops=dict(marker="|", markeredgecolor="black", markersize=14, markeredgewidth=2),
        flierprops=dict(marker="o", markersize=3, markerfacecolor="black", alpha=0.25),
        tick_labels=[c.title() for c in cat_order] if ax is ax1 else [""]*len(cat_order),
    )
    ax.set_xlim(xmin, xmax)
    ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)
    for patch, cat in zip(box["boxes"], cat_order):
        patch.set_facecolor(base_col[cat])
        patch.set_alpha(0.6)
    for median in box["medians"]:
        median.set_color("black")
        median.set_linewidth(1.4)

# small diagonal lines to show the axis break
d = .015
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)
kwargs.update(transform=ax2.transAxes)
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

ax1.set_title("Per-accession counts")
ax2.set_xlabel("Number of AMRs per accession")
plt.tight_layout()
plt.savefig("../data/boxplot_amr_per_accession_kept_brokenaxis.png", dpi=600)

