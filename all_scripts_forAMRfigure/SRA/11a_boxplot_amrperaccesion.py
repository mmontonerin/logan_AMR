import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 
#csv_path = "../data/card_metadata_aro_extended.csv"
csv_path = "../data/card_metadata_aro_dateloc_meta.csv"

chunksize = 10_000_000
per_acc_chunks = []

# Read only needed columns
usecols = ["acc", "metagenome_category", "ARO_ID"]

base_col = {
    "livestock":   "#fac723",
    "human":       "#936fac",
    "wastewater":  "#e95e50",
    "marine":      "#0cb2af",
    "freshwater":  "#a1c65d",
    "soil":        "#f29222",
}

for chunk in pd.read_csv(csv_path, usecols=usecols, chunksize=chunksize, low_memory=False):

    wanted = set(base_col)
    chunk = chunk[chunk["metagenome_category"].isin(wanted)]

    # Count AMRs per accession per category in the chunk
    counts = (
        chunk.groupby(["metagenome_category", "acc"])
             .size()
             .rename("count")
             .reset_index()
    )
    per_acc_chunks.append(counts)

# Concatenate and group again to consolidate duplicates across chunks
per_acc = pd.concat(per_acc_chunks)
per_acc = (
    per_acc.groupby(["metagenome_category", "acc"])["count"]
           .sum()
           .reset_index()
)

# ---------- 2. build a list of count-arrays, one per category ----------
cat_order = (
    per_acc.groupby("metagenome_category")["count"]
           .mean()
           .sort_values()
           .index
)

data_for_plot = [per_acc.loc[per_acc["metagenome_category"] == cat, "count"] for cat in cat_order]

# ------------------------------------------------------------------
# 3A. STANDARD LINEAR-SCALE BOX PLOT
# ------------------------------------------------------------------
# fig, ax = plt.subplots(figsize=(8, 3 + 0.35 * len(cat_order)))

# Banner style plot
fig, ax = plt.subplots(figsize=(18, 5))

box = ax.boxplot(
    data_for_plot,
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
    patch.set_alpha(0.9)

for median in box["medians"]:
    median.set_color("black")
    median.set_linewidth(1.4)

ax.set_xlabel("Number of AMRs per accession")
ax.set_title("Per-accession count distribution")
ax.set_xlim(left=0)
ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.savefig("../data/boxplot_amr_per_accession_kept_linear_bannerstyle.png", dpi=600)
plt.savefig("../data/boxplot_amr_per_accession_kept_linear_bannerstyle.svg")
plt.close()

'''
# ------------------------------------------------------------------
# 3B.  LOG-SCALE VARIANT  (comment-out if you don't need it)
# ------------------------------------------------------------------
# Uncomment this block if you prefer everything on log10(count+1)

data_log = [np.log10(d + 1) for d in data_for_plot]

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
plt.savefig("../data/boxplot_amr_per_accession_kept_log2.png", dpi=600)
plt.savefig("../data/boxplot_amr_per_accession_kept_log2.svg")
plt.close()

# ------------------------------------------------------------------
# 3C.  BROKEN-AXIS VARIANT  (comment-out if you don't need it)
# ------------------------------------------------------------------
# Shows 0–600 AMRs in the top panel and the long tail in the bottom.

lims_main  = (0, 600)    # bulk of the data
lims_tail  = (600, max(map(max, data_for_plot))*1.05)  # outliers

fig, (ax1, ax2) = plt.subplots(
    nrows=2, sharey=True, figsize=(8, 4 + 0.4 * len(cat_order)),
    gridspec_kw=dict(height_ratios=[3, 1])
)
for ax, (xmin, xmax) in zip((ax1, ax2), (lims_main, lims_tail)):
    box = ax.boxplot(
        data_for_plot,
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
plt.savefig("../data/boxplot_amr_per_accession_kept_brokenaxis2.png", dpi=600)
plt.savefig("../data/boxplot_amr_per_accession_kept_brokenaxis2.svg")
plt.close()
'''