import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ---------- 1. read the data ----------
csv_path = Path("../data/full_card_metadata_aro_allfilters_metagenomes_WHOcategories_latlon.csv")
df = pd.read_csv(csv_path, low_memory=False)

# ---------- 2. counts per accession inside each category ----------
per_acc = (
    df.groupby(["metagenome_category", "acc"])
      .size()
      .rename("count")
      .reset_index()
)

# ---------- 3. build a list of count-arrays, one per category ----------
# sort by mean so the y-axis goes from “smallest typical spread” to largest
cat_order = (
    per_acc.groupby("metagenome_category")["count"]
           .mean()
           .sort_values()
           .index
)

data_for_plot = [per_acc.loc[per_acc["metagenome_category"] == cat, "count"] for cat in cat_order]

base_col = {
    "livestock":   "#fac723",
    "human":       "#936fac",
    "wastewater":  "#e95e50",
    "marine":      "#0cb2af",
    "freshwater":  "#a1c65d",
    "soil":        "#f29222",
}

# ---------- 4. plot ----------
fig, ax = plt.subplots(figsize=(8, 4 + 0.4*len(cat_order)))

box = ax.boxplot(
    data_for_plot,
    vert=False,
    patch_artist=True,     # lets us colour the boxes
    tick_labels=[c.title() for c in cat_order],
    showmeans=True,
    meanprops=dict(marker="|", markeredgecolor="black", markersize=14, markeredgewidth=2),
)

# colour each box to match your palette
for patch, cat in zip(box["boxes"], cat_order):
    patch.set_facecolor(base_col.get(cat, "#888888"))
    patch.set_alpha(0.6)

# median line is drawn automatically; tweak style if you like
for median in box["medians"]:
    median.set_color("black")
    median.set_linewidth(1.4)

ax.set_xlabel("Number of AMRs per accession")
ax.set_title("Per-accession count distribution within metagenome categories")
ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)

plt.tight_layout()
plt.savefig("../data/boxplot_amr_per_accession.png", dpi=600)
plt.savefig("../data/boxplot_amr_per_accession.svg")
