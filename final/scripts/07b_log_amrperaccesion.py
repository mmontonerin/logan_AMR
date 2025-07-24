import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------
# file paths
# -------------------------------------------------
#csv_path = "../data/card_metadata_aro_dateloc_meta.csv"
csv_path = "../data/card_metadata_aro_extended.csv"
plot_png = "../data/amr_density_log2enrichment_extendeddata.png"
plot_svg = "../data/amr_density_log2enrichment_extendeddata.svg"

# -------------------------------------------------
# load & preprocess: per-sample AMR counts
# -------------------------------------------------

chunksize = 10_000_000
usecols = ["acc", "metagenome_category", "ARO_ID"]
wanted_categories = ["livestock", "human", "wastewater", "marine", "freshwater", "soil"]

base_col = {
    "livestock":  "#fac723",
    "human":      "#936fac",
    "wastewater": "#e95e50",
    "marine":     "#0cb2af",
    "freshwater": "#a1c65d",
    "soil":       "#f29222",
}

per_acc_chunks = []

for chunk in pd.read_csv(csv_path, usecols=usecols, chunksize=chunksize, low_memory=False):
    chunk = chunk[chunk["metagenome_category"].isin(wanted_categories)]
    counts = (
        chunk.groupby(["metagenome_category", "acc"])
             .size()
             .rename("count")
             .reset_index()
    )
    per_acc_chunks.append(counts)

# Combine all chunks
per_acc = pd.concat(per_acc_chunks)
per_acc = (
    per_acc.groupby(["metagenome_category", "acc"])["count"]
           .sum()
           .reset_index()
)

# -------------------------------------------------
# balanced random sampling and AMR-per-sample mean
# -------------------------------------------------

# 1. Create dictionary of accs per category
groups = per_acc.groupby("metagenome_category")["acc"].unique().to_dict()

# 2. Create lookup: acc → count
acc_to_count = per_acc.set_index("acc")["count"].to_dict()

# 3. Balanced sampling
rng = np.random.default_rng(42)
n_per_cat = min(len(accs) for accs in groups.values())
n_bootstraps = 100

boot_means = {cat: [] for cat in groups}

for _ in range(n_bootstraps):
    for cat, accs in groups.items():
        sampled = rng.choice(accs, size=n_per_cat, replace=False)
        mean_amr = np.mean([acc_to_count[acc] for acc in sampled])
        boot_means[cat].append(mean_amr)

# 4. Mean AMR count per category
mean_per_cat = {cat: np.mean(vals) for cat, vals in boot_means.items()}
overall_mean = np.mean(list(mean_per_cat.values()))

# 5. log₂ enrichment vs. overall mean
log2_enrichment = {
    cat: np.log2(mean / overall_mean)
    for cat, mean in mean_per_cat.items()
}

# -------------------------------------------------
# bar plot of enrichment
# -------------------------------------------------

cats_sorted = sorted(log2_enrichment, key=log2_enrichment.get)
vals = [log2_enrichment[cat] for cat in cats_sorted]
colors = [base_col.get(cat, "#bbbbbb") for cat in cats_sorted]

fig, ax = plt.subplots(figsize=(4, 4))
ax.barh(cats_sorted, vals, color=colors, edgecolor="black", alpha=0.9)
ax.axvline(0, color="black", linewidth=1)

ax.set_xlabel("log₂ enrichment of AMR genes per accession\n"
              "vs. average across metagenome categories", fontsize=11)
ax.set_ylabel("")
ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)

xmax, xmin = max(vals), min(vals)
pad = (xmax - xmin) * 0.02
ax.text(xmin - pad, len(cats_sorted) + 0.3,
        "← fewer AMR genes per accession", ha="left", va="center")
ax.text(xmax + pad, len(cats_sorted) + 0.3,
        "more AMR genes per accession →", ha="right", va="center")

plt.tight_layout()
plt.savefig(plot_png, dpi=600)
plt.savefig(plot_svg)
plt.close()
