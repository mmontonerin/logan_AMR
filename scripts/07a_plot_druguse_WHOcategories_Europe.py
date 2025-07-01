import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# --------------------------------------------------------------
# paths
# --------------------------------------------------------------
csv_file = "../data/Antibiotic_Use_AWaRe.csv"
out_dir  = "../data/aware_use_europe"
os.makedirs(out_dir, exist_ok=True)

# --------------------------------------------------------------
# load & filter
# --------------------------------------------------------------
use_df = (
    pd.read_csv(csv_file, low_memory=False)
      .query("WHORegionCode == 'EUR'")
      .loc[lambda d: d["YEAR"].between(2016, 2022)]
)

if use_df.empty:
    print("No EUR rows between 2016-2022 - nothing to plot.")
    quit()

# --------------------------------------------------------------
# aggregate: DID per (year, AWaReLabel)
# --------------------------------------------------------------
agg = (
    use_df.groupby(["YEAR", "AWaReLabel"], as_index=False)["DID"]
          .sum()
)

# full x-axis 2016–2022 so every plot is same width
years = np.arange(2016, 2023)

# base colours (same hues you used before)
base_colour = {
    "Access":         "#4b8e0d",
    "Watch":          "#d95f02",
    "Reserve":        "#c51b1b",
    "Not classified": "#999999"
}

for cat, colour in base_colour.items():
    # y values (sum DID) – fill missing years with 0
    y = (
        agg.loc[agg["AWaReLabel"] == cat]
           .set_index("YEAR")["DID"]
           .reindex(years, fill_value=0)
           .sort_index()
    )

    fig, ax = plt.subplots(figsize=(10, 3))
    ax.plot(years, y.values, color=colour, linewidth=3)
    ax.set_xlim(years[0], years[-1])
    ax.set_xticks(years)
    ax.set_xlabel("Year")
    ax.set_ylabel("Defined Daily Dose per 1000 inhabitants")
    ax.set_title(f"{cat} antibiotic usage Europe")
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", alpha=0.3, linestyle="--")

    fname = os.path.join(out_dir, f"aware_{cat.lower().replace(' ', '_')}.png")
    fig.tight_layout()
    fig.savefig(fname, dpi=600)
    fig.savefig(fname.replace(".png", ".svg"))
    plt.close()

print("plots written to:", os.path.abspath(out_dir))
