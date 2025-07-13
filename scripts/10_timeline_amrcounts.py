# -----------------------------------------------------------
# 0.  Imports & colour palette
# -----------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pathlib import Path
from scipy.stats import linregress   # for slope, intercept, rvalue

base_col = {
    "livestock":   "#fac723",
    "human":       "#936fac",
    "wastewater":  "#e95e50",
    "marine":      "#0cb2af",
    "freshwater":  "#a1c65d",
    "soil":        "#f29222",
}
wanted = set(base_col)               # only keep these categories

# -----------------------------------------------------------
# 1.  Load, clean, aggregate
# -----------------------------------------------------------
csv_path = Path("../data/full_card_metadata_aro_allfilters_metagenomes_WHOcategories_latlon.csv")
df = pd.read_csv(csv_path, low_memory=False)

# ---- 1 a.  keep rows whose category is one of ours ------------
df["metagenome_category_lc"] = df["metagenome_category"].str.lower()
df = df[df["metagenome_category_lc"].isin(wanted)]

# ---- 1 b.  parse the date & extract year ---------------------
# strip [] if present; coerce errors to NaT and drop missing
df["collection_date_sam"] = (
    df["collection_date_sam"]
      .str.strip("[]")
      .pipe(pd.to_datetime, errors="coerce", utc=False)
)
df = df.dropna(subset=["collection_date_sam"])

# ---- 1 c.  AMR rows per accession ----------------------------
per_acc = (
    df.groupby(["metagenome_category_lc", "acc"])
      .agg(
          date=("collection_date_sam", "min"),   # earliest date for that accession
          count=("acc", "size")                  # AMR rows per accession
      )
      .reset_index(drop=False)
)

per_acc = per_acc[per_acc["date"] >= pd.Timestamp("2010-01-01")]

def dt2num(series):
    """Convert pandas/NumPy datetime to the float format matplotlib uses."""
    return mdates.date2num(pd.to_datetime(series).to_numpy())

# -----------------------------------------------------------
# 2.  One panel per category  (scatter + regression)
# -----------------------------------------------------------
n_cat = len(wanted)
fig, axes = plt.subplots(
    nrows=n_cat, sharex=True, figsize=(8, 2.8 * n_cat),
    constrained_layout=True
)

# make sure plotting order is consistent with earlier scripts
cat_order = sorted(wanted)    # or choose your own order

for ax, cat in zip(axes, cat_order):
    sub = per_acc.loc[per_acc["metagenome_category_lc"] == cat]
    x_dates = pd.to_datetime(sub["date"])
    y_vals  = sub["count"].values

    # scatter on real dates
    ax.scatter(x_dates, y_vals, s=18, alpha=0.6,
               color=base_col[cat], edgecolors="none", label=cat.title())

    # regression in date-space: convert to float first
    x_num = dt2num(x_dates)
    slope, intercept, r, *_ = linregress(x_num, y_vals)

    # regression line (two end points are enough for a straight line)
    span = np.linspace(x_num.min(), x_num.max(), 2)
    ax.plot(
        mdates.num2date(span),
        intercept + slope * span,
        color=base_col[cat], linewidth=1.3, zorder=3,
        label=f"R = {r:.2f}"
    )

    # prettify axis
    ax.set_ylabel("AMRs / accession")
    ax.legend(loc="upper left", frameon=False)
    ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.4)

    # year ticks every 2 years (change interval if you have lot of years)
    ax.xaxis.set_major_locator(mdates.YearLocator(base=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

axes[-1].set_xlabel("Collection year")
fig.suptitle("Trend of AMRs per accession over time")
plt.savefig("../data/scatter_trend_per_category-2010min.png", dpi=600)
plt.close()

# -----------------------------------------------------------
# 3.  All categories in ONE plot
# -----------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)

for cat in cat_order:
    sub = per_acc.loc[per_acc["metagenome_category_lc"] == cat]
    x_dates = pd.to_datetime(sub["date"])
    y_vals  = sub["count"].values

    # scatter
    ax.scatter(x_dates, y_vals, s=15, alpha=0.35,
               color=base_col[cat], edgecolors="none")

    # regression
    x_num = dt2num(x_dates)
    slope, intercept, r, *_ = linregress(x_num, y_vals)
    span = np.linspace(x_num.min(), x_num.max(), 2)
    ax.plot(
        mdates.num2date(span),
        intercept + slope * span,
        color=base_col[cat], linewidth=2,
        label=f"{cat.title()}  (R={r:.2f})"
    )

# axis formatting
ax.set_xlabel("Collection year")
ax.set_ylabel("AMRs / accession")
ax.set_title("AMR trend across all metagenome categories")
ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.4)
ax.xaxis.set_major_locator(mdates.YearLocator(base=2))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax.legend(frameon=False, fontsize=9)
plt.savefig("../data/scatter_trend_all_categories-2010min.png", dpi=600)
