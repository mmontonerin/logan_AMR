import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.proportion import proportions_ztest

"""
Slope-graph of AMR prevalence:
2000-2011  vs  2012-2023
"""

# ---------------------------------------------------------------#
# 0 Prepare data
# ---------------------------------------------------------------#


sns.set_theme(style="ticks")

HITS_CSV   = "../data/full_card_metadata_aro_allfilters_metagenomes2.csv"
META_CSV   = "../data/SRA_metadata_allfilters.csv"
OUT_PNG    = "../data/slope_AMR_2000-11_vs_2012-23.png"
OUT_SVG    = "../data/slope_AMR_2000-11_vs_2012-23.svg"

DRUG_SET = {
    "glycopeptide antibiotic", "tetracycline antibiotic",
    "fluoroquinolone antibiotic", "aminoglycoside antibiotic",
    "penicillin beta-lactam", "glycylcycline", "cephalosporin"
}

SAMPLE_ID = "acc" # SRA sample ID
CONT_COL  = "geo_loc_name_country_continent_calc"
CAT_COL   = "metagenome_category"

ERA_A_YEARS = range(2000, 2012)          # 2000-11
ERA_B_YEARS = range(2012, 2024)          # 2012-23
ERA_LABELS  = ("2000-11", "2012-23")

# ---------------------------------------------------------------#
# 1 · Load both tables
# ---------------------------------------------------------------#
hits_df = pd.read_csv(HITS_CSV, low_memory=False)
meta_df = pd.read_csv(META_CSV, low_memory=False)

# ---------------------------------------------------------------#
# 2 · Common tidy-up
# ---------------------------------------------------------------#
def tidy(df):
    df = df.copy()
    df["collection_date_sam"] = pd.to_datetime(
        df["collection_date_sam"].str.strip("[]"), errors="coerce"
    )
    df = df.dropna(subset=["collection_date_sam", CAT_COL])
    df["year_bin"] = df["collection_date_sam"].dt.year
    return df

hits_df = tidy(hits_df)
meta_df = tidy(meta_df)

# keep only the drug classes of interest in the positive table
hits_df = (
    hits_df.assign(ARO_DrugClass=hits_df["ARO_DrugClass"]
                   .dropna().str.split(";"))
            .explode("ARO_DrugClass")
            .query("ARO_DrugClass in @DRUG_SET")
)

# ---------------------------------------------------------------#
# 3 · Build total- and positive-sample tables
# ---------------------------------------------------------------#
# 3·A  TOTAL samples  (all runs, unique per SAMPLE_ID)
total_samples_df = (
    meta_df.drop_duplicates(subset=[SAMPLE_ID])
           .groupby([CONT_COL, CAT_COL, "year_bin"], as_index=False)
           .agg(total_samples=(SAMPLE_ID, "count"))
)

# 3·B  POSITIVE samples for each drug class
positive_samples_df = (
    hits_df.drop_duplicates(subset=[SAMPLE_ID, "ARO_DrugClass"])
           .groupby(["ARO_DrugClass", CONT_COL, CAT_COL, "year_bin"],
                    as_index=False)
           .agg(positive_samples=(SAMPLE_ID, "count"))
)

# ---------------------------------------------------------------#
# 4 · Merge and compute prevalence
# ---------------------------------------------------------------#
df = positive_samples_df.merge(
        total_samples_df,
        on=[CONT_COL, CAT_COL, "year_bin"],
        how="left", validate="many_to_one"
     )
df["prevalence"] = df["positive_samples"] / df["total_samples"]

# ---------------------------------------------------------------#
# 5 · Collapse to the two eras
# ---------------------------------------------------------------#
def era_label(y):
    if y in ERA_A_YEARS: return ERA_LABELS[0]
    if y in ERA_B_YEARS: return ERA_LABELS[1]
    return np.nan

df["era"] = df["year_bin"].map(era_label)
df = df.dropna(subset=["era"])

era_df = (df.groupby(["ARO_DrugClass", CONT_COL, CAT_COL, "era"],
                     as_index=False)
            .agg(positive=("positive_samples", "sum"),
                 total    =("total_samples",    "sum"))
            .assign(prevalence=lambda d: d.positive / d.total)
)

# Pivot wide: each row has two era columns
wide = era_df.pivot_table(
            index=["ARO_DrugClass", CONT_COL, CAT_COL],
            columns="era",
            values=["positive", "total", "prevalence"]
       ).dropna(axis=0)

def ztest_row(row):
    pos_a, tot_a = row["positive"][ERA_LABELS[0]], row["total"][ERA_LABELS[0]]
    pos_b, tot_b = row["positive"][ERA_LABELS[1]], row["total"][ERA_LABELS[1]]
    stat, p = proportions_ztest([pos_a, pos_b], [tot_a, tot_b])
    return pd.Series({
        "prev_A": row["prevalence"][ERA_LABELS[0]],
        "prev_B": row["prevalence"][ERA_LABELS[1]],
        "delta" : row["prevalence"][ERA_LABELS[1]] -
                  row["prevalence"][ERA_LABELS[0]],
        "p"     : p
    })

slope_tbl = wide.apply(ztest_row, axis=1).reset_index()

# ---------------------------------------------------------------#
# A · compute within-era slopes  (one number per era per group)
# ---------------------------------------------------------------#
from scipy.stats import linregress

def slope_within_era(group):
    """Return slope & p of prevalence ~ year within one era."""
    if group['year_bin'].nunique() < 3:        # too few points
        return pd.Series({'slope': np.nan, 'p': np.nan})
    res = linregress(group['year_bin'], group['prevalence'])
    return pd.Series({'slope': res.slope, 'p': res.pvalue})

# attach era label first
df['era'] = df['year_bin'].map(era_label)
df_era    = df.dropna(subset=['era'])

slopes = (df_era
    .groupby(['era',               # 2000-11  /  2012-23
              'ARO_DrugClass',
              CONT_COL,
              CAT_COL])
    .apply(slope_within_era)
    .reset_index()
)

# ---------------------------------------------------------------#
# B · build the pivot table with multi-index columns
# ---------------------------------------------------------------#
# choose which era you want; loop if you want both in a single run
for era_label in ERA_LABELS:                 # ('2000-11', '2012-23')
    era_slopes = slopes.query("era == @era_label")

    # wide table: rows = drug, columns = (continent, category)
    table = era_slopes.pivot_table(
                index='ARO_DrugClass',
                columns=[CONT_COL, CAT_COL],
                values='slope'
            )

    # -----------------------------------------------------------#
    # C · optional significance mask (p ≥ 0.05 → fade the colour)
    # -----------------------------------------------------------#
    pvals = era_slopes.pivot_table(
                index='ARO_DrugClass',
                columns=[CONT_COL, CAT_COL],
                values='p'
            )
    mask = pvals >= 0.05       # True where slope is NOT significant

    # -----------------------------------------------------------#
    # D · plot as a heat-map
    # -----------------------------------------------------------#
    plt.figure(figsize=(12, 5))
    sns.heatmap(table,
                cmap='vlag',
                center=0,           # white = flat
                linewidths=.4,
                linecolor='grey',
                mask=mask,          # hide non-significant slopes
                cbar_kws=dict(label="Slope of prevalence / year"))
    plt.title(f"Trend inside era {era_label}\n"
              "(red = decline, blue = increase, blank = p≥0.05)")
    plt.tight_layout()
    plt.savefig(f"../data/within_era_trend_{era_label}.png", dpi=300)
    plt.savefig(f"../data/within_era_trend_{era_label}.svg")




'''
# ---------------------------------------------------------------#
# 6 · Draw the slope-graph
# ---------------------------------------------------------------#
def pstar(p):
    return "" if p >= .05 else "*" if p >= .01 else "**" if p >= .001 else "***"

fig, ax = plt.subplots(figsize=(8, 6))

for (_, row) in slope_tbl.iterrows():
    x = [0, 1]
    y = [row["prev_A"], row["prev_B"]]
    colour = "tab:red" if row["delta"] < 0 else "tab:blue"
    ax.plot(x, y, color=colour, linewidth=2,
            marker="o", markersize=6, alpha=.9)
    label = (f"{row[CONT_COL]} • {row['ARO_DrugClass']} "
             f"{pstar(row['p'])}")
    ax.text(1.03, y[1], label, va="center", fontsize=7)

ax.set_xticks([0, 1], ERA_LABELS)
ax.set_xlim(-.05, 1.40)
ax.set_ylim(0, 1.05)
ax.set_ylabel("Prevalence")
ax.set_title("AMR prevalence change\n"
             "blue = increase, red = decrease   "
             "* p<0.05  ** p<0.01  *** p<0.001",
             pad=12)
sns.despine(ax=ax)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=1200)
fig.savefig(OUT_SVG)
print("Saved:", OUT_PNG, "and", OUT_SVG)
'''
