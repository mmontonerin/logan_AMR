import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

# ------------------------------------------------------------------
# 0  CONFIG
# ------------------------------------------------------------------
input_file = "../data/full_card_metadata_aro_allfilters_metagenomes_WHOcategories.csv"
out_dir    = "../data/who_category_timelines_human_europe"
os.makedirs(out_dir, exist_ok=True)

# helper : same as before ------------------------------------------
def date_collection_metagenome(df):
    df = df.copy()
    df["collection_date_sam"] = pd.to_datetime(
        df["collection_date_sam"].str.strip("[]"), errors="coerce"
    )
    df = df.dropna(subset=[
        "collection_date_sam", "metagenome_category", "WHO_categories"
    ])
    df["year_bin"] = df["collection_date_sam"].dt.year
    return df

# ------------------------------------------------------------------
# 1  LOAD + PRE-PROCESS  (add human / Europe filter)
# ------------------------------------------------------------------
df = pd.read_csv(input_file, low_memory=False)
df = date_collection_metagenome(df)

# -------- filter to human metagenomes collected in Europe ----------
df = df[
    (df["metagenome_category"] == "human") &
    (df["geo_loc_name_country_continent_calc"] == "Europe")
]

# if nothing left, exit quietly
if df.empty:
    print("No human-Europe samples found — nothing to plot.")
    quit()

df["quarter_bin"] = df["collection_date_sam"].dt.to_period("Q")

df["_drug_lower"]    = df["ARO_DrugClass"].str.lower().fillna("")
df["_genefam_lower"] = df["AMR_GeneFamily"].str.lower().fillna("")

# ------------------------------------------------------------------
# 2  DEFINE WINDOW 2016-2022 (re-compute after filtering)
# ------------------------------------------------------------------
all_quarters = pd.period_range(df["quarter_bin"].min(),
                               df["quarter_bin"].max(),
                               freq="Q")

start_q = pd.Period("2016Q1", freq="Q")
end_q   = pd.Period("2022Q4", freq="Q")
all_quarters = all_quarters[(all_quarters >= start_q) & (all_quarters <= end_q)]

if all_quarters.empty:
    print("Filtered data contain no quarters between 2016Q1–2022Q4 — nothing to plot.")
    quit()

# ------------------------------------------------------------------
# 3  DRUG CLASS & COLOUR DEFINITIONS
# ------------------------------------------------------------------
drug_classes = {
    "Access": [
        "penicillin beta-lactam", "tetracycline antibiotic",
        "nitrofuran antibiotic", "nitroimidazole antibiotic"
    ],
    "Watch": [
        "fluoroquinolone antibiotic", "macrolide antibiotic",
        "glycopeptide antibiotic", "carbapenem",
        "lincosamide antibiotic", "monobactam"
    ],
    "Reserve": [
        "glycylcycline", "phosphonic acid antibiotic",
        "oxazolidinone antibiotic", "colistin"
    ],
    "Mixed": [
        "aminoglycoside antibiotic", "cephalosporin"
    ],
    "Not classified": [
        "isoniazid-like antibiotic"
    ]
}

category_base = {
    "Access":         "#4b8e0d",  # green
    "Watch":          "#d95f02",  # orange
    "Reserve":        "#c51b1b",  # red
    "Mixed":          "#6a4c93",  # purple
    "Not classified": "#999999"   # grey
}

# ------------------------------------------------------------------
# 4  PLOTS
# ------------------------------------------------------------------
for who_cat, drugs in drug_classes.items():
    base_colour  = category_base[who_cat]
    palette_full = sns.light_palette(base_colour, n_colors=len(drugs) + 3, reverse=True)
    palette      = palette_full[:-1]          # drop palest shade

    x         = np.arange(len(all_quarters))
    year_pos  = [i for i, p in enumerate(all_quarters) if p.quarter == 1]
    year_labs = [str(p.year) for p in all_quarters if p.quarter == 1]

    fig, ax = plt.subplots(figsize=(12, 4))

    for i, drug in enumerate(drugs):
        if drug == "colistin":
            mask = df["_genefam_lower"].str.contains("colistin")
        else:
            mask = df["_drug_lower"].str.contains(drug)

        subset = (
            df.loc[mask, ["quarter_bin", "acc"]]
              .drop_duplicates()                          # one per (quarter, acc)
        )

        timeline = (
            subset.groupby("quarter_bin")["acc"]
                   .nunique()
                   .reindex(all_quarters, fill_value=0)
                   .sort_index()
        )

        ax.plot(x, timeline.values,
                label=drug,
                color=palette[i],
                linewidth=2)

    ax.set_xticks(year_pos)
    ax.set_xticklabels(year_labs, rotation=0)
    ax.set_title(f"{who_cat} AMR-positive samples timeline\nHuman metagenome (Europe)")
    ax.set_ylabel("Unique SRA samples")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=10))
    ax.tick_params(axis="x", rotation=45)
    ax.legend(title="drug class", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()

    fname = os.path.join(out_dir,
                         f"timeline_{who_cat.replace(' ', '_').lower()}_human_europe.png")
    fig.savefig(fname, dpi=600)
    fig.savefig(fname.replace(".png", ".svg"))
    plt.close()

print("plots written to:", os.path.abspath(out_dir))
