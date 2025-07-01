import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


# Inputs and output paths
amr_csv = "../data/full_card_metadata_aro_allfilters_metagenomes2.csv"
sra_csv = "../data/SRA_metadata_allfilters_logan.csv"
out_dir = "../data/who_prevalence_wastewater_europe"
os.makedirs(out_dir, exist_ok=True)

def date_collection_metagenome(df):
    df = df.copy()
    df["collection_date_sam"] = pd.to_datetime(
        df["collection_date_sam"].str.strip("[]"), errors="coerce")
    df = df.dropna(subset=["collection_date_sam",
                           "metagenome_category"])
    df["year_bin"] = df["collection_date_sam"].dt.year
    return df

# ------------------------------------------------------------------
# load & filter: human metagenomes - Europe only
# ------------------------------------------------------------------
amr_df = date_collection_metagenome(pd.read_csv(amr_csv, low_memory=False))
sra_df = date_collection_metagenome(pd.read_csv(sra_csv, low_memory=False))

amr_df = amr_df[
    (amr_df["metagenome_category"] == "wastewater") &
    (amr_df["geo_loc_name_country_continent_calc"] == "Europe")
]
sra_df = sra_df[
    (sra_df["metagenome_category"] == "wastewater") &
    (sra_df["geo_loc_name_country_continent_calc"] == "Europe")
]

if amr_df.empty or sra_df.empty:
    print("No wastewater-Europe rows in one of the tables - nothing to plot.")
    quit()

# ------------------------------------------------------------------
# quarter bins 2016-Q1 … 2022-Q4
# ------------------------------------------------------------------
amr_df["quarter_bin"] = amr_df["collection_date_sam"].dt.to_period("Q")
sra_df["quarter_bin"] = sra_df["collection_date_sam"].dt.to_period("Q")

timeline = pd.period_range("2016Q1", "2022Q4", freq="Q")

# ------------------------------------------------------------------
# WHO drug-class lists  (same as before)
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
    "Access":         "#4b8e0d",
    "Watch":          "#d95f02",
    "Reserve":        "#c51b1b",
    "Mixed":          "#6a4c93",
    "Not classified": "#999999"
}

# lowercase columns for quick matching
amr_df["_drug_lower"]    = amr_df["ARO_DrugClass"].str.lower().fillna("")
amr_df["_genefam_lower"] = amr_df["AMR_GeneFamily"].str.lower().fillna("")

# total samples (denominator) – one per acc per quarter
total_per_q = (
    sra_df.drop_duplicates(subset=["acc", "quarter_bin"])
          .groupby("quarter_bin")["acc"]
          .nunique()
          .reindex(timeline, fill_value=0)
)

# ------------------------------------------------------------------
# build one figure per WHO category
# ------------------------------------------------------------------
sns.set_theme(style="ticks")

for who_cat, drugs in drug_classes.items():
    base = category_base[who_cat]
    shades = sns.light_palette(base, n_colors=len(drugs) + 1, reverse=True)[:-1]

    fig, ax = plt.subplots(figsize=(10, 3))
    x = np.arange(len(timeline))
    year_pos  = [i for i, p in enumerate(timeline) if p.quarter == 1]
    year_lbls = [str(p.year) for p in timeline if p.quarter == 1]

    # -- one line per drug class ------------------------------------
    for col, drug in zip(shades, drugs):
        if drug == "colistin":
            mask = amr_df["_genefam_lower"].str.contains("colistin")
        else:
            mask = amr_df["_drug_lower"].str.contains(drug)

        hits = (
            amr_df.loc[mask, ["quarter_bin", "acc"]]
                  .drop_duplicates()
                  .groupby("quarter_bin")["acc"]
                  .nunique()
                  .reindex(timeline, fill_value=0)
        )

        prevalence = hits / total_per_q.replace(0, np.nan)   # avoid /0
        ax.plot(x, prevalence.values, label=drug, color=col, linewidth=2)

    # -- cosmetics ---------------------------------------------------
    ax.set_ylim(0, 1)
    ax.set_xticks(year_pos)
    ax.set_xticklabels(year_lbls, rotation=0)
    ax.set_ylabel("prevalence (0‒1)")
    ax.set_xlabel("")
    ax.set_title(f"{who_cat} AMR-positive prevalence\nWastewater metagenome (Europe)")
    ax.legend(title="drug class", bbox_to_anchor=(1.02, 1), loc="upper left")
    sns.despine()
    fig.tight_layout()

    fname = os.path.join(out_dir,
                         f"prev_{who_cat.replace(' ', '_').lower()}_wastewater_europe.png")
    fig.savefig(fname, dpi=600)
    fig.savefig(fname.replace(".png", ".svg"))
    plt.close()

print("plots saved to:", os.path.abspath(out_dir))
