import os, math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

amr_csv = "../data/full_card_metadata_aro_allfilters_metagenomes_WHOcategories.csv"
sra_csv = "../data/SRA_metadata_allfilters_logan.csv"
out_dir = "../data/who_radar_mean_prevalence"
os.makedirs(out_dir, exist_ok=True)

# Yearly bins not used, but might be needed later
def date_collection_metagenome(df):
    df = df.copy()
    df["collection_date_sam"] = pd.to_datetime(
        df["collection_date_sam"].str.strip("[]"), errors="coerce")
    
    required = ["collection_date_sam", "metagenome_category", "WHO_categories"]
    present  = [col for col in required if col in df.columns]
    df = df.dropna(subset=present)

    df["year_bin"] = df["collection_date_sam"].dt.year
    return df

# Stay between 2016‒2021 and keep only selected metagenome categories

keep_meta = ["livestock", "human", "wastewater", "marine", "freshwater", "soil"]

amr = date_collection_metagenome(pd.read_csv(amr_csv, low_memory=False))
sra = date_collection_metagenome(pd.read_csv(sra_csv, low_memory=False))

amr = amr.dropna(subset=["WHO_categories"])

mask_time = amr["collection_date_sam"].between("2016-01-01", "2021-12-31")
amr = amr.loc[mask_time & amr["metagenome_category"].isin(keep_meta)].copy()
mask_time = sra["collection_date_sam"].between("2016-01-01", "2022-12-31")
sra = sra.loc[mask_time & sra["metagenome_category"].isin(keep_meta)].copy()

if amr.empty or sra.empty:
    print("No rows after filtering")
    quit()

# lowercase for Access/Watch/Reserve detection
amr["who_class"] = amr["WHO_categories"].str.split(";")

# --------------------------------------------------------------
# explode WHO categories so each row has one WHO label
# --------------------------------------------------------------
amr = amr.explode("who_class").dropna(subset=["who_class"])

# fix naming to match five labels
wanted_labels = ["Access", "Watch", "Reserve", "Mixed", "Not classified"]
amr = amr[amr["who_class"].isin(wanted_labels)]

# --------------------------------------------------------------
# total & positive counts   (unique acc per metagenome category)
# --------------------------------------------------------------
sample_key = "acc"

total_cnt = (
    sra.drop_duplicates([sample_key])
       .groupby("metagenome_category")[sample_key]
       .nunique()
)

positive_cnt = (
    amr.drop_duplicates([sample_key, "who_class"])
       .groupby(["metagenome_category", "who_class"])[sample_key]
       .nunique()
).unstack(fill_value=0)

# align tables & compute prevalence matrix
positive_cnt = positive_cnt.reindex(index=keep_meta, columns=wanted_labels, fill_value=0)
total_cnt    = total_cnt.reindex(keep_meta, fill_value=0)

prevalence = positive_cnt.div(total_cnt, axis=0).fillna(0)

# --------------------------------------------------------------
# RADAR PLOT – all six meta categories on one figure
# --------------------------------------------------------------
n_axes = len(wanted_labels)
angles = np.linspace(0, 2*math.pi, n_axes, endpoint=False).tolist()
angles += angles[:1]                       # close the loop

fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))

# colour palette (same dark bases you used)
base_col = {
    "livestock":   "#fac723",
    "human":       "#936fac",
    "wastewater":  "#e95e50",
    "marine":      "#0cb2af",
    "freshwater":  "#a1c65d",
    "soil":        "#f29222",
}

for meta in keep_meta:
    vals = prevalence.loc[meta].tolist()
    vals += vals[:1]                       # close loop
    ax.plot(angles, vals, label=meta,
            linewidth=2.5, color=base_col.get(meta, "#333333"))
    ax.fill(angles, vals, alpha=0.15, color=base_col.get(meta, "#333333"))

# cosmetic tweaks
ax.set_xticks(angles[:-1])
ax.set_xticklabels(wanted_labels, fontsize=10)
ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_yticklabels(["0.2", "0.4", "0.6", "0.8", "1.0"], fontsize=8)
ax.set_ylim(0, 1.0)
ax.set_title("Average AMR prevalence (2016‒2022)\nWHO classes per metagenome type",
             pad=20)
ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.05))
fig.tight_layout()

outfile = os.path.join(out_dir, "radar_prevalence_human_livestock_etc.png")
fig.savefig(outfile, dpi=600)
fig.savefig(outfile.replace(".png", ".svg"))
plt.close()

print("Radar plot saved to:", outfile)

for yr in range(2016, 2022):               # 2021 inclusive
    amr_yr = amr[amr["collection_date_sam"].dt.year == yr]
    sra_yr = sra[sra["collection_date_sam"].dt.year == yr]

    total_cnt = (
        sra_yr.drop_duplicates(sample_key)
              .groupby("metagenome_category")[sample_key]
              .nunique()
              .reindex(keep_meta, fill_value=0)
    )

    positive_cnt = (
        amr_yr.drop_duplicates([sample_key, "who_class"])
              .groupby(["metagenome_category", "who_class"])[sample_key]
              .nunique()
              .unstack(fill_value=0)
              .reindex(index=keep_meta, columns=wanted_labels, fill_value=0)
    )

    prevalence = positive_cnt.div(total_cnt, axis=0).fillna(0)
    if prevalence.values.sum() == 0:
        print(f"{yr}: no data - skipping")
        continue

    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))
    for meta in keep_meta:
        vals = prevalence.loc[meta].tolist() 
        vals += vals[:1]
        ax.plot(angles, vals, linewidth=2.5, label=meta,
                color=base_col.get(meta, "#333333"))
        ax.fill(angles, vals, alpha=0.15, color=base_col.get(meta, "#333333"))

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(wanted_labels, fontsize=10)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_ylim(0, 1.0)
    ax.set_title(f"AMR prevalence ({yr})", pad=20)
    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.05))
    fig.tight_layout()

    fn = os.path.join(out_dir, f"radar_prevalence_{yr}.png")
    fig.savefig(fn, dpi=600)
    fig.savefig(fn.replace(".png", ".svg"))
    plt.close()

print("Radars saved to:", os.path.abspath(out_dir))