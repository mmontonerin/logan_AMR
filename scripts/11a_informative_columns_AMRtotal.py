from pathlib import Path
from typing import Optional
import pandas as pd

# ───────────────────────────────
# file paths
# ───────────────────────────────
src_csv = Path("../data/card_metadata_aro.csv")
dst_csv = Path("../data/card_metadata_aro_informativecolumns.csv")

# ───────────────────────────────
# metagenome category map
# ───────────────────────────────
categories = {
    "human": [
        "human gut metagenome", "human metagenome", "human oral metagenome",
        "human skin metagenome", "human feces metagenome", "human vaginal metagenome",
        "human nasopharyngeal metagenome", "human lung metagenome",
        "human saliva metagenome", "human reproductive system metagenome",
        "human urinary tract metagenome", "human eye metagenome",
        "human blood metagenome", "human bile metagenome",
        "human tracheal metagenome", "human brain metagenome",
        "human milk metagenome", "human semen metagenome",
        "human skeleton metagenome"
    ],
    "livestock": [
        "bovine gut metagenome", "bovine metagenome", "pig gut metagenome",
        "pig metagenome", "chicken gut metagenome", "chicken metagenome",
        "sheep gut metagenome", "sheep metagenome"
    ],
    "marine": ["marine metagenome", "seawater metagenome"],
    "freshwater": ["freshwater metagenome", "lake water metagenome",
                   "groundwater metagenome"],
    "soil": ["soil metagenome"],
    "wastewater": ["wastewater metagenome"]
}

def get_category(organism: str | float) -> Optional[str]:
    """Return env category name or None."""
    if pd.isna(organism):
        return None
    for cat, values in categories.items():
        if organism in values:
            return cat
    return None

def assign_category(row: pd.Series) -> Optional[str]:
    """Rules 1–3 from your original code."""
    if row["organism_type"] == "Isolate":
        return pd.NA
    cat = get_category(row["organism"])
    return cat if cat is not None else "other"

# ───────────────────────────────
# streaming parameters
# ───────────────────────────────
chunk_size = 1_000_000           # adapt to your machine
needed_cols = ["acc", "organism"]  # read minimal set; we create the rest
acc_seen = set()
first_write = True

for chunk in pd.read_csv(src_csv, dtype=str, usecols=needed_cols,
                         chunksize=chunk_size):

    chunk = chunk.dropna(subset=["acc"])                 # discard rows w/o acc
    chunk = chunk[~chunk["acc"].isin(acc_seen)]          # keep only NEW accs
    acc_seen.update(chunk["acc"])

    # 1. organism_type
    chunk["organism_type"] = chunk["organism"]\
        .str.contains("metagenome", case=False, na=False)\
        .map({True: "Metagenome", False: "Isolate"})

    # 2. metagenome_category
    chunk["metagenome_category"] = chunk.apply(assign_category, axis=1)

    # 3. keep just the three columns in the required order
    out_chunk = chunk[["acc", "organism_type", "metagenome_category"]]

    # 4. append to output
    out_chunk.to_csv(
        dst_csv,
        mode="w" if first_write else "a",
        index=False,
        header=first_write
    )
    first_write = False
