import pandas as pd

df_sra = pd.read_csv("../data/SRA_metadata_before20231211_logan.csv", dtype=str, engine='python')
output_file = "../data/SRA_metadata_before20231211_logan_extended.csv"
#diag_file1 = "../data/SRA_metagenomic_unexpected_counts.csv"
#diag_file2 = "../data/SRA_metagenomic_unexpected_librarysource_counts.csv"

# Define metagenome categories
categories = {
    "human": [
        "human gut metagenome", "human metagenome", "human oral metagenome", "human skin metagenome",
        "human feces metagenome", "human vaginal metagenome", "human nasopharyngeal metagenome", 
        "human lung metagenome", "human saliva metagenome", "human reproductive system metagenome", 
        "human urinary tract metagenome", "human eye metagenome", "human blood metagenome", 
        "human bile metagenome", "human tracheal metagenome", "human brain metagenome", 
        "human milk metagenome", "human semen metagenome", "human skeleton metagenome"
    ],
    "livestock": [
        "bovine gut metagenome", "bovine metagenome", "pig gut metagenome", "pig metagenome", 
        "chicken gut metagenome", "chicken metagenome", "sheep gut metagenome", "sheep metagenome"
    ],
    "marine": ["marine metagenome", "seawater metagenome"],
    "freshwater": ["freshwater metagenome", "lake water metagenome", "groundwater metagenome"],
    "soil": ["soil metagenome"],
    "wastewater": ["wastewater metagenome"]
}

categories_sp = {
    "livestock": [
        "Sus scrofa", "Sus scrofa domesticus", "Sus scrofa affinis", "Bos taurus", "Gallus gallus", "Equus caballus", 
        "Equs caballus", "Ovis aries", "Ovis", "Bos indicus", "Bos mutus", "Bos primigenius", "Bos frontalis",
        "Bos gaurus", "Gallus", "Capra hircus", "Capra aegagrus", "Capra ibex"
    ],
    "human": [
        "Homo sapiens"
    ]  
}

# Function to map organism to metagenome category
def get_category(organism: str) -> str | None:
    if pd.isna(organism):
        return None
    for category, values in categories.items():
        if organism in values:
            return category
    return None

# Function to map organism to metagenome category if organism is sp
def get_category_sp(organism: str) -> str | None:
    if pd.isna(organism):
        return None
    for category, values in categories_sp.items():
        if organism in values:
            return category
    return None


# organism_type  (Metagenome | Isolate)
is_metagenome = (
    df_sra["organism"].str.contains("metagenome", case=False, na=False) |
    df_sra["librarysource"].str.contains("METAGENOMIC|METATRANSCRIPTOMIC", case=False, na=False)
)

df_sra["organism_type"] = is_metagenome.map({True: "Metagenome", False: "Isolate"})


# Assign metagenome_category
# - Logic hierarchy:
#       1. when organism string already contains 'metagenome' → use get_category
#       2. else if librarysource matches the meta(genomic|transcriptomic) → use get_category_sp
#       3. Other metagenome categories → "other"
#    isolates get <NA>

# metagenome_category
def decide_category(row):
    if row["organism_type"] == "Isolate":
        return pd.NA

    # rule 1: organism string already includes 'metagenome'
    cat = get_category(row["organism"])
    if cat is not None:
        return cat

    # rule 2: librarysource categorisation (species names)
    if pd.notna(row["librarysource"]) and pd.Series(row["librarysource"]).str.contains(
        "METAGENOMIC|METATRANSCRIPTOMIC", case=False, regex=True
    ).bool():
        cat = get_category_sp(row["organism"])
        if cat is not None:
            return cat

    # rule 3: default bucket
    return "other"

df_sra["metagenome_category"] = df_sra.apply(decide_category, axis=1)

# Output file
df_sra.to_csv(output_file, index=False)


'''# Diagnostics ─ list & count librarysource values where:
# librarysource contains "METAGENOMIC
mask_unexpected = (
    ~df_sra["organism"].str.contains("metagenome", case=False, na=False)
)

librarysource_counts = (
    df_sra.loc[mask_unexpected, "librarysource"]
          .value_counts(dropna=False)
          .sort_values(ascending=False)
)   
librarysource_df = librarysource_counts.reset_index()
librarysource_df.columns = ["librarysource", "count"]       
# Save diagnostic librarysource counts to CSV
librarysource_df.to_csv(diag_file2, index=False, header=False)

# Diagnostics ─ list & count organism values where:
# librarysource contains "METAGENOMIC"
# BUT organism does *not* contain "metagenome"


mask_unexpected = (
    df_sra["librarysource"].str.contains("METAGENOMIC", case=False, na=False) &
    ~df_sra["organism"].str.contains("metagenome", case=False, na=False)
)

organism_counts = (
    df_sra.loc[mask_unexpected, "organism"]
          .value_counts(dropna=False)
          .sort_values(ascending=False)
)

diagnostic_df = organism_counts.reset_index()

diagnostic_df.columns = ["organism", "count"]

# Save diagnostic organism counts to CSV
diagnostic_df.to_csv(diag_file1, index=False)
'''