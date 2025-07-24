import pandas as pd

sra_db = pd.read_csv("../../data/SRA_metadata_before20231211.csv", dtype=str)
plasmid_db = pd.read_csv("../data/plasmid_data.tsv", sep="\t", dtype=str)

# SRA table has accession numbers in the format "ERR2138710"
# Plasmid table has contig numbers in the format "ERR2138710_20"
plasmid_db["acc"] = plasmid_db["seq_name"].str.split("_").str[0]

# Instead of printing the entire plasmid sequence, I will keep the plasmid length, in case there is some correlations in AMR findings to that
plasmid_db["plasmid_length"] = plasmid_db["sequence"].str.len()

sra_columns_to_keep = [
    "acc",
    "assay_type",
    "bioproject",
    "biosample",
    "organism",
    "librarysource",
    "collection_date_sam",
    "geo_loc_name_country_calc",
    "geo_loc_name_country_continent_calc",
    "releasedate"
]

merged_db = (
    plasmid_db[["acc", "seq_name", "plasmid_length"]]
        .merge(sra_db[sra_columns_to_keep], on="acc", how="left")
)

merged_db.to_csv("../data/plasmids_sra_metadata.csv", index=False)
print("Merged plasmid and SRA metadata saved to 'plasmids_sra_metadata.csv")



