import pandas as pd

csv_df = pd.read_csv("../data/plasmids_sra_metadata.csv", dtype=str)       # has column "seq_name"
tsv_df = pd.read_csv("../data/argnorm_output.tsv", sep="\t", dtype=str)  # has column "Contig id"

# Change column name so that they match
tsv_df = tsv_df.rename(columns={"Contig id": "seq_name"})

amr_metadata_full = pd.merge(
    csv_df,              # left side: the big metadata table
    tsv_df,              # right side: AMR table
    on="seq_name",       # common key
    how="inner"          # keep only the matches (AMR-positive)
)

amr_metadata_full.to_csv("../data/amr_metadata_full.csv", index=False)
print(f"✔ Wrote {len(amr_metadata_full):,} records with AMR details to amr_metadata_full.csv")
