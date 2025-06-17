import pandas as pd

# Only keep assay types that contain most of the genome/gene sequences
assay_types_to_keep = [
    "WGS",
    "WGA",
    "RNA-Seq"
    "Synthetic-Long-Read",
    "WCS",
    "Hi-C",
    "ssRNA-seq",
    "FL-cDNA",
    "EST",
    "CLONE",
    "CLONEEND",
    "POOLCLONE"]

# Define input and output files
input_file = "../data/SRA_metadata_before20231211_date_and_continent_metagenomes.csv"
output_file = "../data/SRA_metadata_before20231211_date_and_continent_metagenomes_assaytype.csv"

df = pd.read_csv(input_file, low_memory=False)

# Filter out rows of data with assay types not in the list
df = df[df['assay_type'].notna() & df['assay_type'].isin(assay_types_to_keep)]

# Save the new table
df.to_csv(output_file, index=False)