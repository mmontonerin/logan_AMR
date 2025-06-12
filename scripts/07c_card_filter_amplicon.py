import pandas as pd

# Define the file paths
input_file = "./data/card_metadata_aro_geolocation_metagenomes_onlydateloc.csv"
output_file = "./data/card_metadata_aro_geolocation_metagenomes_onlydateloc_noamplicon.csv"

df = pd.read_csv(input_file, low_memory=False)

# Filter out rows where assay_type is "AMPLICON"
df = df[df['assay_type'] != 'AMPLICON']

# Save the new table
df.to_csv(output_file, index=False)