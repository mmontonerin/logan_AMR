import pandas as pd

# Define the file paths
input_file = "./data/card_metadata_aro_geolocation_metagenomes_onlydateloc_noamplicon.csv"
output_file = "./data/card_metadata_aro_geolocation_metagenomes_onlydateloc_noamplicon_nouncalculated.csv"

df = pd.read_csv(input_file, low_memory=False)

# Filter out rows where assay_type is "AMPLICON"
df = df[df['geo_loc_name_country_continent_calc'] != 'uncalculated']

# Save the new table
df.to_csv(output_file, index=False)