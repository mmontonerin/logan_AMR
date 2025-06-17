import pandas as pd

# Define input and output files
input_file = "../data/SRA_metadata_before20231211.csv"
output_file = "../data/SRA_metadata_before20231211_date_and_continent.csv"

df = pd.read_csv(input_file, low_memory=False)

# Filter out rows of data published after the date threshold
df = df[df['collection_date_sam'].notna() & df['geo_loc_name_country_continent_calc'].notna()
        & df['collection_date_sam'] != 'uncalculated' & df['geo_loc_name_country_continent_calc'] != 'uncalculated']

# Save the new table
df.to_csv(output_file, index=False)


