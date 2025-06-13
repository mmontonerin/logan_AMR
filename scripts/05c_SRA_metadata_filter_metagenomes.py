import pandas as pd

# Define input and output files
input_file = "../data/SRA_metadata_before20231211_date_and_continent.csv"
output_file = "../data/SRA_metadata_before20231211_date_and_continent_metagenomes.csv"

df = pd.read_csv(input_file, low_memory=False)

# Filter out rows of data published after the date threshold
df = df[df['organism'].notna() & df['organism'].str.contains("metagenome")]

# Save the new table
df.to_csv(output_file, index=False)