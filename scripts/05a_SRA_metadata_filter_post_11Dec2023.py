import pandas as pd

# Define input and output files
input_file = "../data/SRA_metadata.csv"
output_file = "../data/SRA_metadata_before20231211.csv"
date_threshold = "2023-12-11"

df = pd.read_csv(input_file, low_memory=False, parse_dates=["releasedate"])

# Filter out rows of data published after the date threshold
df = df[df['releasedate'] <= date_threshold]

# Save the new table
df.to_csv(output_file, index=False)


