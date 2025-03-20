import pandas as pd

# Define input and output files
metadata_file = "./data/SRA_metadata.csv"
summary_file = "./data/summary_statistics.csv"
date_threshold = "2023-12-11"

# Read CSV in chunks and filter
filtered_count = 0
chunk_size = 1000000

for chunk in pd.read_csv(metadata_file, chunksize=chunk_size, parse_dates=["releasedate"]):
    filtered_chunk = chunk[chunk["releasedate"] <= date_threshold]
    filtered_count += len(filtered_chunk)

# Load summary file and add new column
summary_df = pd.read_csv(summary_file)

# Add the count as a new column
summary_df["filtered_SRA_count"] = filtered_count

# Save updated summary file
summary_df.to_csv(summary_file, index=False)

print(f"Filtered sample count: {filtered_count} added to {summary_file}")
