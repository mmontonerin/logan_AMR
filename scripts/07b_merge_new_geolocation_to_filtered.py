import pandas as pd
import os
from datetime import datetime
import time

def log_progress(message):
    """Helper function to log progress with timestamp"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)

log_progress("Starting data processing")

# Step 1: Load SRA metadata as a dictionary (small file optimization)
log_progress("Loading geolocation data file into memory")
metadata_dict = {}

df = pd.read_csv("./data/biosample_geographical_location.202503.csv", dtype=str)

# Filter out rows that do not contain geolocation information
df = df[df['lat_lon'] != '0101000020E61000000AD7A3703D4A41C06666666666660CC0']

for _, row in df.iterrows():
    metadata_dict[row["accession"]] = {
        "Country": row.get("atribute_name", ""),
        "lat_lon": row.get("attribute_value", ""),
        "WKT": row.get("lat_lon", ""),
        "elevation": row.get("elevation", ""),
        "country_abv": row.get("country", ""),
        "biome": row.get("biome", ""),
        "geoloc_confidence_0-6": row.get("confidence", "")
    }

log_progress("Metadata loading complete")

# Step 2: Process the large alignment file in chunks and merge
output_file = "./data/card_metadata_aro_geolocation_metagenomes_onlydateloc.csv"
first_chunk = True  # To write the header only once
total_processed = 0
start_time = time.time()

for i, chunk in enumerate(pd.read_csv("./data/card_aro_entries_filtered_metagenomes.csv", dtype=str, chunksize=500000)):
    chunk_start = time.time()
    log_progress(f"Processing alignment chunk {i+1}")
    
    # Merge by mapping sample_name to metadata_dict
    chunk = chunk.assign(**{
        "Country": chunk["sample_name"].map(lambda x: metadata_dict.get(x, {}).get("Country", "")),
        "lat_lon": chunk["sample_name"].map(lambda x: metadata_dict.get(x, {}).get("lat_lon", "")),
        "WKT": chunk["sample_name"].map(lambda x: metadata_dict.get(x, {}).get("WKT", "")),
        "elevation": chunk["sample_name"].map(lambda x: metadata_dict.get(x, {}).get("elevation", "")),
        "country_abv": chunk["sample_name"].map(lambda x: metadata_dict.get(x, {}).get("country_abv", "")),
        "biome": chunk["sample_name"].map(lambda x: metadata_dict.get(x, {}).get("biome", "")),
        "geoloc_confidence_0-6": chunk["sample_name"].map(lambda x: metadata_dict.get(x, {}).get("geoloc_confidence_0-6", ""))
    })
    
    # Save merged chunk
    chunk.to_csv(output_file, mode='w' if first_chunk else 'a', index=False, header=first_chunk)
    first_chunk = False  # Ensure only the first chunk writes headers
    
    total_processed += len(chunk)
    log_progress(f"Chunk {i+1} processed in {time.time() - chunk_start:.2f} seconds")

log_progress(f"Processing complete. Total records processed: {total_processed}")