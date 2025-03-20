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
log_progress("Loading SRA metadata file into memory")
start_time = time.time()
metadata_dict = {}
total_metadata_rows = 0

for i, chunk in enumerate(pd.read_csv("./data/SRA_metadata.csv", dtype=str, chunksize=500000)):
    for _, row in chunk.iterrows():
        metadata_dict[row["acc"]] = row.to_dict()  # Store each row as a dictionary
    total_metadata_rows += len(chunk)
    log_progress(f"Processed metadata chunk {i+1} with {len(chunk)} rows. Total rows: {total_metadata_rows}")

metadata_time = time.time() - start_time
log_progress(f"Metadata loading complete. Loaded {total_metadata_rows} records in {metadata_time:.2f} seconds")

# Step 2: Process the large alignment file in chunks and merge
output_file = "./data/card_metadata.csv"
first_chunk = True  # To write the header only once
total_processed = 0
start_time = time.time()

# Get total rows in alignment file for progress reporting
log_progress("Calculating total rows in alignment file...")
total_rows = sum(1 for _ in pd.read_csv("./data/card_alignment_v1.1_contigs.csv", nrows=0, skiprows=lambda x: x > 0))
log_progress(f"Total rows in alignment file: {total_rows}")

for i, chunk in enumerate(pd.read_csv("./data/card_alignment_v1.1_contigs.csv", dtype=str, chunksize=500000)):
    chunk_start = time.time()
    log_progress(f"Processing alignment chunk {i+1} with {len(chunk)} rows")
    
    # Check for missing keys before merging
    missing_keys = [acc for acc in chunk["acc"] if acc not in metadata_dict]
    if missing_keys:
        missing_count = len(missing_keys)
        log_progress(f"Warning: {missing_count} accessions not found in metadata. First few: {missing_keys[:5]}")
    
    # Create a temporary dataframe from metadata for accessions in this chunk
    chunk_acc_list = chunk["acc"].unique().tolist()
    metadata_df = pd.DataFrame([metadata_dict.get(acc, {}) for acc in chunk_acc_list])
    
    if not metadata_df.empty:
        # Ensure we have an index column for the merge
        metadata_df["acc"] = chunk_acc_list
        
        # Check for column overlap to avoid unintentional overwrites
        overlap_columns = set(chunk.columns) & set(metadata_df.columns) - {"acc"}
        if overlap_columns:
            log_progress(f"Warning: Column overlap detected: {overlap_columns}. Columns from metadata will overwrite alignment columns.")
        
        # Perform left merge - without suffixes to keep original column names
        merged_chunk = pd.merge(
            chunk, 
            metadata_df,
            on="acc",
            how="left"
        )
        
        # Append chunk to file
        if first_chunk:
            # Create directory if it doesn't exist
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            
        merged_chunk.to_csv(output_file, index=False, mode="w" if first_chunk else "a", header=first_chunk)
        first_chunk = False  # Ensure header is written only for the first chunk
    
    total_processed += len(chunk)
    chunk_time = time.time() - chunk_start
    overall_progress = (total_processed / total_rows) * 100 if total_rows > 0 else 0
    
    log_progress(f"Chunk {i+1} processed in {chunk_time:.2f} seconds. Overall progress: {overall_progress:.1f}% ({total_processed}/{total_rows})")

total_time = time.time() - start_time
log_progress(f"Merging complete. Processed {total_processed} rows in {total_time:.2f} seconds")
log_progress(f"Results saved to {output_file}")