import pandas as pd
import os
from datetime import datetime
import time
import csv
import re

def log_progress(message):
    """Helper function to log progress with timestamp"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)

log_progress("Starting data processing")

# File paths
aro_index_file = "./data/aro_index.tsv"
card_metadata_file = "./data/card_metadata_filtered.csv"
output_file = "./data/card_metadata_aro.csv"
debug_file = "./data/debug_missing_aros.csv"

# Step 1: Load ARO index as a dictionary
log_progress(f"Loading ARO index file into memory from {aro_index_file}")
start_time = time.time()
aro_dict = {}
aro_ids_in_index = set()
total_aro_rows = 0

# Load the TSV file
with open(aro_index_file, 'r') as tsv_file:
    reader = csv.reader(tsv_file, delimiter='\t')
    headers = next(reader)  # Get the headers
    
    # Find the indices of the columns we need
    try:
        aro_id_idx = headers.index('ARO Accession')
        prot_acc_idx = headers.index('Protein Accession')
        amr_gene_family_idx = headers.index('AMR Gene Family')
        drug_class_idx = headers.index('Drug Class')
        resistance_idx = headers.index('Resistance Mechanism')
        
        log_progress(f"Found all required columns in the ARO index file")
    except ValueError as e:
        log_progress(f"Error finding columns: {e}")
        log_progress(f"Available columns: {headers}")
        raise
    
    # Create the dictionary with the specific columns we need
    for row in reader:
        if len(row) > max(aro_id_idx, prot_acc_idx, amr_gene_family_idx, drug_class_idx, resistance_idx):
            # Store with original key
            original_id = row[aro_id_idx].strip()
            aro_dict[original_id] = {
                'ARO_ProtAccession': row[prot_acc_idx],
                'AMR_GeneFamily': row[amr_gene_family_idx],
                'ARO_DrugClass': row[drug_class_idx],
                'ARO_ResistanceMechanism': row[resistance_idx]
            }
            aro_ids_in_index.add(original_id)
            total_aro_rows += 1

aro_time = time.time() - start_time
log_progress(f"ARO index loading complete. Loaded {total_aro_rows} records in {aro_time:.2f} seconds")
log_progress(f"Number of unique ARO IDs in index: {len(aro_ids_in_index)}")
log_progress(f"First 5 ARO IDs from index: {list(aro_ids_in_index)[:5]}")

# Step 2: Process metadata and collect information about missing ARO IDs
log_progress(f"First pass: collecting information about missing ARO IDs")
missing_aros = {}
total_rows = 0
unique_aros_in_metadata = set()

for i, chunk in enumerate(pd.read_csv(card_metadata_file, dtype=str, chunksize=500000)):
    log_progress(f"Scanning chunk {i+1} with {len(chunk)} rows")
    
    if 'ARO_ID' not in chunk.columns:
        log_progress(f"ERROR: ARO_ID column not found in chunk {i+1}. Columns: {chunk.columns.tolist()}")
        continue
    
    # Clean ARO_ID values (remove any whitespace)
    chunk['ARO_ID'] = chunk['ARO_ID'].astype(str).str.strip()
    
    # Record all unique ARO_IDs in this chunk
    chunk_aros = set(chunk['ARO_ID'].unique())
    unique_aros_in_metadata.update(chunk_aros)
    
    # Check which ones are missing
    for aro_id in chunk_aros:
        if aro_id not in aro_ids_in_index:
            if aro_id not in missing_aros:
                missing_aros[aro_id] = 0
            missing_aros[aro_id] += chunk['ARO_ID'].eq(aro_id).sum()
    
    total_rows += len(chunk)

log_progress(f"Total unique ARO_IDs in metadata: {len(unique_aros_in_metadata)}")
log_progress(f"Total missing unique ARO_IDs: {len(missing_aros)}")

# Save debug info about missing ARO IDs
log_progress(f"Saving debug information about missing ARO IDs to {debug_file}")
missing_df = pd.DataFrame([
    {"ARO_ID": k, "Count": v, "In_Metadata": k in unique_aros_in_metadata, "In_Index": k in aro_ids_in_index}
    for k, v in missing_aros.items()
])
if not missing_df.empty:
    missing_df = missing_df.sort_values("Count", ascending=False)
    os.makedirs(os.path.dirname(debug_file), exist_ok=True)
    missing_df.to_csv(debug_file, index=False)
    log_progress(f"Top 10 missing ARO IDs by frequency: {missing_df.head(10).to_dict('records')}")

# Step 3: Now do the actual merging
log_progress(f"Starting the actual merging process")
log_progress(f"Output will be saved to {output_file}")

first_chunk = True
total_processed = 0
missing_count = 0
start_time = time.time()

for i, chunk in enumerate(pd.read_csv(card_metadata_file, dtype=str, chunksize=500000)):
    chunk_start = time.time()
    log_progress(f"Processing card_metadata chunk {i+1} with {len(chunk)} rows")
    
    if 'ARO_ID' not in chunk.columns:
        log_progress(f"ERROR: ARO_ID column not found in chunk {i+1}. Columns: {chunk.columns.tolist()}")
        continue
        
    # Clean ARO_ID values (remove any whitespace)
    chunk['ARO_ID'] = chunk['ARO_ID'].astype(str).str.strip()
    
    # Create columns for the ARO data we want to add
    chunk['ARO_ProtAccession'] = None
    chunk['AMR_GeneFamily'] = None
    chunk['ARO_DrugClass'] = None
    chunk['ARO_ResistanceMechanism'] = None
    
    # Populate the new columns from our dictionary
    chunk_missing = 0
    for idx, row in chunk.iterrows():
        aro_id = row['ARO_ID']
        
        if aro_id in aro_dict:
            for key, value in aro_dict[aro_id].items():
                chunk.at[idx, key] = value
        else:
            chunk_missing += 1
    
    missing_count += chunk_missing
    if chunk_missing > 0:
        log_progress(f"Warning: {chunk_missing} ARO_IDs not found in ARO index in this chunk.")
    
    # Write the chunk to the output file
    if first_chunk:
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    chunk.to_csv(output_file, index=False, mode="w" if first_chunk else "a", header=first_chunk)
    first_chunk = False  # Ensure header is written only for the first chunk
    
    total_processed += len(chunk)
    chunk_time = time.time() - chunk_start
    
    log_progress(f"Chunk {i+1} processed in {chunk_time:.2f} seconds. Total rows processed so far: {total_processed}")

total_time = time.time() - start_time
log_progress(f"Merging complete. Processed {total_processed} rows in {total_time:.2f} seconds")
log_progress(f"Total ARO_IDs not found: {missing_count} ({(missing_count/total_processed)*100:.2f}% of total)")
log_progress(f"Results saved to {output_file}")
log_progress(f"For detailed analysis of missing ARO IDs, check {debug_file}")