import pandas as pd
import sys
sys.stdout.flush()

# Define input and output file paths
input_file = "./data/card_metadata_filtered.csv"
output_file = "./data/card_metadata_unique.csv"

# Define chunk size
chunk_size = 1000000

# Use a set to track unique 'acc' values
temp_storage = set()

# Process CSV in chunks and write unique rows
with pd.read_csv(input_file, chunksize=chunk_size, low_memory=False) as reader:
    for i, chunk in enumerate(reader):
        print(f"Processing chunk {i + 1}...")
        chunk.drop(columns=["ARO_ID"], inplace=True)
        
        # Remove duplicates within the chunk based on 'acc'
        chunk.drop_duplicates(subset=["acc"], inplace=True)
        
        # Filter out rows with 'acc' values already seen
        unique_chunk = chunk.loc[~chunk["acc"].isin(temp_storage)]
        
        # Update temporary storage
        temp_storage.update(unique_chunk["acc"].tolist())
        
        # Write to file, adding header only for the first chunk
        unique_chunk.to_csv(output_file, mode='w' if i == 0 else 'a', index=False, header=(i == 0))
        print(f"Chunk {i + 1} processed and written to file.")

print("Processing complete.")

