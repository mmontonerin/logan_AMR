import pandas as pd
import sys
sys.stdout.flush()

# Only keep assay types that contain most of the genome/gene sequences
assay_types_to_keep = [
    "WGS",
    "WGA",
    "RNA-Seq"
    "Synthetic-Long-Read",
    "WCS",
    "Hi-C",
    "ssRNA-seq",
    "FL-cDNA",
    "EST",
    "CLONE",
    "CLONEEND",
    "POOLCLONE"]

def process_csv_files(input_file, output_file, chunk_size=30000000):
    # Create empty DataFrame to store filtered results
    filtered_df = pd.DataFrame()

    # Process card-aro-metadata CSV in chunks
    with pd.read_csv(input_file, chunksize=chunk_size, low_memory=False) as reader:
        for i, chunk in enumerate(reader):
            print(f"Processing chunk {i + 1}...")

            # Keep rows of data containing sampling date and location data
            # Filter out rows of data with organisms that are not metagenomes
            # Only keep assay types that contain most of the genome sequences
            filtered_chunk = chunk[
                chunk['collection_date_sam'].notna() &
                (chunk['collection_date_sam'] != 'uncalculated') &
                chunk['geo_loc_name_country_continent_calc'].notna() &
                (chunk['geo_loc_name_country_continent_calc'] != 'uncalculated') &
                chunk['organism'].notna()
            ]

            # Append the filtered chunk to the results DataFrame
            filtered_df = pd.concat([filtered_df, filtered_chunk], ignore_index=True)

            # Print progress
            print(f"Processed chunk {i + 1} with {len(filtered_chunk)} rows, current total count: {len(filtered_df)}")

    # Save the filtered DataFrame to a CSV file
    filtered_df.to_csv(output_file, index=False)
    print(f"Processing complete. Filtered data saved to {output_file}")

if __name__ == "__main__":
    # Define input and output files
    input_file = "../data/card_metadata_aro.csv"
    output_file = "../data/new_full_card_metadata_aro_dateloc_organism.csv"
    
    # chunk size based on number of rows to be read at a time in the CSV file
    # Adapt to available memory
    process_csv_files(input_file, output_file, chunk_size=30000000)