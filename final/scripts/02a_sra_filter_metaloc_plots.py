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

def process_csv_files(input_file, output_file, chunk_size=20000000):
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
                chunk['assay_type'].isin(assay_types_to_keep) &
                (chunk['organism_type'] != 'Isolate')
            ]

            # Append the filtered chunk to the results DataFrame
            filtered_df = pd.concat([filtered_df, filtered_chunk], ignore_index=True)

    # Save the filtered DataFrame to a CSV file
    filtered_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    # Define input and output files
    input_file = "../data/SRA_metadata_before20231211_logan_extended.csv"
    output_file = "../data/SRA_metadata_before20231211_logan_dateloc_meta.csv"
    
    process_csv_files(input_file, output_file, chunk_size=20000000)