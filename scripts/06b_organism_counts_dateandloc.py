import pandas as pd
import sys
sys.stdout.flush()

def process_csv_files(input_file, output_file, full_counts_file, chunk_size=1000000):
    """
    Process large CSV files in chunks to generate organism counts, but only for rows where both
    'collection_date_sam' and 'geo_loc_name_country_continent_calc' are not NA.
    
    Parameters:
    -----------
    input_file : str
        Path to the card_metadata CSV file.
    output_file : str
        Path to save the output CSV file with top 15 organism counts and "other".
    full_counts_file : str
        Path to save the CSV file with the full organism counts.
    chunk_size : int
        Number of rows to process at a time.
    """
    # Initialize a dictionary to store counts of each organism
    organism_counts = {}

    # Process CSV in chunks
    with pd.read_csv(input_file, chunksize=chunk_size, low_memory=False) as reader:
        for i, chunk in enumerate(reader):
            print(f"Processing chunk {i + 1}...")

            # Filter rows where both 'collection_date_sam' and 'geo_loc_name_country_continent_calc' are not NaN
            valid_chunk = chunk[chunk['collection_date_sam'].notna() & chunk['geo_loc_name_country_continent_calc'].notna()]

            # Count occurrences of each "organism" in the filtered chunk
            for organism in valid_chunk['organism'].dropna():
                organism_counts[organism] = organism_counts.get(organism, 0) + 1

            # Print progress
            sys.stdout.flush()
            print(f"Processed chunk {i + 1} with {len(valid_chunk)} rows, current total count: {sum(organism_counts.values())}")

    # Create a DataFrame for full organism counts
    full_counts_df = pd.DataFrame(organism_counts.items(), columns=["organism", "count"])

    # Save the full organism counts to a CSV
    full_counts_df.to_csv(full_counts_file, index=False)

    # Sort organisms by count and get the top 15
    top_organisms = full_counts_df.sort_values(by="count", ascending=False).head(15)

    # Calculate the "other" count (all remaining organisms not in the top 15)
    other_count = full_counts_df.loc[~full_counts_df['organism'].isin(top_organisms['organism']), 'count'].sum()

    # Create a DataFrame for the "other" row
    other_organisms = set(full_counts_df.loc[~full_counts_df['organism'].isin(top_organisms['organism']), 'organism'])
    other_row = pd.DataFrame({"organism": ["other"], "count": [other_count], "unique_in_other": [len(other_organisms)]})
    
    # Combine the top 15 organisms with the "other" row
    top_organisms = pd.concat([top_organisms, other_row], ignore_index=True)

    # Save the top 15 organisms and "other" to the results CSV file
    top_organisms.to_csv(output_file, index=False)

    # Print the final confirmation
    print(f"Processing complete. Table saved to {output_file}")
    print(f"Full counts saved to {full_counts_file}")
    print(top_organisms)

# Example usage
if __name__ == "__main__":
    # Replace with actual file paths
    input_file = "./data/card_metadata_unique.csv"
    output_file = "./data/organism_counts_top15_with_coldate_and_loc.csv"
    full_counts_file = "./data/organism_counts_full_with_coldate_and_loc.csv"
    
    # Adjust chunk size based on available memory
    process_csv_files(input_file, output_file, full_counts_file, chunk_size=1000000)
