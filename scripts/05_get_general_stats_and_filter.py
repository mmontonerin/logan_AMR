import pandas as pd
import sys
sys.stdout.flush()

def process_csv_files(sra_metadata_path, card_metadata_path, chunk_size=100000):
    """
    Process large CSV files in chunks to generate summary statistics and filtered datasets.
    
    Parameters:
    -----------
    sra_metadata_path : str
        Path to SRA_metadata.csv
    card_metadata_path : str
        Path to card_metadata_filtered.csv
    chunk_size : int
        Number of rows to process at a time
    """
    # Initialize counters and sets
    total_sra_count = 0
    unique_accessions = set()
    with_dates = set()
    with_country = set()
    with_continent = set()
    
    # Process SRA_metadata.csv to get total count
    print("Processing SRA_metadata.csv...")
    for chunk in pd.read_csv(sra_metadata_path, chunksize=chunk_size, low_memory=False):
        total_sra_count += len(chunk)
    print(f"Total SRA entries: {total_sra_count}")
    
    # Process card_metadata_filtered.csv
    print("Processing card_metadata_aro.csv...")
    
    # Create files for filtered data
    with_dates_file = open('./data/entries_with_dates.csv', 'w')
    with_dates_continent_file = open('./data/entries_with_dates_continent.csv', 'w')
    
    # Create the iterator outside the loop
    chunks_iterator = pd.read_csv(card_metadata_path, chunksize=chunk_size, low_memory=False)
    

    for i, chunk in enumerate(chunks_iterator):
        # Write headers to filtered files if this is the first chunk
        if i == 0:
            chunk.iloc[0:0].to_csv(with_dates_file, index=False)
            chunk.iloc[0:0].to_csv(with_dates_continent_file, index=False)

        # Process accessions
        acc_values = chunk['acc'].dropna().unique()
        unique_accessions.update(acc_values)
        
        # Identify rows with dates
        date_mask = chunk['collection_date_sam'].notna()
        date_rows = chunk[date_mask]
        date_accs = date_rows['acc'].unique()
        with_dates.update(date_accs)
        
        # Write rows with dates to file
        if not date_rows.empty:
            date_rows.to_csv(with_dates_file, mode='a', header=(i==0), index=False)

        
        # Identify rows with country
        country_mask = chunk['geo_loc_name_country_calc'].notna()
        country_accs = chunk.loc[country_mask, 'acc'].unique()
        with_country.update(country_accs)
        
        # Identify rows with continent
        continent_mask = chunk['geo_loc_name_country_continent_calc'].notna()
        continent_accs = chunk.loc[continent_mask, 'acc'].unique()
        with_continent.update(continent_accs)
        
        # Identify rows with both date and continent
        date_continent_mask = date_mask & continent_mask
        date_continent_rows = chunk[date_continent_mask]
        if not date_continent_rows.empty:
            date_continent_rows.to_csv(with_dates_continent_file, mode='a', header=(i==0), index=False)
        
        print(f"Processed chunk {i+1} with {len(chunk)} rows")
    
    # Close files
    with_dates_file.close()
    with_dates_continent_file.close()
    
    # Calculate combined sets
    with_date_continent = with_dates.intersection(with_continent)
    with_date_country_continent = with_dates.intersection(with_country).intersection(with_continent)

    # Create summary table
    summary_df = pd.DataFrame({
        "Total_SRA": [total_sra_count],
        "Samples_with_AMR": [len(unique_accessions)],
        "col_date_available": [len(with_dates)],
        "country_available": [len(with_country)],
        "continent_available": [len(with_continent)],
        "date_and_continent_available": [len(with_date_continent)],
        "date_country_continent_available": [len(with_date_country_continent)]
    })
    
    # Save summary to CSV
    summary_df.to_csv('./data/summary_statistics.csv', index=False)
    print("Summary statistics saved to ./data/summary_statistics.csv")
    print("Filtered data saved to ./data/with_dates.csv and ./data/with_dates_continent.csv")
    
    return summary_df

# Example usage
if __name__ == "__main__":
    # Replace with actual file paths
    sra_path = "./data/SRA_metadata.csv"
    card_path = "./data/card_metadata_aro.csv"
    
    # Adjust chunk size based on available memory
    result = process_csv_files(sra_path, card_path, chunk_size=100000)
    print(result)