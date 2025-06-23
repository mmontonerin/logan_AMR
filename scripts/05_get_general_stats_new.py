import pandas as pd
import sys
sys.stdout.flush()

def process_csv_files(sra_metadata_path, card_metadata_path, chunk_size=30000000):
    # Initialize counters and sets
    total_sra_count = 0
    unique_accessions = set()
    with_dates = set()
    with_continent = set()
    with_date_continent = set()
    with_metagenome = set()
    with_metagenome_assaytype = set()
    with_metagenomecategory = set()

    # Process SRA_metadata.csv to get total count
    print("Processing SRA_metadata.csv...")
    for chunk in pd.read_csv(sra_metadata_path, chunksize=chunk_size, low_memory=False):
        total_sra_count += len(chunk)
    print(f"Total SRA entries: {total_sra_count}")
    
    # Process card_metadata_filtered.csv
    print("Processing card_metadata_aro.csv...")
    
    # Create the iterator outside the loop
    chunks_iterator = pd.read_csv(card_metadata_path, chunksize=chunk_size, low_memory=False)


    for i, chunk in enumerate(chunks_iterator):
        # Process accessions
        acc_values = chunk['acc'].dropna().unique()
        unique_accessions.update(acc_values)
        
        # Identify rows with dates
        date_mask = chunk['collection_date_sam'].notna() & chunk['collection_date_sam'] != 'uncalculated'
        date_rows = chunk[date_mask]
        date_accs = date_rows['acc'].unique()
        with_dates.update(date_accs)
        
        # Identify rows with continent
        continent_mask = chunk['geo_loc_name_country_continent_calc'].notna() & chunk['geo_loc_name_country_continent_calc'] != 'uncalculated'
        continent_accs = chunk.loc[continent_mask, 'acc'].unique()
        with_continent.update(continent_accs)
        
         # Identify rows with continent and date
        continentdate_mask = continent_mask & date_mask
        continentdate_accs = chunk.loc[continentdate_mask, 'acc'].unique()
        with_date_continent.update(continentdate_accs)

        # Identify rows with metagenome        
        metagenome_mask = chunk['organism'].notna() & chunk['organism'].str.contains("metagenome", case=False) & continentdate_mask
        metagenome_accs = chunk.loc[metagenome_mask, 'acc'].unique()
        with_metagenome.update(metagenome_accs)

        # Identify rows with metagenome and assay type
        metagenome_assaytype_mask = chunk['assay_type'].notna() & \
            chunk['assay_type'].isin([
                "WGS", "WGA", "RNA-Seq", "Synthetic-Long-Read",
                "WCS", "Hi-C", "ssRNA-seq", "FL-cDNA", "EST",
                "CLONE", "CLONEEND", "POOLCLONE"
            ]) & metagenome_mask
        metagenome_assaytype_accs = chunk.loc[metagenome_assaytype_mask, 'acc'].unique()
        with_metagenome_assaytype.update(metagenome_assaytype_accs)

        # Identify rows with metagenome category
        metagenomecategory_mask = chunk['organism'].notna() & chunk['organism'].isin([
            "human gut metagenome", "human metagenome", "human oral metagenome", "human skin metagenome",
            "human feces metagenome", "human vaginal metagenome", "human nasopharyngeal metagenome",
            "human lung metagenome", "human saliva metagenome", "human reproductive system metagenome",
            "human urinary tract metagenome", "human eye metagenome", "human blood metagenome",
            "human bile metagenome", "human tracheal metagenome", "human brain metagenome",
            "human milk metagenome", "human semen metagenome", "human skeleton metagenome",
            "bovine gut metagenome", "bovine metagenome", "pig gut metagenome", "pig metagenome",
            "chicken gut metagenome", "chicken metagenome", "sheep gut metagenome", "sheep metagenome",
            "fish metagenome", "fish gut metagenome", 
            "marine metagenome", "seawater metagenome",
            "freshwater metagenome", "lake water metagenome", "groundwater metagenome",
            "soil metagenome",
            "wastewater metagenome",
            "plant metagenome", "root metagenome", "leaf metagenome"
        ]) & metagenome_assaytype_mask
        metagenomecategory_accs = chunk.loc[metagenomecategory_mask, 'acc'].unique()
        with_metagenomecategory.update(metagenomecategory_accs)

        print(f"Processed chunk {i+1} with {len(chunk)} rows")
    

    # Create summary table
    summary_df = pd.DataFrame({
        "Total_SRA": [total_sra_count],
        "Samples_with_AMR": [len(unique_accessions)],
        "col_date_available": [len(with_dates)],
        "continent_available": [len(with_continent)],
        "date_and_continent_available": [len(with_date_continent)],
        "metagenome_available": [len(with_metagenome)],
        "metagenome_assaytype_available": [len(with_metagenome_assaytype)],
        "metagenome_category_available": [len(with_metagenomecategory)]
    })
    
    # Save summary to CSV
    summary_df.to_csv('../data/summary_statistics_new.csv', index=False)
    
    return summary_df

# Example usage
if __name__ == "__main__":
    # Replace with actual file paths
    sra_path = "../data/SRA_metadata_before20231211.csv"
    card_path = "../data/card_metadata_aro_geolocation.csv"
    
    # Adjust chunk size based on available memory
    result = process_csv_files(sra_path, card_path, chunk_size=30000000)
    print(result)