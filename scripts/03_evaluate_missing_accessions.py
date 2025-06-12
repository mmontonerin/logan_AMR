import csv
import sys

# List of metadata columns
metadata_columns = [
    'assay_type', 'center_name', 'consent', 'experiment', 'sample_name', 
    'instrument', 'librarylayout', 'libraryselection', 'librarysource', 
    'platform', 'sample_acc', 'biosample', 'organism', 'sra_study', 
    'releasedate', 'bioproject', 'mbytes', 'loaddate', 'avgspotlen', 
    'mbases', 'insertsize', 'library_name', 'biosamplemodel_sam', 
    'collection_date_sam', 'geo_loc_name_country_calc', 
    'geo_loc_name_country_continent_calc'
]

def is_line_completely_missing(line, metadata_column_indices):
    """
    Check if all metadata columns are empty/missing for a given line
    """
    # Check if all specified columns are empty or 'NA'
    return all(
        line[idx].strip() in ['', 'NA', 'None'] 
        for idx in metadata_column_indices
    )

def filter_metadata(input_filepath, output_filepath):
    # Tracking variables
    total_rows = 0
    rows_with_data = 0
    completely_missing_rows = 0
    completely_missing_accessions = set()

    # Open input and output files
    with open(input_filepath, 'r') as csvfile_in, \
         open(output_filepath, 'w', newline='') as csvfile_out, \
         open('./data/card_missing_metadata.csv', 'w', newline='') as missing_file:
        
        # Use csv reader and writers for robust parsing
        reader = csv.reader(csvfile_in)
        writer_out = csv.writer(csvfile_out)
        writer_missing = csv.writer(missing_file)
        
        # Read and write header
        header = next(reader)
        writer_out.writerow(header)
        writer_missing.writerow(header)
        
        # Map metadata column names to their indices
        try:
            metadata_column_indices = [
                header.index(col) for col in metadata_columns
            ]
            acc_index = header.index('acc')
        except ValueError as e:
            print(f"Error: Column not found in CSV. {e}")
            sys.exit(1)

        # Process each line
        for line in reader:
            total_rows += 1

            # Check for missing metadata
            is_completely_missing = is_line_completely_missing(line, metadata_column_indices)
            
            if is_completely_missing:
                completely_missing_rows += 1
                completely_missing_accessions.add(line[acc_index])
                writer_missing.writerow(line)
            else:
                # Write rows with some metadata to output file
                writer_out.writerow(line)
                rows_with_data += 1

    # Print results
    print(f"Total rows in original dataset: {total_rows}")
    print(f"Rows with metadata: {rows_with_data}")
    print(f"Rows with ALL metadata missing: {completely_missing_rows}")
    print(f"Unique accessions with ALL metadata missing: {len(completely_missing_accessions)}")
    print(f"\nFiltered metadata saved to: {output_filepath}")
    print(f"Completely missing rows saved to: ./data/card_missing_metadata.csv")

# Run the filtering
filter_metadata(
    "./data/card_metadata.csv", 
    "./data/card_metadata_filtered.csv"
)