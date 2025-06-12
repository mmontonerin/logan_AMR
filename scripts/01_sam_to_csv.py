import re

def calculate_percent_identity(de_value):
    """Calculate percent identity from divergence (de field)"""
    if de_value is not None:
        return (1 - de_value) * 100  # Convert divergence to percent identity
    return None  # Return None if de is missing

def get_field_value(fields, tag):
    """Helper function to extract field values based on the tag"""
    for field in fields:
        if field.startswith(f"{tag}:f:"):
            return float(field.split(":")[2])
    return None

def extract_acc_and_contig_id(field_1):
    """Extract the accession and contig_id from the first field"""
    parts = field_1.split('_')
    acc = parts[0]  # Extract accession
    contig_id = '_'.join(parts[:2])  # Extract contig_id
    return acc, contig_id

def extract_aro_id(field3):
    """Extract full ARO ID (ARO:<number>) from field 3 using regex"""
    match = re.search(r"ARO:[0-9]+", field3)
    if match:
        return match.group(0)  # Return the full matched ARO ID
    return None  # Return None if no ARO ID found

def count_total_lines(file_path):
    """Count the total number of non-header lines in the SAM file"""
    with open(file_path, 'r') as f:
        return sum(1 for line in f if not line.startswith('@'))  # Ignore headers

def parse_sam(file_path, output_csv):
    """Parse the SAM file and extract relevant information to CSV with progress updates"""
    total_lines = count_total_lines(file_path)  # Get total number of non-header lines
    progress_intervals = [int(total_lines * i / 10) for i in range(1, 11)]  # 10%, 20%, ... 100%
    
    processed_lines = 0  # Track progress

    with open(file_path, 'r') as f, open(output_csv, 'w') as out_csv:
        # Write CSV header
        out_csv.write("acc,contig_id,ARO_ID,Identity\n")

        for line in f:
            if line.startswith('@'):
                continue  # Skip header lines
            
            fields = line.strip().split("\t")
            
            # Extract accession and contig_id from the first field
            acc, contig_id = extract_acc_and_contig_id(fields[0])
            
            # Extract other fields
            de_value = get_field_value(fields, 'de')

            # Extract ARO ID from field 3
            aro_id = extract_aro_id(fields[2])

            # Skip the line if required fields are missing
            if None in (de_value, aro_id):
                continue
            
            percent_identity = calculate_percent_identity(de_value)
            
            # Filtering condition
            if percent_identity >= 80:
                out_csv.write(f"{acc},{contig_id},{aro_id},{percent_identity:.2f}\n")

            # Update progress
            processed_lines += 1
            if processed_lines in progress_intervals:
                percent_done = (processed_lines / total_lines) * 100
                print(f"Progress: {percent_done:.0f}% completed...", flush=True)

    print("Processing complete. Results saved to", output_csv)

# Path to SAM file and output CSV file
file_path = './data/card_alignment_v1.1_contigs.sam'
output_csv = './data/card_alignment_v1.1_contigs.csv'

# Call the function to parse the SAM file and write the results to CSV
parse_sam(file_path, output_csv)