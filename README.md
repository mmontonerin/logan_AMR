# Logan AMR

## Installations needed
```
awscli (https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)
zstd (https://github.com/facebook/zstd)
```

## Data needed
#### CARD alignment
Alignment of CARD nucleotide database to Logan v1.1 contigs.

Downloaded using:
```
aws s3 cp s3://serratus-rayan/beetles/logan_feb27_run/minimap2-concat/card_alignment_v1.1_contigs.sam.zst .
zstd -d card_alignment_v1.1_contigs.sam.zst
```
Compressed: 12.6 Gb
Uncompressed: 85 Gb

#### SRA Metadata
General metadata provided by Kristen (SRA_metadata.csv)
Additional geolocation data provided by Alex (bgl_gm4326_gp4326.csv)
Downloaded using:
```
aws s3 cp s3://serratus-rayan/beetles/SRA_metadata.csv.zst .
zstd -d SRA_metadata.csv.zst

wget https://serratus-public.s3.us-east-1.amazonaws.com/tmp/bgl_gm4326_gp4326.csv.gz
gunzip bgl_gm4326_gp4326.csv.gz
```
Compressed SRA metadata: 2 Gb
Uncompressed SRA metadata: 9.6 Gb

Uncompressed geolocation: 7.9 Gb

#### CARD metadata
Downloaded from the CARD database: https://card.mcmaster.ca/download/
Downloaded using:
```
wget https://card.mcmaster.ca/download/0/broadstreet-v3.3.0.tar.bz2
tar -xf broadstreet-v3.3.0.tar.bz2
```

I used the script 01_sam_to_csv.py to create a CSV table from the SAM CARD alignment file to Logan v1.1, containing accession_number, contig_id, ARO_id, alignment_length, %identity
I also selected out those alignments of less than 100 bp and of lower alignment identity threshold than 80%.


## Merge metadata with card-alignment CSV table

I used the script 02_merge_metadata2.py to do that. I also keeps track of progress.
Note, minimum of 45 GB RAM memory needed to store SRA_metadata.csv into memory as dictionary. 
Used an ec2 r7a.2xlarge instance.

```
python 02_merge_metadata2.py > 20250304-merge.log 2>&1 &

# Keep track of progress:
tail -f 20250304-merge.log
```

There are SRA accessions that are present in the CARD alignment file, but missing in the SRA metadata file, because they might have been removed from the databases. 
I used the script 03_evaluate_missing_accessions.py to filter out those with missing metadata, and also get an idea of how many are we filtering:
```
Rows with ALL metadata missing: 828884
Unique accessions with ALL metadata missing: 7463
```

## Merge ARO data with card-alignment-metadata CSV table
I used the script 04_merge_aro_card.py



