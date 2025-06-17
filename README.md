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
python -u 02_merge_metadata2.py > 20250304-merge.log 2>&1 &

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
Using the script 04a_merge_aro_card.py
It also contains a progress tracker. Used an ec2 r7a.2xlarge instance.

```
python -u 04a_merge_aro_card.py > merge_arocard.log 2>&1 &

# Keep track of progress:
tail -f merge_arocard.log
```

## Merge updated geolocation data
Using the script 04b_merge_new_geolocation.py
It also contains a progress tracker. Used an ec2 r7a.2xlarge instance.

```
python -u 04b_merge_new_geolocation.py > merge_geolocation.log 2>&1 &

# Keep track of progress:
tail -f merge_arocard.log
```

## Table filters

#### Filter out the data not containing continent nor date of sampling

```
['collection_date_sam'].notna()
['geo_loc_name_country_continent_calc'].notna()
['collection_date_sam'] != 'uncalculated'
['geo_loc_name_country_continent_calc'] != 'uncalculated'
```

#### Keep only metagenome data

```
['organism'].notna()
['organism'].str.contains("metagenome")
```

#### Only keep assay types that contain most of the genome/gene sequences

```
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

['assay_type'].notna()
['assay_type'].isin(assay_types_to_keep)
```

#### Subcluster metagenome categories
```
categories = {
    "human": [
        "human gut metagenome", "human metagenome", "human oral metagenome", "human skin metagenome",
        "human feces metagenome", "human vaginal metagenome", "human nasopharyngeal metagenome", 
        "human lung metagenome", "human saliva metagenome", "human reproductive system metagenome", 
        "human urinary tract metagenome", "human eye metagenome", "human blood metagenome", 
        "human bile metagenome", "human tracheal metagenome", "human brain metagenome", 
        "human milk metagenome", "human semen metagenome", "human skeleton metagenome"
    ],
    "livestock": [
        "bovine gut metagenome", "bovine metagenome", "pig gut metagenome", "pig metagenome", 
        "chicken gut metagenome", "chicken metagenome", "sheep gut metagenome", "sheep metagenome"
    ],
    "fish": ["fish metagenome", "fish gut metagenome"],
    "marine": ["marine metagenome", "seawater metagenome"],
    "freshwater": ["freshwater metagenome", "lake water metagenome", "groundwater metagenome"],
    "soil": ["soil metagenome"],
    "wastewater": ["wastewater metagenome"],
    "plant": ["plant metagenome", "root metagenome", "leaf metagenome"]
} 
```


## Card-ARO-metadata filter
First overall filter, including date, location, metagenome only, assay_type selection

```
python -u 04c_cardaro_filter.py > cardaro_filter.log 2>&1 &

# Keep track of progress:
tail -f cardaro_filter.log
```

Add cluster classification of metagenome categories
```
python 04d_card_filter_metagenome_categories.py &
```

## SRA metadata filter
Filter out data post 11 December 2023 (the date in which Logan was run), inclusive.
```
python 05a_SRA_metadata_filter_post_11Dec2023.py &
```

Filter out data not containing date and continent. Keep only metagenomes, and selected assay_types
```
python 05b_SRA_metadata_filter_date_and_continent.py &
python 05c_SRA_metadata_filter_metagenomes.py &
python 05d_SRA_metadata_filter_assaytype.py &
```

Add cluster classification of metagenome categories
```
python 05e_SRA_metadata_filter_metagenomecategory.py &
```

# Data accessibility:
### SRA metadata
```
# SRA metadata post 2023-12-11 inclusive
s3://serratus-merce/AMR/SRA_metadata_before20231211.csv.zst 
# SRA metadata filtered
s3://serratus-merce/AMR/SRA_metadata_before20231211_date_and_continent_metagenomes_assaytype.csv.zst
# SRA metadata filtered and metagenomic clusters only
s3://serratus-merce/AMR/SRA_metadata_allfilters.csv.zst  
```
### Card-ARO-metadata table
```
# Complete table of results, including ARO information, and extended geolocation
s3://serratus-merce/AMR/card_metadata_aro_geolocation.csv.zst
# Filtered table
s3://serratus-merce/AMR/full_card_metadata_aro_allfilters.csv.zst
# Filtered table and metagenomic clusters only
s3://serratus-merce/AMR/full_card_metadata_aro_allfilters_metagenomes2.csv.zst
```

