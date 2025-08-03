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


## 2. Merge metadata with card-alignment CSV table

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

## 4. Merge ARO data with card-alignment-metadata CSV table
Using the script 04_merge_aro_card.py
It also contains a progress tracker. Used an ec2 r7a.2xlarge instance.

```
python -u 04_merge_aro_card.py > merge_arocard.log 2>&1 &

# Keep track of progress:
tail -f merge_arocard.log
```

## 5. Extended informative columns on organism type and metagenome category

Add informative columns to both SRA and AMR detected tables
```
# In AMR table
# Create full extended table
python 05b_informative_columns_AMRtotal_extended.py

# Create minimal table to limit size of table for some specific plots:
python 05b_informative_columns_AMRtotal_minimal.py

# In SRA table
python 05c_informative_columns_SRAtotal.py
```

It creates two new columns:
* metagenome_category, which classifies metagenome samples into the more abundant groups
```
# Found in column organism
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
    "marine": ["marine metagenome", "seawater metagenome"],
    "freshwater": ["freshwater metagenome", "lake water metagenome", "groundwater metagenome"],
    "soil": ["soil metagenome"],
    "wastewater": ["wastewater metagenome"]
}

# Found in column librarysource, when organism is named after the host:
categories_sp = {
    "livestock": [
        "Sus scrofa", "Sus scrofa domesticus", "Sus scrofa affinis", "Bos taurus", "Gallus gallus", "Equus caballus", 
        "Equs caballus", "Ovis aries", "Ovis", "Bos indicus", "Bos mutus", "Bos primigenius", "Bos frontalis",
        "Bos gaurus", "Gallus", "Capra hircus", "Capra aegagrus", "Capra ibex"
    ],
    "human": [
        "Homo sapiens"
    ]  
}
```

* organism_type, which identifies samples as Metagenome or Isolate
```
# Isolate to be classified whenever column organism does not contain the word metagenome or column library source does not contain METAGENOMIC / TRANSCRIPTOMIC
```

## 6. Table filters for specific plots

```
# AMR:
python 06a_cardaro_filter_metaloc_plots.py
# SRA:
python 06b_sra_filter_metaloc_plots.py
```

We also add the updated geolocation
```
python 06c_merge_geolocation_to_card_202506.py
```

Filters:
```
# Filter out the data not containing continent nor date of sampling
['collection_date_sam'].notna()
['geo_loc_name_country_continent_calc'].notna()
['collection_date_sam'] != 'uncalculated'
['geo_loc_name_country_continent_calc'] != 'uncalculated'

# Keep only metagenome data
['organism_type'] != 'Isolate'

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

['assay_type'].notna()
['assay_type'].isin(assay_types_to_keep)
```

## Generate plots

Square plot with density of 
* AMR+/AMR- SRA accessions
* AMR+ isolate/metagenome (organism type)
* AMR+ metagenome category

```
python 07_squarify_plot.py
```

log2 Enrichment of AMR+ vs AMR- plasmids found in metagenome categories
Data is randomly subsampled to avoid biases caused by differences in data count

```
python 08_log_enrichment_normalize.py
```

Map plot
```
python 09_map_plot_bubblesize.py
```

Timeline of discovery
```
python 10_rateofdiscovery.py
```

Number of AMR genes detected per accession (in metagenome samples)
Boxplot:
```
python 11_boxplot_amrperaccession.py
```
log2 enrichment:
```
python 11_log_amrperaccession.py
```

# Plasmids

Data used for this analysis was provided by Antonio Camargo, from his pipeline of plasmid discovery from logan contigs.
The AMR detection in plasmids was done with AMRFinderPlus.

### Initial data: 
```
# All AMR genes found in plasmids
s3://serratus-merce/AMR/argnorm_output.tsv.zst
# All plasmids
s3://serratus-merce/AMR/plasmid_data.tsv.zst 
```

### Pipeline
Scripts found in all_scripts_forAMRfigure/plasmids

#### 1. Create tables

Create table of plasmids with selected SRA metadata
```
python 01_create_totalplasmids_table1.py
```

Add SRA metadata to AMR detection table
```
python 01b_amrpositive_table.py
```

Add extended data (geolocation, as well as organism type, and metagenome category columns)
```
# On AMR positive table
python 01c_metagenomecategory_amr_extended.py

# On full plasmids table
python 01c_metagenomecategory_extended.py
```

#### 2. Generate plots

Ring plot with density of 
* AMR+/AMR- plasmids
* AMR+ isolate/metagenome (organism type)
* AMR+ metagenome category
  
```
python 02_ring_plot_plasmids.py
```

Enrichment of AMR+ vs AMR- plasmids found in metagenome categories
Data is randomly subsampled to avoid biases caused by differences in data count

```
python 03_log_enrichment_normalized.py
```

# Data accessibility:
### SRA metadata
```
# SRA metadata post 2023-12-11 inclusive, only data present in logan
s3://serratus-merce/AMR/AMR/SRA_metadata_before20231211_logan.csv.zst
# SRA metadata extended columns
s3://serratus-merce/AMR/AMR/SRA_metadata_before20231211_logan_extended.csv.zst
# SRA metadata filtered and metagenomic only
s3://serratus-merce/AMR/AMR/SRA_metadata_before20231211_logan_dateloc_meta.csv.zst 
```
### Card-ARO-metadata table
```
# Complete AMR table with extended columns
s3://serratus-merce/AMR/AMR/card_metadata_aro_extended.csv.zst
# Only informative columns needed for specific plots
s3://serratus-merce/AMR/AMR/card_metadata_aro_informativecolumns_minimal.csv.zst
# Filtered table, only metagenomics
s3://serratus-merce/AMR/AMR/card_metadata_aro_dateloc_meta.csv.zst
# Filtered table, only metagenomics, and new geolocation
s3://serratus-merce/AMR/AMR/card_metadata_aro_dateloc_meta_latlon.csv.zst
```
### Plasmids datasets
```
# Plasmids data
s3://serratus-merce/AMR/plasmids/plasmid_data.tsv.zst
# AMR detected in plasmids table
s3://serratus-merce/AMR/plasmids/argnorm_output.tsv.zst
# Plasmids data extended columns 
s3://serratus-merce/AMR/plasmids/plasmids_sra_metadata_extended.csv.zst
# AMR+ plasmids extended columns
s3://serratus-merce/AMR/plasmids/amr_metadata_extended.csv.zst
``` 






