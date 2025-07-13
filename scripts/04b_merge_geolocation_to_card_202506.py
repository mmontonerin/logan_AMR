import pandas as pd

# Biosample-geographical-location table 506, new one
geo = (
    pd.read_csv("../data/biosample_geographical_location.202506.csv", dtype=str)
      # discard the null geolocations
      .loc[lambda d: ~d['lat_lon'].isin([
          "0101000020E61000000AD7A3703D4A41C06666666666660CC0",
          "0101000020E6100000EC51B81E85EB0FC052B81E85EB913140"
      ])]
      # rename columns to better names and acc to match card table
      .rename(columns={
          "lat_lon"       : "WTK",
          "attribute_name": "lat_lon",  
          "accession"     : "sample_name"       
      })
      [["sample_name", "lat_lon", "WTK", "biome"]]
)

# CARD metadata table

card = pd.read_csv("../data/full_card_metadata_aro_allfilters_metagenomes_WHOcategories.csv", dtype=str)

# Remove columns from previous geolocation merge
card = card.drop(columns=[c for c in ("lat_lon", "WTK", "biome") if c in card.columns])

# Merge and save

merged = pd.merge(card, geo, on="sample_name", how="left")

merged.to_csv(
    "../data/full_card_metadata_aro_allfilters_metagenomes_WHOcategories_latlon.csv",
    index=False
)
