import pandas as pd

# Read the CSV file
df = pd.read_csv("./data/card_metadata_aro_geolocation_metagenomes_onlydateloc_noamplicon.csv", low_memory=False)

# Count occurrences of each ARO_DrugClass
drug_class_counts = df['ARO_DrugClass'].dropna().str.split(';').explode().value_counts()

# Save to CSV
drug_class_counts.to_csv("./data/forTiscar_drug_class_counts.csv", header=["Count"])
