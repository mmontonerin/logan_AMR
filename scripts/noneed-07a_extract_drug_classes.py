import pandas as pd

# Read the CSV file
df = pd.read_csv("./data/card_aro_entries_filtered_metagenomes.csv", low_memory=False)

# Filter out rows where assay_type is "AMPLICON"
df = df[df['assay_type'] != 'AMPLICON']

# Extract only the first ARO_DrugClass from the semicolon-separated list
df['first_aro_drug_class'] = df['ARO_DrugClass'].apply(lambda x: x.split(';')[0] if isinstance(x, str) else "")

# Print drug classes in a separate file
drug_classes = df['first_aro_drug_class'].unique()
with open("./data/drug_classes_first.txt", "w") as f:
    for drug_class in drug_classes:
        f.write(f"{drug_class}\n")

# Extract only the first ARO_DrugClass from the semicolon-separated list
df['last_aro_drug_class'] = df['ARO_DrugClass'].apply(lambda x: x.split(';')[-1] if isinstance(x, str) else "")

# Print drug classes in a separate file
drug_classes = df['last_aro_drug_class'].unique()
with open("./data/drug_classes_last.txt", "w") as f:
    for drug_class in drug_classes:
        f.write(f"{drug_class}\n")

# Extract only the first ARO_DrugClass from the semicolon-separated list
df['all_aro_drug_class'] = df['ARO_DrugClass']

# Print drug classes in a separate file
drug_classes = df['all_aro_drug_class'].unique()
with open("./data/drug_classes_all.txt", "w") as f:
    for drug_class in drug_classes:
        f.write(f"{drug_class}\n")