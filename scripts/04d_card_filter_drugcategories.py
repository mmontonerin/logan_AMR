import pandas as pd

# Define the file paths
input_file = "../data/full_card_metadata_aro_allfilters_metagenomes2.csv"
output_file = "../data/full_card_metadata_aro_allfilters_metagenomes_WHOcategories.csv"

# Classify drug classes according to WHO categories
drug_classes = {
    "Access": [
        "penicillin beta-lactam", "tetracycline antibiotic", "nitrofuran antibiotic", 
        "nitroimidazole antibiotic"
    ],
    "Watch": [
        "fluoroquinolone antibiotic", "macrolide antibiotic", "glycopeptide antibiotic", 
        "carbapenem", "lincosamide antibiotic", "monobactam"
    ],
    "Reserve": [
        "glycylcycline", "phosphonic acid antibiotic", "oxazolidinone antibiotic"  # Colistin is handled separately because it appears in AMR_GeneFamily instead of ARO_DrugClass
    ],
    "Mixed": [
        "aminoglycoside antibiotic", "cephalosporin"
    ],
    "Not classified": [
        "isoniazid-like antibiotic"
    ]     
}

# Flatten the drug classes into a dictionary for easy lookup
drug_to_cat = {
    drug.lower(): category
    for category, drugs in drug_classes.items()
    for drug in drugs
}

# Input csv DataFrame
df = pd.read_csv(input_file, low_memory=False)

# Decide order when multiple categories are present
_sort_key = {"Access": 0, "Watch": 1, "Reserve": 2, "Mixed": 3, "Not classified": 4}.get

def assign_who_categories(row):
    cats = set()
    # ----- ARO_DrugClass ----------------------------------------------------
    if isinstance(row["ARO_DrugClass"], str):
        for token in row["ARO_DrugClass"].split(";"):
            token = token.strip().lower()
            cat   = drug_to_cat.get(token)
            if cat:
                cats.add(cat)

    # ----- special case: colistin found in AMR_GeneFamily -------------------
    if isinstance(row["AMR_GeneFamily"], str) and "colistin" in row["AMR_GeneFamily"].lower():
        cats.add("Reserve")
    
    if cats:
        # keep stable WHO order and collapse with “;”
        return ";".join(sorted(cats, key=_sort_key))
    return None

df["WHO_categories"] = df.apply(assign_who_categories, axis=1)
df.to_csv(output_file, index=False)