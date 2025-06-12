import pandas as pd

# Define the file paths
input_file = "./data/card_aro_entries_filtered_metagenomes.csv"
output_file = "./data/card_aro_entries_metagenomes_and_drugclass.csv"

df = pd.read_csv(input_file, low_memory=False)

# Define drug categories
action = {
    "cell_wall_synthesis_inhibitors": {
        "penicillin beta-lactam", "cephalosporin", "carbapenem", "glycopeptide antibiotic", "cycloserine-like antibiotic"
    },
    "protein_synthesis_inhibitors": {
        "aminoglycoside antibiotic", "tetracycline antibiotic", "macrolide antibiotic", "lincosamide antibiotic", 
        "phenicol antibiotic", "streptogramin A antibiotic", "streptogramin B antibiotic", "oxazolidinone antibiotic", 
        "elfamycin antibiotic", "pleuromutilin antibiotic", "glycylcycline"
    },
    "dna_rna_synthesis_inhibitors": {
        "fluoroquinolone antibiotic", "rifamycin antibiotic", "diarylquinoline antibiotic", "bicyclomycin-like antibiotic"
    },
    "folic_acid_synthesis_inhibitors": {
        "sulfonamide antibiotic", "diaminopyrimidine antibiotic"
    },
    "membrane_disruptors": {
        "antibacterial free fatty acids", "disinfecting agents and antiseptics"
    },
    "other_mechanisms": {
        "phosphonic acid antibiotic", "pyrazine antibiotic", "nitrofuran antibiotic", "nitroimidazole antibiotic", 
        "mupirocin-like antibiotic", "fusidane antibiotic", "aminocoumarin antibiotic", "isoniazid-like antibiotic", 
        "salicylic acid antibiotic", "polyamine antibiotic", "nucleoside antibiotic"
    }
}

activity = {
    "narrow_spectrum": {
        "penicillin beta-lactam", "vancomycin", "rifamycin antibiotic", "isoniazid-like antibiotic", 
        "fusidane antibiotic", "cycloserine-like antibiotic"
    },
    "broad_spectrum": {
        "tetracycline antibiotic", "fluoroquinolone antibiotic", "macrolide antibiotic", "cephalosporin", 
        "carbapenem", "sulfonamide antibiotic", "phenicol antibiotic", "oxazolidinone antibiotic", 
        "streptogramin A antibiotic", "streptogramin B antibiotic"
    }
}

structure = {
    "beta_lactams": {
        "penicillin beta-lactam", "cephalosporin", "carbapenem"
    },
    "aminoglycosides": {
        "aminoglycoside antibiotic"
    },
    "tetracyclines": {
        "tetracycline antibiotic", "glycylcycline"
    },
    "macrolides": {
        "macrolide antibiotic"
    },
    "lincosamides": {
        "lincosamide antibiotic"
    },
    "phenicols": {
        "phenicol antibiotic"
    },
    "fluoroquinolones": {
        "fluoroquinolone antibiotic"
    },
    "rifamycins": {
        "rifamycin antibiotic"
    },
    "sulfonamides": {
        "sulfonamide antibiotic"
    },
    "diaminopyrimidines": {
        "diaminopyrimidine antibiotic"
    },
    "nitrofurans": {
        "nitrofuran antibiotic"
    },
    "nitroimidazoles": {
        "nitroimidazole antibiotic"
    },
    "glycopeptides": {
        "glycopeptide antibiotic"
    },
    "streptogramins": {
        "streptogramin A antibiotic", "streptogramin B antibiotic"
    },
    "pleuromutilins": {
        "pleuromutilin antibiotic"
    },
    "oxazolidinones": {
        "oxazolidinone antibiotic"
    },
    "aminocoumarins": {
        "aminocoumarin antibiotic"
    },
    "polyamines": {
        "polyamine antibiotic"
    },
    "nucleosides": {
        "nucleoside antibiotic"
    },
    "pyrazines": {
        "pyrazine antibiotic"
    },
    "phosphonic_acids": {
        "phosphonic acid antibiotic"
    },
    "bicyclomycins": {
        "bicyclomycin-like antibiotic"
    },
    "elfamycins": {
        "elfamycin antibiotic"
    },
    "fusidanes": {
        "fusidane antibiotic"
    },
    "isoniazid_like": {
        "isoniazid-like antibiotic"
    },
    "cycloserine_like": {
        "cycloserine-like antibiotic"
    },
    "salicylic_acids": {
        "salicylic acid antibiotic"
    },
    "free_fatty_acids": {
        "antibacterial free fatty acids"
    },
    "disinfectants": {
        "disinfecting agents and antiseptics"
    }
}

# Extract only the first ARO_DrugClass from the semicolon-separated list
df['first_aro_drug_class'] = df['ARO_DrugClass'].apply(lambda x: x.split(';')[0] if isinstance(x, str) else "")

# Function to map first_aro_drug_class to drug category
def get_action(drug):
    for category, values in action.items():
        if drug in values:
            return category
    return None

def get_activity(drug):
    for category, values in activity.items():
        if drug in values:
            return category
    return None

def get_structure(drug):
    for category, values in structure.items():
        if drug in values:
            return category
    return None

df['drug_action'] = df['first_aro_drug_class'].apply(get_action)
df['drug_activity'] = df['first_aro_drug_class'].apply(get_activity)
df['drug_structure'] = df['first_aro_drug_class'].apply(get_structure)

# Save the new table
df.to_csv(output_file, index=False)