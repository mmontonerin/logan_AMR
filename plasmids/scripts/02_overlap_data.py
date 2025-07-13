# 1. Load the tables
import pandas as pd
csv_df = pd.read_csv("../data/plasmids_sra_metadata.csv", dtype=str)       # has column "seq_name"
tsv_df = pd.read_csv("../data/argnorm_output.tsv", sep="\t", dtype=str)  # has column "Contig id"

# 2. Build unique-ID sets
seq_set    = set(csv_df["seq_name"].dropna().unique())
contig_set = set(tsv_df["Contig id"].dropna().unique())

# 3. Quick sanity check / counts
total_seq     = len(seq_set)
total_contig  = len(contig_set)
both_overlap  = len(seq_set & contig_set)
print(f"Plasmids total: {total_seq}")
print(f"AMR-positive: {total_contig}")
print(f"In both: {both_overlap}")

# 4. Venn diagram
from matplotlib import pyplot as plt
from matplotlib_venn import venn2

fig, ax = plt.subplots(figsize=(4,4))

v = venn2([seq_set, contig_set],
          set_labels=("Total plasmids", "AMR-positive plasmids"), ax=ax)

# ── Hide the empty part ───────────────────────────
v.get_patch_by_id('01').set_alpha(0.5)      # area unique to Contig id
v.get_label_by_id('01').set_text('')      # remove the “0” label

plt.tight_layout()
plt.savefig("../data/venn_plasmids_amr.png", dpi=600)
plt.savefig("../data/venn_plasmids_amr.svg")
plt.close()

sizes = [len(contig_set), len(seq_set) - len(contig_set)]
labels = ["AMR-positive plasmids", "Total plasmids"]

fig, ax = plt.subplots(figsize=(4,4), subplot_kw=dict(aspect="equal"))

# Outer ring (total seq_name)
ax.pie([len(seq_set)], radius=1, wedgeprops=dict(width=0.3, edgecolor='w'),
       labels=[f"Total plasmids ({len(seq_set)})"])

# Inner ring (those seq_names that are Contig ids)
ax.pie([len(contig_set)], radius=1-0.3, wedgeprops=dict(width=0.3, edgecolor='w'),
       labels=[f"AMR-positive plasmids ({len(contig_set)})"])

plt.tight_layout()
plt.savefig("../data/ring_plasmids_amr.png", dpi=600)
plt.savefig("../data/ring_plasmids_amr.svg")
plt.close()