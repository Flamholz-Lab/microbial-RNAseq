import pandas as pd
import sys

#inputs
tpm_file  = sys.argv[1]
ann_file  = sys.argv[2]

# Read the gene counts table
tpm_file = pd.read_csv(tpm_file, sep=r"\s+", engine="python")
tpm_file.columns = ["gene_id", "Chr", "Start", "End", "Length", "Reads", "TPM", 
                  "ExonLength", "ExonReads", "ExonTPM", "IntronLength", 
                  "IntronReads", "IntronTPM", "UniqueLength", "UniqueReads", 
                  "UniqueTPM", "UniqueExonLength", "UniqueExonReads", 
                  "UniqueExonTPM", "UniqueIntronLength", "UniqueIntronReads", 
                  "UniqueIntronTPM"]

# Read annotation file - ONLY columns 0-8 (ignore extra columns)
ann_cols = [0,1,2,3,4,5,6,7,8]  # First 9 columns only
ann = pd.read_csv(ann_file, sep=r"\s+", usecols=ann_cols, 
                  header=None, engine="python")

# Read raw lines, parse manually
with open(ann_file) as f:
    lines = [line.strip().split("\t", 8) for line in f if not line.startswith("#")]

ann = pd.DataFrame(lines, columns=[f"col{i}" for i in range(9)])
ann["gene_id"] = ann["col8"].str.extract(r'locus_tag=([^;]+)', expand=False)
ann["Product"] = ann["col8"].str.extract(r'product=([^;\t]+)', expand=False)
ann["Name"] = ann["col8"].str.extract(r'Name=([^;]+)', expand=False)


# Clean Name field
ann["Name"] = ann["Name"].fillna("").str.strip()

# Keep only rows with gene_id and drop duplicates
ann_small = ann[["gene_id", "Name", "Product"]].dropna(subset=["gene_id"])

# Merge counts with annotation on gene_id
merged = tpm_file.merge(ann_small, on="gene_id", how="left")

# Select and reorder final columns
final = merged[["gene_id", "Chr", "Reads", "TPM", "Name", "Product"]]

# Write output
final.to_csv("final_gene_counts.csv", index=False)
