import pandas as pd
import sys

#inputs
count_file  = sys.argv[1]
map_file  = sys.argv[2]
count_bam_name = sys.argv[3]

# Load featureCounts output
counts = pd.read_csv(count_file, sep="\t", comment="#")

# Load gene → organism mapping
mapping = pd.read_csv(map_file, sep="\t",
                      names=["contig","source","feature","start","end",
                             "score","strand","frame","attributes","organism"])
mapping["gene_id"] = mapping["attributes"].str.extract(r'ID=([^;]+)')

# Load ALL Prokka .tsv annotation files and combine
annotation_frames = []
import glob

for tsv in glob.glob("prokka/*/*.tsv"):
    df = pd.read_csv(tsv, sep="\t")
    annotation_frames.append(df)

prokka_annotations = pd.concat(annotation_frames, ignore_index=True)
# Prokka .tsv columns: locus_tag, ftype, length_bp, gene, EC_number, COG, product
prokka_annotations = prokka_annotations[["locus_tag", "gene", "product"]]

# Merge everything together
final = counts.merge(mapping[["gene_id", "organism"]],
                     left_on="Geneid", right_on="gene_id")

final = final.merge(prokka_annotations,
                    left_on="Geneid", right_on="locus_tag",
                    how="left")  # left join so genes with no annotation aren't dropped

# Select and rename final columns
final = final[["Geneid", "gene", "product", "organism", count_bam_name]]
final.columns = ["gene_id", "gene_name", "product", "organism", "count"]

# Optional: flag unannotated genes clearly
final["product"] = final["product"].fillna("hypothetical protein")
final["gene_name"] = final["gene_name"].fillna("unknown")

final.to_csv("final_gene_counts.csv", index=False)