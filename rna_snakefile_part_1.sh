import os
import glob

configfile: "/Users/juliemcdonald/Documents/flamholz_lab/tri_strain_pilot/260616_quintara_RNAseq/260714_4_genome_analysis/rna_snakefile_part_1_config.yaml"

GENOMES = glob.glob(os.path.join(config["genomes"], "*.fa"))
THREADS = config["threads"]
print("Genomes found:", GENOMES)


#To run:
# snakemake --use-conda --jobs 10 --snakefile rna_snakefile_part_1.sh 
# snakemake --snakefile rna_snakefile_part_1.sh --use-conda --jobs 10 --allowed-rules source_mapping --dry-run

rule all:
    input:
        "combined_genomes.fasta",
        "combined_annotation.gff3",
        "gene_to_organism_map.tsv"

# ── Merge genomes into one "barcoded" reference ────────────────────────────────────────────
rule merge:
    input: GENOMES
    output: "combined_genomes.fasta"
    threads: THREADS
    conda: "/Users/juliemcdonald/miniconda3/envs/bakta"
    shell:
        """
        for fa in {input}; do
            SAMPLE=$(basename "${{fa%.fa}}")
            echo "Processing: $fa"
            echo "Sample name: ${{SAMPLE}}"
            
            sed "s/^>/>${{SAMPLE}}_/" "$fa" > "${{SAMPLE}}_renamed.fa"
            echo "Renamed file size: $(wc -l < ${{SAMPLE}}_renamed.fa) lines"
        done

        echo "Renamed files found:"
        ls *_renamed.fa

        cat *_renamed.fa > {output}
        echo "Output size: $(wc -l < {output}) lines"
        """


# ── Annotate combined genome with Bakta ────────────────────────────────────────────
rule bakta:
    input: 
        fasta = "combined_genomes.fasta",
        db = "/Volumes/One_Touch/databases/bakta_db"
    output: "combined_annotation.gff3"
    threads: THREADS
    conda: "/Users/juliemcdonald/miniconda3/envs/bakta"
    shell:
    # can't seem to get this working, not sure why. 
    # ERROR: tRNAscan-SE could not be executed! Please make sure tRNAscan-SE is installed and executable or skip requiring workflow steps via via '--skip-trna'.
        """
        export PERL5LIB="${{PERL5LIB:-}}"
        bakta {input.fasta} \
            --db {input.db} \
            --keep-contig-headers \
            --output bakta \
            --threads {threads}

        cp bakta/combined.gff {output}
        """

# ── Map each gene to its source organism ────────────────────────────────────────────
rule source_mapping:
    input: "combined_annotation.gff3"
    output: "gene_to_organism_map.tsv"
    run:
        with open(input[0]) as f, open(output[0], "w") as out:
            for line in f:
                # skip comment lines and FASTA section
                if line.startswith("#"):
                    continue
                if line.startswith(">"):
                    break
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                contig = parts[0]          # e.g. ec_contig_10
                organism = contig.split("_")[0]  # e.g. ec
                out.write("\t".join(parts) + "\t" + organism + "\n")