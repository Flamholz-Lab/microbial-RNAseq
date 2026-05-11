import os
import glob

configfile: "/Users/juliemcdonald/Documents/flamholz_lab/RNASeq_testing/two_part_snakefile/rna_snakefile_part_1_config.yaml"

GENOMES = config["genomes"]
THREADS = config["threads"]


#To run:
# snakemake --snakefile /Users/juliemcdonald/Documents/flamholz_lab/RNASeq_testing/two_part_snakefile/rna_snakefile_part_1.sh --use-conda --jobs 10

rule all:
    input:
        "combined_genomes.fasta",
        "combined_annotation.gff",
        "gene_to_organism_map.tsv"

# ── Merge genomes into one "barcoded" reference ────────────────────────────────────────────
rule merge:
    input: GENOMES
    output: "combined_genomes.fasta"
    threads: THREADS
    conda: "/Users/juliemcdonald/miniconda3/envs/alignment"
    shell:
        """
        for genome_dir in {input}; do
            for fa_file in $genome_dir/*.fa; do
                prefix=$(basename $fa_file .fa)
                sed 's/^>contig/>&'"${{prefix}}"'_/' "$fa_file" > "$fa_file.tmp"
                mv "$fa_file.tmp" "$fa_file"
            done
        done
        cat {input}/*.fa > {output}
        """

# ── Annotate each individual genome with prokka ─────────────────────────────────────────
rule annotate:
    input: GENOMES
    output:
        gff    = "combined_annotation.gff",
        prokka = directory("prokka")
    threads: THREADS
    conda: "/Users/juliemcdonald/miniconda3/envs/prokka"
    shell:
        """
        rm -rf prokka/
        for genome_dir in {input}; do
            for fa_file in $genome_dir/*.fa; do
                prefix=$(basename $fa_file .fa)
                prokka --outdir "prokka/$prefix" --prefix "$prefix" "$fa_file" --cpus {threads}
            done
        done
        cat prokka/*/*.gff > {output.gff}
        """

# ── Create gene-to-organism map ─────────────────────────────
rule source_mapping:
    input:
        genomes = GENOMES,
        prokka  = "prokka"
    output:
        "gene_to_organism_map.tsv"
    threads: THREADS
    conda: "/Users/juliemcdonald/miniconda3/envs/alignment"
    shell:
        """
        for bin in {input.genomes}/*.fa; do
            name=$(basename $bin .fa)
            grep -v "^#" {input.prokka}/$name/$name.gff | \
            awk -v org="$name" '$3=="CDS" {{print $0"\t"org}}'
        done > {output}
        """
