import os
import glob

configfile: "rna_snakefile_config.yaml"

GENOMES = config["genomes"]
RNA     = config["rna"]
THREADS = config["threads"]

rule all:
    input:
        "combined_genomes.fasta",
        "combined_annotation.gff",
        "aligned.bam",
        "aligned.bam.bai",
        "featureCounts/counts.txt",
        "gene_to_organism_map.tsv",
        "final_gene_counts.csv"

# ── Merge genomes into one "barcoded" reference ────────────────────────────────────────────
rule merge:
    input: GENOMES
    output: "combined_genomes.fasta"
    threads: THREADS
    conda: "/ru-auth/local/home/jmcdonald/miniconda3/envs/alignment"
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
    output: "combined_annotation.gff"
    threads: THREADS
    conda: "/ru-auth/local/home/jmcdonald/miniconda3/envs/prokka"
    shell:
        """
        rm -rf prokka/
        for genome_dir in {input}; do
            for fa_file in $genome_dir/*.fa; do
                prefix=$(basename $fa_file .fa)
                prokka --outdir "prokka/$prefix" --prefix "$prefix" "$fa_file"
            done
        done
        cat prokka/*/*.gff > {output}
        """

# ── Align RNAseq reads to combined genome reference ────────────────────────────
#do you need to deal with gene length here?
rule alignment:
    input:
        genome = "combined_genomes.fasta",
        rna    = RNA
    output:
        bam    = "aligned.bam",
        bai    = "aligned.bam.bai"
    threads: THREADS
    conda: "/ru-auth/local/home/jmcdonald/miniconda3/envs/alignment"
    shell:
        """
        bwa index {input.genome}
        bwa mem -t {threads} {input.genome} {input.rna} | \
        samtools view -b -F 4 | samtools sort -o {output.bam}
        samtools index {output.bam}
        """

# ── Count reads per gene ────────────────────────────────────
rule featureCounts:
    input:
        annotation = "combined_annotation.gff",
        bam        = "aligned.bam"
    output:
        "featureCounts/counts.txt"
    threads: THREADS
    conda: "/ru-auth/local/home/jmcdonald/miniconda3/envs/featureCounts"
    shell:
        """
        featureCounts -a {input.annotation} -o {output} \
            -t CDS -g ID -T {threads} {input.bam}
        """


#add tpm calculation step


# ── Create gene-to-organism map ─────────────────────────────
rule source_mapping:
    input:
        genomes = GENOMES,
        #annot   = "combined_annotation.gff"  # Ensures annotate runs first
    output: 
        "gene_to_organism_map.tsv"
    threads: THREADS
    conda: "/ru-auth/local/home/jmcdonald/miniconda3/envs/alignment"
    shell:
        """
        for bin in {input.genomes}/*.fa; do
            name=$(basename $bin .fa)
            grep -v "^#" prokka/$name/$name.gff | \
            awk -v org="$name" '$3=="CDS" {{print $0"\t"org}}'
        done > {output}
        """

# ── Final table ─────────────────────────────────────────────
rule formatting:
    input:
        counts   = "featureCounts/counts.txt",
        mapping  = "gene_to_organism_map.tsv",
        bam      = "aligned.bam"
    output: "final_gene_counts.csv"
    conda: "/ru-auth/local/home/jmcdonald/miniconda3/envs/alignment"
    shell:
        """
        python /lustre/fs4/flam_lab/scratch/jmcdonald/scripts/gene_count_integration.py {input.counts} {input.mapping} {input.bam}
        """