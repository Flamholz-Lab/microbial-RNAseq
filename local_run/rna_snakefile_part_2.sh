import os
import glob
from pathlib import Path


configfile: "/Users/juliemcdonald/Documents/flamholz_lab/RNASeq_testing/two_part_snakefile/rna_snakefile_part_2_config.yaml"

COMBINED_FASTA = config["combined_fasta"]
COMBINED_ANNO = config["combined_anno"]
RNA     = config["rna"]
GTF_PROCESS_PATH = config["gtf_process_path"]
COUNT_PATH = config["count_path"]
THREADS = config["threads"]
MAP = config["map"]
OUT = Path(config["outdir"])


#To run:
# snakemake --snakefile /Users/juliemcdonald/Documents/flamholz_lab/RNASeq_testing/two_part_snakefile/rna_snakefile_part_2.sh --use-conda --jobs 10

rule all:
    input:
        OUT / "aligned.bam",
        OUT / "aligned.bam.bai",
        OUT / "featureCounts/counts.txt",
        OUT / "final_gene_counts.csv",
        OUT / "combined.gtf",
        OUT / "aligned_genes.out"

# ── Align RNAseq reads to combined genome reference ────────────────────────────
rule alignment:
    input:
        genome = COMBINED_FASTA,
        rna    = RNA
    output:
        bam    = OUT / "aligned.bam",
        bai    = OUT / "aligned.bam.bai"
    threads: THREADS
    conda: "/Users/juliemcdonald/miniconda3/envs/alignment"
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
        annotation = COMBINED_ANNO,
        bam        = OUT / "aligned.bam"
    output: OUT / "featureCounts/counts.txt"
    threads: THREADS
    conda: "/Users/juliemcdonald/miniconda3/envs/featureCounts"
    shell:
        """
        featureCounts -a {input.annotation} -o {output} \
            -t CDS -g ID -T {threads} {input.bam}
        """

# ── GTF processing ────────────────────────────────────
rule GTF_processing:
    input:
        annotation = COMBINED_ANNO,
        path       = GTF_PROCESS_PATH
    output: OUT / "combined.gtf"
    threads: THREADS
    conda: "/Users/juliemcdonald/miniconda3/envs/tpmcalculator"
    shell:
        """
        #Remove extra headers in the prokka GFF
        grep -v "^##sequence-region" {input.annotation} > clean.gff
        
        #Convert GFF to GTF
        
        docker run --rm \
          -v "$PWD":"$PWD" \
          -w "$PWD" \
          quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0 \
          agat_convert_sp_gff2gtf.pl \
            --gff clean.gff \
            --out clean.gtf

        # Replace gene_id with locus_tag in the GTF

        python {input.path} clean.gtf {output}

        #Clean up intermediate files 

        rm clean.gff
        rm clean.gtf
        rm clean.agat.log
        """

# ── Calculate TPMs ────────────────────────────────────
rule tpmcalculator:
    input:
        gtf = OUT / "combined.gtf",
        bam = OUT / "aligned.bam"
    output:
        out = OUT / "aligned_genes.out",
        ent = OUT / "aligned_genes.ent",
        uni = OUT / "aligned_genes.uni"
    threads: THREADS
    conda: "/Users/juliemcdonald/miniconda3/envs/tpmcalculator"
    shell:
        """
        TPMCalculator -g {input.gtf} -b {input.bam}
        mv aligned_genes.out {OUT}/
        mv aligned_genes.ent {OUT}/
        mv aligned_genes.uni {OUT}/
        """

# ── Final table ─────────────────────────────────────────────
rule formatting:
    input:
        tpm   = OUT / "aligned_genes.out",
        mapping  = MAP,
        path     = COUNT_PATH
    output: OUT / "final_gene_counts.csv"
    conda: "/Users/juliemcdonald/miniconda3/envs/alignment"
    shell:
        """
        python {input.path} {input.tpm} {input.mapping}
        mv final_gene_counts.csv {OUT}/
        """