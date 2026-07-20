'''
-------------------------------------------
Pipeline for:
1. Determining rRNA contamination by species in metatranscriptomic sample (requires a directory of individual species genomes)
2. Aligning RNAseq reads to all species genomes (made as a combined fasta in Part 1)
3. Determining per-gene TPM of all genes

To run:
snakemake --snakefile rna_snakefile_part_2.sh --use-conda --jobs 10
-------------------------------------------
'''

import os
import glob
from pathlib import Path

configfile: "/Users/juliemcdonald/Documents/flamholz_lab/tri_strain_pilot/260616_quintara_RNAseq/260720_rRNA_percent/rna_snakefile_part_2_config.yaml"

GENOMES = config["genomes"]
COMBINED_FASTA = config["combined_fasta"]
COMBINED_ANNO = config["combined_anno"]
RNA_F     = config["rna_F"]
RNA_R     = config["rna_R"]
GTF_PROCESS_PATH = config["gtf_process_path"]
COUNT_PATH = config["count_path"]
THREADS = config["threads"]
MAP = config["map"]
OUT = Path(config["outdir"])


rule all:
    input:
        OUT / "rRNA_stats/rRNA_reference.fa",
        OUT / "rRNA_stats/scafstats.txt",
        OUT / "aligned.bam",
        OUT / "aligned.bam.bai",
        OUT / "featureCounts/counts.txt",
        OUT / "final_gene_counts.csv",
        OUT / "combined.gtf",
        OUT / "aligned_genes.out"

# ── Determine rRNA contamination status: preprocessing ────────────────────────
rule barrnap:
    input:
        ref_genome = GENOMES

    output:
        renamed_genomes = directory(OUT / "rRNA_stats/renamed_genomes"),
        rrna_extracted  = directory(OUT / "rRNA_stats/rrna_extracted"),
        rrna_ref        = OUT / "rRNA_stats/rRNA_reference.fa"

    conda: "/Users/juliemcdonald/miniconda3/envs/barrnap"

    shell:
        """
        mkdir -p {output.renamed_genomes} {output.rrna_extracted}

        for f in {input.ref_genome}/*.fa*; do
            name=$(basename "$f")
            name="${{name%.*}}"
            sed "s/^>/>${{name}}_/" "$f" > {output.renamed_genomes}/${{name}}.fa
        done

        for f in {output.renamed_genomes}/*.fa; do
            name=$(basename "$f" .fa)
            barrnap --kingdom "bac" "$f" --outseq {output.rrna_extracted}/${{name}}_rrna.fa
        done

        cat {output.rrna_extracted}/*_rrna.fa > {output.rrna_ref}
        """

# ── Determine rRNA contamination status: bbduk ─────────────────────────────────
rule bbduk:
    input:
        ref = OUT / "rRNA_stats/rRNA_reference.fa",
        renamed_genomes = directory(OUT / "rRNA_stats/renamed_genomes"),
        rna_f  = RNA_F,
        rna_r  = RNA_R

    output:
        scafstats = OUT / "rRNA_stats/scafstats.txt"
    
    conda: "/Users/juliemcdonald/miniconda3/envs/bbduk"

    shell:
        """
        bbduk.sh in1={input.rna_f} in2={input.rna_r} \
            ref={input.ref} \
            k=31 \
            hdist=1 \
            scafstats={output.scafstats}

        rm -r {input.renamed_genomes}
        """


# ── Align RNAseq reads to combined genome reference ────────────────────────────
rule alignment:
    input:
        genome = COMBINED_FASTA,
        rna_f    = RNA_F,
        rna_r  = RNA_R

    output:
        bam    = OUT / "aligned.bam",
        bai    = OUT / "aligned.bam.bai"

    threads: THREADS
    conda: "/Users/juliemcdonald/miniconda3/envs/alignment"

    shell:
        """
        bwa index {input.genome}
        bwa mem -t {threads} {input.genome} {input.rna_f} {input.rna_r} | \
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
        featureCounts -p -a {input.annotation} -o {output} \
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
