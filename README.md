# microbial-RNAseq
RNA sequencing analysis pipelines

## Inputs
This repository contains a snakemake pipeline that takes a directory of assembled genomes (in .fa format) and RNAseq reads (in .fastq format) as input and runs through genome annotation, RNA normalization, mapping, and counting steps. It also contains a jupyter notebook to visualize outputs (`/genome-RNA-scatter.ipynb`).

To run, you’ll need three files in your working directory or path:

`rna_snakefile` \
`rna_snakefile_config.yaml` \
`rna_snakefile_run.sh` \
`gene_count_integration.py` 

First edit the config file to reflect your directory containing genomes and your RNA .fastq file. The actual snakefile (/rna_snakefile.sh) should not need to be edited. Note that the snakefile utilizes conda environments that I’ve created under my profile in HPC storage. TBD if others in the lab can access those. If not, you will want to make your own conda environments that have the right packages.

To run on the HPC use:

`sbatch rna_snakefile_run.sh`

The pipeline should run quickly, about 10 min or so but scales with the number of genomes you need to annotate. It should be able to run on your own device.

## Outputs
A formatted .csv file with the gene name, function, source organism, and RNAseq counts: \
`final_gene_counts.csv/` 

In `example_outputs/`, I ran the snakefile on two different sample E. coli transcriptome experiments and aligned to C. necator, E. coli, and P. aeruginosa genomes. Obviously there were limited to no hits to C. necator or P. aeruginosa. I generated final_counts sheets for each condition, then merged and analyzed them into `merged_gene_counts.csv` using `genome-RNA-scatter.ipynb`. I used plotly in that notebook to generate `wt_vs_lim.html`, which is an interactive map of the nutrient-limited vs. WT transcriptomes.

With real data I'll plot multiple genome-resolved samples at once to look at regulation more closely. 

## Notes and musings
What count normalization to use for sequencing data? 
Options (https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/):
  1. Raw counts -- how many reads aligned to an mRNA. Shouldn't use this because it doesn't account for gene size.
  2. FPKM/RPKM -- Per million scaling factor for sequencing depth (reads per million), and gene length (reads per million per kilobase
  3. TPM -- transcripts per million. Normalize by gene length first, then sequencing depth so that comparison between samples is easier.
TPM seems to be the gold standard, so I'm using that. So we are measuring a relative abundance of the mRNA to some absolute mRNA content that we don't know (how could we find out?).

We expect relative abundance of a given gene's mRNA to change by some scaling factor (e.g. the slope of condition 1 vs condition 2), and that scaling factor is roughly proportional to the growth rates in condition 1 vs. 2 (Keren et al. 2013). When the relative abundance of a gene's mRNA deviates signficantly (note to self: how do we derive significance here?) from said line, we suspect some biological relevance to nutrient starvation.

I seem to be getting mismatched hits (wrong organism RNA -> wrong organism genome) from very small "genes" that are annotated by Prokka. I wonder if I should set a minimum gene length. Otherwise can go through unexpected hits manually to see if they make sense. 

