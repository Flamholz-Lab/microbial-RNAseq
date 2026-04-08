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

The pipeline should run quickly, about 10 min or so. It can easily run on your own device.

## Outputs
A formatted .csv file with the gene name, function, source organism, and RNAseq counts: \
`final_gene_counts.csv/` 

In `example_outputs/`, I ran the snakefile on two different sample E. coli transcriptome experiments and aligned to C. necator, E. coli, and P. aeruginosa genomes. Obviously there were limited to no hits to C. necator or P. aeruginosa. I generated final_counts sheets for each condition, then merged and analyzed them into `merged_gene_counts.csv` using `genome-RNA-scatter.ipynb`. I used plotly in that notebook to generate `wt_vs_lim.html`, which is an interactive map of the nutrient-limited vs. WT transcriptomes.

With real data I'll plot multiple genome-resolved samples at once to look at regulation more closely. 
