#!/bin/bash

# ==============================================
# C. psittaci Metagenomic Analysis Pipeline
# Version: 1.1
# Description: This pipeline processes metagenomic sequencing data
#              for C. psittaci identification and characterization
# Update: 20260328
# ==============================================

# NOTE：pleaes note that this script is an example pipeline for reproducing the analysis. All environment variables, file paths, and software parameters must be adapted to your specific computing environment and data structure. The actual execution may require adjustments based on your local configuration and data characteristics.


# ==============================================
# C. psittaci Core Genome Phylogenetic Analysis
# Method: Prokka v1.14.6 -> Panaroo v1.2.8 (strict mode, core_threshold=1.0) -> IQ-TREE 2.2.0 with GTR+G+F
# ==============================================

# 1. Annotating 117 genomes with Prokka
for genome in genomes/*.fasta; do
    sample=$(basename $genome .fasta)
    prokka --outdir annotation/${sample} --prefix ${sample} ${genome}
done

# 2. Detecting core genes with panaroo
panaroo \
    -i annotation/*/*.gff \
    -o panaroo_output \
    --mode strict \
    --core_threshold 1.0 \    
    --aligner mafft \
    --remove-invalid \
    --threads 8

# 3. Building phylogenetic tree with IQ-TREE..."
iqtree2 \
    -s panaroo_output/core_gene_alignment.aln \
    -m GTR+G+F \
    -bb 1000 \
    -alrt 1000 \
    -nt AUTO \
    -pre cpsittaci_panaroo_tree
# bb -> UFBoot 1000 times
# -alrt -> SH-aLRT times
