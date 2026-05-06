#!/bin/bash

# Core SNP analysis and Minimum Spanning Tree (MST) for C. psittaci
# Steps: Snippy -> Gubbins -> filter by recombination proportion -> snp-dist -> GrapeTree


# ==============================================
# C. psittaci Metagenomic Analysis Pipeline
# Version: 1.1
# Description: This pipeline processes metagenomic sequencing data
#              for C. psittaci identification and characterization
# Update: 20260506
# ==============================================

# NOTE：pleaes note that this script is an example pipeline for reproducing the analysis. All environment variables, file paths, and software parameters must be adapted to your specific computing environment and data structure. The actual execution may require adjustments based on your local configuration and data characteristics.



# ---------- 1. Run Snippy for all genomes ----------
awk '{print $1"\t./genomes/"$1".fasta"}' sample_list.txt > input.tab
snippy-multi input.tab --ref reference/GCF_000204255.1.fasta --cpus 8 > runme.sh #  Generate batch script for multiple samples
bash runme.sh

# ---------- 2. Generate full core alignment and remove recombination ----------
snippy-clean_full_aln core.full.aln > clean.full.aln # Clean the full alignment by replacing non-standard characters with N
run_gubbins.py -p gubbins clean.full.aln # Identify and filter recombination with Gubbins

# ---------- 3. Extract core SNPs only (polymorphic sites) ----------
snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln 

# ---------- 4. Calculate pairwise SNP distance matrix ----------
snp-dists -c -b clean.core.aln > snp_matrix.csv

# ---------- 5. Build Minimum Spanning Tree (MST) using GrapeTree's MSTreeV2 ----------
grapetree \
    --method MSTreeV2 \
    --distance-matrix snp_matrix.csv \
    --output-tree mst.nwk \
    --output-graphics mst.pdf   # optional

