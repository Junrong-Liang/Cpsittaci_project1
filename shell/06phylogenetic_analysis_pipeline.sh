#!/bin/bash

# ==============================================
# C. psittaci Metagenomic Analysis Pipeline
# Version: 1.0
# Description: This pipeline processes metagenomic sequencing data
#              for C. psittaci identification and characterization
# ==============================================

# NOTE：pleaes note that this script is an example pipeline for reproducing the analysis. All environment variables, file paths, and software parameters must be adapted to your specific computing environment and data structure. The actual execution may require adjustments based on your local configuration and data characteristics.


# ==============================================
# C. psittaci Core Genome Phylogenetic Analysis
# ==============================================

# 1.  Annotating 117 genomes with Prokka
for genome in genomes/*.fasta; do
    sample=$(basename $genome .fasta)
    prokka --outdir annotation/${sample} --prefix ${sample} ${genome}
done

# 2. Step 2: Detecting core genes with panaroo
panaroo \
    -i annotation/*/*.gff \
    -o panaroo_output \
    --mode strict \
    --core_threshold 1.0 \    
    --aligner mafft \
    --remove-invalid \
    --threads 8

# 3. Step 3: Building phylogenetic tree with IQ-TREE..."
iqtree2 \
    -s panaroo_output/core_gene_alignment.aln \
    -m GTR+G+F \
    -bb 1000 \
    -alrt 1000 \
    -nt AUTO \
    -pre cpsittaci_panaroo_tree

# 4. Step 4: Calling core SNPs with Snippy..."
# 4.1 Prepare input file for snippy-multi
# sample_list.txt file format：
# path/to/sample1.fasta
# path/to/sample2.fasta
# path/to/sample3.fasta
awk '{print $1"\t./genomes/"$1".fasta"}' sample_list.txt > input.tab
 
# 4.2 Generate batch script for multiple samples
snippy-multi input.tab --ref reference/GCF_000204255.1.fasta --cpus 8 > runme.sh

# 4.3 Run Snippy analysis for all samples
bash runme.sh

# 4.4 Process core SNP alignment
snippy-clean_full_aln core.full.aln > clean.full.aln # Clean the full alignment by replacing non-standard characters with N

run_gubbins.py -p gubbins clean.full.aln # Identify and filter recombination with Gubbins

snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln # Extract polymorphic sites from filtered alignment

FastTree -gtr -nt clean.core.aln > clean.core.tree # Build phylogenetic tree from core SNPs using FastTree

