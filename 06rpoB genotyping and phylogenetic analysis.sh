#!/bin/bash

# ==============================================
# C. psittaci Metagenomic Analysis Pipeline
# Version: 1.0
# Description: This pipeline processes metagenomic sequencing data
#              for C. psittaci identification and characterization
# ==============================================

# NOTE：pleaes note that this script is an example pipeline for reproducing the analysis. All environment variables, file paths, and software parameters must be adapted to your specific computing environment and data structure. The actual execution may require adjustments based on your local configuration and data characteristics.


# Set parameters
INPUT_DIR="./raw_data"            # Directory containing raw FASTQ files
OUTPUT_DIR="./analysis_results"   # Output directory
REFERENCE_DIR="./references"      # Directory containing reference files
NUM_THREADS=8                     # Number of threads for parallel processing
SAMPLE_ID=$1                      # Sample ID passed as first argument

# Create output directories
mkdir -p ${OUTPUT_DIR}/${SAMPLE_ID}/{quality_control,host_removed,alignment,assembly,virulence,amr,variants,t3ss,phylogeny}
mkdir -p ${OUTPUT_DIR}/logs


# ==============================================
# STEP 10: RPOB CONSENSUS SEQUENCE EXTRACTION
# ==============================================
echo "Step 10: rpoB consensus sequence extraction for ${SAMPLE_ID}"

# Align reads to rpoB reference
bwa mem \
    -t ${NUM_THREADS} \
    ${REFERENCE_DIR}/rpoB_ref.fasta \
    ${OUTPUT_DIR}/${SAMPLE_ID}/alignment/${SAMPLE_ID}_cp_candidate.fastq.gz \
    2>&1 | samtools view -bS - \
    | samtools sort -o ${OUTPUT_DIR}/${SAMPLE_ID}/phylogeny/${SAMPLE_ID}_rpoB_aligned.bam

# Generate consensus sequence
samtools mpileup \
    -A \
    -d 1000 \
    -B \
    -Q 0 \
    -f ${REFERENCE_DIR}/rpoB_ref.fasta \
    ${OUTPUT_DIR}/${SAMPLE_ID}/phylogeny/${SAMPLE_ID}_rpoB_aligned.bam \
    | ivar consensus \
    -p ${OUTPUT_DIR}/${SAMPLE_ID}/phylogeny/${SAMPLE_ID}_rpoB_consensus \
    -n N \
    -m 1 \
    -t 0.5

# Multiple sequence alignment
cat ${REFERENCE_DIR}/rpoB_reference_sequences.fasta \
    ${OUTPUT_DIR}/${SAMPLE_ID}/phylogeny/${SAMPLE_ID}_rpoB_consensus.fa \
    > ${OUTPUT_DIR}/${SAMPLE_ID}/phylogeny/${SAMPLE_ID}_rpoB_all.fasta

mafft \
    --auto \
    --thread ${NUM_THREADS} \
    ${OUTPUT_DIR}/${SAMPLE_ID}/phylogeny/${SAMPLE_ID}_rpoB_all.fasta \
    > ${OUTPUT_DIR}/${SAMPLE_ID}/phylogeny/${SAMPLE_ID}_rpoB_aligned.fasta

# Phylogenetic tree construction
fasttree \
    -nt \
    -gtr \
    ${OUTPUT_DIR}/${SAMPLE_ID}/phylogeny/${SAMPLE_ID}_rpoB_aligned.fasta \
    > ${OUTPUT_DIR}/${SAMPLE_ID}/phylogeny/${SAMPLE_ID}_rpoB_tree.newick