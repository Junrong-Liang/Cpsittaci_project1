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
# STEP 1: QUALITY CONTROL WITH FASTP
# ==============================================
echo "Step 1: Quality control for ${SAMPLE_ID}"
fastp \
    --in1 ${INPUT_DIR}/${SAMPLE_ID}_R1.fastq.gz \
    --out1 ${OUTPUT_DIR}/${SAMPLE_ID}/quality_control/${SAMPLE_ID}_clean_R1.fastq.gz \
    --json ${OUTPUT_DIR}/${SAMPLE_ID}/quality_control/${SAMPLE_ID}_fastp.json \
    --html ${OUTPUT_DIR}/${SAMPLE_ID}/quality_control/${SAMPLE_ID}_fastp.html \
    --cut_mean_quality 20 \
    --n_base_limit 5 \
    --length_required 34 \
    --cut_tail \
    --dont_eval_duplication \
    --trim_tail1 1 \
    --disable_trim_poly_g \
    --low_complexity_filter \
    --complexity_threshold 7 \
    --thread ${NUM_THREADS} \
    --detect_adapter_for_pe \
    --correction \
    --overrepresentation_analysis \
    2>&1 | tee ${OUTPUT_DIR}/logs/${SAMPLE_ID}_step1_fastp.log
