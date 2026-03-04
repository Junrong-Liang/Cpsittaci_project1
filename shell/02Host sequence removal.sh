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
# STEP 2: HOST SEQUENCE REMOVAL
# ==============================================
echo "Step 2: Host sequence removal for ${SAMPLE_ID}"

# Step 2.1: Kraken2 classification
echo "Step 2.1: Kraken2 classification"
kraken2 \
    --db /path/to/kraken2_standard_db \
    --threads ${NUM_THREADS} \
    --unclassified-out ${OUTPUT_DIR}/${SAMPLE_ID}/host_removed/${SAMPLE_ID}_unclassified#.fastq \
    --classified-out ${OUTPUT_DIR}/${SAMPLE_ID}/host_removed/${SAMPLE_ID}_classified#.fastq \
    --report ${OUTPUT_DIR}/${SAMPLE_ID}/host_removed/${SAMPLE_ID}_kraken2_report.txt \
    --output ${OUTPUT_DIR}/${SAMPLE_ID}/host_removed/${SAMPLE_ID}_kraken2_output.txt \
    --use-names \
    --gzip-compressed \
    ${OUTPUT_DIR}/${SAMPLE_ID}/quality_control/${SAMPLE_ID}_clean_R1.fastq.gz \
    2>&1 | tee ${OUTPUT_DIR}/logs/${SAMPLE_ID}_step2_1_kraken2.log

# Step 2.2: BWA alignment to hg38 host genome
echo "Step 2.2: BWA alignment to hg38 host genome"
# Index human reference genome if not already indexed
if [ ! -f "${REFERENCE_DIR}/hg38.fasta.bwt" ]; then
    bwa index ${REFERENCE_DIR}/hg38.fasta
fi

# Align to host genome
bwa mem \
    -t ${NUM_THREADS} \
    ${REFERENCE_DIR}/hg38.fasta \
    ${OUTPUT_DIR}/${SAMPLE_ID}/host_removed/${SAMPLE_ID}_unclassified_1.fastq \
    2>&1 | samtools view -bS - \
    | samtools sort -o ${OUTPUT_DIR}/${SAMPLE_ID}/host_removed/${SAMPLE_ID}_host_aligned.bam

# Filter out host reads with coverage ≥80% and similarity ≥90%
samtools view -F 4 ${OUTPUT_DIR}/${SAMPLE_ID}/host_removed/${SAMPLE_ID}_host_aligned.bam \
    | awk '$5 >= 0.9*256 {print}' \
    | awk 'BEGIN{OFS="\t"} {split($6, cigar, /[A-Z]/); 
         total_len=0; for(i in cigar) total_len+=cigar[i]; 
         if(total_len/length($10) >= 0.8) print $1}' \
    > ${OUTPUT_DIR}/${SAMPLE_ID}/host_removed/${SAMPLE_ID}_host_reads.txt

# Extract non-host reads
seqtk subseq \
    ${OUTPUT_DIR}/${SAMPLE_ID}/host_removed/${SAMPLE_ID}_unclassified_1.fastq \
    ${OUTPUT_DIR}/${SAMPLE_ID}/host_removed/${SAMPLE_ID}_host_reads.txt \
    -r \
    > ${OUTPUT_DIR}/${SAMPLE_ID}/host_removed/${SAMPLE_ID}_nonhost.fastq

gzip ${OUTPUT_DIR}/${SAMPLE_ID}/host_removed/${SAMPLE_ID}_nonhost.fastq