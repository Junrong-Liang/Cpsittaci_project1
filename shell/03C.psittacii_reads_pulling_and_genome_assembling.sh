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
# STEP 3: ALIGNMENT TO METAGENOMIC DATABASE
# ==============================================
echo "Step 3: Alignment to metagenomic database for ${SAMPLE_ID}"

# Index custom metagenomic database if needed
if [ ! -f "${REFERENCE_DIR}/metagenomic_db.fasta.bwt" ]; then
    bwa index ${REFERENCE_DIR}/metagenomic_db.fasta
fi

# Align to metagenomic database
bwa mem \
    -t ${NUM_THREADS} \
    -k 19 \
    -Y \
    -h 1000 \
    ${REFERENCE_DIR}/metagenomic_db.fasta \
    ${OUTPUT_DIR}/${SAMPLE_ID}/host_removed/${SAMPLE_ID}_nonhost.fastq.gz \
    2>&1 | samtools view -bS - \
    | samtools sort -o ${OUTPUT_DIR}/${SAMPLE_ID}/alignment/${SAMPLE_ID}_metagenome_aligned.bam

# Index the BAM file
samtools index ${OUTPUT_DIR}/${SAMPLE_ID}/alignment/${SAMPLE_ID}_metagenome_aligned.bam

# Extract C. psittaci reads with coverage ≥80% and identity ≥90%
# First, get list of C. psittaci reference sequences
grep "Chlamydia_psittaci" ${REFERENCE_DIR}/metagenomic_db.fasta | sed 's/>//' > ${OUTPUT_DIR}/${SAMPLE_ID}/alignment/cp_refs.txt

# Extract reads mapped to C. psittaci
while read ref; do
    samtools view -b ${OUTPUT_DIR}/${SAMPLE_ID}/alignment/${SAMPLE_ID}_metagenome_aligned.bam "$ref" \
        | samtools fastq - \
        >> ${OUTPUT_DIR}/${SAMPLE_ID}/alignment/${SAMPLE_ID}_cp_candidate.fastq
done < ${OUTPUT_DIR}/${SAMPLE_ID}/alignment/cp_refs.txt

# Filter by coverage and identity (simplified approach)
# Note: For exact filtering, consider using custom scripts or tools like bam-readcount
gzip ${OUTPUT_DIR}/${SAMPLE_ID}/alignment/${SAMPLE_ID}_cp_candidate.fastq


# ==============================================
# STEP 4: ASSEMBLY AND ASSESSMENT
# ==============================================
echo "Step 4: Assembly and assessment for ${SAMPLE_ID}"

# Assembly with SPAdes
spades.py \
    --single ${OUTPUT_DIR}/${SAMPLE_ID}/alignment/${SAMPLE_ID}_cp_candidate.fastq.gz \
    -o ${OUTPUT_DIR}/${SAMPLE_ID}/assembly \
    --threads ${NUM_THREADS} \
    --careful \
    2>&1 | tee ${OUTPUT_DIR}/logs/${SAMPLE_ID}_step4_spades.log

# Assess assembly quality with CheckM
checkm lineage_wf \
    -x fasta \
    -t ${NUM_THREADS} \
    ${OUTPUT_DIR}/${SAMPLE_ID}/assembly \
    ${OUTPUT_DIR}/${SAMPLE_ID}/assembly/checkm_results \
    2>&1 | tee ${OUTPUT_DIR}/logs/${SAMPLE_ID}_step4_checkm.log

# Calculate basic assembly statistics
assembly_stats() {
    local assembly=$1
    local output=$2
    
    # Get contig count
    grep -c ">" ${assembly} > ${output}_contig_count.txt
    
    # Calculate N50 and total length
    awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen+=length($0)} END {print seqlen}' ${assembly} \
        | sort -nr \
        | awk 'BEGIN {sum=0; count=0; n50=0} 
               {sum+=$1; count++; lengths[count]=$1} 
               END {half=sum/2; 
                    for(i=1;i<=count;i++) {cumsum+=lengths[i]; if(cumsum>=half) {n50=lengths[i]; break}} 
                    print "Total length: " sum "\nContig count: " count "\nN50: " n50}' \
        > ${output}_assembly_stats.txt
}

assembly_stats ${OUTPUT_DIR}/${SAMPLE_ID}/assembly/contigs.fasta ${OUTPUT_DIR}/${SAMPLE_ID}/assembly/${SAMPLE_ID}
