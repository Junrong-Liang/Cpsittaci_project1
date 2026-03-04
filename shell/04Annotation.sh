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
# STEP 5: VIRULENCE GENE ANNOTATION
# ==============================================
echo "Step 5: Virulence gene annotation for ${SAMPLE_ID}"

# Gene prediction with Pyrodigal
pyrodigal \
    -i ${OUTPUT_DIR}/${SAMPLE_ID}/assembly/contigs.fasta \
    -o ${OUTPUT_DIR}/${SAMPLE_ID}/virulence/${SAMPLE_ID}_genes.fna \
    -a ${OUTPUT_DIR}/${SAMPLE_ID}/virulence/${SAMPLE_ID}_proteins.faa \
    -d ${OUTPUT_DIR}/${SAMPLE_ID}/virulence/${SAMPLE_ID}_genes.ffn \
    --min-gene 90 \
    --max-gene 6000 \
    2>&1 | tee ${OUTPUT_DIR}/logs/${SAMPLE_ID}_step5_pyrodigal.log

# BLAST against VFDB
makeblastdb \
    -in ${REFERENCE_DIR}/VFDB.fasta \
    -dbtype nucl \
    -out ${REFERENCE_DIR}/VFDB

blastn \
    -query ${OUTPUT_DIR}/${SAMPLE_ID}/virulence/${SAMPLE_ID}_genes.ffn \
    -db ${REFERENCE_DIR}/VFDB \
    -out ${OUTPUT_DIR}/${SAMPLE_ID}/virulence/${SAMPLE_ID}_vfdb_blast.txt \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
    -evalue 1e-5 \
    -num_threads ${NUM_THREADS} \
    2>&1 | tee ${OUTPUT_DIR}/logs/${SAMPLE_ID}_step5_blastn.log

# ==============================================
# STEP 6: ANTIMICROBIAL RESISTANCE GENE IDENTIFICATION
# ==============================================
echo "Step 6: AMR gene identification for ${SAMPLE_ID}"

rgi main \
    -i ${OUTPUT_DIR}/${SAMPLE_ID}/assembly/contigs.fasta \
    -o ${OUTPUT_DIR}/${SAMPLE_ID}/amr/${SAMPLE_ID}_rgi \
    -t contig \
    -a BLAST \
    -n ${NUM_THREADS} \
    --clean \
    2>&1 | tee ${OUTPUT_DIR}/logs/${SAMPLE_ID}_step6_rgi.log

# ==============================================
# STEP 7: VARIANT CALLING
# ==============================================
echo "Step 7: Variant calling for ${SAMPLE_ID}"

# For NCBI samples (use contigs as input)
if [[ ${SAMPLE_ID} == *"NCBI"* ]]; then
    snippy \
        --ref ${REFERENCE_DIR}/Chlamydia_psittaci_6BC.fasta \
        --contigs ${OUTPUT_DIR}/${SAMPLE_ID}/assembly/contigs.fasta \
        --outdir ${OUTPUT_DIR}/${SAMPLE_ID}/variants \
        --prefix ${SAMPLE_ID} \
        --cpus ${NUM_THREADS} \
        2>&1 | tee ${OUTPUT_DIR}/logs/${SAMPLE_ID}_step7_snippy.log
# For CDC samples (use reads as input)
elif [[ ${SAMPLE_ID} == *"CDC"* ]]; then
    snippy \
        --ref ${REFERENCE_DIR}/Chlamydia_psittaci_6BC.fasta \
        --pe1 ${OUTPUT_DIR}/${SAMPLE_ID}/alignment/${SAMPLE_ID}_cp_candidate.fastq.gz \
        --outdir ${OUTPUT_DIR}/${SAMPLE_ID}/variants \
        --prefix ${SAMPLE_ID} \
        --cpus ${NUM_THREADS} \
        2>&1 | tee ${OUTPUT_DIR}/logs/${SAMPLE_ID}_step7_snippy.log
fi

# ==============================================
# STEP 8: T3SS EFFECTOR PROTEIN PREDICTION
# ==============================================
echo "Step 8: T3SS effector prediction for ${SAMPLE_ID}"

java -jar /path/to/TTSS_STD_2.0.2_src_all.jar \
    -m TTSS_STD_2.0.2_src_all.jar \
    -t selective \
    -i ${OUTPUT_DIR}/${SAMPLE_ID}/virulence/${SAMPLE_ID}_proteins.faa \
    -o ${OUTPUT_DIR}/${SAMPLE_ID}/t3ss/${SAMPLE_ID}_t3ss_results.txt \
    2>&1 | tee ${OUTPUT_DIR}/logs/${SAMPLE_ID}_step8_effectivet3.log
