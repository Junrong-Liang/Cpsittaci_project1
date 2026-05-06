## Description
This project provides scripts and associated files to reproduce key steps of a bioinformatics pipeline for analyzing the epidemiology and phylodynamics of *Chlamydia psittaci*.
The pipeline covers: quality control of raw sequencing reads, removal of host genetic information, homology analysis, genome annotation, phylogenetic inference, and genotyping.

## Dependancies
* **beast v1.10.4** – Bayesian evolutionary analysis sampling trees; used for time‑scaled phylodynamic inference.
* **bwa v0.7.17‑r1188** Burrows‑Wheeler Aligner for mapping reads to a reference genome.
* **checkm v1.1.3** – Assesses the quality of genome bins recovered from metagenomes (completeness, contamination).
* **fastp 0.23.1** – FASTQ preprocessor for quality filtering, adapter trimming, and read correction.
* **gubbins v3.0.0** – Identifies and removes recombinant regions from core genome alignments.
* **iqtree 2.2.0** – Compute maximum‑likelihood phylogenetic tree inference with model selection.
* **kraken2 v2.1.2** – Taxonomic classification of sequencing reads using exact k‑mer matching.
* **panaroo v1.2.8** – Pangenome analysis pipeline that identifies core and accessory genes from genome annotations.
* **samtools v1.7** – Manipulates SAM/BAM alignment files (sorting, indexing, variant calling utilities).
* **spades v3.15.2** – De novo genome assembler for metagenomic and isolate data.
* **prokka v1.14.6** – Rapid prokaryotic genome annotation (coding sequences, rRNAs, tRNAs).
* **quast (latest)** – Quality assessment tool for genome assemblies (contiguity, misassemblies, gene content).
* **snp‑dist v0.7.0** – Computes pairwise SNP distances from a core genome alignment.
* **GrapeTree v1.5.0** – Visualizes and builds minimum spanning trees (MSTreeV2) from SNP distance matrices.
