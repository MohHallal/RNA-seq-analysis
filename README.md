# RNA-Seq Analysis Pipeline for Single-End Reads

This README file provides detailed instructions for setting up and executing the RNA-Seq analysis pipeline, which is specifically designed for single-end reads. The pipeline integrates several bioinformatics tools including Kallisto, Salmon, RSEM, DESeq2, and discoSNP, and generates outputs such as quantification, PCA plots, and significant gene lists.
# Important Notice

This pipeline has been developed and tested on the IFB cluster. We cannot guarantee that it will work exactly the same way on other platforms due to differences in system configurations, software versions, and environment settings.

# Pipeline Overview

The pipeline involves the following major steps:

    Indexing the reference transcriptome.
    Quantifying gene expression.
    Performing PCA analysis on quantification results.
    Performing differential expression analysis using DESeq2.
    Identifying SNPs using discoSNP and filtering VCF files.

# Prerequisites

Ensure the following modules are available and load them before running the pipeline:

    Snakemake
    Conda
    Kallisto
    Salmon
    RSEM
    Bowtie2
    Samtools
    discoSNP
    R

# Directory Structure

    data/: Directory containing the input FASTQ files (.fastq.gz).
    metadata/: Directory containing metadata files for DESeq2 and PCA.
    scripts/: Directory containing R scripts for PCA and DESeq2 analyses.
    config.yaml: Configuration file with parameters and resources.
    env.yaml: Conda environment file specifying the required packages.
    Snakefile: Snakemake file defining the workflow.

# Configuration

The config.yaml file includes parameters and resource allocations for different tools. Ensure that the sample names in the metadata files match the names of the read files without the .fastq.gz extension.
# Execution Instructions

    1- Prepare Data and Metadata:
        Place your FASTQ files in the data/ directory.
        Create metadata files (metadata_deseq2.tsv and metadata_PCA.tsv) and place them in the metadata/ directory.

    2- Load Necessary Modules:
module load snakemake
module load conda
    3- Run the Pipeline:
Use the following command to execute the pipeline. Adjust memory and CPU settings based on your available resources.
srun --mem=32000 snakemake --use-conda --cores 32
# Metadata Files

    metadata_deseq2.tsv: This file must include each sample in a separate line along with its corresponding condition, with tab separation between the sample name and condition, without column names.
    metadata_PCA.tsv: This file must include columns (sample, individu, origine, and name) even if they do not correspond to your context.

# Notes

    The pipeline is designed for single-end reads of RNA-Seq data.
    Ensure that the sample names in the metadata tables match the names of the read files without the .fastq.gz extension.
    It is important to use the same order of samples in the metadata files and the samples list in the config.yaml file.
    The pipeline generates output directories and log files for each step, which are useful for tracking progress and troubleshooting.

By following the instructions in this README, you will be able to execute the RNA-Seq analysis pipeline and obtain quantification, PCA, differential expression results, and SNPs from your single-end RNA-Seq data.
