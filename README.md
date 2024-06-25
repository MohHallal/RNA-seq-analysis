# RNA-seq and Variant Calling Pipeline

This repository contains a comprehensive bioinformatics pipeline for RNA-seq data analysis and variant calling. The pipeline is implemented using Snakemake and includes steps for indexing, quantification, principal component analysis (PCA), differential expression analysis, and variant calling.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Contributing](#contributing)
- [License](#license)

## Overview

The pipeline includes the following key steps:

1. **Indexing and Quantification:**
   - Kallisto
   - Salmon
   - RSEM

2. **Principal Component Analysis (PCA):**
   - Kallisto
   - Salmon
   - DiscoSNP++

3. **Differential Expression Analysis:**
   - Kallisto
   - Salmon

4. **Variant Calling:**
   - DiscoSNP++

## Installation

To run this pipeline, you need to have Conda installed. Follow the steps below to set up the environment:

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/rna-seq-variant-calling-pipeline.git
   cd rna-seq-variant-calling-pipeline
   
Create and activate the Conda environment: 

conda env create -f env.yaml
conda activate rna-seq-pipeline

# Configuration

The pipeline requires a configuration file (config.yaml) to specify various parameters. Default values are precised in the example configuration file.
    
# Usage

To run the pipeline, use the following command:

snakemake --cores <number_of_cores>

You can also specify additional options such as --use-conda to ensure that the Conda environment is used for each rule.
Pipeline Steps
# Indexing and Quantification

    Kallisto:
        index_kallisto: Creates an index for Kallisto using a reference FASTA file.
        count_kallisto_single: Quantifies RNA-seq data using Kallisto for single-end reads.

    Salmon:
        index_salmon: Creates an index for Salmon using a reference FASTA file.
        quant_salmon: Quantifies RNA-seq data using Salmon.

    RSEM:
        RSEM_index: Prepares a reference for RSEM using Bowtie2.
        RSEM_quant: Quantifies RNA-seq data using RSEM.

# Principal Component Analysis (PCA)

    Kallisto:
        PCA_kallisto: Performs PCA on Kallisto quantification results and generates plots.

    Salmon:
        PCA_salmon: Performs PCA on Salmon quantification results and generates plots.

    DiscoSNP++:
        PCA_discosnp: Performs PCA on variant calling results from DiscoSNP++ and generates plots.
# Differential Expression Analysis

    Kallisto:
        DESEq2_kallisto: Performs differential expression analysis using DESeq2 on Kallisto quantification results.

    Salmon:
        DESEq2_salmon: Performs differential expression analysis using DESeq2 on Salmon quantification results.

# Variant Calling

    DiscoSNP++:
        pop_list_discosnp: Generates a list of sample files for DiscoSNP++.
        discosnp: Runs DiscoSNP++ for variant calling.
        filter_vcf: Filters the VCF file generated by DiscoSNP++.

# Contributing

Contributions are welcome! Please open an issue or submit a pull request for any improvements or bug fixes.
# License

This project is licensed under the MIT License. See the LICENSE file for details.
