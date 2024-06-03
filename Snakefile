configfile: "config.yaml"
rule all:
    input:
        "PCA_discosnp/PC1_PC2_discosnp.png",
       	"PCA_discosnp/PC3_PC4_discosnp.png"
        "discosnp/snp_filter_3.recode.vcf"
        "discosnp/discoRes_k_31_c_3_D_100_P_3_b_0_coherent.vcf"
#        expand("RSEM/quant/{sample}_quant/{sample}_quant",sample=config["samples"])
rule index_kallisto:
    input:
        config["fasta_ref"]
    output:
        "kallisto/index/kallisto_index"
    conda:
        "env.yaml"
    resources:
        mem_gb=config["kallisto_index_mem_gb"]
    threads: config["kallisto_index_threads"]
    log:
        "kallisto/index/kallisto_index.log"
    shell:
        """
        mkdir -p kallisto
        kallisto index {input} -i {output}
        """
rule count_kallisto_single:
    input:
        index="kallisto/index/kallisto_index",
        sample="data/{sample}.fastq.gz"
    output:
        quant_file=directory("kallisto/quant/{sample}_quant")
    conda:
        "env.yaml"
    resources:
        mem_gb=config["kallisto_quant_mem_gb"]
    threads: config["kallisto_quant_threads"]
    log:
        "kallisto/quant/{sample}_quant/kallisto_quant.log"
    shell:
        """
        mkdir -p {output.quant_file}
        kallisto quant {input.sample} -i {input.index} -o {output.quant_file} --single -l {config[l]} -s {config[s]} -t {threads}
        """
rule PCA_kallisto:
    input:
        expand("kallisto/quant/{sample}_quant", sample=config["samples"]),
        "metadata/metadata_PCA.tsv",
        "scripts/PCA_kallisto.R"
    output:
        "PCA_expression_kallisto/PC1_PC2_kallisto_expression.png",
        "PCA_expression_kallisto/PC3_PC4_kallisto_expression.png"    
    resources:
        mem_gb=config["pca_kallisto_mem_gb"]
    threads: config["pca_kallisto_threads"]
    conda:
        "env.yaml"
    shell:
        """
       	module load r
        mkdir -p PCA_expression_kallisto
        Rscript scripts/PCA_kallisto.R metadata/metadata_PCA.tsv 
        """
rule DESEq2_kallisto:
    input:
        expand("kallisto/quant/{sample}_quant", sample=config["samples"]),
        "metadata/metadata_deseq2.tsv"
    output:
        "DESEq2_kallisto/sign_genes_pvalue_0.05.csv"
    resources:
        mem_gb=config["deseq2_kallisto_mem_gb"]
    threads: config["deseq2_kallisto_threads"]
    shell:
        """
        module load r
        mkdir -p DESEq2_kallisto
        Rscript scripts/DESEq2_kallisto.R metadata/metadata_deseq2.tsv
        """
rule index_salmon:
    input:
        config["fasta_ref"]
    output:
        directory("salmon/index/")
    conda:
        "env.yaml"
    resources:
        mem_gb=config["salmon_index_mem_gb"]
    threads: config["salmon_index_threads"]
    shell:
        """
        mkdir -p salmon
        salmon index -t {input} -i {output} -k {config[k_salmon]} -p {config[salmon_index_threads]}
        """
rule quant_salmon:
    input:
        index="salmon/index/",
        sample="data/{sample}.fastq.gz"
    output:
        quant_file=directory("salmon/quant/{sample}_quant")
    conda:
        "env.yaml"
    resources:
        mem_gb=config["salmon_quant_mem_gb"]
    threads: config["salmon_quant_threads"]
    log:
        "salmon/quant/{sample}_quant/salmon_quant.log"
    shell:
        """
        mkdir -p {output.quant_file}
        salmon quant -i {input.index} -l {config[libtype_salmon]} -r {input.sample} -o {output.quant_file} -p {config[salmon_quant_threads]}
        """
rule PCA_salmon:
    input:
        expand("salmon/quant/{sample}_quant", sample=config["samples"]),
        "metadata/metadata_PCA.tsv"
    output:
        "PCA_expression_salmon/PC1_PC2_salmon_expression.png",
        "PCA_expression_salmon/PC3_PC4_salmon_expression.png"
    resources:
        mem_gb=config["pca_salmon_mem_gb"]
    threads: config["pca_salmon_threads"]
    conda:
        "env.yaml"
    shell:
        """
        module load r
        mkdir -p PCA_expression_salmon
        Rscript scripts/PCA_salmon.R metadata/metadata_PCA.tsv
        """
rule DESEq2_salmon:
    input:
        expand("salmon/quant/{sample}_quant", sample=config["samples"]),
        "metadata/metadata_deseq2.tsv"
    output:
        "DESEq2_salmon/sign_genes_pvalue_0.05.csv"
    resources:
        mem_gb=config["deseq2_salmon_mem_gb"]
    threads: config["deseq2_salmon_threads"]
    shell:
        """
        module load r
        mkdir -p DESEq2_salmon
        Rscript scripts/DESEq2_salmon.R metadata/metadata_deseq2.tsv
        """
rule RSEM_index:
    input:
        config["fasta_ref"]
    output:
        "RSEM/index"
    conda:
        "env.yaml"
    resources:
        mem_gb=config["RSEM_index_mem_gb"]
    threads: config["RSEM_index_threads"]
    log:
        "RSEM/index/rsem_index.log"
    shell:
        """
        module load rsem
        module load bowtie2
        module load samtools
        mkdir -p RSEM/index
        rsem-prepare-reference {input} {output}/quant --bowtie2 -p {threads}
        """
rule RSEM_quant:
    input:
        index="RSEM/index",
        sample="data/{sample}.fastq.gz"
    output:
        quant_file=directory("RSEM/quant/{sample}_quant/{sample}_quant")
    conda:
        "env.yaml"
    resources:
        mem_gb=config["rsem_quant_mem_gb"]
    threads: config["rsem_quant_threads"]
    log:
        "RSEM/quant/{sample}_quant/rsem_quant.log"
    shell:
        """
        module load rsem
        module load bowtie2
        module load bowtie
        module load samtools
        mkdir lib
        cd lib
        ln -s /shared/ifbstor1/software/miniconda/envs/samtools-1.13/lib/libcrypto.so libcrypto.so.1.0.0
        export LD_LIBRARY_PATH=$(pwd) 
        samtools --version
        cd ..
        mkdir -p {output.quant_file}
        rsem-calculate-expression {input.sample} {input.index}/index {output} --bowtie2 -p {threads} 
        rm -r lib
        """
rule pop_list_discosnp:
    output:
        "discosnp/pop_list.txt"
    run:
        import subprocess
        import os
        shell("mkdir -p discosnp")
        snakefile_dir = os.path.abspath(os.getcwd())
        data_dir = os.path.join(snakefile_dir, "data")
        samples=config["samples"]
        with open(output[0], "w") as f:
            for sample in samples:
                full_path = os.path.join(data_dir, sample + ".fastq.gz")
                f.write(full_path + "\n")
rule discosnp:
    output:
        "discosnp/discoRes_k_31_c_3_D_100_P_3_b_0_coherent.vcf"
    input:
        "discosnp/pop_list.txt"
    resources:
         mem_gb=config["discosnp_mem_gb"]
    threads: config["discosnp_threads"]
    conda:
        "env.yaml"
    shell:
        """
        module load discosnp/2.4.3
        run_discoSnp++.sh -r {input} -g {config[fasta_ref]}
        mv discoRes* discosnp/
        """
rule filter_vcf:
    output:
        "discosnp/snp_filter_3.recode.vcf"
    input:
        "discosnp/discoRes_k_31_c_3_D_100_P_3_b_0_coherent.vcf"
    conda:
        "env.yaml"
    shell:
        """
        vcftools --vcf discosnp/discoRes_k_31_c_3_D_100_P_3_b_0_coherent.vcf --max-missing 0.90 --min-alleles 2 --max-alleles 2 --recode --out discosnp/snp_filter_1
        vcftools --vcf discosnp/snp_filter_1.recode.vcf --max-missing 1 --min-alleles 2 --max-alleles 2 --recode --out discosnp/snp_filter_2
        vcftools --vcf discosnp/snp_filter_2.recode.vcf --thin 1000 --recode --out discosnp/snp_filter_3
        vcftools --vcf discosnp/snp_filter_3.recode.vcf --missing-indv
        mv out.imiss discosnp/
        """
rule PCA_discosnp:
    output:
        "PCA_discosnp/PC1_PC2_discosnp.png",
        "PCA_discosnp/PC3_PC4_discosnp.png"
    input:
        "discosnp/snp_filter_3.recode.vcf",
        "metadata/metadata_PCA.tsv"
    shell:
        """
        module load r
        mkdir -p PCA_discosnp
        Rscript scripts/PCA_discosnp.R metadata/metadata_PCA.tsv
        """

