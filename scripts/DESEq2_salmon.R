library("DESeq2")

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
# Check if input file path is provided
if (length(args) < 1) {
  stop("Metadata file path is missing.")
}
gorpol_origin <- read.table(input_file, header = FALSE, row.names = 1)
colnames(gorpol_origin) <- c("condition")
rownames(gorpol_origin) <- paste0(rownames(gorpol_origin), "_quant")

samples <- list.files(path = "salmon/quant",full.names = TRUE)
sf_files <- list.files(samples, pattern = "\\.sf$", full.names = TRUE)

data_frames <- list()
sample_origins <- list()
for (i in 1:length(sf_files)) {
  data_frames[[i]] <- read.table(sf_files[i], header = TRUE)
}
gorpol_reads <- do.call(cbind, lapply(data_frames, function(x) x$NumReads))

row.names(gorpol_reads) <- data_frames[[1]]$Name
colnames_gorpol <- list()
for (path in sf_files) {
  sample <- strsplit(path[[1]],"/")[[1]][3]
  colnames_gorpol <- c(colnames_gorpol,as.character(sample)) 
}

#colnames_gorpol <- lapply(colnames_gorpol, function(x) substr(x, 1, nchar(x) - 3))
colnames(gorpol_reads) <- colnames_gorpol
head(round(gorpol_reads))
print(gorpol_origin)
dds <- DESeqDataSetFromMatrix(countData = round(gorpol_reads), colData = gorpol_origin, design =~condition)

dds_gorpol <- DESeq(dds)

results_gorpol <- results(dds_gorpol)
summary_gorpol <- summary(results_gorpol)

results_gorpol_pval <- results_gorpol[order(results_gorpol$pvalue),]
head(results_gorpol_pval)

significant_genes <- results_gorpol_pval[!is.na(results_gorpol_pval$padj) & results_gorpol_pval$padj < 0.05, ]

sink("DESEq2_salmon/deseq2_summary.txt")
paste("The number of DGE found with pvalue of 0.05: ",nrow(significant_genes))
write.csv(results_gorpol_pval, "DESEq2_salmon/sign_genes_pvalue_0.05.csv")
summary(results_gorpol)
print("")
print(mcols(results_gorpol)$description)
sink()



