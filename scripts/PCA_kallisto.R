if(!require(tidyverse)){
    install.packages("tidyverse")
}
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
metadata <- read.csv(input_file, sep = "\t")
# Check if input file path is provided
if (length(args) < 1) {
  stop("Metadata file path is missing.")
}
samples <- list.files(path = "kallisto/quant/",full.names = TRUE)
sf_files <- list.files(samples, pattern = "\\.tsv$", full.names = TRUE)

data_frames <- list()
sample_origins <- list()
for (i in 1:length(sf_files)) {
  data_frames[[i]] <- read.table(sf_files[i], header = TRUE)
}
gorpol_reads <- do.call(cbind, lapply(data_frames, function(x) x$est_counts))
# La correction des noms de colonnes et des lignes
rownames(gorpol_reads) <- data_frames[[1]]$target_id

colnames_gorpol <- list()
for (path in sf_files) {
  sample <- strsplit(path[[1]],"/")[[1]][4]
  colnames_gorpol <- c(colnames_gorpol,as.character(sample)) 
}
colnames_gorpol <- lapply(colnames_gorpol, function(x) substr(x, 1, nchar(x) - 6))

colnames(gorpol_reads) <- colnames_gorpol
gorpol_reads <- t(gorpol_reads)
zero_cols <- colSums(gorpol_reads == 0) == nrow(gorpol_reads)
gorpol_reads <- gorpol_reads[, !zero_cols]

pca_result <- prcomp(gorpol_reads, scale. = TRUE)
pc_scores <- pca_result$x

pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

# print the result
data <-pc_scores[,1:5]

data <- merge(data, metadata, by="sample", all = "TRUE")

gg1 <- ggplot(data, aes(x = PC1, y = PC2, label = name)) +
  geom_point(size = 10,aes(shape = origine, color = individu)) +
  geom_text(hjust = 0, vjust = 0, size = 6) +  # ajuster l'alignement du texte près des points
  labs(x = "PC1", y = "PC2") +  # titres des axes
  theme_minimal() +
  theme(legend.position = "top",
         legend.text = element_text(size = 18),  # Adjust legend text size
        legend.title = element_text(size = 18),
        plot.background = element_rect(fill = "white"))

gg2 <- ggplot(data, aes(x = PC3, y = PC4, label = name)) +
  geom_point(size = 10,aes(shape = origine, color = individu)) +
  geom_text(hjust = 0, vjust = 0, size = 6) +  # ajuster l'alignement du texte près des points
  labs(x = "PC3", y = "PC4") +  # titres des axes
  theme_minimal() +
  theme(legend.position = "top",
         legend.text = element_text(size = 18),  # Adjust legend text size
        legend.title = element_text(size = 18),
        plot.background = element_rect(fill = "white"))

ggsave("PCA_expression_kallisto/PC1_PC2_kallisto_expression.png", gg1, width = 10, height = 8, units = "in")
ggsave("PCA_expression_kallisto/PC3_PC4_kallisto_expression.png", gg2, width = 10, height = 8, units = "in")

sink("PCA_expression_kallisto/pca_summary.txt")
print(data)
print(summary(pca_result))
sink()
