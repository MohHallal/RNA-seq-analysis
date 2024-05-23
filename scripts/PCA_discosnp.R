# PCA on local computer

library(ape)
library(ade4)
library(adegenet)
library(pegas)
library(RColorBrewer)
library(ggsci)
library(ggplot2)
library(dplyr)

# import individuals list and dataset

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
metadata <- read.csv(input_file, sep = "\t")
# Check if input file path is provided
if (length(args) < 1) {
  stop("Metadata file path is missing.")
}

popsGorpol <- read.csv("metadata/metadata_deseq2.tsv", header = FALSE, sep = "")
Gorpol_vcf <- read.vcf("discosnp/snp_filter_3.recode.vcf", from = 1, to = 60000)
Gorpol_ind <- loci2genind(Gorpol_vcf)
# define pops in the genind file and check order
pop(Gorpol_ind) <- popsGorpol[,2]
Gorpol_ind$pop

# perform PCA

# no need to scale if no missing data
pcaGorpol <- dudi.pca(Gorpol_ind, scale=FALSE, nf = 10, scannf = FALSE)

pc_scores <- pcaGorpol$li
pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")
# print the result
data <-pc_scores[,1:5]
data$sample <- popsGorpol[,1]
data <- merge(data, metadata, by="sample", all = "TRUE")

gg1 <- ggplot(data, aes(x = Axis1, y = Axis2, label = name)) +
  geom_point(size = 10,aes(shape = origine, color = individu)) +
  geom_text(hjust = 0, vjust = 0, size = 6) +  # ajuster l'alignement du texte près des points
  labs(x = "PC1", y = "PC2") +  # titres des axes
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 18),  # Adjust legend text size
        legend.title = element_text(size = 18),
        plot.background = element_rect(fill = "white"))

gg2 <- ggplot(data, aes(x = Axis3, y = Axis4, label = name)) +
  geom_point(size = 10,aes(shape = origine, color = individu)) +
  geom_text(hjust = 0, vjust = 0, size = 6) +  # ajuster l'alignement du texte près des points
  labs(x = "PC3", y = "PC4") +  # titres des axes
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 18),  # Adjust legend text size
        legend.title = element_text(size = 18),
        plot.background = element_rect(fill = "white"))
ggsave("PCA_discosnp/PC1_PC2_discosnp.png", gg1, width = 10, height = 8, units = "in")
ggsave("PCA_discosnp/PC3_PC4_discosnp.png", gg2, width = 10, height = 8, units = "in")
sink("PCA_discosnp/PCA_discosnp_summary.txt")
summary(pcaGorpol)
data
sink()
