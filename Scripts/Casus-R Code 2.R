setwd("C:/Users/Mahmut Otman/OneDrive/Documenten/R-Studio/Transcriptomics-Reuma")
print(getwd())


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("Rsubread", "dplyr", "Rsamtools", "readr"), ask = FALSE)


library(readr)
library(dplyr)
library(Rsamtools)
library(Rsubread)

# BAM-bestanden in vector opslaan
all_samples <- c(
  "normal1.BAM", "normal2.BAM", "normal3.BAM", "normal4.BAM",
  "RA1.BAM", "RA2.BAM", "RA3.BAM", "RA4.BAM"
)

# featureCounts uitvoeren (paired-end reads)
count_matrix <- featureCounts(
  files = all_samples,
  annot.ext = "Homo_sapiens.GRCh38.114.gtf.gz",  # GTF annotatiebestand
  isPairedEnd = TRUE,                            # TRUE omdat BAM's paired-end zijn
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)


counts <- count_matrix$counts
annotation <- count_matrix$annotation


colnames(counts) <- c("normal1", "normal2", "normal3", "normal4", "RA1", "RA2", "RA3", "RA4")


head(annotation)
head(counts)
str(count_matrix)

# Opslaan van count matrix als CSV bestand
write.csv(counts, "bewerkt_countmatrix.csv")


summary(rowSums(counts))

print(count_matrix$stat)

print(countBam("normal1.BAM")$records)


counts_sorted <- counts[order(rowSums(counts), decreasing = TRUE), ]
head(counts_sorted)

