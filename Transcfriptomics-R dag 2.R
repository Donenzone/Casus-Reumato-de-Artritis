# Werkdirectory instellen
setwd("C:/Users/Mahmut Otman/OneDrive/Documenten/R-Studio/Transcriptomics-Reuma")
print(getwd())  # Controleren of working directory goed staat

# Packages installeren indien nodig
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("Rsubread", "dplyr", "Rsamtools", "readr"), ask = FALSE)

# Packages laden
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

# Count matrix en annotatie ophalen
counts <- count_matrix$counts
annotation <- count_matrix$annotation

# Kolomnamen aanpassen voor leesbaarheid
colnames(counts) <- c("normal1", "normal2", "normal3", "normal4", "RA1", "RA2", "RA3", "RA4")

# Eerste rijen en structuur van resultaten bekijken
head(annotation)
head(counts)
str(count_matrix)

# Opslaan van count matrix als CSV bestand
write.csv(counts, "bewerkt_countmatrix.csv")

# Samenvatting van total counts per gen (over alle samples)
summary(rowSums(counts))

# Overzicht tellingen van featureCounts
print(count_matrix$stat)

# Check aantal reads in één BAM bestand
print(countBam("normal1.BAM")$records)

# Sorteer count matrix op totaal aantal counts per gen (hoog naar laag) voor screenshot
counts_sorted <- counts[order(rowSums(counts), decreasing = TRUE), ]
head(counts_sorted)

