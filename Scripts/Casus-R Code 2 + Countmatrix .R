# Zet de werkdirectory waar je data zich bevindt
setwd("C:/Users/Mahmut Otman/OneDrive/Documenten/R-Studio/Transcriptomics-Reuma")
print(getwd())  # Check of het werkdirectory correct is ingesteld

# Laad benodigde libraries
library(readr)   # Voor het inlezen van tab-delimited bestanden (zoals .txt of .tsv)
library(dplyr)   # Voor data-manipulatie (optioneel in dit script)

# -----------------------------
# Inlezen van count matrix
# -----------------------------

# Lees de tellingsmatrix in die eerder gegenereerd is met featureCounts
# row.names = 1: gebruik de eerste kolom (gene-ID's) als rij-namen
# check.names = FALSE: voorkomt dat R automatisch kolomnamen aanpast (zoals "." in plaats van "-")
counts <- read.delim("count_matrix.txt", row.names = 1, check.names = FALSE)

# Afronden van waarden: featureCounts levert soms decimalen (bijv. bij gemiddeldes); deze worden afgerond
counts <- round(counts)

# Kolomnamen toewijzen voor duidelijkheid (per sample)
colnames(counts) <- c("normal1", "normal2", "normal3", "normal4", "RA1", "RA2", "RA3", "RA4")

# -----------------------------
# Verkenning van de data
# -----------------------------

# Bekijk de eerste paar rijen van de matrix (gen-ID’s met counts)
head(counts)

# Bekijk de structuur van het object (aantal genen, type matrix, etc.)
str(counts)

# Samenvatting van het totaal aantal reads per gen (rowSums)
# Handig voor het identificeren van genen met lage expressie
summary(rowSums(counts))

# Sorteer de genen op totaal aantal reads (hoogste expressie bovenaan)
counts_sorted <- counts[order(rowSums(counts), decreasing = TRUE), ]
head(counts_sorted)  # Bekijk de top-genen met hoogste totale read count

# -----------------------------
# Opslaan voor latere analyse
# -----------------------------

# Sla de (onbewerkte, afgeronde) count matrix op als CSV-bestand
write.csv(counts, "bewerkt_countmatrix.csv")
