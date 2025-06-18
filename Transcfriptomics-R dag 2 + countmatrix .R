# Werkdirectory instellen
setwd("C:/Users/Mahmut Otman/OneDrive/Documenten/R-Studio/Transcriptomics-Reuma")
print(getwd())  # Controleren of je in de juiste map zit

# Packages laden
library(readr)
library(dplyr)

# Countmatrix inladen die al is gemaakt (vanuit WC2 of gegeven bestand)
counts <- read.delim("count_matrix.txt", row.names = 1, check.names = FALSE)

counts <- round(counts)

# Kolomnamen vervangen door sample-namen (zoals eerder gedaan in WC2)
colnames(counts) <- c("normal1", "normal2", "normal3", "normal4", "RA1", "RA2", "RA3", "RA4")

# Structuur inspecteren
head(counts)
str(counts)

# Samenvatting total counts per gen (voor filtering of controle)
summary(rowSums(counts))

# Genen sorteren op totaal aantal counts (hoog naar laag)
counts_sorted <- counts[order(rowSums(counts), decreasing = TRUE), ]
head(counts_sorted)

# Opslaan van aangepaste countmatrix (optioneel)
write.csv(counts, "bewerkt_countmatrix.csv")

