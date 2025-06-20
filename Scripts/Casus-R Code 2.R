# Zet de werkdirectory waar al je bestanden staan
setwd("C:/Users/Mahmut Otman/OneDrive/Documenten/R-Studio/Transcriptomics-Reuma")
print(getwd())  # Print om te controleren of het pad goed staat

# Installeer BiocManager als die nog niet geïnstalleerd is
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Installeer benodigde Bioconductor- en CRAN-packages zonder bevestiging (ask = FALSE)
BiocManager::install(c("Rsubread", "dplyr", "Rsamtools", "readr"), ask = FALSE)

# Laad de benodigde libraries
library(readr)       # Voor data in- en uitlezen
library(dplyr)       # Voor data-manipulatie
library(Rsamtools)   # Voor BAM-bestand verwerking
library(Rsubread)    # Voor read alignment en quantificatie

# ------------------------
# Reads tellen per gen
# ------------------------

# Geef de namen van alle BAM-bestanden op (gemaakt in vorige stap)
all_samples <- c(
  "normal1.BAM", "normal2.BAM", "normal3.BAM", "normal4.BAM",  # Controles
  "RA1.BAM", "RA2.BAM", "RA3.BAM", "RA4.BAM"                    # Reuma-patiënten
)

# Voer featureCounts uit om reads per gen te tellen
count_matrix <- featureCounts(
  files = all_samples,                              # BAM-bestanden
  annot.ext = "Homo_sapiens.GRCh38.114.gtf.gz",     # GTF-genannotatiebestand
  isPairedEnd = TRUE,                               # Aangeven dat reads paired-end zijn
  isGTFAnnotationFile = TRUE,                       # Bestandsformaat is GTF
  GTF.attrType = "gene_id",                         # Gebruik gene_id als unieke identifier
  useMetaFeatures = TRUE                            # Tel op gen-niveau (alle exons worden samengevoegd)
)

# Haal count data (tellingen) en annotaties (gen-informatie) apart uit het resultaat
counts <- count_matrix$counts         # Matrix: rijen = genen, kolommen = samples
annotation <- count_matrix$annotation # Gen-coördinaten zoals chromosoom, start, eind, etc.

# Geef herkenbare kolomnamen aan de matrix (voor overzicht en analyse)
colnames(counts) <- c("normal1", "normal2", "normal3", "normal4", "RA1", "RA2", "RA3", "RA4")

# Bekijk de eerste paar rijen van de annotatie en de countmatrix
head(annotation)
head(counts)

# Bekijk de structuur van het hele count_matrix-object
str(count_matrix)

# ------------------------
# Opslaan & Inspectie
# ------------------------

# Sla de countmatrix op als CSV-bestand voor latere analyse (bijv. in DESeq2)
write.csv(counts, "bewerkt_countmatrix.csv")

# Geef een samenvatting van het totaal aantal reads per gen (rowSums)
summary(rowSums(counts))  # Kan gebruikt worden om filtering te doen op lage expressie

# Bekijk statistieken van de featureCounts-aanroep, zoals:
# - hoeveel reads zijn gemapt
# - hoeveel zijn gebruikt
# - hoeveel zijn genegeerd
print(count_matrix$stat)

# Tel het aantal reads in één specifiek BAM-bestand
print(countBam("normal1.BAM")$records)

# Sorteer de genen op totaal aantal reads (hoog naar laag)
counts_sorted <- counts[order(rowSums(counts), decreasing = TRUE), ]
head(counts_sorted)  # Bekijk de top-expressie genen
