# Zet de werkdirectory naar de projectmap
setwd("C:/Users/Mahmut Otman/OneDrive/Documenten/R-Studio/Transcriptomics-Reuma/")
getwd()  # Controleer of de werkdirectory correct is ingesteld

# Installeer benodigde Bioconductor-pakketten als deze nog niet aanwezig zijn
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Installeer transcriptomics-gerelateerde packages
BiocManager::install("Rsubread")    # Voor het alignen en tellen van reads
BiocManager::install("DESeq2")      # Voor differentiële expressieanalyse
BiocManager::install("Rsamtools")   # Voor manipulatie van BAM-bestanden

# Laad de geïnstalleerde packages
library(Rsubread)
library(DESeq2)
library(Rsamtools)

# Bekijk de handleidingen van de gebruikte packages
browseVignettes('Rsubread')
browseVignettes('DESeq2')

# Instellen van de index voor de referentiegenoom
indexSplit = TRUE
buildindex(
  basename = 'human_genome',  # Basisnaam van de gegenereerde indexbestanden
  reference = 'Homo_sapiens.GRCh38.dna.toplevel.fa/Homo_sapiens.GRCh38.dna.toplevel.fa',  # Pad naar referentie FASTA
  memory = 12000,             # Geheugenlimiet in MB
  indexSplit = TRUE           # Splits de index voor snellere toegang
)

list.files()  # Controleer welke bestanden beschikbaar zijn

# === Aligneren van reads ===
# Normale (controle) monsters
align.normal1 <- align(index = "human_genome",
                       readfile1 = "Data_RA_raw/Data_RA_raw/SRR4785819_1_subset40k.fastq",
                       readfile2 = "Data_RA_raw/Data_RA_raw/SRR4785819_2_subset40k.fastq",
                       output_file = "normal1.BAM")

align.normal2 <- align(index = "human_genome",
                       readfile1 = "Data_RA_raw/Data_RA_raw/SRR4785820_1_subset40k.fastq",
                       readfile2 = "Data_RA_raw/Data_RA_raw/SRR4785820_2_subset40k.fastq",
                       output_file = "normal2.BAM")

align.normal3 <- align(index = "human_genome",
                       readfile1 = "Data_RA_raw/Data_RA_raw/SRR4785828_1_subset40k.fastq",
                       readfile2 = "Data_RA_raw/Data_RA_raw/SRR4785828_2_subset40k.fastq",
                       output_file = "normal3.BAM")

align.normal4 <- align(index = "human_genome",
                       readfile1 = "Data_RA_raw/Data_RA_raw/SRR4785831_1_subset40k.fastq",
                       readfile2 = "Data_RA_raw/Data_RA_raw/SRR4785831_2_subset40k.fastq",
                       output_file = "normal4.BAM")

# Reumatoïde artritis (RA) monsters
RA1 <- align(index = "human_genome",
             readfile1 = "Data_RA_raw/Data_RA_raw/SRR4785979_1_subset40k.fastq",
             readfile2 = "Data_RA_raw/Data_RA_raw/SRR4785979_2_subset40k.fastq",
             output_file = "RA1.BAM")

RA2 <- align(index = "human_genome",
             readfile1 = "Data_RA_raw/Data_RA_raw/SRR4785980_1_subset40k.fastq",
             readfile2 = "Data_RA_raw/Data_RA_raw/SRR4785980_2_subset40k.fastq",
             output_file = "RA2.BAM")

RA3 <- align(index = "human_genome",
             readfile1 = "Data_RA_raw/Data_RA_raw/SRR4785986_1_subset40k.fastq",
             readfile2 = "Data_RA_raw/Data_RA_raw/SRR4785986_2_subset40k.fastq",
             output_file = "RA3.BAM")

RA4 <- align(index = "human_genome",
             readfile1 = "Data_RA_raw/Data_RA_raw/SRR4785988_1_subset40k.fastq",
             readfile2 = "Data_RA_raw/Data_RA_raw/SRR4785988_2_subset40k.fastq",
             output_file = "RA4.BAM")

# Controleer of de indexbestanden correct zijn aangemaakt
list.files(pattern = "human_genome")

# === Sorteren van BAM-bestanden ===
samples <- c("normal1", "normal2", "normal3", "normal4",  # Controlemonsters
             "RA1", "RA2", "RA3", "RA4")                   # RA-monsters

# Sorteren van alle BAM-bestanden
lapply(samples, function(s) {
  sortBam(file = paste0(s, ".BAM"), destination = paste0(s, ".sorted"))
})

# Indexeren van de gesorteerde BAM-bestanden
lapply(samples, function(s) {
  indexBam(paste0(s, ".sorted.bam"))
})

# Indexeren van de referentie FASTA (vereist voor sommige downstream analyses)
indexFa("Homo_sapiens.GRCh38.dna.toplevel.fa/Homo_sapiens.GRCh38.dna.toplevel.fa")

# === Metadata aanmaken ===
# Maak een data frame met informatie over de samples
sample_info <- data.frame(
  Sample_ID = c("SRR4785819", "SRR4785820", "SRR4785828", "SRR4785831",
                "SRR4785979", "SRR4785980", "SRR4785986", "SRR4785988"),
  Leeftijd = c(31, 15, 31, 42, 54, 66, 60, 59),
  Geslacht = rep("female", 8),
  Groep = c(rep("Control", 4), rep("Rheumatoid arthritis", 4)),
  Diagnose_status = c(rep("Geen RA", 4), rep("RA >12 maanden", 4)),
  ACPA_status = c(rep("Negatief", 4), rep("Positief", 4)),
  Biopttype = rep("Synoviumbiopt", 8)
)

# Schrijf metadata weg naar een CSV-bestand in juiste mapstructuur
dir.create("data/metadata/samples", recursive = TRUE, showWarnings = FALSE)
write.csv(sample_info, file = "data/metadata/samples/sample_info.csv", row.names = FALSE)

# Laad en bekijk de opgeslagen metadata
sample_info_loaded <- read.csv("data/metadata/samples/sample_info.csv")
head(sample_info_loaded)
