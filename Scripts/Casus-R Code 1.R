setwd("C:/Users/Mahmut Otman/OneDrive/Documenten/R-Studio/Transcriptomics-Reuma/")

getwd()


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsubread")
BiocManager::install("DESeq2")
BiocManager::install("Rsamtools")


library(Rsubread)
library(DESeq2)
library(Rsamtools)

browseVignettes('Rsubread')
browseVignettes('DESeq2')

indexSplit = TRUE

buildindex(
  basename = 'human_genome',
  reference = 'Homo_sapiens.GRCh38.dna.toplevel.fa/Homo_sapiens.GRCh38.dna.toplevel.fa',
  memory = 12000,
  indexSplit = TRUE)

list.files()

# Normal samples
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

# RA samples
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


list.files(pattern = "human_genome")


samples <- c("normal1", "normal2", "normal3", "normal4", # Normaal
             "RA1", "RA2", "RA3", "RA4") # RA
lapply(samples, function(s) {
  sortBam(file = paste0(s, ".BAM"), destination = paste0(s, ".sorted"))
})

# Indexeren van gesorteerde BAM-bestanden
lapply(samples, function(s) {
  indexBam(paste0(s, ".sorted.bam"))
})


indexFa("Homo_sapiens.GRCh38.dna.toplevel.fa/Homo_sapiens.GRCh38.dna.toplevel.fa")


