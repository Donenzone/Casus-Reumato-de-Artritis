setwd("C:/Users/Mahmut Otman/OneDrive/Documenten/R-Studio/Transcriptomics-Reuma")
print(getwd())


library(readr)
library(dplyr)


counts <- read.delim("count_matrix.txt", row.names = 1, check.names = FALSE)
counts <- round(counts)


colnames(counts) <- c("normal1", "normal2", "normal3", "normal4", "RA1", "RA2", "RA3", "RA4")


head(counts)
str(counts)

summary(rowSums(counts))


counts_sorted <- counts[order(rowSums(counts), decreasing = TRUE), ]
head(counts_sorted)


write.csv(counts, "bewerkt_countmatrix.csv")

