counts <- read.delim("count_matrix.txt", row.names = 1, check.names = FALSE)
counts <- round(counts)
colnames(counts) <- c("normal1", "normal2", "normal3", "normal4", "RA1", "RA2", "RA3", "RA4")

treatment <- factor(c(rep("normal", 4), rep("RA", 4)))
treatment_table <- data.frame(treatment = treatment)
rownames(treatment_table) <- colnames(counts)

# Instal packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("pathview")

library(DESeq2)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(pathview)

# Maak DESeqDataSet aan
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = treatment_table,
                              design = ~ treatment)

# Voer analyse uit
dds <- DESeq(dds)
resultaten <- results(dds)

# Controleer of er weer rownames zijn
head(rownames(resultaten))


# Resultaten opslaan in een bestand
#Bij het opslaan van je tabel kan je opnieuw je pad instellen met `setwd()` of het gehele pad waar je de tabel wilt opslaan opgeven in de code.
write.table(resultaten, file = 'ResultatenWC3.csv', row.names = TRUE, col.names = TRUE)

sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)

hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]

head(laagste_p_waarde)

EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj')

# Alternatieve plot zonder p-waarde cutoff (alle genen zichtbaar)
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0)

dev.copy(png, 'VolcanoplotWC.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()


# Zorg dat je genen (rownames(resultaten)) Entrez-ID's krijgen
gene_info <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = rownames(resultaten),
                                   keytype = "SYMBOL", # of "ENSEMBL" indien van toepassing
                                   columns = c("ENTREZID", "SYMBOL"))


# Voeg EntrezID toe aan je resultaten
resultaten$symbol <- rownames(resultaten)
resultaten <- merge(as.data.frame(resultaten), gene_info, by.x = "symbol", by.y = "SYMBOL")

# Verwijder NA’s en dupliceerde ENTREZID’s
resultaten <- resultaten[!is.na(resultaten$ENTREZID), ]
resultaten <- resultaten[!duplicated(resultaten$ENTREZID), ]

# Maak gene_vector: named vector met log2FC per EntrezID
gene_vector <- resultaten$log2FoldChange
names(gene_vector) <- resultaten$ENTREZID


pathview(
  gene.data = gene_vector,
  pathway.id = "hsa05323",   # KEGG ID voor Rheumatoid arthritis
  species = "hsa",           # Humane species
  gene.idtype = "ENTREZ",    # Omdat gene_vector EntrezID's bevat
  limit = list(gene = 5)     # Log2FC bereik tussen -5 en 5
)

str(resultaten)
head(resultaten)

head(rownames(resultaten))

________________________________________________________________________________
Go-Analyse

library(goseq)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)

# Voorbeeld dataset 'resultaten' moet Entrez IDs en padj bevatten
# significantie cutoff definiëren
significant_genes <- resultaten$padj < 0.05 & abs(resultaten$log2FoldChange) > 1

# Named vector maken (1 = DEG, 0 = niet DEG)
gene_vector <- integer(nrow(resultaten))
names(gene_vector) <- resultaten$ENTREZID
gene_vector[significant_genes] <- 1

# Verwijder genen zonder geldige EntrezID
gene_vector <- gene_vector[!is.na(names(gene_vector)) & names(gene_vector) != ""]

# Exons per gen ophalen
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
exons <- exonsBy(txdb, by = "gene")

# Per gen exonlengte berekenen
gene_lengths <- sapply(exons, function(x) sum(width(reduce(x))))

# Filter op genen in gene_vector
gene_lengths <- gene_lengths[names(gene_lengths) %in% names(gene_vector)]

# Zet vectors in dezelfde volgorde
gene_vector <- gene_vector[names(gene_lengths)]

# Bereken PWF
pwf <- nullp(gene_vector, bias.data = gene_lengths)

# GO annotaties ophalen
all_genes <- names(gene_vector)
go_annotations <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = all_genes,
                                        columns = "GO",
                                        keytype = "ENTREZID")

# Filter op Biological Process
go_annotations <- go_annotations[go_annotations$ONTOLOGY == "BP", ]

# Gen-naar-GO mapping
gene2cat <- split(go_annotations$GO, go_annotations$ENTREZID)

# GOseq uitvoeren
goseq_results <- goseq(pwf, gene2cat = gene2cat)

# Resultaten bekijken
head(goseq_results)

# Top 10 visualiseren
sig_go <- goseq_results %>%
  filter(over_represented_pvalue < 0.05) %>%
  arrange(over_represented_pvalue) %>%
  head(10)

ggplot(sig_go, aes(x = reorder(term, -over_represented_pvalue), y = -log10(over_represented_pvalue))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "GO term", y = "-log10(p-value)", title = "Top 10 overgerepresenteerde GO termen") +
  theme_minimal()




# Stel: je hebt al je goseq_results (uit goseq) klaarliggen
# Filter top 10 op laagste p-waarde (over_represented_pvalue)
top10 <- goseq_results %>%
  filter(!is.na(over_represented_pvalue)) %>%
  arrange(over_represented_pvalue) %>%
  head(10) %>%
  mutate(hitsPerc = numDEInCat * 100 / numInCat)

# Haal GO termen (beschrijvingen) erbij met GO.db
top10$term <- Term(GOTERM[top10$category])

# Maak dotplot
ggplot(top10, aes(x = hitsPerc, y = term,
                  colour = over_represented_pvalue,
                  size = numDEInCat)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", trans = "log10") +
  expand_limits(x = 0) +
  labs(x = "Hits (%)",
       y = "GO term",
       colour = "p value",
       size = "Count",
       title = "Top 10 Overrepresented GO terms") +
  theme_minimal()
