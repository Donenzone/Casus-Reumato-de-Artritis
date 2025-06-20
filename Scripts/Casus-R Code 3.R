# === 1. Data inladen en voorbereiden ===
# Lees de count matrix in en rond de waarden af naar hele getallen
counts <- read.delim("count_matrix.txt", row.names = 1, check.names = FALSE)
counts <- round(counts)

# Hernoem de kolomnamen voor duidelijkheid (4 normale monsters, 4 RA monsters)
colnames(counts) <- c("normal1", "normal2", "normal3", "normal4", "RA1", "RA2", "RA3", "RA4")

# Maak een factorvector voor de condities (groepsindeling)
treatment <- factor(c(rep("normal", 4), rep("RA", 4)))

# Zet de treatment vector om naar een data frame en koppel aan samples
treatment_table <- data.frame(treatment = treatment)
rownames(treatment_table) <- colnames(counts)

# === 2. Installeer en laad benodigde packages (indien nog niet aanwezig) ===
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Installeer Bioconductor packages
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("pathview")

# Laad de libraries in
library(DESeq2)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(pathview)

# === 3. Differentiële genexpressie analyse ===
# Maak een DESeqDataSet object aan met de count data en groepsindeling
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = treatment_table,
                              design = ~ treatment)

# Voer de DESeq-analyse uit (normalisatie + statistiek)
dds <- DESeq(dds)

# Haal resultaten op
resultaten <- results(dds)

# === 4. Visualisatie: MA-plot ===
plotMA(resultaten, ylim = c(-5, 5))  # MA-plot toont gemiddelde expressie vs. verandering

# Sla MA-plot op als PNG
png("data/processed/deseq2_results/MAplot.png", width = 6, height = 6, units = "in", res = 300)
plotMA(resultaten, ylim = c(-5, 5))
dev.off()

# === 5. Opslaan van resultaten ===
# Maak map aan als die nog niet bestaat
dir.create("data/processed/deseq2_results", recursive = TRUE, showWarnings = FALSE)

# Sla DESeqDataSet object op
saveRDS(dds, file = "data/processed/deseq2_results/dds.rds")

# Bekijk de eerste genen
head(rownames(resultaten))

# Schrijf resultaten weg naar CSV
write.table(resultaten, file = 'ResultatenWC3.csv', row.names = TRUE, col.names = TRUE)

# Tel significante genen met log2FC > 1 en padj < 0.05 (up/down gereguleerd)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)

# Sorteer resultaten op log2FC en p-waarde
hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]
head(laagste_p_waarde)

# === 6. Volcano plot ===
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj')

# Alternatieve volcano plot zonder drempelwaarden (alle genen zichtbaar)
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0)

# Sla Volcano plot op
dev.copy(png, 'VolcanoplotWC.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

# === 7. Annotatie van genen met EntrezID (nodig voor KEGG & GO analyses) ===
gene_info <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = rownames(resultaten),
                                   keytype = "SYMBOL",
                                   columns = c("ENTREZID", "SYMBOL"))

# Voeg EntrezID toe aan resultaten
resultaten$symbol <- rownames(resultaten)
resultaten <- merge(as.data.frame(resultaten), gene_info, by.x = "symbol", by.y = "SYMBOL")

# Verwijder genen zonder EntrezID of met dubbele ID's
resultaten <- resultaten[!is.na(resultaten$ENTREZID), ]
resultaten <- resultaten[!duplicated(resultaten$ENTREZID), ]

# Maak vector aan voor pathway analyse
gene_vector <- resultaten$log2FoldChange
names(gene_vector) <- resultaten$ENTREZID

# === 8. Pathway visualisatie met KEGG Pathview ===
pathview(
  gene.data = gene_vector,
  pathway.id = "hsa05323",   # KEGG pathway ID voor reumatoïde artritis
  species = "hsa",           # Humane soort
  gene.idtype = "ENTREZ",
  limit = list(gene = 5)     # Schaal: log2FC van -5 tot +5
)

# Bekijk structuur van resultaten
str(resultaten)
head(resultaten)
head(rownames(resultaten))

# === 9. GO-enrichment analyse met goseq ===
# Laad aanvullende packages
library(goseq)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)

# Selecteer significante differentieel geëxprimeerde genen (log2FC > 1, padj < 0.05)
significant_genes <- resultaten$padj < 0.05 & abs(resultaten$log2FoldChange) > 1

# Maak binaire vector van genen (1 = DEG, 0 = niet DEG)
gene_vector <- integer(nrow(resultaten))
names(gene_vector) <- resultaten$ENTREZID
gene_vector[significant_genes] <- 1

# Verwijder NA's of lege namen
gene_vector <- gene_vector[!is.na(names(gene_vector)) & names(gene_vector) != ""]

# Bepaal genlengtes met behulp van exoninformatie
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
exons <- exonsBy(txdb, by = "gene")
gene_lengths <- sapply(exons, function(x) sum(width(reduce(x))))
gene_lengths <- gene_lengths[names(gene_lengths) %in% names(gene_vector)]

# Zorg dat volgorde van genen in lengte en vector gelijk is
gene_vector <- gene_vector[names(gene_lengths)]

# Bereken Probability Weighting Function (PWF) voor bias-correctie
pwf <- nullp(gene_vector, bias.data = gene_lengths)

# Haal GO-annotaties op
all_genes <- names(gene_vector)
go_annotations <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = all_genes,
                                        columns = "GO",
                                        keytype = "ENTREZID")

# Filter alleen op GO termen binnen 'Biological Process' (BP)
go_annotations <- go_annotations[go_annotations$ONTOLOGY == "BP", ]

# Zet GO termen om in lijststructuur voor goseq
gene2cat <- split(go_annotations$GO, go_annotations$ENTREZID)

# Voer GO-analyse uit
goseq_results <- goseq(pwf, gene2cat = gene2cat)
head(goseq_results)

# === 10. Visualisatie van GO-resultaten ===
# Selecteer top 10 significante GO-termen
sig_go <- goseq_results %>%
  filter(over_represented_pvalue < 0.05) %>%
  arrange(over_represented_pvalue) %>%
  head(10)

# Maak barplot van top 10 GO-termen
ggplot(sig_go, aes(x = reorder(term, -over_represented_pvalue), y = -log10(over_represented_pvalue))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "GO term", y = "-log10(p-value)", title = "Top 10 overgerepresenteerde GO termen") +
  theme_minimal()

# === 11. Alternatieve visualisatie: Dotplot ===
top10 <- goseq_results %>%
  filter(!is.na(over_represented_pvalue)) %>%
  arrange(over_represented_pvalue) %>%
  head(10) %>%
  mutate(hitsPerc = numDEInCat * 100 / numInCat)

# Voeg GO termen toe
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
