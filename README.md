# Casus Reumatoïde Artritis - Transcriptomics Analyse

## Repository Structuur

| Map / Bestand          | Inhoud                                                                                      |
|------------------------|---------------------------------------------------------------------------------------------|
| [`Data/Raw/`](Data/Raw/)              | Ruwe of originele inputdata (bijv. `Count_matrix.txt`)                                                |
| [`Data/Processed/`](Data/Processed/)  | Bewerkt of gefilterd materiaal (bijv. `Bewerkt_countmatrix.csv`)                                       |
| [`Scripts/`](Scripts/)                | R-scripts die de analyse uitvoeren                                                                     |
| [`Results/`](Results/)                | Visualisaties en outputbestanden, zoals `Volcanoplot.png`, `GO_BP_resultaten.csv`                      |
| [`Bronnen/`](Bronnen/)                | Wetenschappelijke literatuur en casusmateriaal                                                         |
| [`Assets/`](Assets/)                  | Flowschema & Workflow verslag                                                                          |
| [`Data_stewardship/`](Data_stewardship/) | Documentatie over toepassing van de competentie *beheren*                                           |
| [`README.md`](README.md)              | Dit bestand                                                                                             |


---

## Inleiding

Reumatoïde artritis (RA) is een systemische auto-immuunziekte waarbij het immuunsysteem het eigen lichaam aanvalt, in het bijzonder het synoviaal weefsel in gewrichten. Dit leidt tot synovitis, een ontsteking van het gewrichtsslijmvlies, en veroorzaakt chronische pijn en gewrichtsschade. Hoewel de precieze oorzaak onbekend is, spelen genetische factoren, omgevingsinvloeden en een ontregeld immuunsysteem een rol (Gabriel, 2001; Radu & Bungau, 2021). Vroege diagnose en behandeling zijn essentieel om onomkeerbare schade te beperken.

In deze casus is RNA-sequencing uitgevoerd op synoviumbiopten van vier ACPA-positieve RA-patiënten en vier gezonde controles. Het doel is om genexpressieverschillen te analyseren, ziekte-gerelateerde genen te identificeren, en biologische processen en pathways te koppelen aan RA met behulp van transcriptomics.

De gebruikte bronnen zijn te vinden in de map [`Bronnen/`](Bronnen/).

---

## Methode

De RNA-seq data van 8 samples (4 RA, 4 controles) zijn eerst uitgelijnd op het humane referentiegenoom (GRCh38) met behulp van het `Rsubread`-pakket in R. Gen-tellingen zijn gegenereerd met `featureCounts`, op basis van een GTF-annotatiebestand.

Na filtering van laag-exprimerende genen is een differentiële genexpressie-analyse uitgevoerd met `DESeq2`. Data zijn genormaliseerd en significante genen zijn geselecteerd op basis van een aangepaste p-waarde < 0.05 en |log2 fold change| > 1.

De differentieel geëxprimeerde genen (DEGs) zijn geannoteerd en functioneel geanalyseerd:

- KEGG pathway-analyse met `pathview`
- GO-verrijkingsanalyse met `goseq` (inclusief correctie voor genlengte-bias)

Het volledige script staat in [`Scripts/`](Scripts/).  
Een overzicht van de workflow is te vinden in [`Assets/Workflow_Flowschema.png`](Assets/Workflow_Flowschema.png).

---

## Resultaten

De analyse identificeerde meerdere genen met significant veranderde expressie tussen RA-patiënten en gezonde controles. Genen met sterke expressieveranderingen zijn betrokken bij ontstekings- en immuunprocessen, passend bij de bekende pathofysiologie van RA.

Belangrijkste bevindingen:

- Volcano plot: [`Results/Volcanoplot.png`](Results/Volcanoplot.png)
- KEGG pathway `hsa05323` (RA): visualisatie in [`Results/hsa05323.pathview.png`](Results/hsa05323.pathview.png)
- Top 10 verrijkte GO-biologische processen, waaronder "immune response" en "inflammatory response", in [`Results/GO_BP_resultaten.csv`](Results/GO_BP_resultaten.csv) en [`Results/GO_plot_zoom_png.png`](Results/GO_plot_zoom_png.png)

Deze resultaten illustreren hoe transcriptomics belangrijke ziekteprocessen kan blootleggen.

---

## Conclusie

De transcriptomics-analyse biedt inzicht in de moleculaire mechanismen van RA. De geïdentificeerde differentieel geëxprimeerde genen en verrijkte pathways ondersteunen de centrale rol van ontsteking en immuunactivatie.

De resultaten bevestigen eerder gepubliceerde bevindingen en illustreren hoe RNA-sequencing een krachtig hulpmiddel is voor het identificeren van biomarkers en therapeutische targets. Verder onderzoek kan zich richten op integratie met proteomics, validatie-experimenten en grotere patiëntcohorten.

De opbouw en structuur van dit project via GitHub draagt bij aan reproduceerbaarheid en transparantie van het onderzoek.

---

## Data Stewardship & Beheer

Deze repository is opgezet met aandacht voor goed data- en scriptbeheer. Voor uitleg over hoe de competentie *Beheren* is toegepast, zie de map [`Data_stewardship/`](Data_stewardship/). Hierin zijn documenten opgenomen over:

- het organiseren van data- en scripts volgens FAIR-principes
- het toepassen van versiebeheer met GitHub
- het vastleggen van projectinformatie op een reproduceerbare manier

Deze README en de bijbehorende folderstructuur vormen samen een reproduceerbare analyse-omgeving waarin transparantie, vindbaarheid en herbruikbaarheid centraal staan.

---

## Vragen of opmerkingen?

Open een issue of neem contact op via deze GitHub-pagina.
