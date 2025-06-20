Inhoud en Structuur van de Repository
Map / Bestand	Inhoud
data/raw/	Ruwe of originele inputdata, zoals count_matrix.txt
data/processed/	Bewerkt of gefilterd materiaal, bijvoorbeeld bewerkt_countmatrix.csv
scripts/	Alle R-scripts die de analyse uitvoeren
results/	Visualisaties zoals volcano plots, resultaten in CSV-formaat, GO-analyses
sources/	Wetenschappelijke literatuur en casusmateriaal
assets/	Afbeeldingen voor README of verslag
data_stewardship/	Beschrijving van de toepassing van de competentie beheren
README.md	Dit bestand met een overzicht van de inhoud en toelichting van het project

Inleiding
Reumatoïde artritis (RA) is een systemische auto-immuunziekte waarbij het immuunsysteem het eigen lichaam aanvalt, met name het synoviaal weefsel in de gewrichten. Dit leidt tot synovitis, een ontsteking van het gewrichtsslijmvlies, die op den duur chronische pijn en gewrichtsschade veroorzaakt. De precieze oorzaak van RA is nog niet volledig bekend, maar een combinatie van genetische aanleg, omgevingsfactoren en een dysfunctie van het immuunsysteem speelt een belangrijke rol (Gabriel, 2001; Radu & Bungau, 2021). Vroege diagnose en adequate behandeling zijn cruciaal om blijvende schade te beperken.

In deze casus is RNA-sequencing toegepast op synoviumbiopten van vier RA-patiënten (positief voor ACPA) en vier gezonde controles. Het doel is om verschillen in genexpressie tussen de groepen te analyseren, genen te identificeren die mogelijk bijdragen aan de ziekte, en de betrokken biologische processen en pathways in kaart te brengen. Met behulp van transcriptomics wordt zo meer inzicht verkregen in de moleculaire mechanismen achter RA.

De gebruikte wetenschappelijke bronnen en casusmateriaal zijn terug te vinden in de map /sources/.

Methode
De RNA-seq data van de 8 samples (4 RA, 4 controles) zijn eerst uitgelijnd op het humane referentiegenoom (GRCh38) met behulp van het R-pakket Rsubread. Vervolgens zijn gen-tellingen bepaald met featureCounts, gebruikmakend van het bijbehorende GTF-annotatiebestand.

Na de pre-processing is de differentiële genexpressie-analyse uitgevoerd met DESeq2. Hierbij werden laag-exprimerende genen verwijderd en de data genormaliseerd om verschillen in sequencingdiepte te corrigeren. Voor elk gen zijn log2 fold changes en p-waarden berekend, waarbij genen met een adjusted p-waarde < 0.05 en een absolute log2 fold change > 1 als significant zijn aangemerkt.

De significant differentieel geëxprimeerde genen (DEGs) zijn geannoteerd met EntrezID’s en gebruikt voor functionele analyses. Pathway-analyse is gedaan met KEGG Pathview, waarmee veranderingen in genexpressie binnen bekende biologische routes zijn gevisualiseerd. Daarnaast is een Gene Ontology (GO) verrijkingsanalyse uitgevoerd met het R-pakket goseq, dat bias door genlengte corrigeert.

Alle scripts voor de analyse, inclusief commentaar en workflow, zijn beschikbaar in de map /scripts/. Een overzichtelijke flowchart van de analyseworkflow is te vinden in /assets/Flowschema_Analyse.png.

Resultaten
De analyse toonde een aantal genen met significant verschillende expressie tussen RA-patiënten en gezonde controles. Diverse genen met hoge log2 fold changes zijn betrokken bij ontstekings- en immuunprocessen, wat overeenkomt met de bekende pathofysiologie van RA.

De volcano plot (opgeslagen als /results/VolcanoplotWC.png) visualiseert het onderscheid tussen genen met verhoogde en verlaagde expressie. De KEGG-pathway voor reumatoïde artritis (hsa05323) toont duidelijke expressieveranderingen in meerdere pathway-genen. GO-analyse identificeerde de top 10 verrijkte biologische processen, waaronder “immune response” en “inflammatory response”, welke zijn weergegeven in een dotplot (/results/GOplot.png) en in tabelvorm in /results/top10_GO.csv.

Deze bevindingen bevestigen dat transcriptomics een krachtig instrument is om ziektegerelateerde genen en mechanismen te ontrafelen.

Conclusie
De transcriptomics-analyse van synoviumweefsel van RA-patiënten levert belangrijke inzichten op in de moleculaire processen die aan deze ziekte ten grondslag liggen. Door vergelijking met gezonde controles zijn meerdere differentieel geëxprimeerde genen en relevante biologische pathways geïdentificeerd, vooral gerelateerd aan ontsteking en immuunregulatie.

Deze resultaten bevestigen bestaande kennis over RA en illustreren de waarde van transcriptomics voor het ontrafelen van ziekteprocessen op moleculair niveau. Voor toekomstig onderzoek wordt aanbevolen om de genexpressiegegevens te combineren met proteomics en klinische data, en de analyse uit te breiden naar grotere patiëntengroepen en longitudinale datasets (bijvoorbeeld pre- en post-behandeling).

Daarnaast benadrukt dit project het belang van zorgvuldig databeheer en versiecontrole, zoals geïmplementeerd met GitHub. Door een gestructureerde mappenindeling en documentatie is de volledige analyse reproduceerbaar en transparant.

