# Casus Reumatoïde Artritis – Transcriptomics Analyse

## Repository Structuur

| Map / Bestand       | Inhoud                                                                                  |
|---------------------|------------------------------------------------------------------------------------------|
| [Data/Raw/](Data/Raw/)             | Ruwe of originele inputdata (bijv. `Count_matrix.txt`)                                      |
| [Data/Processed/](Data/Processed/) | Bewerkt of gefilterd materiaal (bijv. `Bewerkt_countmatrix.csv`)                            |
| [Scripts/](Scripts/)               | R-scripts die de analyse uitvoeren                                                          |
| [Results/](Results/)               | Visualisaties en outputbestanden, zoals `VolcanoplotWC.png`, `top10_GO.csv`                |
| [Bronnen/](Bronnen/)               | Wetenschappelijke literatuur en casusmateriaal                                              |
| [Assets/](Assets/)                 | Flowschema & workflow-verslag                                                               |
| [Data_stewardship/](Data_stewardship/) | Documentatie over toepassing van de competentie *Beheren*                               |
| README.md               | Dit bestand                                                                                 |

---

## Inleiding

Reumatoïde artritis (RA) is een chronische auto-immuunziekte die vooral het gewrichtsslijmvlies aantast, wat leidt tot ontsteking, pijn en uiteindelijk gewrichtsschade. Deze aandoening treft wereldwijd miljoenen mensen en heeft een grote impact op kwaliteit van leven en maatschappelijke kosten. Hoewel de exacte oorzaak van RA onbekend is, wordt aangenomen dat een combinatie van genetische factoren, omgevingsinvloeden en een ontregeld immuunsysteem een rol speelt. In het synoviaal weefsel van aangetaste gewrichten zijn immuuncellen actief die ontstekingsmediatoren produceren en daardoor schade veroorzaken (Gabriel, 2001; Radu & Bungau, 2021). Vroege diagnostiek en behandeling zijn essentieel om progressie te beperken.

In deze studie is RNA-sequencing gebruikt om genexpressiepatronen in synoviumbiopten van vier RA-patiënten (ACPA-positief) en vier gezonde controles te onderzoeken. Door verschillen in genexpressie te analyseren kunnen genen en biologische processen geïdentificeerd worden die betrokken zijn bij het ziekteproces. Het doel is om een beter moleculair begrip van RA te krijgen en mogelijke biomarkers of therapeutische targets aan te wijzen. Deze transcriptomics-analyse kan zo de kennis over RA verbeteren en bijdragen aan gerichte behandelingen.

Alle gebruikte literatuur en data zijn opgenomen in de map [Bronnen/](Bronnen/).

---

## Methode

De RNA-seq data van acht samples (vier RA-patiënten en vier controles) zijn geanalyseerd volgens een gestandaardiseerd bioinformatica-protocol. Eerst werden de ruwe sequencelijnen uitgelijnd op het menselijke referentiegenoom GRCh38 met het R-pakket Rsubread. Vervolgens is met featureCounts het aantal reads per gen bepaald, gebruikmakend van een GTF-annotatiebestand.

Daarna is een kwaliteitscontrole uitgevoerd waarbij genen met lage expressie zijn verwijderd om ruis te minimaliseren. De genexpressiegegevens zijn genormaliseerd met het DESeq2-pakket om verschillen in sequencing diepte te corrigeren. Hierna is een differentiële genexpressie-analyse gedaan met de criteria: een adjusted p-waarde kleiner dan 0,05 en een absolute log2 fold change groter dan 1 om significante genen te identificeren.

De gevonden differentieel geëxprimeerde genen (DEGs) zijn geannoteerd en vervolgens onderworpen aan functionele analyses. Hiervoor is gebruikgemaakt van KEGG-pathwayanalyse via het pakket Pathview en Gene Ontology (GO) verrijkingsanalyse met goseq. Deze laatste houdt rekening met bias zoals genlengte.

Alle R-scripts voor deze stappen zijn te vinden in de map [Scripts/](Scripts/).  
Een visueel overzicht van de workflow is opgenomen in [Assets/Workflow_Flowschema](Assets
/Workflow_Flowschema).

---

## Resultaten

De differentiële genexpressie-analyse tussen synovium van RA-patiënten en gezonde controles identificeerde diverse genen met significante veranderingen in expressie. Veel van deze genen zijn betrokken bij immuun- en ontstekingsprocessen, zoals cytokines, chemokines en genen die een rol spelen in de activatie van immuuncellen. Dit sluit aan bij de bekende pathofysiologie van RA, waarin ontsteking centraal staat.

Een volcano plot visualiseert de spreiding van genen op basis van significatie en expressieverandering [Results/VolcanoplotWC.png](Resultaten/Deseq2_results/Volcanoplot.png.png).
erder werd een KEGG-pathwayanalyse uitgevoerd, waarbij onder andere de pathway ‘reumatoïde artritis’ (hsa05323) sterk verrijkt bleek. Deze resultaten zijn zichtbaar in [Results/hsa05323.pathview.png](Resultaten/Pathway_analysis_results/hsa05323.pathview.png).
Gene Ontology-analyse toonde verrijking in biologische processen zoals immuunrespons, ontsteking, en celactivatie. De top 10 GO-processen zijn gedocumenteerd in [Results/top10_GO.csv](Resultaten/Pathway_analysis_results/GO_BP_resultaten.csv) en visueel samengevat in [Results/GOplot.png](Resultaten/Pathway_analysis_results/GO_plot_zoom_png.png).
Deze bevindingen bevestigen het belang van immuun- en ontstekingsmechanismen bij RA en onderstrepen het nut van transcriptomics voor het ontrafelen van ziekteprocessen.

---

## Conclusie

Deze transcriptomics-analyse van synoviumweefsel van RA-patiënten geeft waardevolle inzichten in de moleculaire veranderingen die bijdragen aan de ziekte. De significant differentieel geëxprimeerde genen wijzen duidelijk op een centrale rol voor immuunactivatie en ontsteking. De gevonden pathways en biologische processen sluiten goed aan bij de huidige kennis over RA-pathogenese.

Deze studie bevestigt dat RNA-sequencing een krachtig hulpmiddel is om complexe ziekteprocessen te bestuderen en potentiële biomarkers te identificeren. De data benadrukken het belang van vroege ontstekingsremming en kunnen richting geven aan de ontwikkeling van gerichte therapieën.

Voor toekomstig onderzoek wordt aanbevolen om deze transcriptomics-resultaten te combineren met proteomics en grotere, longitudinale datasets om ziekteprogressie beter te begrijpen. Daarnaast is gestructureerd beheer van data en scripts, bijvoorbeeld via GitHub, essentieel voor transparantie en reproduceerbaarheid van onderzoek.

Dit project illustreert niet alleen de wetenschappelijke waarde van transcriptomics bij RA, maar ook het belang van goede data stewardship en versiebeheer voor betrouwbaar wetenschappelijk onderzoek.

---

## Data Stewardship & Beheer

Voor een toelichting op het beheer van data en het toepassen van versiecontrole met GitHub, zie de bestanden in [Data_stewardship/](Data_stewardship/).

Deze repository bevat alle data, scripts, en documentatie om de analyse te reproduceren en te begrijpen.

---

## Vragen of opmerkingen?

Open gerust een issue of neem contact op via GitHub.
