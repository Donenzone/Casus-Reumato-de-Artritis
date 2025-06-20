# Casus Reumatoïde Artritis - Transcriptomics Analyse

## Repository Structuur

| Map / Bestand         | Inhoud                                                                                      |
|----------------------|---------------------------------------------------------------------------------------------|
| `Data/Raw/`          | Ruwe of originele inputdata (bijv. `Count_matrix.txt`)                                    |
| `Data/Processed/`    | Bewerkt of gefilterd materiaal (bijv. `Bewerkt_countmatrix.csv`)                           |
| `Scripts/`           | R-scripts die de analyse uitvoeren                                                         |
| `Results/`           | Visualisaties en outputbestanden, zoals `VolcanoplotWC.png`, `top10_GO.csv`                 |
| `Bronnen/`           | Wetenschappelijke literatuur en casusmateriaal                                            |
| `Assets/`            | Flowschema & Workflow verslag                                                |
| `Data_stewardship/`  | Documentatie over toepassing van de competentie *beheren*                                |
| `README.md`          | Dit bestand                                                                                |

---

## Inleiding  
Reumatoïde artritis (RA) is een systemische auto-immuunziekte waarbij het immuunsysteem het eigen lichaam aanvalt, vooral het synoviaal weefsel in gewrichten. Dit veroorzaakt synovitis, een ontsteking van het gewrichtsslijmvlies, wat leidt tot chronische pijn en schade. Hoewel de precieze oorzaak onbekend is, spelen genetische aanleg, omgevingsfactoren en een ontregeld immuunsysteem een belangrijke rol (Gabriel, 2001; Radu & Bungau, 2021). Vroege diagnose en behandeling zijn essentieel om schade te beperken.  

In deze casus is RNA-sequencing toegepast op synoviumbiopten van vier RA-patiënten (ACPA-positief) en vier gezonde controles. Het doel is genexpressieverschillen te analyseren, ziekte-gerelateerde genen te identificeren en biologische processen en pathways te koppelen aan RA. Transcriptomics verschaft inzicht in de moleculaire mechanismen van deze aandoening.  

De gebruikte bronnen zijn te vinden in de map [`Bronnen/`](Bronnen/).

---

## Methode  
De RNA-seq data van 8 samples (4 RA, 4 controles) zijn uitgelijnd op het humane referentiegenoom (GRCh38) met het R-pakket **Rsubread**. Met **featureCounts** zijn gen-tellingen gegenereerd op basis van een GTF-annotatiebestand.  

Na pre-processing is een differentiële genexpressie-analyse uitgevoerd met **DESeq2**. Laag-exprimerende genen zijn verwijderd, data genormaliseerd, en significante genen geïdentificeerd met criteria: adjusted p-waarde < 0.05 en |log2 fold change| > 1.  

De differentieel geëxprimeerde genen (DEGs) zijn geannoteerd en gebruikt voor functionele analyses:  
- KEGG-pathwayanalyse met **Pathview**  
- Gene Ontology (GO) verrijkingsanalyse met **goseq**, rekening houdend met genlengte-bias  

Het volledige script is beschikbaar in [`Scripts/`](Scripts/). Een overzicht van de workflow is te vinden in [`Workflow_Flowschema`](Assets/Workflow_Flowschema).

---

## Resultaten  
De analyse identificeerde meerdere genen met significant veranderde expressie tussen RA-patiënten en controles. Veel genen met hoge log2 fold changes zijn betrokken bij ontstekings- en immuunprocessen, conform bekende RA-pathologie.  

- De volcano plot is te bekijken in [`Volcanoplot.png`](Resultaten/Deseq2_results/Volcanoplot.png.png).  
- KEGG-pathway voor reumatoïde artritis (hsa05323) toont duidelijke expressieveranderingen.  
- Top 10 GO-processen, waaronder "immune response" en "inflammatory response", staan in [`GO_BP_resultaten.csv`](Resultaten/Pathway_analysis_results/GO_BP_resultaten.csv) en zijn gevisualiseerd in [`GOplot.png`](Results/GOplot.png).  

Deze bevindingen onderstrepen de kracht van transcriptomics om ziektegerelateerde mechanismen te identificeren.

---

## Conclusie  
De transcriptomics-analyse biedt waardevolle inzichten in de moleculaire basis van RA. Differentieel geëxprimeerde genen en pathways wijzen op belangrijke rollen van ontsteking en immuunregulatie.  

De resultaten bevestigen bestaande kennis en tonen de geschiktheid van transcriptomics voor het ontrafelen van ziekteprocessen. Toekomstig onderzoek kan zich richten op proteomics en grotere, longitudinale datasets.  

Dit project benadrukt tevens het belang van gestructureerd data- en scriptbeheer, zoals toegepast met GitHub, voor reproduceerbaarheid en transparantie.

---

## Data Stewardship & Beheer  
Voor een toelichting op het beheer van data en het toepassen van versiecontrole met GitHub, zie de bestanden in [`Data_stewardship/`](Data_stewardship/).

---

*Deze repository bevat alle data, scripts, en documentatie om de analyse te reproduceren en te begrijpen.*

---

**Vragen of opmerkingen?** Open gerust een issue of neem contact op via GitHub.

