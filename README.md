Inhoud / structuur
Map / Bestand	Inhoud
data/raw/	Ruwe of originele inputdata (zoals count_matrix.txt)
data/processed/	Bewerkt of gefilterd materiaal zoals bewerkt_countmatrix.csv
scripts/	Alle R-scripts die de analyse uitvoeren
resultaten/	Volcano plot, resultaten CSV’s, GO-analyses
bronnen/	Wetenschappelijke literatuur en casusmateriaal
assets/	Afbeeldingen voor README of verslag
data_stewardship/	Beschrijving van hoe je de competentie beheren hebt toegepast
README.md	Dit bestand



1. Inleiding

Reumatoïde artritis (RA) is een systemische auto-immuunziekte waarbij het immuunsysteem het eigen lichaam aanvalt, met name het synoviaal weefsel van gewrichten. Kenmerkend is synovitis: een ontsteking van het gewrichtsslijmvlies die leidt tot chronische pijn en gewrichtsschade. Hoewel de exacte oorzaak onbekend is, speelt een combinatie van genetische aanleg, omgevingsfactoren en een ontspoord immuunsysteem een rol (Gabriel, 2001; Radu & Bungau, 2021). Vroege diagnose en behandeling zijn essentieel om schade te beperken.

In deze casus is RNA-sequencing toegepast op synoviumbiopten van vier RA-patiënten (positief voor ACPA) en vier gezonde controles. Het doel is om genexpressieverschillen te analyseren, betrokken genen te identificeren, en biologische processen en pathways te koppelen aan RA. Door gebruik te maken van transcriptomics-technieken wordt inzicht verkregen in de onderliggende moleculaire mechanismen.

De gebruikte bronnen zijn terug te vinden in de map [`/bronnen`](./bronnen), waaronder wetenschappelijke artikelen en handleidingen uit de werkcolleges.


 2. Methode

De analyse werd uitgevoerd in R. Allereerst werden de sequencing reads van de 8 samples (4 RA, 4 controles) uitgelijnd op het humane referentiegenoom (GRCh38) met behulp van het `Rsubread`-pakket. Daarna werden tellingen per gen gegenereerd met `featureCounts`, waarbij het bijbehorende GTF-annotatiebestand werd gebruikt.

Na pre-processing is een differentiële genexpressie-analyse uitgevoerd met het `DESeq2`-pakket. Hierbij werd voor elk gen het verschil in expressie tussen de RA-groep en controle berekend, inclusief log2FoldChange en p-waarden. Genen met een p-adj < 0.05 en |log2FC| > 1 werden als significant beschouwd.

Vervolgens werd pathway-analyse uitgevoerd met het `pathview`-pakket, waarbij log2FoldChange-waarden werden gekoppeld aan KEGG-pathways. Daarnaast werd met `goseq` een Gene Ontology (GO) analyse gedaan om overgerepresenteerde biologische processen te identificeren.

Scripts, data en outputs zijn ondergebracht in de mappen `/scripts`, `/data`, en `/resultaten`, zodat de volledige analyse reproduceerbaar is. Zie ook het bestand [`data_stewardship.md`](./data_stewardship/data_stewardship.md) voor uitleg over databeheer.


3. Resultaten

De differentiële genexpressie-analyse leverde meerdere genen op met een significant verschil tussen RA-patiënten en gezonde personen. Enkele genen met de hoogste log2FoldChange-waarden zijn betrokken bij ontstekingsprocessen en immuunregulatie, wat in lijn is met de bekende pathologie van RA. De resultaten zijn opgeslagen in [`ResultatenWC3.csv`](./resultaten/ResultatenWC3.csv).

De volcano plot (zie [`VolcanoplotWC.png`](./resultaten/VolcanoplotWC.png)) toont de verdeling van log2FoldChange versus de p-waarde. Hieruit blijkt dat er genen zijn met sterk verhoogde of verlaagde expressie bij RA.

Via `pathview` werd de KEGG-pathway voor reumatoïde artritis (hsa05323) gevisualiseerd. Meerdere genen in deze pathway vertoonden duidelijke expressieveranderingen. In de GO-analyse (via `goseq`) werden de top 10 overgerepresenteerde processen geïdentificeerd, waaronder "immune response" en "inflammatory response". Deze resultaten zijn gevisualiseerd in een dotplot (`GOplot.png`) en in tabelformaat opgeslagen in `top10_GO.csv`.

Deze bevindingen ondersteunen het gebruik van transcriptomics om ziektegerelateerde genen en mechanismen te identificeren.



4. Conclusie

De transcriptomics-analyse van synoviumweefsel heeft waardevolle inzichten opgeleverd in de moleculaire processen die een rol spelen bij reumatoïde artritis. Door het vergelijken van genexpressie tussen RA-patiënten en gezonde controles zijn meerdere differentieel geëxprimeerde genen en bijbehorende pathways geïdentificeerd. Met name ontstekingsgerelateerde genen en immuunprocessen kwamen sterk naar voren in de analyses.

Deze resultaten bevestigen bekende kenmerken van RA en laten zien dat transcriptomics geschikt is om ziekteprocessen op cellulair niveau te ontrafelen. Daarnaast tonen de KEGG- en GO-analyses aan dat betrokken genen zich groeperen in biologisch relevante netwerken.

Aanbevolen wordt om vervolgonderzoek te doen op eiwitniveau (proteomics) om te bevestigen welke genen functioneel actief zijn in RA. Ook zou uitbreiding naar grotere datasets of tijdreeksen (bijv. voor en na behandeling) nog meer inzicht kunnen bieden.

Tot slot toont dit project het belang aan van zorgvuldig data- en scriptbeheer: dankzij de heldere structuur en versiecontrole via GitHub is de gehele analyse reproduceerbaar en inzichtelijk voor anderen.




