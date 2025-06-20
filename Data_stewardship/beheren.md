 ============================================
 Efficiënt en gestructureerd databeheer
 ============================================

 Het belang van databeheer

 - Reproduceerbaarheid van het onderzoek, doordat anderen (of jijzelf op een later moment)
   exact kunnen nagaan welke stappen zijn genomen.
 - Voorkomen van dataverlies of verwarring, doordat alles logisch en eenduidig wordt opgeslagen.
 - Soepelere samenwerking tussen teamleden, met duidelijke afspraken over bestandslocaties en -namen.
 - Waarborging van de integriteit van resultaten, doordat transparant is hoe analyses zijn uitgevoerd.

 Onsamenhangende of ongeorganiseerde werkwijze verhoogt de kans op fouten, duplicatie van bestanden
 en onherstelbaar dataverlies.

 --------------------------------------------
 Mappenstructuur
 --------------------------------------------

 Tijdens het project is de volgende mappenstructuur gebruikt:

 /data
 ├── raw/                Originele, onbewerkte bestanden
 │   ├── fastq/
 │   ├── count_matrices/
 │   └── other/
 ├── processed/          Bewerkte data of output van analyses
 │   ├── deseq2_results/
 │   └── pathway_analysis/
 ├── metadata/           Sample-informatie en experimentele details
 │   ├── samples/
 │   └── experiment_info/
 /scripts/               Alle R-scripts die analyses uitvoeren
 /figures/               Grafieken en visualisaties
 /docs/                  Verslaglegging of notities
 README.md               Projectoverzicht en toelichting

 Deze structuur helpt overzicht te behouden en voorkomt dat data of resultaten per ongeluk
 worden overschreven of verplaatst.

 --------------------------------------------
 Naamgeving van bestanden
 --------------------------------------------

 - Gebruik van kleine letters (lowercase)
 - Gebruik underscores (_) in plaats van spaties
 - Korte maar duidelijke beschrijvingen, bv:
   count_matrix_raw.txt, volcanoplot_RA_vs_Control.png

 Dit voorkomt verwarring en maakt bestanden makkelijk vindbaar, zeker voor geautomatiseerde scripts.

 --------------------------------------------
 Versiebeheer (Git & GitHub)
 --------------------------------------------

 - Wijzigingen aan scripts of data systematisch bijhouden
 - Teruggaan naar eerdere versies bij fouten
 - Samenwerken zonder overschrijven van elkaars werk
 - Aantonen wanneer en waarom aanpassingen zijn gedaan

 Voorbeelden van commits tijdens dit project:
   - Eerste upload van ruwe data
   - Toevoegen van DESeq2-analyse
   - Genereren van volcano plot en export als PNG

 --------------------------------------------
 Documentatie in scripts (voorbeeld R)
 --------------------------------------------
 === 8. Pathway visualisatie met KEGG Pathview ===
 pathview(
   gene.data = gene_vector,
   pathway.id = "hsa05323",     KEGG pathway ID voor reumatoïde artritis
   species = "hsa",             Humane soort
   gene.idtype = "ENTREZ",
   limit = list(gene = 5)       Schaal: log2FC van -5 tot +5
 )

 Goede documentatie maakt het mogelijk dat anderen (of jijzelf) de code begrijpen,
 aanpassen en hergebruiken. Het voorkomt fouten door onduidelijkheid of verkeerde aannames.

 --------------------------------------------
 Veilig omgaan met gevoelige data
 --------------------------------------------

 - Anonimiseren of pseudonimiseren van persoonsgegevens
 - Opslaan op versleutelde schijven of beveiligde servers
 - Gebruik van gesloten repositories en toegangsbeheer
 - Naleving van AVG (Algemene Verordening Gegevensbescherming)

 Zorgvuldige omgang met gevoelige data is wettelijk en ethisch verplicht.

 --------------------------------------------
 Open data en het publiceren van datasets
 --------------------------------------------

 - Maakt controle en herhaling van analyses mogelijk
 - Stimuleert hergebruik en vervolgonderzoek
 - Vergemakkelijkt samenwerking binnen en buiten het vakgebied

 Bij publicatie is het essentieel dat data goed gestructureerd, gedocumenteerd en begrijpelijk is.
