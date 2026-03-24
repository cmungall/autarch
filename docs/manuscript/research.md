
============================================================
---
provider: falcon
model: Edison Scientific Literature
cached: false
start_time: '2026-03-10T20:39:09.305386'
end_time: '2026-03-10T21:08:11.240986'
duration_seconds: 1741.94
provider_config:
  timeout: null
  max_retries: 3
  parameters:
    allowed_domains: []
    temperature: 0.1
citation_count: 52
---

## Question

research the history of reaction and enzyme classification, covering the enzyne commission through to GO and databases like MetaCyc; is enzyme classification important? who uses. Through to modern databases. Also any approaches to automate reaction classification

## Output

Question: You are an expert researcher providing comprehensive, well-cited information.

Provide detailed information focusing on:
1. Key concepts and definitions with current understanding
2. Recent developments and latest research (prioritize 2023-2024 sources)
3. Current applications and real-world implementations
4. Expert opinions and analysis from authoritative sources
5. Relevant statistics and data from recent studies

Format as a comprehensive research report with proper citations. Include URLs and publication dates where available.
Always prioritize recent, authoritative sources and provide specific citations for all major claims.

research the history of reaction and enzyme classification, covering the enzyne commission through to GO and databases like MetaCyc; is enzyme classification important? who uses. Through to modern databases. Also any approaches to automate reaction classification

# Research report: History and importance of enzyme/reaction classification (EC → GO → modern databases) and automation of reaction classification

## Executive summary
Enzyme/reaction classification began as a response to ambiguous and inconsistent enzyme naming and evolved into the Enzyme Commission (EC) system, a stable reaction-centric identifier scheme that remains foundational for curation, genome annotation, metabolic modeling, and pathway engineering. (tipton2000historyofthe pages 1-2, tipton2000historyofthe pages 2-3, tipton2000historyofthe pages 3-5)

Modern resources increasingly separate (i) **reaction identity** (e.g., Rhea) from (ii) **protein entities** (UniProt) and (iii) **pathway context** (Reactome/MetaCyc/EcoCyc), and they interoperate through controlled identifiers such as **EC**, **GO molecular function**, **Rhea reaction IDs**, and **ChEBI** chemical entities. (bansal2022rheathereaction pages 1-2, bateman2023uniprottheuniversal pages 2-3, milacic2024thereactomepathway pages 5-6)

Since 2023–2024, enzyme/reaction classification has also become a central target for **automation** using transformers and other ML approaches that predict EC numbers from protein sequences or directly from reaction SMILES, and for “reaction-to-enzyme” tools used in synthetic biology and biocatalysis. (kim2023functionalannotationof pages 1-2, qian2024ageneralmodel pages 1-2, shi2024remeanintegrated pages 1-2)

| Period / Year | Milestone | Impact & Significance | Key Source |
|---|---|---|---|
| **1833** | Discovery of diastase & naming convention | Payen and Persoz identified diastase and proposed the suffix '-ase' for enzymes. | (tipton2000historyofthe pages 1-2) |
| **1876** | Term 'Enzyme' coined | W.F. Kühne introduced the term 'enzyme'. | (tipton2000historyofthe pages 1-2) |
| **1956** | International Commission on Enzymes | Established by Marcel Florkin/IUB to standardize nomenclature amidst confusion from trivial names. | (mcdonald2014fifty‐fiveyearsof pages 2-4, tipton2000historyofthe pages 3-5) |
| **1961** | First IUB Report (Enzyme List) | Publication of the first official Enzyme List with the 4-digit EC number system (Class.Subclass.Sub-subclass.Serial). | (mcdonald2014fifty‐fiveyearsof pages 2-4, tipton2000historyofthe pages 2-3) |
| **1965–1992** | Printed Editions | Successive updates (1965, 1972, 1978, 1984, 1992) managed the growing list, eventually exceeding print limitations. | (mcdonald2014fifty‐fiveyearsof pages 2-4, tipton2000historyofthe pages 6-7) |
| **2000s** | Shift to Online (ExplorEnz) | Moved from print to MySQL-based web curation (ExplorEnz) to handle rapid data growth and formatting. | (mcdonald2014fifty‐fiveyearsof pages 2-4) |
| **2019–2021** | UniProt adopts Rhea | UniProtKB adopted Rhea as the reference vocabulary for enzyme annotation; reaction coverage grew from 6,654 to 9,294 (68% of Rhea). | (bansal2022rheathereaction pages 1-2) |
| **2023** | EcoCyc & Flux Models | EcoCyc update highlighted its role in generating steady-state metabolic flux models and organizing data via sophisticated ontology. | (karp2023theecocycdatabase pages 1-2) |
| **2023** | DeepECtransformer | Transformer-based model predicted EC numbers for 464 unannotated *E. coli* genes; outperformed previous methods. | (kim2023functionalannotationof pages 1-2) |
| **2024** | BEC-Pred (ML) | BERT-based model for EC prediction from reaction SMILES achieved 91.6% accuracy, outperforming sequence/graph methods. | (qian2024ageneralmodel pages 1-2) |
| **2024** | ReactZyme | Framed enzyme-reaction prediction as a retrieval task to annotate enzymes by reaction rather than family. | (hua2024reactzymeabenchmark pages 1-2) |


*Table: A historical overview of enzyme classification milestones from the 19th century to modern AI-driven annotation and database integration.*

| Resource (URL) | Scope & Curation Focus | Primary Identifiers | Recent Statistics (2021–2024) | Typical Applications |
|---|---|---|---|---|
| **UniProt**<br>(uniprot.org) | Global protein sequence & functional annotation. | Rhea (reactions), EC, GO, ChEBI. | 251M total records; ~10,540 Rhea reactions annotating ~24.8M sequences (2023). (bateman2023uniprottheuniversal pages 2-3, xie2024deeplearningin pages 2-3) | Sequence annotation, functional prediction, genome annotation. |
| **Rhea**<br>(rhea-db.org) | Expert-curated biochemical reactions; UniProt's reference vocabulary. | Rhea IDs, ChEBI, EC, GO, UniProt. | 13,673 unique reactions; 11,886 participants (2021/22). (bansal2022rheathereaction pages 1-1, bansal2022rheathereaction pages 5-6) | Reaction standards, metabolic modeling, enzyme annotation. |
| **Reactome**<br>(reactome.org) | Human pathways & processes as ordered molecular transformations. | Reactome IDs, Rhea, ChEBI, UniProt. | ~16,002 reactions; ~2,600 pathways; 11,630 human genes (2024/26). (ragueneau2026thereactomeknowledgebase pages 1-2, presern2023enzymedatabasesin pages 9-11) | Pathway enrichment, variant effect analysis, drug targeting. |
| **MetaCyc**<br>(metacyc.org) | Curated metabolic pathways from all life domains; experimental focus. | EC, M-numbers, BioCyc IDs. | ~3,153 pathways; ~19,020 reactions (2024). (xie2024deeplearningin pages 2-3, presern2023enzymedatabasesin pages 9-11) | Metabolic reconstruction (PathoLogic), reference knowledgebase. |
| **KEGG**<br>(genome.jp/kegg) | Integrated database of genomes, pathways, and biological chemicals. | KO (Orthology), R-numbers, EC. | >20,000 reactions (2024); ~12,000 reactions (2023). (xie2024deeplearningin pages 2-3, presern2023enzymedatabasesin pages 9-11) | Pathway mapping, orthology annotation, de novo pathway design. |
| **EcoCyc**<br>(ecocyc.org) | Model organism (*E. coli* K-12) genome & machinery. | BioCyc IDs, EC, GO. | Gold-standard curation; generates steady-state flux models for every release (2023). (karp2023theecocycdatabase pages 1-2) | Metabolic engineering, flux balance analysis, operon prediction. |
| **BRENDA**<br>(brenda-enzymes.org) | Comprehensive enzyme functional & kinetic data. | EC numbers, UniProt. | ~177,000 $K_m$ and ~87,000 $k_{cat}$ values; ~157,000 raw kinetic entries (2023/25). (wei2025findingthedark pages 7-8, presern2023enzymedatabasesin pages 8-9) | Kinetic modeling, enzyme characterization, biocatalysis. |
| **SABIO-RK**<br>(sabio.h-its.org) | Structured kinetic data with experimental conditions. | Reactome/KEGG/UniProt links. | ~51,000 $K_m$ and ~28,000 $k_{cat}$ entries (2023). (presern2023enzymedatabasesin pages 8-9) | Systems biology modeling, precise kinetic analysis. |


*Table: A comparison of key databases for enzyme and reaction data, highlighting their identifiers, scale, and primary use cases based on 2021–2024 literature.*

## 1) Key concepts and definitions (current understanding)

### 1.1 What is being classified: enzymes vs reactions
The EC system classifies **enzyme-catalyzed reactions** (i.e., the overall chemical transformation), not protein families or mechanisms; consequently, distinct proteins (including isoenzymes) can share an EC number if they catalyze the same overall reaction. (tipton2000historyofthe pages 3-5, mcdonald2014fifty‐fiveyearsof pages 7-9)

A persistent limitation of a purely reaction-based system is that it may group enzymes with different mechanisms and sequences under one EC number and can obscure multifunctionality or alternative reactions. (tipton2000historyofthe pages 3-5, mcdonald2014fifty‐fiveyearsof pages 6-7)

### 1.2 EC numbers
The EC system uses a **four-level hierarchical code** (e.g., 3.1.1.1): first digit = broad class; second/third digits provide finer reaction detail; fourth digit encodes substrate specificity/serial identity within a sub-subclass. (tipton2000historyofthe pages 3-5, qian2024ageneralmodel pages 2-3)

EC numbers are non-reusable; misclassifications and obsoleted activities remain as deleted/transferred entries, supporting stability but leaving historical gaps. (tipton2000historyofthe pages 3-5, mcdonald2014fifty‐fiveyearsof pages 7-9)

### 1.3 Gene Ontology (GO) molecular function for catalytic activity
GO provides an ontology-based representation of gene/protein functions. For enzymes, GO’s molecular function branch contains catalytic activity terms, and its organization is reported to “largely follow” the EC hierarchy. (garapati2026comprehensiveannotationof pages 1-2)

GO has recently undertaken explicit alignment of its **catalytic activity branch (GO:0003824)** with **EC** and **Rhea** to improve precision and interoperability of catalytic function terms with reaction definitions. (g2026thegeneontology pages 5-7)

### 1.4 Reaction knowledgebases and chemical ontologies (Rhea and ChEBI)
Rhea is an expert-curated biochemical reaction knowledgebase that represents reaction participants using **ChEBI** and serves as the reference vocabulary for enzyme/transporter reaction annotation in UniProtKB. (bansal2022rheathereaction pages 1-2)

Rhea provides cross-references and identifier mapping across many resources (KEGG/Reactome/MetaCyc/EcoCyc; EC; GO molecular function), enabling reaction normalization across databases. (bansal2022rheathereaction pages 5-6)

## 2) History of enzyme nomenclature and classification

### 2.1 Early nomenclature milestones
The “-ase” naming convention traces back to Payen and Persoz’s work on diastase (1833), and the term “enzyme” was coined by Kühne (1876). (tipton2000historyofthe pages 1-2)

### 2.2 Drivers of standardization
Enzyme names historically suffered from ambiguity (multiple names per enzyme; same name used for different enzymes; arcane/trivial names), motivating formal nomenclature and classification efforts. (tipton2000historyofthe pages 1-2, mcdonald2014fifty‐fiveyearsof pages 1-2)

### 2.3 Establishment of the EC system
The International Commission on Enzymes was established in 1956 to address nomenclature/classification, leading to an official enzyme list published in 1961 and successive printed updates through 1992. (mcdonald2014fifty‐fiveyearsof pages 2-4, tipton2000historyofthe pages 3-5)

The field’s shift from printed lists to online systems (e.g., ExplorEnz) was driven by rapid growth and impracticality of maintaining printed volumes (e.g., >2000 A4 pages), as well as the need for richer search and curation workflows. (mcdonald2014fifty‐fiveyearsof pages 2-4)

### 2.4 Persistent conceptual tensions
Key limitations recognized over decades include: grouping by overall reaction despite mechanistic diversity; lack of reaction direction; challenges in representing multi-step or multifunctional enzymes; and constraints that EC assignment requires direct experimental evidence (sequence similarity alone is insufficient). (tipton2000historyofthe pages 3-5, mcdonald2014fifty‐fiveyearsof pages 7-9)

## 3) Modern databases and how classification is used today (2023–2024 emphasis)

### 3.1 UniProtKB (2023) and standardized reaction annotation
UniProtKB now includes large-scale Rhea-based reaction annotations: **10,540 Rhea reactions** are linked to **24,842,646 UniProtKB sequence records** (including **226,101 reviewed Swiss-Prot** records). (bateman2023uniprottheuniversal pages 2-3)

UniProt also deploys automated annotation systems that can include catalytic activity and EC predictions, such as ARBA (rules-based, self-training) producing 27,338 rules (release 2022_03), aiming to reduce “uncharacterized protein” labels at scale. (bateman2023uniprottheuniversal pages 4-5)

### 3.2 Rhea as an interoperability hub (reaction identity)
Rhea is positioned explicitly as a bridge between UniProt (proteins) and ChEBI (small molecules), and as a standard vocabulary for UniProtKB enzyme/transporter reaction annotations. (bansal2022rheathereaction pages 5-6, bansal2022rheathereaction pages 1-2)

Quantitative interoperability/coverage signals include (i) growth of UniProt-linked Rhea reactions from 6,654 (Oct 2019) to 9,294 (June 2021), about **68%** of Rhea at that time, and (ii) Rhea release 119 (June 2021) containing **13,673 unique reactions**, **11,886 unique participants**, and **15,500 literature references**. (bansal2022rheathereaction pages 1-2, bansal2022rheathereaction pages 1-1)

A Rhea reaction page can cross-reference UniProtKB proteins, EC numbers, and GO terms, illustrating how reaction, enzyme classification, and ontology annotations are connected in practice. (bansal2022rheathereaction media 1bea24f1)

### 3.3 Reactome (2024): pathways as ordered reaction networks
Reactome describes human biology as an **ordered network of molecular transformations**, serving both as an archive and as a tool for analyzing gene expression and somatic mutation catalogs (e.g., tumor mutation sets). (milacic2024thereactomepathway pages 1-2)

Reactome integrates external standards: it uses UniProt proteins and ChEBI small molecules and takes balanced chemical reaction equations from Rhea; GO alignment is used for controlled vocabularies to support interoperability. (milacic2024thereactomepathway pages 5-6)

### 3.4 EcoCyc (2023) and BioCyc/MetaCyc ecosystem
EcoCyc is a curated knowledgebase for *E. coli* K-12 MG1655, with pages for gene products, metabolites, reactions, operons, and pathways. Each release generates an executable steady-state metabolic flux model to predict fluxes and growth under gene knockouts and nutrient conditions, supporting both experimentalists and modelers. (karp2023theecocycdatabase pages 1-2)

MetaCyc, a curated multi-organism metabolic pathway resource, is reported (in 2023–2024 reviews) at roughly **~3,100–3,153 pathways** and **~18,500–19,020 reactions**; KEGG and Reactome are reported at **~12,000 to >20,000 reactions** (KEGG) and **~15,000 reactions** (Reactome). (presern2023enzymedatabasesin pages 9-11, xie2024deeplearningin pages 2-3)

### 3.5 Kinetic and enzyme-property databases (BRENDA; SABIO-RK)
BRENDA and SABIO-RK illustrate why EC/GO/reaction identifiers matter beyond naming: kinetic data must connect a protein to a reaction/substrate context.

A 2023 enzyme-database review reports approximate kinetic entry counts: BRENDA with **~177,000 Km** and **~87,000 kcat** values; SABIO-RK with **~51,000 Km** and **~28,000 kcat** values, with SABIO-RK emphasizing structured rate laws and experimental conditions and using ontologies (ChEBI, GO, etc.) for standardization. (presern2023enzymedatabasesin pages 8-9)

## 4) Is enzyme classification important? Who uses it?

### 4.1 Why EC/GO/reaction identifiers remain essential
Despite known limitations, EC numbers provide stable, widely recognized identifiers that support unambiguous linking across resources and publications and provide a scaffold for curated updates and public review. (tipton2000historyofthe pages 3-5, mcdonald2014fifty‐fiveyearsof pages 7-9)

Modern reaction normalization (Rhea+ChEBI) and pathway knowledgebases (Reactome, MetaCyc/EcoCyc) depend on stable function/reaction identifiers to connect proteins, chemicals, and pathways into computable networks. (milacic2024thereactomepathway pages 5-6, bansal2022rheathereaction pages 5-6)

A key practical problem remains coverage: many reactions lack EC assignment because full enzyme characterization is required, leaving sparse and imbalanced labeled datasets—especially in secondary metabolism—motivating both curation and automation efforts. (qian2024ageneralmodel pages 2-3)

### 4.2 Primary user communities and real-world implementations
* **Biocurators and database teams** rely on EC/GO/Rhea mappings to curate standardized knowledge and to maintain interoperability across UniProt, GO, Reactome, and pathway databases. (bansal2022rheathereaction pages 1-2, milacic2024thereactomepathway pages 5-6)
* **Genome annotation pipelines** use EC and GO as core function labels; GO’s TreeGrafter/InterProScan pipeline is used for computational GO inference in UniProt proteomes and NCBI RefSeq eukaryotic genomes. (g2026thegeneontology pages 5-7)
* **Systems biologists and metabolic modelers** use curated reaction/pathway knowledgebases; EcoCyc explicitly provides executable flux models per release and is used as a gold standard for method development. (karp2023theecocycdatabase pages 1-2)
* **Biocatalysis and synthetic biology practitioners** use reaction-to-enzyme mining tools that rely on reaction similarity, EC expansion, and kinetic/condition prediction to choose candidate enzymes for non-natural reactions. (shi2024remeanintegrated pages 1-2, shi2024remeanintegrated pages 2-4)

## 5) Recent developments (prioritizing 2023–2024)

### 5.1 Alignment of GO catalytic activity with EC/Rhea (ontology modernization)
GO reports a concerted alignment of catalytic activity terms with EC and Rhea, including updates to >2000 cross-references, obsoletion of >1000 redundant/out-of-scope terms, and creation/maintenance of a mapping set with **>9,400 Rhea and EC terms mapped to GO terms**. (g2026thegeneontology pages 5-7)

### 5.2 UniProt–Rhea scale-up and FAIR chemical querying
UniProt’s Rhea-based reaction annotation at tens of millions of sequences, along with ChEBI-driven ligand FAIRification efforts, indicates a shift toward chemically grounded, computable enzyme annotations for large-scale protein–ligand and function analyses. (bateman2023uniprottheuniversal pages 2-3)

### 5.3 Reaction-to-enzyme platforms for engineering (REME, 2024)
REME (NAR 2024) exemplifies real-world implementation: it uses atom-to-atom mapping, atom type change identification, and reaction similarity (RDKit/RXNFP/DRFP), and expands candidate enzymes via EC numbers and sequence homology beyond BRENDA-derived sets; it then evaluates candidates with deep-learning predictors of substrate specificity and kinetics/conditions. (shi2024remeanintegrated pages 1-2, shi2024remeanintegrated pages 2-4)

## 6) Automating reaction and enzyme classification

### 6.1 Sequence → EC prediction (transformers and protein language models)
DeepECtransformer (Nature Communications, 2023) uses transformer layers to predict EC numbers and evaluated on a curated test set of **2,013,612 enzyme sequences** spanning **5,360 EC numbers**; it predicted EC numbers for **464 previously unannotated *E. coli* genes** and experimentally validated three predicted activities. (kim2023functionalannotationof pages 1-2)

This work also highlights a major practical constraint: strong class imbalance affects performance, with F1 correlating with training-set size (Spearman 0.6872, p<0.001). (kim2023functionalannotationof pages 1-2)

### 6.2 Reaction SMILES → EC prediction (reaction-first automation)
BEC-Pred (Journal of Cheminformatics, 2024) predicts EC numbers directly from substrate/product SMILES using a BERT-based approach and reports **91.6% accuracy**, outperforming comparator sequence- and graph-based methods by **5.5%** and improving F1 by ~6%. (qian2024ageneralmodel pages 1-2)

The same work emphasizes why automated reaction classification is still difficult: many reactions lack EC numbers and datasets are skewed, especially for secondary metabolism, limiting deep learning for underrepresented classes. (qian2024ageneralmodel pages 2-3)

### 6.3 Benchmarks and reframing as retrieval
CARE (2024) and ReactZyme (2024) reflect a broader trend toward benchmarking and framing enzyme annotation as classification and retrieval, using data derived from Swiss-Prot and Rhea and noting persistent gaps in functional annotation at genome scale. (yang2024careabenchmark pages 1-2, hua2024reactzymeabenchmark pages 1-2)

## 7) Limitations and open gaps

* **Coverage gaps (unassigned EC numbers)**: many reactions remain outside EC due to evidentiary requirements and incomplete characterization, limiting both curated knowledge and ML training sets. (qian2024ageneralmodel pages 2-3, mcdonald2014fifty‐fiveyearsof pages 7-9)
* **Granularity mismatch**: EC is reaction-centric and may not capture mechanistic differences, organism-specific context, or multifunctionality cleanly; GO may be broader and sometimes insufficiently specific for reaction-level needs. (tipton2000historyofthe pages 3-5, hua2024reactzymeabenchmark pages 1-2)
* **Identifier harmonization** remains a practical barrier in integrative workflows; resources like Rhea’s mapping service explicitly target this need by accepting EC and GO MF identifiers and mapping them across major pathway databases. (bansal2022rheathereaction pages 5-6)

## References (URLs and publication dates)

Key sources used prominently in this report:

* Tipton KF, Boyce S. *History of the enzyme nomenclature system*. **Bioinformatics**. 2000-01. https://doi.org/10.1093/bioinformatics/16.1.34 (tipton2000historyofthe pages 1-2)
* McDonald AG, Tipton KF. *Fifty-five years of enzyme classification: advances and difficulties*. **FEBS J**. 2014-01. https://doi.org/10.1111/febs.12530 (mcdonald2014fifty‐fiveyearsof pages 2-4)
* Bateman A et al. *UniProt: the Universal Protein Knowledgebase in 2023*. **Nucleic Acids Res**. 2023-11. https://doi.org/10.1093/nar/gkac1052 (bateman2023uniprottheuniversal pages 2-3)
* Bansal P et al. *Rhea, the reaction knowledgebase in 2022*. **Nucleic Acids Res**. 2022-11. https://doi.org/10.1093/nar/gkab1016 (bansal2022rheathereaction pages 1-2)
* Milacic M et al. *The Reactome Pathway Knowledgebase 2024*. **Nucleic Acids Res**. 2024-11. https://doi.org/10.1093/nar/gkad1025 (milacic2024thereactomepathway pages 5-6)
* Karp PD et al. *The EcoCyc Database (2023)*. **EcoSal Plus**. 2023-12. https://doi.org/10.1128/ecosalplus.esp-0002-2023 (karp2023theecocycdatabase pages 1-2)
* Kim GB et al. *Functional annotation of enzyme-encoding genes using deep learning with transformer layers*. **Nat Commun**. 2023-11. https://doi.org/10.1038/s41467-023-43216-z (kim2023functionalannotationof pages 1-2)
* Qian W et al. *A general model for predicting enzyme functions based on enzymatic reactions*. **J Cheminform**. 2024-03. https://doi.org/10.1186/s13321-024-00827-y (qian2024ageneralmodel pages 1-2)
* Shi Z et al. *REME: an integrated platform for reaction enzyme mining and evaluation*. **Nucleic Acids Res**. 2024-05. https://doi.org/10.1093/nar/gkae405 (shi2024remeanintegrated pages 1-2)

Note on GO 2023–2024 sourcing: a GO “knowledgebase 2023” NAR update could not be retrieved via the available search results in this run; the most direct evidence on EC/Rhea alignment came from the GO NAR 2026 update, which explicitly describes alignment work that references the 2023 update. (g2026thegeneontology pages 5-7)

References

1. (tipton2000historyofthe pages 1-2): Keith Tipton and Sinéad Boyce. History of the enzyme nomenclature system. Bioinformatics, 16 1:34-40, Jan 2000. URL: https://doi.org/10.1093/bioinformatics/16.1.34, doi:10.1093/bioinformatics/16.1.34. This article has 209 citations and is from a highest quality peer-reviewed journal.

2. (tipton2000historyofthe pages 2-3): Keith Tipton and Sinéad Boyce. History of the enzyme nomenclature system. Bioinformatics, 16 1:34-40, Jan 2000. URL: https://doi.org/10.1093/bioinformatics/16.1.34, doi:10.1093/bioinformatics/16.1.34. This article has 209 citations and is from a highest quality peer-reviewed journal.

3. (tipton2000historyofthe pages 3-5): Keith Tipton and Sinéad Boyce. History of the enzyme nomenclature system. Bioinformatics, 16 1:34-40, Jan 2000. URL: https://doi.org/10.1093/bioinformatics/16.1.34, doi:10.1093/bioinformatics/16.1.34. This article has 209 citations and is from a highest quality peer-reviewed journal.

4. (bansal2022rheathereaction pages 1-2): Parit Bansal, Anne Morgat, Kristian B Axelsen, Venkatesh Muthukrishnan, Elisabeth Coudert, Lucila Aimo, Nevila Hyka-Nouspikel, Elisabeth Gasteiger, Arnaud Kerhornou, Teresa Batista Neto, Monica Pozzato, Marie-Claude Blatter, Alex Ignatchenko, Nicole Redaschi, and Alan Bridge. Rhea, the reaction knowledgebase in 2022. Nucleic Acids Research, 50:D693-D700, Nov 2022. URL: https://doi.org/10.1093/nar/gkab1016, doi:10.1093/nar/gkab1016. This article has 275 citations and is from a highest quality peer-reviewed journal.

5. (bateman2023uniprottheuniversal pages 2-3): A. Bateman, M. Martin, S. Orchard, M. Magrane, Shadab Ahmad, E. Alpi, E. Bowler-Barnett, R. Britto, Hema Bye-A-Jee, Austra Cukura, Paul Denny, Tunca Dogan, Thankgod Ebenezer, Jun Fan, Penelope Garmiri, Leonardo Jose da Costa Gonzales, E. Hatton-Ellis, Abdulrahman Hussein, A. Ignatchenko, Giuseppe Insana, Rizwan Ishtiaq, Vishal Joshi, D. Jyothi, Swaathi Kandasaamy, A. Lock, Aurélien Luciani, Marija Lugarić, Jie Luo, Yvonne Lussi, Alistair MacDougall, F. Madeira, Mahdi Mahmoudy, Alok Mishra, Katie Moulang, A. Nightingale, Sangya Pundir, G. Qi, Shriya Raj, P. Raposo, Daniel L Rice, Rabie Saidi, Rafael Santos, Elena Speretta, J. Stephenson, Prabhat Totoo, Edward Turner, N. Tyagi, Preethi Vasudev, Kate Warner, Xavier Watkins, Rossana Zaru, H. Zellner, A. Bridge, L. Aimo, Ghislaine Argoud-Puy, A. Auchincloss, K. Axelsen, Parit Bansal, Delphine Baratin, Teresa M Batista Neto, M. Blatter, Jerven T. Bolleman, E. Boutet, L. Breuza, B. Gil, Cristina Casals-Casas, Kamal Chikh Echioukh, E. Coudert, Béatrice A. Cuche, Edouard de Castro, A. Estreicher, M. Famiglietti, M. Feuermann, E. Gasteiger, P. Gaudet, S. Gehant, V. Gerritsen, A. Gos, N. Gruaz, C. Hulo, Nevila Hyka-Nouspikel, F. Jungo, A. Kerhornou, Philippe le Mercier, D. Lieberherr, P. Masson, A. Morgat, Venkatesh Muthukrishnan, S. Paesano, I. Pedruzzi, S. Pilbout, L. Pourcel, S. Poux, Monica Pozzato, Manuela Pruess, Nicole Redaschi, C. Rivoire, Christian J. A. Sigrist, K. Sonesson, S. Sundaram, Cathy H. Wu, C. Arighi, L. Arminski, Chuming Chen, Yongxing Chen, Hongzhan Huang, K. Laiho, P. McGarvey, D. Natale, K. Ross, C. R. Vinayaka, Qinghua Wang, Yuqi Wang, and Jian Zhang. Uniprot: the universal protein knowledgebase in 2023. Nucleic Acids Research, 51:D523-D531, Nov 2023. URL: https://doi.org/10.1093/nar/gkac1052, doi:10.1093/nar/gkac1052. This article has 6107 citations and is from a highest quality peer-reviewed journal.

6. (milacic2024thereactomepathway pages 5-6): Marija Milacic, Deidre Beavers, Patrick Conley, Chuqiao Gong, Marc Gillespie, Johannes Griss, Robin Haw, Bijay Jassal, Lisa Matthews, Bruce May, Robert Petryszak, Eliot Ragueneau, Karen Rothfels, Cristoffer Sevilla, Veronica Shamovsky, Ralf Stephan, Krishna Tiwari, Thawfeek Varusai, Joel Weiser, Adam Wright, Guanming Wu, Lincoln Stein, Henning Hermjakob, and Peter D’Eustachio. The reactome pathway knowledgebase 2024. Nucleic Acids Research, 52:D672-D678, Nov 2024. URL: https://doi.org/10.1093/nar/gkad1025, doi:10.1093/nar/gkad1025. This article has 1279 citations and is from a highest quality peer-reviewed journal.

7. (kim2023functionalannotationof pages 1-2): Gi Bae Kim, Ji Yeon Kim, Jong An Lee, Charles J. Norsigian, Bernhard O. Palsson, and Sang Yup Lee. Functional annotation of enzyme-encoding genes using deep learning with transformer layers. Nature Communications, Nov 2023. URL: https://doi.org/10.1038/s41467-023-43216-z, doi:10.1038/s41467-023-43216-z. This article has 109 citations and is from a highest quality peer-reviewed journal.

8. (qian2024ageneralmodel pages 1-2): Wenjia Qian, Xiaorui Wang, Yu Kang, Peichen Pan, Tingjun Hou, and Chang-Yu Hsieh. A general model for predicting enzyme functions based on enzymatic reactions. Journal of Cheminformatics, Mar 2024. URL: https://doi.org/10.1186/s13321-024-00827-y, doi:10.1186/s13321-024-00827-y. This article has 19 citations and is from a peer-reviewed journal.

9. (shi2024remeanintegrated pages 1-2): Zhenkun Shi, Dehang Wang, Yang Li, Rui Deng, Jiawei Lin, Cui Liu, Haoran Li, Ruoyu Wang, Muqiang Zhao, Zhitao Mao, Qianqian Yuan, Xiaoping Liao, and Hongwu Ma. Reme: an integrated platform for reaction enzyme mining and evaluation. Nucleic Acids Research, 52:W299-W305, May 2024. URL: https://doi.org/10.1093/nar/gkae405, doi:10.1093/nar/gkae405. This article has 12 citations and is from a highest quality peer-reviewed journal.

10. (mcdonald2014fifty‐fiveyearsof pages 2-4): Andrew G. McDonald and Keith F. Tipton. Fifty‐five years of enzyme classification: advances and difficulties. The FEBS Journal, 281:583-592, Jan 2014. URL: https://doi.org/10.1111/febs.12530, doi:10.1111/febs.12530. This article has 170 citations.

11. (tipton2000historyofthe pages 6-7): Keith Tipton and Sinéad Boyce. History of the enzyme nomenclature system. Bioinformatics, 16 1:34-40, Jan 2000. URL: https://doi.org/10.1093/bioinformatics/16.1.34, doi:10.1093/bioinformatics/16.1.34. This article has 209 citations and is from a highest quality peer-reviewed journal.

12. (karp2023theecocycdatabase pages 1-2): Peter D. Karp, Suzanne Paley, Ron Caspi, Anamika Kothari, Markus Krummenacker, Peter E. Midford, Lisa R. Moore, Pallavi Subhraveti, Socorro Gama-Castro, Victor H. Tierrafria, Paloma Lara, Luis Muñiz-Rascado, César Bonavides-Martinez, Alberto Santos-Zavaleta, Amanda Mackie, Gwanggyu Sun, Travis A. Ahn-Horst, Heejo Choi, Markus W. Covert, Julio Collado-Vides, and Ian Paulsen. The ecocyc database (2023). EcoSal Plus, Dec 2023. URL: https://doi.org/10.1128/ecosalplus.esp-0002-2023, doi:10.1128/ecosalplus.esp-0002-2023. This article has 94 citations.

13. (hua2024reactzymeabenchmark pages 1-2): Chenqing Hua, Bozitao Zhong, Sitao Luan, Liang Hong, Guy Wolf, Doina Precup, and Shuangjia Zheng. Reactzyme: a benchmark for enzyme-reaction prediction. ArXiv, Aug 2024. URL: https://doi.org/10.48550/arxiv.2408.13659, doi:10.48550/arxiv.2408.13659. This article has 31 citations.

14. (xie2024deeplearningin pages 2-3): Xueying Xie, Lin Gui, Baixue Qiao, Guohua Wang, Shan Huang, Yuming Zhao, and Shanwen Sun. Deep learning in template-free de novo biosynthetic pathway design of natural products. Briefings in Bioinformatics, Sep 2024. URL: https://doi.org/10.1093/bib/bbae495, doi:10.1093/bib/bbae495. This article has 12 citations and is from a domain leading peer-reviewed journal.

15. (bansal2022rheathereaction pages 1-1): Parit Bansal, Anne Morgat, Kristian B Axelsen, Venkatesh Muthukrishnan, Elisabeth Coudert, Lucila Aimo, Nevila Hyka-Nouspikel, Elisabeth Gasteiger, Arnaud Kerhornou, Teresa Batista Neto, Monica Pozzato, Marie-Claude Blatter, Alex Ignatchenko, Nicole Redaschi, and Alan Bridge. Rhea, the reaction knowledgebase in 2022. Nucleic Acids Research, 50:D693-D700, Nov 2022. URL: https://doi.org/10.1093/nar/gkab1016, doi:10.1093/nar/gkab1016. This article has 275 citations and is from a highest quality peer-reviewed journal.

16. (bansal2022rheathereaction pages 5-6): Parit Bansal, Anne Morgat, Kristian B Axelsen, Venkatesh Muthukrishnan, Elisabeth Coudert, Lucila Aimo, Nevila Hyka-Nouspikel, Elisabeth Gasteiger, Arnaud Kerhornou, Teresa Batista Neto, Monica Pozzato, Marie-Claude Blatter, Alex Ignatchenko, Nicole Redaschi, and Alan Bridge. Rhea, the reaction knowledgebase in 2022. Nucleic Acids Research, 50:D693-D700, Nov 2022. URL: https://doi.org/10.1093/nar/gkab1016, doi:10.1093/nar/gkab1016. This article has 275 citations and is from a highest quality peer-reviewed journal.

17. (ragueneau2026thereactomeknowledgebase pages 1-2): Eliot Ragueneau, Chuqiao Gong, Pierre Sinquin, Cristoffer Sevilla, Deidre Beavers, Alexander Grentner, Johannes Griss, Gregory F J Hogue, Nancy T Li, Lisa Matthews, Bruce May, Marija Milacic, Helia Mohammadi, Robert Petryszak, Karen Rothfels, Veronica Shamovsky, Ralf Stephan, Krishna Tiwari, Joel Weiser, Adam Wright, Marc Gillespie, Guanming Wu, Lincoln Stein, Henning Hermjakob, and Peter D’Eustachio. The reactome knowledgebase 2026. Nucleic Acids Research, 54:D673-D681, Nov 2026. URL: https://doi.org/10.1093/nar/gkaf1223, doi:10.1093/nar/gkaf1223. This article has 6 citations and is from a highest quality peer-reviewed journal.

18. (presern2023enzymedatabasesin pages 9-11): Uroš Prešern and Marko Goličnik. Enzyme databases in the era of omics and artificial intelligence. International Journal of Molecular Sciences, 24:16918, Nov 2023. URL: https://doi.org/10.3390/ijms242316918, doi:10.3390/ijms242316918. This article has 22 citations.

19. (wei2025findingthedark pages 7-8): Galen Wei, Xinchun Ran, Runeem Ai-Abssi, and Zhongyue Yang. Finding the dark matter: large language model‐based enzyme kinetic data extractor and its validation. Protein Science : A Publication of the Protein Society, Aug 2025. URL: https://doi.org/10.1002/pro.70251, doi:10.1002/pro.70251. This article has 2 citations.

20. (presern2023enzymedatabasesin pages 8-9): Uroš Prešern and Marko Goličnik. Enzyme databases in the era of omics and artificial intelligence. International Journal of Molecular Sciences, 24:16918, Nov 2023. URL: https://doi.org/10.3390/ijms242316918, doi:10.3390/ijms242316918. This article has 22 citations.

21. (mcdonald2014fifty‐fiveyearsof pages 7-9): Andrew G. McDonald and Keith F. Tipton. Fifty‐five years of enzyme classification: advances and difficulties. The FEBS Journal, 281:583-592, Jan 2014. URL: https://doi.org/10.1111/febs.12530, doi:10.1111/febs.12530. This article has 170 citations.

22. (mcdonald2014fifty‐fiveyearsof pages 6-7): Andrew G. McDonald and Keith F. Tipton. Fifty‐five years of enzyme classification: advances and difficulties. The FEBS Journal, 281:583-592, Jan 2014. URL: https://doi.org/10.1111/febs.12530, doi:10.1111/febs.12530. This article has 170 citations.

23. (qian2024ageneralmodel pages 2-3): Wenjia Qian, Xiaorui Wang, Yu Kang, Peichen Pan, Tingjun Hou, and Chang-Yu Hsieh. A general model for predicting enzyme functions based on enzymatic reactions. Journal of Cheminformatics, Mar 2024. URL: https://doi.org/10.1186/s13321-024-00827-y, doi:10.1186/s13321-024-00827-y. This article has 19 citations and is from a peer-reviewed journal.

24. (garapati2026comprehensiveannotationof pages 1-2): Phani V. Garapati, Rossana Zaru, H. Attrill, Gilberto dos Santos, J. Goodman, Jim Thurmond, and Steven J. Marygold. Comprehensive annotation of the enzymes of drosophila melanogaster. G3: Genes | Genomes | Genetics, Dec 2026. URL: https://doi.org/10.1093/g3journal/jkaf294, doi:10.1093/g3journal/jkaf294. This article has 0 citations.

25. (g2026thegeneontology pages 5-7): Suzi A James P Seth J Michael Dustin Marc Pascale Nomi Aleksander Balhoff Carbon Cherry Ebert Feuermann G, Suzi A Aleksander, J. Balhoff, S. Carbon, J. M. Cherry, Dustin Ebert, M. Feuermann, Pascale Gaudet, Nomi L. Harris, David P. Hill, Patrick Kalita, Raymond Lee, H. Mi, Sierra A T Moxon, C. Mungall, A. Muruganujan, Tremayne Mushayahama, Paul W. Sternberg, P. D. Thomas, K. V. Van Auken, E. Wong, V. Wood, J. Ramsey, D. Siegele, Rex L. Chisholm, Robert H. T. Dodson, Petra Fey, M. C. Aspromonte, María Victoria Nugnes, Ximena Aixa Castro Naser, Silvio C. E. Tosatto, Michelle Giglio, Suvarna Nadendla, Giulia Antonazzo, H. Attrill, Nicholas H. Brown, Gilberto dos Santos, Steven J. Marygold, Katja Röper, V. Strelets, Christopher J. Tabone, Jim Thurmond, Pinglei Zhou, Rossana Zaru, R. Lovering, Colin Logie, Daqing Chen, A. Naba, K. Christie, L. Corbani, L. Ni, D. Sitnikov, Cynthia L. Smith, James Seager, L. Cooper, J. Elser, Pankaj Jaiswal, Parul Gupta, Sushma Naithani, Pascal Carme, Kim Rutherford, J. D. De Pons, M. Dwinell, G. Hayman, M. Kaldunski, A. Kwitek, S. Laulederkind, M. Tutaj, M. Vedi, Shur-Jen Wang, Peter D’Eustachio, L. Aimo, K. Axelsen, Alan J. Bridge, Nevila Hyka-Nouspikel, A. Morgat, Gene Goldbold, S. Engel, S. Miyasato, R. Nash, G. Sherlock, S. Weng, E. Bakker, T. Berardini, L. Reiser, A. Auchincloss, Ghislaine Argoud-Puy, M. Blatter, E. Boutet, L. Breuza, Cristina Casals-Casas, E. Coudert, A. Estreicher, M. Famiglietti, A. Gos, N. Gruaz-Gumowski, C. Hulo, F. Jungo, P. L. Mercier, D. Lieberherr, P. Masson, I. Pedruzzi, L. Pourcel, S. Poux, C. Rivoire, S. Sundaram, Alex Bateman, A. Adesina, E. Bowler-Barnett, D. Carpentier, Paul Denny, Alexandr Ignatchenko, Rizwan Ishtiaq, Antonia Lock, Yvonne Lussi, M. Magrane, María J. Martín, Sandra Orchard, Pedro Raposo, Elena Speretta, Nidhi Tyagi, Nadya Urakova, Kate Warner, Conny Wing-Heng Yu, J. Chan, Stavros Diamantakis, Mark Quinton-Tulloch, Daniela Raciti, Malcolm Fisher, Christina James-Zorn, V. Ponferrada, Aaron Zorn, D. Howe, S. Ramachandran, Leyla Ruzicka, and M. Westerfield. The gene ontology knowledgebase in 2026. Nucleic Acids Research, 54:D1779-D1792, Dec 2026. URL: https://doi.org/10.1093/nar/gkaf1292, doi:10.1093/nar/gkaf1292. This article has 1026 citations and is from a highest quality peer-reviewed journal.

26. (mcdonald2014fifty‐fiveyearsof pages 1-2): Andrew G. McDonald and Keith F. Tipton. Fifty‐five years of enzyme classification: advances and difficulties. The FEBS Journal, 281:583-592, Jan 2014. URL: https://doi.org/10.1111/febs.12530, doi:10.1111/febs.12530. This article has 170 citations.

27. (bateman2023uniprottheuniversal pages 4-5): A. Bateman, M. Martin, S. Orchard, M. Magrane, Shadab Ahmad, E. Alpi, E. Bowler-Barnett, R. Britto, Hema Bye-A-Jee, Austra Cukura, Paul Denny, Tunca Dogan, Thankgod Ebenezer, Jun Fan, Penelope Garmiri, Leonardo Jose da Costa Gonzales, E. Hatton-Ellis, Abdulrahman Hussein, A. Ignatchenko, Giuseppe Insana, Rizwan Ishtiaq, Vishal Joshi, D. Jyothi, Swaathi Kandasaamy, A. Lock, Aurélien Luciani, Marija Lugarić, Jie Luo, Yvonne Lussi, Alistair MacDougall, F. Madeira, Mahdi Mahmoudy, Alok Mishra, Katie Moulang, A. Nightingale, Sangya Pundir, G. Qi, Shriya Raj, P. Raposo, Daniel L Rice, Rabie Saidi, Rafael Santos, Elena Speretta, J. Stephenson, Prabhat Totoo, Edward Turner, N. Tyagi, Preethi Vasudev, Kate Warner, Xavier Watkins, Rossana Zaru, H. Zellner, A. Bridge, L. Aimo, Ghislaine Argoud-Puy, A. Auchincloss, K. Axelsen, Parit Bansal, Delphine Baratin, Teresa M Batista Neto, M. Blatter, Jerven T. Bolleman, E. Boutet, L. Breuza, B. Gil, Cristina Casals-Casas, Kamal Chikh Echioukh, E. Coudert, Béatrice A. Cuche, Edouard de Castro, A. Estreicher, M. Famiglietti, M. Feuermann, E. Gasteiger, P. Gaudet, S. Gehant, V. Gerritsen, A. Gos, N. Gruaz, C. Hulo, Nevila Hyka-Nouspikel, F. Jungo, A. Kerhornou, Philippe le Mercier, D. Lieberherr, P. Masson, A. Morgat, Venkatesh Muthukrishnan, S. Paesano, I. Pedruzzi, S. Pilbout, L. Pourcel, S. Poux, Monica Pozzato, Manuela Pruess, Nicole Redaschi, C. Rivoire, Christian J. A. Sigrist, K. Sonesson, S. Sundaram, Cathy H. Wu, C. Arighi, L. Arminski, Chuming Chen, Yongxing Chen, Hongzhan Huang, K. Laiho, P. McGarvey, D. Natale, K. Ross, C. R. Vinayaka, Qinghua Wang, Yuqi Wang, and Jian Zhang. Uniprot: the universal protein knowledgebase in 2023. Nucleic Acids Research, 51:D523-D531, Nov 2023. URL: https://doi.org/10.1093/nar/gkac1052, doi:10.1093/nar/gkac1052. This article has 6107 citations and is from a highest quality peer-reviewed journal.

28. (bansal2022rheathereaction media 1bea24f1): Parit Bansal, Anne Morgat, Kristian B Axelsen, Venkatesh Muthukrishnan, Elisabeth Coudert, Lucila Aimo, Nevila Hyka-Nouspikel, Elisabeth Gasteiger, Arnaud Kerhornou, Teresa Batista Neto, Monica Pozzato, Marie-Claude Blatter, Alex Ignatchenko, Nicole Redaschi, and Alan Bridge. Rhea, the reaction knowledgebase in 2022. Nucleic Acids Research, 50:D693-D700, Nov 2022. URL: https://doi.org/10.1093/nar/gkab1016, doi:10.1093/nar/gkab1016. This article has 275 citations and is from a highest quality peer-reviewed journal.

29. (milacic2024thereactomepathway pages 1-2): Marija Milacic, Deidre Beavers, Patrick Conley, Chuqiao Gong, Marc Gillespie, Johannes Griss, Robin Haw, Bijay Jassal, Lisa Matthews, Bruce May, Robert Petryszak, Eliot Ragueneau, Karen Rothfels, Cristoffer Sevilla, Veronica Shamovsky, Ralf Stephan, Krishna Tiwari, Thawfeek Varusai, Joel Weiser, Adam Wright, Guanming Wu, Lincoln Stein, Henning Hermjakob, and Peter D’Eustachio. The reactome pathway knowledgebase 2024. Nucleic Acids Research, 52:D672-D678, Nov 2024. URL: https://doi.org/10.1093/nar/gkad1025, doi:10.1093/nar/gkad1025. This article has 1279 citations and is from a highest quality peer-reviewed journal.

30. (shi2024remeanintegrated pages 2-4): Zhenkun Shi, Dehang Wang, Yang Li, Rui Deng, Jiawei Lin, Cui Liu, Haoran Li, Ruoyu Wang, Muqiang Zhao, Zhitao Mao, Qianqian Yuan, Xiaoping Liao, and Hongwu Ma. Reme: an integrated platform for reaction enzyme mining and evaluation. Nucleic Acids Research, 52:W299-W305, May 2024. URL: https://doi.org/10.1093/nar/gkae405, doi:10.1093/nar/gkae405. This article has 12 citations and is from a highest quality peer-reviewed journal.

31. (yang2024careabenchmark pages 1-2): Jason Yang, Ariane Mora, Shengchao Liu, Bruce J. Wittmann, Anima Anandkumar, Frances H. Arnold, and Yisong Yue. Care: a benchmark suite for the classification and retrieval of enzymes. Preprint, Jan 2024. URL: https://doi.org/10.48550/arxiv.2406.15669, doi:10.48550/arxiv.2406.15669. This article has 26 citations.

## Citations

1. tipton2000historyofthe pages 1-2
2. bansal2022rheathereaction pages 1-2
3. karp2023theecocycdatabase pages 1-2
4. kim2023functionalannotationof pages 1-2
5. qian2024ageneralmodel pages 1-2
6. hua2024reactzymeabenchmark pages 1-2
7. presern2023enzymedatabasesin pages 8-9
8. garapati2026comprehensiveannotationof pages 1-2
9. g2026thegeneontology pages 5-7
10. bansal2022rheathereaction pages 5-6
11. bateman2023uniprottheuniversal pages 2-3
12. bateman2023uniprottheuniversal pages 4-5
13. milacic2024thereactomepathway pages 1-2
14. milacic2024thereactomepathway pages 5-6
15. qian2024ageneralmodel pages 2-3
16. shi2024remeanintegrated pages 1-2
17. tipton2000historyofthe pages 2-3
18. tipton2000historyofthe pages 3-5
19. tipton2000historyofthe pages 6-7
20. xie2024deeplearningin pages 2-3
21. bansal2022rheathereaction pages 1-1
22. ragueneau2026thereactomeknowledgebase pages 1-2
23. presern2023enzymedatabasesin pages 9-11
24. wei2025findingthedark pages 7-8
25. shi2024remeanintegrated pages 2-4
26. yang2024careabenchmark pages 1-2
27. https://doi.org/10.1093/bioinformatics/16.1.34
28. https://doi.org/10.1111/febs.12530
29. https://doi.org/10.1093/nar/gkac1052
30. https://doi.org/10.1093/nar/gkab1016
31. https://doi.org/10.1093/nar/gkad1025
32. https://doi.org/10.1128/ecosalplus.esp-0002-2023
33. https://doi.org/10.1038/s41467-023-43216-z
34. https://doi.org/10.1186/s13321-024-00827-y
35. https://doi.org/10.1093/nar/gkae405
36. https://doi.org/10.1093/bioinformatics/16.1.34,
37. https://doi.org/10.1093/nar/gkab1016,
38. https://doi.org/10.1093/nar/gkac1052,
39. https://doi.org/10.1093/nar/gkad1025,
40. https://doi.org/10.1038/s41467-023-43216-z,
41. https://doi.org/10.1186/s13321-024-00827-y,
42. https://doi.org/10.1093/nar/gkae405,
43. https://doi.org/10.1111/febs.12530,
44. https://doi.org/10.1128/ecosalplus.esp-0002-2023,
45. https://doi.org/10.48550/arxiv.2408.13659,
46. https://doi.org/10.1093/bib/bbae495,
47. https://doi.org/10.1093/nar/gkaf1223,
48. https://doi.org/10.3390/ijms242316918,
49. https://doi.org/10.1002/pro.70251,
50. https://doi.org/10.1093/g3journal/jkaf294,
51. https://doi.org/10.1093/nar/gkaf1292,
52. https://doi.org/10.48550/arxiv.2406.15669,
