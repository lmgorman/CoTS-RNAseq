# CoTS-RNAseq
Bioinformatics analysis of RNA-seq data collected from corals (_Acropora hyacinthus_ and _Porites_ sp.) being eaten by CoTS in Mo'orea, French Polynesia during the 2022-26 CoTS outbreak.

Four treatments:

_Acropora hyacinthus_ eaten: 497R, 568R, 410R, 549R, 414R, 321R A, 581R, 370R

_Acropora hyacinthus_ control: 419R, 331R A, 336R, 512R, 380R, 571R, 586R, 468R

_Porites_ sp. eaten: 225R, 235R, 43R, 211R, 218R, 34R, 227R, 16R

_Porites_ sp. control: 61R, 86R, 236R, 244R, 82R, 71R, 253R, 76R

# Māuruuru roa
As guests, we would like to give thanks to the land, water and resources of Polynesia, in particular Mo’orea, which allowed us to carry out this study. We would also like to thank the traditional land owners of the past, present and future - Māuruuru roa.

# **Sample collection methodology**

Porites sp. were collected from the lagoon at Ta'ahiamanu (-17.495210°, -149.851340°), Mo’orea, French Polynesia, across a 6-day period (30 January - 4 February 2024), between 9am-12pm. Similarly, Acropora hyacinthus specimens were collected from the outer reef (due to no individuals being eaten in the lagoon) adjacent to the Plage publique de Temae (-17.501364, -149.756490), Mo’orea, French Polynesia across a 5 week period (12 March - 24 April 2024), between 9am-12pm. 

A. hyacinthus and Porites sp. that had been preyed on by CoTS (“eaten”) and ones that had not been preyed on (“control”) were chosen (Figure 1). Control colonies were chosen based on colonies that were not being actively eaten by CoTS and also showed no observable historic CoTS feeding scars (Figure 1A, C). Predation was determined if a CoTS had its stomach everted over the colony, which was checked by inverting the CoTS (Figure 2). For the colonies being eaten, coral samples were excised from where the CoTS had been in contact with the coral colony. Coral biopsies were taken with stainless steel coral bone cutters, cleaned by wiping in seawater and rinsing with ethanol. The biopsies were placed in 1.5ml screwcap tubes containing 0.75ml of DNA/RNA shield (Zymo, Cat number: R1100-50). Colonies (eaten and control) were sampled at the same times (9am-12pm) successively, alternating between an eaten and control sample, to avoid any changes in gene expression due to diurnal fluctuations. Each coral colony represented one biological replicate and altogether, 15 biological replicates of eaten and control colonies of each species (Porites sp. and A. hyacinthus) were taken (30 colonies of each genera in total), the number of samples sent for sequencing is detailed below. 




<img width="1217" height="886" alt="image" src="https://github.com/user-attachments/assets/2432de55-0f61-4a06-b590-2b78fa5c042d" />



**Figure 1.** _Porites_ sp. (**A**, **B**) and _Acropora hyacinthus_ (**C**, **D**) individuals during the 2022/24 CoTS outbreak in Mo'orea, French Polynesia. Control colonies for each species were chosen based on no observable CoTS feeding scars (**A**, **C**). Colonies being fed upon by CoTS (**B**, **D**) were chosen that were actively consumed by CoTS during sampling and harboured both an observable feeding scar (skeleton devoid of coral tissue (red arrow)) in addition to harbouring areas with remaining coral tissue (pink arrow) for sampling. 


# RNA extraction
Prior to extraction, Acropora hyacinthus tissues were already dissociated from the skeleton and were thus only mixed by vortexing. For Porites sp. tissue was still remaining on the coral skeleton and therefore, the sample was first pulverised with a metal blade to access fully and then mixed by vortexing. Aliquots from each sample (300-1500 µL) were then transferred to a new tube and RNA was extracted using the Zymo RNA/ DNA quick prep mini kit (Cat number: D7001) following manufacturers instructions. Samples of total RNA with a 260/280 ratio between 1.8-2.2 measured on an N60 nanospectrophotometer (IMPLEN, Munich, Germany) and containing two RNA peaks (18S and 28S) on the Agilent Bioanalyser were sent for sequencing. Altogether, eight biological replicates for each treatment (n = 2 treatments - eaten versus control) for each species were sent for sequencing.

RNA sequencing and library construction
RNA sequencing and library construction were carried out at MGX-Montpellier GenomiX (1). Ribosomal RNA was first removed by Illumina Ribo-Zero (Cat number: 20040526) to purify each sample and leave mRNA. Libraries were prepared via Polyadenylated RNA selection with polyT capture beads using Illumina's Truseq stranded mRNA kit (Cat number: 20020594) according to manufacturer's instructions. Libraries were validated by DNA quantification for fragment size and concentration using a Fragment Bioanalyzer (Agilent) and qPCR (Roche), respectively. The resulting RNA was sequenced 2 x 150pb on one lane of a Novaseq S4 (Illumina).

# Bioinformatics
Raw RNA sequence reads were run through fastaQC (v.0.12.1, Babraham Institute, 2023) to check read quality metrics. One sample of pair-end reads in the Porites sp. eaten treatment had an abnormally high number of reads in the sequence library (> 610 million) and seqTK software v.1.4 (Li, 2023) was used to randomly subset to 90 million reads to mirror the average number of reads in the other samples. Pair-end reads from all samples were then trimmed using fastp v.0.23.2. (Chen, 2023) with the following parameters: a qualified base quality value of 20; an unqualified base limit of 10%; window size of 5 and a quality threshold of 20 for the window. Trimmed reads were then assessed for their quality using fastaQC and multiQC (v.1.12, Ewels et al., 2016). Acropora hyacinthus samples were then mapped to the A. hyacinthus genome (Lopez-Nadam et al., 2023) using the STAR software (v.2.7.11.b, Dobin et al., 2013). Porites cannot be reliably identified to species visually (Forsman et al., 2009), therefore Porites species samples were mapped to both an Australian P. lutea genome (Robbins et al., 2019) and a French Polynesia P. evermanni genome (Planes et al., 2019). All samples from Porites species showed higher mapping to the P. evermanni genome. Thus, the P. evermanni genome was chosen for subsequent analysis. RNA Seq reads from A. hyacinthus and Porites sp. were then aligned to the corresponding genome using STAR (v.2.7.11.b). 

Samples were checked for outliers in terms of the number of total mapped reads, as well as for global gene expression using Principle Coordinate Analysis (PCA) plots. Samples with a low percent of total mapped reads (<5%) and showing increased dispersion on PCA plots were removed from the datasets, corresponding to one eaten A. hyacinthus (sample ID: 497R) and four eaten Porites sp. (sample IDs: 225R, 34R, 43R and 16R) replicates being removed from further analysis (Table S1). For A. hyacinthus, seven and eight biological replicates were used for eaten and control colonies, respectively. Whereas for Porites sp. four and eight biological replicates were used for eaten and control colonies, respectively. Expression tables and gene count matrices were formed using StringTie v3.0.0 (Shumate et al., 2022), resulting in 27,110 genes for A. hyacinthus and 40,389 for Porites sp. Genes with 0 counts were then removed from the dataset and, samples were filtered using pOverA to retain genes present in 44% of samples in A. hyacinthus (proportion of samples in smallest treatment group) and 34% of samples in Porites sp. with transcript counts >10, resulting in 18,258 total genes for A. hyacinthus and 22,214 for Porites spp.

# Analysis
Genes differentially expressed in “eaten” colonies were identified using differential gene expression analysis. Gene counts were normalized using a variance stabilized transformation in the DESeq2 package in R (Love et al., 2014). A permutational analysis of variance (PERMANOVA) was then used to test for the effects of predation status for all genes that passed filtering thresholds (18,258 and 22,214 genes for A. hyacinthus and Porites sp., respectively) followed by a permutational analysis of dispersion (PERMDISP) to test for variation in dispersion in the vegan package in R (Oksanen et al., 2025). A Wald test in DESeq2 was then used to compute the differentially expressed genes between eaten and control colonies.

Prior to functional enrichment, the general feature format (gff) file of the A. hyacinthus genome was formatted and cleaned using the AGAT toolkit v1.4.1 plugin in python (Dainat, 2022). Functional enrichment of the A. hyacinthus genome was performed as in Conn et al. (2025) using funannotate v.1.8.17 (Palmer & Stajich, 2020). Briefly, the A. hyacinthus genome was soft-masked using the mask function in funannotate v.1.8.17 (Palmer & Stajich, 2020). The soft-masked A. hyacinthus genome was then annotated by firstly assigning gene function and ontology using the eggnog mapper v.2.1.12 (Cantalapiedra et al., 2021) and InterProScan v.5.73-104.0 (Jones et al., 2014). The outputs of these were then inputted into funannotate wrapper to run alongside funannotate. Funannotate then searches the PFAM (Paysan-Lafosse et al., 2025), CAZy (Drula et al., 2022), UniProt (The UniProt Consortium, 2025) and GO (Ashburner et al., 2000; The Gene Ontology Consortium, 2026) databases and assigns all the results to each gene. Due to the format of the P. evermanni genome, an alternative method was used to annotate the functional enrichment. The Blast2GO suite (Götz et al., 2008) in OmicsBox (v3.5; OmicsBox, 2019) was used to functionally annotate the P. evermanni genome. These ontologies were then appended to the differentially expressed genes in each dataset using the topGO package in R (Alexa & Rahnenfuhrer, 2025). A node size of 5 was chosen with an ontology of biological process (“BP”) and an annotation setting of annFUN.gene2GO. A weighted Fisher exact test with P < 0.05 was used to find significant GO terms between eaten and control colonies. 

# Pure venom extraction 
To confirm whether the genes that were expressed during the transcriptomics were in fact venom-related, pure venom was extracted from _Porites_ sp. and _Acropora hyacinthus_ colonies. Pure venom from _Porites_ sp. colonies were taken from colonies at the same location as transcriptomic samples (Ta'ahiamanu, Moorea; -17.495210°, -149.851340°), whereas _A. hyacinthus_ pure venom was taken from colonies at a slightly different location in Moorea (-17.480528°, -149.84791667°). Pure venom was extracted during October 2024 (9-10 months after transcriptomic samples were taken). Small coral nubbins (_n_ = 5) around ~2.5cm in length were taken from _Porites_ sp. and _Acropora hyacinthus_ and placed into 96% EtOH for 30 seconds, allowing nematocysts to successfully discharge venom. Extra care was taken to make sure minimal mechanical pressure was applied to the nubbins prior to placing in ethanol and nubbins were only touched with forceps. This procedure was repeated five times for each species to get a total of 5 biological replicates _per_ species. Samples were then transported at room temperature to CNRS, Montpellier where they underwent proteomic analysis.

# Data Availability
Raw sequences can be found at NCBI under the BioProject accession: PRJNA1230057. 

All coral collections were conducted under permits issued by the French Polynesian authorities (Arrêté n° 3091/MPR/DRM and CITES n° FR1998700269-E), within the framework of a LabEx CORAIL grant “AcantCorVenin” to SCM & LG and Pacific Funds grants “COTS Pacifique” to SCM, LG, MB & HP. Ethical permits were granted by CNRS Animal Experimentation, R-13-CNRS-F1-16 to Yann Lacube and ANZCCART ComPass Animal Welfare Training certificate to SCM. 

# Bioinformatics Analysis happening:
See here - https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/_posts/2025-02-22-COTS-Project-Bioinformatics-Analysis.md 


See here - https://github.com/lmgorman/CoTS-RNAseq/blob/bee3172fc9aed27c553984132a5a2e36353e6cb4/LMG_Lab_Notebook.md

# References
1. MGX-Montpellier GenomiX, Univ. Montpellier, CNRS, INSERM, Montpellier, France
