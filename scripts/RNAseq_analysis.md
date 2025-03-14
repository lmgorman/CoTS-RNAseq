---
Title: "RNAseq Differential Gene Expression Analysis"
Authors: "Lucy M Gorman, Ariana S Huffmyer"
Date: "3/13/2025"
Output: html_document
Editor_options: 
  chunk_output_type: console
---

Based on code from Ariana S Huffmyer:
https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/rna-seq_DEG.Rmd
https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/rna-seq_functional_enrichment_topGO.Rmd

Functional enrichment of _Acropora hyacinthus_ and _Porites_ sp. eaten versus non-eaten by CoTS DEGs.  

# Set up 

Before loading R make sure you click "Run as Administrator"

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE) #Set Strings to character
install.packages("BiocManager", repos = "https://cloud.r-project.org")
```

```{r}
if ("tidyverse" %in% rownames(installed.packages()) == 'FALSE') install.packages('tidyverse') 
if ("genefilter" %in% rownames(installed.packages()) == 'FALSE') install.packages('genefilter') 
if ("DESeq2" %in% rownames(installed.packages()) == 'FALSE') install.packages('DESeq2') 
if ("RColorBrewer" %in% rownames(installed.packages()) == 'FALSE') install.packages('RColorBrewer') 
if ("WGCNA" %in% rownames(installed.packages()) == 'FALSE') install.packages('WGCNA') 
if ("flashClust" %in% rownames(installed.packages()) == 'FALSE') install.packages('flashClust') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("ComplexHeatmap" %in% rownames(installed.packages()) == 'FALSE') install.packages('ComplexHeatmap') 
if ("goseq" %in% rownames(installed.packages()) == 'FALSE') install.packages('goseq') 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("clusterProfiler" %in% rownames(installed.packages()) == 'FALSE') install.packages('clusterProfiler') 
if ("pheatmap" %in% rownames(installed.packages()) == 'FALSE') install.packages('pheatmap') 
if ("magrittr" %in% rownames(installed.packages()) == 'FALSE') install.packages('magrittr') 
if ("rtracklayer" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("rtracklayer")
if ("GenomicRanges" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("GenomicRanges")
if ("plyranges" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("plyranges")
if ("GSEABase" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("GSEABase")
if ("GOSim" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("GOSim")
if ("stats" %in% rownames(installed.packages()) == 'FALSE') install.packages("stats")
if ("ggdendro" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("ggdendro")
if ("GO.db" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("GO.db")
if ("rrvgo" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("rrvgo")

#BiocManager::install("org.Ce.eg.db", force=TRUE) #install if needed 

#BiocManager::install("topGO")
#BiocManager::install("biomaRt")
#BiocManager::install("Rgraphviz")

library("tidyverse")
library("genefilter")
library("DESeq2")
library("RColorBrewer")
library("WGCNA")
library("flashClust")
library("gridExtra")
library("ComplexHeatmap")
library("goseq")
library("dplyr")
library("clusterProfiler")
library("pheatmap")
library("magrittr")
library("rtracklayer")
library("GenomicRanges")
library("plyranges")
library("GSEABase")
#library("GOSim")
library("stats")
library("ggdendro")
library("GO.db")
library("rrvgo")
library("cowplot")

library("topGO")
library("biomaRt")
library("Rgraphviz")
```
# Read in data files and filtering data

```
#Load metadata sheet with sample name and eaten vs control information
metadata <- read.csv("data/rna_seq/sample_rnaseq_metadata.csv", header = TRUE, sep = ",")%>%dplyr::select(sample, eaten_vs_control)
metadata$code<-paste0(metadata$eaten_vs_control)
head(metadata)
```
```
#Load A. hyacinthus gene count matrix generated from stringtie
gcount <- as.data.frame(read.csv("data/rna_seq/Acropora_gene_count_matrix.csv", row.names="gene_id"), colClasses = double)
head(gcount)
```
Check that there are no genes with 0 counts across all samples. 

```{r}
dim(gcount) 

gcount<-gcount %>%
     mutate(Total = rowSums(.[, 1:54]))%>%
    filter(!Total==0)%>%
    dplyr::select(!Total)

dim(gcount)
```
We started with X genes. which is the total number in the annotation. About X genes that had 0's across all samples (not detected) with X remaining.

Conduct data filtering, this includes:  

*pOverA*: Specifying the minimum count for a proportion of samples for each gene. Here, we are using a pOverA of 0.05. This is because we have 16 samples with a minimum of n=8 samples  per group. Therefore, we will accept genes that are present in 8/16 = 0.50 of the samples because we expect different expression by treatment group (eaten x non-eaten). We are further setting the minimum count of genes to 10, such that 50% of the samples must have a gene count of >10 in order for the gene to remain in the data set.  

Filter in the package "genefilter". Pre-filtering our dataset to reduce the memory size dataframe, increase the speed of the transformation and testing functions, and improve quality of statistical analysis by removing low-coverage counts. Removed counts could represent outliers in the data and removing these improves sensitivity of statistical tests.

```{r}
filt <- filterfun(pOverA(0.1,10))

#create filter for the counts data
gfilt <- genefilter(gcount, filt)

#identify genes to keep by count filter
gkeep <- gcount[gfilt,]

#identify gene lists
gn.keep <- rownames(gkeep)

#gene count data filtered in PoverA, P percent of the samples have counts over A
gcount_filt <- as.data.frame(gcount[which(rownames(gcount) %in% gn.keep),])

#How many rows do we have before and after filtering?
nrow(gcount) #Before
nrow(gcount_filt) #After
```

Before filtering, we had approximately X genes. After filtering for pOverA, we have approximately X genes. This indicates that there were X genes present in X% of samples at <10 counts per gene. X/X = 0.X, so I selected X% filter because I want to keep genes that are uniquely expressed in one treatment group (n=X out of X).  

Prepare sample names. Keep R and any numbers following R in the colnames of gcount_filt.  

```{r}
# Function to extract 'R' followed by numbers from column names
extract_R_numbers <- function(col_names) {
  sub("R(\\d+).*", "\\1", col_names)
}

# Apply the function to all column names
new_col_names <- sapply(names(gcount_filt), extract_R_numbers)

# Rename the columns
names(gcount_filt) <- paste0("R", new_col_names)

# Output the modified dataframe
names(gcount_filt)
length(names(gcount_filt))
```

Column names are now correct and there are X as expected.  

Print the sample names in the metadata file. 

```{r}
metadata$sample
```
