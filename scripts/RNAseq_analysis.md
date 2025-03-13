---
title: "RNAseq Differential Gene Expression Analysis"
author: "Lucy M Gorman"
date: "3/13/2025"
output: html_document
editor_options: 
  chunk_output_type: console
---

Based on code from Ariana S Huffmyer:
https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/rna-seq_DEG.Rmd
https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/rna-seq_functional_enrichment_topGO.Rmd

Functional enrichment of _Acropora hyacinthus_ and _Porites_ sp. eaten versus non-eaten by CoTS DEGs.  

# Set up 

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE) #Set Strings to character
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
#if ("GOSim" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("GOSim")
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
#Load A. hyacinthus gene count matrix generated from stringtie
```
gcount <- as.data.frame(read.csv("data/rna_seq/Acropora_gene_count_matrix.csv", row.names="gene_id"), colClasses = double)
head(gcount)
```
# Read in data files 
