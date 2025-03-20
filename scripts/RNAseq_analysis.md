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
```
knitr::opts_chunk$set(echo = TRUE)
options(stringsAsFactors = FALSE)
if ("tidyverse" %in% rownames(installed.packages()) == 'FALSE') install.packages('tidyverse') 
if ("genefilter" %in% rownames(installed.packages()) == 'FALSE') install.packages('genefilter') 
if ("RColorBrewer" %in% rownames(installed.packages()) == 'FALSE') install.packages('RColorBrewer') 
if ("WGCNA" %in% rownames(installed.packages()) == 'FALSE') install.packages('WGCNA') 
if ("flashClust" %in% rownames(installed.packages()) == 'FALSE') install.packages('flashClust') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("ComplexHeatmap" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("ComplexHeatmap")
if ("DESeq2" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("DESeq2")
if ("goseq" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("goseq") 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("clusterProfiler" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("clusterProfiler") 
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
BiocManager::install("org.Ce.eg.db", force=TRUE) #install if needed 
BiocManager::install("topGO")
BiocManager::install("biomaRt")
BiocManager::install("Rgraphviz")
BiocManager::install("EnhancedVolcano")
```

```
library("DESeq2")
library("tidyverse")
library("genefilter")
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
library("GOSim")
library("stats")
library("ggdendro")
library("GO.db")
library("rrvgo")
library("cowplot")
library("topGO")
library("biomaRt")
library("Rgraphviz")
library("EnhancedVolcano")
```

#Load metadata sheet with sample name and eaten vs control information
```
metadata <- read.csv("D:/RNAseq/sample_rnaseq_metadata_gorman.csv", header = TRUE, sep = ",")%>%dplyr::select(sample, eatenvscontrol, code)
metadata$code<-paste0(metadata$code)
head(metadata)
```
#Load A. hyacinthus gene count matrix generated from stringtie
```
gcount <- as.data.frame(read.csv("D:/RNAseq/acropora_gene_count_matrix.csv", row.names="gene_id"), colClasses = double)
head(gcount)
```
#Check that there are no genes with 0 counts across all samples (remember only 16 columns since we are just looking at Acropora!). 
```
dim(gcount) 
gcount<-gcount %>%
  mutate(Total = rowSums(.[, 1:16]))%>%
  filter(!Total==0)%>%
  dplyr::select(!Total)
dim(gcount)
```
```
filt <- filterfun(pOverA(0.1,10))
#create filter for the counts data
gfilt <- genefilter(gcount, filt)
```
#identify genes to keep by count filter
```
gkeep <- gcount[gfilt,]

#identify gene lists
gn.keep <- rownames(gkeep)
```
#gene count data filtered in PoverA, P percent of the samples have counts over A
```
gcount_filt <- as.data.frame(gcount[which(rownames(gcount) %in% gn.keep),])
```
#How many rows do we have before and after filtering?
```
nrow(gcount) #Before
nrow(gcount_filt) #After
# Function to extract 'R' followed by numbers from column names
extract_R_numbers <- function(col_names) {
  sub("R(\\d+).*", "\\1", col_names)
}
```

# Apply the function to all column names
```
new_col_names <- sapply(names(gcount_filt), extract_R_numbers)
```

# Rename the columns - do I need to do this?
```
#names(gcount_filt) <- paste0("R", new_col_names)
# Output the modified dataframe
names(gcount_filt)
length(names(gcount_filt))
#Column names are now correct and there are X as expected.
```

# Print the sample names in the metadata file. 
```
metadata$sample<-as.factor(metadata$sample)
metadata$eatenvscontrol<-as.factor(metadata$eatenvscontrol)
metadata$code<-as.factor(metadata$code)
```

# Set levels of factors. 
```
metadata$eatenvscontrol<-factor(metadata$eatenvscontrol, levels=c("eaten", "control"))
metadata$code<-factor(metadata$code, levels=c("Ahya_eaten", "Ahya_control"))
```
```
#Make sure the metadata and the columns in the gene count matrix are all the same.  
metadata$sample
colnames(gcount_filt)

list<-colnames(gcount_filt)
list<-as.factor(list)

#Make sure metadata$sample and column names are identical
metadata$sample <- tolower(trimws(metadata$sample))
colnames(gcount_filt) <- tolower(trimws(colnames(gcount_filt)))
identical(metadata$sample, colnames(gcount_filt))

metadata$sample<-as.factor(metadata$sample)
```
##Metadata 
```
## Re-order the levels - havent included this as it was messing my code up!
##metadata$sample <- factor(as.character(metadata$sample), levels=list)

# Re-order the data.frame
metadata_ordered <- metadata[order(metadata$sample),]
metadata_ordered$sample

head(metadata_ordered)

metadata_ordered$sample
colnames(gcount_filt)
head(gcount_filt)

gcount_filt <- gcount_filt[, order(colnames(gcount_filt))]
identical(metadata_ordered$sample, colnames(gcount_filt))
```
```
## Check unique sample names from metadata
unique(metadata_ordered$sample)
# Check unique column names from gcount_filt
unique(colnames(gcount_filt))
# Clean up metadata sample names
metadata_ordered$sample <- tolower(trimws(metadata_ordered$sample))
```

# Clean up gcount column names
```
colnames(gcount_filt) <- tolower(trimws(colnames(gcount_filt)))
identical(metadata_ordered$sample, colnames(gcount_filt))

gdds <- DESeqDataSetFromMatrix(countData = gcount_filt,
                               colData = metadata,
                               design = ~eatenvscontrol)
```
```
SF.gdds <- estimateSizeFactors(gdds) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 for us to use vst
print(sizeFactors(SF.gdds)) #View size factors

all(sizeFactors(SF.gdds)) < 4
```
All size factors are less than 4, so we can use VST transformation. 
```
gvst <- vst(gdds, blind=FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
head(assay(gvst), 3) #view transformed gene count data for the first three genes in the dataset.

   trim.321ra.gtf trim.331ra.gtf trim.336r.gtf trim.370r.gtf
Ahyacinthus26340       6.502666       6.867042      6.502666      6.502666
Ahyacinthus26341       7.217751       6.722764      6.502666      6.502666
Ahyacinthus26342       8.458109       6.980634      6.502666      6.502666
                 trim.380r.gtf trim.410r.gtf trim.414r.gtf trim.419r.gtf
Ahyacinthus26340      7.076405      6.502666      6.502666      7.518834
Ahyacinthus26341      6.502666      6.957145      7.237757      6.502666
Ahyacinthus26342      6.502666      7.385119      8.021587      6.502666
                 trim.468r.gtf trim.497r.gtf trim.512r.gtf trim.549r.gtf
Ahyacinthus26340      7.460916      8.548814      7.703149      6.502666
Ahyacinthus26341      6.502666      6.502666      6.950481      6.502666
Ahyacinthus26342      6.502666      6.502666      7.265970      6.502666
                 trim.568r.gtf trim.571r.gtf trim.581r.gtf trim.586r.gtf
Ahyacinthus26340      6.502666      6.942077      6.502666      6.502666
Ahyacinthus26341      6.502666      6.502666      6.502666      6.502666
Ahyacinthus26342      6.502666      7.039807      6.502666      6.502666

dim(gvst)
[1] 22719    16
```

#Plot a heatmap to sample distances
```
gsampleDists <- dist(t(assay(gvst))) #calculate distance matix
gsampleDistMatrix <- as.matrix(gsampleDists) #distance matrix
rownames(gsampleDistMatrix) <- colnames(gvst) #assign row names
colnames(gsampleDistMatrix) <- NULL #assign col names
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

pht<-pheatmap(gsampleDistMatrix, #plot matrix
              clustering_distance_rows=gsampleDists, #cluster rows
              clustering_distance_cols=gsampleDists, #cluster columns
              col=colors) #set colors

save_pheatmap_pdf(pht, "D:/RNAseq/sample_distances.pdf")
dev.off()
```

# Conduct PERMANOVA and PCA for all genes
```
#Export data for PERMANOVA test.  
test<-t(assay(gvst)) #export as matrix
test<-as.data.frame(test)

#add category columns
test$sample<-rownames(test)
test$eatenvscontrol<-metadata$eatenvscontrol[match(test$sample, metadata$sample)]


#Build PERMANOVA model for all genes in the dataset.  
install.packages('vegan') 
install.packages('factoextra') 
library(vegan)
library(factoextra)

dim(test)
scaled_test <-prcomp(test[c(1:22719)], scale=TRUE, center=TRUE)
fviz_eig(scaled_test)

# scale data
vegan <- scale(test[c(1:22719)])

# PerMANOVA 
permanova<-adonis2(vegan ~ eatenvscontrol, data = test, method='eu')
permanova

[1] Output
adonis2(formula = vegan ~ eatenvscontrol, data = test, method = "eu")
         Df SumOfSqs      R2      F Pr(>F)   
Model     1    44001 0.12912 2.0756  0.005 **
Residual 14   296784 0.87088                 
Total    15   340785 1.00000                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```


## Examine PCA and sample distances of all genes  
```
#Plot a PCA of samples by eatenvscontrol phenotype for all genes.        

gPCAdata <- plotPCA(gvst, intgroup = c("eatenvscontrol"), returnData=TRUE)
percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data

allgenesfilt_PCA <- ggplot(gPCAdata, aes(PC1, PC2, color=eatenvscontrol)) + 
  geom_point(size=3) + 
  ggtitle(expression(bold("All genes (22719 genes)")))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(name="Eaten vs Control", values=c("control"="gray", "eaten"="orange"))+
  geom_text(x=0, y=10, label="p[eatenvscontrol]=0.001", color="black")+
  theme_classic() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines 
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()); allgenesfilt_PCA

ggsave("D:/RNAseq/pca.pdf", allgenesfilt_PCA, width=8, height=4)
ggsave("D:/RNAseq/pca.jpeg", allgenesfilt_PCA, width=8, height=4)

#Add a plot with elipses.  

allgenes_ellipse<-allgenesfilt_PCA + stat_ellipse();allgenes_ellipse
ggsave("D:/RNAseq/pca_ellipse.pdf", allgenes_ellipse, width=9, height=5)
ggsave("D:/RNAseq/pca_ellipse.jpeg", allgenes_ellipse, width=9, height=5)
```

# Run DEG analysis of eatenvscontrol   
```
#Use the likelihood ratio approach based on guidance from DESeq2: "DESeq2 offers two kinds of hypothesis tests: the Wald test, where we use the estimated standard error of a log2 fold change to test if it is equal to zero, and the likelihood ratio test (LRT). The LRT examines two models for the counts, a full model with a certain number of terms and a reduced model, in which some of the terms of the full model are removed. The test determines if the increased likelihood of the data using the extra terms in the full model is more than expected if those extra terms are truly zero.

#The LRT is therefore useful for testing multiple terms at once, for example testing 3 or more levels of a factor at once, or all interactions between two variables. The LRT for count data is conceptually similar to an analysis of variance (ANOVA) calculation in linear regression, except that in the case of the Negative Binomial GLM, we use an analysis of deviance (ANODEV), where the deviance captures the difference in likelihood between a full and a reduced model." 

## 1. View DEG's by eatenvscontrol 

#Run Wald test.  
DEG_T <- DESeq(gdds, test="Wald")


#View results by eatenvscontrol effects. Create dataframe for each comparison and add a column for the specific contrast.  
resultsNames(DEG_T)
#this shows the comparisons made (shown similar to lm output in terms of interactions)

eatenvscontrol<-results(DEG_T, contrast=c("eatenvscontrol","control","eaten"))
eatenvscontrol <- as.data.frame(subset(eatenvscontrol, padj<0.01))
eatenvscontrol <- as.data.frame(subset(eatenvscontrol, abs(log2FoldChange)>1.5))
eatenvscontrol$contrast <- c("eatenvscontrol")
eatenvscontrol$gene <- rownames(eatenvscontrol)
rownames(eatenvscontrol) <- NULL

dim(eatenvscontrol)
# [1] 470   8
```

#There are 470 genes differentially expressed between
```
#Combine this list of DEGs into one dataframe. 
DEG_eatenvscontrol<-rbind(eatenvscontrol)
dim(DEG_eatenvscontrol)
length(unique(DEG_eatenvscontrol$gene))

#There are 470 total genes and X unique genes, indicating that some are shared in the contrasts. 
#Get a VST normalized gene count matrix for these DEGs. 

DEG_eatenvscontrol_vst <- gdds[unique(DEG_eatenvscontrol$gene)]
dim(DEG_eatenvscontrol_vst)

DEG_eatenvscontrol_vst <- varianceStabilizingTransformation(DEG_eatenvscontrol_vst) 
```

### Upset plot of genes between species at later R session - Acropora vs Porites
##Need more than one comparison to look at!
#Create an upset diagram showing genes unique and shared
```
#BiocManager::install("VennDetail")
#library("VennDetail")
#list1<-DEG_eatenvscontrol%>%
#filter(contrast=="species")
#venn <- venndetail(list(list1))
#upset_eatenvscontrol<-plot(venn, type = "venn")
#pdf("D:/RNAseq/upset_eatenvscontrol.pdf", width=8, height=6)
#par(mar = c(5, 5, 5, 5))
#plot(ven, type = "upset")
#dev.off()
```

### Plot a volcano plot for these comparisons.   
```
class(DEG_eatenvscontrol)
DEG_eatenvscontrol <- DESeq(gdds)
#Eaten vs control
res <- results(DEG_eatenvscontrol, contrast=c("eatenvscontrol","control","eaten"))

volcano_eatenvscontrol<-EnhancedVolcano(res,
                                     lab = rownames(res),
                                     x = 'log2FoldChange',
                                     y = 'pvalue',
                                     title = 'eatenvscontrol',
                                     pCutoff = 0.01,
                                     FCcutoff = 1.5,
                                     pointSize = 3.0,
                                     labSize = 0.0, 
                                     legendPosition = 'none')



ggsave("D:/RNAseq/volcano_eatenvscontrol.pdf", plot=volcano_eatenvscontrol, width=18, height=6)
ggsave("D:/RNAseq/volcano_eatenvscontrol.png", plot=volcano_eatenvscontrol, width=18, height=6)
```


(https://github.com/user-attachments/assets/88333f73-3b6e-4a16-998b-269ba2b6a551)
[acropora_gene_count_matrix.csv](https://github.com/user-attachments/files/19368199/acropora_gene_count_matrix.csv)
![volcano_eatenvscontrol](https://github.com/user-attachments/assets/6ec7d590-3a7b-4746-a385-68cb831ae3f9)
[volcano_eatenvscontrol.pdf](https://github.com/user-attachments/files/19368198/volcano_eatenvscontrol.pdf)
[sample_rnaseq_metadata_gorman.csv](https://github.com/user-attachments/files/19368197/sample_rnaseq_metadata_gorman.csv)
[sample_distances.pdf](https://github.com/user-attachments/files/19368196/sample_distances.pdf)
[pca_ellipse.pdf](https://github.com/user-attachments/files/19368194/pca_ellipse.pdf)
![pca_ellipse](https://github.com/user-attachments/assets/b4cc9046-0dad-48bb-9778-1d46b812ec51)
[pca.pdf](https://github.com/user-attachments/files/19368193/pca.pdf)

