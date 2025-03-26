#Remove 497R as it only mapped 2%, remove column from gcount matrix
#Remove 497R from metadata
#Rerun


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
if ("viridis" %in% rownames(installed.packages()) == 'FALSE') install.packages("viridis")
if ("topGO" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("topGO")
if ("biomaRt" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("biomaRt")
if ("Rgraphviz" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("Rgraphviz")
if ("EnhancedVolcano" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("EnhancedVolcano")
if ("vegan" %in% rownames(installed.packages()) == 'FALSE') install.packages("vegan") 
if ("factoextra" %in% rownames(installed.packages()) == 'FALSE') install.packages("factoextra") 

library("vegan")
library("factoextra")
library("viridis")
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

#Load metadata sheet with sample name and eaten vs control information
metadata <- read.csv("D:/RNAseq/sample_rnaseq_metadata_gorman.csv", header = TRUE, sep = ",")%>%dplyr::select(sample, eatenvscontrol, code)
metadata$code<-paste0(metadata$code)
head(metadata)
#Load A. hyacinthus gene count matrix generated from stringtie
gcount <- as.data.frame(read.csv("D:/RNAseq/acropora_gene_count_matrix.csv", row.names="gene_id"), colClasses = double)
head(gcount)
#Check that there are no genes with 0 counts across all samples (remember only 16 columns since we are just looking at Acropora!). 
dim(gcount)
    #[1] 27110    16
gcount<-gcount %>%
  mutate(Total = rowSums(.[, 1:16]))%>%
  filter(!Total==0)%>%
  dplyr::select(!Total)
dim(gcount)
    #[1] 25287    16

filt <- filterfun(pOverA(0.5,10))
#create filter for the counts data
gfilt <- genefilter(gcount, filt)
help(genefilter)
#identify genes to keep by count filter
gkeep <- gcount[gfilt,]

#identify gene lists
gn.keep <- rownames(gkeep)

#gene count data filtered in PoverA, P percent of the samples have counts over A
gcount_filt <- as.data.frame(gcount[which(rownames(gcount) %in% gn.keep),])

#How many rows do we have before and after filtering?
nrow(gcount) #Before
    #[1] 25287
nrow(gcount_filt) #After
    #[1] 17634
# Function to extract 'R' followed by numbers from column names
extract_R_numbers <- function(col_names) {
  sub("R(\\d+).*", "\\1", col_names)
}

# Apply the function to all column names
new_col_names <- sapply(names(gcount_filt), extract_R_numbers)

# Rename the columns - do I need to do this?
#names(gcount_filt) <- paste0("R", new_col_names)
# Output the modified dataframe
names(gcount_filt)
length(names(gcount_filt))
    #[1] 16
#Column names are now correct and there are 16 as expected.


# Print the sample names in the metadata file. 
metadata$sample<-as.factor(metadata$sample)
metadata$eatenvscontrol<-as.factor(metadata$eatenvscontrol)
metadata$code<-as.factor(metadata$code)

# Set levels of factors. 
metadata$eatenvscontrol<-factor(metadata$eatenvscontrol, levels=c("eaten", "control"))
metadata$code<-factor(metadata$code, levels=c("Ahya_eaten", "Ahya_control"))

#Make sure the metadata and the columns in the gene count matrix are all the same.  
metadata$sample
colnames(gcount_filt)

list<-colnames(gcount_filt)
list<-as.factor(list)

#Make sure metadata$sample and column names are identical
metadata$sample <- tolower(trimws(metadata$sample))
colnames(gcount_filt) <- tolower(trimws(colnames(gcount_filt)))
identical(metadata$sample, colnames(gcount_filt))
    #[1] TRUE
# Re-order the data.frame
metadata$sample<-as.factor(metadata$sample)
metadata_ordered <- metadata[order(metadata$sample),]
metadata_ordered$sample

head(metadata_ordered)

metadata_ordered$sample
colnames(gcount_filt)
head(gcount_filt)
#Check names are identical
gcount_filt <- gcount_filt[, order(colnames(gcount_filt))]
identical(metadata_ordered$sample, colnames(gcount_filt))
    #[1] FALSE
## Check unique sample names from metadata
unique(metadata_ordered$sample)
# Check unique column names from gcount_filt
unique(colnames(gcount_filt))
# Clean up metadata sample names
metadata_ordered$sample <- tolower(trimws(metadata_ordered$sample))

# Clean up gcount column names
colnames(gcount_filt) <- tolower(trimws(colnames(gcount_filt)))
identical(metadata_ordered$sample, colnames(gcount_filt))
    #[1] TRUE
#Make sure that the reference level for estimating effect size is control
metadata$eatenvscontrol <- relevel(metadata$eatenvscontrol, ref = "control")

gdds <- DESeqDataSetFromMatrix(countData = gcount_filt,
                               colData = metadata,
                               design = ~eatenvscontrol)

SF.gdds <- estimateSizeFactors(gdds) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 for us to use vst
print(sizeFactors(SF.gdds)) #View size factors

all(sizeFactors(SF.gdds)) < 4
    #[1] TRUE

#All size factors are less than 4, so we can use VST transformation.  
gvst <- vst(gdds, blind=FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
head(assay(gvst), 3) #view transformed gene count data for the first three genes in the dataset. 

dim(gvst)
    #[1] 17634    16

#Plot a heatmap to sample distances
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
## Need a diagnol line in the heatmap -shows  metadata is mapping correctly/ calling name of samples correctly
## not mapping correctly to sample name = no diagnol line 

# Conduct PERMANOVA and PCA for all genes
#Export data for PERMANOVA test.  
test<-t(assay(gvst)) #export as matrix
test<-as.data.frame(test)

#add category columns
test$sample<-rownames(test)
test$eatenvscontrol<-metadata$eatenvscontrol[match(test$sample, metadata$sample)]

dim(test)
    #[1]    16 17636
scaled_test <-prcomp(test[c(1:17634)], scale=TRUE, center=TRUE)
fviz_eig(scaled_test)

# scale data
vegan <- scale(test[c(1:17634)])

# PerMANOVA 
permanova<-adonis2(vegan ~ eatenvscontrol, data = test, method='eu')
permanova
#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999
#adonis2(formula = vegan ~ eatenvscontrol, data = test, method = "eu")
#Df SumOfSqs      R2      F Pr(>F)    
#Model     1    38615 0.14599 2.3932  0.001 ***
#  Residual 14   225895 0.85401                  
#Total    15   264510 1.00000                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## Examine PCA and sample distances of all genes  
#Plot a PCA of samples by eatenvscontrol phenotype for all genes.        

gPCAdata <- plotPCA(gvst, intgroup = c("eatenvscontrol"), returnData=TRUE)
percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data

allgenesfilt_PCA <- ggplot(gPCAdata, aes(PC1, PC2, color=eatenvscontrol)) + 
  geom_point(size=3) + 
  ggtitle(expression(bold("All genes (17634 genes)")))+
  geom_text(aes(x = PC1, y = PC2, label = name), vjust = -0.5) +
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
##Run analysis with and without data point PCA -> DEGS 

####

### Run DEG analysis of eatenvscontrol   
#Use the likelihood ratio approach based on guidance from DESeq2: "DESeq2 offers two kinds of hypothesis tests: the Wald test, where we use the estimated standard error of a log2 fold change to test if it is equal to zero, and the likelihood ratio test (LRT). The LRT examines two models for the counts, a full model with a certain number of terms and a reduced model, in which some of the terms of the full model are removed. The test determines if the increased likelihood of the data using the extra terms in the full model is more than expected if those extra terms are truly zero.

#The LRT is therefore useful for testing multiple terms at once, for example testing 3 or more levels of a factor at once, or all interactions between two variables. The LRT for count data is conceptually similar to an analysis of variance (ANOVA) calculation in linear regression, except that in the case of the Negative Binomial GLM, we use an analysis of deviance (ANODEV), where the deviance captures the difference in likelihood between a full and a reduced model." 
## 1. View DEG's by eatenvscontrol 

#Run Wald test.  
DEG_T <- DESeq(gdds, test="Wald")

#View results by eatenvscontrol effects. Create dataframe for each comparison and add a column for the specific contrast.  
resultsNames(DEG_T)
#this shows the comparisons made (shown similar to lm output in terms of interactions)
head(DEG_T)
eatenvscontrol<-results(DEG_T, contrast=c("eatenvscontrol","control","eaten"))
eatenvscontrol <- as.data.frame(subset(eatenvscontrol, padj<0.01))
eatenvscontrol <- as.data.frame(subset(eatenvscontrol, abs(log2FoldChange)>1.5))
eatenvscontrol$contrast <- c("eatenvscontrol")
eatenvscontrol$gene <- rownames(eatenvscontrol)
rownames(eatenvscontrol) <- NULL

#Ask Ariana how to append the column names from DEG_T to the gene results in the eatensvscontrol dataframe


dim(eatenvscontrol)
# [1] 481   8


#There are 481 genes differentially expressed between

#Combine this list of DEGs into one dataframe. 
DEG_eatenvscontrol<-rbind(eatenvscontrol)
DEG_eatenvscontrol

# Extract the 'sample' column from eatenvscontrol
sample_column <- eatenvscontrol$sample

# Add the 'sample' column to DEG_eatenvscontrol
DEG_eatenvscontrol$sample <- sample_column

# Check the result
head(DEG_eatenvscontrol)

dim(DEG_eatenvscontrol)
length(unique(DEG_eatenvscontrol$gene))

#Get a VST normalized gene count matrix for these DEGs. 

DEG_eatenvscontrol_vst <- gdds[unique(DEG_eatenvscontrol$gene)]
dim(DEG_eatenvscontrol_vst)

DEG_eatenvscontrol_vst <- varianceStabilizingTransformation(DEG_eatenvscontrol_vst) 

### Upset plot of genes between species at later R session - Acropora vs Porites
##Need more than one comparison to look at!
#Create an upset diagram showing genes unique and shared
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

### Plot a volcano plot for these comparisons.   
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


### Plot a heatmap for genes that are differentially expressed between eatenvscontrol.
#Set themes and metadata.  
colnames(colData(gdds))
nrow(colData(gdds))
ncol(gcount_filt)
df <- as.data.frame(colData(gdds)[, "eatenvscontrol"])
str(df)
##error here
legend_names_col = colnames(assay(DEG_eatenvscontrol_vst))

df$eatenvscontrol <- factor(as.character(df$eatenvscontrol), levels=c("eaten", "control"))

df<-df%>%
  arrange(eatenvscontrol)

ann_colors = list(
  eatenvscontrol = c("control" = "lightblue","eaten"= "red2"), 
)

pdf("D:/RNAseq/DEG_eatenvscontrol.pdf", width=15, height=20)
pheatmap(assay(DEG_eatenvscontrol_vst), cluster_rows=TRUE, show_rownames=FALSE, color=inferno(10), show_colnames=TRUE, fontsize_row=3, scale="row", cluster_cols=TRUE, annotation_col=df, annotation_colors = ann_colors, labels_col=legend_names_col)
dev.off()

### Run a PERMANOVA and plot a PCA for temperature comparisons.  
#Conduct PERMANOVA by eatenvscontrol
#Export data for PERMANOVA test.  
test_DEG_eatenvscontrol<-t(assay(DEG_eatenvscontrol_vst)) #export as matrix
test_DEG_eatenvscontrol<-as.data.frame(test_DEG_eatenvscontrol)

#add category columns
test_DEG_eatenvscontrol$sample<-rownames(test_DEG_eatenvscontrol)
test_DEG_eatenvscontrol$eatenvscontrol<-metadata$eatenvscontrol[match(test_DEG_eatenvscontrol$sample, metadata$sample)]


#Build PERMANOVA model for DEG's that are different between temperatures.  
dim(test_DEG_eatenvscontrol)

scaled_test_DEG_eatenvscontrol <-prcomp(test_DEG_eatenvscontrol[c(1:481)], scale=TRUE, center=TRUE)
fviz_eig(scaled_test_DEG_eatenvscontrol)

# scale data
vegan_DEG_eatenvscontrol <- scale(test_DEG_eatenvscontrol[c(1:481)])

# PerMANOVA 
permanova_DEG<-adonis2(vegan_DEG_eatenvscontrol ~ eatenvscontrol, data = test_DEG_eatenvscontrol, method='eu')
permanova_DEG

#Df SumOfSqs      R2      F Pr(>F)    
#Model     1   1543.8 0.21898 3.9252  0.001 ***
# Residual 14   5506.2 0.78102                  
#Total    15   7050.0 1.00000                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#About X% of variance is explained by first PC followed by <X% on following PCs. 

#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#Plot a PCA of DEG's for eatenvscontrol.        
gPCAdata_DEG_eatenvscontrol <- plotPCA(DEG_eatenvscontrol_vst, intgroup = c("eatenvscontrol"), returnData=TRUE)
percentVar_DEG_eatenvscontrol <- round(100*attr(gPCAdata_DEG_eatenvscontrol, "percentVar")) #plot PCA of samples with all data

DEG_PCA_eatenvscontrol <- ggplot(gPCAdata_DEG_eatenvscontrol, aes(PC1, PC2, color=eatenvscontrol)) + 
  geom_point(size=3) + 
  ggtitle(expression(bold("DEGs: Eaten vs Control (470 genes)")))+
  xlab(paste0("PC1: ",percentVar_DEG_eatenvscontrol[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_DEG_eatenvscontrol[2],"% variance")) +
  scale_color_manual(name="eatenvscontrol", values=c("control"="gray", "eaten"="orange"))+
  geom_text(x=0, y=-8, label="p[eatenvscontrol]=0.001", color="black")+
  theme_classic() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines 
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()); DEG_PCA_eatenvscontrol

ggsave("D:/RNAseq/DEG_PCA_eatenvscontrol.pdf", DEG_PCA_eatenvscontrol, width=8, height=4)
ggsave("D:/RNAseq/DEG_PCA_eatenvscontrol.png", DEG_PCA_eatenvscontrol, width=8, height=4)


#From this we can see very clear separation by eaten vs control with high variance explained by PC1.

### Plot a heatmap for genes that are differentially expressed between eatenvscontrol 

#Set themes and metadata.  
legend_names_col = colnames(assay(DEG_eatenvscontrol_vst))
df <- as.data.frame(colData(gdds)[, "eatenvscontrol", drop = FALSE])

df$eatenvscontrol <- factor(as.character(df$eatenvscontrol), levels=c("eaten", "control"))
df<-df%>%
  arrange(eatenvscontrol)

ann_colors = list(
  eatenvscontrol = c("control" = "lightblue","eaten"= "red2")
)


#Plot heatmap.  

pdf("D:/RNAseq/DEG_heatmap_eatenvscontrol_acro.pdf", width=15, height=20)
DEG_eatenvscontrol_acro <- pheatmap(assay(DEG_eatenvscontrol_vst), cluster_rows=TRUE, show_rownames=FALSE, color=inferno(10), show_colnames=TRUE, fontsize_row=3, scale="row", cluster_cols=TRUE, annotation_col=df, annotation_colors = ann_colors, labels_col=legend_names_col)
dev.off()
ggsave("D:/RNAseq/DEG_heatmap_eatenvscontrol_acro2.pdf", DEG_eatenvscontrol_acro, width=8, height=4)

### Plot a volcano plot for DEGs these comparisons.   
class(DEG_eatenvscontrol_vst)
DEG_eatenvscontrol_470 <- DESeq(gdds)
#Eaten vs control
res <- results(DEG_eatenvscontrol_vst, contrast=c("eatenvscontrol","control","eaten"))

volcano_eatenvscontrol_470<-EnhancedVolcano(res,
                                        lab = rownames(res),
                                        x = 'log2FoldChange',
                                        y = 'pvalue',
                                        title = 'eatenvscontrol',
                                        pCutoff = 0.01,
                                        FCcutoff = 1.5,
                                        pointSize = 3.0,
                                        labSize = 0.0, 
                                        legendPosition = 'none')



ggsave("D:/RNAseq/volcano_DEG_eatenvscontrol_470.pdf", plot=volcano_eatenvscontrol_470, width=18, height=6)
ggsave("D:/RNAseq/volcano_DEG_eatenvscontrol_470.png", plot=volcano_eatenvscontrol_470, width=18, height=6)

## Save lists of gene counts and DEG information 
#DEG_acro_list: List of DEG's and fold change/pvalues (eatenvscontrol)
#DEG_eatenvscontrol_vst: Count matrix of genes from DEG list 

deg_acro_list<-eatenvscontrol%>%
  dplyr::select(gene, everything())
write_csv(deg_acro_list, "D:/RNAseq/deg_acro_list.csv")

dim(DEG_eatenvscontrol_vst)
deg_acro_counts<-t(assay(DEG_eatenvscontrol_vst))
deg_acro_counts<-as.data.frame(deg_acro_counts)
dim(deg_acro_counts)
deg_acro_counts$sample<-rownames(deg_acro_counts)
deg_acro_counts<-deg_acro_counts%>%
  dplyr::select(sample, everything())
write_csv(deg_acro_counts, "D:/RNAseq/deg_acro_counts.csv")
