knitr::opts_chunk$set(echo = TRUE)
options(stringsAsFactors = FALSE)
BiocManager::install("org.Ce.eg.db", force=TRUE) #install if needed 
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
metadata_por <- read.csv("D:/RNAseq/metadata_porites.csv", header = TRUE, sep = ",")%>%dplyr::select(sample, eatenvscontrol, code)
metadata_por$code<-paste0(metadata_por$code)
head(metadata_por)

#Load Porites gene count matrix generated from stringtie
gcount_por <- as.data.frame(read.csv("D:/RNAseq/porites_gene_count_matrix.csv", row.names="gene_id"), colClasses = double)
head(gcount_por)
#Check that there are no genes with 0 counts across all samples (remember only 16 columns since we are just looking at Porites). 
dim(gcount_por) 
   #[1] 40389    16
gcount_por<-gcount_por %>%
  mutate(Total = rowSums(.[, 1:16]))%>%
  filter(!Total==0)%>%
  dplyr::select(!Total)
dim(gcount_por)
  #[1] 35599    16

#Create filter for the counts data
filt_por <- filterfun(pOverA(0.5,10))
gfilt_por <- genefilter(gcount_por, filt_por)
#identify genes to keep by count filter
gkeep_por <- gcount_por[gfilt_por,]
#identify gene lists
gn.keep_por <- rownames(gkeep_por)
#gene count data filtered in PoverA, P percent of the samples have counts over A
gcount_filt_por <- as.data.frame(gcount_por[which(rownames(gcount_por) %in% gn.keep_por),])

#How many rows do we have before and after filtering?
nrow(gcount_por) #Before
    #[1] 35599
nrow(gcount_filt_por) #After
    #[1] 16279

# Function to extract 'R' followed by numbers from column names
extract_R_numbers_por <- function(col_names) {
  sub("R(\\d+).*", "\\1", col_names)
}

# Apply the function to all column names
new_col_names_por <- sapply(names(gcount_filt_por), extract_R_numbers_por)

# Rename the columns - do I need to do this?
# names(gcount_filt_por) <- paste0("R", new_col_names_por)
# Output the modified dataframe
names(gcount_filt_por)
length(names(gcount_filt_por))
    #[1] 16
#Column names are now correct and there are 16 as expected.

# Print the sample names in the metadata_por file. 
metadata_por$sample<-as.factor(metadata_por$sample)
metadata_por$eatenvscontrol<-as.factor(metadata_por$eatenvscontrol)
metadata_por$code<-as.factor(metadata_por$code)

# Set levels of factors. 
metadata_por$eatenvscontrol<-factor(metadata_por$eatenvscontrol, levels=c("eaten", "control"))
metadata_por$code<-factor(metadata_por$code, levels=c("Por_eaten", "Por_control"))

#Make sure the metadata and the columns in the gene count matrix are all the same.  
metadata_por$sample
colnames(gcount_filt_por)
list_por<-colnames(gcount_filt_por)
list_por<-as.factor(list_por)

#Make sure metadata$sample and column names are identical
metadata_por$sample <- tolower(trimws(metadata_por$sample))
colnames(gcount_filt_por) <- tolower(trimws(colnames(gcount_filt_por)))
identical(metadata_por$sample, colnames(gcount_filt_por))
    #[1] TRUE

## Re-order the levels - havent included the second part as it was messing my code up!
metadata_por$sample<-as.factor(metadata_por$sample)
##metadata_por$sample <- factor(as.character(metadata_por$sample), levels=list)

# Re-order the data.frame
metadata_por_ordered <- metadata_por[order(metadata_por$sample),]
metadata_por_ordered$sample

head(metadata_por_ordered)

metadata_por_ordered$sample
colnames(gcount_filt_por)
head(gcount_filt_por)

gcount_filt_por <- gcount_filt_por[, order(colnames(gcount_filt_por))]

# Remove all spaces from the 'sample' column
metadata_por_ordered$sample <- gsub("\\s+", "", metadata_por_ordered$sample)
# Remove all spaces from the column names of 'gcount_filt_por'
colnames(gcount_filt_por) <- gsub("\\s+", "", colnames(gcount_filt_por))
#Check names are the same
identical(metadata_por_ordered$sample, colnames(gcount_filt_por))
    #[1] TRUE

## Check unique sample names from metadata
unique(metadata_por_ordered$sample)
# Check unique column names from gcount_filt_por
unique(colnames(gcount_filt_por))
# Clean up metadata sample names
metadata_por_ordered$sample <- tolower(trimws(metadata_por_ordered$sample))

# Clean up gcount column names
colnames(gcount_filt_por) <- tolower(trimws(colnames(gcount_filt_por)))
identical(metadata_por_ordered$sample, colnames(gcount_filt_por))
    #[1] TRUE
#Make sure that the reference level for estimating effect size is control
metadata_por$eatenvscontrol <- relevel(metadata_por$eatenvscontrol, ref = "control")
gdds_por <- DESeqDataSetFromMatrix(countData = gcount_filt_por,
                               colData = metadata_por,
                               design = ~eatenvscontrol)

SF.gdds_por <- estimateSizeFactors(gdds_por) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 for us to use vst
print(sizeFactors(SF.gdds_por)) #View size factors

all(sizeFactors(SF.gdds_por)) < 4
    #[1] TRUE (but look like we have 11 and 8 when we print the factors - ask Ariana)

#All size factors are less than 4, so we can use VST transformation.  
gvst_por <- vst(gdds_por, blind=FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
head(assay(gvst_por), 3) #view transformed gene count data for the first three genes in the dataset. 

dim(gvst_por)
    #[1] 16279    16

#Plot a heatmap to sample distances
gsampleDists_por <- dist(t(assay(gvst_por))) #calculate distance matix
gsampleDistMatrix_por <- as.matrix(gsampleDists_por) #distance matrix
rownames(gsampleDistMatrix_por) <- colnames(gvst_por) #assign row names
colnames(gsampleDistMatrix_por) <- NULL #assign col names
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

pht_por<-pheatmap(gsampleDistMatrix_por, #plot matrix
              clustering_distance_rows=gsampleDists_por, #cluster rows
              clustering_distance_cols=gsampleDists_por, #cluster columns
              col=colors) #set colors

save_pheatmap_pdf(pht, "D:/RNAseq/sample_distances_por.pdf")
dev.off()

### Conduct PERMANOVA and PCA for all genes
#Export data for PERMANOVA test.  
test_por<-t(assay(gvst_por)) #export as matrix
test_por<-as.data.frame(test_por)

#add category columns
test_por$sample<-rownames(test_por)
test_por$eatenvscontrol<-metadata_por$eatenvscontrol[match(test_por$sample, metadata_por$sample)]

#Build PERMANOVA model for all genes in the dataset.  
dim(test_por)
    #[1]    16 16281 (remember two of the values are for the right measurement are just names so actually = 29555)
scaled_test_por <-prcomp(test_por[c(1:16279)], scale=TRUE, center=TRUE)
fviz_eig(scaled_test_por)

# scale data
vegan_por <- scale(test_por[c(1:16279)])

# PerMANOVA 
permanova_por<-adonis2(vegan_por ~ eatenvscontrol, data = test_por, method='eu')
permanova_por
#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999
#adonis2(formula = vegan_por ~ eatenvscontrol, data = test_por, method = "eu")
#Df SumOfSqs     R2      F Pr(>F)  
#Model     1    38800 0.1589 2.6448  0.065 .
#Residual 14   205385 0.8411                
#Total    15   244185 1.0000                
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Examine PCA and sample distances of all genes  

#Plot a PCA of samples by eatenvscontrol phenotype for all genes.        
gPCAdata_por <- plotPCA(gvst_por, intgroup = c("eatenvscontrol"), returnData=TRUE)
percentVar_por <- round(100*attr(gPCAdata_por, "percentVar")) #plot PCA of samples with all data

allgenesfilt_PCA_por <- ggplot(gPCAdata_por, aes(PC1, PC2, color=eatenvscontrol)) + 
  geom_point(size=3) + 
  ggtitle(expression(bold("All genes (16279 genes)")))+
  geom_text(aes(x = PC1, y = PC2, label = name), vjust = -0.5) +
  xlab(paste0("PC1: ",percentVar_por[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_por[2],"% variance")) +
  scale_color_manual(name="Eaten vs Control", values=c("control"="gray", "eaten"="orange"))+
  geom_text(x=0, y=10, label="p[eatenvscontrol]=0.065", color="black")+
  theme_classic() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines 
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()); allgenesfilt_PCA_por

ggsave("D:/RNAseq/pca_por.pdf", allgenesfilt_PCA_por, width=8, height=4)
ggsave("D:/RNAseq/pca_por.jpeg", allgenesfilt_PCA_por, width=8, height=4)

#Add a plot with elipses.  
allgenes_ellipse_por<-allgenesfilt_PCA_por + stat_ellipse();allgenes_ellipse_por
ggsave("D:/RNAseq/pca_ellipse_por.pdf", allgenes_ellipse_por, width=9, height=5)
ggsave("D:/RNAseq/pca_ellipse_por.jpeg", allgenes_ellipse_por, width=9, height=5)
###### Remove left cluster as it is clusering only based on read depth (225R, 34R, 43R, 16R)





### Run DEG analysis of eatenvscontrol   

#Use the likelihood ratio approach based on guidance from DESeq2: "DESeq2 offers two kinds of hypothesis tests: the Wald test, where we use the estimated standard error of a log2 fold change to test if it is equal to zero, and the likelihood ratio test (LRT). The LRT examines two models for the counts, a full model with a certain number of terms and a reduced model, in which some of the terms of the full model are removed. The test determines if the increased likelihood of the data using the extra terms in the full model is more than expected if those extra terms are truly zero.
#The LRT is therefore useful for testing multiple terms at once, for example testing 3 or more levels of a factor at once, or all interactions between two variables. The LRT for count data is conceptually similar to an analysis of variance (ANOVA) calculation in linear regression, except that in the case of the Negative Binomial GLM, we use an analysis of deviance (ANODEV), where the deviance captures the difference in likelihood between a full and a reduced model." 

## 1. View DEG's by eatenvscontrol 
#Run Wald test.  
DEG_T_por <- DESeq(gdds_por, test="Wald")
#View results by eatenvscontrol effects. Create dataframe for each comparison and add a column for the specific contrast.  
resultsNames(DEG_T_por)
#this shows the comparisons made (shown similar to lm output in terms of interactions)
eatenvscontrol_por<-results(DEG_T_por, contrast=c("eatenvscontrol","control","eaten"))
eatenvscontrol_por <- as.data.frame(subset(eatenvscontrol_por, padj<0.01))
eatenvscontrol_por <- as.data.frame(subset(eatenvscontrol_por, abs(log2FoldChange)>1.5))
eatenvscontrol_por$contrast <- c("eatenvscontrol_por")
eatenvscontrol_por$gene <- rownames(eatenvscontrol_por)
rownames(eatenvscontrol_por) <- NULL
dim(eatenvscontrol_por)
# [1] 3   8

#Combine this list of DEGs into one dataframe. 
DEG_eatenvscontrol_por<-rbind(eatenvscontrol_por)
dim(DEG_eatenvscontrol_por)
length(unique(DEG_eatenvscontrol_por$gene))

#There are 30 total genes and 8 unique genes, indicating that some are shared in the contrasts. 
#Get a VST normalized gene count matrix for these DEGs. 
DEG_eatenvscontrol_vst_por <- gdds_por[unique(DEG_eatenvscontrol_por$gene)]
dim(DEG_eatenvscontrol_vst_por)

# Now apply variance stabilizing transformation
DEG_eatenvscontrol_vst_por <- varianceStabilizingTransformation(gcount_filt_por_filtered)


DEG_eatenvscontrol_vst_por <- varianceStabilizingTransformation(DEG_eatenvscontrol_vst_por) 
#Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
# every gene contains at least one zero, cannot compute log geometric means


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
class(DEG_eatenvscontrol_por)
DEG_eatenvscontrol_por <- DESeq(gdds_por)
#Eaten vs control
res_por <- results(DEG_eatenvscontrol_por, contrast=c("eatenvscontrol","control","eaten"))

volcano_eatenvscontrol_por<-EnhancedVolcano(res_por,
                                        lab = rownames(res_por),
                                        x = 'log2FoldChange',
                                        y = 'pvalue',
                                        title = 'eatenvscontrol',
                                        pCutoff = 0.01,
                                        FCcutoff = 1.5,
                                        pointSize = 3.0,
                                        labSize = 0.0, 
                                        legendPosition = 'none')
ggsave("D:/RNAseq/volcano_eatenvscontrol_por.pdf", plot=volcano_eatenvscontrol_por, width=18, height=6)
ggsave("D:/RNAseq/volcano_eatenvscontrol_por.png", plot=volcano_eatenvscontrol_por, width=18, height=6)

### Plot a heatmap for genes that are differentially expressed between eatenvscontrol.
#Set themes and metadata.  
colnames(colData(gdds_por))
nrow(colData(gdds_por))
ncol(gcount_filt_por)
legend_names_col = colnames(assay(DEG_eatenvscontrol_vst_por))

#try again
# Extract the 'eatenvscontrol' column as a vector from colData(gdds_por)
df_por <- data.frame(eatenvscontrol = colData(gdds_por)[, "eatenvscontrol"])

# Convert the 'eatenvscontrol' column to a factor with levels 'eaten' and 'control'
df_por$eatenvscontrol <- factor(as.character(df_por$eatenvscontrol), levels = c("eaten", "control"))

# Check the structure of df_por after conversion
str(df_por)

df_por<-df_por%>%
  arrange(eatenvscontrol)
ann_colors = list(
  eatenvscontrol = c("control" = "lightblue","eaten"= "red2") 
)

pdf("D:/RNAseq/DEG_eatenvscontrol_por.pdf", width=15, height=20)
porDEGS_heatmap <- pheatmap(assay(DEG_eatenvscontrol_vst_por), cluster_rows=TRUE, show_rownames=FALSE, color=inferno(10), show_colnames=TRUE, fontsize_row=3, scale="row", cluster_cols=TRUE, annotation_col=df, annotation_colors = ann_colors, labels_col=legend_names_col)
dev.off()
ggsave("D:/RNAseq/porDEGS_heatmap.pdf", plot=porDEGS_heatmap, width=18, height=6)

### Run a PERMANOVA and plot a PCA for temperature comparisons.  
#Conduct PERMANOVA by eatenvscontrol
#Export data for PERMANOVA test.  
test_DEG_eatenvscontrol_por<-t(assay(DEG_eatenvscontrol_vst_por)) #export as matrix
test_DEG_eatenvscontrol_por<-as.data.frame(test_DEG_eatenvscontrol_por)

#add category columns
test_DEG_eatenvscontrol_por$sample<-rownames(test_DEG_eatenvscontrol_por)
test_DEG_eatenvscontrol_por$eatenvscontrol<-metadata_por$eatenvscontrol[match(test_DEG_eatenvscontrol_por$sample, metadata_por$sample)]


#Build PERMANOVA model for DEG's that are different between temperatures.  
dim(test_DEG_eatenvscontrol_por)
scaled_test_DEG_eatenvscontrol_por <-prcomp(test_DEG_eatenvscontrol_por[c(1:30)], scale=TRUE, center=TRUE)
fviz_eig(scaled_test_DEG_eatenvscontrol_por)

# scale data
vegan_DEG_eatenvscontrol_por <- scale(test_DEG_eatenvscontrol_por[c(1:30)])

# PerMANOVA 
permanova_DEG<-adonis2(vegan_DEG_eatenvscontrol_por ~ eatenvscontrol, data = test_DEG_eatenvscontrol_por, method='eu')
permanova_DEG

#  Df SumOfSqs      R2      F Pr(>F)    
#Model     1    78.18 0.17372 2.9435  0.001 ***
#  Residual 14   371.82 0.82628                  
#Total    15   450.00 1.00000                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###Plot a PCA of DEG's for eatenvscontrol.   
#Plot PCA doesn't work as DEG_eatenvscontrol_vst_por doesn't work (line 307)
gPCAdata_DEG_eatenvscontrol_por <- plotPCA(DEG_eatenvscontrol_vst_por, intgroup = c("eatenvscontrol"), returnData=TRUE)
percentVar_DEG_eatenvscontrol_por <- round(100*attr(gPCAdata_DEG_eatenvscontrol_por, "percentVar")) #plot PCA of samples with all data

DEG_PCA_eatenvscontrol_por <- ggplot(gPCAdata_DEG_eatenvscontrol_por, aes(PC1, PC2, color=eatenvscontrol)) + 
  geom_point(size=3) + 
  ggtitle(expression(bold("DEGs:Porites Eaten vs Control (30 genes)")))+
  xlab(paste0("PC1: ",percentVar_DEG_eatenvscontrol_por[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_DEG_eatenvscontrol_por[2],"% variance")) +
  scale_color_manual(name="eatenvscontrol", values=c("control"="gray", "eaten"="orange"))+
  geom_text(x=0, y=-8, label="p[eatenvscontrol]=0.001", color="black")+
  theme_classic() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines 
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()); DEG_PCA_eatenvscontrol_por

ggsave("D:/RNAseq/DEG_PCA_eatenvscontrol_por.pdf", DEG_PCA_eatenvscontrol_por, width=8, height=4)
ggsave("D:/RNAseq/DEG_PCA_eatenvscontrol_por.png", DEG_PCA_eatenvscontrol_por, width=8, height=4)