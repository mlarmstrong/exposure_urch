##DOSE EXPOSURE STUDY SCRIPT ###
#Dose Exposure--Gene Expression Scripts
#website with tutorial information: 
#https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#rich-visualization-and-reporting-of-results

#Madison Armstrong 
#Last Modified 7/7/2025


setwd("~/Desktop/purp development/Dose_exposure")

#Libraries
#if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install(c("DESeq2","edgeR","arrayQualityMetrics", "variancePartition"))

# URL format: 
package_url <- "https://bioconductor.org/packages/tximport/"

# Install from the specified URL
install.packages(package_url, repos = NULL, type = "source")

library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(tximport)
library(variancePartition)
library(pheatmap)
library(ashr)

#figure theme to keep things consistent
theme_box <- function(base_size = 11, base_family = '') {
  theme_classic() %+replace% 
    theme(text = element_text(size = 20), 
          panel.border = element_rect(color = "black", fill = NA, size = 1)) 
}

#create a vector of filenames, reading in a table that contains the sample IDs
dir <- "~/Desktop/purp dev data/Dose_exposure/Salmon_quant"
files_list <- list.files(dir)
files_list <- files_list[-27]
files <- file.path(dir,files_list, "quant.sf")
all(file.exists(files))
names <- gsub("_quant|_USP.*", "", files_list)
names(files) <- names

## Import quantification tables 
#Read in matrix of RSEM expected read counts 
#all_quant.sf is normalized... we don't want that!!
all.quant <- read.delim("Dose_exposure/Salmon_quant/all_quant.sf")
tx2gene<-data.frame(all.quant[,1], all.quant[,1])
colnames(tx2gene)<-c("TXNAME", "GENEID")

data <- tximport(files, type="salmon", tx2gene=tx2gene)

metadata <- read.table(file.path("Dose_exposure/RNA_data/RNAdose_metadata.txt"), header=TRUE, row.names=1)

#make sure things are in the correct order
all(rownames(metadata) %in% colnames(data))
all(rownames(metadata) == colnames(data))

dds <- DESeqDataSetFromTximport(data,
                              colData = metadata,
                              design = ~ Matepair + Treatment) 

dds

#pre-filtering step: using same as Leslie
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# identify genes that pass expression cutoff
isexpr <- rowSums(fpm(dds) > 1) >= 0.5 * ncol(dds)

quantLog <- log2(fpm(dds)[isexpr, ] + 1)
# Check sample names alignment
all(colnames(quantLog) == rownames(metadata))  # Should return TRUE

#varpar analysis to understand what variation is driving the pcas
#https://www.bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/variancePartition.html

# Specify variables to consider
# Treatment is a fixed effect, main thing I am investigating
#Matepair and Stage are all categorical and fixed effects too 
form <- ~ Treatment + Matepair + Stage

# Fit model and extract results
# 1) fit linear mixed model on gene expression
# If categorical variables are specified, a linear mixed model is used
# If all variables are modeled as fixed effects, a linear model is used
# each entry in results is a regression model fit on a single gene
# 2) extract variance fractions from each model fit
# for each gene, returns fraction of variation attributable
#       to each variable
# Interpretation: the variance explained by each variables
# after correcting for all other variables
# Note that geneExpr can either be a matrix, and EList output by voom() in the limma package,or an ExpressionSet

# Re-run variance partitioning
varPart <- fitExtractVarPartModel(quantLog, form, metadata)

# sort variables (i.e. columns) by median fraction of variance explained
vp <- sortCols(varPart)

# Bar plot of variance fractions for the first 10 genes
plotPercentBars(vp[1:10, ], colors=my_colors)

# violin plot of contribution of each variable to total variance
plotVarPart(vp)
C <- canCorPairs(form, metadata)
# Plot correlation matrix between all pairs of variables
plotCorrMatrix(C)

##SEPARATE OUT STAGE FIRST####
####1. GASTRULA####
metadata.gast <- subset(metadata, Stage == "Gastrula")

#subset data
files_list_gast<-paste(rownames(metadata.gast))
files.gast <- file.path(dir,files_list_gast, "quant.sf")
all(file.exists(files.gast))
names.gast <- gsub("_quant|_USP.*", "", files_list_gast)
names(files.gast) <- names.gast
data.gast <- tximport(files.gast, type="salmon", tx2gene=tx2gene)

#With the count matrix, cts, and the sample information, coldata, we can construct a DESeqDataSet:
dds.gast <- DESeqDataSetFromTximport(data.gast,
                              colData = metadata.gast,
                              design = ~ Matepair + Treatment) 

dds.gast 

#pre-filtering step: using same as Leslie
keep <- rowSums(counts(dds.gast)) >= 10
dds.gast <- dds.gast[keep,]

#variance stabilizing transformation-- extracting transformed values
vsd.gast <- vst(dds.gast, blind=FALSE)

vsd.gast$Treatment <- factor(as.character(vsd.gast$Treatment), 
                        levels=c("C", "NP100", "NP500", "NP1000"))


###2. PLOTTING PCAS-GASTRULA####
pcaData<-
  plotPCA(vsd.gast, intgroup=c("Matepair", "Treatment"),  pcsToUse=1:2, returnData=TRUE) #using ntop=500 top features by variance
percentVar <- round(100 * attr(pcaData, "percentVar"))

#pca.gast<-
ggplot(pcaData, aes(x=PC1, y=PC2, color=Matepair, shape=Treatment)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  labs(title= "A. Gastrula")+
  scale_color_brewer(palette="Dark2")+
  theme(legend.position="right") + theme_box()

#linear models
##PC1
MPT1<-lm(PC1~Matepair+Treatment, data=pcaData)
##PC2
MPT3<-lm(PC2~Matepair+Treatment, data=pcaData)

summary(MPT1) #PC1 comparisons: MP3 & MP4 sig diff
summary (MPT2) #MP3 and NP500 sig diff
anova(MPT3) #matepair and treatment


#### Identify Differentially Expressed GENES ####  

#By default, DESeq2 perform pair-wise comparison of the first and the last variable
#in the experimental design variables and provide a result table
dose.gast <- DESeq(dds.gast)
# Report: Number of transcripts
num_transcripts <- nrow(dose.gast)
num_transcripts #10805

resultsNames(dose.gast) #comparing treatments to C and matepairs 3 vs 1 and 4 vs 1

#####control to treatment contrasts ####
###########NP100 vs C #####
res100<-results(dose.gast, contrast=c("Treatment", "C", "NP100"))
res100 #gives readout of padj, log2FoldChange and other stuff
summary(res100) 

# Subset colData to only include samples from the relevant contrast
subset_res100 <- colData(dose.gast)$Treatment %in% c("C", "NP100")
ordered_res100 <-res100[order(res100$padj),] #ordered by lowest pvalue

#plot top genes--
par(mfrow=c(3,5))
for(i in 1:15) {
  plotCounts(dose.gast[,subset_res100], gene=rownames(ordered_res100[i,]), intgroup="Treatment")
}

#top DIFFERENTIALLY EXPRESSED genes C vs NP100 heatmap
#filter by padj=0.05 and then extract logfoldchange
res100_filtered <- res100[!is.na(res100$padj) & res100$padj < 0.05, ]

# Subset the DESeq2 object to include only "Control" and "NP100"
dds_gast_100 <- dose.gast[, dose.gast$Treatment %in% c("C", "NP100")]
dds_gast_100$Treatment <- droplevels(dds_gast_100$Treatment)

# Remove genes with zero counts across all samples (to prevent vst() errors)
dds_gast_100 <- dds_gast_100[rowSums(counts(dds_gast_100)) > 0, ]

# Perform Variance Stabilizing Transformation (VST)
dose.gast.vsd1 <- vst(dds_gast_100, blind=FALSE)

# Order by absolute log2 fold change (biggest up/downregulation)
#extract log2foldchange from those genes... and insert a C=0 column
res100_top <- res100_filtered[order(abs(res100_filtered$log2FoldChange), decreasing=TRUE), ]

# Select the top 50 most differentially expressed genes (list)
res100_top_genes <- rownames(res100_top)[1:50]

# Extract the transformed expression matrix for the top genes
res100_top_genes_matrix_var <- assay(dose.gast.vsd1)[res100_top_genes, ]

df <- as.data.frame(colData(dose.gast.vsd1)[,c("Treatment","Matepair")])
treatment_factor <- as.numeric(factor(df$Treatment)) #make it sort by treatment for easier comparisons
dist_matrix <- dist(treatment_factor)

pheatmap(res100_top_genes_matrix_var, cluster_rows=TRUE, show_rownames=F,
         cluster_cols=TRUE, annotation_col=df, 
         scale="row,clustering_distance_cols=dist_matrix")


#log-fold change diffs MA-plot
#res100L<-lfcShrink(dose.gast, contrast=c("Treatment", "C", "NP100"), type="ashr")
#res100L
#plotMA(res100L, ylim=c(-2,2)) #MA-plot to show differences in log-fold change expression levels...

############NP500 vs C #####
res500<-results(dose.gast, contrast=c("Treatment", "C", "NP500"))
summary(res500) 
res500

# Subset colData to only include samples from the relevant contrast
subset_res500 <- colData(dose.gast)$Treatment %in% c("C", "NP500")
ordered_res500 <-res500[order(res500$padj),] #lowest pval

#plot top genes--
par(mfrow=c(3,5))
for(i in 1:15) {
  plotCounts(dose.gast[,subset_res500], gene=rownames(ordered_res500[i,]), intgroup="Treatment")
}

#top DIFFERENTIALLY EXPRESSED genes C vs NP heatmap
#filter by padj=0.05 and then extract logfoldchange
res500_filtered <- res500[!is.na(res500$padj) & res500$padj < 0.05, ]

# Subset the DESeq2 object to include only "Control" and "NP500"
dds_gast_500 <- dose.gast[, dose.gast$Treatment %in% c("C", "NP500")]
dds_gast_500$Treatment <- droplevels(dds_gast_500$Treatment)
table(colData(dds_gast_500)$Treatment) #verify you are pulling the right data
# Remove genes with zero counts across all samples (to prevent vst() errors)
dds_gast_500 <- dds_gast_500[rowSums(counts(dds_gast_500)) > 0, ]

# Perform Variance Stabilizing Transformation (VST)
dose.gast.vsd2 <- vst(dds_gast_500, blind=FALSE)

# Order by absolute log2 fold change (biggest up/downregulation)
#extract log2foldchange from those genes... and insert a C=0 column
res500_top <- res500_filtered[order(abs(res500_filtered$log2FoldChange), decreasing=TRUE), ]

# Select the top 100 most differentially expressed genes (list)
res500_top_genes <- rownames(res500_top)[1:50]

# Extract the transformed expression matrix for the top genes
res500_top_genes_matrix_var <- assay(dose.gast.vsd2)[res500_top_genes, ]

df <- as.data.frame(colData(dose.gast.vsd2)[,c("Treatment","Matepair")])
treatment_factor <- as.numeric(factor(df$Treatment)) #make it sort by treatment for easier comparisons
dist_matrix <- dist(treatment_factor)

pheatmap(res500_top_genes_matrix_var, cluster_rows=TRUE, show_rownames=F,
         cluster_cols=TRUE, annotation_col=df, 
         scale="row",clustering_distance_cols=dist_matrix)



############NP1000 vs C ####
res1000<-results(dose.gast, contrast=c("Treatment", "C", "NP1000"))
summary(res1000) 

# Subset colData to only include samples from the relevant contrast
subset_res1000 <- colData(dose.gast)$Treatment %in% c("C", "NP1000")

ordered_res1000 <-res1000[order(res1000$padj),]

#plot top genes--
par(mfrow=c(3,5))
for(i in 1:15) {
  plotCounts(dose.gast[,subset_res1000], gene=rownames(ordered_res1000[i,]), intgroup="Treatment")
}

#top DIFFERENTIALLY EXPRESSED genes C vs NP heatmap
#filter by padj=0.05 and then extract logfoldchange
res1000_filtered <- res1000[!is.na(res1000$padj) & res1000$padj < 0.05, ]

# Subset the DESeq2 object to include only "Control" and "NP100"
dds_gast_1000 <- dose.gast[, dose.gast$Treatment %in% c("C", "NP1000")]
dds_gast_1000$Treatment <- droplevels(dds_gast_1000$Treatment)
table(colData(dds_gast_1000)$Treatment) #verify you are pulling the right data
# Remove genes with zero counts across all samples (to prevent vst() errors)
dds_gast_1000 <- dds_gast_1000[rowSums(counts(dds_gast_1000)) > 0, ]

# Perform Variance Stabilizing Transformation (VST)
dose.gast.vsd3 <- vst(dds_gast_1000, blind=FALSE)

# Order by absolute log2 fold change (biggest up/downregulation)
#extract log2foldchange from those genes... and insert a C=0 column
res1000_top <- res1000_filtered[order(abs(res1000_filtered$log2FoldChange), decreasing=TRUE), ]

# Select the top 100 most differentially expressed genes (list)
res1000_top_genes <- rownames(res1000_top)[1:50]

# Extract the transformed expression matrix for the top genes
res1000_top_genes_matrix_var <- assay(dose.gast.vsd3)[res1000_top_genes, ]

df <- as.data.frame(colData(dose.gast.vsd3)[,c("Treatment","Matepair")])
treatment_factor <- as.numeric(factor(df$Treatment)) #make it sort by treatment for easier comparisons
dist_matrix <- dist(treatment_factor)

pheatmap(res1000_top_genes_matrix_var, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, 
         scale="row",clustering_distance_cols=dist_matrix)

#####treatment to treatment contrasts####
#NP100 vs NP1000
res1<-results(dose.gast, contrast=c("Treatment", "NP100", "NP1000"))
summary(res1) 


#NP500 vs NP1000
res5<-results(dose.gast, contrast=c("Treatment", "NP500", "NP1000"))
summary(res5) 

res15<-results(dose.gast, contrast=c("Treatment", "NP100", "NP500"))
summary(res15)

###3. LRT Gastrula--DEGs all treatments######
# The full model was specified previously with the `design = ~ sampletype`:

# Likelihood ratio test
#this tests if any of the predictors, treatment or mate pair, signifciantly contribute to gene expression
dds.gast_lrt <- DESeq(dds.gast, test="LRT", reduced = ~ Matepair) 
#null=only MP, LRT tests where treatment significantly improves fit
# Extract results
res.gast_LRT <- results(dds.gast_lrt)
res.gast_LRT
#####SAVE FOR TOPGO
write.csv(res.gast_LRT,file="Dose_exposure/RNA_data/results.gast_LRT.csv", row.names=TRUE)

# Subset the LRT results to return genes with padj < 0.05
#The p-values are determined solely by the difference in deviance between the ‘full’ and ‘reduced’ model formula
#(not log2 fold changes)
sig_res.gast_LRT <- res.gast_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
# Get sig gene lists
sigLRT_gast.genes <- sig_res.gast_LRT %>% 
  pull(gene)
length(sigLRT_gast.genes) #171 genes woo

# Perform Variance Stabilizing Transformation (VST)
dds.gast_lrt.vsd <- vst(dds.gast_lrt, blind=FALSE)

# Extract the transformed expression matrix for the top genes
sigLRT_gast.genes_matrix_var <- assay(dds.gast_lrt.vsd)[sigLRT_gast.genes, ]

df <- as.data.frame(colData(dds.gast_lrt.vsd)[,c("Treatment","Matepair")])
#reorder everything
df$Treatment <- factor(df$Treatment, levels = c("C", "NP100", "NP500", "NP1000")) 
df$Matepair <- factor(df$Matepair, levels = c("MP1", "MP2", "MP3", "MP4")) 
ordered_columns <- order(df$Treatment, df$Matepair)
sigLRT_gast.genes_matrix_var <- sigLRT_gast.genes_matrix_var[, ordered_columns]
df <- df[ordered_columns, ]  # Ensure annotation dataframe is in the same order

treatment_factor <- as.numeric(factor(df$Treatment)) #make it sort by treatment for easier comparisons
dist_matrix <- dist(treatment_factor)

#colors
annot_colors=list(
  Matepair=c("MP1"="#1b9e77", "MP2"="#D95F02", "MP3"="#7570b3", "MP4"="#e7298a"),
  Treatment=c("C"= "#A6CEE3" , "NP100"= "#1F78B4", "NP500"= "#33A02C", "NP1000"="#006D2C")
)

library(viridis)
gast.heat<-
  pheatmap(sigLRT_gast.genes_matrix_var, cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df, 
         annotation_colors=annot_colors,
         scale="row",clustering_distance_cols=dist_matrix,
         cutree_cols=4,
         color = magma(100))
  
ggsave("gast.heatmap.png",gast.heat, width=20, height=16, units = "cm") 


##4. PLUTEUS####
#Read in matrix of RSEM expected read counts 
metadata.plut <- subset(metadata, Stage == "Pluteus")
#subset data
files_list_plut<-paste(rownames(metadata.plut))
files.plut <- file.path(dir,files_list_plut, "quant.sf")
all(file.exists(files.plut))
names.plut <- gsub("_quant|_USP.*", "", files_list_plut)
names(files.plut) <- names.plut
data.plut <- tximport(files.plut, type="salmon", tx2gene=tx2gene)


#With the count matrix, cts, and the sample information, coldata, we can construct a DESeqDataSet:
dds.plut <- DESeqDataSetFromTximport(data.plut,
                                   colData = metadata.plut,
                                   design = ~ Matepair + Treatment) 

dds.plut 

#pre-filtering step: using same as Leslie
keep <- rowSums(counts(dds.plut)) >= 10
dds.plut <- dds.plut[keep,]

#variance stabilizing transformation-- extracting transformed values
vsd.plut <- vst(dds.plut, blind=FALSE)

vsd.plut$Treatment <- factor(as.character(vsd.plut$Treatment), 
                             levels=c("C", "NP100", "NP500", "NP1000"))

###5. PLOTTING PCAS-Pluteus####
pcaData1<-
  plotPCA(vsd.plut, intgroup=c("Matepair", "Treatment"),  
          pcsToUse=1:2, returnData=TRUE) #using ntop=500 top features by variance
percentVar2 <- round(100 * attr(pcaData1, "percentVar"))


#pca.plut<-
ggplot(pcaData1, aes(x=PC1, y=PC2, color=Matepair, shape=Treatment)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar2[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar2[2],"% variance")) + 
  labs(title= "B. Pluteus")+
  scale_color_brewer(palette="Dark2")+
  theme(legend.position="right")+ theme_box()

######6. PCA figures for paper####
library(gridExtra)
#pc1:2 gast and plut and with MP2 dropped
fullpca12 <- 
pca.gast + pca.plut + plot_layout(guides="collect")
ggsave("fullpca12.png",fullpca12, width=30, height=18, units = "cm")


#PC~linear models for pluteus: rerun code above to ensure correct axes
pMPT1<-lm(PC1~Matepair+Treatment, data=pcaData1)
##PC2
pMPT3<-lm(PC2~Matepair+Treatment, data=pcaData1)
summary(pMPT) #PC1 comparisons: MP3 & MP4 sig diff (like Gastrula)
summary (pMPT2) #MP2&3 sig diff (different than Gastrula)
anova(pMPT3) #matepair highly sig

#### Identify Differentially Expressed GENES ####  

#By default, DESeq2 perform pair-wise comparison of the first and the last variable
#in the experimental design variables and provide a result table
dose.plut <- DESeq(dds.plut)
# Report: Number of transcripts
num_transcriptsp <- nrow(dose.plut)
num_transcriptsp #14,567

resultsNames(dose.plut) #comparing treatments to C and matepairs 3 vs 1 and 4 vs 1

# Subsetting dds object and performing DE analysis for each condition
#~100 DEGs for different comparsions for pluteus...

#####control to treatment comparison####
############NP100 vs C ####
resp100<-results(dose.plut, contrast=c("Treatment", "C", "NP100"))
summary(resp100) 
resp100 #check comparison is correct

# Subset colData to only include samples from the relevant contrast
subset_resp500 <- colData(dose.plut)$Treatment %in% c("C", "NP100")

ordered_resp100 <-resp100[order(resp100$padj),]

#plot top genes--
par(mfrow=c(3,5))
for(i in 1:15) {
  plotCounts(dose.plut[,subset_resp100], gene=rownames(ordered_resp100[i,]), intgroup="Treatment")
}

#top DIFFERENTIALLY EXPRESSED genes C vs NP heatmap
#filter by padj=0.05 and then extract logfoldchange
resp100_filtered <- resp100[!is.na(resp100$padj) & resp100$padj < 0.05, ]

# Subset the DESeq2 object to include only "Control" and "NP100"
dds_plut_100 <- dose.plut[, dose.plut$Treatment %in% c("C", "NP100")]
dds_plut_100$Treatment <- droplevels(dds_plut_100$Treatment)
table(colData(dds_plut_100)$Treatment) #verify you are pulling the right data
# Remove genes with zero counts across all samples (to prevent vst() errors)
dds_plut_100 <- dds_plut_100[rowSums(counts(dds_plut_100)) > 0, ]

# Perform Variance Stabilizing Transformation (VST)
dose.plut.vsd1 <- vst(dds_plut_100, blind=FALSE)

# Order by absolute log2 fold change (biggest up/downregulation)
#extract log2foldchange from those genes... and insert a C=0 column
resp100_top <- resp100_filtered[order(abs(resp100_filtered$log2FoldChange), decreasing=TRUE), ]

# Select the top 100 most differentially expressed genes (list)
resp100_top_genes <- rownames(resp100_top)[1:50]

# Extract the transformed expression matrix for the top genes
resp100_top_genes_matrix_var <- assay(dose.plut.vsd1)[resp100_top_genes, ]

df <- as.data.frame(colData(dose.plut.vsd1)[,c("Treatment","Matepair")])
treatment_factor <- as.numeric(factor(df$Treatment)) #make it sort by treatment for easier comparisons
dist_matrix <- dist(treatment_factor)

pheatmap(resp100_top_genes_matrix_var, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, 
         scale="row",clustering_distance_cols=dist_matrix)

############NP500 vs C ####
resp500<-results(dose.plut, contrast=c("Treatment", "C", "NP500"))
summary(resp500) 
# Subset colData to only include samples from the relevant contrast
subset_resp500 <- colData(dose.plut)$Treatment %in% c("C", "NP500")

ordered_resp500 <-resp500[order(resp500$padj),]

#plot top genes–
par(mfrow=c(3,5))
for(i in 1:15) {
  plotCounts(dose.plut[,subset_resp500], gene=rownames(ordered_resp500[i,]), intgroup="Treatment")
}

#top DIFFERENTIALLY EXPRESSED genes C vs NP heatmap
#filter by padj=0.05 and then extract logfoldchange
resp500_filtered <- resp500[!is.na(resp500$padj) & resp500$padj < 0.05, ]

# Subset the DESeq2 object to include only "Control" and "NP500"
dds_plut_500 <- dose.plut[, dose.plut$Treatment %in% c("C", "NP500")]
dds_plut_500$Treatment <- droplevels(dds_plut_500$Treatment)
table(colData(dds_plut_500)$Treatment) #verify you are pulling the right data
# Remove genes with zero counts across all samples (to prevent vst() errors)
dds_plut_500 <- dds_plut_500[rowSums(counts(dds_plut_500)) > 0, ]

# Perform Variance Stabilizing Transformation (VST)
dose.plut.vsd2 <- vst(dds_plut_500, blind=FALSE)

# Order by absolute log2 fold change (biggest up/downregulation)
#extract log2foldchange from those genes... and insert a C=0 column
resp500_top <- resp500_filtered[order(abs(resp500_filtered$log2FoldChange), decreasing=TRUE), ]

# Select the top 100 most differentially expressed genes (list)
resp500_top_genes <- rownames(resp500_top)[1:50]

# Extract the transformed expression matrix for the top genes
resp500_top_genes_matrix_var <- assay(dose.plut.vsd2)[resp500_top_genes, ]

df <- as.data.frame(colData(dose.plut.vsd2)[,c("Treatment","Matepair")])
treatment_factor <- as.numeric(factor(df$Treatment)) #make it sort by treatment for easier comparisons
dist_matrix <- dist(treatment_factor)

pheatmap(resp500_top_genes_matrix_var, cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df, 
         scale="row",clustering_distance_cols=dist_matrix)


#######NP1000 vs C#######
resp1000<-results(dose.plut, contrast=c("Treatment", "C", "NP1000"))
summary(resp1000) 

# Subset colData to only include samples from the relevant contrast
subset_resp1000 <- colData(dose.plut)$Treatment %in% c("C", "NP1000")
ordered_resp1000 <-resp1000[order(resp1000$padj),]

#plot top genes–
par(mfrow=c(3,5))
for(i in 1:15) {
  plotCounts(dose.plut[,subset_resp1000], gene=rownames(ordered_resp1000[i,]), intgroup="Treatment")
}

#top DIFFERENTIALLY EXPRESSED genes C vs NP heatmap
#filter by padj=0.05 and then extract logfoldchange
resp1000_filtered <- resp1000[!is.na(resp1000$padj) & resp1000$padj < 0.05, ]

# Subset the DESeq2 object to include only "Control" and "NP1000"
dds_plut_1000 <- dose.plut[, dose.plut$Treatment %in% c("C", "NP1000")]
dds_plut_1000$Treatment <- droplevels(dds_plut_1000$Treatment)
table(colData(dds_plut_1000)$Treatment) #verify you are pulling the right data
# Remove genes with zero counts across all samples (to prevent vst() errors)
dds_plut_1000 <- dds_plut_1000[rowSums(counts(dds_plut_1000)) > 0, ]

# Perform Variance Stabilizing Transformation (VST)
dose.plut.vsd3 <- vst(dds_plut_1000, blind=FALSE)

# Order by absolute log2 fold change (biggest up/downregulation)
#extract log2foldchange from those genes... and insert a C=0 column
resp1000_top <- resp1000_filtered[order(abs(resp1000_filtered$log2FoldChange), decreasing=TRUE), ]

# Select the top 100 most differentially expressed genes (list)
resp1000_top_genes <- rownames(resp1000_top)[1:50]

# Extract the transformed expression matrix for the top genes
resp1000_top_genes_matrix_var <- assay(dose.plut.vsd3)[resp1000_top_genes, ]

df <- as.data.frame(colData(dose.plut.vsd3)[,c("Treatment","Matepair")])
treatment_factor <- as.numeric(factor(df$Treatment)) #make it sort by treatment for easier comparisons
dist_matrix <- dist(treatment_factor)

pheatmap(resp1000_top_genes_matrix_var, cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df, 
         scale="row",clustering_distance_cols=dist_matrix)

#####treatment to treatment#####
#NP100 vs NP1000
resp1<-results(dose.plut, contrast=c("Treatment", "NP100", "NP1000"))
summary(resp1) 

resp5<-results(dose.plut, contrast=c("Treatment", "NP500", "NP1000"))
summary(resp5) 

resp15<-results(dose.plut, contrast=c("Treatment", "NP100", "NP500"))
summary(resp15) 

###7. LRT Pluteus--DEGs all treatments######
# The full model was specified previously with the `design = ~ sampletype`:

# Likelihood ratio test
#this tests if any of the predictors, treatment or mate pair, signifciantly contribute to gene expression
dds.plut_lrt <- DESeq(dds.plut, test="LRT", reduced = ~ Matepair) 
#null=only MP, LRT tests where treatment signficantly improves fit

# Extract results
res.plut_LRT <- results(dds.plut_lrt)
res.plut_LRT
#####SAVE FOR TOPGO
write.csv(res.plut_LRT,file="Dose_exposure/RNA_data/results.plut_LRT.csv", row.names=TRUE)


# Subset the LRT results to return genes with padj < 0.05
#The p-values are determined solely by the difference in deviance between the ‘full’ and ‘reduced’ model formula
#(not log2 fold changes)
sig_res.plut_LRT <- res.plut_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
# Get sig gene lists
sig_res.plut_LRT.genes <- sig_res.plut_LRT %>% 
  pull(gene)
length(sig_res.plut_LRT.genes) #51 genes woo

# Perform Variance Stabilizing Transformation (VST)
dds.plut_lrt.vsd <- vst(dds.plut_lrt, blind=FALSE)

# Extract the transformed expression matrix for the top genes
sigLRT_plut.genes_matrix_var <- assay(dds.plut_lrt.vsd)[sig_res.plut_LRT.genes, ]

df1 <- as.data.frame(colData(dds.plut_lrt.vsd)[,c("Treatment","Matepair")])
#reorder everything
df1$Treatment <- factor(df1$Treatment, levels = c("C", "NP100", "NP500", "NP1000")) 
df1$Matepair <- factor(df1$Matepair, levels = c("MP1", "MP2", "MP3", "MP4")) 
ordered_columns <- order(df1$Treatment, df1$Matepair)
sigLRT_plut.genes_matrix_var <- sigLRT_plut.genes_matrix_var[, ordered_columns]
df1 <- df1[ordered_columns, ]  # Ensure annotation dataframe is in the same order

treatment_factor1 <- as.numeric(factor(df1$Treatment)) #make it sort by treatment for easier comparisons
dist_matrix1 <- dist(treatment_factor1)

#colors
annot_colors=list(
  Matepair=c("MP1"="#1b9e77", "MP2"="#D95F02", "MP3"="#7570b3", "MP4"="#e7298a"),
  Treatment=c("C"= "#A6CEE3" , "NP100"= "#1F78B4", "NP500"= "#33A02C", "NP1000"="#006D2C")
)


plut.heat<-
  pheatmap(sigLRT_plut.genes_matrix_var, cluster_rows=TRUE, show_rownames=FALSE,
           cluster_cols=TRUE, annotation_col=df1, 
           annotation_colors=annot_colors,
           scale="row",clustering_distance_cols=dist_matrix1,
           cutree_cols=4,border_color=NA,
           color = magma(100))

ggsave("plut.heatmap.png",plut.heat, width=20, height=16, units = "cm") 
