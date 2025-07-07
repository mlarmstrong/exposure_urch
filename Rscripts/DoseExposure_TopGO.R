##DOSE EXPOSURE STUDY SCRIPT ###
#TopGO script
#Madison Armstrong 
#Last Modified 7/7/2025


setwd("~/Desktop/purp dev data/Dose_exposure")
#Libraries
if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("biomaRt","topGO", "stringr"))

library(biomaRt)
##https://www.ensembl.org/info/data/biomart/index.html
library(topGO)
library(tidyverse)
library(stringr)
library(dplyr)

#read in files####
#read in .csv from LRT 

###gastrula####
gast.results<-read.csv("Dose_exposure/RNA_data/results.gast_LRT.csv")
parname <- "padj"  # Corrected column name
# Define significance threshold (typically 0.05 for DEGs)
significance_threshold <- 0.05  
# Create a new column to indicate significant genes (1) and non-significant genes (0)
gast.results$significant <- ifelse(gast.results[[parname]] < significance_threshold, 1, 0)
parsig.gast = which(gast.results$significant==1)
length(parsig.gast) #171 yay it worked
gast.results.sig <- gast.results[!is.na(gast.results$significant) & gast.results$significant != 0, ] #get rid of NAs & 0s

###pluteus####
plut.results<-read.csv("Dose_exposure/RNA_data/results.plut_LRT.csv")
parnames <- "padj"  # Corrected column name
# Define significance threshold (typically 0.05 for DEGs)
significance_threshold <- 0.05  
# Create a new column to indicate significant genes (1) and non-significant genes (0)
plut.results$significant <- ifelse(plut.results[[parnames]] < significance_threshold, 1, 0)
parsig.plut = which(plut.results$significant==1)
length(parsig.plut) #51 yay it worked
plut.results.sig <- plut.results[!is.na(plut.results$significant) & plut.results$significant != 0, ] #get rid of NAs & 0s

#Get annotations from biomaRt####
metazoa_mart <- useMart(biomart = 'metazoa_mart',
                        host = "https://metazoa.ensembl.org/")
Spurp_mart <- useDataset(mart = metazoa_mart, dataset = "spurpuratus_eg_gene")
listAttributes(Spurp_mart)
Spurp_genes <- getBM(attributes = c("ensembl_gene_id", "description", "chromosome_name", "start_position", "end_position",
                                    "external_gene_name", "gene_biotype", "go_id", "name_1006", "definition_1006", 
                                    "namespace_1003","ensembl_transcript_id"),
                     filters = c("with_geneid"),
                     values = TRUE,
                     mart = Spurp_mart)
head(Spurp_genes)
write.csv(Spurp_genes,"Dose_exposure/RNA_data/Spurp_genes2.csv")

### Create GO Key for all transcripts in the whole genome
uniqTranscripts <- unique(Spurp_genes$ensembl_transcript_id)  # Get unique transcript IDs
gos <- data.frame(ensembl_id = uniqTranscripts, GO=NA)  # Initialize data frame

for (t in uniqTranscripts) {
  sub <- subset(Spurp_genes, ensembl_transcript_id == t)  # Subset data for each transcript
  gostring <- paste(sub$go_id, collapse=",")  # Concatenate GO terms
  gos$GO[match(t, gos$ensembl_id)] <- gostring # Assign GO terms to the transcript
}

write.table(gos,"Dose_exposure/RNA_data/Spurp_GOmap.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

#use S.purp dataset to pull gene_info to fill in gene info
#ensembl_gene_id + description
Spurp_genes_sub<- Spurp_genes[c("ensembl_gene_id", "description", "ensembl_transcript_id")]

#pull out versions for transcripts from gast.results.sig
gast.results.sig.sep<-read.csv("Dose_exposure/RNA_data/gast.results.sig.full.csv")
gast.genes<-left_join(gast.results.sig.sep, Spurp_genes_sub, by=join_by(transcript==ensembl_transcript_id))
gast.genes <- gast.genes %>%
  distinct(transcript, .keep_all = TRUE)
write.csv(gast.genes,"Dose_exposure/RNA_data/gast.results.sig.full.genes.csv")

#pull out versions for transcripts from plut.results.sig
plut.results.sig.sep<-read.csv("Dose_exposure/RNA_data/plut.results.sig.full.csv")
plut.genes<-left_join(plut.results.sig.sep, Spurp_genes_sub, by=join_by(transcript==ensembl_transcript_id))
plut.genes <- plut.genes %>%
  distinct(transcript, .keep_all = TRUE)
write.csv(plut.genes, "Dose_exposure/RNA_data/plut.results.sig.full.genes.csv")

#compare gast and plut datasets to see overlap
shared.genes<-inner_join(gast.genes, plut.genes, by="transcript") #x is gast, y is plut
#remove unwanted columns, like duplicate columns
shared.genes<-shared.genes %>% select(-c(X.1.x, significant.x, significant.y, ensembl_gene_id.y, description.y))
write.csv(shared.genes, "Dose_exposure/RNA_data/shared.genes.csv")

#GO enrichment for DEGs####
geneID2GO <- readMappings(file="Dose_exposure/RNA_data/Spurp_GOmap.txt", sep="\t", IDsep=",")

# Read in set of all genes and make 'background' gene set
all.genes.gast <- gast.results$X #all gastrula genes
all.genes.plut <- plut.results$X #all pluteus genes

####Gastrula####
candgenes.gast <- gast.results$X[parsig.gast] #93 genes so things look good
myIG.gast <-factor(as.numeric(all.genes.gast%in%candgenes.gast))
myIG.gast <- as.numeric(as.character(myIG.gast)) #make sure it is 0s & 1s
names(myIG.gast) <- all.genes.gast
names(myIG.gast) <- sub("\\..*", "", names(myIG.gast))  # Remove everything after the first dot
topGO_selectFun <- function(x) return(x == 1)  # Select genes labeled as "1"

parname = "Gastrula"

##now separate out by category
##BP
GOdata.gast <- new("topGOdata",ontology="BP",allGenes=myIG.gast, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,
              geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata.gast,test.stat)
res <- GenTable(GOdata.gast,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("Dose_exposure/RNA_data/",parname,".GO_BP.csv",sep=""))

##MF
GOdata <- new("topGOdata",ontology="MF",allGenes=myIG.gast, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("Dose_exposure/RNA_data/",parname,".GO_MF.csv",sep=""))

##CC
GOdata <- new("topGOdata",ontology="CC",allGenes=myIG.gast, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("Dose_exposure/RNA_data/",parname,".GO_CC.csv",sep=""))

####Pluteus####
candgenes.plut <- plut.results$X[parsig.plut] #51 genes so things look good
myIG.plut <-factor(as.numeric(all.genes.plut%in%candgenes.plut))
myIG.plut <- as.numeric(as.character(myIG.plut)) #make sure it is 0s & 1s
names(myIG.plut) <- all.genes.plut
names(myIG.plut) <- sub("\\..*", "", names(myIG.plut))  # Remove everything after the first dot
topGO_selectFun <- function(x) return(x == 1)  # Select genes labeled as "1"

parname2 = "Pluteus"

##now separate out by category
##BP
GOdata.plut <- new("topGOdata",ontology="BP",allGenes=myIG.plut, nodeSize=10,
                   annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,
                   geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata.plut,test.stat)
res <- GenTable(GOdata.plut,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("Dose_exposure/RNA_data/",parname2,".GO_BP.csv",sep=""))

##MF
GOdata <- new("topGOdata",ontology="MF",allGenes=myIG.plut, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("Dose_exposure/RNA_data/",parname2,".GO_MF.csv",sep=""))

##CC
GOdata <- new("topGOdata",ontology="CC",allGenes=myIG.plut, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO,  geneSelectionFun=topGO_selectFun)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>2,]
write.csv(filt,paste("Dose_exposure/RNA_data/",parname2,".GO_CC.csv",sep=""))

