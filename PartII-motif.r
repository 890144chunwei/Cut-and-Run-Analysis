library("DESeq2")
library(ggplot2)
library(pheatmap)
library(dplyr)

CaR_SRCAP <- read.delim("~/Desktop/CaR/FC/SRCAP_ckit_combined_tss.txt", comment.char="#")
CaR_SRCAP$Location <- paste(CaR_SRCAP$Chr,CaR_SRCAP$Start,sep="_")
CaR_SRCAP_fc <- CaR_SRCAP[,7:10]
row.names(CaR_SRCAP_fc) <- CaR_SRCAP$Location

info_CaR_SRCAP <- data.frame(X = c("WT_ckit_SRCAP1","WT_ckit_SRCAP2","Mut_ckit_SRCAP1","Mut_ckit_SRCAP2"), 
                            Condition = c("WT","WT","Mut","Mut") )
dds_CaR_SRCAP <- DESeqDataSetFromMatrix(countData = CaR_SRCAP_fc, colData = info_CaR_SRCAP, design = ~Condition)

keep_CaR_SRCAP <- rowSums2(counts(dds_CaR_SRCAP)) >= 150
dds_CaR_SRCAP <- dds_CaR_SRCAP[keep_CaR_SRCAP,]
DE_CaR_SRCAP <- DESeq(dds_CaR_SRCAP)
plotMA(results(DE_CaR_SRCAP, contrast = c('Condition','Mut','WT'), alpha = 0.05), ylim =c(-5,5),main="CaR_SRCAP: Mut vs WT")

nCounts_CaR_SRCAP <- counts(DE_CaR_SRCAP ,normalized = T)

write.csv(nCounts_CaR_SRCAP,"~/Desktop/CaR/CaR_SRCAP_normCnt.txt")
Cnt_CaR_SRCAP <- read.csv("~/Desktop/CaR/CaR_SRCAP_normCnt.txt", row.names=1)

#PCA
tnormCounts_CaR_SRCAP <- t(nCounts_CaR_SRCAP)
pca_CaR_SRCAP <- prcomp(tnormCounts_CaR_SRCAP, scale. = TRUE)
screeplot(pca_CaR_SRCAP, type = "l", main = "Screeplot for ATAC_germ")
abline(h=2700, col = 'red', lty =2)

summary(pca_CaR_SRCAP)
pca_CaR_SRCAP$x
pca_CaR_SRCAP <- data.frame(info_CaR_SRCAP[,2], pca_CaR_SRCAP$x[,1:2])
colnames(pca_CaR_SRCAP) <- c('Condition', 'PC1', 'PC2')
ggplot(pca_CaR_SRCAP, aes(PC1, PC2, group=Condition)) + geom_point(size=6, aes(shape = Condition, color = Condition))+ 
  scale_shape_manual(values = c(15,16,15,16))+
  scale_color_manual(values = c("red3", "blue4"))+
  labs(title="CaR_SRCAP",x= "PC1", y= "PC2") + 
  theme(axis.text.y= element_text(size =16), axis.text.x= element_text(size =16), axis.title.x.bottom = element_text(size=20),axis.title.y = element_text(size=20))+
  xlim(-100,70) + ylim(-60, 80)

res_CaR_SRCAP <- results(DE_CaR_SRCAP , contrast = c("Condition", "Mut", "WT"), alpha = 0.05, cooksCutoff = FALSE)

res_CaR_SRCAP <- res_CaR_SRCAP[order(res_CaR_SRCAP$pvalue),]
write.csv(res_CaR_SRCAP, "~/Desktop/CaR/CaR_SRCAP_pval.txt")

Pval_CaR_SRCAP <- read.csv("~/Desktop/CaR/CaR_SRCAP_pval.txt", row.names = 1)
row.names(Pval_CaR_SRCAP) <- Pval_CaR_SRCAP$Location
row.names(CaR_SRCAP) <- CaR_SRCAP$Location

Pval_CaR_SRCAP <- na.omit(Pval_CaR_SRCAP)
Pval_CaR_SRCAP <- merge(Pval_CaR_SRCAP, CaR_SRCAP ,by=0)
Pval_CaR_SRCAP <- Pval_CaR_SRCAP[,-1]
Pval_CaR_SRCAP <- Pval_CaR_SRCAP[order(Pval_CaR_SRCAP$pvalue),]
write.csv(Pval_CaR_SRCAP, "~/Desktop/CaR/CaR_SRCAP_pval.txt")

Pval_CaR_SRCAP <- read.csv("~/Desktop/CaR/CaR_SRCAP_pval.txt", row.names = 1)

BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
BiocManager::install("org.Mm.eg.db", character.only = TRUE)
library("org.Mm.eg.db", character.only = TRUE)

CaR_SRCAP_norm <- merge(Cnt_CaR_SRCAP, CaR_SRCAP,by=0)
CaR_SRCAP_norm <- CaR_SRCAP_norm[!duplicated(CaR_SRCAP_norm$Entrez.ID),]
row.names(CaR_SRCAP_norm)<- CaR_SRCAP_norm$Entrez.ID
CaR_SRCAP_overlap_all <- merge(CaR_SRCAP_norm, Overlap,by=c('Entrez.ID','Entrez.ID'))
CaR_SRCAP_overlap_DDR <- merge(CaR_SRCAP_norm, Overlap_DDR,by=c('Entrez.ID','Entrez.ID'))
CaR_SRCAP_overlap_HSC <- merge(CaR_SRCAP_norm, Overlap_HSC, by=c('Entrez.ID','Entrez.ID'))
CaR_SRCAP_overlap_Remodel <- merge(CaR_SRCAP_norm, Overlap_Remodel, by=c('Entrez.ID','Entrez.ID'))
write.csv(CaR_SRCAP_overlap_all, "~/Desktop/CaR/CaR_SRCAP_overlap_all.csv")
write.csv(CaR_SRCAP_overlap_DDR, "~/Desktop/CaR/CaR_SRCAP_overlap_DDR.csv")
write.csv(CaR_SRCAP_overlap_HSC, "~/Desktop/CaR/CaR_SRCAP_overlap_HSC.csv")
write.csv(CaR_SRCAP_overlap_Remodel, "~/Desktop/CaR/CaR_SRCAP_overlap_Remodel.csv")

BiocManager::install("ChIPseeker")
library("ChIPseeker")
library(clusterProfiler)
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

promoter <- getPromoters(TxDb=TxDb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(ATAC_germ_wt_peak, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")

DNArepair_genelist <- read.delim2("~/Dropbox/Mac_Desktop/DNArepair_genelist_full.txt")
DNArepair_genelist<-DNArepair_genelist[,2:3]
DNArepair_genelist<-unique(DNArepair_genelist)
row.names(DNArepair_genelist)<-DNArepair_genelist$converted_alias
HSC.genelist <- read.csv("~/Dropbox/Mac_Desktop/HSC gene list_full.csv")
HSC.genelist<-HSC.genelist[,2:3]
HSC.genelist<-unique(HSC.genelist)
row.names(HSC.genelist)<-HSC.genelist$converted_alias
Remodeling_genelist <- read.delim2("~/Dropbox/Mac_Desktop/Remodeling_genelist_full.txt")
Remodeling_genelist<-Remodeling_genelist[,2:3]
Remodeling_genelist<-unique(Remodeling_genelist)
row.names(Remodeling_genelist)<-Remodeling_genelist$converted_alias
write.csv(Overlap_DDR, "~/Desktop/ATAC-seq/Overlap_DDR.csv",row.names = FALSE, col.names = FALSE)
write.csv(Overlap_HSC, "~/Desktop/ATAC-seq/Overlap_HSC.csv",row.names = FALSE, col.names = FALSE)
write.csv(Overlap_Remodel, "~/Desktop/ATAC-seq/Overlap_Remodel.csv",row.names = FALSE, col.names = FALSE)

BiocManager::install("memes")
library("memes")
BiocManager::install("plyranges")
library("plyranges")
BiocManager::install("GenomicRanges")
library("GenomicRanges")
BiocManager::install("IRanges")
library("IRanges")
install.packages("magrittr")
library("magrittr")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene

Pval_CaR_H2AZ<-Pval_CaR_H2AZ[order(Pval_CaR_H2AZ$baseMean, decreasing =T),]
Pval_CaR_H2AZ_f <- Pval_CaR_H2AZ[1:200,]
Pval_CaR_SRCAP<-Pval_CaR_SRCAP[order(Pval_CaR_SRCAP$baseMean, decreasing =T),]
Pval_CaR_SRCAP_f <- Pval_CaR_SRCAP[1:200,]
Combined_peaks <- rbind(Pval_CaR_H2AZ_f[,8:10],Pval_CaR_SRCAP_f[,8:10])
ms.genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
Combined_peaks_g <- makeGRangesFromDataFrame(Combined_peaks, keep.extra.columns = TRUE)
Combined_peaks_s <- resize(Combined_peaks_g,100,"center")
Combined_peaks_s <- get_sequence(Combined_peaks_s, ms.genome)
Combined_peaks_dreme <- runDreme(Combined_peaks_s, control = "shuffle", seed =3, e=0.001,nmotifs = 3)
library(universalmotif)
combined_peaks_list <- to_list(Combined_peaks_dreme)
view_motifs(combined_peaks_list)
plot(combined_peaks_list[1])
