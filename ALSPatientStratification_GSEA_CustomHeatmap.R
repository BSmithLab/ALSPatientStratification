##Custom Gene Set Enrichment Heatmaps
#Manual Heatmaps

#Written By: Jarrett Eshima
#For: Dr. Barbara Smith Lab
#Date: September 27th 2021

#Note: this code is repetitive, each pathway considered is separated by hashes

##############################################################################################################################################################
# Heatmap from GSEA Rank Metric Score
##############################################################################################################################################################
library(RColorBrewer)
library(gplots)
##############################################################################################################################################################

#Genes downregulated in Alzheimer's

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/GSEA/FINAL/ManualHeatmaps/AlzDown")

TE = read.csv("TE_Alz_Heatmap.csv") #File generated from GSEA 
OX = read.csv("OX_Alz_Heatmap.csv") #File generated from GSEA 
GLIA = read.csv("GLIA_Alz_Heatmap.csv") #File generated from GSEA 
HC = read.csv("HC_Alz_Heatmap.csv") #File generated from GSEA 

Alz_heatmapD = data.frame(matrix(NA,nrow=nrow(TE),ncol=4))
rownames(Alz_heatmapD) = TE$SYMBOL
colnames(Alz_heatmapD) = c("HC","GLIA","TE","OX")

TEscore = as.numeric(unlist(TE$TE.Score))
OXscore = as.numeric(unlist(OX$OX.Score))
GLIAscore = as.numeric(unlist(GLIA$GLIA.Score))
HCscore = as.numeric(unlist(HC$HC.Score))

Alz_heatmapD$TE = TEscore
Alz_heatmapD$OX = OXscore
Alz_heatmapD$GLIA = GLIAscore
Alz_heatmapD$HC = HCscore


Alz_heatmapD2 = matrix(as.numeric(unlist(Alz_heatmapD),ncol = 4))
tmp = matrix(Alz_heatmapD2,ncol=4)
typeof(tmp)
rownames(tmp) = TE$SYMBOL
colnames(tmp)= c("Control","ALS-GLIA","ALS-TE","ALS-OX")


mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))
heatmap.2(tmp,scale="none",col=mypalette,trace="none",density.info="none",key.xlab = "GSEA Rank Metric Score",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "ALS Subtype",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 0.7,cexCol = 1.25,Colv = F)


##############################################################################################################################################################

#Genes unregulated in Alzheimer's

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/GSEA/FINAL/ManualHeatmaps/AlzUp")

TE = read.csv("TE_Alz_Heatmap.csv") #File generated from GSEA 
OX = read.csv("OX_Alz_Heatmap.csv") #File generated from GSEA 
GLIA = read.csv("GLIA_Alz_Heatmap.csv") #File generated from GSEA 
HC = read.csv("HC_Alz_Heatmap.csv") #File generated from GSEA 

Alz_heatmapU = data.frame(matrix(NA,nrow=nrow(TE),ncol=4))
rownames(Alz_heatmapU) = TE$SYMBOL
colnames(Alz_heatmapU) = c("HC","GLIA","TE","OX")

TEscore = as.numeric(unlist(TE$TE.Score))
OXscore = as.numeric(unlist(OX$OX.Score))
GLIAscore = as.numeric(unlist(GLIA$GLIA.Score))
HCscore = as.numeric(unlist(HC$HC.Score))

Alz_heatmapU$TE = TEscore
Alz_heatmapU$OX = OXscore
Alz_heatmapU$GLIA = GLIAscore
Alz_heatmapU$HC = HCscore


Alz_heatmapU2 = matrix(as.numeric(unlist(Alz_heatmapU),ncol = 4))
tmp = matrix(Alz_heatmapU2,ncol=4)
typeof(tmp)
rownames(tmp) = TE$SYMBOL
colnames(tmp)= c("Control","ALS-GLIA","ALS-TE","ALS-OX")


mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))
heatmap.2(tmp,scale="none",col=mypalette,trace="none",density.info="none",key.xlab = "GSEA Rank Metric Score",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "ALS Subtype",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 0.55,cexCol = 1.25,Colv = F)

##############################################################################################################################################################

#GPCR Ligand Binding

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/GSEA/FINAL/ManualHeatmaps/GPCR")

TE = read.csv("TE_GPCR_Heatmap.csv") #File generated from GSEA 
OX = read.csv("OX_GPCR_Heatmap.csv") #File generated from GSEA 
GLIA = read.csv("GLIA_GPCR_Heatmap.csv") #File generated from GSEA 
HC = read.csv("HC_GPCR_Heatmap.csv") #File generated from GSEA 

GPCR_heatmap = data.frame(matrix(NA,nrow=nrow(TE),ncol=4))
rownames(GPCR_heatmap) = TE$SYMBOL
colnames(GPCR_heatmap) = c("HC","GLIA","TE","OX")

TEscore = as.numeric(unlist(TE$TE.Score))
OXscore = as.numeric(unlist(OX$OX.Score))
GLIAscore = as.numeric(unlist(GLIA$GLIA.Score))
HCscore = as.numeric(unlist(HC$HC.Score))

GPCR_heatmap$TE = TEscore
GPCR_heatmap$OX = OXscore
GPCR_heatmap$GLIA = GLIAscore
GPCR_heatmap$HC = HCscore


GPCR_heatmap2 = matrix(as.numeric(unlist(GPCR_heatmap),ncol = 4))
tmp = matrix(GPCR_heatmap2,ncol=4)
typeof(tmp)
rownames(tmp) = TE$SYMBOL
colnames(tmp)= c("Control","ALS-GLIA","ALS-TE","ALS-OX")


mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))
heatmap.2(tmp,scale="none",col=mypalette,trace="none",density.info="none",key.xlab = "GSEA Rank Metric Score",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "ALS Subtype",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 0.5,cexCol = 1.25,Colv = F)


##############################################################################################################################################################

#Adaptive Immune System

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/GSEA/FINAL/ManualHeatmaps/Adaptive Immune")

TE = read.csv("TE_Adaptive_Heatmap.csv") #File generated from GSEA 
OX = read.csv("OX_Adaptive_Heatmap.csv") #File generated from GSEA 
GLIA = read.csv("GLIA_Adaptive_Heatmap.csv") #File generated from GSEA 
HC = read.csv("HC_Adaptive_Heatmap.csv") #File generated from GSEA 

Immune_heatmap = data.frame(matrix(NA,nrow=nrow(TE),ncol=4))
rownames(Immune_heatmap) = TE$SYMBOL
colnames(Immune_heatmap) = c("HC","GLIA","TE","OX")

TEscore = as.numeric(unlist(TE$TE.Score))
OXscore = as.numeric(unlist(OX$OX.Score))
GLIAscore = as.numeric(unlist(GLIA$GLIA.Score))
HCscore = as.numeric(unlist(HC$HC.Score))

Immune_heatmap$TE = TEscore
Immune_heatmap$OX = OXscore
Immune_heatmap$GLIA = GLIAscore
Immune_heatmap$HC = HCscore


Immune_heatmap2 = matrix(as.numeric(unlist(Immune_heatmap),ncol = 4))
tmp = matrix(Immune_heatmap2,ncol=4)
typeof(tmp)
rownames(tmp) = TE$SYMBOL
colnames(tmp)= c("Control","ALS-GLIA","ALS-TE","ALS-OX")


mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))
heatmap.2(tmp,scale="none",col=mypalette,trace="none",density.info="none",key.xlab = "GSEA Rank Metric Score",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "ALS Subtype",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 0.65,cexCol = 1.25,Colv = F)

##############################################################################################################################################################

#Immunoregulatory Interactions

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/GSEA/FINAL/ManualHeatmaps/Immunoregulatory")

TE = read.csv("TE_IReg_Heatmap.csv") #File generated from GSEA 
OX = read.csv("OX_IReg_Heatmap.csv") #File generated from GSEA 
GLIA = read.csv("GLIA_IReg_Heatmap.csv") #File generated from GSEA 
HC = read.csv("HC_IReg_Heatmap.csv") #File generated from GSEA 

ImmuneR_heatmap = data.frame(matrix(NA,nrow=nrow(TE),ncol=4))
rownames(ImmuneR_heatmap) = TE$SYMBOL
colnames(ImmuneR_heatmap) = c("HC","GLIA","TE","OX")

TEscore = as.numeric(unlist(TE$TE.Score))
OXscore = as.numeric(unlist(OX$OX.Score))
GLIAscore = as.numeric(unlist(GLIA$GLIA.Score))
HCscore = as.numeric(unlist(HC$HC.Score))

ImmuneR_heatmap$TE = TEscore
ImmuneR_heatmap$OX = OXscore
ImmuneR_heatmap$GLIA = GLIAscore
ImmuneR_heatmap$HC = HCscore


ImmuneR_heatmap2 = matrix(as.numeric(unlist(ImmuneR_heatmap),ncol = 4))
tmp = matrix(ImmuneR_heatmap2,ncol=4)
typeof(tmp)
rownames(tmp) = TE$SYMBOL
colnames(tmp)= c("Control","ALS-GLIA","ALS-TE","ALS-OX")


mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))
heatmap.2(tmp,scale="none",col=mypalette,trace="none",density.info="none",key.xlab = "GSEA Rank Metric Score",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "ALS Subtype",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 1,cexCol = 1.25,Colv = F)


##############################################################################################################################################################

#Diseases of Glycosylation

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/GSEA/FINAL/ManualHeatmaps/Glycosylation")

TE = read.csv("TE_Glyco_Heatmap.csv") #File generated from GSEA 
OX = read.csv("OX_Glyco_Heatmap.csv") #File generated from GSEA 
GLIA = read.csv("GLIA_Glyco_Heatmap.csv") #File generated from GSEA 
HC = read.csv("HC_Glyco_Heatmap.csv") #File generated from GSEA 

Glyco_heatmap = data.frame(matrix(NA,nrow=nrow(TE),ncol=4))
rownames(Glyco_heatmap) = TE$SYMBOL
colnames(Glyco_heatmap) = c("HC","GLIA","TE","OX")

TEscore = as.numeric(unlist(TE$TE.Score))
OXscore = as.numeric(unlist(OX$OX.Score))
GLIAscore = as.numeric(unlist(GLIA$GLIA.Score))
HCscore = as.numeric(unlist(HC$HC.Score))

Glyco_heatmap$TE = TEscore
Glyco_heatmap$OX = OXscore
Glyco_heatmap$GLIA = GLIAscore
Glyco_heatmap$HC = HCscore


Glyco_heatmap2 = matrix(as.numeric(unlist(Glyco_heatmap),ncol = 4))
tmp = matrix(Glyco_heatmap2,ncol=4)
typeof(tmp)
rownames(tmp) = TE$SYMBOL
colnames(tmp)= c("Control","ALS-GLIA","ALS-TE","ALS-OX")


mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))
heatmap.2(tmp,scale="none",col=mypalette,trace="none",density.info="none",key.xlab = "GSEA Rank Metric Score",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "ALS Subtype",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 1.1,cexCol = 1.25,Colv = F)


##############################################################################################################################################################

#PTM

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/GSEA/FINAL/ManualHeatmaps/PTM")

TE = read.csv("TE_PTM_Heatmap.csv") #File generated from GSEA 
OX = read.csv("OX_PTM_Heatmap.csv") #File generated from GSEA 
GLIA = read.csv("GLIA_PTM_Heatmap.csv") #File generated from GSEA 
HC = read.csv("HC_PTM_Heatmap.csv") #File generated from GSEA 

PTM_heatmap = data.frame(matrix(NA,nrow=nrow(TE),ncol=4))
rownames(PTM_heatmap) = TE$SYMBOL
colnames(PTM_heatmap) = c("HC","GLIA","TE","OX")

TEscore = as.numeric(unlist(TE$TE.Score))
OXscore = as.numeric(unlist(OX$OX.Score))
GLIAscore = as.numeric(unlist(GLIA$GLIA.Score))
HCscore = as.numeric(unlist(HC$HC.Score))

PTM_heatmap$TE = TEscore
PTM_heatmap$OX = OXscore
PTM_heatmap$GLIA = GLIAscore
PTM_heatmap$HC = HCscore


PTM_heatmap2 = matrix(as.numeric(unlist(PTM_heatmap),ncol = 4))
tmp = matrix(PTM_heatmap2,ncol=4)
typeof(tmp)
rownames(tmp) = TE$SYMBOL
colnames(tmp)= c("Control","ALS-GLIA","ALS-TE","ALS-OX")


mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))
heatmap.2(tmp,scale="none",col=mypalette,trace="none",density.info="none",key.xlab = "GSEA Rank Metric Score",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "ALS Subtype",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 0.7,cexCol = 1.25,Colv = F)


##############################################################################################################################################################

#ECM

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/GSEA/FINAL/ManualHeatmaps/ECM")

TE = read.csv("TE_ECM_Heatmap.csv") #File generated from GSEA 
OX = read.csv("OX_ECM_Heatmap.csv") #File generated from GSEA 
GLIA = read.csv("GLIA_ECM_Heatmap.csv") #File generated from GSEA 
HC = read.csv("HC_ECM_Heatmap.csv") #File generated from GSEA 

ECM_heatmap = data.frame(matrix(NA,nrow=nrow(TE),ncol=4))
rownames(ECM_heatmap) = TE$SYMBOL
colnames(ECM_heatmap) = c("HC","GLIA","TE","OX")

TEscore = as.numeric(unlist(TE$TE.Score))
OXscore = as.numeric(unlist(OX$OX.Score))
GLIAscore = as.numeric(unlist(GLIA$GLIA.Score))
HCscore = as.numeric(unlist(HC$HC.Score))

ECM_heatmap$TE = TEscore
ECM_heatmap$OX = OXscore
ECM_heatmap$GLIA = GLIAscore
ECM_heatmap$HC = HCscore


ECM_heatmap2 = matrix(as.numeric(unlist(ECM_heatmap),ncol = 4))
tmp = matrix(ECM_heatmap2,ncol=4)
typeof(tmp)
rownames(tmp) = TE$SYMBOL
colnames(tmp)= c("Control","ALS-GLIA","ALS-TE","ALS-OX")


mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))
heatmap.2(tmp,scale="none",col=mypalette,trace="none",density.info="none",key.xlab = "GSEA Rank Metric Score",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "ALS Subtype",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 0.9,cexCol = 1.25,Colv = F)


##############################################################################################################################################################

#Tyrosine Kinase Signaling

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/GSEA/FINAL/ManualHeatmaps/TyrosineKinase")

TE = read.csv("TE_TK_Heatmap.csv") #File generated from GSEA 
OX = read.csv("OX_TK_Heatmap.csv") #File generated from GSEA 
GLIA = read.csv("GLIA_TK_Heatmap.csv") #File generated from GSEA 
HC = read.csv("HC_TK_Heatmap.csv") #File generated from GSEA 

TK_heatmap = data.frame(matrix(NA,nrow=nrow(TE),ncol=4))
rownames(TK_heatmap) = TE$SYMBOL
colnames(TK_heatmap) = c("HC","GLIA","TE","OX")

TEscore = as.numeric(unlist(TE$TE.Score))
OXscore = as.numeric(unlist(OX$OX.Score))
GLIAscore = as.numeric(unlist(GLIA$GLIA.Score))
HCscore = as.numeric(unlist(HC$HC.Score))

TK_heatmap$TE = TEscore
TK_heatmap$OX = OXscore
TK_heatmap$GLIA = GLIAscore
TK_heatmap$HC = HCscore


TK_heatmap2 = matrix(as.numeric(unlist(TK_heatmap),ncol = 4))
tmp = matrix(TK_heatmap2,ncol=4)
typeof(tmp)
rownames(tmp) = TE$SYMBOL
colnames(tmp)= c("Control","ALS-GLIA","ALS-TE","ALS-OX")


mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))
heatmap.2(tmp,scale="none",col=mypalette,trace="none",density.info="none",key.xlab = "GSEA Rank Metric Score",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "ALS Subtype",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 0.9,cexCol = 1.25,Colv = F)

##############################################################################################################################################################

#GTPase 

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/GSEA/FINAL/ManualHeatmaps/GTPase")

TE = read.csv("TE_GTPase_Heatmap.csv") #File generated from GSEA 
OX = read.csv("OX_GTPase_Heatmap.csv") #File generated from GSEA 
GLIA = read.csv("GLIA_GTPase_Heatmap.csv") #File generated from GSEA 
HC = read.csv("HC_GTPase_Heatmap.csv") #File generated from GSEA 

GTP_heatmap = data.frame(matrix(NA,nrow=nrow(TE),ncol=4))
rownames(GTP_heatmap) = TE$SYMBOL
colnames(GTP_heatmap) = c("HC","GLIA","TE","OX")

TEscore = as.numeric(unlist(TE$TE.Score))
OXscore = as.numeric(unlist(OX$OX.Score))
GLIAscore = as.numeric(unlist(GLIA$GLIA.Score))
HCscore = as.numeric(unlist(HC$HC.Score))

GTP_heatmap$TE = TEscore
GTP_heatmap$OX = OXscore
GTP_heatmap$GLIA = GLIAscore
GTP_heatmap$HC = HCscore


GTP_heatmap2 = matrix(as.numeric(unlist(GTP_heatmap),ncol = 4))
tmp = matrix(GTP_heatmap2,ncol=4)
typeof(tmp)
rownames(tmp) = TE$SYMBOL
colnames(tmp)= c("Control","ALS-GLIA","ALS-TE","ALS-OX")


mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))
heatmap.2(tmp,scale="none",col=mypalette,trace="none",density.info="none",key.xlab = "GSEA Rank Metric Score",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "ALS Subtype",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 1.1,cexCol = 1.25,Colv = F)

##############################################################################################################################################################

#Synaptic Transmission

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/GSEA/FINAL/ManualHeatmaps/SynapticTransmission")

TE = read.csv("TE_Syn_Heatmap.csv") #File generated from GSEA 
OX = read.csv("OX_Syn_Heatmap.csv") #File generated from GSEA 
GLIA = read.csv("GLIA_Syn_Heatmap.csv") #File generated from GSEA 
HC = read.csv("HC_Syn_Heatmap.csv") #File generated from GSEA 

Syn_heatmap = data.frame(matrix(NA,nrow=nrow(TE),ncol=4))
rownames(Syn_heatmap) = TE$SYMBOL
colnames(Syn_heatmap) = c("HC","GLIA","TE","OX")

TEscore = as.numeric(unlist(TE$TE.Score))
OXscore = as.numeric(unlist(OX$OX.Score))
GLIAscore = as.numeric(unlist(GLIA$GLIA.Score))
HCscore = as.numeric(unlist(HC$HC.Score))

Syn_heatmap$TE = TEscore
Syn_heatmap$OX = OXscore
Syn_heatmap$GLIA = GLIAscore
Syn_heatmap$HC = HCscore


Syn_heatmap2 = matrix(as.numeric(unlist(Syn_heatmap),ncol = 4))
tmp = matrix(Syn_heatmap2,ncol=4)
typeof(tmp)
rownames(tmp) = TE$SYMBOL
colnames(tmp)= c("Control","ALS-GLIA","ALS-TE","ALS-OX")


mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))
heatmap.2(tmp,scale="none",col=mypalette,trace="none",density.info="none",key.xlab = "GSEA Rank Metric Score",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "ALS Subtype",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 1.2,cexCol = 1.25,Colv = F)


##############################################################################################################################################################

#Pol II

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/GSEA/FINAL/ManualHeatmaps/Pol2")

TE = read.csv("TE_Pol_Heatmap.csv") #File generated from GSEA 
OX = read.csv("OX_Pol_Heatmap.csv") #File generated from GSEA 
GLIA = read.csv("GLIA_Pol_Heatmap.csv") #File generated from GSEA 
HC = read.csv("HC_Pol_Heatmap.csv") #File generated from GSEA 

Pol2_heatmap = data.frame(matrix(NA,nrow=nrow(TE),ncol=4))
rownames(Pol2_heatmap) = TE$SYMBOL
colnames(Pol2_heatmap) = c("HC","GLIA","TE","OX")

TEscore = as.numeric(unlist(TE$TE.Score))
OXscore = as.numeric(unlist(OX$OX.Score))
GLIAscore = as.numeric(unlist(GLIA$GLIA.Score))
HCscore = as.numeric(unlist(HC$HC.Score))

Pol2_heatmap$TE = TEscore
Pol2_heatmap$OX = OXscore
Pol2_heatmap$GLIA = GLIAscore
Pol2_heatmap$HC = HCscore


Pol2_heatmap2 = matrix(as.numeric(unlist(Pol2_heatmap),ncol = 4))
tmp = matrix(Pol2_heatmap2,ncol=4)
typeof(tmp)
rownames(tmp) = TE$SYMBOL
colnames(tmp)= c("Control","ALS-GLIA","ALS-TE","ALS-OX")


mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))
heatmap.2(tmp,scale="none",col=mypalette,trace="none",density.info="none",key.xlab = "GSEA Rank Metric Score",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "ALS Subtype",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 1.1,cexCol = 1.25,Colv = F)

##############################################################################################################################################################

#Small Molecule Transport

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/GSEA/FINAL/ManualHeatmaps/Small Molecule Transport")

TE = read.csv("TE_Small_Heatmap.csv") #File generated from GSEA 
OX = read.csv("OX_Small_Heatmap.csv") #File generated from GSEA 
GLIA = read.csv("GLIA_Small_Heatmap.csv") #File generated from GSEA 
HC = read.csv("HC_Small_Heatmap.csv") #File generated from GSEA 

Small_heatmap = data.frame(matrix(NA,nrow=nrow(TE),ncol=4))
rownames(Small_heatmap) = TE$SYMBOL
colnames(Small_heatmap) = c("HC","GLIA","TE","OX")

TEscore = as.numeric(unlist(TE$TE.Score))
OXscore = as.numeric(unlist(OX$OX.Score))
GLIAscore = as.numeric(unlist(GLIA$GLIA.Score))
HCscore = as.numeric(unlist(HC$HC.Score))

Small_heatmap$TE = TEscore
Small_heatmap$OX = OXscore
Small_heatmap$GLIA = GLIAscore
Small_heatmap$HC = HCscore


Small_heatmap2 = matrix(as.numeric(unlist(Small_heatmap),ncol = 4))
tmp = matrix(Small_heatmap2,ncol=4)
typeof(tmp)
rownames(tmp) = TE$SYMBOL
colnames(tmp)= c("Control","ALS-GLIA","ALS-TE","ALS-OX")


mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))
heatmap.2(tmp,scale="none",col=mypalette,trace="none",density.info="none",key.xlab = "GSEA Rank Metric Score",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "ALS Subtype",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 0.55,cexCol = 1.25,Colv = F)


##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
##Necessary function

#Convert to VST Z-score
zscore = function(x){
  Z = rep(NA,length(x))
  for(i in 1:length(x)){
    mu = mean(x)
    sigma = sd(x)
    Z[i] = (x[i]-mu)/sigma
  }
  return(Z)
}

##############################################################################################################################################################
# Heatmap from VST transformed Counts (z-score normalized)
##############################################################################################################################################################

#Pathway feature lists
AlzDown_featurelist = rownames(Alz_heatmapD)
AlzUp_featurelist = rownames(Alz_heatmapU)
TK_featurelist = rownames(TK_heatmap)
ECM_featurelist = rownames(ECM_heatmap)
PTM_featurelist = rownames(PTM_heatmap)
Immune_featurelist = rownames(Immune_heatmap)
Glyco_featurelist = rownames(Glyco_heatmap)
Pol2_featurelist = rownames(Pol2_heatmap)
Small_featurelist = rownames(Small_heatmap)
IReg_featurelist = rownames(ImmuneR_heatmap)
GPCR_featurelist = rownames(GPCR_heatmap)
GTPase_featurelist = rownames(GTP_heatmap)
Synapse_featurelist = rownames(Syn_heatmap)

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/GSEA")
Counts = read.csv("ALS_HC_ClassifierMatrix.csv")
rownames(Counts) = Counts[,1]
Counts = Counts[,-1]

#TE feature list
hist(nchar(rownames(Counts)))
TEindex = which(nchar(rownames(Counts)) > 25)
TE_featurelist = rownames(Counts)[TEindex]

#Reorder the columns so that all Control samples are first, then Glia, then TE, then OX
subtypes = read.csv("ALS451_coldata_SUBTYPES.csv")

TEcounts = Gliacounts = Oxcounts = HCcounts = matrix(NA,nrow = nrow(Counts),ncol = ncol(Counts))
colnames(TEcounts) = colnames(Gliacounts) = colnames(Oxcounts) = colnames(HCcounts) = colnames(Counts)
rownames(TEcounts) = rownames(Gliacounts) = rownames(Oxcounts) = rownames(HCcounts) = rownames(Counts)

#Parse count data by ALS subtype
for(i in 1:nrow(subtypes)){
  
  if(subtypes$Subtype[i] == "TE"){
    
    tag = subtypes$Subject[i]
    tmpindex = which(colnames(Counts) == tag)
    tmpcounts = Counts[,tmpindex]
    
    
    for(j in 1:ncol(TEcounts)){
      if(colnames(TEcounts)[j] == tag){
        TEcounts[,j] = tmpcounts
      }
    }
    
  }
  
  if(subtypes$Subtype[i] == "OX"){
    
    tag = subtypes$Subject[i]
    tmpindex = which(colnames(Counts) == tag)
    tmpcounts = Counts[,tmpindex]
    
    
    for(j in 1:ncol(Oxcounts)){
      if(colnames(Oxcounts)[j] == tag){
        Oxcounts[,j] = tmpcounts
      }
    }
    
  }
  
  if(subtypes$Subtype[i] == "GLIA"){
    
    tag = subtypes$Subject[i]
    tmpindex = which(colnames(Counts) == tag)
    tmpcounts = Counts[,tmpindex]
    
    
    for(j in 1:ncol(Gliacounts)){
      if(colnames(Gliacounts)[j] == tag){
        Gliacounts[,j] = tmpcounts
      }
    }
    
  }
  
}

#Clean up
table(subtypes$Subtype)

empty1 = which(is.na(colSums(TEcounts)))
empty2 = which(is.na(colSums(Oxcounts)))
empty3 = which(is.na(colSums(Gliacounts)))

TEcounts = TEcounts[,-empty1] #correct dimension
Oxcounts = Oxcounts[,-empty2] #correct dimension
Gliacounts = Gliacounts[,-empty3] #correct dimension

ALSnames = c(colnames(TEcounts),colnames(Oxcounts),colnames(Gliacounts))
HCcounts = Counts[,! colnames(Counts) %in% ALSnames] #correct dimension

#Reordered
FinalCounts = cbind(HCcounts,Gliacounts,TEcounts,Oxcounts)

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Pub Data/10-4-21")

##############################################################################################################################################################

#Transposable Elements

TE_heatmap_counts = FinalCounts[rownames(FinalCounts) %in% TE_featurelist,]
tmp = as.numeric(unlist(TE_heatmap_counts))
TE_heatmap_counts_num = matrix(tmp,nrow = nrow(TE_heatmap_counts))

colnames(TE_heatmap_counts_num) = colnames(TE_heatmap_counts)
rownames(TE_heatmap_counts_num) = rownames(TE_heatmap_counts)

TE_heatmap_zscore = apply(TE_heatmap_counts_num,1,zscore)

#transpose
TE_heatmap_zscore = t(TE_heatmap_zscore)
colnames(TE_heatmap_zscore) = colnames(TE_heatmap_counts_num)

mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))

pdf("GSEA_TE_Heatmap_Subtypes.pdf",width = 24, height = 18)
heatmap.2(TE_heatmap_zscore,scale="none",Colv = F,col=mypalette,trace="none",density.info="none",key.xlab = "Z-Score Counts",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "Sample",ylab="Transposable Element",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 0.2,cexCol = 0.15)
dev.off()

#Get the TE dendrogram rowname list for figure
tmp = heatmap.2(TE_heatmap_zscore,scale="none",Colv = F,col=mypalette,trace="none",density.info="none",key.xlab = "Z-Score Counts",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "Sample",ylab="Transposable Element",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 0.2,cexCol = 0.15)
tmp2 = tmp$rowInd
TErn = TE_featurelist[tmp2]
TErn_top2bottom = rev(TErn)
write.csv(TErn_top2bottom,"GSEA_TE_Heatmap_Subtypes_Rownames.csv")
##############################################################################################################################################################

#ECM

ECM_heatmap_counts = FinalCounts[rownames(FinalCounts) %in% ECM_featurelist,]
tmp = as.numeric(unlist(ECM_heatmap_counts))
ECM_heatmap_counts_num = matrix(tmp,nrow = nrow(ECM_heatmap_counts))

colnames(ECM_heatmap_counts_num) = colnames(ECM_heatmap_counts)
rownames(ECM_heatmap_counts_num) = rownames(ECM_heatmap_counts)

ECM_heatmap_zscore = apply(ECM_heatmap_counts_num,1,zscore)

#transpose
ECM_heatmap_zscore = t(ECM_heatmap_zscore)
colnames(ECM_heatmap_zscore) = colnames(ECM_heatmap_counts_num)

mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))

pdf("GSEA_ECM_Heatmap_Subtypes.pdf",width = 24, height = 18)
heatmap.2(ECM_heatmap_zscore,scale="none",Colv = F,col=mypalette,trace="none",density.info="none",key.xlab = "Z-Score Counts",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "Sample",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 1.1,cexCol = 0.15)
dev.off()
##############################################################################################################################################################

#PTM

PTM_heatmap_counts = FinalCounts[rownames(FinalCounts) %in% PTM_featurelist,]
tmp = as.numeric(unlist(PTM_heatmap_counts))
PTM_heatmap_counts_num = matrix(tmp,nrow = nrow(PTM_heatmap_counts))

colnames(PTM_heatmap_counts_num) = colnames(PTM_heatmap_counts)
rownames(PTM_heatmap_counts_num) = rownames(PTM_heatmap_counts)

PTM_heatmap_zscore = apply(PTM_heatmap_counts_num,1,zscore)

#transpose
PTM_heatmap_zscore = t(PTM_heatmap_zscore)
colnames(PTM_heatmap_zscore) = colnames(PTM_heatmap_counts_num)

mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))

pdf("GSEA_PTM_Heatmap_Subtypes.pdf",width = 24, height = 18)
heatmap.2(PTM_heatmap_zscore,scale="none",Colv = F,col=mypalette,trace="none",density.info="none",key.xlab = "Z-Score Counts",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "Sample",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 1.1,cexCol = 0.15)
dev.off()
##############################################################################################################################################################

#Tyrosine Kinase

TK_heatmap_counts = FinalCounts[rownames(FinalCounts) %in% TK_featurelist,]
tmp = as.numeric(unlist(TK_heatmap_counts))
TK_heatmap_counts_num = matrix(tmp,nrow = nrow(TK_heatmap_counts))

colnames(TK_heatmap_counts_num) = colnames(TK_heatmap_counts)
rownames(TK_heatmap_counts_num) = rownames(TK_heatmap_counts)

TK_heatmap_zscore = apply(TK_heatmap_counts_num,1,zscore)

#transpose
TK_heatmap_zscore = t(TK_heatmap_zscore)
colnames(TK_heatmap_zscore) = colnames(TK_heatmap_counts_num)

mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))

pdf("GSEA_TK_Heatmap_Subtypes.pdf",width = 24, height = 18)
heatmap.2(TK_heatmap_zscore,scale="none",Colv = F,col=mypalette,trace="none",density.info="none",key.xlab = "Z-Score Counts",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "Sample",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 1.1,cexCol = 0.15)
dev.off()
##############################################################################################################################################################

#Adaptive Immune System

Immune_heatmap_counts = FinalCounts[rownames(FinalCounts) %in% Immune_featurelist,]
tmp = as.numeric(unlist(Immune_heatmap_counts))
Immune_heatmap_counts_num = matrix(tmp,nrow = nrow(Immune_heatmap_counts))

colnames(Immune_heatmap_counts_num) = colnames(Immune_heatmap_counts)
rownames(Immune_heatmap_counts_num) = rownames(Immune_heatmap_counts)

Immune_heatmap_zscore = apply(Immune_heatmap_counts_num,1,zscore)

#transpose
Immune_heatmap_zscore = t(Immune_heatmap_zscore)
colnames(Immune_heatmap_zscore) = colnames(Immune_heatmap_counts_num)

mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))

pdf("GSEA_AdaptiveImmune_Heatmap_Subtypes.pdf",width = 24, height = 18)
heatmap.2(Immune_heatmap_zscore,scale="none",Colv = F,col=mypalette,trace="none",density.info="none",key.xlab = "Z-Score Counts",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "Sample",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 0.7,cexCol = 0.15)
dev.off()
##############################################################################################################################################################

#Downregulated Alzheimer's Genes

AD_heatmap_counts = FinalCounts[rownames(FinalCounts) %in% AlzDown_featurelist,]
tmp = as.numeric(unlist(AD_heatmap_counts))
AD_heatmap_counts_num = matrix(tmp,nrow = nrow(AD_heatmap_counts))

colnames(AD_heatmap_counts_num) = colnames(AD_heatmap_counts)
rownames(AD_heatmap_counts_num) = rownames(AD_heatmap_counts)

AD_heatmap_zscore = apply(AD_heatmap_counts_num,1,zscore)

#transpose
AD_heatmap_zscore = t(AD_heatmap_zscore)
colnames(AD_heatmap_zscore) = colnames(AD_heatmap_counts_num)

mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))

pdf("GSEA_AlzDown_Heatmap_Subtypes.pdf",width = 24, height = 18)
heatmap.2(AD_heatmap_zscore,scale="none",Colv = F,col=mypalette,trace="none",density.info="none",key.xlab = "Z-Score Counts",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "Sample",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 0.7,cexCol = 0.15)
dev.off()
##############################################################################################################################################################

#Upregulated Alzheimer's Genes

AU_heatmap_counts = FinalCounts[rownames(FinalCounts) %in% AlzUp_featurelist,]
tmp = as.numeric(unlist(AU_heatmap_counts))
AU_heatmap_counts_num = matrix(tmp,nrow = nrow(AU_heatmap_counts))

colnames(AU_heatmap_counts_num) = colnames(AU_heatmap_counts)
rownames(AU_heatmap_counts_num) = rownames(AU_heatmap_counts)

AU_heatmap_zscore = apply(AU_heatmap_counts_num,1,zscore)

#transpose
AU_heatmap_zscore = t(AU_heatmap_zscore)
colnames(AU_heatmap_zscore) = colnames(AU_heatmap_counts_num)

mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))

pdf("GSEA_AlzUp_Heatmap_Subtypes.pdf",width = 24, height = 18)
heatmap.2(AU_heatmap_zscore,scale="none",Colv = F,col=mypalette,trace="none",density.info="none",key.xlab = "Z-Score Counts",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "Sample",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 0.6,cexCol = 0.15)
dev.off()
##############################################################################################################################################################

#Diseases of Glycosylation Genes

Glyco_heatmap_counts = FinalCounts[rownames(FinalCounts) %in% Glyco_featurelist,]
tmp = as.numeric(unlist(Glyco_heatmap_counts))
Glyco_heatmap_counts_num = matrix(tmp,nrow = nrow(Glyco_heatmap_counts))

colnames(Glyco_heatmap_counts_num) = colnames(Glyco_heatmap_counts)
rownames(Glyco_heatmap_counts_num) = rownames(Glyco_heatmap_counts)

Glyco_heatmap_zscore = apply(Glyco_heatmap_counts_num,1,zscore)

#transpose
Glyco_heatmap_zscore = t(Glyco_heatmap_zscore)
colnames(Glyco_heatmap_zscore) = colnames(Glyco_heatmap_counts_num)

mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))

pdf("GSEA_Glyco_Heatmap_Subtypes.pdf",width = 24, height = 18)
heatmap.2(Glyco_heatmap_zscore,scale="none",Colv = F,col=mypalette,trace="none",density.info="none",key.xlab = "Z-Score Counts",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "Sample",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 1.1,cexCol = 0.15)
dev.off()
##############################################################################################################################################################

#RNA Pol II

Pol2_heatmap_counts = FinalCounts[rownames(FinalCounts) %in% Pol2_featurelist,]
tmp = as.numeric(unlist(Pol2_heatmap_counts))
Pol2_heatmap_counts_num = matrix(tmp,nrow = nrow(Pol2_heatmap_counts))

colnames(Pol2_heatmap_counts_num) = colnames(Pol2_heatmap_counts)
rownames(Pol2_heatmap_counts_num) = rownames(Pol2_heatmap_counts)

Pol2_heatmap_zscore = apply(Pol2_heatmap_counts_num,1,zscore)

#transpose
Pol2_heatmap_zscore = t(Pol2_heatmap_zscore)
colnames(Pol2_heatmap_zscore) = colnames(Pol2_heatmap_counts_num)

mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))

pdf("GSEA_Pol2_Heatmap_Subtypes.pdf",width = 24, height = 18)
heatmap.2(Pol2_heatmap_zscore,scale="none",Colv = F,col=mypalette,trace="none",density.info="none",key.xlab = "Z-Score Counts",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "Sample",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 1.1,cexCol = 0.15)
dev.off()
##############################################################################################################################################################

#Small Molecule Transport

Small_heatmap_counts = FinalCounts[rownames(FinalCounts) %in% Small_featurelist,]
tmp = as.numeric(unlist(Small_heatmap_counts))
Small_heatmap_counts_num = matrix(tmp,nrow = nrow(Small_heatmap_counts))

colnames(Small_heatmap_counts_num) = colnames(Small_heatmap_counts)
rownames(Small_heatmap_counts_num) = rownames(Small_heatmap_counts)

Small_heatmap_zscore = apply(Small_heatmap_counts_num,1,zscore)

#transpose
Small_heatmap_zscore = t(Small_heatmap_zscore)
colnames(Small_heatmap_zscore) = colnames(Small_heatmap_counts_num)

mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))

pdf("GSEA_SmallMolTransport_Heatmap_Subtypes.pdf",width = 24, height = 18)
heatmap.2(Small_heatmap_zscore,scale="none",Colv = F,col=mypalette,trace="none",density.info="none",key.xlab = "Z-Score Counts",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "Sample",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 0.6,cexCol = 0.15)
dev.off()

##############################################################################################################################################################

#ImmunoRegulatory

IReg_heatmap_counts = FinalCounts[rownames(FinalCounts) %in% IReg_featurelist,]
tmp = as.numeric(unlist(IReg_heatmap_counts))
IReg_heatmap_counts_num = matrix(tmp,nrow = nrow(IReg_heatmap_counts))

colnames(IReg_heatmap_counts_num) = colnames(IReg_heatmap_counts)
rownames(IReg_heatmap_counts_num) = rownames(IReg_heatmap_counts)

IReg_heatmap_zscore = apply(IReg_heatmap_counts_num,1,zscore)

#transpose
IReg_heatmap_zscore = t(IReg_heatmap_zscore)
colnames(IReg_heatmap_zscore) = colnames(IReg_heatmap_counts_num)

mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))

pdf("GSEA_IReg_Heatmap_Subtypes.pdf",width = 24, height = 18)
heatmap.2(IReg_heatmap_zscore,scale="none",Colv = F,col=mypalette,trace="none",density.info="none",key.xlab = "Z-Score Counts",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "Sample",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 1.1,cexCol = 0.15)
dev.off()

##############################################################################################################################################################

#GPCR

GPCR_heatmap_counts = FinalCounts[rownames(FinalCounts) %in% GPCR_featurelist,]
tmp = as.numeric(unlist(GPCR_heatmap_counts))
GPCR_heatmap_counts_num = matrix(tmp,nrow = nrow(GPCR_heatmap_counts))

colnames(GPCR_heatmap_counts_num) = colnames(GPCR_heatmap_counts)
rownames(GPCR_heatmap_counts_num) = rownames(GPCR_heatmap_counts)

GPCR_heatmap_zscore = apply(GPCR_heatmap_counts_num,1,zscore)

#transpose
GPCR_heatmap_zscore = t(GPCR_heatmap_zscore)
colnames(GPCR_heatmap_zscore) = colnames(GPCR_heatmap_counts_num)

mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))

pdf("GSEA_GPCR_Heatmap_Subtypes.pdf",width = 24, height = 18)
heatmap.2(GPCR_heatmap_zscore,scale="none",Colv = F,col=mypalette,trace="none",density.info="none",key.xlab = "Z-Score Counts",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "Sample",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 1.1,cexCol = 0.15)
dev.off()

##############################################################################################################################################################

#GTPase

GTPase_heatmap_counts = FinalCounts[rownames(FinalCounts) %in% GTPase_featurelist,]
tmp = as.numeric(unlist(GTPase_heatmap_counts))
GTPase_heatmap_counts_num = matrix(tmp,nrow = nrow(GTPase_heatmap_counts))

colnames(GTPase_heatmap_counts_num) = colnames(GTPase_heatmap_counts)
rownames(GTPase_heatmap_counts_num) = rownames(GTPase_heatmap_counts)

GTPase_heatmap_zscore = apply(GTPase_heatmap_counts_num,1,zscore)

#transpose
GTPase_heatmap_zscore = t(GTPase_heatmap_zscore)
colnames(GTPase_heatmap_zscore) = colnames(GTPase_heatmap_counts_num)

mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))

pdf("GSEA_GTPase_Heatmap_Subtypes.pdf",width = 24, height = 18)
heatmap.2(GTPase_heatmap_zscore,scale="none",Colv = F,col=mypalette,trace="none",density.info="none",key.xlab = "Z-Score Counts",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "Sample",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 1.1,cexCol = 0.15)
dev.off()

##############################################################################################################################################################

#Synaptic Transmission

Syn_heatmap_counts = FinalCounts[rownames(FinalCounts) %in% Synapse_featurelist,]
tmp = as.numeric(unlist(Syn_heatmap_counts))
Syn_heatmap_counts_num = matrix(tmp,nrow = nrow(Syn_heatmap_counts))

colnames(Syn_heatmap_counts_num) = colnames(Syn_heatmap_counts)
rownames(Syn_heatmap_counts_num) = rownames(Syn_heatmap_counts)

Syn_heatmap_zscore = apply(Syn_heatmap_counts_num,1,zscore)

#transpose
Syn_heatmap_zscore = t(Syn_heatmap_zscore)
colnames(Syn_heatmap_zscore) = colnames(Syn_heatmap_counts_num)

mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))

pdf("GSEA_SynapticTransmission_Heatmap_Subtypes.pdf",width = 24, height = 18)
heatmap.2(Syn_heatmap_zscore,scale="none",Colv = F,col=mypalette,trace="none",density.info="none",key.xlab = "Z-Score Counts",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'median'),xlab = "Sample",ylab="Gene",labRow = rownames(tmp),labCol = colnames(tmp),margins=c(8,8),cexRow = 1.1,cexCol = 0.15)
dev.off()