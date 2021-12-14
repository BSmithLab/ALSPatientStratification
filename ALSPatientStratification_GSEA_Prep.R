#Short script to obtain gene expression file ready for GSEA analysis
#GMT/GMX files can be obtained from GSEA Molecular Signatures Database
#Custom gene sets for TREM2-dependent and TREM2-indepedent disease associated microglia obtained from [1]
#Custom gene set for SOD1 microglia obtained from [2]
#.CLS phenotype file generated from subtypes assigned in previous script

#[1] Keren-Shaul, Hadas, et al. "A unique microglia type associated with restricting development of Alzheimer’s disease." Cell 169.7 (2017): 1276-1290.
#[2] Chiu, Isaac M., et al. "A neurodegeneration-specific gene-expression signature of acutely isolated microglia from an amyotrophic lateral sclerosis mouse model." Cell reports 4.2 (2013): 385-401.

#Written By Jarrett Eshima
#Dr. Barbara Smith Lab
#Sept 3 2021

wd = "C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/GSEA"
setwd(wd)

##################################################################################################################

#Remove TE and ENSEMBL Genes for Enrichment

dat = read.csv("AllSubjects_ClassifierMatrix.csv") #Generated in this script
rownames(dat) = dat[,1]
dat = dat[,-1]

rn = which(nchar(rownames(dat))>20) #use feature character number to easily filter for TEs
TE_Free = dat[-rn,]

elab = which(substr(rownames(TE_Free),1,4) == "ENSG")

ENSMBLFree = TE_Free[-elab,]


write.csv(ENSMBLFree,"AllSubjects_Classifier_GSEA_NoENSG_NoTE.csv")

#Manually convert to tab separated

##################################################################################################################
##Supplemental code to run a simple hierarchical clustering analysis on the TE features ONLY

#First, check to see if locus-specific TEs can be collapsed down to the subfamily level
##################################################################################################################

ClassMat = read.csv("ALS_HC_ClassifierMatrix.csv") #Generated in previous script
rownames(ClassMat) = ClassMat[,1]

hist(nchar(ClassMat[,1])) #Transcripts with very long names are TEs
abline(v=20,col="red") #character length threshold

ClassMat2 = ClassMat[,-1]

#Simple Hierarchical Clustering of TE features
charlen = table(nchar(rownames(ClassMat2))>20)

TE_only = data.frame(matrix(NA,nrow = charlen[[2]],ncol=ncol(ClassMat2)))
rn = which(nchar(rownames(ClassMat2))>20)

rownames(TE_only) = rownames(ClassMat2)[rn]
colnames(TE_only) = colnames(ClassMat2)

for(i in 1:nrow(ClassMat2)){
  for(j in 1:nrow(TE_only)){
    
    if(rownames(ClassMat2)[i] == rownames(TE_only)[j]){
      TE_only[j,] = ClassMat2[i,]
    }
    
  }
}

tmp = as.numeric(unlist(TE_only))
TE_only_num = matrix(tmp,nrow = nrow(TE_only),ncol = ncol(TE_only))
rownames(TE_only_num) = rownames(TE_only)
colnames(TE_only_num) = colnames(TE_only)

TEclust = hclust(dist(TE_only),method = "average")

#Manually check hierarchical clustering to see if TEs cluster by subfamily
pdf(file= "TE_Hierarchical_Clustering.pdf", width = 24, height = 18)
plot(TEclust,cex=0.3)
dev.off()

pdf(file= "TE_Hierarchical_Clustering_Heatmap.pdf", width = 28, height = 18)
mypalette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(18))
heatmap.2(TE_only_num,scale="row",Colv = F,col=mypalette,trace="none",density.info="none",key.xlab = "Z-Score Counts",dendrogram = "row",hclustfun = function(x) hclust(x,method = 'average'),xlab = "Sample",ylab="Transposable Element",labRow = rownames(TE_only_num),labCol = colnames(TE_only_num),margins=c(8,8),cexRow = 0.2,cexCol = 0.5)
dev.off()

#TEs do not appear to cluster by subfamily, enrichment of locus-specific TEs not performed in GSEA-
#bias would be introduced by collapsing locus specific TE names to TE subfamilies
#A more robust analysis can be performed using WGCNA (supplemental code)

##################################################################################################
