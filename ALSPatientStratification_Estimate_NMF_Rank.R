############################################  Estimate Rank for NMF Clustering  #################################################

#Notes written by: Jarrett Eshima
#Date: April, 2021
#For: Use by the Dr. Barbara Smith Lab at Arizona State University

###############################################################################################################################
#User Notes:
#This should be run before running NMF 'in earnest'.
#This code will allow you to estimate the matrix rank (number of clusters, k)
#Once you identify the proper rank, proceed to the SAKE GUI to run 'real' NMF analysis

#This is much more costly than the code appears, each rank (number of clusters) is computed 'nrun' times
#Each run gets a different initialization values for matrices W and H
#The cophenetic correlation coefficient and other metrics (such as dispersion) will inform you on the reproducibility of your clusters/groups (indirectly) 

#It is recommended to set nrun to at least 30-50 to get a reasonable estimate for performance metrics like coph corr
#Source: Burnet 2010 (https://www.pnas.org/content/pnas/101/12/4164.full.pdf)
#Source: Hutchins 2008 (https://academic.oup.com/bioinformatics/article/24/23/2684/180798?login=true)

# Rank Estimation References: 
# https://cran.r-project.org/web/packages/NMF/vignettes/NMF-vignette.pdf

###############################################################################################################################
#Load dependent packages
library(Biobase)
library(destiny)
library(NMF)

###############################################################################################################################

#NovaSeq Matrix
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2")
data = read.csv("FINAL_NOVASEQ_ALS451_MAD10k_10-4-21_SYMBOL.csv") #File generated in previous script
rownames(data) = data[,1]
data = data[,-1]

tdata = t(data)
tdata = data.frame(tdata)

estrank_nova = as.ExpressionSet(tdata)
dim(estrank_nova)

minrank = 2
maxrank = 6

if(requireNamespace("Biobase",quietly = T)){
  estim.r.nova = nmf(estrank_nova,minrank:maxrank,nrun=50,method = "nsNMF") #Don't include a seed argument to get the default method, "Random"
}

plot(estim.r.nova)


#Check Overfitting
#par(mar = c(6,2,2,2))
sizeGrWindow(12,9)
par(mar = c(0,4,2,0))
if(requireNamespace("Biobase", quietly=TRUE)){
  consensusmap(estim.r.nova, annCol=estrank_nova, labCol=NA, labRow=NA)
}

#Finalized Data
#save.image("GSE153960_NMF_RankCalculation_10-4-21.RData")
#load("GSE153960_NMF_RankCalculation_10-4-21.RData")

###############################################################################################################################

#HiSeq Matrix
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2")
data = read.csv("FINAL_HISEQ_ALS451_MAD10k_10-4-21_SYMBOL.csv") #File generated in previous script
rownames(data) = data[,1]
data = data[,-1]

tdata = t(data)
tdata = data.frame(tdata)

estrank_hi = as.ExpressionSet(tdata)

minrank = 2
maxrank = 6

if(requireNamespace("Biobase",quietly = T)){
  estim.r.hi = nmf(estrank_hi,minrank:maxrank,nrun=50,method = "nsNMF") #Don't include a seed argument to get the default method, "Random"
}

plot(estim.r.hi)


#Check Overfitting
if(requireNamespace("Biobase", quietly=TRUE)){
  consensusmap(estim.r.hi, annCol=estrank_hi, labCol=NA, labRow=NA)
}


#save.image("GSE153960_NMF_RankCalculation_10-4-21.RData")
#load("GSE153960_NMF_RankCalculation_10-4-21.RData")

