#######################################  Unsupervised Clustering Preparation  #########################################################
#Written by: Jarrett Eshima
#Date: April, 2021
#For: Use by the Dr. Barbara Smith Lab at Arizona State University

#Note: This short script is designed to clean up a very small number of "duplicated" gene names in the MAD10k matrix
#This script is entirely manual

#######################################################################################################################################

#HISEQ
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2")
expression = read.csv("HISEQ_ALS451_MAD10k50_10-4-21_SYMBOL.csv") #File generated in previous script

#ID duplicated Gene Names
rownames(expression) = expression[,1]

#3 genes are duplicated: 'CSF2RA', 'P2RY8', 'TBC1D26' 
#2 Genuine Duplicates: 'CSF2RA', 'P2RY8'

#Check that gene is genuinely duplicated
which(expression[,1] == "CSF2RA") #duplicate
expression[4841:4842,]
which(expression[,1] == "P2RY8") #duplicate
expression[671:672,]
which(expression[,1] == "TBC1D26") #not duplicated
expression[2769,2:5]
expression[4245,2:5]

#Use code above to verify the gene is actually duplicated (exact same count values)
duplicatedgenes = c(4842,672)
expression= expression[-duplicatedgenes,]

rownames(expression) = make.unique(expression[,1])
expression = expression[,-1]
rownames(expression) = make.unique(gsub('\\..*','',rownames(expression)))

fexpression = expression

fexpression = fexpression[-10001:-100048,] #Remove the extra genes included for duplicate removal

file = "FINAL_HISEQ_ALS451_MAD10k_10-4-21_SYMBOL.csv"
write.csv(fexpression,file)
file = "FINAL_HISEQ_ALS451_MAD10k_10-4-21_SYMBOL.txt"
write.table(fexpression,file,sep="\t",quote=F,row.names=T,col.names=T)


####################################################################################################################################
####################################################################################################################################


#NOVASEQ
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2")
expression = read.csv("NOVASEQ_ALS451_MAD10k50_10-4-21_SYMBOL.csv") #File generated in previous script


#ID duplicated Gene Names
rownames(expression) = expression[,1]

#8 genes are duplicated: 'ASMTL-AS1', 'CSF2RA', 'DHRSX-IT1', 'IL9R', 'P2RY8', 'RNA5-8S5', 'TBC1D26', 'WASH6P'
#13 Genuine Duplicates: 'ASMTL-AS1', 'CSF2RA', 'DHRSX-IT1', 'IL9R', 'P2RY8', 'RNA5-8S5' (X7), 'WASH6P'

#Check that gene is genuinely duplicated
which(expression[,1] == "ASMTL-AS1") #duplicate
expression[7025:7026,]
which(expression[,1] == "CSF2RA") #duplicate
expression[4896:4897,]
which(expression[,1] == "DHRSX-IT1") #duplicate
expression[7362:7363,]
which(expression[,1] == "IL9R") #duplicate
expression[3848:3849,]
which(expression[,1] == "P2RY8") #duplicate
expression[4849:4850,]
which(expression[,1] == "RNA5-8S5") #duplicated (X7)
expression[710:717,10:15]
which(expression[,1] == "TBC1D26") #not duplicated
expression[2223,2:4]
expression[2640,2:4]
which(expression[,1] == "WASH6P") #duplicated
expression[7931:7932,]

#Use code above to verify the gene is actually duplicated (exact same count values)
duplicatedgenes = c(7026,4897,7363,3849,4850,711,712,713,714,715,716,717,7932)
expression= expression[-duplicatedgenes,]

rownames(expression) = make.unique(expression[,1])
expression = expression[,-1]
rownames(expression) = make.unique(gsub('\\..*','',rownames(expression)))

fexpression = expression

fexpression = fexpression[-10001:-100037,] #Remove the extra genes included for duplicate removal

file = "FINAL_NOVASEQ_ALS451_MAD10k_10-4-21_SYMBOL.csv"
write.csv(fexpression,file)
file = "FINAL_NOVASEQ_ALS451_MAD10k_10-4-21_SYMBOL.txt"
write.table(fexpression,file,sep="\t",quote=F,row.names=T,col.names=T)
