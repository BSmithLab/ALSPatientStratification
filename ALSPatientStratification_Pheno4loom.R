#Short script to combine sample Phenotype Data and subtype labels into file for Loom generation

#Written By: Jarrett Eshima
#For: Dr. Barbara Smith Lab

load("D:/Jarrett/Research/Summer 2021/SQuIRE Count ALS/FINAL/DESeq2Completed_10-4-21.RData") #File generated in previous script
ALS_Pheno = ALScoldata_451_TE 

#You can also read in the ColData file available in Prudencio et al. supplemental files
ALS_Pheno = read.csv("ColData.csv")

blank = rep(NA,nrow(ALS_Pheno))

ALS_Pheno = cbind(ALS_Pheno,blank)
ALS_Pheno = cbind(ALS_Pheno,blank)

colnames(ALS_Pheno)[9:10] = c("Subtype","Factor")

#NovaSEQ Robust Subtype Assignment

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2/NovaSeq")
NovaSubtype = read.csv("NovaSeq_RobustSubtypeAssignment_10-11-21_majority.csv") #File generated in previous script
rownames(NovaSubtype) = NovaSubtype[,1]
NovaSubtype = NovaSubtype[,-1]

#HiSEQ Robust Subtype Assignment

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2/HiSeq")
HiSubtype = read.csv("HiSeq_RobustSubtypeAssignment_10-11-21_majority.csv") #File generated in previous script
rownames(HiSubtype) = HiSubtype[,1]
HiSubtype = HiSubtype[,-1]

for(i in 1:nrow(ALS_Pheno)){
  for(j in 1:ncol(NovaSubtype)){
    if(rownames(ALS_Pheno)[i] == colnames(NovaSubtype)[j]){
      ALS_Pheno$Factor[i] = NovaSubtype[13,j]
    }
  }
}

for(i in 1:nrow(ALS_Pheno)){
  for(j in 1:ncol(HiSubtype)){
    if(rownames(ALS_Pheno)[i] == colnames(HiSubtype)[j]){
      ALS_Pheno$Factor[i] = HiSubtype[13,j]
    }
  }
}

#TE = 1
#OX = 2
#GLIA = 3

for(i in 1:nrow(ALS_Pheno)){
  if(ALS_Pheno$Factor[i] == 1){
    ALS_Pheno$Subtype[i] = "TE"
  }else if(ALS_Pheno$Factor[i] == 2){
    ALS_Pheno$Subtype[i] = "OX"
  }else{
    ALS_Pheno$Subtype[i] = "GLIA"
  }
}

table(ALS_Pheno$Subtype)


##### Extra code to have samples grouped by sequencing platform ######

plat = table(ALS_Pheno$sequencing_platform)

NovaPheno = data.frame(matrix(NA,nrow=plat[[2]],ncol=ncol(ALS_Pheno)))
HiPheno = data.frame(matrix(NA,nrow=plat[[1]],ncol=ncol(ALS_Pheno)))
count1 = count2 = 1
for(i in 1:nrow(ALS_Pheno)){
  if(ALS_Pheno$sequencing_platform[i] == "HiSeq"){
    HiPheno[count1,] = ALS_Pheno[i,]
    rownames(HiPheno)[count1] = rownames(ALS_Pheno)[i]
    count1 = count1+1
  }else if (ALS_Pheno$sequencing_platform[i] == "NovaSeq"){
    NovaPheno[count2,] = ALS_Pheno[i,]
    rownames(NovaPheno)[count2] = rownames(ALS_Pheno)[i]
    count2 = count2+1
  }
}

colnames(NovaPheno) = colnames(ALS_Pheno)
colnames(HiPheno) = colnames(ALS_Pheno)

setwd('C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2')
write.csv(NovaPheno,"Nova_ALS_coldata_SUBTYPES.csv")
write.csv(HiPheno,"Hi_ALS_coldata_SUBTYPES.csv")


#Combined phenotype file
Final_Pheno = rbind(NovaPheno,HiPheno) #NovaSeq samples first
colnames(Final_Pheno) = colnames(ALS_Pheno)
setwd('C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2')
write.csv(Final_Pheno,"ALS451_coldata_SUBTYPES.csv")
