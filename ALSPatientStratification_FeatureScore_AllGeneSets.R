##Use NMF 'feature scores' to determine the subset of MAD genes to use for supervised classification and GSEA

#This script generates all feature sets considered for development of supervised classification models
#The final feature set can be generated using ALSPatientStratification_FeatureScore_v2.R

#Written By: Jarrett Eshima
#For: Dr. Barbara Smith Lab

## Read in "totfeatures" files. These are obtained from SAKE by downloading all feature scores following NMF clustering.

NumReplicates = 11

#Read in the data
for(i in 1:NumReplicates){
  setwd(paste("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2/NovaSeq/",seq(1,NumReplicates,1),sep="")[i])
  tmp = read.csv('totfeatures.csv') #File generated in previous script
  rownames(tmp) = tmp$GeneCard
  tmp = tmp[,-1]
  repname = paste("NovaRep",i,sep="")
  assign(repname,tmp)
}

NumReplicates = 11

for(i in 1:NumReplicates){
  setwd(paste("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2/HiSeq/",seq(1,NumReplicates,1),sep="")[i])
  tmp = read.csv('totfeatures.csv') #File generated in previous script
  rownames(tmp) = tmp$GeneCard
  tmp = tmp[,-1]
  repname = paste("HiRep",i,sep="")
  assign(repname,tmp)
}


#subset out data based on averaged feature score across 11 runs
#both top x features for each subtype and top x features overall

#NovaSeq
AverageFeatureScoreNova = matrix(NA,nrow=nrow(NovaRep1),ncol = 2)
rownames(AverageFeatureScoreNova) = rownames(NovaRep1)
colnames(AverageFeatureScoreNova) = c("AverageScore","Rep1Group")
AverageFeatureScoreNova = data.frame(AverageFeatureScoreNova)

for(i in 1:nrow(AverageFeatureScoreNova)){
  r1 = which(rownames(NovaRep1) == rownames(AverageFeatureScoreNova)[i])
  r2 = which(rownames(NovaRep2) == rownames(AverageFeatureScoreNova)[i])
  r3 = which(rownames(NovaRep3) == rownames(AverageFeatureScoreNova)[i])
  r4 = which(rownames(NovaRep4) == rownames(AverageFeatureScoreNova)[i])
  r5 = which(rownames(NovaRep5) == rownames(AverageFeatureScoreNova)[i])
  r6 = which(rownames(NovaRep6) == rownames(AverageFeatureScoreNova)[i])
  r7 = which(rownames(NovaRep7) == rownames(AverageFeatureScoreNova)[i])
  r8 = which(rownames(NovaRep8) == rownames(AverageFeatureScoreNova)[i])
  r9 = which(rownames(NovaRep9) == rownames(AverageFeatureScoreNova)[i])
  r10 = which(rownames(NovaRep10) == rownames(AverageFeatureScoreNova)[i])
  r11 = which(rownames(NovaRep11) == rownames(AverageFeatureScoreNova)[i])

  tmp1 = mean(c(NovaRep1$featureScore[r1],NovaRep2$featureScore[r2],NovaRep3$featureScore[r3],NovaRep4$featureScore[r4],NovaRep5$featureScore[r5],NovaRep6$featureScore[r6],NovaRep7$featureScore[r7],NovaRep8$featureScore[r8],NovaRep9$featureScore[r9],NovaRep10$featureScore[r10]),NovaRep11$featureScore[r11])
  AverageFeatureScoreNova[i,1] = tmp1 #the average feature score across NMF replicates
  AverageFeatureScoreNova[i,2] = NovaRep1$Group[r1] #the NMF cluster of the feature from the first run only

  if((i %% 500) == 0) cat("% Done:",i/nrow(AverageFeatureScoreNova)*100,"\n")
}

orderNova = AverageFeatureScoreNova[order(AverageFeatureScoreNova$AverageScore,decreasing = T),]

head(orderNova)


#HiSeq
AverageFeatureScoreHi = matrix(NA,nrow=nrow(HiRep1),ncol = 2)
rownames(AverageFeatureScoreHi) = rownames(HiRep1)
colnames(AverageFeatureScoreHi) = c("AverageScore","Rep1Group")
AverageFeatureScoreHi = data.frame(AverageFeatureScoreHi)

for(i in 1:nrow(AverageFeatureScoreHi)){
  r1 = which(rownames(HiRep1) == rownames(AverageFeatureScoreHi)[i])
  r2 = which(rownames(HiRep2) == rownames(AverageFeatureScoreHi)[i])
  r3 = which(rownames(HiRep3) == rownames(AverageFeatureScoreHi)[i])
  r4 = which(rownames(HiRep4) == rownames(AverageFeatureScoreHi)[i])
  r5 = which(rownames(HiRep5) == rownames(AverageFeatureScoreHi)[i])
  r6 = which(rownames(HiRep6) == rownames(AverageFeatureScoreHi)[i])
  r7 = which(rownames(HiRep7) == rownames(AverageFeatureScoreHi)[i])
  r8 = which(rownames(HiRep8) == rownames(AverageFeatureScoreHi)[i])
  r9 = which(rownames(HiRep9) == rownames(AverageFeatureScoreHi)[i])
  r10 = which(rownames(HiRep10) == rownames(AverageFeatureScoreHi)[i])
  r11 = which(rownames(HiRep11) == rownames(AverageFeatureScoreHi)[i])
  
  tmp1 = mean(c(HiRep1$featureScore[r1],HiRep2$featureScore[r2],HiRep3$featureScore[r3],HiRep4$featureScore[r4],HiRep5$featureScore[r5],HiRep6$featureScore[r6],HiRep7$featureScore[r7],HiRep8$featureScore[r8],HiRep9$featureScore[r9],HiRep10$featureScore[r10],HiRep11$featureScore[r11]))
  AverageFeatureScoreHi[i,1] = tmp1
  AverageFeatureScoreHi[i,2] = HiRep1$Group[r1]
  
  if((i %% 500) == 0) cat("% Done:",i/nrow(AverageFeatureScoreHi)*100,"\n")
}

orderHi = AverageFeatureScoreHi[order(AverageFeatureScoreHi$AverageScore,decreasing = T),]

head(orderHi)



#Top features overall - NovaSeq
Ntop250ALL = rownames(orderNova)[1:250]
Ntop500ALL = rownames(orderNova)[1:500]
Ntop1000ALL = rownames(orderNova)[1:1000]
Ntop2000ALL = rownames(orderNova)[1:2000]
Ntop5000ALL = rownames(orderNova)[1:5000]
#N10kALL = rownames(orderNova)

#Top features overall - HiSeq
Htop250ALL = rownames(orderHi)[1:250]
Htop500ALL = rownames(orderHi)[1:500]
Htop1000ALL = rownames(orderHi)[1:1000]
Htop2000ALL = rownames(orderHi)[1:2000]
Htop5000ALL = rownames(orderHi)[1:5000]
#H10kALL = rownames(orderHi)


#Top features by group - NovaSeq

count1=count2=count3=1
G1 = G2 = G3 = data.frame(matrix(NA,nrow(orderNova),ncol(orderNova)))
colnames(G1)=colnames(G2)=colnames(G3)=colnames(orderNova)

for(i in 1:nrow(orderNova)){
  if(orderNova$Rep1Group[i] == 1){
    G1[count1,] = orderNova[i,]
    tmprn1 = rownames(orderNova)[i]
    rownames(G1)[count1] = tmprn1
    count1=count1+1
  }else if(orderNova$Rep1Group[i] == 2){
    G2[count2,] = orderNova[i,]
    tmprn2 = rownames(orderNova)[i]
    rownames(G2)[count2] = tmprn2
    count2 = count2+1
  }else{
    G3[count3,] = orderNova[i,]
    tmprn3 = rownames(orderNova)[i]
    rownames(G3)[count3] = tmprn3
    count3 = count3+1
  }
}

remi1 = which(is.na(G1$AverageScore))
remi2 = which(is.na(G2$AverageScore))
remi3 = which(is.na(G3$AverageScore))
G1 = G1[-remi1,]
G2 = G2[-remi2,]
G3 = G3[-remi3,]

Ntop99group = c(rownames(G1)[1:33],rownames(G2)[1:33],rownames(G3)[1:33])
Ntop250group = c(rownames(G1)[1:250],rownames(G2)[1:250],rownames(G3)[1:250])
Ntop500group = c(rownames(G1)[1:500],rownames(G2)[1:500],rownames(G3)[1:500])
Ntop1000group = c(rownames(G1)[1:1000],rownames(G2)[1:1000],rownames(G3)[1:1000])
Ntop2000group = c(rownames(G1)[1:2000],rownames(G2)[1:2000],rownames(G3)[1:2000])
#Ntop5000group = c(rownames(G1)[1:5000],rownames(G2)[1:5000],rownames(G3)[1:5000]) #This is basically all the features

#Top features by group - HiSeq

count1=count2=count3=1
G1 = G2 = G3 = data.frame(matrix(NA,nrow(orderHi),ncol(orderHi)))
colnames(G1)=colnames(G2)=colnames(G3)=colnames(orderHi)

for(i in 1:nrow(orderHi)){
  if(orderHi$Rep1Group[i] == 1){
    G1[count1,] = orderHi[i,]
    tmprn1 = rownames(orderHi)[i]
    rownames(G1)[count1] = tmprn1
    count1=count1+1
  }else if(orderHi$Rep1Group[i] == 2){
    G2[count2,] = orderHi[i,]
    tmprn2 = rownames(orderHi)[i]
    rownames(G2)[count2] = tmprn2
    count2 = count2+1
  }else{
    G3[count3,] = orderHi[i,]
    tmprn3 = rownames(orderHi)[i]
    rownames(G3)[count3] = tmprn3
    count3 = count3+1
  }
}

remi1 = which(is.na(G1$AverageScore))
remi2 = which(is.na(G2$AverageScore))
remi3 = which(is.na(G3$AverageScore))
G1 = G1[-remi1,]
G2 = G2[-remi2,]
G3 = G3[-remi3,]

Htop99group = c(rownames(G1)[1:33],rownames(G2)[1:33],rownames(G3)[1:33])
Htop250group = c(rownames(G1)[1:250],rownames(G2)[1:250],rownames(G3)[1:250])
Htop500group = c(rownames(G1)[1:500],rownames(G2)[1:500],rownames(G3)[1:500])
Htop1000group = c(rownames(G1)[1:1000],rownames(G2)[1:1000],rownames(G3)[1:1000])
Htop2000group = c(rownames(G1)[1:2000],rownames(G2)[1:2000],rownames(G3)[1:2000]) 
#Htop5000group = c(rownames(G1)[1:5000],rownames(G2)[1:5000],rownames(G3)[1:5000])

### Combine NovaSeq and HiSeq feature sets

#Combined All
Ctop250ALL = names(table(c(Ntop250ALL,Htop250ALL)))
Ctop500ALL = names(table(c(Ntop500ALL,Htop500ALL)))
Ctop1000ALL = names(table(c(Ntop1000ALL,Htop1000ALL)))
Ctop2000ALL = names(table(c(Ntop2000ALL,Htop2000ALL)))
Ctop5000ALL = names(table(c(Ntop5000ALL,Htop5000ALL)))

#Combined Group
Ctop99group = names(table(c(Ntop99group,Htop99group)))
Ctop250group = names(table(c(Ntop250group,Htop250group)))
Ctop500group = names(table(c(Ntop500group,Htop500group)))
Ctop1000group = names(table(c(Ntop1000group,Htop1000group)))
Ctop2000group = names(table(c(Ntop2000group,Htop2000group)))
#Ctop5000group = names(table(c(Ntop5000group,Htop5000group)))

#Overlapping All
Otop250ALL = Ntop250ALL[Ntop250ALL %in% Htop250ALL]
Otop500ALL = Ntop500ALL[Ntop500ALL %in% Htop500ALL]
Otop1000ALL = Ntop1000ALL[Ntop1000ALL %in% Htop1000ALL]
Otop2000ALL = Ntop2000ALL[Ntop2000ALL %in% Htop2000ALL]
Otop5000ALL = Ntop5000ALL[Ntop5000ALL %in% Htop5000ALL]

#Overlapping Group
Otop99group = Ntop99group[Ntop99group %in% Htop99group]
Otop250group = Ntop250group[Ntop250group %in% Htop250group]
Otop500group = Ntop500group[Ntop500group %in% Htop500group]
Otop1000group = Ntop1000group[Ntop1000group %in% Htop1000group]
Otop2000group = Ntop2000group[Ntop2000group %in% Htop2000group]
#Otop5000group = Ntop5000group[Ntop5000group %in% Htop5000group]


## GLOBAL OVERLAP
#Top X features based on average score across sequencing platforms - ALL
alloverlappinggenes = orderNova[rownames(orderNova) %in% rownames(orderHi),]
alloverlappinggenes = rownames(alloverlappinggenes)

NovaOverlap = orderNova[rownames(orderNova) %in% alloverlappinggenes,]
HiOverlap = orderHi[rownames(orderHi) %in% alloverlappinggenes,]

COverlap = matrix(NA,nrow=nrow(NovaOverlap),ncol = ncol(NovaOverlap))
rownames(COverlap) = rownames(NovaOverlap) ; colnames(COverlap) = colnames(NovaOverlap)

for(i in 1:nrow(COverlap)){
  NovaIndex = which(rownames(COverlap)[i] == rownames(NovaOverlap))
  HiIndex = which(rownames(COverlap)[i] == rownames(HiOverlap))
  
  COverlap[i,1] = mean(c(NovaOverlap$AverageScore[NovaIndex],HiOverlap$AverageScore[HiIndex]))
}

COverlap = data.frame(COverlap)

orderCOverlapAll = COverlap[order(COverlap$AverageScore,decreasing = T),]

COtop500ALL = rownames(orderCOverlapAll)[1:500]
COtop1000ALL = rownames(orderCOverlapAll)[1:1000]
COtop2000ALL = rownames(orderCOverlapAll)[1:2000]
COtop2500ALL = rownames(orderCOverlapAll)[1:2500]
COtop5000ALL = rownames(orderCOverlapAll)[1:5000]

#####################################################################################################################################################################################
#ENSMBL LUT for VST Files (both VST files required because MAD10k and other processed files contain platform-specific genes)

#Get ENSEMBL ID and Symbol LUT for novaseq and hiseq platforms (and combine)
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2/PreProcessing")

mad10000.ens = read.csv("NOVASEQ_ALS451_VST_Gene_TE_10-4-21.csv",header = T) #File generated in previous script
rownames(mad10000.ens) = mad10000.ens[,1]
mad10000.ens = mad10000.ens[,-1]
egenes = rownames(mad10000.ens)
newEG = sub("\\..*","",egenes) #Keep text to the left of the dot
ensembl_version = "https://dec2016.archive.ensembl.org" #Updated: https://jan2020.archive.ensembl.org
species="human"
ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host=ensembl_version)
gene_positions_nova <- biomaRt::getBM(filters = 'ensembl_gene_id',attributes=c('ensembl_gene_id','hgnc_symbol'), values = newEG, mart = ensembl)

mad10000.ens = read.csv("HISEQ_ALS451_VST_Gene_TE_10-4-21.csv",header=T) #File generated in previous script
rownames(mad10000.ens) = mad10000.ens[,1]
mad10000.ens = mad10000.ens[,-1]
egenes = rownames(mad10000.ens)
newEG = sub("\\..*","",egenes) #Keep text to the left of the dot
ensembl_version = "https://dec2016.archive.ensembl.org" #Updated: https://jan2020.archive.ensembl.org
species="human"
ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host=ensembl_version)
gene_positions_hi <- biomaRt::getBM(filters = 'ensembl_gene_id',attributes=c('ensembl_gene_id','hgnc_symbol'), values = newEG, mart = ensembl)

ke = which(! gene_positions_nova$ensembl_gene_id %in% gene_positions_hi$ensembl_gene_id) #Only finds the genes found in the HiSeq Gene Symbol LUT but not the NovaSeq Gene Symbol LUT
ke2 = gene_positions_nova[ke,]
Combined_gene_positions = rbind(gene_positions_hi,ke2)


#####################################################################################################################################################################################

#Parse VST Datasets for Classifier Development
#Completed: HISEQ
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2/PreProcessing")
filen = "NOVASEQ_ALS451_VST_Gene_TE_10-4-21.csv"
vstcounts = read.csv(filen) #File generated in previous script
rownames(vstcounts) = vstcounts[,1]
vstcounts = vstcounts[,-1]


#Covert Ensembl rownames to gene symbol rownames (10,000+ names X 10,000+ names = 100 million if statements)
blank = rep(NA,nrow(vstcounts))
vstsymbol = cbind(blank,vstcounts)
vstsymbol[,1] = rownames(vstcounts)

#This loop takes a very long time, it is recommended to save the global environment after completion
DONE = F
for(i in 1:nrow(vstsymbol)){
  for(j in 1:nrow(Combined_gene_positions)){
    if(DONE == F){
      if(sub("\\..*","",vstsymbol[i,1]) == Combined_gene_positions$ensembl_gene_id[j]){
        vstsymbol[i,1] = Combined_gene_positions$hgnc_symbol[j]
        DONE = T
      }
    }
  }
  DONE = F
  if((i %% 100) == 0) cat("% Done:",i/nrow(vstsymbol)*100,"\n")
}

#Replace missing gene symbol with original ENSEMBL ID
for(i in 1:nrow(vstsymbol)){
  if(vstsymbol[i,1] == ""){
    vstsymbol[i,1] = rownames(vstsymbol)[i]
  }
}


finalvstsymbol = vstsymbol

rownames(finalvstsymbol) = make.unique(finalvstsymbol[,1])
finalvstsymbol = finalvstsymbol[,-1]
rownames(finalvstsymbol) = make.unique(gsub('\\..*','',rownames(finalvstsymbol)))
finalvstsymbol =data.frame(finalvstsymbol)


wd = 'C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2/'
setwd(wd)

#FullNovaSymbol = finalvstsymbol
#write.csv(FullNovaSymbol,"Full_Nova_Matrix_VST_8-23-21.csv")
#FullHiSymbol = finalvstsymbol
#write.csv(FullHiSymbol,"Full_Hiseq_Matrix_VST_8-23-21.csv")


#save.image("SAKEfeatures2ClassifierMatrix_v2.RData")
load("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2/SAKEfeatures2ClassifierMatrix_v2.RData")


##Correct all feature set lists for excel autoformat driven errors (MARCH3 --> 3-Mar)
myfeaturesets = c("COtop500ALL","COtop1000ALL","COtop2000ALL","COtop2500ALL","COtop5000ALL", "Ctop250ALL","Ctop500ALL","Ctop1000ALL","Ctop2000ALL","Ctop5000ALL","Ctop99group","Ctop250group","Ctop500group","Ctop1000group","Ctop2000group","Otop250ALL","Otop500ALL","Otop1000ALL","Otop2000ALL","Otop5000ALL","Otop99group","Otop250group","Otop500group","Otop1000group","Otop2000group")
excelgenelist = c("1-Sep","2-Sep","3-Sep","4-Sep","5-Sep","6-Sep","7-Sep","8-Sep","9-Sep","10-Sep","11-Sep","12-Sep","1-Mar","2-Mar","3-Mar","4-Mar","5-Mar","6-Mar","7-Mar","8-Mar","9-Mar","10-Mar","11-Mar","12-Mar")
truegenelist = c("SEPT1","SEPT2","SEPT3","SEPT4","SEPT5","SEPT6","SEPT7","SEPT8","SEPT9","SEPT10","SEPT11","SEPT12","MARCH1","MARCH2","MARCH3","MARCH4","MARCH5","MARCH6","MARCH7","MARCH8","MARCH9","MARCH10","MARCH11","MARCH12")

ExcelerrorLUT = data.frame(cbind(excelgenelist,truegenelist))

for(i in 1:length(myfeaturesets)){
  
  tmp = get(myfeaturesets[i])
  
  for(j in 1:length(tmp)){
    for(k in 1:nrow(ExcelerrorLUT)){
      
      if(tmp[j] == ExcelerrorLUT$excelgenelist[k]){
        
        tmp[j] = ExcelerrorLUT$truegenelist[k]
        
      }
      
    }
  }
  nam = myfeaturesets[i]
  assign(paste(nam),tmp)
  if((i %% 1) == 0) cat("Feature Set Completed:",myfeaturesets[i],"\n")
} 


#FullHiSymbol2 = FullHiSymbol[rownames(FullHiSymbol) %in% rownames(FullNovaSymbol),]
#FullNovaSymbol2 = FullNovaSymbol[rownames(FullNovaSymbol) %in% rownames(FullHiSymbol),]

#Rownames must be in the same order - TRUE for version 2 matrices
#table(rownames(FullNovaSymbol2) == rownames(FullHiSymbol2))


#A handful of SAKE-important genes for one platform are sex-dependent genes on the other platform.
#Examples: CDH1, VAMP8, IGFLR1
#Solution: Genes that are sex-dependent on either platform are excluded from the supervised classification analysis and GSEA

#Remove platform-specific sex-dependent genes - NOVASEQ
MAD10k = FullNovaSymbol

platdif1 = which(! Ctop5000ALL %in% rownames(MAD10k))
Ctop5000ALL2 = Ctop5000ALL[-platdif1]
platdif2 = which(! Ctop500ALL %in% rownames(MAD10k))
Ctop500ALL2 = Ctop500ALL[-platdif2]
platdif3 = which(! Ctop1000ALL %in% rownames(MAD10k))
Ctop1000ALL2 = Ctop1000ALL[-platdif3]
platdif4 = which(! Ctop2000ALL %in% rownames(MAD10k))
Ctop2000ALL2 = Ctop2000ALL[-platdif4]

#platdif5 = which(!Ctop99group %in% rownames(MAD10k)) #No missing genes
#Ctop99group2 = Ctop99group[-platdif5]
platdif6 = which(! Ctop250group %in% rownames(MAD10k))
Ctop250group2 = Ctop250group[-platdif6]
platdif7 = which(! Ctop500group %in% rownames(MAD10k))
Ctop500group2 = Ctop500group[-platdif7]
platdif8 = which(! Ctop1000group %in% rownames(MAD10k))
Ctop1000group2 = Ctop1000group[-platdif8]
platdif9 = which(! Ctop2000group %in% rownames(MAD10k))
Ctop2000group2 = Ctop2000group[-platdif9]


#Remove platform-specific sex-dependent genes - HISEQ
MAD10k = FullHiSymbol

h1 = which(! Ctop5000ALL2 %in% rownames(MAD10k))
Ctop5000ALL3 = Ctop5000ALL2[-h1]
h2 = which(! Ctop500ALL2 %in% rownames(MAD10k))
Ctop500ALL3 = Ctop500ALL2[-h2]
h3 = which(! Ctop1000ALL2 %in% rownames(MAD10k))
Ctop1000ALL3 = Ctop1000ALL2[-h3]
h4 = which(! Ctop2000ALL2 %in% rownames(MAD10k))
Ctop2000ALL3 = Ctop2000ALL2[-h4]

#h5 = which(!Ctop99group %in% rownames(MAD10k)) #No missing genes
#Ctop99group3 = Ctop99group2[-h5]
h6 = as.numeric(which(! Ctop250group2 %in% rownames(MAD10k)))
Ctop250group3 = Ctop250group2[-h6]
h7 = which(! Ctop500group2 %in% rownames(MAD10k))
Ctop500group3 = Ctop500group2[-h7]
h8 = which(! Ctop1000group2 %in% rownames(MAD10k))
Ctop1000group3 = Ctop1000group2[-h8]
h9 = which(! Ctop2000group2 %in% rownames(MAD10k))
Ctop2000group3 = Ctop2000group2[-h9]



##NovaSeq
#Combined All (duplicates removed)
MAD10k = FullNovaSymbol
CA500N = MAD10k[rownames(MAD10k)%in% Ctop500ALL3,]
CA1000N = MAD10k[rownames(MAD10k) %in% Ctop1000ALL3,]
CA2000N = MAD10k[rownames(MAD10k) %in% Ctop2000ALL3,]
CA5000N =  MAD10k[rownames(MAD10k) %in% Ctop5000ALL3,]
#Combined Group (duplicates removed)
CG99N = MAD10k[rownames(MAD10k) %in% Ctop99group,]
CG250N = MAD10k[rownames(MAD10k) %in% Ctop250group3,]
CG500N = MAD10k[rownames(MAD10k) %in% Ctop500group3,]
CG1000N = MAD10k[rownames(MAD10k) %in% Ctop1000group3,]
CG2000N = MAD10k[rownames(MAD10k) %in% Ctop2000group3,]
#Overlapping All (overlapping genes in a specific feature set)
#OA250N = MAD10k[rownames(MAD10k) %in% Otop250ALL,] #Too few features
OA500N = MAD10k[rownames(MAD10k) %in% Otop500ALL,]
OA1000N = MAD10k[rownames(MAD10k) %in% Otop1000ALL,]
OA2000N = MAD10k[rownames(MAD10k) %in% Otop2000ALL,]
OA5000N = MAD10k[rownames(MAD10k) %in% Otop5000ALL,]
#Overlapping Group (overlapping genes in a specific feature set)
#OG99N = MAD10k[rownames(MAD10k) %in% Otop99group,] #Too few features
OG250N = MAD10k[rownames(MAD10k) %in% Otop250group,]
OG500N = MAD10k[rownames(MAD10k) %in% Otop500group,]
OG1000N = MAD10k[rownames(MAD10k) %in% Otop1000group,]
OG2000N = MAD10k[rownames(MAD10k) %in% Otop2000group,]
#Global Overlap
GO500N = MAD10k[rownames(MAD10k) %in% COtop500ALL,]
GO1000N = MAD10k[rownames(MAD10k) %in% COtop1000ALL,]
GO2000N = MAD10k[rownames(MAD10k) %in% COtop2000ALL,]
GO2500N = MAD10k[rownames(MAD10k) %in% COtop2500ALL,]
GO5000N = MAD10k[rownames(MAD10k) %in% COtop5000ALL,]


##Hiseq
#Combined All (duplicates removed)
MAD10k = FullHiSymbol
CA500H = MAD10k[rownames(MAD10k) %in% Ctop500ALL3,]
CA1000H = MAD10k[rownames(MAD10k) %in% Ctop1000ALL3,]
CA2000H = MAD10k[rownames(MAD10k) %in% Ctop2000ALL3,]
CA5000H =  MAD10k[rownames(MAD10k) %in% Ctop5000ALL3,]
#Combined Group (duplicates removed)
CG99H = MAD10k[rownames(MAD10k) %in% Ctop99group,]
CG250H = MAD10k[rownames(MAD10k) %in% Ctop250group3,]
CG500H = MAD10k[rownames(MAD10k) %in% Ctop500group3,]
CG1000H = MAD10k[rownames(MAD10k) %in% Ctop1000group3,]
CG2000H = MAD10k[rownames(MAD10k) %in% Ctop2000group3,]
#Overlapping All (overlapping genes in a specific feature set)
#OA250H = MAD10k[rownames(MAD10k) %in% Otop250ALL,] #Too few features
OA500H = MAD10k[rownames(MAD10k) %in% Otop500ALL,]
OA1000H = MAD10k[rownames(MAD10k) %in% Otop1000ALL,]
OA2000H = MAD10k[rownames(MAD10k) %in% Otop2000ALL,]
OA5000H = MAD10k[rownames(MAD10k) %in% Otop5000ALL,]
#Overlapping Group (overlapping genes in a specific feature set)
#OG99H = MAD10k[rownames(MAD10k) %in% Otop99group,] #Too few features
OG250H = MAD10k[rownames(MAD10k) %in% Otop250group,]
OG500H = MAD10k[rownames(MAD10k) %in% Otop500group,]
OG1000H = MAD10k[rownames(MAD10k) %in% Otop1000group,]
OG2000H = MAD10k[rownames(MAD10k) %in% Otop2000group,]
#Global Overlap
GO500H = MAD10k[rownames(MAD10k) %in% COtop500ALL,]
GO1000H = MAD10k[rownames(MAD10k) %in% COtop1000ALL,]
GO2000H = MAD10k[rownames(MAD10k) %in% COtop2000ALL,]
GO2500H = MAD10k[rownames(MAD10k) %in% COtop2500ALL,]
GO5000H = MAD10k[rownames(MAD10k) %in% COtop5000ALL,]


#Check that you can column bind the two platform VST count matrices
table(rownames(OG2000H) == rownames(OG2000N))
table(rownames(GO5000H) == rownames(GO5000N))


#Final Classifier Matrices
CA500 = cbind(CA500N,CA500H)
CA1000 = cbind(CA1000N,CA1000H)
CA2000 = cbind(CA2000N,CA2000H)
CA5000 =  cbind(CA5000N,CA5000H)
#Combined Group (duplicates removed)
CG99 = cbind(CG99N,CG99H)
CG250 = cbind(CG250N,CG250H)
CG500 = cbind(CG500N,CG500H)
CG1000 = cbind(CG1000N,CG1000H)
CG2000 = cbind(CG2000N,CG2000H)
#Overlapping All (overlapping genes in a specific feature set)
#OA250 = cbind(OA250N,OA250H)
OA500 = cbind(OA500N,OA500H)
OA1000 = cbind(OA1000N,OA1000H)
OA2000 = cbind(OA2000N,OA2000H)
OA5000 = cbind(OA5000N,OA5000H)
#Overlapping Group (overlapping genes in a specific feature set)
#OG99 = cbind(OG99N,OG99H)
OG250 = cbind(OG250N,OG250H)
OG500 = cbind(OG500N,OG500H)
OG1000 = cbind(OG1000N,OG1000H)
OG2000 = cbind(OG2000N,OG2000H)
#Global Overlap
GO500 = cbind(GO500N,GO500H)
GO1000 = cbind(GO1000N,GO1000H)
GO2000 = cbind(GO2000N,GO2000H)
GO2500 = cbind(GO2500N,GO2500H)
GO5000 = cbind(GO5000N,GO5000H)


#Manually save each feature set 
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Supervised Classification/451")

filen2 = "Supervised_Classifier_GlobalOverlap_5000.csv"
write.csv(GO5000,filen2) 

#save.image("SAKEfeatures2ClassifierMatrix_Final.RData")

#Load Finalized SAKE Workspace
load("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/GSE153960/3k/SAKEfeatures2ClassifierMatrix_Final.RData")


