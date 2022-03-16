##Use NMF 'feature scores' to determine the subset of MAD genes to use for supervised classification and GSEA

#This script generates the final feature set used in enrichment, network construction, and univariate analysis

#Written By: Jarrett Eshima
#For: Dr. Barbara Smith Lab

## Read in "totfeatures" files. These are obtained from SAKE by downloading all feature scores following NMF clustering.

NumReplicates = 11

#Read in the data
for(i in 1:NumReplicates){
  setwd(paste("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2/NovaSeq/",seq(1,NumReplicates,1),sep="")[i])
  tmp = read.csv('totfeatures.csv') #Generated from SAKE NMF output
  rownames(tmp) = tmp$GeneCard
  tmp = tmp[,-1]
  repname = paste("NovaRep",i,sep="")
  assign(repname,tmp)
}

NumReplicates = 11

for(i in 1:NumReplicates){
  setwd(paste("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2/HiSeq/",seq(1,NumReplicates,1),sep="")[i])
  tmp = read.csv('totfeatures.csv') #Generated from SAKE NMF output
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
Ntop1000ALL = rownames(orderNova)[1:1000]

#Top features overall - HiSeq
Htop1000ALL = rownames(orderHi)[1:1000]


### Combine NovaSeq and HiSeq feature sets

#Combined All
Ctop1000ALL = names(table(c(Ntop1000ALL,Htop1000ALL)))


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

##Parse VST Datasets for Classifier Development (run this section for both sequencing platform cohorts)
##NovaSeq
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2/PreProcessing")
filen = "NOVASEQ_ALS451_VST_Gene_TE_10-4-21.csv" #File generated in previous script
vstcounts = read.csv(filen)
rownames(vstcounts) = vstcounts[,1]
vstcounts = vstcounts[,-1]


#Covert Ensembl rownames to gene symbol rownames (10,000+ names X 10,000+ names = 100 million if statements)
blank = rep(NA,nrow(vstcounts))
vstsymbol = cbind(blank,vstcounts)
vstsymbol[,1] = rownames(vstcounts)

#This loop takes a very long time (6-12 hrs), it is recommended to save the global environment after completion (or vectorize for speed)
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
FullNovaSymbol = finalvstsymbol


##HISEQ
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2/PreProcessing")
filen = "HISEQ_ALS451_VST_Gene_TE_10-4-21.csv"
vstcounts = read.csv(filen) #File generated in previous script
rownames(vstcounts) = vstcounts[,1]
vstcounts = vstcounts[,-1]


#Covert Ensembl rownames to gene symbol rownames
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
FullHiSymbol = finalvstsymbol


#wd = 'C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Enrichment/'
#setwd(wd)
#save.image("SAKEfeatures4Enrichment.RData")
#load("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Enrichment/SAKEfeatures4Enrichment.RData")


##Correct all feature set lists for excel autoformat driven errors (MARCH3 --> 3-Mar)
myfeaturesets = "Ctop1000ALL"
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



##Remove platform-specific sex-dependent genes (~20 total)
#NOVASEQ
MAD10k = FullNovaSymbol
platdif = which(! Ctop1000ALL %in% rownames(MAD10k))
Ctop1000ALL2 = Ctop1000ALL[-platdif]
#HISEQ
MAD10k = FullHiSymbol
h1 = which(! Ctop1000ALL2 %in% rownames(MAD10k))
Ctop1000ALL3 = Ctop1000ALL2[-h1]


##NovaSeq
MAD10k = FullNovaSymbol
CA1000N = MAD10k[rownames(MAD10k) %in% Ctop1000ALL3,]

##Hiseq
MAD10k = FullHiSymbol
CA1000H = MAD10k[rownames(MAD10k) %in% Ctop1000ALL3,]


#Check that you can column bind the two platform VST count matrices
table(rownames(CA1000N) == rownames(CA1000H))


#Final Enrichment Matrix
CA1000 = cbind(CA1000N,CA1000H)


#Manually save feature set 
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Enrichment")

filen2 = "Supervised_Classifier_CombinedPlatform_All_1000.csv"
write.csv(CA1000,filen2) 
