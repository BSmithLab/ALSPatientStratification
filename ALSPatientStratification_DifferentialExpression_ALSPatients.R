#######################################  IDENTIFY SEX-DEPENDENT GENES using DESeq2  #################################################

#Written by: Jarrett Eshima
#Date: April, 2021
#For: Use by the Dr. Barbara Smith Lab at Arizona State University

#Note: This code is not fully automated, as it has been designed to run in individual but related parts.
#Some parts can be skipped, depending on the downstream analysis you are interested in.

## Study Meta Information --
# n = 439 donors/subjects
# nsamples = 1659 RNA-seq files
# ntissues = 11
# GEO Series: GSE153960
# GEO SuperSeries: GSE137810

# n = 473 tissue-matched ALS patient samples
# n = 451 fully processed samples
# n = 140 from GSE124439
# n = 311 from GSE153960 
#######################################################################################################################################

#Load Libraries
library(dplyr)
library(Biobase)
library(limma)
library(DESeq2)
library(biomaRt)

########################################################################################################################################
#############################################  PART 1: LOAD DATA   #####################################################################
########################################################################################################################################
#SPECIFIC Gene Expression Omnibus (GEO) Study - GSE153960
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153960
#https://www.jci.org/articles/view/139741


setwd('C:/Users/jeshima/Documents/Smith Lab/Fall 2020/ALS/GSE153960')
first = read.table(gzfile("GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020.txt.gz"),sep="\t") #Available online at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153960
dim(first)


GSE15_RawCounts = first
GSE15_RawCounts2 = GSE15_RawCounts
rownames(GSE15_RawCounts) = first[,1]
GSE15_RawCounts = GSE15_RawCounts[,-1]

ENS_RawCounts = GSE15_RawCounts

########################################################################################################################################
#####################################  PART 2: CONVERT ENSEMBLs   ######################################################################
########################################################################################################################################
#Convert ENSEMBL IDs to Gene Symbols
#Important Recommendation: Only the UTY and XIST ENSEMBL IDs are needed to complete Part 2
#You can wait until after filtering by median absolute deviation to convert Ensembl IDs to Gene Symbols (much quicker)
#Some small gene naming issues can occur if you use the Gene Symbol matrix generated in this Part of the code

ensembl.genes = c("ENSG00000183878","ENSG00000229807") #UTY, XIST - respectively 
newEG = sub("\\..*","",ensembl.genes) #Keep text to the left of the dot

ensembl_version = "https://dec2016.archive.ensembl.org"
species="human"
ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host=ensembl_version)
gene_positions <- biomaRt::getBM(filters = 'ensembl_gene_id',attributes=c('ensembl_gene_id','hgnc_symbol'), values = newEG, mart = ensembl)


#Can you copy the rownames over? All values must be true
table(sub("\\..*","",rownames(GSE15_RawCounts)) == gene_positions$ensembl_gene_id)


DONE = F
for(i in 1:nrow(GSE15_RawCounts2)){
  for(j in 1:nrow(gene_positions)){
    if(DONE == F){
      if(sub("\\..*","",GSE15_RawCounts2[i,1]) == gene_positions$ensembl_gene_id[j]){
        GSE15_RawCounts2[i,1] = gene_positions$hgnc_symbol[j]
        DONE = T
      }
    }
  }
  DONE = F
  if((i %% 1000) == 0) cat("% Done:",i/nrow(GSE15_RawCounts2)*100,"\n")
}

#Replace missing gene symbol with original ENSEMBL ID
for(i in 1:nrow(GSE15_RawCounts2)){
  if(GSE15_RawCounts2[i,1] == ""){
    GSE15_RawCounts2[i,1] = rownames(GSE15_RawCounts)[i]
  }
}


uniquern = make.names(GSE15_RawCounts2[,1],unique = T)
rownames(GSE15_RawCounts2) = uniquern
GSE15_RawCounts2 = GSE15_RawCounts2[,-1]

#Dataset to move forward with
CleanCountData = GSE15_RawCounts2

#ENSEMBL GENE ID Dataset
ENS_RawCounts

#Order rownames alphabetically
CleanCountData = CleanCountData[order(row.names(CleanCountData)),]

########################################################################################################################################
#####################################  PART 3: ASSIGN SUBJECT SEX   ####################################################################
########################################################################################################################################

### Subset the columns such that only Frontal and Motor postmortem tissue is used (tissue match GSE124439)

#Identify gender using expression of XY Chromosome genes (UTY and XIST)
#UTY - Y Chromosome - Male
#XIST - X Chromosome - Female

MALE = which(rownames(CleanCountData) == "UTY")
tmp = CleanCountData[MALE,]
FEMALE = which(rownames(CleanCountData) == "XIST")
tmp2 = CleanCountData[FEMALE,]

#Visualize cutoff
UTY = as.numeric(tmp)
hist(UTY,25)
abline(v=120,col="red")
XIST = as.numeric(tmp2)
hist(XIST,75)
abline(v=3000,col="red")

#UTY count cutoff: 120
#XIST count cutoff: 3000

SubjectSex = matrix(NA,nrow=ncol(CleanCountData),ncol=1)
rownames(SubjectSex) = colnames(CleanCountData)

for(i in 1:ncol(CleanCountData)){
  if(tmp[i] > 120 && tmp2[i] < 3000){
    SubjectSex[i,1] = "Male"
  }else if(tmp[i] < 120 && tmp2[i] > 3000){
    SubjectSex[i,1] = "Female"
  }
}


#Coldata matrix for DESeq2
setwd('D:/Jarrett/Research/Fall 2020/ALS/GSE153960')
clinicaldata = read.csv("ColData.csv") #Available in Prudencio et al. supplemental files
rownames(clinicaldata) = clinicaldata[,1]
clinicaldata = clinicaldata[,-1]

coldata = clinicaldata
coldata = cbind(coldata,SubjectSex)

#Check DESeq matrix requirements (all lines should be TRUE)
all(rownames(coldata) == colnames(CleanCountData))
ncol(CleanCountData)==nrow(coldata)
CleanCountData = data.matrix(CleanCountData)
all(is.numeric(CleanCountData))


#Save Outputs
#setwd("D:/Jarrett/Research/Summer 2021/SQuIRE Count ALS/FINAL")
#save.image("CleanCountData_StartPoint_10-1-21.RData")
#load("CleanCountData_StartPoint_10-1-21.RData")

load("D:/Jarrett/Research/Summer 2021/SQuIRE Count ALS/FINAL/CleanCountData_StartPoint_10-1-21.RData")

########################################################################################################################################
#####################################  PART 4: FILTER ALS PATIENTS   ###################################################################
########################################################################################################################################

#Inclusion Criteria: 
#1) ALS Spectrum MND or Pre-fALS diagnosis/disease group
#2) Frontal or Motor Cortex postmortem tissue              

#Disease Groups: ALS-TDP, ALS/FTLD, FTLD-TDP, ALS-SOD1, ALS/AD, Control, FTLD-TAU, FTLD-FUS

ALSSubjIndex = which(substr(coldata$disease_group,1,3) == "ALS")
manualadd = c(218) #This subject was labeled ALS in GSE124439 but received a different label in GSE153960
ALSSubjIndex = c(ALSSubjIndex,manualadd)

ALSCountData = CleanCountData[,ALSSubjIndex]
dim(ALSCountData)
ALScoldata = coldata[ALSSubjIndex,]

ENSCountData = ENS_RawCounts[,ALSSubjIndex]
ENScoldata = coldata[ALSSubjIndex,] #Not necessary but including for clarity during analysis


# Restrict dataset further to contain only the GSE124439 matched tissues (frontal and motor cortex)
T1 = which(ALScoldata$tissue == "Frontal Cortex")
T2 = which(ALScoldata$tissue == "Lateral Motor Cortex")
T3 = which(ALScoldata$tissue == "Medial Motor Cortex")
T4 = which(ALScoldata$tissue == "Other Motor Cortex")

TissueIndex = c(T1,T2,T3,T4)
ALSCountData = ALSCountData[,TissueIndex]
dim(ALSCountData)
ALScoldata = ALScoldata[TissueIndex,]
dim(ALScoldata)

ENSCountData = ENSCountData[,TissueIndex]
ENScoldata = ENScoldata[TissueIndex,]

#Check that DESeq2 will work - All outputs must be true
all(rownames(ALScoldata) == colnames(ALSCountData))
ncol(ALSCountData)==nrow(ALScoldata)
ALSCountData = data.matrix(ALSCountData)
all(is.numeric(ALSCountData))

########################################################################################################################################
##############################  PART 5: ADD TRANSPOSABLE ELEMENT COUNTS  ###############################################################
########################################################################################################################################
### Add on TE Counts

##Read in TE File
#The TE Counts file (TECountsF) is generated through SQuIRE_Count_PostProcessing.R and SQuIRE_Count_Meta.R scripts
setwd("D:/Jarrett/Research/Summer 2021/SQuIRE Count ALS/FINAL")
TECountsF = read.csv("TECounts_HGND_9-29-21.csv")
rownames(TECountsF) = TECountsF[,1]
TECountsF = TECountsF[,-1]

## For the Gene Symbol Matrix

#Clean up colnames so that ALScoldata rowname syntax matches TECounts colname syntax 
colnames(TECountsF) = gsub("-",".",colnames(TECountsF))
#Only include the subjects for which TE quantification was successfully completed
ACD = ALSCountData[,colnames(ALSCountData) %in% colnames(TECountsF)]

#Add in rows and columns for TE entries
TEBlank = matrix(NA,nrow=nrow(TECountsF),ncol=ncol(ACD))
rownames(TEBlank) = rownames(TECountsF)
tmpACD = rbind(ACD,TEBlank)
si = nrow(tmpACD) - nrow(TEBlank) +1
ei = nrow(tmpACD)

for(i in 1:ncol(tmpACD)){
  for(j in 1:ncol(TECountsF)){
    if(colnames(tmpACD)[i] == colnames(TECountsF)[j]){
      tmpACD[si:ei,i] = TECountsF[,j]
    }
  }
}
dim(tmpACD)

write.csv(tmpACD,"GSE153960_HGNCGene_Plus_TE_Counts_10-3-21.csv")

#Phenotype Data
tmpAcoldata = ALScoldata[rownames(ALScoldata) %in% colnames(TECountsF),]

## For the ENSEMBL Gene Matrix

#Only include the subjects for which TE quantification was successfully completed
ECD = ENSCountData[,colnames(ENSCountData) %in% colnames(TECountsF)]

#Add in rows and columns for TE entries
TEBlank = matrix(NA,nrow=nrow(TECountsF),ncol=ncol(ECD))
rownames(TEBlank) = rownames(TECountsF)
tmpECD = rbind(ECD,TEBlank)
si = nrow(tmpECD) - nrow(TEBlank) +1
ei = nrow(tmpECD)

for(i in 1:ncol(tmpECD)){
  for(j in 1:ncol(TECountsF)){
    if(colnames(tmpECD)[i] == colnames(TECountsF)[j]){
      tmpECD[si:ei,i] = TECountsF[,j]
    }
  }
}

write.csv(tmpECD,"GSE153960_ENSMBLGene_Plus_TE_Counts_10-4-21.csv")

tmpEcoldata = ENScoldata[rownames(ENScoldata) %in% colnames(TECountsF),]

## Clean up Phenotype Data (use GEO repository metadata file to fill in missing sequencing platform entries)
#Read in MetaData
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping")
Meta = read.csv("GSE153960_MetaData.txt") #Available online at: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA644618&o=acc_s%3Aa

#ALS Patients Only
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/TE Quant/FINAL")
SRRs = read.csv("SRR_ID_List.csv") #File included in Supplemental Tables

index1 = rep(NA,nrow(Meta))
count = 1
for(i in 1:nrow(Meta)){
  for(j in 1:nrow(SRRs)){
    if(Meta$Run[i] == SRRs$SRR.ID[j])
      index1[count] = i
    count = count+1
  }
}

index1 = index1[!is.na(index1)]

MetaData = Meta[index1,] 

CGND_IDs = MetaData$sample_id_alt
CGND_IDs = gsub("-",".",CGND_IDs)
SequencingPlatform_IDs = MetaData$Instrument

Convert = matrix(cbind(CGND_IDs,SequencingPlatform_IDs),ncol=2)

for(i in 1:nrow(Convert)){
  for(j in 1:nrow(tmpAcoldata)){
    if(Convert[i,1] == rownames(tmpAcoldata)[j]){
      tmpAcoldata$sequencing_platform[j] = Convert[i,2]
    }
  }
}

for(i in 1:nrow(Convert)){
  for(j in 1:nrow(tmpEcoldata)){
    if(Convert[i,1] == rownames(tmpEcoldata)[j]){
      tmpEcoldata$sequencing_platform[j] = Convert[i,2]
    }
  }
}

ALScoldata_451_TE = tmpAcoldata
ALSCountData_451_TE = tmpACD

ENScoldata_451_TE = tmpEcoldata
ENSCountData_451_TE = tmpECD


dim(ALScoldata_451_TE);dim(ALSCountData_451_TE)
dim(ENScoldata_451_TE);dim(ENSCountData_451_TE)

#Check - should all be true
table(rownames(ALScoldata_451_TE) == colnames(ALSCountData_451_TE))
table(rownames(ENScoldata_451_TE) == colnames(ENSCountData_451_TE))
table(colnames(ENSCountData_451_TE) == colnames(ALSCountData_451_TE))


#Save Outputs
#setwd("D:/Jarrett/Research/Summer 2021/SQuIRE Count ALS/FINAL")
#save.image("CleanCountData_StartPoint_10-1-21.RData")
#load("CleanCountData_StartPoint_10-1-21.RData")

load("D:/Jarrett/Research/Summer 2021/SQuIRE Count ALS/FINAL/CleanCountData_StartPoint_10-1-21.RData")


########################################################################################################################################
#####################################  PART 6:  Subset Patient Samples  ################################################################
########################################################################################################################################
###Not Necessary but included for convenience of further analysis 

#This code will split the n = 451 patient dataset into the samples analyzed exclusively in GSE153960 (n=331)
#and the samples analyzed in both the reference study GSE124439 and GSE153960 (n=140)

#6 unique ALS samples in GSE122439 (CGND-HRA-): 00296, 00297, 00529, 00569, 00792, 01218
#Remaining 142 ALS samples are used in both studies (2 had missing paired end FASTA files - 140 total)

#Remove duplicate subjects
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/TE Quant/FINAL")
TamSamples = read.csv("TamSubset.csv"); TamSamples = TamSamples[,-1] #This file is included in the Supplemental Tables

subsamp = match(TamSamples,colnames(ALSCountData_451_TE))
subsamp = subsamp[!is.na(subsamp)]

#Check that correct columns are found
colnames(ALSCountData_451_TE)[subsamp]

#For Gene Symbol Matrix

ALSCountData_311_TE = ALSCountData_451_TE[,-subsamp]
ALSCountData_140_TE = ALSCountData_451_TE[,subsamp]

ALScoldata_311_TE = ALScoldata_451_TE[-subsamp,]
ALScoldata_140_TE = ALScoldata_451_TE[subsamp,]

dim(ALSCountData_311_TE); dim(ALSCountData_140_TE)
dim(ALScoldata_311_TE); dim(ALScoldata_140_TE)

#For ENSEMBL Gene Matrix
ENSCountData_311_TE = ENSCountData_451_TE[,-subsamp]
ENSCountData_140_TE = ENSCountData_451_TE[,subsamp]

ENScoldata_311_TE = ENScoldata_451_TE[-subsamp,]
ENScoldata_140_TE = ENScoldata_451_TE[subsamp,]

dim(ENSCountData_311_TE); dim(ENSCountData_140_TE)
dim(ENScoldata_311_TE); dim(ENScoldata_140_TE)

########################################################################################################################################
#####################################  PART 7:  Review Data Frames  ####################################################################
########################################################################################################################################
## Review Datasets

#Gene Symbol - Count Data
dim(ALSCountData_451_TE) #All ALS Subjects
dim(ALSCountData_311_TE) #Only unique ALS subjects from GSE153960
dim(ALSCountData_140_TE) #Only ALS subjects common between GSE124439 and GSE153960 (very similar to Tam et al dataset)

#Clinical Data
dim(ALScoldata_451_TE)
dim(ALScoldata_311_TE)
dim(ALScoldata_140_TE)

#ENSEMBL ID - Count Data
dim(ENSCountData_451_TE)
dim(ENSCountData_311_TE)
dim(ENSCountData_140_TE)

#(Redundant) ENSEMBL Clinical Data
dim(ENScoldata_451_TE)
dim(ENScoldata_311_TE)
dim(ENScoldata_140_TE)

########################################################################################################################################
##############################  PART 8: DESeq2 Dependent Gene Removal  #################################################################
########################################################################################################################################
library(DESeq2)

setwd('C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2')

################### FOR DATASETS CONTAINING MORE THAN ONE SEQUENCING PLATFORM ##################################

#Over 17k genes are significantly differentially expressed when considering sequencing platform
#To avoid this platform-dependent bias, patients are separated into two groups based on sequencing platform
#Authors of GSE153960 took a similar approach to their analysis.

for(i in 1:nrow(ENScoldata_451_TE)){
  if(ENScoldata_451_TE$sequencing_platform[i] == "Illumina HiSeq 2500"){
    ENScoldata_451_TE$sequencing_platform[i] = "HiSeq"
  }else{
    ENScoldata_451_TE$sequencing_platform[i] = "NovaSeq"
  }
}

catt = table(ENScoldata_451_TE$sequencing_platform)
catt

HiSeqCounts = data.frame(matrix(NA,nrow(ENSCountData_451_TE),ncol=catt[[1]]))
NovaSeqCounts = data.frame(matrix(NA,nrow(ENSCountData_451_TE),ncol=catt[[2]]))
HiSeqCol = data.frame(matrix(NA,nrow = catt[[1]],ncol(ENScoldata_451_TE)))
NovaSeqCol = data.frame(matrix(NA,nrow = catt[[2]],ncol(ENScoldata_451_TE)))
count=count1=1
for(i in 1:nrow(ENScoldata_451_TE)){
  if(ENScoldata_451_TE$sequencing_platform[i] == "HiSeq"){
    HiSeqCounts[,count] = ENSCountData_451_TE[,i]
    HiSeqCol[count,] = ENScoldata_451_TE[i,]
    count = count+1
  }else{
    NovaSeqCounts[,count1] = ENSCountData_451_TE[,i]
    NovaSeqCol[count1,] = ENScoldata_451_TE[i,]
    count1=count1+1
  }
}

rownames(HiSeqCounts) = rownames(NovaSeqCounts) = rownames(ENSCountData_451_TE)
HiSeqi = which(ENScoldata_451_TE$sequencing_platform == "HiSeq")
NovaSeqi = which(ENScoldata_451_TE$sequencing_platform == "NovaSeq")
hnames = colnames(ENSCountData_451_TE)[HiSeqi]
nnames = colnames(ENSCountData_451_TE)[NovaSeqi]
colnames(HiSeqCounts) = hnames
colnames(NovaSeqCounts) = nnames

rownames(HiSeqCol) = hnames
rownames(NovaSeqCol) = nnames
colnames(HiSeqCol) = colnames(NovaSeqCol) = colnames(ENScoldata_451_TE)

table(rownames(HiSeqCol) == colnames(HiSeqCounts))
table(rownames(NovaSeqCol) == colnames(NovaSeqCounts))

###############################################################################################################
#Save the platform-separated raw count files for downstream univariate analysis

HSfn = "HiSeq_RawCounts_ENSG_ALS451.csv"
write.csv(HiSeqCounts,HSfn)
NSfn = "NovaSeq_RawCounts_ENSG_ALS451.csv"
write.csv(NovaSeqCounts,NSfn)

###############################################################################################################
#Sex-dependent genes identified after separating subjects based on sequencing platform

###HiSeq Cohort

#Rounding non-integer counts recommended by RSEM authors for DESeq differential expression
rALSCountData = round(HiSeqCounts,0)

dds = DESeqDataSetFromMatrix(countData = rALSCountData, colData = HiSeqCol, design= ~ SubjectSex, tidy=F)
dds$SubjectSex = relevel(dds$SubjectSex,ref = "Male")
dseq = DESeq(dds,betaPrior=T)
res = results(dseq)
sig = res[! is.na(res$padj) & res$padj<0.05,]
write.csv(sig,"HISEQ_ALS451_Subj_SexSigGenes_10-4-21.csv")

vsd = varianceStabilizingTransformation(dseq)
vstcounts = assay(vsd)
filtdata = vstcounts[! (rownames(vstcounts) %in% rownames(sig)),]

file = "HISEQ_ALS451_VST_Gene_TE_10-4-21.csv"
#write.table(filtdata,file,sep="\t",quote=F,row.names=T,col.names=T)
write.csv(filtdata,file)

madval = apply(as.matrix(filtdata),1,mad)
madval = madval[order(madval,decreasing=T)]

#Important Note: There are three "duplicated" gene names in the HiSEQ MAD10k gene set
#Two are genuine duplicates and are removed in ALSPatientStratification_ClusterPrep.R
mad10000 = filtdata[match(names(madval)[1:10000],rownames(filtdata)),] #This gives slightly less than 10k genes
file = "HISEQ_451_ALS_MAD10k_GeneENSEMBL_TE_10-4-21.csv"
write.table(mad10000,file,sep="\t",quote=F,row.names=T,col.names=T)
write.csv(mad10000,file)

#For duplicate rownames (<10 total) - grab more than you need an manually adjust in Excel
mad10050 = filtdata[match(names(madval)[1:10050],rownames(filtdata)),]
file = "HISEQ_451_ALS_MAD10k50_GeneENSEMBL_TE_10-4-21.csv"
write.csv(mad10050,file)

#For convenience, code for top 5k genes by mad is also included
# mad5000 = filtdata[match(names(madval)[1:5000],rownames(filtdata)),]
# file = "HISEQ_451_ALS_MAD5k_Gene_TE.txt"
# write.table(mad5000,file,sep="\t",quote=F,row.names=T,col.names=T)

###############################################################################################################

###NovaSeq Cohort
rENSCountData = round(NovaSeqCounts,0)
ddsens = DESeqDataSetFromMatrix(countData = rENSCountData, colData = NovaSeqCol, design= ~ SubjectSex, tidy=F)
ddsens$SubjectSex = relevel(ddsens$SubjectSex,ref = "Male")
dseqens = DESeq(ddsens,betaPrior=T)
resens = results(dseqens)
sigens = resens[! is.na(resens$padj) & resens$padj<0.05,]
write.csv(sigens,"NOVASEQ_ALS451_Subj_SexSigGenes_10-4-21.csv")

vsdens = varianceStabilizingTransformation(dseqens)
vstcountsens = assay(vsdens)
fData.ens = vstcountsens[! (rownames(vstcountsens) %in% rownames(sigens)),]

file = "NOVASEQ_ALS451_VST_Gene_TE_10-4-21.csv"
#write.table(fData.ens,file,sep="\t",quote=F,row.names=T,col.names=T)
write.csv(fData.ens,file)

madVal.ens = apply(as.matrix(fData.ens),1,mad)
madVal.ens = madVal.ens[order(madVal.ens,decreasing=T)]

# mad5000.ens = fData.ens[match(names(madVal.ens)[1:5000],rownames(fData.ens)),]
# file = "NOVASEQ_451_ALS_MAD5k_GeneENSEMBL_TE_8-17-21.txt"
# write.table(mad5000.ens,file,sep="\t",quote=F,row.names=T,col.names=T)

#Important Note: There are 14 "duplicated" gene names in the NovaSEQ MAD10k gene set
#13 are genuine duplicates and are removed in ALSPatientStratification_ClusterPrep.R
mad10000.ens = fData.ens[match(names(madVal.ens)[1:10000],rownames(fData.ens)),] #This gives slightly less than 10k genes
file = "NOVASEQ_451_ALS_MAD10k_GeneENSEMBL_TE_10-4-21.txt"
write.table(mad10000.ens,file,sep="\t",quote=F,row.names=T,col.names=T)
write.csv(mad10000.ens,file)

#For duplicate rownames (<15 total) - grab more than you need
mad10050.ens = fData.ens[match(names(madVal.ens)[1:10050],rownames(fData.ens)),]
file = "NOVASEQ_451_ALS_MAD10k50_GeneENSEMBL_TE_10-4-21.csv"
write.csv(mad10050,file)


########################################################################################################################################
##############################  PART 9: Convert ENSEMBL IDs  ###########################################################################
########################################################################################################################################

#Convenient code to convert MAD 5k and MAD 10k ENSEMBLs to Symbols (much quicker than converting all genes up front)
#It is recommended to convert ENSEMBL IDs to Gene Symbols for SAKE unsupervised clustering, however it is not necessary

egenes = rownames(mad10050.ens)
newEG = sub("\\..*","",egenes) #Keep text to the left of the dot

ensembl_version = "https://dec2016.archive.ensembl.org" #Updated: https://jan2020.archive.ensembl.org
species="human"
ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host=ensembl_version)
gene_positions <- biomaRt::getBM(filters = 'ensembl_gene_id',attributes=c('ensembl_gene_id','hgnc_symbol'), values = newEG, mart = ensembl)


#Covert Ensembl rownames to gene symbol rownames (10,000+ names X 10,000+ names = 100 million if statements)
blank = rep(NA,nrow(mad10050.ens))
mad10000.ens1 = cbind(blank,mad10050.ens)
mad10000.ens1[,1] = rownames(mad10050.ens)

DONE = F
for(i in 1:nrow(mad10000.ens1)){
  for(j in 1:nrow(gene_positions)){
    if(DONE == F){
      if(sub("\\..*","",mad10000.ens1[i,1]) == gene_positions$ensembl_gene_id[j]){
        mad10000.ens1[i,1] = gene_positions$hgnc_symbol[j]
        DONE = T
      }
    }
  }
  DONE = F
  if((i %% 100) == 0) cat("% Done:",i/nrow(mad10000.ens1)*100,"\n")
}

#Replace missing gene symbol with original ENSEMBL ID
for(i in 1:nrow(mad10000.ens1)){
  if(mad10000.ens1[i,1] == ""){
    mad10000.ens1[i,1] = rownames(mad10000.ens1)[i]
  }
}

finalmad10000 = mad10000.ens1
rownames(finalmad10000) = finalmad10000[,1]
finalmad10000 = finalmad10000[,-1]

file = "NOVASEQ_ALS451_MAD10k50_10-4-21_SYMBOL.csv"
write.csv(finalmad10000,file)
file = "NOVASEQ_ALS451_MAD10k50_10-4-21_SYMBOL.txt"
write.table(finalmad10000,file,sep="\t",quote=F,row.names=T,col.names=T)


#Save Outputs
#setwd("D:/Jarrett/Research/Summer 2021/SQuIRE Count ALS/FINAL")
#save.image("DESeq2Completed_10-4-21.RData")
load("D:/Jarrett/Research/Summer 2021/SQuIRE Count ALS/FINAL/DESeq2Completed_10-4-21.RData")
