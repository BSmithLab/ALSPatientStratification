#######################################  IDENTIFY-DEPENDENT GENES using DESeq2  #################################################

#Written by: Jarrett Eshima
#Date: April, 2021
#For: Use by the Dr. Barbara Smith Lab at Arizona State University

## Description: This code is designed to read in control subject data and filter the count matrix according to important genes identified through 
# unsupervised clustering and supervised classification. Those steps should be completed before running this code.
# Processed control subject count information is utilized in the univariate analysis and gene set enrichment steps.

#######################################################################################################################################

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

#New Working Directory
setwd('C:/Users/jeshima/Documents/Smith Lab/Fall 2020/ALS/GSE153960')
setwd('D:/Jarrett/Research/Fall 2020/ALS/GSE153960')
first = read.table(gzfile("GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020.txt.gz"),sep="\t") #Available online at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153960
dim(first)

GSE15_RawCounts = first
rownames(GSE15_RawCounts) = first[,1]
GSE15_RawCounts = GSE15_RawCounts[,-1]
GSE15_RawCounts[1:5,1:10]


#Read in control names
setwd('C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/Controls')
OMNDlist = read.table("Control_SRR_OtherNeurologicalDisorders_Final83_commasep.txt",sep=",",header = F) #File included in Supplemental Tables
HClist = read.table("Control_SRR_NonNeurologicalControl_Final95_commasep.txt",sep=",",header = F) #File included in Supplemental Tables

OMNDlist = unlist(list(OMNDlist))
HClist = unlist(list(HClist))

#Read in MetaData file to convert SRR IDs to CGND-HRA IDs
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping")
Meta = read.csv("GSE153960_MetaData.txt") #Available online at: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA644618&o=acc_s%3Aa

OMND_CGND = rep(NA,length(OMNDlist))
for(i in 1:length(OMNDlist)){
  for(j in 1:nrow(Meta)){
    if(OMNDlist[i] == Meta$Run[j]){
      OMND_CGND[i] = Meta$sample_id_alt[j]
    }
  }
}
OMND_CGND = OMND_CGND[!is.na(OMND_CGND)]

HC_CGND = rep(NA,length(HClist))
for(i in 1:length(HClist)){
  for(j in 1:nrow(Meta)){
    if(HClist[i] == Meta$Run[j]){
      HC_CGND[i] = Meta$sample_id_alt[j]
    }
  }
}
HC_CGND = HC_CGND[!is.na(HC_CGND)]

OMND_subj = OMND_CGND
OMND_subj = gsub("-",".",OMND_subj)
HC_subj = HC_CGND
HC_subj = gsub("-",".",HC_subj)

ALSCountData_451_TE = read.csv("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/ALS_FullCountMatrix.csv") #This file is provided in the supplemental tables

ALS_subj = colnames(ALSCountData_451_TE)[-1]
Control_subj = c(OMND_subj,HC_subj)
All_subj = c(ALS_subj,Control_subj)

GSE153960_AllSubjects = GSE15_RawCounts[,colnames(GSE15_RawCounts) %in% All_subj]
GSE153960_Controls = GSE15_RawCounts[,colnames(GSE15_RawCounts) %in% Control_subj]

#There are sample IDs included in the meta data file that were not run / deposited into the repository (pandemic limitations)
missing_controls = All_subj[! All_subj %in% colnames(GSE15_RawCounts)]

#Missing HC
missing_HC = missing_controls[missing_controls %in% HC_subj]
#Missing OMND
missing_OMND = missing_controls[missing_controls %in% OMND_subj]

#Save points - this script can take awhile to get through
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS")
save.image("GSE153960_ALS_plus_Controls.RData")
#load("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/GSE153960_ALS_plus_Controls.RData")

#########################################################################################################################################################
# Convert ENSMBL ID to HGNC SYMBOL to identify control subject sex
# Build the phenotype dataframe
#########################################################################################################################################################


egenes = c("ENSG00000183878","ENSG00000229807") #UTY, XIST - respectively 
ensembl_version = "https://dec2016.archive.ensembl.org" 
species="human"
ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host=ensembl_version)
gene_positions <- biomaRt::getBM(filters = 'ensembl_gene_id',attributes=c('ensembl_gene_id','hgnc_symbol'), values = egenes, mart = ensembl)

GSE153960_Controls_SYM = GSE153960_Controls
blank = rep(NA,nrow(GSE153960_Controls_SYM))
GSE153960_Controls_SYM = cbind(blank,GSE153960_Controls_SYM)
GSE153960_Controls_SYM[,1] = rownames(GSE153960_Controls_SYM)

DONE = F
for(i in 1:nrow(GSE153960_Controls_SYM)){
  for(j in 1:nrow(gene_positions)){
    if(DONE == F){
      if(sub("\\..*","",GSE153960_Controls_SYM[i,1]) == gene_positions$ensembl_gene_id[j]){
        GSE153960_Controls_SYM[i,1] = gene_positions$hgnc_symbol[j]
        DONE = T
      }
    }
  }
  DONE = F
  if((i %% 1000) == 0) cat("% Done:",i/nrow(GSE153960_Controls_SYM)*100,"\n")
}

#Replace missing gene symbol with original ENSEMBL ID
for(i in 1:nrow(GSE153960_Controls_SYM)){
  if(GSE153960_Controls_SYM[i,1] == ""){
    GSE153960_Controls_SYM[i,1] = rownames(GSE153960_Controls_SYM)[i]
  }
}

uniquern = make.names(GSE153960_Controls_SYM[,1],unique = T)
rownames(GSE153960_Controls_SYM) = uniquern
GSE153960_Controls_SYM = GSE153960_Controls_SYM[,-1]

#Dataset to move forward with
CleanCountData_Controls = GSE153960_Controls_SYM

#Order rownames alphabetically
CleanCountData_Controls = CleanCountData_Controls[order(row.names(CleanCountData_Controls)),]

#Save points - this script can take awhile to get through
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS")
save.image("GSE153960_ALS_plus_Controls.RData")
#load("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/GSE153960_ALS_plus_Controls.RData")


######### Subset the columns such that only Frontal and Motor postmortem tissue is used (tissue match GSE124439)

#Identify gender using expression of XY Chromosome genes (UTY and XIST)
#UTY - Y Chromosome - Male marker
#XIST - X Chromosome - Female marker

MALE = which(rownames(CleanCountData_Controls) == "UTY")
tmp = CleanCountData_Controls[MALE,]
FEMALE = which(rownames(CleanCountData_Controls) == "XIST")
tmp2 = CleanCountData_Controls[FEMALE,]

#Visualize
UTY = as.numeric(tmp)
hist(UTY,25)
abline(v=120,col="red") #count cutoff
XIST = as.numeric(tmp2)
hist(XIST,75)
abline(v=3000,col="red") #count cutoff

#Loop through subjects (columns) to determine sex
#Note that the ColData matrix and Raw Count Matrix have subjects in the same order

SubjectSex = matrix(NA,nrow=ncol(CleanCountData_Controls),ncol=1)
rownames(SubjectSex) = colnames(CleanCountData_Controls)

for(i in 1:ncol(GSE153960_Controls)){
  if(tmp[i] > 100 && tmp2[i] < 3000){
    SubjectSex[i,1] = "Male"
  }else if(tmp[i] < 100 && tmp2[i] > 3000){
    SubjectSex[i,1] = "Female"
  }
}

#Coldata matrix for DESeq2
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS")
clinicaldata = read.csv("ColData.csv") #Available in Prudencio et al. supplemental files
rownames(clinicaldata) = clinicaldata[,1]
clinicaldata = clinicaldata[,-1]

coldata = clinicaldata

#Filter coldata for only the control subjects meeting the inclusion criteria
finalcoldata = coldata[rownames(coldata) %in% rownames(SubjectSex),]

finalcoldata = cbind(finalcoldata,SubjectSex)

#Check DESeq matrix requirements (all lines should be TRUE)
all(rownames(finalcoldata) == colnames(CleanCountData_Controls))
ncol(CleanCountData_Controls)==nrow(finalcoldata)
CleanCountData_Controls = data.matrix(CleanCountData_Controls)
all(is.numeric(CleanCountData_Controls))

#write.csv(finalcoldata,"Control_coldata_SUBTYPES_10-21-21.csv")

#Save Outputs
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS")
save.image("GSE153960_ALS_plus_Controls.RData")
#load("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/GSE153960_ALS_plus_Controls.RData")


#########################################################################################################################################################
# Use ALS TE count matrix to clean up control TE files
#########################################################################################################################################################

TEKEEP = read.csv("ALS_TECounts.csv") #File included in Supplemental Tables
rownames(TEKEEP) = TEKEEP[,1]
TEKEEP = TEKEEP[,-1]

TEKEEP2 = rownames(TEKEEP)

## Healthy Controls
HealthyControlCounts = CleanCountData_Controls[,colnames(CleanCountData_Controls) %in% HC_subj]
## Other Neurological Disease Controls
OMNDControlCounts = CleanCountData_Controls[,colnames(CleanCountData_Controls) %in% OMND_subj]

#Generate containers for TE counts
HC_TECounts = matrix(NA,nrow=length(TEKEEP2),ncol=ncol(HealthyControlCounts))
rownames(HC_TECounts) = TEKEEP2
colnames(HC_TECounts) = colnames(HealthyControlCounts)

OND_TECounts = matrix(NA,nrow=length(TEKEEP2),ncol=ncol(OMNDControlCounts))
rownames(OND_TECounts) = TEKEEP2
colnames(OND_TECounts) = colnames(OMNDControlCounts)


#Run lines 240-380 for both the Healthy Control TE counts and Other Neurological Disease TE counts (Only necessary if TE count .txt files are stored in separate folders)
wd = "D:/Jarrett/Research/Summer 2021/SQuIRE Count OMND"
setwd(wd)
allfiles = list.files(path = wd,pattern = ".txt")
allfilelength = nchar(allfiles)

#Identify substring handle
handle = "SRR"
subfile = sub(handle,"",allfiles)
accessions = sub("\\_.*","",subfile)
clean_accessions = names(table(accessions))


#Get character lengths for substr arguments
AL = nchar(accessions)
HL = as.numeric(nchar(handle))

#Loop through files and grab the index to organize
count1=count2=count3=count4 = 1
TEfiles = subFfiles = refGenefiles = abundfiles = rep(NA,length(accessions))
for(i in 1:length(allfiles)){
  if(substr(allfiles[i],HL+AL[i]+1,allfilelength[i]) == "_TEcounts.txt"){
    TEfiles[count1] = i
    count1 = count1+1
  }else if(substr(allfiles[i],HL+AL[i]+1,allfilelength[i]) == "_subFcounts.txt"){
    subFfiles[count2] = i
    count2 = count2+1
  }else if(substr(allfiles[i],HL+AL[i]+1,allfilelength[i]) == "_refGenecounts.txt"){
    refGenefiles[count3] = i
    count3 = count3+1
  }else if(substr(allfiles[i],HL+AL[i]+1,allfilelength[i]) == "_abund.txt"){
    abundfiles[count4] = i
    count4 = count4+1
  }
}

#Remove NA values
TEI = TEfiles[!is.na(TEfiles)]
subFI = subFfiles[!is.na(subFfiles)]
refGeneI = refGenefiles[!is.na(refGenefiles)]
abundI = abundfiles[!is.na(abundfiles)]


#Read in Subjects - This can take awhile depending on the file number and file size
#It is not recommended to load more than 100 TE count files unless your computer has a large RAM resource

pseudoc = 1
for(i in 1:length(TEI)){
  nam = paste("HC",clean_accessions[i],sep = "") #You can change the object handle for OND and HC groups but not necessary due to unique accessions
  index = TEI[i]
  assign(nam,read.table(allfiles[index],header = T))
  if((pseudoc %% 5) == 0) cat("TE Read-In %:",pseudoc/(length(TEI))*100,"\n")
  pseudoc = pseudoc+1
}


for(i in 1:length(TEI)){
  tmp = get(paste("HC",clean_accessions[i],sep = ""))

  tmp2 = tmp[tmp$TE_ID %in% TEKEEP2,]
  assign(paste("HC",clean_accessions[i],sep = ""),tmp2)
  if((i %% 1) == 0) cat("Subject Completed:",i,"\n")
}


#Filter the unnecessary columns
for(i in 1:length(TEI)){
  tmp = get(paste("HC",clean_accessions[i],sep = ""))
  tmpn = tmp$TE_ID
  tmpc = tmp$tot_counts
  filt = matrix(c(tmpn,tmpc),ncol=2)
  colnames(filt) = c("TE_ID","Counts")
  filt = data.frame(filt)
  assign(paste("HC",clean_accessions[i],sep = ""),filt)
  if((i %% 1) == 0) cat("Subject Completed:",i,"\n")
}


#Give zero counts for the missing TEs
thresh = length(TEKEEP2)

for(i in 1:length(TEI)){
  tmp = get(paste("HC",clean_accessions[i],sep = ""))

  if(nrow(tmp) < thresh){

    missingTEs = TEKEEP2[! TEKEEP2 %in% tmp$TE_ID ]
    myrep = rep(0,length(missingTEs))
    tmp2 = c(missingTEs,myrep)
    missingTEs2 = matrix(tmp2,nrow=length(missingTEs),ncol=2)
    colnames(missingTEs2) = colnames(tmp)

    tmp = rbind(tmp,missingTEs2)
    assign(paste("HC",clean_accessions[i],sep = ""),tmp)
    if((i %% 1) == 0) cat("Subject Completed:",i,"\n")

  }

}


for(i in 1:length(TEI)){
  tmp = get(paste("HC",clean_accessions[i],sep = ""))
  
  srr = paste("SRR",clean_accessions[i],sep="")
  
  for(j in 1:nrow(Meta)){
    if(srr == Meta$Run[j]){
      CGND = Meta$sample_id_alt[j]
    }
  }
  
  blank = rep(CGND,nrow(tmp))
  
  tmp2 = cbind(tmp,blank)
  colnames(tmp2) = c("TE_ID","Counts","CGND")
  
  assign(paste("HC",clean_accessions[i],sep = ""),tmp2)
  
}

## FOR HC Subjects
for(i in 1:length(TEI)){
  tmp = get(paste("HC",clean_accessions[i],sep = ""))
  tmp$CGND = gsub("-",".",tmp$CGND)

  tmpindex = which(colnames(HC_TECounts) == tmp$CGND[1])

  if(length(tmpindex)>0){
  for(j in 1:nrow(tmp)){
    for(k in 1:nrow(HC_TECounts)){

      if(tmp$TE_ID[j] == rownames(HC_TECounts)[k]){

        HC_TECounts[k,tmpindex] = tmp$Counts[j]

      }

    }
  }
  }
  if((i %% 1) == 0) cat("Subject Completed:",i,"\n")
}

## FOR OND Subjects
for(i in 1:length(TEI)){
  tmp = get(paste("HC",clean_accessions[i],sep = ""))
  tmp$CGND = gsub("-",".",tmp$CGND)
  
  tmpindex = which(colnames(OND_TECounts) == tmp$CGND[1])
  
  if(length(tmpindex)>0){
    for(j in 1:nrow(tmp)){
      for(k in 1:nrow(OND_TECounts)){
      
      if(tmp$TE_ID[j] == rownames(OND_TECounts)[k]){
        
        OND_TECounts[k,tmpindex] = tmp$Counts[j]
        
      }
      
      }
    }
  }
  if((i %% 1) == 0) cat("Subject Completed:",i,"\n")
}

write.csv(OND_TECounts,"OND_TECounts.csv")


#Can you use rbind? Must all be true
table(colnames(OND_TECounts) == colnames(OMNDControlCounts))
table(colnames(HC_TECounts) == colnames(HealthyControlCounts))


HC_TECountscopy = HC_TECounts
OND_TECountscopy = OND_TECounts

#Convert from character to numeric
HC_TECountscopy = matrix(as.numeric(HC_TECountscopy),nrow=nrow(HC_TECountscopy),ncol=ncol(HC_TECountscopy))
OND_TECountscopy = matrix(as.numeric(OND_TECountscopy),nrow=nrow(OND_TECountscopy),ncol=ncol(OND_TECountscopy))
colnames(HC_TECountscopy) = colnames(HC_TECounts); rownames(HC_TECountscopy) = rownames(HC_TECounts)
colnames(OND_TECountscopy) = colnames(OND_TECounts); rownames(OND_TECountscopy) = rownames(OND_TECounts)

FinalHealthyControlCounts = rbind(HealthyControlCounts,HC_TECountscopy)
FinalONDControlCounts = rbind(OMNDControlCounts,OND_TECountscopy)

dim(HealthyControlCounts);dim(HC_TECountscopy);dim(FinalHealthyControlCounts)
dim(OMNDControlCounts);dim(OND_TECountscopy);dim(FinalONDControlCounts)

#Save count matrices
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/Controls")
write.csv(FinalHealthyControlCounts,"HealthyControl_GeneTEcounts.csv")
write.csv(FinalONDControlCounts,"OMNDControl_GeneTEcounts.csv")

table(rownames(FinalHealthyControlCounts) == rownames(FinalONDControlCounts))
FinalControlCounts = cbind(FinalHealthyControlCounts,FinalONDControlCounts)
write.csv(FinalControlCounts,"Control_GeneTEcounts.csv")


#Parse Finalcoldata into HC and OND
HC_coldata = finalcoldata[rownames(finalcoldata) %in% colnames(FinalHealthyControlCounts),]
OND_coldata = finalcoldata[rownames(finalcoldata) %in% colnames(FinalONDControlCounts),]
Control_coldata = finalcoldata[rownames(finalcoldata) %in% colnames(FinalControlCounts),]

#Check - should all be true
table(rownames(HC_coldata) == colnames(FinalHealthyControlCounts))
table(rownames(OND_coldata) == colnames(FinalONDControlCounts))

#Save Outputs
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS")
#save.image("GSE153960_ALS_plus_Controls_v3.RData")
load("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/Controls/GSE153960_ALS_plus_Controls_v3.RData")

###############################################################################################

dim(FinalControlCounts)   
dim(FinalHealthyControlCounts)
dim(FinalONDControlCounts) 

dim(finalcoldata)
dim(HC_coldata)
dim(OND_coldata)

ind = nrow(finalcoldata)
finalcoldata$sequencing_platform[ind] = "NovaSeq" #Obtained from Prudencio et al supplemental files
ind = nrow(HC_coldata)
HC_coldata$sequencing_platform[ind] = "NovaSeq"
###############################################################################################

########################################################################################################################################
#DESeq2 Dependent Gene Removal 
########################################################################################################################################
library(DESeq2)

setwd('C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/Controls')

#Rounding non-integer counts recommended by RSEM authors for DESeq differential expression
#OND cohort run exclusively on NovaSeq platform
rENSCountData = round(FinalControlCounts,0)

################### FOR DATASETS CONTAINING MORE THAN ONE SEQUENCING PLATFORM ##################################

#Separate control cohort by sequencing platform - same as ALS cohort

catt = table(finalcoldata$sequencing_platform)

HiSeqCounts = data.frame(matrix(NA,nrow(FinalControlCounts),ncol=catt[[1]]))
NovaSeqCounts = data.frame(matrix(NA,nrow(FinalControlCounts),ncol=catt[[1]]))
HiSeqCol = data.frame(matrix(NA,nrow = catt[[1]],ncol(finalcoldata)))
NovaSeqCol = data.frame(matrix(NA,nrow = catt[[1]],ncol(finalcoldata)))
count=count1=1
for(i in 1:nrow(finalcoldata)){
  if(finalcoldata$sequencing_platform[i] == "HiSeq 2500"){
    HiSeqCounts[,count] = FinalControlCounts[,i]
    HiSeqCol[count,] = finalcoldata[i,]
    rownames(HiSeqCol)[count] = rownames(finalcoldata)[i]
    count = count+1
  }else{
    NovaSeqCounts[,count1] = FinalControlCounts[,i]
    NovaSeqCol[count1,] = finalcoldata[i,]
    rownames(NovaSeqCol)[count1] = rownames(finalcoldata)[i]
    count1=count1+1
  }
}

rownames(HiSeqCounts) = rownames(NovaSeqCounts) = rownames(FinalControlCounts)
HiSeqi = which(finalcoldata$sequencing_platform == "HiSeq 2500")
NovaSeqi = which(finalcoldata$sequencing_platform == "NovaSeq")
hnames = rownames(finalcoldata)[HiSeqi]
nnames = rownames(finalcoldata)[NovaSeqi]
colnames(HiSeqCounts) = hnames
colnames(NovaSeqCounts) = nnames

colnames(HiSeqCol) = colnames(NovaSeqCol) = colnames(finalcoldata)

table(colnames(NovaSeqCounts) == rownames(NovaSeqCol)) #quick check for DESeq2
table(colnames(HiSeqCounts) == rownames(HiSeqCol)) #quick check for DESeq2

###############################################################################################################
#Sex-dependent genes identified after separating subjects based on sequencing platform
#Because the supervised classifier feature set is used to filter the control subjects, removing sex-dependent genes is not necessary

#HiSeq
rENSCountData = round(HiSeqCounts,0)
ddsens = DESeqDataSetFromMatrix(countData = rENSCountData, colData = HiSeqCol, design= ~ SubjectSex, tidy=F)
ddsens$SubjectSex = relevel(ddsens$SubjectSex,ref = "Male")
dseqens = DESeq(ddsens,betaPrior=T)
resens = results(dseqens)
sigens = resens[! is.na(resens$padj) & resens$padj<0.05,]
write.csv(sigens,"HISEQ_AllControl_SexSigGenes_ENSG_10-25-21.csv")

vsdens = varianceStabilizingTransformation(dseqens)
vstcountsens = assay(vsdens)
fData.ens = vstcountsens[! (rownames(vstcountsens) %in% rownames(sigens)),]

file = "HISEQ_AllControl_VST_ENSG_10-25-21.csv"
write.csv(vstcountsens,file)
# You can save the VST count matrix with sex-dependent genes removed, however this may cause issues depending on the final feature set selected
# file = "HISEQ_HealthyControl_VST_Symbol_SexDepGeneRemoval_10-25-21.csv"
# write.csv(fData.ens,file)


#NovaSeq
rENSCountData = round(NovaSeqCounts,0)
ddsens = DESeqDataSetFromMatrix(countData = rENSCountData, colData = NovaSeqCol, design= ~ SubjectSex, tidy=F)
ddsens$SubjectSex = relevel(ddsens$SubjectSex,ref = "Male")
dseqens = DESeq(ddsens,betaPrior=T)
resens = results(dseqens)
sigens = resens[! is.na(resens$padj) & resens$padj<0.05,]
write.csv(sigens,"NOVASEQ_AllControl_SexSigGenes_ENSG_10-25-21.csv")

vsdens = varianceStabilizingTransformation(dseqens)
vstcountsens = assay(vsdens)
fData.ens = vstcountsens[! (rownames(vstcountsens) %in% rownames(sigens)),]

file = "NOVASEQ_AllControl_VST_ENSG_10-25-21.csv"
write.csv(vstcountsens,file)

# file = "NOVASEQ_OND_VST_Symbol_SexDepGeneRemoval_10-25-21.csv"
# write.csv(fData.ens,file)


###############################################################################################################
#Build preliminary count matrix for GSEA (run GSEA_prep.R after this to generate final GSEA compatible count file)

ClassifierMatrix = read.csv("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/GSEA/Supervised_Classifier_CombinedPlatform_All_1000.csv") #File generated in previous script

rownames(ClassifierMatrix) = ClassifierMatrix[,1]
ClassifierMatrix = ClassifierMatrix[,-1]

GSEA_genes = rownames(ClassifierMatrix)

HISEQ_VST = read.csv("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/Controls/HISEQ_AllControl_VST_ENSG_10-25-21.csv") #File generated earlier in this script
NOVASEQ_VST = read.csv("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/Controls/NOVASEQ_AllControl_VST_ENSG_10-25-21.csv") #File generated earlier in this script
rownames(HISEQ_VST) = HISEQ_VST[,1]
rownames(NOVASEQ_VST) = NOVASEQ_VST[,1]
HISEQ_VST = HISEQ_VST[,-1]
NOVASEQ_VST = NOVASEQ_VST[,-1]

##Convert ENSEMBL ID to Gene Symbol for the feature subset selected

egenes = GSEA_genes
ensembl_version = "https://dec2016.archive.ensembl.org" 
species="human"
ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host=ensembl_version)
gene_positions <- biomaRt::getBM(filters = 'hgnc_symbol',attributes=c('ensembl_gene_id','hgnc_symbol'), values = egenes, mart = ensembl)

table(substr(egenes,1,4) == "ENSG")
table(substr(egenes,1,3) == "chr")

NOVASEQ_VST_SYM = NOVASEQ_VST
blank = rep(NA,nrow(NOVASEQ_VST_SYM))
NOVASEQ_VST_SYM = cbind(blank,NOVASEQ_VST_SYM)
NOVASEQ_VST_SYM[,1] = rownames(NOVASEQ_VST_SYM)

#This loop takes approximately 20-30 minutes 
DONE = F
for(i in 1:nrow(NOVASEQ_VST_SYM)){
  for(j in 1:nrow(gene_positions)){
    if(DONE == F){
      if(sub("\\..*","",NOVASEQ_VST_SYM[i,1]) == gene_positions$ensembl_gene_id[j]){
        NOVASEQ_VST_SYM[i,1] = gene_positions$hgnc_symbol[j]
        DONE = T
      }
    }
  }
  DONE = F
  if((i %% 1000) == 0) cat("% Done:",i/nrow(NOVASEQ_VST_SYM)*100,"\n")
}

#Replace missing gene symbol with original ENSEMBL ID
for(i in 1:nrow(NOVASEQ_VST_SYM)){
  if(NOVASEQ_VST_SYM[i,1] == ""){
    NOVASEQ_VST_SYM[i,1] = rownames(NOVASEQ_VST_SYM)[i]
  }
}


rownames(NOVASEQ_VST_SYM) = make.unique(NOVASEQ_VST_SYM[,1])
NOVASEQ_VST_SYM = NOVASEQ_VST_SYM[,-1]
rownames(NOVASEQ_VST_SYM) = make.unique(gsub('\\..*','',rownames(NOVASEQ_VST_SYM)))

CleanCountData_Controls_nova = NOVASEQ_VST_SYM[rownames(NOVASEQ_VST_SYM) %in% GSEA_genes,] #Correct dimensions

HISEQ_VST_SYM = HISEQ_VST
blank = rep(NA,nrow(HISEQ_VST_SYM))
HISEQ_VST_SYM = cbind(blank,HISEQ_VST_SYM)
HISEQ_VST_SYM[,1] = rownames(HISEQ_VST_SYM)

DONE = F
for(i in 1:nrow(HISEQ_VST_SYM)){
  for(j in 1:nrow(gene_positions)){
    if(DONE == F){
      if(sub("\\..*","",HISEQ_VST_SYM[i,1]) == gene_positions$ensembl_gene_id[j]){
        HISEQ_VST_SYM[i,1] = gene_positions$hgnc_symbol[j]
        DONE = T
      }
    }
  }
  DONE = F
  if((i %% 1000) == 0) cat("% Done:",i/nrow(HISEQ_VST_SYM)*100,"\n")
}

#Replace missing gene symbol with original ENSEMBL ID
for(i in 1:nrow(HISEQ_VST_SYM)){
  if(HISEQ_VST_SYM[i,1] == ""){
    HISEQ_VST_SYM[i,1] = rownames(HISEQ_VST_SYM)[i]
  }
}


rownames(HISEQ_VST_SYM) = make.unique(HISEQ_VST_SYM[,1])
HISEQ_VST_SYM = HISEQ_VST_SYM[,-1]
rownames(HISEQ_VST_SYM) = make.unique(gsub('\\..*','',rownames(HISEQ_VST_SYM)))


CleanCountData_Controls_hiseq = HISEQ_VST_SYM[rownames(HISEQ_VST_SYM) %in% GSEA_genes,] #Correct dimensions

#Can you column bind the hiseq and novaseq cohorts? (in retrospect it would have been faster to do this first)
table(rownames(CleanCountData_Controls_nova) == rownames(CleanCountData_Controls_hiseq))
CleanCountData_Controls = cbind(CleanCountData_Controls_nova,CleanCountData_Controls_hiseq)

#Order rownames alphabetically
#CleanCountData_Controls = CleanCountData_Controls[order(row.names(CleanCountData_Controls)),]

CleanCountData_HC = CleanCountData_Controls[,colnames(CleanCountData_Controls) %in% rownames(HC_coldata)]
CleanCountData_OND = CleanCountData_Controls[,colnames(CleanCountData_Controls) %in% rownames(OND_coldata)]

write.csv(CleanCountData_Controls,"Control_SupervisedClassifier_Matrix_10-25-21.csv")


#Convert periods to dashes 
colnames(CleanCountData_HC) = gsub("\\.","-",colnames(CleanCountData_HC))
colnames(CleanCountData_OND)  = gsub("\\.","-",colnames(CleanCountData_OND))
colnames(CleanCountData_Controls)  = gsub("\\.","-",colnames(CleanCountData_Controls))
colnames(ClassifierMatrix) = gsub("\\.","-",colnames(ClassifierMatrix))


##Combined ALS and Control Classifier Matrices for GSEA
table(rownames(ClassifierMatrix) == rownames(CleanCountData_Controls))
table(rownames(ClassifierMatrix) == rownames(CleanCountData_HC))
table(rownames(ClassifierMatrix) == rownames(CleanCountData_OND))

#All ALS subjects and all control subjects 
AllSubjectsClassifier = cbind(ClassifierMatrix,CleanCountData_Controls)
write.csv(AllSubjectsClassifier,"AllSubjects_ClassifierMatrix.csv")

#All ALS subjects and healthy control subjects
ALS_HC_Classifier = cbind(ClassifierMatrix,CleanCountData_HC)
write.csv(ALS_HC_Classifier,"ALS_HC_ClassifierMatrix.csv")

#All ALS subjects and other neurological disease subjects
ALS_OND_Classifier = cbind(ClassifierMatrix,CleanCountData_OND)
write.csv(ALS_OND_Classifier,"ALS_OND_ClassifierMatrix.csv")

#Save Outputs
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS")
#save.image("GSE153960_ALS_plus_Controls_v3.RData")
load("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/Controls/GSE153960_ALS_plus_Controls_v3.RData")
