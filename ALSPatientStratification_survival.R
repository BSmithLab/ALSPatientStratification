##Survival analysis for ALS subtypes

#Written By: Jarrett Eshima
#For: Dr. Barbara Smith Lab
#Date: June 25th 2021

library(survival)
library(survminer)
library(ranger)
library(ggplot2)
library(dplyr)

##########################################################################################################################################
################################    GSE153960 Survival (n = 451)     #####################################################################
##########################################################################################################################################

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping")
ALSMetaData = read.csv("Eshima_CleanMetaData_ALS_Patients.csv") #File generated in supplemental R script (metadata extractor)
ControlMetaData = read.csv("Eshima_CleanMetaData_Control_Patients.csv") #File generated in supplemental R script (metadata extractor)
EshimaSubjects = c(ALSMetaData$sample_id_alt,ControlMetaData$sample_id_alt)

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Survival")
RosettaStone = read.csv("GEO_Collaborator_PatientData.csv") #Provided by NYGC ALS Consortium, by request
EshimaRosettaStone = RosettaStone[RosettaStone$ExternalSampleId %in% EshimaSubjects,]

#Parse into ALS and Control Rosetta Stones
for(i in 1:nrow(EshimaRosettaStone)){
  if(EshimaRosettaStone$Subject.Group[i] == "Non-Neurological Control"){
    EshimaRosettaStone$Subject.Group[i] = "HC"
  }
}

for(i in 1:nrow(EshimaRosettaStone)){
  if(EshimaRosettaStone$Subject.Group[i] == "Other Neurological Disorders"){
    EshimaRosettaStone$Subject.Group[i] = "OMND"
  }
}

table(EshimaRosettaStone$Subject.Group) #93 HC, 42 OMND, 451 ALS -- GOOD

OI = which(EshimaRosettaStone$Subject.Group == "OMND")
OMNDRosettaStone = EshimaRosettaStone[OI,]
HCI = which(EshimaRosettaStone$Subject.Group == "HC")
HCRosettaStone = EshimaRosettaStone[HCI,]
ALSI = c(OI,HCI)
ALSRosettaStone = EshimaRosettaStone[-ALSI,]
ControlRosettaStone = rbind(HCRosettaStone,OMNDRosettaStone)


#Subtype Labels
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Survival")
mycoldata = read.csv("ALS451_coldata_SUBTYPES.csv") #manually add the "Subject" header to the first column ; file generated in previous script
subtypes = data.frame(mycoldata$Subtype)
rownames(subtypes) = gsub("\\.","-",mycoldata$Subject) 
colnames(subtypes) = "ESubtype"
head(subtypes)

SurvivalPatientData = ALSRosettaStone
tmp = c(1,11,12,15,16,17,19,21,20,22)
SurvivalPatientData = SurvivalPatientData[,-tmp]
colnames(SurvivalPatientData) = c("CGND","ID","Sex","Group","Subcategory","SiteofOnset","SiteExtra","AgeofOnset","AgeofDeath","Duration","Mutation","DiseaseGroup")

table(SurvivalPatientData$ID)
UniqueSubjects = names(table(SurvivalPatientData$ID))

SurvivalMat = data.frame(matrix(NA,nrow=length(UniqueSubjects),ncol=ncol(SurvivalPatientData)))
rownames(SurvivalMat) = UniqueSubjects
colnames(SurvivalMat) = colnames(SurvivalPatientData)

tmp = data.frame(matrix(NA,nrow(SurvivalMat),ncol=3))
colnames(tmp) = c("RNAseq1","RNAseq2","RNAseq3")
rownames(tmp) = rownames(SurvivalMat)

count=1
for(i in 1:nrow(SurvivalMat)){
  for(j in 1:nrow(SurvivalPatientData)){
    if(rownames(SurvivalMat)[i] == SurvivalPatientData$ID[j]){
      SurvivalMat[i,] = SurvivalPatientData[j,]
      tmp[i,count] = SurvivalPatientData$CGND[j]
      count = count+1
    }
  }
  count = 1
}

SurvivalMat2 = cbind(SurvivalMat,tmp)

sublabs = c("ESubtype1","ESubtype2","ESubtype3","Final")
blanks = matrix(NA,nrow(SurvivalMat),ncol=length(sublabs))
colnames(blanks) = sublabs

FinalSurvivalMat = cbind(SurvivalMat2,blanks)


for(i in 1:nrow(FinalSurvivalMat)){
  for(j in 1:nrow(subtypes)){
    if(FinalSurvivalMat$RNAseq1[i] == rownames(subtypes)[j]){
      FinalSurvivalMat$ESubtype1[i] = subtypes$ESubtype[j]
    }
    if(is.na(FinalSurvivalMat$RNAseq2[i]) == F){
      if(FinalSurvivalMat$RNAseq2[i] == rownames(subtypes)[j]){
        FinalSurvivalMat$ESubtype2[i] = subtypes$ESubtype[j]
      }
    }
    
    if(is.na(FinalSurvivalMat$RNAseq3[i]) == F){
      if(FinalSurvivalMat$RNAseq3[i] == rownames(subtypes)[j]){
        FinalSurvivalMat$ESubtype3[i] = subtypes$ESubtype[j]
      }
    }
  }
}


#This is the majority agreement approach
for(i in 1:nrow(FinalSurvivalMat)){
  if(length( which(table(c(FinalSurvivalMat$ESubtype1[i],FinalSurvivalMat$ESubtype2[i],FinalSurvivalMat$ESubtype3[i])) == max(table(c(FinalSurvivalMat$ESubtype1[i],FinalSurvivalMat$ESubtype2[i],FinalSurvivalMat$ESubtype3[i]))))) == 1){
    FinalSurvivalMat$Final[i] = names(which(table(c(FinalSurvivalMat$ESubtype1[i],FinalSurvivalMat$ESubtype2[i],FinalSurvivalMat$ESubtype3[i])) == max(table(c(FinalSurvivalMat$ESubtype1[i],FinalSurvivalMat$ESubtype2[i],FinalSurvivalMat$ESubtype3[i])))))
  }else{
    FinalSurvivalMat$Final[i] = "Discordant"
  }
}


#Clean up missing entries
for(i in 1:nrow(FinalSurvivalMat)){
  if(FinalSurvivalMat$AgeofOnset[i] == "Unknown"){
    FinalSurvivalMat$AgeofOnset[i] = NA
  }
  
  if(FinalSurvivalMat$Duration[i] == "Unknown"){
    FinalSurvivalMat$Duration[i] = NA
  }
}

a = which(FinalSurvivalMat$Duration == "Not Applicable")
FinalSurvivalMat$Duration[a] = NA
b = which(FinalSurvivalMat$AgeofOnset == "Not Applicable")
FinalSurvivalMat$AgeofOnset[b] = NA

#Add status

status = rep(NA,nrow(FinalSurvivalMat))
for(i in 1:nrow(FinalSurvivalMat)){
  if(is.na(FinalSurvivalMat$AgeofOnset[i])){
    status[i] = 0
  }else if (is.na(FinalSurvivalMat$Duration[i])){
    status[i] = 0
  }else if(is.na(FinalSurvivalMat$AgeofDeath[i])){
    status[i] = 0
  }else{
    status[i]= 1
  }
}


FinalSurvival = cbind(FinalSurvivalMat,status)
colnames(FinalSurvival)[10] = "time"
FinalSurvival$time = as.numeric(FinalSurvival$time)

km = with(FinalSurvival,Surv(time,status))
km_fit = survfit(Surv(time, status) ~ Final, data=FinalSurvival)
surv_pvalue(km_fit)
print(km_fit)
summary(km_fit, times = c(1,30,60,90*(1:10)))

par(mar = c(5, 5, 3, 2))
plot(km_fit,col = c("darkorchid2","goldenrod","navy","firebrick"),main="ALS Subtype Kaplan-Meier Survival",xlab="Disease Duration (months)",ylab="Survival Probability",lwd=2.5,cex.axis = 1.5,cex.lab=1.5,cex.main=1.5,xaxt="n")
timeseq = seq(0,156,12)
axis(side = 1, at = timeseq,labels = T,tick = T,cex.axis = 1.5)
legend(132,1.02,legend=c("ALS-Glia","ALS-Ox","ALS-TD","Discordant"),lty = 1,lwd=3.5,col=c("goldenrod","navy","firebrick","darkorchid2"))
text(10,0.2,"p = 0.023",cex=1.5)

segments(0,0.5,42,0.5,lty="dashed",lwd=2.5)
segments(28,0,28,0.5,lty="dashed",lwd=2.5)
segments(36,0,36,0.5,lty="dashed",lwd=2.5)
segments(42,0,42,0.5,lty="dashed",lwd=2.5)


#############################################################################################################################

## Considering only subjects with agreement between motor and frontal cortex samples ("Discordant" removed)

#Remove Discordant
tmpind = which(FinalSurvival$Final == 'Discordant') #Exclude discordant Subjects
CleanSurvivalMat = FinalSurvival[-tmpind,]
tail(CleanSurvivalMat)

#Add status

status = rep(NA,nrow(CleanSurvivalMat))
for(i in 1:nrow(CleanSurvivalMat)){
  if(is.na(CleanSurvivalMat$AgeofOnset[i])){
    status[i] = 0
  }else if (is.na(CleanSurvivalMat$time[i])){
    status[i] = 0
  }else if(is.na(CleanSurvivalMat$AgeofDeath[i])){
    status[i] = 0
  }else{
    status[i]= 1
  }
}

CleanSurvivalMat = cbind(CleanSurvivalMat,status)
CleanSurvivalMat$time = as.numeric(CleanSurvivalMat$time)

par(mar = c(5, 5, 4, 3))
km = with(CleanSurvivalMat,Surv(time,status))
km_fit = survfit(Surv(time, status) ~ Final, data=CleanSurvivalMat)
surv_pvalue(km_fit)
print(km_fit)
summary(km_fit, times = c(1,30,60,90*(1:10)))
plot(km_fit,col = c("goldenrod","navy","firebrick"),main="ALS Subtype Kaplan-Meier Survival",xlab="Disease Duration (months)",ylab="Survival Probability",lwd=2.5,cex.axis = 1.5,cex.lab=1.5,cex.main=1.5,xaxt="n")
timeseq = seq(0,156,12)
axis(side = 1, at = timeseq,labels = T,tick = T,cex.axis = 1.5)
legend(135,1.02,legend=c("ALS-Glia","ALS-Ox","ALS-TD"),lty = 1,lwd=3.5,col=c("goldenrod","navy","firebrick"))
text(10,0.2,"p = 0.0095",cex=1.5)

segments(0,0.5,42,0.5,lty="dashed",lwd=2.5)
segments(28,0,28,0.5,lty="dashed",lwd=2.5)
segments(36,0,36,0.5,lty="dashed",lwd=2.5)
segments(42,0,42,0.5,lty="dashed",lwd=2.5)


#############################################################################################################################

#Pairwise survival comparisons - No Discordant

tmp = which(CleanSurvivalMat$Final == "GLIA")
CleanSurvivalMat_NoGlia = CleanSurvivalMat[-tmp,]

tmp = which(CleanSurvivalMat$Final == "OX")
CleanSurvivalMat_NoOx = CleanSurvivalMat[-tmp,]

tmp = which(CleanSurvivalMat$Final == "TE")
CleanSurvivalMat_NoTD = CleanSurvivalMat[-tmp,]


par(mar = c(5, 5, 4, 3))
km = with(CleanSurvivalMat_NoGlia,Surv(time,status))
km_fit = survfit(Surv(time, status) ~ Final, data=CleanSurvivalMat_NoGlia)
surv_pvalue(km_fit)
print(km_fit)


par(mar = c(5, 5, 4, 3))
km = with(CleanSurvivalMat_NoOx,Surv(time,status))
km_fit = survfit(Surv(time, status) ~ Final, data=CleanSurvivalMat_NoOx)
surv_pvalue(km_fit)
print(km_fit)


par(mar = c(5, 5, 4, 3))
km = with(CleanSurvivalMat_NoTD,Surv(time,status))
km_fit = survfit(Surv(time, status) ~ Final, data=CleanSurvivalMat_NoTD)
surv_pvalue(km_fit)
print(km_fit)


#Pairwise survival comparisons - With Discordant

tmp = which(FinalSurvival$Final == "GLIA")
CleanSurvivalMat_MinGD = FinalSurvival[-tmp,]
tmp = which(CleanSurvivalMat_MinGD$Final == "Discordant")
CleanSurvivalMat_MinGD = CleanSurvivalMat_MinGD[-tmp,]

tmp = which(FinalSurvival$Final == "OX")
CleanSurvivalMat_MinGO = FinalSurvival[-tmp,]
tmp = which(CleanSurvivalMat_MinGO$Final == "GLIA")
CleanSurvivalMat_MinGO = CleanSurvivalMat_MinGO[-tmp,]

tmp = which(FinalSurvival$Final == "TE")
CleanSurvivalMat_MinTD = FinalSurvival[-tmp,]
tmp = which(CleanSurvivalMat_MinTD$Final == "Discordant")
CleanSurvivalMat_MinTD = CleanSurvivalMat_MinTD[-tmp,]

tmp = which(FinalSurvival$Final == "TE")
CleanSurvivalMat_MinTO = FinalSurvival[-tmp,]
tmp = which(CleanSurvivalMat_MinTO$Final == "OX")
CleanSurvivalMat_MinTO = CleanSurvivalMat_MinTO[-tmp,]

tmp = which(FinalSurvival$Final == "TE")
CleanSurvivalMat_MinTG = FinalSurvival[-tmp,]
tmp = which(CleanSurvivalMat_MinTG$Final == "GLIA")
CleanSurvivalMat_MinTG = CleanSurvivalMat_MinTG[-tmp,]

tmp = which(FinalSurvival$Final == "Discordant")
CleanSurvivalMat_MinDO = FinalSurvival[-tmp,]
tmp = which(CleanSurvivalMat_MinDO$Final == "OX")
CleanSurvivalMat_MinDO = CleanSurvivalMat_MinDO[-tmp,]

#Ox vs TD - p=0.3
par(mar = c(5, 5, 4, 3))
km = with(CleanSurvivalMat_MinGD,Surv(time,status))
km_fit = survfit(Surv(time, status) ~ Final, data=CleanSurvivalMat_MinGD)
surv_pvalue(km_fit)
print(km_fit)

#TD vs Discordant - p=0.25
par(mar = c(5, 5, 4, 3))
km = with(CleanSurvivalMat_MinGO,Surv(time,status))
km_fit = survfit(Surv(time, status) ~ Final, data=CleanSurvivalMat_MinGO)
surv_pvalue(km_fit)
print(km_fit)

#Ox vs Glia - p=0.015
par(mar = c(5, 5, 4, 3))
km = with(CleanSurvivalMat_MinTD,Surv(time,status))
km_fit = survfit(Surv(time, status) ~ Final, data=CleanSurvivalMat_MinTD)
surv_pvalue(km_fit)
print(km_fit)

#Glia vs Discordant - p=0.091
par(mar = c(5, 5, 4, 3))
km = with(CleanSurvivalMat_MinTO,Surv(time,status))
km_fit = survfit(Surv(time, status) ~ Final, data=CleanSurvivalMat_MinTO)
surv_pvalue(km_fit)
print(km_fit)

#Ox vs Discordant - p=0.76
par(mar = c(5, 5, 4, 3))
km = with(CleanSurvivalMat_MinTG,Surv(time,status))
km_fit = survfit(Surv(time, status) ~ Final, data=CleanSurvivalMat_MinTG)
surv_pvalue(km_fit)
print(km_fit)

#Glia vs TD - p=0.0043
par(mar = c(5, 5, 4, 3))
km = with(CleanSurvivalMat_MinDO,Surv(time,status))
km_fit = survfit(Surv(time, status) ~ Final, data=CleanSurvivalMat_MinDO)
surv_pvalue(km_fit)
print(km_fit)
