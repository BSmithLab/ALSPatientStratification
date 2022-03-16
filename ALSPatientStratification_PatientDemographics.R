##Patient Demographics

#Written By: Jarrett Eshima
#For: Dr. Barbara Smith Lab
#Date: September 2021


setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping")
ALSMetaData = read.csv("Eshima_CleanMetaData_ALS_Patients.csv") #File generated in previous script
HCMetaData = read.csv("Eshima_CleanMetaData_HC_10-15-21.csv") #File generated in previous script
ONDMetaData = read.csv("Eshima_CleanMetaData_OND_10-15-21.csv") #File generated in previous script
EshimaSubjects = c(ALSMetaData$sample_id_alt,HCMetaData$sample_id_alt,ONDMetaData$sample_id_alt)

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Survival")
Demographics = read.csv("GEO_Collaborator_PatientData.csv") #File provided by NYGC ALS Consortium, by request
EshimaDemographics = Demographics[Demographics$ExternalSampleId %in% EshimaSubjects,]

#Parse into ALS and Control Rosetta Stones
for(i in 1:nrow(EshimaDemographics)){
  if(EshimaDemographics$Subject.Group[i] == "Non-Neurological Control"){
    EshimaDemographics$Subject.Group[i] = "HC"
  }
}

for(i in 1:nrow(EshimaDemographics)){
  if(EshimaDemographics$Subject.Group[i] == "Other Neurological Disorders"){
    EshimaDemographics$Subject.Group[i] = "OMND"
  }
}


OI = which(EshimaDemographics$Subject.Group == "OMND")
OMNDDemographics = EshimaDemographics[OI,]
HCI = which(EshimaDemographics$Subject.Group == "HC")
HCDemographics = EshimaDemographics[HCI,]
ALSI = c(OI,HCI)
ALSDemographics = EshimaDemographics[-ALSI,]
ControlDemographics = rbind(HCDemographics,OMNDDemographics)

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Survival")
mycoldata = read.csv("ALS451_coldata_SUBTYPES.csv") #File generated in previous script
subtypes = data.frame(mycoldata$Subtype)
rownames(subtypes) = mycoldata$Subject
colnames(subtypes) = "ESubtype"
rownames(subtypes) = gsub("\\.","-",rownames(subtypes))
head(subtypes)

#############################################################################################################
#ALS Patient Demographics

PatientData = ALSDemographics

table(PatientData$ExternalSubjectId)
UniqueSubjects = names(table(PatientData$ExternalSubjectId))

tmpALS = data.frame(matrix(NA,nrow=length(UniqueSubjects),ncol=ncol(PatientData)))
rownames(tmpALS) = UniqueSubjects
colnames(tmpALS) = colnames(PatientData)

tmp = data.frame(matrix(NA,nrow(tmpALS),ncol=3))
colnames(tmp) = c("RNAseq1","RNAseq2","RNAseq3")
rownames(tmp) = rownames(tmpALS)

count=1
for(i in 1:nrow(tmpALS)){
  for(j in 1:nrow(PatientData)){
    if(rownames(tmpALS)[i] == PatientData$ExternalSubjectId[j]){
      tmpALS[i,] = PatientData[j,]
      tmp[i,count] = PatientData$ExternalSampleId[j]
      count = count+1
    }
  }
  count = 1
}

FinalALSDemographics = cbind(tmpALS,tmp)
ALSDemographics_SampleConverted = FinalALSDemographics

for(i in 1:nrow(ALSDemographics_SampleConverted)){
  for(j in 1:nrow(ALSDemographics)){
    
    if(! is.na(ALSDemographics_SampleConverted$RNAseq1[i])){
      if(ALSDemographics_SampleConverted$RNAseq1[i] == ALSDemographics$ExternalSampleId[j]){
        ALSDemographics_SampleConverted$RNAseq1[i] = ALSDemographics$Sample.Source[j]
      }
    }
    
    if(! is.na(ALSDemographics_SampleConverted$RNAseq2[i])){
      if(ALSDemographics_SampleConverted$RNAseq2[i] == ALSDemographics$ExternalSampleId[j]){
        ALSDemographics_SampleConverted$RNAseq2[i] = ALSDemographics$Sample.Source[j]
      }
    }
    
    if(! is.na(ALSDemographics_SampleConverted$RNAseq3[i])){
      if(ALSDemographics_SampleConverted$RNAseq3[i] == ALSDemographics$ExternalSampleId[j]){
        ALSDemographics_SampleConverted$RNAseq3[i] = ALSDemographics$Sample.Source[j]
      }
    }
    
  }
  if((i %% 10) == 0) cat("% Done:",i/nrow(ALSDemographics_SampleConverted)*100,"\n")
}

table(ALSDemographics_SampleConverted$Sex)

sublabs = c("ESubtype1","ESubtype2","ESubtype3","Final")
blanks = matrix(NA,nrow(FinalALSDemographics),ncol=length(sublabs))
colnames(blanks) = sublabs

FinalSurvivalMat = cbind(FinalALSDemographics,blanks)


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


table(FinalSurvivalMat$Final)

G = which(FinalSurvivalMat$Final == "GLIA")
O = which(FinalSurvivalMat$Final == "OX")
TD = which(FinalSurvivalMat$Final == "TE")
D = which(FinalSurvivalMat$Final == "Discordant")
G = FinalSurvivalMat[G,]
O = FinalSurvivalMat[O,]
TD = FinalSurvivalMat[TD,]
D = FinalSurvivalMat[D,]

#Age of onset - summary metrics
tmp = O$Age.at.Symptom.Onset
tmp = as.numeric(tmp)
tmp = tmp[!is.na(tmp)]
mean(tmp)
sd(tmp)/sqrt(length(tmp))

#Age of death - summary metrics
tmp = G$Age.at.Death
tmp = as.numeric(tmp)
tmp = tmp[!is.na(tmp)]
mean(tmp)
sd(tmp)/sqrt(length(tmp))

#Site of onset - summary metrics
tmp = O$Site.of.Motor.Onset
table(tmp)

#FTLD Comorbidity - summary metrics
#Site of onset
tmp = D$disease_group
table(tmp)

#Disease Duration - summary metrics
tmp = G$Disease.Duration.in.Months
tmp = as.numeric(tmp)
tmp = tmp[!is.na(tmp)]
mean(tmp)
sd(tmp)/sqrt(length(tmp))


### Across all ALS Patients
#Age of onset
tmp = FinalALSDemographics$Age.at.Symptom.Onset
tmp = as.numeric(tmp)
tmp = tmp[!is.na(tmp)]
mean(tmp)
sd(tmp)/sqrt(length(tmp))

#Age of death
tmp = FinalALSDemographics$Age.at.Death
tmp = as.numeric(tmp)
tmp = tmp[!is.na(tmp)]
mean(tmp)
sd(tmp)/sqrt(length(tmp))

#Disease Duration
tmp = FinalALSDemographics$Disease.Duration.in.Months
tmp = as.numeric(tmp)
tmp = tmp[!is.na(tmp)]
mean(tmp)
sd(tmp)/sqrt(length(tmp))



############################################################################################################################################
################################ Plot clinical parameters ##################################################################################
############################################################################################################################################

#This section considers clinical parameters without discordant subjects
#Discordant subjects are considered in the next section of this script

################################################# Age of onset

TEaoo = as.numeric(TE$Age.at.Symptom.Onset)
TEaoo = TEaoo[!is.na(TEaoo)]
OXaoo = as.numeric(O$Age.at.Symptom.Onset)
OXaoo = OXaoo[!is.na(OXaoo)]
GLaoo = as.numeric(G$Age.at.Symptom.Onset)
GLaoo = GLaoo[!is.na(GLaoo)]

par(mar = c(5, 5, 3, 2))
boxplot(GLaoo,OXaoo,TEaoo,xaxt="n",main=c("Subtype Age of Onset"),cex.axis = 1.5,cex.main=1.5,col=c("goldenrod1","navy","firebrick"),xlab="Subtype",ylab="Age of Onset (years)",cex.lab=1.5,pch=20,cex=1.4)
axis(at=1:3,side=1,labels=c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=1.5)

#Quick ANOVA to check for post-hoc tests
Age = FinalSurvivalMat$Age.at.Symptom.Onset
ST = FinalSurvivalMat$Final
anovadat = data.frame(cbind(Age,ST))
anovadat$ST = factor(anovadat$ST,ordered = T)
levels(anovadat$ST)

tmpindex = rep(NA,nrow(anovadat))
count=1
for(i in 1:nrow(anovadat)){
  if(anovadat$Age[i] == "Unknown"){
    tmpindex[count] = i
    count = count+1
  }else if(anovadat$Age[i] == "Not Applicable"){
    tmpindex[count] = i
    count = count+1
  }else if(anovadat$ST[i] == "Discordant"){
    tmpindex[count] = i
    count = count+1
  }
}
tmpindex = tmpindex[!is.na(tmpindex)]
NoDisc = anovadat[-tmpindex,]

oneway = aov(Age~ST,data=NoDisc)
summary(oneway) #suggests no significant differences in age of onset by subtype

#Post hoc (not necessary)
t.test(OXaoo,GLaoo) #p = 0.2
t.test(TEaoo,GLaoo) #p = 0.85
t.test(TEaoo,OXaoo) #p = 0.25



################################################# Site of onset (categories: Bulbar, Limb, Other)

#There is extra code at the end of the script to consider all categories: Axial, Axial-Limb, Axial-Bulbar, Bulbar, Limb, Bulbar-Limb, Generalized, Unknown
SiteofOnset = FinalSurvivalMat
for(i in 1:nrow(SiteofOnset)){
  if(SiteofOnset$Site.of.Motor.Onset[i] != "Limb" && SiteofOnset$Site.of.Motor.Onset[i] != "Bulbar"){
    SiteofOnset$Site.of.Motor.Onset[i] = "Other"
  }
}


tmpdata = SiteofOnset

ite = which(tmpdata$Final == "TE")
iox = which(tmpdata$Final == "OX")
igl = which(tmpdata$Final == "GLIA")

tmpTE = tmpdata[ite,]
tmpOX = tmpdata[iox,]
tmpGL = tmpdata[igl,]

TEsoo = tmpTE$Site.of.Motor.Onset
OXsoo = tmpOX$Site.of.Motor.Onset
GLsoo = tmpGL$Site.of.Motor.Onset

a = table(TEsoo)
b = table(OXsoo)
c = table(GLsoo)

barplot(a,main="ALS-TE Site of Onset",names.arg = c("Bulbar","Limb","Other"),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="firebrick",ylim=c(0,60))
barplot(b,main="ALS-Ox Site of Onset",names.arg = c("Bulbar","Limb","Other"),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="navy",ylim=c(0,60))
barplot(c,main="ALS-Glia Site of Onset",names.arg = c("Bulbar","Limb","Other"),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="goldenrod",ylim=c(0,60))

#Chi-squared test of independence (categorical/freq data)
CT = matrix(c(a,b,c),nrow = 3,ncol = 3)
rownames(CT)= c("Bulbar","Limb","Other")
colnames(CT)= c("ALS-TE","ALS-OX","ALS-GLIA")
CT = t(CT)

chisq.test(CT) # p = 0.33 
#Do not reject the null that Site of Onset is independent of ALS Subtype - interesting

#Chi-squared goodness of fit test (Bulbar:Limb 1:2?)
chisq.test(c(11,20),p=c(1/3,2/3)) #Glia 
chisq.test(c(23,51),p=c(1/3,2/3)) #Ox
chisq.test(c(17,35),p=c(1/3,2/3)) #TE
#results suggest that a 1:2 ratio of bulbar:limb onset is a good fit for all subtypes

#Extra code to look at just Limb and Bulbar sites of symptom onset

# ind1 = which(SiteofOnset$Site.of.Motor.Onset == "Limb")
# ind2 = which(SiteofOnset$Site.of.Motor.Onset == "Bulbar")
# ind = c(ind1,ind2)
# tmpdata = SiteofOnset[ind,]
# ite = which(tmpdata$Final == "TE")
# iox = which(tmpdata$Final == "OX")
# igl = which(tmpdata$Final == "GLIA")
# 
# tmpTE = tmpdata[ite,]
# tmpOX = tmpdata[iox,]
# tmpGL = tmpdata[igl,]
# 
# TEsoo = tmpTE$Site.of.Motor.Onset
# OXsoo = tmpOX$Site.of.Motor.Onset
# GLsoo = tmpGL$Site.of.Motor.Onset
# 
# a = table(TEsoo)
# b = table(OXsoo)
# c = table(GLsoo)
# barplot(a,main="ALS-TE Site of Onset",names.arg = c("Bulbar","Limb"),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="firebrick",ylim=c(0,60))
# barplot(b,main="ALS-Ox Site of Onset",names.arg = c("Bulbar","Limb"),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="navy",ylim=c(0,60))
# barplot(c,main="ALS-Glia Site of Onset",names.arg = c("Bulbar","Limb"),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="goldenrod",ylim=c(0,60))



################################################# Age of death 

TEaod = as.numeric(TE$Age.at.Death)
TEaod = TEaod[!is.na(TEaod)]
OXaod = as.numeric(O$Age.at.Death)
OXaod = OXaod[!is.na(OXaod)]
GLaod = as.numeric(G$Age.at.Death)
GLaod = GLaod[!is.na(GLaod)]

par(mar = c(5, 5, 3, 2))
boxplot(GLaod,OXaod,TEaod,xaxt="n",main=c("Subtype Age of Death"),cex.axis = 1.5,cex.main=1.5,col=c("goldenrod1","navy","firebrick"),xlab="Subtype",ylab="Age of Death (years)",cex.lab=1.5,pch=20,cex=1.4)
axis(at=1:3,side=1,labels=c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=1.5)

#Quick ANOVA to check for post-hoc tests
Age = FinalSurvivalMat$Age.at.Death
ST = FinalSurvivalMat$Final
anovadat = data.frame(cbind(Age,ST))
anovadat$ST = factor(anovadat$ST,ordered = T)
levels(anovadat$ST)

tmpindex = rep(NA,nrow(anovadat))
count=1
for(i in 1:nrow(anovadat)){
  if(anovadat$Age[i] == "Unknown"){
    tmpindex[count] = i
    count = count+1
  }else if(anovadat$Age[i] == "Not Applicable"){
    tmpindex[count] = i
    count = count+1
  }else if(anovadat$ST[i] == "Discordant"){
    tmpindex[count] = i
    count = count+1
  }
}
tmpindex = tmpindex[!is.na(tmpindex)]
NoDisc2 = anovadat[-tmpindex,]

oneway = aov(Age~ST,data=NoDisc2)
summary(oneway) #suggests no significant differences in age of onset by subtype

#Post hoc (not necessary)
t.test(OXaod,GLaod) #p = 0.27
t.test(TEaod,GLaod) #p = 0.77
t.test(TEaod,OXaod) #p = 0.11


################################################# FTD Comorbidity
G.FTD = names(table(G$Subject.Group.Subcategory)) #3 and 5 are FTD
blank = rep(NA,nrow(G))
for(j in 1:nrow(G)){
  if(G$Subject.Group.Subcategory[j] == G.FTD[3]){
    blank[j] = "Positive"
  }else if(G$Subject.Group.Subcategory[j] == G.FTD[5]){
    blank[j] = "Positive"
  }else{
    blank[j] = "Negative"
  }
}
G2 = cbind(G,blank)
nc = ncol(G2)
colnames(G2)[nc] = "FTD.Comorbidity"


O.FTD = names(table(O$Subject.Group.Subcategory))
blank = rep(NA,nrow(O))
for(j in 1:nrow(O)){
  if(O$Subject.Group.Subcategory[j] == O.FTD[6]){
    blank[j] = "Positive"
  }else if(O$Subject.Group.Subcategory[j] == O.FTD[7]){
    blank[j] = "Positive"
  }else{
    blank[j] = "Negative"
  }
}
O2 = cbind(O,blank)
nc = ncol(O2)
colnames(O2)[nc] = "FTD.Comorbidity"


TE.FTD = names(table(TE$Subject.Group.Subcategory))
blank = rep(NA,nrow(TE))
for(j in 1:nrow(TE)){
  if(TE$Subject.Group.Subcategory[j] == TE.FTD[6]){
    blank[j] = "Positive"
  }else{
    blank[j] = "Negative"
  }
}
TE2 = cbind(TE,blank)
nc = ncol(TE2)
colnames(TE2)[nc] = "FTD.Comorbidity"


x = table(G2$FTD.Comorbidity)
y = table(O2$FTD.Comorbidity)
z = table(TE2$FTD.Comorbidity)

#Chi-squared test for independence
CT = matrix(c(x,y,z),nrow = 2,ncol = 3)
rownames(CT)= c("Negative","Positive")
colnames(CT)= c("ALS-Glia","ALS-Ox","ALS-TD")
CT = t(CT)

chisq.test(CT) # p = 0.594
#Do not reject the null that FTD comorbidity is independent of ALS-subtype

visdata.g = x[[2]]/(x[[1]]+x[[2]])
visdata.o = y[[2]]/(y[[1]]+y[[2]])
visdata.t = z[[2]]/(z[[1]]+z[[2]])
vd = c(visdata.g,visdata.o,visdata.t)

barplot(vd*100,main="ALS Subtype FTLD Comorbidity",names.arg = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5,xlab="Subtype",ylab="FTLD Comorbidity (%)",col=c("goldenrod1","navy","firebrick"),ylim = c(0,25))

###############################################################################################################################
################################## Including Discordant #######################################################################
###############################################################################################################################

################################################# Age of onset

TEaoo = as.numeric(TE$Age.at.Symptom.Onset)
TEaoo = TEaoo[!is.na(TEaoo)]
OXaoo = as.numeric(O$Age.at.Symptom.Onset)
OXaoo = OXaoo[!is.na(OXaoo)]
GLaoo = as.numeric(G$Age.at.Symptom.Onset)
GLaoo = GLaoo[!is.na(GLaoo)]
DIaoo = as.numeric(D$Age.at.Symptom.Onset)
DIaoo = DIaoo[!is.na(DIaoo)]

par(mar = c(5,5,4,3))
boxplot(GLaoo,OXaoo,TEaoo,DIaoo,xaxt="n",main=c("Subtype Age of Onset"),cex.axis = 1.5,cex.main=1.5,col=c("goldenrod1","navy","firebrick","darkorchid2"),xlab="Subtype",ylab="Age of Onset (years)",cex.lab=1.5)
axis(at=1:4,side=1,labels=c("ALS-Glia","ALS-Ox","ALS-TD","Discordant"),cex.axis=1.25)

#Quick ANOVA to check for post-hoc tests
Age = FinalSurvivalMat$Age.at.Symptom.Onset
ST = FinalSurvivalMat$Final
anovadat = data.frame(cbind(Age,ST))
anovadat$ST = factor(anovadat$ST,ordered = T)
levels(anovadat$ST)

tmpindex = rep(NA,nrow(anovadat))
count=1
for(i in 1:nrow(anovadat)){
  if(anovadat$Age[i] == "Unknown"){
    tmpindex[count] = i
    count = count+1
  }else if(anovadat$Age[i] == "Not Applicable"){
    tmpindex[count] = i
    count = count+1
  }
}
tmpindex = tmpindex[!is.na(tmpindex)]
WithDisc = anovadat[-tmpindex,]

oneway = aov(Age~ST,data=WithDisc)
summary(oneway) #suggests no significant differences in age of onset by subtype

#Post hoc (not necessary)
t.test(OXaoo,GLaoo) #p = 0.2
t.test(TEaoo,GLaoo) #p = 0.85
t.test(TEaoo,OXaoo) #p = 0.25
t.test(OXaoo,DIaoo) #p = 0.82
t.test(TEaoo,DIaoo) #p = 0.46
t.test(GLaoo,DIaoo) #p = 0.38




################################################# Site of onset (categories: Bulbar, Limb, Other)

SiteofOnset = FinalSurvivalMat
for(i in 1:nrow(SiteofOnset)){
  if(SiteofOnset$Site.of.Motor.Onset[i] != "Limb" && SiteofOnset$Site.of.Motor.Onset[i] != "Bulbar"){
    SiteofOnset$Site.of.Motor.Onset[i] = "Other"
  }
}


tmpdata = SiteofOnset

ite = which(tmpdata$Final == "TE")
iox = which(tmpdata$Final == "OX")
igl = which(tmpdata$Final == "GLIA")
idi = which(tmpdata$Final == "Discordant")

tmpTE = tmpdata[ite,]
tmpOX = tmpdata[iox,]
tmpGL = tmpdata[igl,]
tmpDI = tmpdata[idi,]

TEsoo = tmpTE$Site.of.Motor.Onset
OXsoo = tmpOX$Site.of.Motor.Onset
GLsoo = tmpGL$Site.of.Motor.Onset
DIsoo = tmpDI$Site.of.Motor.Onset

a = table(TEsoo)
b = table(OXsoo)
c = table(GLsoo)
d = table(DIsoo)

barplot(a,main="ALS-TE Site of Onset",names.arg = c("Bulbar","Limb","Other"),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="firebrick",ylim=c(0,60))
barplot(b,main="ALS-Ox Site of Onset",names.arg = c("Bulbar","Limb","Other"),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="navy",ylim=c(0,60))
barplot(c,main="ALS-Glia Site of Onset",names.arg = c("Bulbar","Limb","Other"),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="goldenrod",ylim=c(0,60))
barplot(d,main="ALS-Discordant Site of Onset",names.arg = c("Bulbar","Limb","Other"),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="darkorchid2",ylim=c(0,60))

#Chi-squared test (categorical/freq data)
CT = matrix(c(a,b,c,d),nrow = 3,ncol = 4)
rownames(CT)= c("Bulbar","Limb","Other")
colnames(CT)= c("ALS-TE","ALS-OX","ALS-GLIA","ALS-Discordant")
CT = t(CT)

chisq.test(CT) # p = 0.563 

#Chi-squared goodness of fit test (Bulbar:Limb 1:2?)
chisq.test(c(11,20),p=c(1/3,2/3)) #Glia 
chisq.test(c(23,51),p=c(1/3,2/3)) #Ox
chisq.test(c(17,35),p=c(1/3,2/3)) #TE
chisq.test(c(8,19),p=c(1/3,2/3)) #Discordant
#results suggest that a 1:2 ratio of bulbar:limb onset is a good fit for all subtypes and discordant subjects



################################################# Age of death 

TEaod = as.numeric(TE$Age.at.Death)
TEaod = TEaod[!is.na(TEaod)]
OXaod = as.numeric(O$Age.at.Death)
OXaod = OXaod[!is.na(OXaod)]
GLaod = as.numeric(G$Age.at.Death)
GLaod = GLaod[!is.na(GLaod)]
DIaod = as.numeric(D$Age.at.Death)
DIaod = DIaod[!is.na(DIaod)]

par(mar = c(5,5,4,3))
boxplot(GLaod,OXaod,TEaod,DIaod,xaxt="n",main=c("Subtype Age of Death"),cex.axis = 1.5,cex.main=1.5,col=c("goldenrod1","navy","firebrick","darkorchid2"),xlab="Subtype",ylab="Age of Death",cex.lab=1.5)
axis(at=1:4,side=1,labels=c("ALS-Glia","ALS-Ox","ALS-TD","Discordant"),cex.axis=1.25)

#Quick ANOVA to check for post-hoc tests
Age = FinalSurvivalMat$Age.at.Death
ST = FinalSurvivalMat$Final
anovadat = data.frame(cbind(Age,ST))
anovadat$ST = factor(anovadat$ST,ordered = T)
levels(anovadat$ST)

oneway = aov(Age~ST,data=anovadat)
summary(oneway) #suggests no significant differences in age of onset by subtype

#Post hoc (not necessary)
t.test(OXaod,GLaod) #p = 0.27
t.test(TEaod,GLaod) #p = 0.77
t.test(TEaod,OXaod) #p = 0.11
t.test(OXaod,DIaod) #p = 0.98
t.test(TEaod,DIaod) #p = 0.21
t.test(GLaod,DIaod) #p = 0.37



################################################# FTD Comorbidity
D.FTD = names(table(D$Subject.Group.Subcategory))
blank = rep(NA,nrow(D))
for(j in 1:nrow(D)){
  if(D$Subject.Group.Subcategory[j] == D.FTD[4]){
    blank[j] = "Positive"
  }else{
    blank[j] = "Negative"
  }
}
D2 = cbind(D,blank)
nc = ncol(D2)
colnames(D2)[nc] = "FTD.Comorbidity"

x = table(G2$FTD.Comorbidity)
y = table(O2$FTD.Comorbidity)
z = table(TE2$FTD.Comorbidity)

z2 = table(D2$FTD.Comorbidity)

#Chi-squared test for independence
CT = matrix(c(x,y,z,z2),nrow = 2,ncol = 4)
rownames(CT)= c("Negative","Positive")
colnames(CT)= c("ALS-TE","ALS-OX","ALS-GLIA","ALS-Discordant")
CT = t(CT)

chisq.test(CT) # p = 0.791

visdata.g = x[[2]]/(x[[1]]+x[[2]])
visdata.o = y[[2]]/(y[[1]]+y[[2]])
visdata.t = z[[2]]/(z[[1]]+z[[2]])
visdata.d = z2[[2]]/(z2[[1]]+z2[[2]])
vd = c(visdata.g,visdata.o,visdata.t,visdata.d)

barplot(vd,main="ALS Subtype FTLD Comorbidity",names.arg = c("ALS-Glia","ALS-Ox","ALS-TD","Discordant"),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="FTD Comorbidity (%)",col=c("goldenrod1","navy","firebrick","darkorchid2"),ylim = c(0,0.25))

###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
#HC Patient Demographics

PatientData = HCDemographics

table(PatientData$ExternalSubjectId)
UniqueSubjects = names(table(PatientData$ExternalSubjectId))

tmpHC = data.frame(matrix(NA,nrow=length(UniqueSubjects),ncol=ncol(PatientData)))
rownames(tmpHC) = UniqueSubjects
colnames(tmpHC) = colnames(PatientData)

tmp = data.frame(matrix(NA,nrow(tmpHC),ncol=3))
colnames(tmp) = c("RNAseq1","RNAseq2","RNAseq3")
rownames(tmp) = rownames(tmpHC)

count=1
for(i in 1:nrow(tmpHC)){
  for(j in 1:nrow(PatientData)){
    if(rownames(tmpHC)[i] == PatientData$ExternalSubjectId[j]){
      tmpHC[i,] = PatientData[j,]
      tmp[i,count] = PatientData$ExternalSampleId[j]
      count = count+1
    }
  }
  count = 1
}

FinalHCDemographics = cbind(tmpHC,tmp)
HCDemographics_SampleConverted = FinalHCDemographics

for(i in 1:nrow(HCDemographics_SampleConverted)){
  for(j in 1:nrow(HCDemographics)){
    
    if(! is.na(HCDemographics_SampleConverted$RNAseq1[i])){
      if(HCDemographics_SampleConverted$RNAseq1[i] == HCDemographics$ExternalSampleId[j]){
        HCDemographics_SampleConverted$RNAseq1[i] = HCDemographics$Sample.Source[j]
      }
    }
    
    if(! is.na(HCDemographics_SampleConverted$RNAseq2[i])){
      if(HCDemographics_SampleConverted$RNAseq2[i] == HCDemographics$ExternalSampleId[j]){
        HCDemographics_SampleConverted$RNAseq2[i] = HCDemographics$Sample.Source[j]
      }
    }
    
    if(! is.na(HCDemographics_SampleConverted$RNAseq3[i])){
      if(HCDemographics_SampleConverted$RNAseq3[i] == HCDemographics$ExternalSampleId[j]){
        HCDemographics_SampleConverted$RNAseq3[i] = HCDemographics$Sample.Source[j]
      }
    }
    
  }
  if((i %% 10) == 0) cat("% Done:",i/nrow(HCDemographics_SampleConverted)*100,"\n")
}

#Sex
table(HCDemographics_SampleConverted$Sex)
#Tissue Site
tmpcontainer = c(HCDemographics_SampleConverted$RNAseq1,HCDemographics_SampleConverted$RNAseq2,HCDemographics_SampleConverted$RNAseq3)
table(tmpcontainer)
#Age of Death
index = which(HCDemographics_SampleConverted$Age.at.Death =="90 or Older")
HCDemographics_SampleConverted$Age.at.Death[index] = 90
mean(as.numeric(HCDemographics_SampleConverted$Age.at.Death))
sd(as.numeric(HCDemographics_SampleConverted$Age.at.Death))
#FTLD not applicable for HC group

###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
#OND Patient Demographics

PatientData = OMNDDemographics

table(PatientData$ExternalSubjectId)
UniqueSubjects = names(table(PatientData$ExternalSubjectId))

tmpOMND = data.frame(matrix(NA,nrow=length(UniqueSubjects),ncol=ncol(PatientData)))
rownames(tmpOMND) = UniqueSubjects
colnames(tmpOMND) = colnames(PatientData)

tmp = data.frame(matrix(NA,nrow(tmpOMND),ncol=3))
colnames(tmp) = c("RNAseq1","RNAseq2","RNAseq3")
rownames(tmp) = rownames(tmpOMND)

count=1
for(i in 1:nrow(tmpOMND)){
  for(j in 1:nrow(PatientData)){
    if(rownames(tmpOMND)[i] == PatientData$ExternalSubjectId[j]){
      tmpOMND[i,] = PatientData[j,]
      tmp[i,count] = PatientData$ExternalSampleId[j]
      count = count+1
    }
  }
  count = 1
}

FinalOMNDDemographics = cbind(tmpOMND,tmp)
ONDDemographics_SampleConverted = FinalOMNDDemographics

for(i in 1:nrow(ONDDemographics_SampleConverted)){
  for(j in 1:nrow(FinalOMNDDemographics)){
    
    if(! is.na(ONDDemographics_SampleConverted$RNAseq1[i])){
      if(ONDDemographics_SampleConverted$RNAseq1[i] == FinalOMNDDemographics$ExternalSampleId[j]){
        ONDDemographics_SampleConverted$RNAseq1[i] = FinalOMNDDemographics$Sample.Source[j]
      }
    }
    
    if(! is.na(ONDDemographics_SampleConverted$RNAseq2[i])){
      if(ONDDemographics_SampleConverted$RNAseq2[i] == FinalOMNDDemographics$ExternalSampleId[j]){
        ONDDemographics_SampleConverted$RNAseq2[i] = FinalOMNDDemographics$Sample.Source[j]
      }
    }
    
    if(! is.na(ONDDemographics_SampleConverted$RNAseq3[i])){
      if(ONDDemographics_SampleConverted$RNAseq3[i] == FinalOMNDDemographics$ExternalSampleId[j]){
        ONDDemographics_SampleConverted$RNAseq3[i] = FinalOMNDDemographics$Sample.Source[j]
      }
    }
    
  }
  if((i %% 10) == 0) cat("% Done:",i/nrow(HCDemographics_SampleConverted)*100,"\n")
}

#Sex
table(ONDDemographics_SampleConverted$Sex)
#Tissue Site
tmpcontainer = c(ONDDemographics_SampleConverted$RNAseq1,ONDDemographics_SampleConverted$RNAseq2,ONDDemographics_SampleConverted$RNAseq3)
table(tmpcontainer)
#Age at Death
mean(as.numeric(ONDDemographics_SampleConverted$Age.at.Death))
sd(as.numeric(ONDDemographics_SampleConverted$Age.at.Death))
#Disease Duration - Not Available
#mean(as.numeric(ONDDemographics_SampleConverted$Disease.Duration.in.Months))
#sd(as.numeric(ONDDemographics_SampleConverted$Disease.Duration.in.Months))

###############################################################################################################################
###############################################################################################################################

################################################# ALS Site of onset (All categories)
#Without Discordant
SiteofOnset = FinalSurvivalMat
for(i in 1:nrow(SiteofOnset)){
  if(SiteofOnset$Site.of.Motor.Onset[i] == "Not Applicable"){
    SiteofOnset$Site.of.Motor.Onset[i] = "Unknown"
  }
}

tmpdata = SiteofOnset

ite = which(tmpdata$Final == "TE")
iox = which(tmpdata$Final == "OX")
igl = which(tmpdata$Final == "GLIA")

tmpTE = tmpdata[ite,]
tmpOX = tmpdata[iox,]
tmpGL = tmpdata[igl,]

TEsoo = tmpTE$Site.of.Motor.Onset
OXsoo = tmpOX$Site.of.Motor.Onset
GLsoo = tmpGL$Site.of.Motor.Onset

a = table(TEsoo)
b = table(OXsoo)
c = table(GLsoo)

barplot(a,main="ALS-TE Site of Onset",names.arg = names(a),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="firebrick",ylim=c(0,60))
barplot(b,main="ALS-Ox Site of Onset",names.arg = names(b),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="navy",ylim=c(0,60))
barplot(c,main="ALS-Glia Site of Onset",names.arg = names(c),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="goldenrod",ylim=c(0,60))

#Chi-squared test of independence
Categories = names(table(names(c(a,b,c))))
CT = matrix(NA,nrow=3,ncol = length(Categories))
rownames(CT) = c("ALS-TE","ALS-OX","ALS-GLIA")
colnames(CT) = Categories

for(i in 1:length(a)){
  for(j in 1:ncol(CT)){
    if(names(a)[i] == colnames(CT)[j]){
      CT[1,j] = a[[i]]
    }
  }
}
for(i in 1:length(b)){
  for(j in 1:ncol(CT)){
    if(names(b)[i] == colnames(CT)[j]){
      CT[2,j] = b[[i]]
    }
  }
}
for(i in 1:length(c)){
  for(j in 1:ncol(CT)){
    if(names(c)[i] == colnames(CT)[j]){
      CT[3,j] = c[[i]]
    }
  }
}

for(i in 1:nrow(CT)){
  for(j in 1:ncol(CT)){
    if(is.na(CT[i,j])){
      CT[i,j] = 0
    }
  }
}


chisq.test(CT) # p = 0.745 
#Do not reject the null that Site of Onset is independent of ALS Subtype



#With Discordant
SiteofOnset = FinalSurvivalMat
for(i in 1:nrow(SiteofOnset)){
  if(SiteofOnset$Site.of.Motor.Onset[i] == "Not Applicable"){
    SiteofOnset$Site.of.Motor.Onset[i] = "Unknown"
  }
}

tmpdata = SiteofOnset

ite = which(tmpdata$Final == "TE")
iox = which(tmpdata$Final == "OX")
igl = which(tmpdata$Final == "GLIA")
idi = which(tmpdata$Final == "Discordant")

tmpTE = tmpdata[ite,]
tmpOX = tmpdata[iox,]
tmpGL = tmpdata[igl,]
tmpDI = tmpdata[idi,]

TEsoo = tmpTE$Site.of.Motor.Onset
OXsoo = tmpOX$Site.of.Motor.Onset
GLsoo = tmpGL$Site.of.Motor.Onset
DIsoo = tmpDI$Site.of.Motor.Onset

a = table(TEsoo)
b = table(OXsoo)
c = table(GLsoo)
d = table(DIsoo)

barplot(a,main="ALS-TE Site of Onset",names.arg = names(a),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="firebrick",ylim=c(0,60))
barplot(b,main="ALS-Ox Site of Onset",names.arg = names(b),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="navy",ylim=c(0,60))
barplot(c,main="ALS-Glia Site of Onset",names.arg = names(c),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="goldenrod",ylim=c(0,60))
barplot(d,main="ALS-Discordant Site of Onset",names.arg = names(d),cex.axis = 1.25,cex.names = 1.25,cex.lab=1.5,xlab="Subtype",ylab="Number of Patients",col="darkorchid2",ylim=c(0,60))

#Chi-squared test of independence
Categories = names(table(names(c(a,b,c,d))))
CT = matrix(NA,nrow=4,ncol = length(Categories))
rownames(CT) = c("ALS-TE","ALS-OX","ALS-GLIA","ALS-Discordant")
colnames(CT) = Categories

for(i in 1:length(a)){
  for(j in 1:ncol(CT)){
    if(names(a)[i] == colnames(CT)[j]){
      CT[1,j] = a[[i]]
    }
  }
}
for(i in 1:length(b)){
  for(j in 1:ncol(CT)){
    if(names(b)[i] == colnames(CT)[j]){
      CT[2,j] = b[[i]]
    }
  }
}
for(i in 1:length(c)){
  for(j in 1:ncol(CT)){
    if(names(c)[i] == colnames(CT)[j]){
      CT[3,j] = c[[i]]
    }
  }
}
for(i in 1:length(d)){
  for(j in 1:ncol(CT)){
    if(names(d)[i] == colnames(CT)[j]){
      CT[4,j] = d[[i]]
    }
  }
}

for(i in 1:nrow(CT)){
  for(j in 1:ncol(CT)){
    if(is.na(CT[i,j])){
      CT[i,j] = 0
    }
  }
}


chisq.test(CT) # p = 0.822
#Do not reject the null that Site of Onset is independent of ALS Subtype
