### ALS Classification - Patel et al. Subtype Scoring Method
# Reference: https://www.science.org/doi/full/10.1126/science.1254257

#Jarrett Eshima
#Dr. Barbara Smith Lab

#Load in dependent libraries
library(scatterplot3d)
library(plot3D)
library(rgl)

#####################################################################################################################################################################################

#Load in cleaned/normalized count matrix (this file is generated in the UnivariateAnalysis.R script)
load("D:/Jarrett/Research/Fall 2021/ProgBM/SOD1/ALSPatientStratification_UnivariateDatasets_SOD1.RData")

#Read in phenotype data
wd2 = "C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Supervised Classification/ALS_Patel"
setwd(wd2)
Pheno = read.csv("ALS451_coldata_SUBTYPES.csv")
SampleList = Pheno$Subject

############# Step 1: Define subtype predictors and visualize preliminary subtype scores ###############################

#Purple Eigengene - Glia; Magenta Eigengene - TD; Turquoise Eigengene - Ox
GliaPredictors = read.csv("GliaFeatures_Purple.csv") #Derived from table S8
GliaPredictors = unlist(GliaPredictors)
table(table(GliaPredictors)>1)#duplicate check
OxPredictors = read.csv("OxFeatures_Turq.csv") #Derived from table S8
OxPredictors = unlist(OxPredictors)
table(table(OxPredictors)>1)#duplicate check
TDPredictors = read.csv("TDFeatures_Mag.csv") #Derived from table S8
TDPredictors = unlist(TDPredictors)
table(table(TDPredictors)>1)#duplicate check

#Features used in more than one subtype? - No
GliaPredictors[GliaPredictors %in% OxPredictors]
GliaPredictors[GliaPredictors %in% TDPredictors]
OxPredictors[OxPredictors %in% TDPredictors]

GliaAverage = NormCounts[rownames(NormCounts) %in% GliaPredictors,]
OxAverage = NormCounts[rownames(NormCounts) %in% OxPredictors,]
TDAverage = NormCounts[rownames(NormCounts) %in% TDPredictors,]

GliaAverage = data.frame(GliaAverage);OxAverage = data.frame(OxAverage);TDAverage = data.frame(TDAverage)

ALSGlia = GliaAverage[,colnames(GliaAverage) %in% SampleList]
ALSOx = OxAverage[,colnames(OxAverage) %in% SampleList]
ALSTD = TDAverage[,colnames(TDAverage) %in% SampleList]

ALSGlia = data.frame(t(ALSGlia))
ALSOx = data.frame(t(ALSOx))
ALSTD = data.frame(t(ALSTD))

AverageGliaScore = rowMeans(ALSGlia,na.rm=T)
AverageOxScore = rowMeans(ALSOx,na.rm=T)
AverageTDScore = rowMeans(ALSTD,na.rm=T)

#Average Expression of 1690 features
CleanNormCounts = NormCounts[,colnames(NormCounts) %in% SampleList] #1690 features
CleanNormCounts = t(CleanNormCounts)
dim(CleanNormCounts)
table(rownames(ALSGlia) == rownames(CleanNormCounts))
AverageExpression = rowMeans(CleanNormCounts,na.rm = T)

#Patel et al. defined subtype scores
SubtypeScore_Glia = AverageGliaScore - AverageExpression
SubtypeScore_Ox = AverageOxScore - AverageExpression
SubtypeScore_TD = AverageTDScore - AverageExpression

SubtypeScores = data.frame(matrix(NA,length(SubtypeScore_Glia),3))
rownames(SubtypeScores) = SampleList
colnames(SubtypeScores) = c("Glia","Ox","TD")

SubtypeScores$Glia = SubtypeScore_Glia
SubtypeScores$Ox = SubtypeScore_Ox
SubtypeScores$TD = SubtypeScore_TD

#Assign colors based on unsupervised clustering results
subtypecolors = shapes = rep(NA,length(SampleList))
for(i in 1:length(SampleList)){
  
  tmp = SampleList[i]
  ind = which(Pheno$Subject == tmp)
  tmp2 = Pheno$Subtype[ind]
  
  if(tmp2 == "TE"){
    subtypecolors[i] = "firebrick"
    shapes[i] = 15
  }else if(tmp2 == "OX"){
    subtypecolors[i] = "navy"
    shapes[i] = 16
  }else{
    subtypecolors[i] = "goldenrod"
    shapes[i] = 17
  }
  
}

#Visualize subtype scoring method
plot3d(SubtypeScores[,1],SubtypeScores[,2],SubtypeScores[,3],xlab = "Glia Score",ylab = "Ox Score",zlab="TD Score",col=subtypecolors,size = 7)

#################################### Step 2: Define 5% cutoff ###################################################################################

set.seed(1234)
times = 100
npredictorsG = length(GliaPredictors)
npredictorsO = length(OxPredictors)
npredictorsT = length(TDPredictors)


RandomPredictors_Glia = data.frame(matrix(NA,times,npredictorsG))
RandomPredictors_Ox = data.frame(matrix(NA,times,npredictorsO))
RandomPredictors_TD = data.frame(matrix(NA,times,npredictorsT))
FeatureList = rownames(NormCounts)

# #Using random features from 1690 set - All thresholds the same (roughly)
# for(i in 1:times){
#   Gindices = sample(1:nrow(NormCounts),npredictorsG,replace = F)
#   RandomPredictors_Glia[i,] = FeatureList[Gindices]
#   Oindices = sample(1:nrow(NormCounts),npredictorsO,replace = F)
#   RandomPredictors_Ox[i,] = FeatureList[Oindices]
#   Tindices = sample(1:nrow(NormCounts),npredictorsT,replace = F)
#   RandomPredictors_TD[i,] = FeatureList[Tindices]
# }

#Using sampling with replacement from subtype predictor lists ("5% cutoff")
for(i in 1:times){
  Gindices = sample(1:length(GliaPredictors),npredictorsG,replace = T)
  RandomPredictors_Glia[i,] = GliaPredictors[Gindices]
  Oindices = sample(1:length(OxPredictors),npredictorsO,replace = T)
  RandomPredictors_Ox[i,] = OxPredictors[Oindices]
  Tindices = sample(1:length(TDPredictors),npredictorsT,replace = T)
  RandomPredictors_TD[i,] = TDPredictors[Tindices]
}


#Iterate over 100 different sets of subtype predictors (according to Patel et al.)
ScoreSamplingG = ScoreSamplingO = ScoreSamplingT = data.frame(matrix(NA,length(SampleList),times))
rownames(ScoreSamplingG) = rownames(ScoreSamplingO) = rownames(ScoreSamplingT) = SampleList
colnames(ScoreSamplingG) = colnames(ScoreSamplingO) = colnames(ScoreSamplingT) = paste("Rep",seq(1,times,1),sep="")

for(i in 1:times){
  
  tmp = RandomPredictors_Glia[i,]
  RSAverage = NormCounts[rownames(NormCounts) %in% tmp,]
  RSAverage = data.frame(RSAverage)
  RSAverage = RSAverage[,colnames(RSAverage) %in% SampleList]
  RSAverage = data.frame(t(RSAverage))
  RandomScore = rowMeans(RSAverage,na.rm=T)
  SubtypeScore_cutoff = RandomScore - AverageExpression
  values = unname(SubtypeScore_cutoff)
  ScoreSamplingG[,i] = values
  
  tmp = RandomPredictors_Ox[i,]
  RSAverage = NormCounts[rownames(NormCounts) %in% tmp,]
  RSAverage = data.frame(RSAverage)
  RSAverage = RSAverage[,colnames(RSAverage) %in% SampleList]
  RSAverage = data.frame(t(RSAverage))
  RandomScore = rowMeans(RSAverage,na.rm=T)
  SubtypeScore_cutoff = RandomScore - AverageExpression
  values = unname(SubtypeScore_cutoff)
  ScoreSamplingO[,i] = values
  
  tmp = RandomPredictors_TD[i,]
  RSAverage = NormCounts[rownames(NormCounts) %in% tmp,]
  RSAverage = data.frame(RSAverage)
  RSAverage = RSAverage[,colnames(RSAverage) %in% SampleList]
  RSAverage = data.frame(t(RSAverage))
  RandomScore = rowMeans(RSAverage,na.rm=T)
  SubtypeScore_cutoff = RandomScore - AverageExpression
  values = unname(SubtypeScore_cutoff)
  ScoreSamplingT[,i] = values
  
}

#Optional code to check normality for use of z distribution in selection of 5% cutoff
# normcheck = sample(1:times,9,replace=F) #Approximately normally distributed for most replicates/iterations
# par(mfrow=c(3,3))
# for(i in 1:length(normcheck)){
#   tmp = normcheck[i]
#   hist(ScoreSamplingG[,tmp])
# }
# for(i in 1:length(normcheck)){
#   tmp = normcheck[i]
#   hist(ScoreSamplingO[,tmp])
# }
# for(i in 1:length(normcheck)){
#   tmp = normcheck[i]
#   hist(ScoreSamplingT[,tmp])
# }
# par(mfrow=c(1,1))

ScoreSamplingG2 = matrix(unlist(ScoreSamplingG),nrow(ScoreSamplingG),ncol(ScoreSamplingG))
rownames(ScoreSamplingG2) = rownames(ScoreSamplingG)
colnames(ScoreSamplingG2) = colnames(ScoreSamplingG)
ScoreSamplingO2 = matrix(unlist(ScoreSamplingO),nrow(ScoreSamplingO),ncol(ScoreSamplingO))
rownames(ScoreSamplingO2) = rownames(ScoreSamplingO)
colnames(ScoreSamplingO2) = colnames(ScoreSamplingO)
ScoreSamplingT2 = matrix(unlist(ScoreSamplingT),nrow(ScoreSamplingT),ncol(ScoreSamplingT))
rownames(ScoreSamplingT2) = rownames(ScoreSamplingT)
colnames(ScoreSamplingT2) = colnames(ScoreSamplingT)


#Average Across Predictor Sets
repthreshG = repthreshO = repthreshT = rep(NA,nrow(ScoreSamplingG))
for(i in 1:nrow(ScoreSamplingG2)){
  
  tmpmean = mean(ScoreSamplingG2[i,])
  tmpsigma = sd(ScoreSamplingG2[i,])
  repthreshG[i] = tmpmean + 1.645*tmpsigma
  
  tmpmean = mean(ScoreSamplingO2[i,])
  tmpsigma = sd(ScoreSamplingO2[i,])
  repthreshO[i] = tmpmean + 1.645*tmpsigma
  
  tmpmean = mean(ScoreSamplingT2[i,])
  tmpsigma = sd(ScoreSamplingT2[i,])
  repthreshT[i] = tmpmean + 1.645*tmpsigma
}

#Thresholds values are with respect to unsupervised clustering subtype proportions
table(Pheno$Subtype)
Gtmp = repthreshG[order(repthreshG,decreasing = T)][84] #84 of 451 are in the Glia cluster
Otmp = repthreshO[order(repthreshO,decreasing = T)][239] #239 of 451 are in the Ox cluster
Ttmp = repthreshT[order(repthreshT,decreasing = T)][128] #128 of 451 are in the TD cluster

#Visualize thresholds
par(mfrow=c(1,3))
hist(repthreshG); abline(v = Gtmp,col="red")
hist(repthreshO); abline(v = Otmp,col="red")
hist(repthreshT); abline(v = Ttmp,col="red")

#Initialize threshold reference container
SubtypeThresholds = data.frame(matrix(NA,nrow(ScoreSamplingG2),4))
colnames(SubtypeThresholds) = c("SampleID","GThresh","OThresh","TDThresh")
SubtypeThresholds$SampleID = rownames(ScoreSamplingG2)
SubtypeThresholds$GThresh = Gtmp
SubtypeThresholds$OThresh = Otmp
SubtypeThresholds$TDThresh = Ttmp


#Pre-calculate average Expression of 1690 features
CleanNormCounts = NormCounts[,colnames(NormCounts) %in% SampleList] #1690 features
CleanNormCounts = t(CleanNormCounts)
table(rownames(ALSGlia) == rownames(CleanNormCounts))
AverageExpression = rowMeans(CleanNormCounts,na.rm = T)

#Filter out control samples
ALSCounts = NormCounts[,colnames(NormCounts) %in% SampleList]
table(colnames(ALSCounts) == SubtypeThresholds$SampleID)

############################### Step 3: Calculate subtype scores ###################################################################################

iter = 1000 #n iterations in bootstrap
GliaScores = OxScores = TDScores = data.frame(matrix(NA,nrow(SubtypeThresholds),iter)) #containers
rownames(GliaScores) = rownames(OxScores) = rownames(TDScores) = SubtypeThresholds$SampleID

for(i in 1:iter){
  sampGP = sample(GliaPredictors,npredictorsG,replace = T)
  sampOP = sample(OxPredictors,npredictorsO,replace = T)
  sampTP = sample(TDPredictors,npredictorsT,replace = T)
  
  GliaAverage = data.frame(matrix(NA,nrow = length(sampGP),ncol = ncol(ALSCounts)))
  OxAverage = data.frame(matrix(NA,nrow = length(sampOP),ncol = ncol(ALSCounts)))
  TDAverage = data.frame(matrix(NA,nrow = length(sampTP),ncol = ncol(ALSCounts)))
  
  for(j in 1:length(sampGP)){
    index = which(rownames(ALSCounts) == sampGP[j])
    GliaAverage[j,] = ALSCounts[index,]
  }
  
  for(k in 1:length(sampOP)){
    index = which(rownames(ALSCounts) == sampOP[k])
    OxAverage[k,] = ALSCounts[index,]
  }
  
  for(l in 1:length(sampTP)){
    index = which(rownames(ALSCounts) == sampTP[l])
    TDAverage[l,] = ALSCounts[index,]
  }
  
  ALSGlia = data.frame(t(GliaAverage))
  ALSOx = data.frame(t(OxAverage))
  ALSTD = data.frame(t(TDAverage))
  
  AverageGliaScore = rowMeans(ALSGlia,na.rm=T)
  AverageOxScore = rowMeans(ALSOx,na.rm=T)
  AverageTDScore = rowMeans(ALSTD,na.rm=T)
  
  #Patel et al. defined subtype scores
  SubtypeScore_Glia = AverageGliaScore - AverageExpression
  SubtypeScore_Ox = AverageOxScore - AverageExpression
  SubtypeScore_TD = AverageTDScore - AverageExpression
  
  GliaScores[,i] = SubtypeScore_Glia
  OxScores[,i] = SubtypeScore_Ox
  TDScores[,i] = SubtypeScore_TD
  
  if((i %% 10) == 0) cat("% Done:",i/iter*100,"\n")
}

#Assess sample expression in the context of subtype thresholds
table(rownames(GliaScores) == SubtypeThresholds$SampleID) #Must be all true
BootClassification = data.frame(matrix(NA,nrow(SubtypeThresholds),3))
colnames(BootClassification) = c("GliaScore","OxScore","TDScore")
rownames(BootClassification) = SubtypeThresholds$SampleID
for(i in 1:length(AverageExpression)){
  GliaCounter = OxCounter = TDCounter = 0
  for(j in 1:iter){
    
    if(GliaScores[i,j] > SubtypeThresholds$GThresh[i]){
      GliaCounter = GliaCounter + 1
    }
    
    if(OxScores[i,j] > SubtypeThresholds$OThresh[i]){
      OxCounter = OxCounter + 1
    }
    
    if(TDScores[i,j] > SubtypeThresholds$TDThresh[i]){
      TDCounter = TDCounter + 1
    }
    
  }
  
  BootClassification$GliaScore[i] = GliaCounter/iter
  BootClassification$OxScore[i] = OxCounter/iter
  BootClassification$TDScore[i] = TDCounter/iter
  if((i %% 10) == 0) cat("% Done:",i/length(AverageExpression)*100,"\n")
}

plot3d(BootClassification$GliaScore,BootClassification$OxScore,BootClassification$TDScore,xlab = "Glia Score",ylab = "Ox Score",zlab="TD Score",col=subtypecolors,size = 7)

############################### Step 4: Establish color code ###################################################################################

#Color Code for Scoring-based Classification (similar to Patel et al.)
classcolors = rep(NA,nrow(BootClassification))
for(i in 1:length(classcolors)){
  
  if(BootClassification$GliaScore[i] > 0.5){
    if(BootClassification$OxScore[i] <= 1 && BootClassification$OxScore[i] > 0.4){
      classcolors[i] = "chartreuse2"
    }else if(BootClassification$TDScore[i] <= 1 && BootClassification$TDScore[i] > 0.4){
      classcolors[i] = "darkorange2"
    }else if(BootClassification$OxScore[i] < 0.4 && BootClassification$TDScore[i] < 0.4){  
      classcolors[i] = "goldenrod"
    }else{
      classcolors[i] = "gray50"
    }
  }else if(BootClassification$OxScore[i] > 0.5){
    if(BootClassification$GliaScore[i] <= 1 && BootClassification$GliaScore[i] > 0.4){
      classcolors[i] = "chartreuse2"
    }else if(BootClassification$TDScore[i] <= 1 && BootClassification$TDScore[i] > 0.4){
      classcolors[i] = "darkorchid2"
    }else if(BootClassification$GliaScore[i] < 0.4 && BootClassification$TDScore[i] < 0.4){
      classcolors[i] = "navy"
    }else{
      classcolors[i] = "gray50"
    }
  }else if(BootClassification$TDScore[i] > 0.5){
    if(BootClassification$GliaScore[i] <= 1 && BootClassification$GliaScore[i] > 0.4){
      classcolors[i] = "darkorange2"
    }else if(BootClassification$OxScore[i] <= 1  && BootClassification$OxScore[i] > 0.4){
      classcolors[i] = "darkorchid2"
    }else if(BootClassification$GliaScore[i] < 0.4 && BootClassification$OxScore[i] < 0.4){
      classcolors[i] = "firebrick"
    }else{
      classcolors[i] = "gray50"
    }
  }else if(BootClassification$GliaScore[i] <= 1 && BootClassification$GliaScore[i] > 0.4){
    if(BootClassification$OxScore[i] <= 1 && BootClassification$OxScore[i] > 0.4){
      classcolors[i] = "chartreuse2"
    }else if(BootClassification$TDScore[i] <= 1 && BootClassification$TDScore[i] > 0.4){
      classcolors[i] = "darkorange2"
    }else{
      classcolors[i] = "gray50"
    }
  }else if(BootClassification$OxScore[i] <= 1 && BootClassification$OxScore[i] > 0.4){
    if(BootClassification$GliaScore[i] <= 1 && BootClassification$GliaScore[i] > 0.4){
      classcolors[i] = "chartreuse2"
    }else if(BootClassification$TDScore[i] <= 1 && BootClassification$TDScore[i] > 0.4){
      classcolors[i] = "darkorchid2"
    }else{
      classcolors[i] = "gray50"
    }
  }else if(BootClassification$TDScore[i] <= 1 && BootClassification$TDScore[i] > 0.4){
    if(BootClassification$GliaScore[i] <= 1 && BootClassification$GliaScore[i] > 0.4){
      classcolors[i] = "darkorange2"
    }else if(BootClassification$OxScore[i] <= 1  && BootClassification$OxScore[i] > 0.4){
      classcolors[i] = "darkorchid2"
    }else{
      classcolors[i] = "gray50"
    }
  }else{ 
    classcolors[i] = "gray50"
  }
  
}

############################################## Step 5: Plot ###################################################################################


#Point Style - classcolors
plot3d(BootClassification$GliaScore,BootClassification$TDScore,BootClassification$OxScore,xlab = "",zlab = "",ylab="",col=classcolors,pch=20,size=7,axes=F)
#Sphere Style - classcolors
plot3d(BootClassification$GliaScore,BootClassification$TDScore,BootClassification$OxScore,type = "s",radius = 0.02,xlab = "",zlab = "",ylab="",col=classcolors,pch=20,size=7,axes=F)
#Sphere Style - subtypecolors
plot3d(BootClassification$GliaScore,BootClassification$TDScore,BootClassification$OxScore,type = "s",radius = 0.02,xlab = "",zlab = "",ylab="",col=subtypecolors,pch=20,size=7,axes=F)


######## FINAL PLOT
#Plot Subtype Scores
plot3d(BootClassification$GliaScore,BootClassification$TDScore,BootClassification$OxScore,xlab = "",zlab = "",ylab="",col=classcolors,pch=20,size=7,axes=F)
#Add bounding grid
axis3d("x",pos = c(0,0,0),tick = F,labels = F)
axis3d("y",pos = c(0,0,0),tick = F,labels = F)
axis3d("z",pos = c(0,0,0),tick = F,labels = F)
axis3d("x",pos = c(1,1,0),tick = F,labels = F)
axis3d("x",pos = c(1,0,1),tick = F,labels = F)
axis3d("y",pos = c(1,1,0),tick = F,labels = F)
axis3d("y",pos = c(0,1,1),tick = F,labels = F)
segments3d(x=c(1,1),y=c(0,0),z=c(0,1))
segments3d(x=c(0,0),y=c(1,1),z=c(0,1))
#Overlay Clustering information 
pch3d(BootClassification$GliaScore,BootClassification$TDScore,BootClassification$OxScore,pch = 21, bg = classcolors,color = subtypecolors,cex=0.1)
#Save Plot
rgl.postscript('Patel_ALSSubtypeScore_3DScatter_Eigen2.pdf',fmt = 'pdf')


#Empty Grid for Legend
plot3d()
axis3d("x",pos = c(0,0,0),tick = F,labels = F)
axis3d("y",pos = c(0,0,0),tick = F,labels = F)
axis3d("z",pos = c(0,0,0),tick = F,labels = F)
axis3d("x",pos = c(1,1,0),tick = F,labels = F)
axis3d("x",pos = c(1,0,1),tick = F,labels = F)
axis3d("y",pos = c(1,1,0),tick = F,labels = F)
axis3d("y",pos = c(0,1,1),tick = F,labels = F)
segments3d(x=c(1,1),y=c(0,0),z=c(0,1))
segments3d(x=c(0,0),y=c(1,1),z=c(0,1))


rgl.postscript('Patel_ALSSubtypeScore_3DScatter_Grid.pdf',fmt = 'pdf')
      
######################## Step 6: Clustering Visualization ###################################################################################

par(mfrow=c(1,1))

#Sample order is the same for all three
classcolors
subtypecolors
BootClassification

table(rownames(BootClassification) == Pheno$Subject[-remsub])

ClustSubtype = ALSPheno$Subtype[-remsub]

for(i in 1:length(ClustSubtype)){
  
  if(ClustSubtype[i] == "GLIA"){
    ClustSubtype[i] = "ALS-Glia"
  }else if(ClustSubtype[i] == "OX"){
    ClustSubtype[i] = "ALS-Ox"
  }else if(ClustSubtype[i] == "TE"){
    ClustSubtype[i] = "ALS-TD"
  }
  
}

#ALS-Glia / ALS-TD Hybrids
GliaTDind = which(classcolors == "darkorange2")
GliaTDScores = BootClassification[GliaTDind,]
GliaTDSamples = rownames(GliaTDScores)
GliaTDClust = ClustSubtype[GliaTDind]
GliaTDbar = table(factor(GliaTDClust,levels = c("ALS-Glia","ALS-Ox","ALS-TD")))
barplot(GliaTDbar,col=c("goldenrod1","navy","firebrick"),main = "ALS-Glia + ALS-TD Hybrids: Clustering Subtypes",ylab="Frequency",ylim=c(0,16),cex.axis = 1.5,cex = 1.5,cex.main = 1.5,cex.lab=1.5)

#ALS-Glia / ALS-Ox Hybrids
GliaOxind = which(classcolors == "chartreuse2")
GliaOxScores = BootClassification[GliaOxind,]
GliaOxSamples = rownames(GliaOxScores)
GliaOxClust = ClustSubtype[GliaOxind]
GliaOxbar = table(factor(GliaOxClust,levels = c("ALS-Glia","ALS-Ox","ALS-TD")))
barplot(GliaOxbar,col=c("goldenrod1","navy","firebrick"),main = "ALS-Glia + ALS-Ox Hybrids: Clustering Subtypes",ylab="Frequency",ylim=c(0,16),cex.axis = 1.5,cex = 1.5,cex.main = 1.5,cex.lab=1.5)


#ALS-Glia
Gliaind = which(classcolors == "goldenrod")
GliaSamples = rownames(BootClassification)[Gliaind]
GliaClust = ClustSubtype[Gliaind]
Gliabar = table(factor(GliaClust,levels = c("ALS-Glia","ALS-Ox","ALS-TD")))
barplot(Gliabar,col=c("goldenrod1","navy","firebrick"),main = "ALS-Glia: Clustering Subtypes",ylab="Frequency",ylim=c(0,16),cex.axis = 1.5,cex = 1.5,cex.main = 1.5,cex.lab=1.5)
table(GliaClust)

#ALS-Ox
Oxind = which(classcolors == "navy")
OxSamples = rownames(BootClassification)[Oxind]
OxClust = ClustSubtype[Oxind]
Oxbar = table(factor(OxClust,levels = c("ALS-Glia","ALS-Ox","ALS-TD")))
barplot(Oxbar,col=c("goldenrod1","navy","firebrick"),main = "ALS-Ox: Clustering Subtypes",ylab="Frequency",ylim=c(0,120),cex.axis = 1.5,cex = 1.5,cex.main = 1.5,cex.lab=1.5)
table(OxClust)

#ALS-TD
TDind = which(classcolors == "firebrick")
TDSamples = rownames(BootClassification)[TDind]
TDClust = ClustSubtype[TDind]
TDbar = table(factor(TDClust,levels = c("ALS-Glia","ALS-Ox","ALS-TD")))
barplot(TDbar,col=c("goldenrod1","navy","firebrick"),main = "ALS-TD: Clustering Subtypes",ylab="Frequency",ylim=c(0,60),cex.axis = 1.5,cex = 1.5,cex.main = 1.5,cex.lab=1.5)
table(TDClust)

#Unclassified
Uncind = which(classcolors == "gray50")
UncSamples = rownames(BootClassification)[Uncind]
UncClust = ClustSubtype[Uncind]
Uncbar = table(factor(UncClust,levels = c("ALS-Glia","ALS-Ox","ALS-TD")))
barplot(Uncbar,col=c("goldenrod1","navy","firebrick"),main = "Unclassified: Clustering Subtypes",ylab="Frequency",ylim=c(0,120),cex.axis = 1.5,cex = 1.5,cex.main = 1.5,cex.lab=1.5)


################################ Step 7: Clinical Assessment ###################################################################################

### Survival
library(survival)
library(survminer)
#Survival in Tam cohort
setwd("C:/Users/jeshima/Desktop/ALS Manuscript/Nature/Communications Submission/Peer Review")

SurvivalDat = read.csv("Table_S10.csv")

HybridSurvival = SurvivalDat
HybridSurvival$RNAseqSample1 = gsub("-","\\.",HybridSurvival$RNAseqSample1)
HybridSurvival$RNAseqSample2 = gsub("-","\\.",HybridSurvival$RNAseqSample2)
HybridSurvival$RNAseqSample3 = gsub("-","\\.",HybridSurvival$RNAseqSample3)


#ALS-Glia / ALS-TD
HybridSamples = GliaTDSamples

hybpatients1 = hybpatients2 = hybpatients3 = rep(NA,nrow(HybridSurvival))
count1 = count2 = count3 = 1
for(i in 1:nrow(HybridSurvival)){
  
  for(j in 1:length(HybridSamples)){
    if(HybridSurvival$RNAseqSample1[i]==HybridSamples[j]){
      hybpatients1[count1] = i
      count1 = count1+1
    }
    
    if(! is.na(HybridSurvival$RNAseqSample2[i])){
      if(HybridSurvival$RNAseqSample2[i] == HybridSamples[j]){
        hybpatients2[count2] = i
        count2 = count2+1
      }
    }
    
    if(! is.na(HybridSurvival$RNAseqSample3[i])){
      if(HybridSurvival$RNAseqSample3[i] == HybridSamples[j]){
        hybpatients3[count3] = i
        count3 = count3+1
      }
    }
  }
  
}

GliaTDIndex = as.numeric(names(table(c(hybpatients1,hybpatients2,hybpatients3))))

GliaTDSurvival = HybridSurvival[GliaTDIndex,]
mean(GliaTDSurvival$Duration)
sd(GliaTDSurvival$Duration)/sqrt(length(GliaTDSurvival$Duration))


for(i in 1:length(GliaTDSurvival$FinalSubtype)){
  
  if(GliaTDSurvival$FinalSubtype[i] == "GLIA"){
    GliaTDSurvival$FinalSubtype[i] = "ALS-Glia"
  }else if(GliaTDSurvival$FinalSubtype[i] == "OX"){
    GliaTDSurvival$FinalSubtype[i] = "ALS-Ox"
  }else if(GliaTDSurvival$FinalSubtype[i] == "TE"){
    GliaTDSurvival$FinalSubtype[i] = "ALS-TD"
  }
}

#Patient-wise subtypes
PatientGliaTD = table(factor(GliaTDSurvival$FinalSubtype,levels = c("ALS-Glia","ALS-Ox","ALS-TD","Discordant")))
barplot(PatientGliaTD,col=c("goldenrod1","navy","firebrick","darkorchid2"),main = "ALS-Glia + ALS-TD Hybrid: Patient Subtypes",ylab="Frequency",ylim=c(0,12),cex.axis = 1.5,cex = 1.5,cex.main = 1.5,cex.lab=1.5)

HybridSurvival$FinalSubtype[GliaTDIndex] = "ALS-Glia + ALS-TD Hybrid" #Patients with 1 or more hybrid samples are assigned the hybrid subtype

HybridSamples = GliaOxSamples

hybpatients1 = hybpatients2 = hybpatients3 = rep(NA,nrow(HybridSurvival))
count1 = count2 = count3 = 1
for(i in 1:nrow(HybridSurvival)){
  
  for(j in 1:length(HybridSamples)){
    if(HybridSurvival$RNAseqSample1[i]==HybridSamples[j]){
      hybpatients1[count1] = i
      count1 = count1+1
    }
    
    if(! is.na(HybridSurvival$RNAseqSample2[i])){
      if(HybridSurvival$RNAseqSample2[i] == HybridSamples[j]){
        hybpatients2[count2] = i
        count2 = count2+1
      }
    }
    
    if(! is.na(HybridSurvival$RNAseqSample3[i])){
      if(HybridSurvival$RNAseqSample3[i] == HybridSamples[j]){
        hybpatients3[count3] = i
        count3 = count3+1
      }
    }
  }
  
}

GliaOxIndex = as.numeric(names(table(c(hybpatients1,hybpatients2,hybpatients3))))

GliaOxSurvival = HybridSurvival[GliaOxIndex,]
mean(GliaOxSurvival$Duration)
sd(GliaOxSurvival$Duration)/sqrt(length(GliaOxSurvival$Duration))


for(i in 1:length(GliaOxSurvival$FinalSubtype)){
  
  if(GliaOxSurvival$FinalSubtype[i] == "GLIA"){
    GliaOxSurvival$FinalSubtype[i] = "ALS-Glia"
  }else if(GliaOxSurvival$FinalSubtype[i] == "OX"){
    GliaOxSurvival$FinalSubtype[i] = "ALS-Ox"
  }else if(GliaOxSurvival$FinalSubtype[i] == "TE"){
    GliaOxSurvival$FinalSubtype[i] = "ALS-TD"
  }
}

#Patient-wise subtypes
PatientGliaOx = table(factor(GliaOxSurvival$FinalSubtype,levels = c("ALS-Glia","ALS-Ox","ALS-TD","Discordant")))
barplot(PatientGliaOx,col=c("goldenrod1","navy","firebrick","darkorchid2"),main = "ALS-Glia + ALS-Ox Hybrid: Patient Subtypes",ylab="Frequency",ylim=c(0,12),cex.axis = 1.5,cex = 1.5,cex.main = 1.5,cex.lab=1.5)

HybridSurvival$FinalSubtype[GliaOxIndex] = "ALS-Glia + ALS-Ox Hybrid" #Patients with 1 or more hybrid samples are assigned the hybrid subtype

for(i in 1:length(HybridSurvival$FinalSubtype)){
  
  if(HybridSurvival$FinalSubtype[i] == "GLIA"){
    HybridSurvival$FinalSubtype[i] = "ALS-Glia"
  }else if(HybridSurvival$FinalSubtype[i] == "OX"){
    HybridSurvival$FinalSubtype[i] = "ALS-Ox"
  }else if(HybridSurvival$FinalSubtype[i] == "TE"){
    HybridSurvival$FinalSubtype[i] = "ALS-TD"
  }
}

#With Discordant

#Add status
status = rep(NA,nrow(HybridSurvival))
for(i in 1:nrow(HybridSurvival)){
  if(is.na(HybridSurvival$AgeofOnset[i])){
    status[i] = 0
  }else if (is.na(HybridSurvival$Duration[i])){
    status[i] = 0
  }else if(is.na(HybridSurvival$AgeofDeath[i])){
    status[i] = 0
  }else{
    status[i]= 1
  }
}

FinalSurvival = cbind(HybridSurvival,status)
FinalSurvival$Duration = as.numeric(FinalSurvival$Duration)

km = with(FinalSurvival,Surv(Duration,status))
km_fit = survfit(Surv(Duration, status) ~ FinalSubtype, data=FinalSurvival)
surv_pvalue(km_fit)
print(km_fit)
summary(km_fit, times = c(1,30,60,90*(1:10)))

par(mar = c(5, 5, 3, 2))
plot(km_fit,col = c("goldenrod1","chartreuse2","darkorange2","navy","firebrick","gray50"),main="Hybrid Subtype Survival",ylim = c(-0.02,1.03),xlab="Disease Duration (months)",ylab="Survival Probability",lwd=2.5,cex.axis = 1.5,cex.lab=1.5,cex.main=1.5,xaxt="n")
timeseq = seq(0,156,12)
axis(side = 1, at = timeseq,labels = T,tick = T,cex.axis = 1.5)
legend(124,0.99,legend=c("ALS-Glia","ALS-Glia + ALS-Ox","ALS-Glia + ALS-TD","ALS-Ox","ALS-TD","Discordant"),lty = 1,lwd=3.5,col=c("goldenrod1","chartreuse2","darkorange2","navy","firebrick","gray50"))



#Without Discordant

remind = which(HybridSurvival$FinalSubtype == "Discordant")
HybridSurvivalnodisc = HybridSurvival[-remind,]

status = rep(NA,nrow(HybridSurvivalnodisc))
for(i in 1:nrow(HybridSurvivalnodisc)){
  if(is.na(HybridSurvivalnodisc$AgeofOnset[i])){
    status[i] = 0
  }else if (is.na(HybridSurvivalnodisc$Duration[i])){
    status[i] = 0
  }else if(is.na(HybridSurvivalnodisc$AgeofDeath[i])){
    status[i] = 0
  }else{
    status[i]= 1
  }
}

FinalSurvival = cbind(HybridSurvivalnodisc,status)
FinalSurvival$Duration = as.numeric(FinalSurvival$Duration)

km = with(FinalSurvival,Surv(Duration,status))
km_fit = survfit(Surv(Duration, status) ~ FinalSubtype, data=FinalSurvival)
surv_pvalue(km_fit)
print(km_fit)
summary(km_fit, times = c(1,30,60,90*(1:10)))

par(mar = c(5, 5, 3, 2))
plot(km_fit,col = c("goldenrod1","chartreuse2","darkorange2","navy","firebrick"),main="Hybrid Subtype Survival",xlab="Disease Duration (months)",ylab="Survival Probability",lwd=2.5,cex.axis = 1.5,cex.lab=1.5,cex.main=1.5,xaxt="n")
timeseq = seq(0,156,12)
axis(side = 1, at = timeseq,labels = T,tick = T,cex.axis = 1.5)
legend(124,0.99,legend=c("ALS-Glia","ALS-Glia + ALS-Ox","ALS-Glia + ALS-TD","ALS-Ox","ALS-TD"),lty = 1,lwd=3.5,col=c("goldenrod1","chartreuse2","darkorange2","navy","firebrick"))

################################################ Other Clinical Parameters

table(HybridSurvival$FinalSubtype)
TDind = which(HybridSurvival$FinalSubtype == "ALS-TD")
Oxind = which(HybridSurvival$FinalSubtype == "ALS-Ox")
Gliaind = which(HybridSurvival$FinalSubtype == "ALS-Glia")
Dind = which(HybridSurvival$FinalSubtype == "Discordant")
GTind = which(HybridSurvival$FinalSubtype == "ALS-Glia + ALS-TD Hybrid")
GOind = which(HybridSurvival$FinalSubtype == "ALS-Glia + ALS-Ox Hybrid")

######################################### Age of onset

TDaoo = as.numeric(HybridSurvival$AgeofOnset[TDind])
TDaoo = TDaoo[!is.na(TDaoo)]
OXaoo = as.numeric(HybridSurvival$AgeofOnset[Oxind])
OXaoo = OXaoo[!is.na(OXaoo)]
GLIAaoo = as.numeric(HybridSurvival$AgeofOnset[Gliaind])
GLIAaoo = GLIAaoo[!is.na(GLIAaoo)]
GTaoo = as.numeric(HybridSurvival$AgeofOnset[GTind])
GTaoo = GTaoo[!is.na(GTaoo)]
GOaoo = as.numeric(HybridSurvival$AgeofOnset[GOind])
GOaoo = GOaoo[!is.na(GOaoo)]
Discaoo = as.numeric(HybridSurvival$AgeofOnset[Dind])
Discaoo = Discaoo[!is.na(Discaoo)]

par(mar = c(5, 5, 3, 2))
boxplot(GLIAaoo,OXaoo,TDaoo,Discaoo,GTaoo,GOaoo,xaxt="n",ylim = c(25,85),main=c("Hybrid Subtype Age of Onset"),cex.axis = 1.5,cex.main=1.5,col=c("goldenrod1","navy","firebrick","gray50","darkorange2","chartreuse2"),xlab="Subtype",ylab="Age of Onset (years)",cex.lab=1.5,pch=20,cex=1.4)
axis(at=1:6,side=1,labels=c("ALS-Glia","ALS-Ox","ALS-TD","Discordant","Glia-TD Hybrid","Glia-Ox Hybrid"),cex.axis=1.2)

#Quick ANOVA to check for post-hoc tests
Age = HybridSurvival$AgeofOnset
ST = HybridSurvival$FinalSubtype
anovadat = data.frame(cbind(Age,ST))
anovadat$ST = factor(anovadat$ST,ordered = T)
levels(anovadat$ST)

anovadat = anovadat[-which(is.na(anovadat$Age)),]

oneway = aov(Age~ST,data=anovadat)
summary(oneway) #suggests no significant differences in age of onset by subtype

#Post hoc (p < 0.05)
t.test(GOaoo,OXaoo) #limited number of Glia-Ox hybrids reduces value of this analysis
t.test(GOaoo,GLIAaoo) #limited number of Glia-Ox hybrids reduces value of this analysis
t.test(GOaoo,TDaoo) #limited number of Glia-Ox hybrids reduces value of this analysis
t.test(GOaoo,GTaoo) #limited number of Glia-Ox hybrids reduces value of this analysis
t.test(GOaoo,Discaoo) #limited number of Glia-Ox hybrids reduces value of this analysis


########################################### Age of death 

TDaod = as.numeric(HybridSurvival$AgeofDeath[TDind])
TDaod = TDaod[!is.na(TDaod)]
OXaod = as.numeric(HybridSurvival$AgeofDeath[Oxind])
OXaod = OXaod[!is.na(OXaod)]
GLIAaod = as.numeric(HybridSurvival$AgeofDeath[Gliaind])
GLIAaod = GLIAaod[!is.na(GLIAaod)]
GTaod = as.numeric(HybridSurvival$AgeofDeath[GTind])
GTaod = GTaod[!is.na(GTaod)]
GOaod = as.numeric(HybridSurvival$AgeofDeath[GOind])
GOaod = GOaod[!is.na(GOaod)]
Discaod = as.numeric(HybridSurvival$AgeofDeath[Dind])
Discaod = Discaod[!is.na(Discaod)]

par(mar = c(5, 5, 3, 2))
boxplot(GLIAaod,OXaod,TDaod,Discaod,GTaod,GOaod,xaxt="n",ylim = c(30,90),main=c("Hybrid Subtype Age of Death"),cex.axis = 1.5,cex.main=1.5,col=c("goldenrod1","navy","firebrick","gray50","darkorange2","chartreuse2"),xlab="Subtype",ylab="Age of Death (years)",cex.lab=1.5,pch=20,cex=1.4)
axis(at=1:6,side=1,labels=c("ALS-Glia","ALS-Ox","ALS-TD","Discordant","Glia-TD Hybrid","Glia-Ox Hybrid"),cex.axis=1.2)

#Quick ANOVA to check for post-hoc tests
Age = HybridSurvival$AgeofDeath
ST = HybridSurvival$FinalSubtype
anovadat = data.frame(cbind(Age,ST))
anovadat$ST = factor(anovadat$ST,ordered = T)
levels(anovadat$ST)

anovadat = anovadat[-which(is.na(anovadat$Age)),]

oneway = aov(Age~ST,data=anovadat)
summary(oneway) #suggests no significant differences in age of onset by subtype

#Post hoc (confounded with age of onset; p < 0.05)
t.test(GOaod,GLIAaod) #limited number of Glia-Ox hybrids reduces value of this analysis
t.test(GOaod,OXaod) #limited number of Glia-Ox hybrids reduces value of this analysis
t.test(GOaod,TDaod) #limited number of Glia-Ox hybrids reduces value of this analysis
t.test(GOaod,Discaod) #limited number of Glia-Ox hybrids reduces value of this analysis
t.test(GOaod,GTaod) #limited number of Glia-Ox hybrids reduces value of this analysis


################################ Site of onset (categories: Bulbar, Limb, Other)

GLIAsoo = HybridSurvival$SiteofOnset[Gliaind]
OXsoo = HybridSurvival$SiteofOnset[Oxind]
TDsoo = HybridSurvival$SiteofOnset[TDind]
Discsoo = HybridSurvival$SiteofOnset[Dind]
GTsoo = HybridSurvival$SiteofOnset[GTind]
GOsoo = HybridSurvival$SiteofOnset[GOind]

table(HybridSurvival$SiteofOnset)

excelplot = data.frame(table(factor(GLIAsoo,levels = c("Axial","Axial and Bulbar","Axial and Limb","Bulbar","Bulbar and Limb","Generalized","Limb","Not Applicable","Unknown"))),table(factor(OXsoo,levels = c("Axial","Axial and Bulbar","Axial and Limb","Bulbar","Bulbar and Limb","Generalized","Limb","Not Applicable","Unknown"))),table(factor(TDsoo,levels = c("Axial","Axial and Bulbar","Axial and Limb","Bulbar","Bulbar and Limb","Generalized","Limb","Not Applicable","Unknown"))),table(factor(Discsoo,levels = c("Axial","Axial and Bulbar","Axial and Limb","Bulbar","Bulbar and Limb","Generalized","Limb","Not Applicable","Unknown"))),table(factor(GTsoo,levels = c("Axial","Axial and Bulbar","Axial and Limb","Bulbar","Bulbar and Limb","Generalized","Limb","Not Applicable","Unknown"))),table(factor(GOsoo,levels = c("Axial","Axial and Bulbar","Axial and Limb","Bulbar","Bulbar and Limb","Generalized","Limb","Not Applicable","Unknown"))))
colnames(excelplot) = c("ALS-Glia","Freq","ALS-Ox","Freq","ALS-TD","Freq","Discordant","Freq","Glia-TD","Freq","Glia-Ox","Freq")

#write.csv(excelplot,"Hybrid_SiteofOnset_Table_PeerReview.csv")

################################ FTLD Comorbidity

Glia.FTLD = HybridSurvival$Subcategory[Gliaind]
Ox.FTLD = HybridSurvival$Subcategory[Oxind]
TD.FTLD = HybridSurvival$Subcategory[TDind]
Disc.FTLD = HybridSurvival$Subcategory[Dind]
GT.FTLD = HybridSurvival$Subcategory[GTind]
GO.FTLD = HybridSurvival$Subcategory[GOind]


GCats = names(table(HybridSurvival$Subcategory[Gliaind])) #3 and 5 are FTD
blank = rep(NA,length(Gliaind))
for(j in 1:length(Glia.FTLD)){
  if(Glia.FTLD[j] == GCats[3]){
    blank[j] = "Positive"
  }else if(Glia.FTLD[j] == GCats[5]){
    blank[j] = "Positive"
  }else{
    blank[j] = "Negative"
  }
}

Glia.FTLD.binary = blank


OCats = names(table(HybridSurvival$Subcategory[Oxind])) #5 and 6 are FTD
blank = rep(NA,length(Oxind))
for(j in 1:length(Ox.FTLD)){
  if(Ox.FTLD[j] == OCats[5]){
    blank[j] = "Positive"
  }else if(Ox.FTLD[j] == OCats[6]){
    blank[j] = "Positive"
  }else{
    blank[j] = "Negative"
  }
}

Ox.FTLD.binary = blank


TCats = names(table(HybridSurvival$Subcategory[TDind])) #6 is FTD
blank = rep(NA,length(TDind))
for(j in 1:length(TD.FTLD)){
  if(TD.FTLD[j] == TCats[6]){
    blank[j] = "Positive"
  }else{
    blank[j] = "Negative"
  }
}

TD.FTLD.binary = blank


DCats = names(table(HybridSurvival$Subcategory[Dind])) #4 is FTD
blank = rep(NA,length(Dind))
for(j in 1:length(Disc.FTLD)){
  if(Disc.FTLD[j] == DCats[4]){
    blank[j] = "Positive"
  }else{
    blank[j] = "Negative"
  }
}

Disc.FTLD.binary = blank


GTCats = names(table(HybridSurvival$Subcategory[GTind])) #3 is FTD
blank = rep(NA,length(GTind))
for(j in 1:length(GT.FTLD)){
  if(GT.FTLD[j] == GTCats[3]){
    blank[j] = "Positive"
  }else{
    blank[j] = "Negative"
  }
}

GT.FTLD.binary = blank


GOCats = names(table(HybridSurvival$Subcategory[GOind])) #No FTLD patients
GO.FTLD.binary = rep("Negative",length(GOind))


#Hybrid FTLD Comorbidity

Glia.CoMorb = length(which(Glia.FTLD.binary == "Positive"))/length(Glia.FTLD.binary)
Ox.CoMorb = length(which(Ox.FTLD.binary == "Positive"))/length(Ox.FTLD.binary)
TD.CoMorb = length(which(TD.FTLD.binary == "Positive"))/length(TD.FTLD.binary)
Disc.CoMorb = length(which(Disc.FTLD.binary == "Positive"))/length(Disc.FTLD.binary)
GT.CoMorb = length(which(GT.FTLD.binary == "Positive"))/length(GT.FTLD.binary)
GO.CoMorb = length(which(GO.FTLD.binary == "Positive"))/length(GO.FTLD.binary)

CoMorbData = c(Glia.CoMorb,Ox.CoMorb,TD.CoMorb,Disc.CoMorb,GT.CoMorb,GO.CoMorb)
barplot(CoMorbData*100,main="Hybrid Subtype FTLD Comorbidity",names.arg = c("ALS-Glia","ALS-Ox","ALS-TD","Discordant","Glia-TD Hybrid","Glia-Ox Hybrid"),cex.axis = 1.4,cex.names = 1.3, cex.lab=1.4,col=c("goldenrod1","navy","firebrick","gray50","darkorange2","chartreuse2"),ylim=c(0,25),ylab = "FTLD Comorbidity (%)",cex.main = 1.5)

