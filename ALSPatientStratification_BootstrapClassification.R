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
