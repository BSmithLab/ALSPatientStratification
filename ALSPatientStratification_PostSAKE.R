#Post-SAKE Univariate Statistical Analysis

#Written By: Jarrett Eshima
#For: Dr. Barbara Smith Lab
#1-9-2021

###############################################################################################################################################################

##Read in SAKE Cluster Data

#HiSeq
#setwd('C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects/HiSeq/RobustSubtypeAssignment/11/11')
#NovaSeq
setwd('C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2/HiSeq/11/11')

clusterexpression = read.csv('Original_plus_NMF.csv') #File generated from SAKE NMF output (save results from each run)

###############################################################################################################################################################
##########################################  PARSE NMF Clusters  ###############################################################################################
###############################################################################################################################################################

#Specify the number of clusters used in SAKE (max = 7)
k=3

#Parse NMF Clusters
samplenames = colnames(clusterexpression)[-1]
Cluster1=Cluster2=Cluster3=Cluster4=Cluster5=Cluster6=Cluster7 = matrix(rep(0,nrow(clusterexpression)*ncol(clusterexpression)),nrow=nrow(clusterexpression))
count=count2=count3=count4=count5=count6=count7=1
ClusterList = rep(NA,length(samplenames))

if(k == 2){
  for(i in 1:length(samplenames)){
    if(substr(samplenames[i],1,4) == "NMF1"){
      Cluster1[,count] = clusterexpression[,i+1]
      count = count+1
      ClusterList[i] = 1
    }else if(substr(samplenames[i],1,4) == "NMF2"){
      Cluster2[,count2] = clusterexpression[,i+1]
      count2 = count2+1
      ClusterList[i] = 2
    }
  }
}else if(k == 3){
  for(i in 1:length(samplenames)){
    if(substr(samplenames[i],1,4) == "NMF1"){
      Cluster1[,count] = clusterexpression[,i+1]
      count = count+1
      ClusterList[i] = 1
    }else if(substr(samplenames[i],1,4) == "NMF2"){
      Cluster2[,count2] = clusterexpression[,i+1]
      count2 = count2+1
      ClusterList[i] = 2
    }else if(substr(samplenames[i],1,4) == "NMF3"){
      Cluster3[,count3] = clusterexpression[,i+1]
      count3 = count3+1
      ClusterList[i] = 3
    }
  }
}else if(k == 4){
  for(i in 1:length(samplenames)){
    if(substr(samplenames[i],1,4) == "NMF1"){
      Cluster1[,count] = clusterexpression[,i+1]
      count = count+1
      ClusterList[i] = 1
    }else if(substr(samplenames[i],1,4) == "NMF2"){
      Cluster2[,count2] = clusterexpression[,i+1]
      count2 = count2+1
      ClusterList[i] = 2
    }else if(substr(samplenames[i],1,4) == "NMF3"){
      Cluster3[,count3] = clusterexpression[,i+1]
      count3 = count3+1
      ClusterList[i] = 3
    }else if(substr(samplenames[i],1,4) == "NMF4"){
      Cluster4[,count4] = clusterexpression[,i+1]
      count4 = count4+1
      ClusterList[i] = 4
    }
  }
}else if(k == 5){
  for(i in 1:length(samplenames)){
    if(substr(samplenames[i],1,4) == "NMF1"){
      Cluster1[,count] = clusterexpression[,i+1]
      count = count+1
      ClusterList[i] = 1
    }else if(substr(samplenames[i],1,4) == "NMF2"){
      Cluster2[,count2] = clusterexpression[,i+1]
      count2 = count2+1
      ClusterList[i] = 2
    }else if(substr(samplenames[i],1,4) == "NMF3"){
      Cluster3[,count3] = clusterexpression[,i+1]
      count3 = count3+1
      ClusterList[i] = 3
    }else if(substr(samplenames[i],1,4) == "NMF4"){
      Cluster4[,count4] = clusterexpression[,i+1]
      count4 = count4+1
      ClusterList[i] = 4
    }else if(substr(samplenames[i],1,4) == "NMF5"){
      Cluster5[,count5] = clusterexpression[,i+1]
      count5 = count5+1
      ClusterList[i] = 5
    }
  }
}else if(k == 6){
  for(i in 1:length(samplenames)){
    if(substr(samplenames[i],1,4) == "NMF1"){
      Cluster1[,count] = clusterexpression[,i+1]
      count = count+1
      ClusterList[i] = 1
    }else if(substr(samplenames[i],1,4) == "NMF2"){
      Cluster2[,count2] = clusterexpression[,i+1]
      count2 = count2+1
      ClusterList[i] = 2
    }else if(substr(samplenames[i],1,4) == "NMF3"){
      Cluster3[,count3] = clusterexpression[,i+1]
      count3 = count3+1
      ClusterList[i] = 3
    }else if(substr(samplenames[i],1,4) == "NMF4"){
      Cluster4[,count4] = clusterexpression[,i+1]
      count4 = count4+1
      ClusterList[i] = 4
    }else if(substr(samplenames[i],1,4) == "NMF5"){
      Cluster5[,count5] = clusterexpression[,i+1]
      count5 = count5+1
      ClusterList[i] = 5
    }else if(substr(samplenames[i],1,4) == "NMF6"){
      Cluster6[,count6] = clusterexpression[,i+1]
      count6 = count6+1
      ClusterList[i] = 6
    }
  }
}else if(k == 7){
  for(i in 1:length(samplenames)){
    if(substr(samplenames[i],1,4) == "NMF1"){
      Cluster1[,count] = clusterexpression[,i+1]
      count = count+1
      ClusterList[i] = 1
    }else if(substr(samplenames[i],1,4) == "NMF2"){
      Cluster2[,count2] = clusterexpression[,i+1]
      count2 = count2+1
      ClusterList[i] = 2
    }else if(substr(samplenames[i],1,4) == "NMF3"){
      Cluster3[,count3] = clusterexpression[,i+1]
      count3 = count3+1
      ClusterList[i] = 3
    }else if(substr(samplenames[i],1,4) == "NMF4"){
      Cluster4[,count4] = clusterexpression[,i+1]
      count4 = count4+1
      ClusterList[i] = 4
    }else if(substr(samplenames[i],1,4) == "NMF5"){
      Cluster5[,count5] = clusterexpression[,i+1]
      count5 = count5+1
      ClusterList[i] = 5
    }else if(substr(samplenames[i],1,4) == "NMF6"){
      Cluster6[,count6] = clusterexpression[,i+1]
      count6 = count6+1
      ClusterList[i] = 6
    }else if(substr(samplenames[i],1,4) == "NMF7"){
      Cluster7[,count7] = clusterexpression[,i+1]
      count7 = count7+1
      ClusterList[i] = 7
    }
  }
}


#Clean up
index = which(colSums(Cluster1)==0)
Cluster1 = Cluster1[,-index]
index = which(colSums(Cluster2)==0)
Cluster2 = Cluster2[,-index]
index = which(colSums(Cluster3)==0)
Cluster3 = Cluster3[,-index]
index = which(colSums(Cluster4)==0)
Cluster4 = Cluster4[,-index]
index = which(colSums(Cluster5)==0)
Cluster5 = Cluster5[,-index]
index = which(colSums(Cluster6)==0)
Cluster6 = Cluster6[,-index]
index = which(colSums(Cluster7)==0)
Cluster7 = Cluster7[,-index]
dim(Cluster1);dim(Cluster2);dim(Cluster3);dim(Cluster4);dim(Cluster5);dim(Cluster6);dim(Cluster7)

ClusterList = matrix(ClusterList,nrow=length(ClusterList))
rownames(ClusterList) = samplenames
colnames(ClusterList) = "Subtype"



###############################################################################################################################################################
###################################  Majority Subtype Assignment  #############################################################################################
###############################################################################################################################################################

## Generate Containers (only run once)

numrep = 11
HiSubtypeMat =  matrix(NA,nrow=numrep,ncol=nrow(ClusterList))
rownames(HiSubtypeMat) = paste("Replicate",seq(1,numrep,1),sep="")
colnames(HiSubtypeMat) = samplenames

numrep = 11
NovaSubtypeMat = matrix(NA,nrow=numrep,ncol=nrow(ClusterList))
rownames(NovaSubtypeMat) = paste("Replicate",seq(1,numrep,1),sep="")
colnames(NovaSubtypeMat) = samplenames

#################################################################################

##HiSeq Samples
itr = 11
NMF1 = "GLIA"
NMF2 = "OX"
NMF3 = "TE"

MajoritySubtype = ClusterList

for(i in 1:nrow(MajoritySubtype)){
  if(MajoritySubtype[i,] == 1){
    MajoritySubtype[i,] = NMF1
  }else if(MajoritySubtype[i,] == 2){
    MajoritySubtype[i,] = NMF2
  }else if(MajoritySubtype[i,] == 3){
    MajoritySubtype[i,] = NMF3
  }
}

HiSubtypeMat[itr,] = MajoritySubtype
#Load next SAKE result file, repeat for number of replicates


setwd('C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2/HiSeq')
file = "HiSeq_RobustSubtypeAssignment_10-11-21_Tiebreaker.csv"
write.csv(HiSubtypeMat,file)

#After completing all 10 or 11
Major2 = rep(NA,ncol(HiSubtypeMat))
for(i in 1:ncol(HiSubtypeMat)){
  tmp = table(HiSubtypeMat[,i])
  locate = as.numeric(which(tmp>5))
  if(length(locate) == 1){
    Major2[i] = names(tmp[locate])
  }
}

#If all values are false, tiebreaker round is not required.
table(is.na(Major2))

FinalHiSubMat = rbind(HiSubtypeMat,Major2)

file = "HiSeq_RobustSubtypeAssignment_10-11-21_majority.csv"
write.csv(FinalHiSubMat,file)



##NovaSeq Samples
itr = 11 #change this value each iteration
NMF1 = "GLIA" #Use the SAKE Enrichment plot to determine subtype identity across SAKE runs
NMF2 = "TE"
NMF3 = "OX"

MajoritySubtype = ClusterList

for(i in 1:nrow(MajoritySubtype)){
  if(MajoritySubtype[i,] == 1){
    MajoritySubtype[i,] = NMF1
  }else if(MajoritySubtype[i,] == 2){
    MajoritySubtype[i,] = NMF2
  }else if(MajoritySubtype[i,] == 3){
    MajoritySubtype[i,] = NMF3
  }
}

NovaSubtypeMat[itr,] = MajoritySubtype
#Load next SAKE result file, repeat for number of replicates


setwd('C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2/NovaSeq')
file = "NovaSeq_RobustSubtypeAssignment_10-11-21_Tiebreaker.csv"
write.csv(NovaSubtypeMat,file)

Major = rep(NA,ncol(NovaSubtypeMat))
for(i in 1:ncol(NovaSubtypeMat)){
  tmp = table(NovaSubtypeMat[,i])
  locate = as.numeric(which(tmp>5))
  if(length(locate) == 1){
    Major[i] = names(tmp[locate])
  }
}

table(is.na(Major))

FinalNovaSubMat = rbind(NovaSubtypeMat,Major)

file = "NovaSeq_RobustSubtypeAssignment_10-11-21_majority.csv"
write.csv(FinalNovaSubMat,file)
