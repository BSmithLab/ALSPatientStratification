## Short script to assign subtype labels to discriminatory features

#NovaSeq Rep - NMF breakdown
#NMF1 = Glia
#NMF2 = TE
#NMF3 = Ox

wd = "C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2/NovaSeq/11"
setwd(wd)
NovaFeatures = read.csv("totfeatures.csv")

for(i in 1:nrow(NovaFeatures)){
  if(NovaFeatures$Group[i] == 1){
    NovaFeatures$Group[i] = "Glia"
  }else if(NovaFeatures$Group[i] == 2){
    NovaFeatures$Group[i] = "TD"
  }else{
    NovaFeatures$Group[i] = "Ox"
  }
}

NovaFeatures = NovaFeatures[,-2]

#HiSeq Rep - NMF breakdown
#NMF1 = Glia
#NMF2 = Ox
#NMF3 = TE

wd = "C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/AllSubjects2/HiSeq/11"
setwd(wd)
HiFeatures = read.csv("totfeatures.csv")

for(i in 1:nrow(HiFeatures)){
  if(HiFeatures$Group[i] == 1){
    HiFeatures$Group[i] = "Glia"
  }else if(HiFeatures$Group[i] == 2){
    HiFeatures$Group[i] = "Ox"
  }else{
    HiFeatures$Group[i] = "TD"
  }
}

HiFeatures = HiFeatures[,-2]


wd = "C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Prognostic Biomarkers"
setwd(wd)

# Classifier = read.csv("Supervised_Classifier_CombinedPlatform_All_1000.csv")
# rownames(Classifier) = Classifier[,1]
# Classifier = Classifier[,-1]

#FinalFeature = data.frame(matrix(NA,nrow=nrow(Classifier),ncol=20))
#rownames(FinalFeature) = rownames(Classifier)
#colnames(FinalFeature) = c(paste("NovaSeq",seq(1,10,1),sep=""),paste("HiSeq",seq(1,10,1),sep=""))


CleanNovaFeatures = NovaFeatures[NovaFeatures$GeneCard %in% rownames(Classifier),]

#These genes are missing due to the combination of features between sequencing platform cohorts
#test = rownames(Classifier)[! rownames(Classifier) %in% NovaFeatures$GeneCard]

CleanHiFeatures = HiFeatures[HiFeatures$GeneCard %in% rownames(Classifier),]

DiscFeature = data.frame(matrix(NA,nrow=nrow(Classifier),ncol=3))
colnames(DiscFeature) = c("Transcript","NovaSeqSubtype","HiSeqSubtype")

DiscFeature$Transcript = rownames(Classifier)

for(i in 1:nrow(DiscFeature)){
  for(j in 1:nrow(CleanNovaFeatures)){
    
    if(DiscFeature$Transcript[i] == CleanNovaFeatures$GeneCard[j]){
      
      DiscFeature$NovaSeqSubtype[i] = CleanNovaFeatures$Group[j]
      
    }
    
  }
}


for(i in 1:nrow(DiscFeature)){
  for(j in 1:nrow(CleanHiFeatures)){
    
    if(DiscFeature$Transcript[i] == CleanHiFeatures$GeneCard[j]){
      
      DiscFeature$HiSeqSubtype[i] = CleanHiFeatures$Group[j]
      
    }
    
  }
}

#Completed Replicate: 1,2,3,4,5,6,7,8,9,10
#FinalFeature$NovaSeq10 = DiscFeature$NovaSeqSubtype
#FinalFeature$HiSeq10 = DiscFeature$HiSeqSubtype

NovaSeq11 = DiscFeature$NovaSeqSubtype
HiSeq11 = DiscFeature$HiSeqSubtype

FinalFeature2 = cbind(FinalFeature,HiSeq11,NovaSeq11)

write.csv(FinalFeature2,"FeatureList_MajoritySubtypeAssignment.csv")

###############################################################################################################################################################################################################################################################

MSA = read.csv("FeatureList_MajoritySubtypeAssignment.csv")
rownames(MSA)= MSA[,1]
MSA = MSA[,-1]

MSA_mat = data.frame(matrix(NA,nrow(MSA),ncol=3))

colnames(MSA_mat) = c("feature","NovaSeq","HiSeq")
MSA_mat$feature = rownames(MSA)

for(i in 1:nrow(MSA_mat)){
  
  tmp = c(unlist(MSA[i,])[[1]],unlist(MSA[i,])[[2]],unlist(MSA[i,])[[3]],unlist(MSA[i,])[[4]],unlist(MSA[i,])[[5]],unlist(MSA[i,])[[6]],unlist(MSA[i,])[[7]],unlist(MSA[i,])[[8]],unlist(MSA[i,])[[9]],unlist(MSA[i,])[[10]],unlist(MSA[i,])[[11]])
  tmp2 = c(unlist(MSA[i,])[[12]],unlist(MSA[i,])[[13]],unlist(MSA[i,])[[14]],unlist(MSA[i,])[[15]],unlist(MSA[i,])[[16]],unlist(MSA[i,])[[17]],unlist(MSA[i,])[[18]],unlist(MSA[i,])[[19]],unlist(MSA[i,])[[20]],unlist(MSA[i,])[[21]],unlist(MSA[i,])[[22]])
  novamaj = table(tmp)
  himaj = table(tmp2)
  
  if(length(novamaj) != 0){
    MSA_mat$NovaSeq[i] = names(which(novamaj == max(novamaj)))
  }
  
  if(length(himaj) != 0)
  MSA_mat$HiSeq[i] = names(which(himaj == max(himaj)))
  
}

write.csv(MSA_mat,"ClassifierFeatureSet1681_SubtypeAssignment.csv")

###############################################################################################################################################################################################################################################################

NovaSigGenes = read.csv("NOVASEQ_ALS_HC_SubtypeSigGenes_11-5-21.csv")
blank = rep(NA,nrow(NovaSigGenes))
NovaSigGenes = cbind(NovaSigGenes,blank)

HiSigGenes = read.csv("HISEQ_ALS_HC_SubtypeSigGenes_11-5-21.csv")
blank = rep(NA,nrow(HiSigGenes))
HiSigGenes = cbind(HiSigGenes,blank)

for(i in 1:nrow(NovaSigGenes)){
  for(j in 1:nrow(MSA_mat)){
    
    if(NovaSigGenes$X[i] == MSA_mat$feature[j]){
      NovaSigGenes$blank[i] = MSA_mat$NovaSeq[j]
    }
    
  }
}

for(i in 1:nrow(HiSigGenes)){
  for(j in 1:nrow(MSA_mat)){
    
    if(HiSigGenes$X[i] == MSA_mat$feature[j]){
      HiSigGenes$blank[i] = MSA_mat$HiSeq[j]
    }
    
  }
}

write.csv(NovaSigGenes,"NOVASEQ_ALS_HC_SubtypeSigGenes_11-5-21_v2.csv")
write.csv(HiSigGenes,"HISEQ_ALS_HC_SubtypeSigGenes_11-5-21_v2.csv")

###############################################################################################################################################################################################################################################################

PlatformSigGenes = read.csv("BothPlatform_ALL_SubtypeSigGenes_11-5-21.csv")
blank = rep(NA,nrow(PlatformSigGenes))
blank2 = rep(NA,nrow(PlatformSigGenes))
PlatformSigGenes = cbind(PlatformSigGenes,blank,blank2)

for(i in 1:nrow(PlatformSigGenes)){
  for(j in 1:nrow(MSA_mat)){
    
    if(PlatformSigGenes$X[i] == MSA_mat$feature[j]){
      PlatformSigGenes$blank[i] = MSA_mat$NovaSeq[j]
    }
    if(PlatformSigGenes$X[i] == MSA_mat$feature[j]){
      PlatformSigGenes$blank2[i] = MSA_mat$HiSeq[j]
    }
    
  }
}

write.csv(PlatformSigGenes,"BothPlatform_ALL_SubtypeSigGenes_11-5-21_v2.csv")
