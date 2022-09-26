

# Differential Expression to assign subtype-specific features - for Enrichr-based enrichment

load("D:/Jarrett/Research/Nat Comm Reviews/RData/ALSPatientStratification_UnivariateDatasets_RINSite_PeerReview.RData")

GliaIDs = FullPheno_sr$Subject[which(FullPheno_sr$Subtype == "GLIA")]
OxIDs = FullPheno_sr$Subject[which(FullPheno_sr$Subtype == "OX")]
TDIDs = FullPheno_sr$Subject[which(FullPheno_sr$Subtype == "TE")]
HCIDs = FullPheno_sr$Subject[which(FullPheno_sr$Subtype == "Control")]
FTLDIDs = FullPheno_sr$Subject[which(FullPheno_sr$Subtype == "OND")]

nfeat = as.numeric(length(Gliaavg))

####################################################################################################################################

NormCounts2 = matrix(as.numeric(unlist(NormCounts)),nrow = nrow(NormCounts),ncol=ncol(NormCounts))
colnames(NormCounts2) = colnames(NormCounts)
rownames(NormCounts2) = rownames(NormCounts)


Glia = NormCounts2[,colnames(NormCounts2) %in% GliaIDs]
Ox = NormCounts2[,colnames(NormCounts2) %in% OxIDs]
TD = NormCounts2[,colnames(NormCounts2) %in% TDIDs]
HC = NormCounts2[,colnames(NormCounts2) %in% HCIDs]
FTLD = NormCounts2[,colnames(NormCounts2) %in% FTLDIDs]

Gliamed = Oxmed = TDmed = HCmed = FTLDmed = rep(NA,nfeat)

for(i in 1:nfeat){
  
  Gliamed[i] = median(Glia[i,],na.rm = T)
  Oxmed[i] = median(Ox[i,],na.rm = T)
  TDmed[i] = median(TD[i,],na.rm = T)
  HCmed[i] = median(HC[i,],na.rm = T)
  FTLDmed[i] = median(FTLD[i,],na.rm = T)
  
}

names(Gliamed) = rownames(NormCounts2)
names(Oxmed) = rownames(NormCounts2)
names(TDmed) = rownames(NormCounts2)
names(HCmed) = rownames(NormCounts2)
names(FTLDmed) = rownames(NormCounts2)

GliaFeats_Up_Median = OxFeats_Up_Median = TDFeats_Up_Median = HCFeats_Up_Median = FTLDFeats_Up_Median = GliaFeats_Down_Median = OxFeats_Down_Median = TDFeats_Down_Median = HCFeats_Down_Median = FTLDFeats_Down_Median = rep(NA,nfeat)

#Median-based approach
for(i in 1:nfeat){
  
  submed = c(Gliamed[i],Oxmed[i],TDmed[i],HCmed[i],FTLDmed[i])
  
  if(Gliamed[i] == max(submed)){
    
    GliaFeats_Up_Median[i] = names(Gliamed[i])
    
  }else if(Oxmed[i] == max(submed)){
    
    OxFeats_Up_Median[i] = names(Oxmed[i])
    
  }else if(TDmed[i] == max(submed)){
    
    TDFeats_Up_Median[i] = names(TDmed[i])
    
  }
  
  if(Gliamed[i] == min(submed)){
    
    GliaFeats_Down_Median[i] = names(Gliamed[i])
    
  }else if(Oxmed[i] == min(submed)){
    
    OxFeats_Down_Median[i] = names(Oxmed[i])
    
  }else if(TDmed[i] == min(submed)){
    
    TDFeats_Down_Median[i] = names(TDmed[i])
    
  }
  
}


GliaFeats_Up_Median = GliaFeats_Up[!is.na(GliaFeats_Up)]
GliaFeats_Down_Median = GliaFeats_Down[!is.na(GliaFeats_Down)]
OxFeats_Up_Median = OxFeats_Up[!is.na(OxFeats_Up)]
OxFeats_Down_Median = OxFeats_Down[!is.na(OxFeats_Down)]
TDFeats_Up_Median = TDFeats_Up[!is.na(TDFeats_Up)]
TDFeats_Down_Median = TDFeats_Down[!is.na(TDFeats_Down)]


hist(nchar(GliaFeats_Up_Median))
a = which(nchar(GliaFeats_Up_Median) >20)
b = which(substr(GliaFeats_Up_Median,1,4) == "ENSG")
c = c(a,b)
GliaFeats_Up_Median = GliaFeats_Up_Median[-c]

hist(nchar(GliaFeats_Down_Median))
a = which(nchar(GliaFeats_Down_Median) >20)
b = which(substr(GliaFeats_Down_Median,1,4) == "ENSG")
c = c(a,b)
GliaFeats_Down_Median = GliaFeats_Down_Median[-c]

hist(nchar(OxFeats_Up_Median))
a = which(nchar(OxFeats_Up_Median) >20)
b = which(substr(OxFeats_Up_Median,1,4) == "ENSG")
c = c(a,b)
OxFeats_Up_Median = OxFeats_Up_Median[-c]

hist(nchar(OxFeats_Down_Median))
a = which(nchar(OxFeats_Down_Median) >20)
b = which(substr(OxFeats_Down_Median,1,4) == "ENSG")
c = c(a,b)
OxFeats_Down_Median = OxFeats_Down_Median[-c]

hist(nchar(TDFeats_Up_Median))
a = which(nchar(TDFeats_Up_Median) >20)
b = which(substr(TDFeats_Up_Median,1,4) == "ENSG")
c = c(a,b)
TDFeats_Up_Median = TDFeats_Up_Median[-c]

hist(nchar(TDFeats_Down_Median))
a = which(nchar(TDFeats_Down_Median) >20)
b = which(substr(TDFeats_Down_Median,1,4) == "ENSG")
c = c(a,b)
TDFeats_Down_Median = TDFeats_Down_Median[-c]


setwd("D:/Jarrett/Research/Nat Comm Reviews/Enrichment/Enrichr")
write.csv(GliaFeats_Up_Median,"Glia_Upregulated_Median.csv")
write.csv(GliaFeats_Down_Median,"Glia_Downregulated_Median.csv")
write.csv(OxFeats_Up_Median,"Ox_Upregulated_Median.csv")
write.csv(OxFeats_Down_Median,"Ox_Downregulated_Median.csv")
write.csv(TDFeats_Up_Median,"TD_Upregulated_Median.csv")
write.csv(TDFeats_Down_Median,"TD_Downregulated_Median.csv")
