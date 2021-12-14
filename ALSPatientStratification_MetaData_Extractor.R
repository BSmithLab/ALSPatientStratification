#MetaData File Obtained from: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA644618&o=acc_s%3Aa

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping")
Meta = read.csv("GSE153960_MetaData.txt")

###ALS Patients Only
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Pub Data/Supplemental Files")
SRRs = read.csv("SRR_IDs.csv")
colnames(SRRs) = SRRs[1,]
SRRs = SRRs[-1,]
SRRs[1:5,]

keepsrr = SRRs$ALS

index1 = rep(NA,nrow(Meta))
count = 1
for(i in 1:nrow(Meta)){
  for(j in 1:length(keepsrr)){
    if(Meta$Run[i] == keepsrr[j])
      index1[count] = i
    count = count+1
  }
}

index1 = index1[!is.na(index1)]

EshimaMetaData = Meta[index1,] 

#Run these two lines to remove less important metadata information
#RC = c(2,3,4,7,8,9,10,11,12,19,23,26)
#TrimmedEshimaMD = EshimaMetaData[,-RC]


setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping")
write.csv(EshimaMetaData,"Eshima_CleanMetaData_ALS_Patients.csv")

#How many ALS PATIENTS (not individual samples)
length(table(EshimaMetaData$Subject_ID))



###Control Subjects
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Clustering w 451ALS/Controls")

HC = SRRs[,2]
OMND = SRRs[,3]

for(i in 1:length(HC)){
  if(HC[i] == ""){
    HC[i] = NA
  }
}

HC = HC[!is.na(HC)]

for(i in 1:length(OMND)){
  if(OMND[i] == ""){
    OMND[i] = NA
  }
}

OMND = OMND[!is.na(OMND)]

Controls = c(unlist(HC),unlist(OMND))


index2 = rep(NA,nrow(Meta))
count = 1
for(i in 1:nrow(Meta)){
  for(j in 1:length(Controls)){
    if(Meta$Run[i] == Controls[j])
      index2[count] = i
      count = count+1
  }
}

index2 = index2[!is.na(index2)]

EshimaMetaData_Controls = Meta[index2,] 


index3 = rep(NA,nrow(Meta))
count = 1
for(i in 1:nrow(Meta)){
  for(j in 1:length(HC)){
    if(Meta$Run[i] == HC[j])
      index3[count] = i
    count = count+1
  }
}

index3 = index3[!is.na(index3)]

EshimaMetaData_HC = Meta[index3,] 


index4 = rep(NA,nrow(Meta))
count = 1
for(i in 1:nrow(Meta)){
  for(j in 1:length(OMND)){
    if(Meta$Run[i] == OMND[j])
      index4[count] = i
    count = count+1
  }
}

index4 = index4[!is.na(index4)]

EshimaMetaData_OMND = Meta[index4,] 


setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping")
write.csv(EshimaMetaData_HC,"Eshima_CleanMetaData_HC.csv")
write.csv(EshimaMetaData_OMND,"Eshima_CleanMetaData_OMND.csv")
write.csv(EshimaMetaData_Controls,"Eshima_CleanMetaData_Control_Patients.csv")





# Generates a condensed metadata table (Subject metadata instead of Sample metadata) - some sample specific information is lost
tmpsrr = EshimaMetaData$Run
tmpMeta = Meta[Meta$Run %in% tmpsrr,]

rolling = NA
dup = names(which(table(tmpMeta$Subject_ID)>1))
for(j in 1:length(dup)){
  tmp = which(tmpMeta$Subject_ID == dup[j])
  
  if(length(tmp) != 0){
    n = length(tmp)-1
    s = sample(tmp,n,replace = F)
  }
  
  rolling = c(rolling,s)
}
rolling = rolling[!is.na(rolling)]
rollingMeta = tmpMeta[-rolling,]
rollingCGND = rollingMeta$sample_id_alt
rollingCGND = gsub("-",".",rollingCGND)
ALScoldata_Subjects = ALScoldata_451_TE[rownames(ALScoldata_451_TE) %in% rollingCGND,]
table(ALScoldata_Subjects$SubjectSex)