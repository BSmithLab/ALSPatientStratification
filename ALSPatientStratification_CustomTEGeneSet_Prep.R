#Post-SAKE GSEA Preparation

#Written By: Jarrett Eshima
#For: Dr. Barbara Smith Lab
#7/29/21

#This script is designed to facilitate gene set enrichment analysis (GSEA, Broad Institute).
#Custom Transposable Element Gene Set limited to TE subfamilies rather than individual TEs


##REPBASE ANNOTATION: https://www.girinst.org/downloads/
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/GSEA/CustomTE")
RepBase = read.table("humrep.txt",header = F,fill = T)
head(RepBase)

Nsubfamilies = table(table(which(RepBase[,1] == "ID")))

HumanTEs = rep(NA,Nsubfamilies[[1]])
count = 1
for(i in 1:nrow(RepBase)){
  if(RepBase[i,1] == "ID"){
    HumanTEs[count] = RepBase[i,2]
    count = count+1
  }
}

which(table(HumanTEs)>1) #No duplicates :)

write.csv(HumanTEs,"RepBase_HumanTransposableElementSubfamilies.csv")

##RepBase does not contain all of the subfamilies annotated by SQuIRE/RepeatMasker
#Potentially a different verision of the file was obtained... this is no problem, supplement missing subfamilies using SQuIRE names
#Then filter out duplicates

squireTEs.i = which(substr(rownames(ALSCountData_451_TE),1,3) == "chr")
squireTEs = ALSCountData_451_TE[squireTEs.i,]

squireTEs2 = rownames(squireTEs)

missingsubfam = rep(NA,length(squireTEs2))
for(i in 1:length(squireTEs2)){
  tmp = squireTEs2[i]
  tmp = sub("\\:.*","",tmp)
  tmp2 = strsplit(tmp,"_")
  missingsubfam[i] = tmp2[[1]][4]
}

AllTEs = c(HumanTEs,missingsubfam)

TEList.all= names(table(AllTEs))

write.csv(TEList.all,"RepBase_RepeatMasker_HumanTransposableElementSubfamilies_ALL.csv")
