#############################  SQuIRE Transposable Element RNA-seq Post-Processing Part 2  ####################################################

#Written by: Jarrett Eshima
#Date: June 4th, 2021
#For: Use by the Dr. Barbara Smith Lab at Arizona State University

## Description: This code is intended to help with the organization of SQuIRE post-processing count .RData files
# This code will provide a single output matrix of TPM counts as rows and RNA-seq samples as columns.
# Only the TEs found in all samples/subjects will be considered for downstream analysis.

#Important Note: This code is not fully automated. User is required to determine the appropriate number of "chunks" to load (dependent on system mem)

###############################################################################################################################################
############ USER PARAMETERS

#Set working directory where .RData files are stored
setwd("D:/Jarrett/Research/Summer 2021/SQuIRE Count ALS/FINAL") #These paths should be the same

#File Handle - the lettering scheme that precedes the integer and .RData file type
handle = "SQuIRE_Post_Chunk"

#Number of "chunks" / .RData files
nchunk = 3

#Chunk Length (number of subjects in each chunk - make sure to adjust for the last chunk)
clength = 50

#If your last chunk is not the same length as other chunks - update arguments below
allequal = F #change this to false if all chunks do not have the same file number
totalfiles = 451

#Name output Count Matrix (leave the .csv part)
savefile = "SQuIRE_TEs_GSE153960_9-28-21.csv"

#Select directory for output file
outputdir = "C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/TE Quant/FINAL"

###############################################################################################################################################
### Read in .RDATA files (from Post-Processing Script) and identify overlapping features

#Piecemeal Approach - .RData batches of 3
myrdata = paste(handle,seq(7,9,1),".RData",sep = "")

#Load in all .RData files
for(y in 1:nchunk){
  load(myrdata[y])
  if((y %% 1) == 0) cat(".RData Chunk Read-In:",y,"\n")
}

#I have 9 chunks, so I am breaking it up into 1-3 and 4-6 and 7-9
#Comment the section corresponding to the chunks not being run

# tmp1 = get(paste("avalanche",clength*1,sep=""))
# tmp2 = get(paste("avalanche",clength*2,sep=""))
# tmp3 = get(paste("avalanche",clength*3,sep=""))
# TE_KEEP = tmp1[tmp1$TE_ID %in% tmp2$TE_ID,]
# TE_KEEP = TE_KEEP[TE_KEEP$TE_ID %in% tmp3$TE_ID,]
# write.csv(TE_KEEP,"TE_KEEP_Chunks1to3.csv")

# tmp1 = get(paste("avalanche",clength*4,sep=""))
# tmp2 = get(paste("avalanche",clength*5,sep=""))
# tmp3 = get(paste("avalanche",clength*6,sep=""))
# TE_KEEP = tmp1[tmp1$TE_ID %in% tmp2$TE_ID,]
# TE_KEEP = TE_KEEP[TE_KEEP$TE_ID %in% tmp3$TE_ID,]
# write.csv(TE_KEEP,"TE_KEEP_Chunks4to6.csv")


tmp1 = get(paste("avalanche",clength*7,sep=""))
tmp2 = get(paste("avalanche",clength*8,sep=""))
tmp3 = get(paste("avalanche",clength*9,sep=""))
TE_KEEP = tmp1[tmp1$TE_ID %in% tmp2$TE_ID,]
TE_KEEP = TE_KEEP[TE_KEEP$TE_ID %in% tmp3$TE_ID,]
write.csv(TE_KEEP,"TE_KEEP_Chunks7to9.csv")


#You can use this part to remove subjects already completed
# samp_memhelp = paste("Subject",seq(1,clength*3,1),sep="")
# av_memhelp = paste("avalanche",seq(1,clength*3,5),sep="")
# rm(list = samp_memhelp)
# rm(list = av_memhelp)

###############################################################################################################################################
#Combine each "half" to get the final list of TEs found in all subjects (minumum 1 count, 99 quality score)
setwd("D:/Jarrett/Research/Summer 2021/SQuIRE Count ALS/FINAL")

half1 = read.csv("TE_KEEP_Chunks1to3.csv")
half2 = read.csv("TE_KEEP_Chunks4to6.csv")
half3 = read.csv("TE_KEEP_Chunks7to9.csv")

FinalTEs = half1[half1$TE_ID %in% half2$TE_ID,]
FinalTEs = FinalTEs[FinalTEs$TE_ID %in% half3$TE_ID,]

IDs = FinalTEs$TE_ID

##############################################################################################################################################
#Generate a count matrix for the transposable elements
setwd("D:/Jarrett/Research/Summer 2021/SQuIRE Count ALS") #These paths should be the same
wd = "D:/Jarrett/Research/Summer 2021/SQuIRE Count ALS" #These paths should be the same
handle = "SRR"

allfiles = list.files(path = wd,pattern = ".txt")

allfilelength = nchar(allfiles)

#Identify substring handle
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

SubjNames = allfiles[TEI]
SubjNames = substr(SubjNames,1,11)

TECounts = matrix(NA,nrow=length(IDs),ncol=length(SubjNames))
rownames(TECounts) = IDs
colnames(TECounts) = SubjNames


for(i in 1:length(TEI)){

  index = TEI[i]
  nam = substr(allfiles[index],1,11)
  assign(nam,read.table(allfiles[index],header = T))
  
  tmp = get(nam)
  tmp = tmp[tmp$TE_ID %in% IDs,]
  
  for(j in 1:nrow(TECounts)){
    for(k in 1:nrow(tmp)){
      
      if(rownames(TECounts)[j] == tmp$TE_ID[k]){
        if(colnames(TECounts)[i] == nam){ #not necessary but just to be sure
          TECounts[j,i] = tmp$tot_counts[k]
        }
        
        
      }
      
    }
  }
  
  if((i %% 1) == 0) cat("TE Counts Completed for Subject:",i,"of",length(TEI),"\n")
  rm(list=nam)
}



##############################################################################################################################################
#TE Count matrix combined with gene count matrix in the next script (ALSPatientStratification_DifferentialExpression_ALSPatients.R)

###############################################################################################################################################

##Clean up accession names to alternative ID (CGND-HRA)

#Convert SRR to HGND
setwd("D:/Jarrett/Research/Summer 2021/SQuIRE Count ALS/Post")
lut = read.table("MetaData_GSE153960.txt",sep = ",",fill = T,header = T) #This file can be obtained from NCBI GEO repository (Accession: PRJNA644618)
cleanlut = lut[lut$Run %in% SubjNames,]
write.csv(cleanlut,"GSE153960_CleanMetaData.csv")


CGND_IDs = cleanlut$sample_id_alt
SRR_IDs = cleanlut$Run

Convert = matrix(cbind(SRR_IDs,CGND_IDs),ncol=2)
TECountsF = TECounts

for(i in 1:nrow(Convert)){
  for(j in 1:ncol(TECounts)){
    if(Convert[i,1] == colnames(TECounts)[j]){
      colnames(TECountsF)[j] = Convert[i,2]
    }
  }
}

setwd("D:/Jarrett/Research/Summer 2021/SQuIRE Count ALS/FINAL")
#write.csv(TECountsF,"TECounts_HGND_9-29-21.csv")
#save.image(file="GSE153960_TECounts_9-29-21.RData")
load("GSE153960_TECounts_6-13-21.RData")
