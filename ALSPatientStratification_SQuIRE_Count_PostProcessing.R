#######################################  SQuIRE Transposable Element RNA-seq Post-Processing  #################################################

#Written by: Jarrett Eshima
#Date: June 1st, 2021
#For: Dr. Barbara Smith Lab at Arizona State University

## Description: This code is intended to help with the post-processing of SQuIRE TE Count Output .txt files.
# This post-processing code can be run automatically (source) after changing the "User Parameters" Section below.
# If you encounter errors running the code from source, run the lines individually
# Note: This computation is extremely resource intensive, especially RAM, since each file is ~0.25Gb. Not feasible to run batches of >100 or so subjects.

# You should clear the global environment every time this code runs successfully... Otherwise .RData files will accumulate and get too large
# The SQuIRE output .txt files should be placed in a clean folder with no other files (folders are OK)
#######################################################################################################################################
############ USER PARAMETERS

#Add repositories for dependent packages
#setRepositories() #Access BioConductor and CRAN extras (R-Forge and Omegahat also available) Paste: 1 2 3 4 5

run = TRUE

#Set working directory where count files are stored
setwd("D:/Jarrett/Research/Summer 2021/SQuIRE Count ALS") #These paths should be the same
wd = "D:/Jarrett/Research/Summer 2021/SQuIRE Count ALS" #These paths should be the same

#File Handle - the lettering scheme that precedes the accession number (unique sample identifier for GEO)
handle = "SRR"

#Minimum score threshold to keep TE
qualityscore = 99 #1% false positives ("TEs with few uniquely aligned reads may be prone to misrepresentation" - Refer to SQuIRE paper)

#Name output RData file - Note: it is important that your name is trailed by a single integer which indicates which run it came from.
#Naming your .RData files in this way will allow you to run SQuIRE_Count_Meta.R without errors.
file = "SQuIRE_Post_Chunk9.RData"

#Select directory for output
outputdir = "D:/Jarrett/Research/Summer 2021/SQuIRE Count ALS/FINAL"

#Manually specify the files you want to run...
chunk_start = 401
chunk_stop = 451

#Does each chunk have the same number of samples? If not, the last chunk will be set to have an alternative number of samples. 
lastchunkequal = F


#######################################################################################################################################
############ SUBSET THE DATA

#Get file names in the directory
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


#######################################################################################################################################
############ READ IN THE DATA - Lower RAM requirement but slower
if(run == TRUE){

  if(chunk_stop-chunk_start > 75 && chunk_stop-chunk_start < 150){
    cat("There may be too many files in this batch...")
  }

  if(chunk_stop-chunk_start > 150 && chunk_stop-chunk_start < 400){
    cat("Warning: Recommend lowering batch size")
  }

  if(chunk_stop-chunk_start > 400){
    cat("Error: Lower batch size")
    break
  }


#Read in Subjects for Chunk
for(i in chunk_start:chunk_stop){
  if(i == chunk_start){
    pseudoc = 1
  }
  nam = paste("Subject",i,sep = "")
  index = TEI[i]
  assign(nam,read.table(allfiles[index],header = T))
  if((pseudoc %% 5) == 0) cat("TE Read-In %:",pseudoc/(chunk_stop-chunk_start+1)*100,"\n")
  pseudoc = pseudoc+1
}

#Filter by Quality Score
for(i in chunk_start:chunk_stop){
  tmp = get(paste("Subject",i,sep=""))
  KEEP = rep(NA,nrow(tmp))
  count = 1
  for(k in 1:nrow(tmp)){
    if(tmp$score[k] >= qualityscore){
      KEEP[count] = k
      count = count+1
      if((k %% 100000) == 0) cat("TE Score Filtering In Progress, % Done:",round(k/nrow(tmp)*100,0),"\n")
    }
  }
  KEEP = KEEP[!is.na(KEEP)]
  assign(paste("Subject",i,sep=""),tmp[KEEP,])
  if((i %% 1) == 0) cat("TEs Filtered by Score for Subject:",i,"\n")
}

#TPM Normalization
for(i in chunk_start:chunk_stop){
    tmp = get(paste("Subject",i,sep=""))
    nr = nrow(tmp)
    nc = ncol(tmp)
  #Perform TPM Normalization
    blank = matrix(NA,nrow=nr,ncol = 1)
    tmp = cbind(tmp,blank)
    colnames(tmp)[nc+1] = "TPM"

    #FPKM direct to TPM
    scalefactor = 1e6/sum(tmp$fpkm)
    tmp$TPM = tmp$fpkm*scalefactor

    assign(paste("Subject",i,sep=""),tmp)
    if((i %% 1) == 0) cat("TE Counts Normalized for Subject:",i,"\n")
}

#Restrict TEs to those present in all subjects in a chunk
#A secondary script: SQuIRE_Count_Meta.R is required to compare names across Chunks

#Overlapping Name Avalanche

for(i in chunk_start:chunk_stop){
  if(i < chunk_stop){
    if(i == chunk_start){
      temp1 = get(paste("Subject",i,sep=""))
      temp2 = get(paste("Subject",i+1,sep=""))
      snow = temp1[temp1$TE_ID %in% temp2$TE_ID,]
    }else{
      temp1 = snow
      temp2 = get(paste("Subject",i+1,sep=""))
      snow = temp1[temp1$TE_ID %in% temp2$TE_ID,]
      #assign(paste("avalanche",count,sep=""),snowball)
      #count = count+1
      if((i %% 1) == 0) cat("TE Avalanche for Subjects:",chunk_start,"-",i+1,"\n")
    }
  }
  if((i %% 5) == 0) assign(paste("avalanche",i,sep=""),snow)
}

if(lastchunkequal == F && i == length(TEI)){
  assign(paste("avalanche",i,sep=""),snow)
}

rm(tmp,temp2,temp1)
setwd(outputdir)
save.image(file=file)
cat("TE Post-Processing Completed Successfully for subjects:",chunk_start,"-",chunk_stop)
}



