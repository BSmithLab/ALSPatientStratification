#Short script to combine sample Phenotype Data and subtype labels into file for Loom generation

#Written By: Jarrett Eshima
#For: Dr. Barbara Smith Lab


#Read in ALS coldata file
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Supervised Classification/Control")
ALS_Pheno = read.csv("ALS451_coldata_SUBTYPES.csv")
All_Pheno = read.csv("ColData.csv")

#Add three blank columns (sex, subtype, factor) to match ncol of ALS phenotype file
blank = rep(NA,nrow(All_Pheno))
All_Pheno = cbind(All_Pheno,blank,blank,blank)
colnames(All_Pheno) = colnames(ALS_Pheno)


#Control Sex can be determined from XIST and UTY expression as performed in previous scripts
#For the sake of time, subject sex will be read in from the NYGC ALS Consortium clinical data file
FullClinical = read.csv("GEO_Collaborator_PatientData.csv")
FullClinical$ExternalSampleId = gsub("-",".",FullClinical$ExternalSampleId)


#Read in the GSEA matrix
FeatureMat = read.csv("AllSubjects_ClassifierMatrix.csv")

Controlids = colnames(FeatureMat)[! colnames(FeatureMat) %in% ALS_Pheno$Subject]
Controlids = Controlids[-1]

#Combining the phenotype information in this way allows you to keep the sample order the same as the classifier dataframe
Control_Pheno = data.frame(matrix(NA,nrow = length(Controlids),ncol = ncol(ALS_Pheno)))
colnames(Control_Pheno) = colnames(ALS_Pheno)
Control_Pheno$Subject = Controlids

for(i in 1:nrow(Control_Pheno)){
  for(j in 1:nrow(All_Pheno)){
    
    if(Control_Pheno$Subject[i] == All_Pheno$Subject[j]){
      
      Control_Pheno[i,] = All_Pheno[j,]
      
    }
    
  }
}

Control_Pheno$Subtype = "Control"
Control_Pheno$Factor = 4


for(i in 1:nrow(Control_Pheno)){
  for(j in 1:nrow(FullClinical)){
    
    if(Control_Pheno$Subject[i] == FullClinical$ExternalSampleId[j]){
      
      Control_Pheno$SubjectSex[i] = FullClinical$Sex[j]
      
    }
    
  }
}

ALS_Control_Pheno = rbind(ALS_Pheno,Control_Pheno)


setwd('C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Supervised Classification/Control')
write.csv(ALS_Control_Pheno,"AllSubjects_coldata_SUBTYPES.csv")

#Note the cohorts including healthy control donors contain one "missing" value for sequencing platform
#Sample CGND.HRA.01976 was analyzed on a NovaSeq instrument 
#Manually corrected this value in excel or similar editor
