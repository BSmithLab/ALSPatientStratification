#Read/Generate Loom file for compatibility with Plaisier Lab [1] Supervised Classification Python Scripts

#Written By: Jarrett Eshima
#For: Dr. Barbara Smith Lab

#[1] https://github.com/plaisier-lab/U5_hNSC_Neural_G0
###########################################################################################################################
##This section only needs to be run once

#Install LoomR
#library(devtools)

#devtools::install_github(repo = "hhoeflin/hdf5r")
#devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")

###########################################################################################################################

library(loomR)


#Generate Loom File
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Supervised Classification/Control")
ExpData = read.csv("AllSubjects_ClassifierMatrix.csv") #File generated in previous script
rownames(ExpData) = ExpData[,1]
ExpData = ExpData[,-1]

#setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Supervised Classification")
loomR::create("AllSubjects_CA1000.loom",ExpData)

#Read in Pheno Data and other relevant files
Pheno = read.csv("AllSubjects_coldata_SUBTYPES.csv") #File generated in previous script
colnames(Pheno)[1] = "ID"

#The next line must only be true
table(Pheno$ID == colnames(ExpData))

#This line is run each time you want to interact with the loom object (after closing the object previously)
lfile = connect(filename = "AllSubjects_CA1000.loom",mode = "r+")

#Convert Pheno Data to Loom file
nsub=nrow(Pheno)
ngene = as.numeric(nrow(ExpData))
lfile$add.row.attribute(list(Selected = rep(1,ngene)),overwrite = T)
lfile[["row_attrs"]]

#lfile$add.col.attribute(list(IDs = Pheno$ID),overwrite = T) #This is the subject metadata/phenotype info
lfile$add.col.attribute(list(Sex = Pheno$SubjectSex),overwrite = T)
lfile$add.col.attribute(list(Tissue = Pheno$tissue),overwrite = T)
lfile$add.col.attribute(list(clust_ID = Pheno$Subtype),overwrite = T)
lfile$add.col.attribute(list(Disease = Pheno$disease_group),overwrite = T)
lfile$add.col.attribute(list(Platform = Pheno$sequencing_platform),overwrite = T)
lfile$add.col.attribute(list(Factor = Pheno$Factor),overwrite = T)
lfile$add.col.attribute(list(SubID = Pheno$ID),overwrite = T)

#Check that everything was added to the loom file
lfile[["col_attrs"]]


#Use this line to properly close out the loom file and avoid corrupting base files
lfile$close_all()
