##Weighted Gene Correlation Network Analysis (WGCNA)
#For application with the feature set used for GSEA

#Written By: Jarrett Eshima
#For: Dr. Barbara Smith Lab

#Reference: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html

library(WGCNA)
########################################################################################################################

wd = "C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/GSEA"
setwd(wd)

ClassMat = read.csv("ALS_HC_ClassifierMatrix.csv") #Generated in previous script
rownames(ClassMat) = ClassMat[,1]

hist(nchar(ClassMat[,1])) #Transcripts with very long names are TEs
abline(v=20,col="red") #character length threshold

ClassMat2 = ClassMat[,-1]

#Simple Hierarchical Clustering of TE features
charlen = table(nchar(rownames(ClassMat2))>20)

TE_only = data.frame(matrix(NA,nrow = charlen[[2]],ncol=ncol(ClassMat2)))
rn = which(nchar(rownames(ClassMat2))>20)

rownames(TE_only) = rownames(ClassMat2)[rn]
colnames(TE_only) = colnames(ClassMat2)

for(i in 1:nrow(ClassMat2)){
  for(j in 1:nrow(TE_only)){
    
    if(rownames(ClassMat2)[i] == rownames(TE_only)[j]){
      TE_only[j,] = ClassMat2[i,]
    }
    
  }
}

tmp = as.numeric(unlist(TE_only))
TE_only_num = matrix(tmp,nrow = nrow(TE_only),ncol = ncol(TE_only))
rownames(TE_only_num) = rownames(TE_only)
colnames(TE_only_num) = colnames(TE_only)
TE_only_num = data.frame(TE_only_num)

datExpr = TE_only_num
datExpr = t(datExpr)
dim(datExpr)
rownames(datExpr) = gsub("\\.","-",rownames(datExpr))
datExpr = data.frame(datExpr)

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Network")

#Read in phenotype/clinical data
datTraits = read.csv("GEO Collaborator Metadata.csv") #File provided by NYGC ALS Consortium

cleanTraits = datTraits[datTraits$ExternalSampleId %in% rownames(datExpr),]

#Remove non-numeric phenotype data
keep = c(9,10,13,20)
finalTraits = cleanTraits[,keep]
rownames(finalTraits) = cleanTraits$ExternalSampleId



gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK


if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

dim(datExpr)

finaldatExpr = datExpr

#Convert to traits numeric
tmp = as.numeric(unlist(finalTraits))
finalTraits2 = matrix(tmp,nrow=nrow(finalTraits),ncol = ncol(finalTraits))
rownames(finalTraits2) = rownames(finalTraits); colnames(finalTraits2) = colnames(finalTraits)


#Convert phenotype data to z-score for heatmap
finalTraits3 = finalTraits2
for(i in 1:ncol(finalTraits2)){
  
  mu = mean(finalTraits2[,i],na.rm = T)
  sdev = sd(finalTraits2[,i],na.rm = T)
  
  tmp = (finalTraits2[,i]-mu)/sdev
  finalTraits3[,i] = tmp
}


#Reorder
pheno = matrix(NA,nrow(finalTraits3),ncol(finalTraits3))
rownames(pheno) = rownames(finaldatExpr)
colnames(pheno) = colnames(finalTraits3)

for(i in 1:nrow(finalTraits3)){
  for(j in 1:nrow(pheno)){
    
    if(rownames(pheno)[j] == rownames(finalTraits3)[i]){
      
      pheno[j,] = finalTraits3[i,]
      
    }
    
  }
}

subtypepheno = pheno[,1:3] 
subtypepheno = data.frame(subtypepheno)
names(subtypepheno) = c("Age of Onset","Age of Death","Disease Duration")

# Re-cluster samples
sampleTree2 = hclust(dist(finaldatExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(subtypepheno, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,cex.dendroLabels = 0.4,#marAll = c(6,6,6,6),
                    groupLabels = names(subtypepheno), 
                    main = "Sample dendrogram and trait heatmap")

########################################################################################################################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 10, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(finaldatExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#Power = 14
softPower = 14
adjacency = adjacency(finaldatExpr, power = softPower)
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "TE clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)



# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 10;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "TE dendrogram and module colors")


# Calculate eigengenes
MEList = moduleEigengenes(finaldatExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes (TEs)",
     xlab = "", sub = "")



MEDissThres = 0.05
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(finaldatExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs


sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()



# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

########################################################################################################################
#Get EigenGenes for TE expression only

nEig = 11
Eig1 = names(finaldatExpr)[moduleColors=="green"]
Eig2 = names(finaldatExpr)[moduleColors=="blue"]
Eig3 = names(finaldatExpr)[moduleColors=="turquoise"]
Eig4 = names(finaldatExpr)[moduleColors=="black"]
Eig5 = names(finaldatExpr)[moduleColors=="brown"]
Eig6 = names(finaldatExpr)[moduleColors=="magenta"]
Eig7 = names(finaldatExpr)[moduleColors=="pink"]
Eig8 = names(finaldatExpr)[moduleColors=="yellow"]
Eig9 = names(finaldatExpr)[moduleColors=="red"]
Eig10 = names(finaldatExpr)[moduleColors=="purple"]
Eig11 = names(finaldatExpr)[moduleColors=="grey"]

numrow = max(c(length(Eig1),length(Eig2),length(Eig3),length(Eig4),length(Eig5),length(Eig6),length(Eig7),length(Eig8),length(Eig9),length(Eig10),length(Eig11)))

EigenTE = matrix(NA,nrow=numrow,ncol=nEig)
colnames(EigenTE) = c("MEgreen","MEblue","MEturquoise","MEblack","MEbrown","MEmagenta","MEpink","MEyellow","MEred","MEpurple","MEgrey")
EigenTE = data.frame(EigenTE)

for(i in 1:nEig){
  tmp = get(paste("Eig",i,sep=""))
  a = length(tmp)
  b = numrow - a
  
  tmp = c(tmp,rep(NA,b))
  assign(paste("Eig",i,sep=""),tmp)
}


EigenTE$MEgreen = Eig1
EigenTE$MEblue = Eig2
EigenTE$MEturquoise = Eig3
EigenTE$MEblack = Eig4
EigenTE$MEbrown = Eig5
EigenTE$MEmagenta = Eig6
EigenTE$MEpink = Eig7
EigenTE$MEyellow = Eig8
EigenTE$MEred = Eig9
EigenTE$MEpurple = Eig10
EigenTE$MEgrey = Eig11

write.csv(EigenTE,"WGCNA_EigenTE_Table_minsize10.csv")

