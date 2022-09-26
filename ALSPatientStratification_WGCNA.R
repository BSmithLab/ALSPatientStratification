##Weighted Gene Co-Expression Network Analysis (WGCNA)
#For application with the feature set used for GSEA

#Written By: Jarrett Eshima
#For: Dr. Barbara Smith Lab

#Key Reference: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
#Citation: Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. BMC bioinformatics, 9(1), 1-13.

library(WGCNA)
########################################################################################################################

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Network")

#Read in expression data
datExpr = read.csv("Table_S5.csv") #File generated in previous script
rownames(datExpr) = datExpr[,1]
datExpr = datExpr[,-1]
datExpr = t(datExpr)
dim(datExpr)
rownames(datExpr) = gsub("\\.","-",rownames(datExpr))
datExpr = data.frame(datExpr)

#Read in phenotype/clinical data
datTraits = read.csv("GEO Collaborator Metadata.csv") #File provided by NYGC ALS Consortium

cleanTraits = datTraits[datTraits$ExternalSampleId %in% rownames(datExpr),]

#Remove non-numeric phenotype data
keep = c(9,10,13,20)
finalTraits = cleanTraits[,keep]
rownames(finalTraits) = cleanTraits$ExternalSampleId


finaldatExpr = datExpr


#Convert to numeric
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

########################################################################################################################

gsg = goodSamplesGenes(finaldatExpr, verbose = 3);
gsg$allOK


if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(finaldatExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(finaldatExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  finaldatExpr = finaldatExpr[gsg$goodSamples, gsg$goodGenes]
}

dim(finaldatExpr)

sampleTree = hclust(dist(finaldatExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


# Re-cluster samples
sampleTree2 = hclust(dist(finaldatExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(subtypepheno, signed = T);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,cex.dendroLabels = 0.4,marAll = c(6,6,6,6),
                    groupLabels = names(subtypepheno), 
                    main = "Sample dendrogram and trait heatmap")

########################################################################################################################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 11, to=20, by=1))
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


#Power = 13
softPower = 13
adjacency = adjacency(finaldatExpr, power = softPower)
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)



# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 25;
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
                    main = "Gene dendrogram and module colors")


# Calculate eigengenes
MEList = moduleEigengenes(finaldatExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
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
# Save module colors and labels for use in subsequent parts
#save(MEs, moduleLabels, moduleColors, geneTree, file = "ALS-networkConstruction-stepByStep.RData")
#load("WGCNA_ALS_Subtype.RData")


# Define numbers of genes and samples
nGenes = ncol(finaldatExpr)
nSamples = nrow(finaldatExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(finaldatExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, subtypepheno, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Add the number of individual genes that comprise the module eigengene 
MEs.n = rep(NA,length(names(MEs0)))
for(i in 1:length(names(MEs0))){
  MEs.n[i] = paste(names(MEs0)[i]," (",table(dynamicColors)[[i]],")",sep="")
}

#Reorder the names
newMEs = rep(NA,length(MEs.n))
for(i in 1:length(MEs)){
  
  char = nchar(names(MEs))[i]
  MEcolor = names(MEs)[i]
  for(j in 1:length(MEs.n)){
  
    if(substr(MEs.n[j],1,char) == MEcolor){
      if(is.na(newMEs[i])){
        newMEs[i] = MEs.n[j]
      }
    }
  }
  
}


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), " (",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 11, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = c("Age of Onset","Age of Death","Disease Duration"),
               yLabels = names(MEs),
               ySymbols = newMEs,
               colorLabels = FALSE,
               colors = rev(greenWhiteRed(50)),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.lab = 1.1,
               cex.text = 1.3,
               zlim = c(-1,1))
               #main = paste("Module-trait relationships"))

#write.csv(MEs,"ALS_EigenGeneMatrix.csv")


#################################################################################################################
#Custom WGCNA heatmaps

moduleTraitCor
textMatrix
WGCNAmat = data.frame(moduleTraitCor)
WGCNAp =  signif(moduleTraitPvalue, 1)
WGCNAp2 = - log10(WGCNAp)

library(ggplot2)

customheat = data.frame(matrix(NA,nrow(WGCNAmat)*ncol(WGCNAmat),3))
colnames(customheat) = c("Module","Response","Value")
plotorder = c(rep("MEpink",3),rep("MEred",3),rep("MEtan",3),rep("MEturquoise",3),rep("MEbrown",3),rep("MEgreen",3),rep("MEpurple",3),rep("MEgrey",3),rep("MEmagenta",3),rep("MEyellow",3),rep("MEblue",3),rep("MEsalmon",3),rep("MEblack",3),rep("MEgreenyellow",3))
customheat$Module = plotorder
customheat$Response = rep(colnames(WGCNAmat),nrow(WGCNAmat))

#Fill plot matrix
for(i in 1:nrow(customheat)){
  
  for(j in 1:nrow(WGCNAmat)){
    
    for(k in 1:ncol(WGCNAmat)){
      
      if(customheat$Module[i] == rownames(WGCNAmat)[j] && customheat$Response[i] == colnames(WGCNAmat)[k]){
        
        customheat$Value[i] = WGCNAmat[j,k]
        
      }
      
    }
    
    
  }
  
}

#Clean up Response variable names
for(i in 1:nrow(customheat)){
  
  if(customheat$Response[i] == "Age.of.Onset"){
    customheat$Response[i] = "Age of Onset"
  }else if(customheat$Response[i] == "Age.of.Death"){
    customheat$Response[i] = "Age of Death"
  }else if(customheat$Response[i] == "Disease.Duration"){
    customheat$Response[i] = "Disease Duration"
  }
  
}

customheat$Module = as.character(customheat$Module)
customheat$Module = factor(customheat$Module,levels=unique(customheat$Module))

p = ggplot(customheat,aes(Module,Response)) + geom_tile(aes(fill=Value),color="black")
p = p+scale_fill_gradient2(low="blue",mid="white",high="red",limits=c(-0.2,0.2))
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "white"))
p


customheatp = data.frame(matrix(NA,nrow(WGCNAp2)*ncol(WGCNAp2),3))
colnames(customheatp) = c("Module","Response","Value")
customheatp$Module = plotorder
customheatp$Response = rep(colnames(WGCNAp2),nrow(WGCNAp2))

#Fill plot matrix
for(i in 1:nrow(customheatp)){
  
  for(j in 1:nrow(WGCNAp2)){
    
    for(k in 1:ncol(WGCNAp2)){
      
      if(customheatp$Module[i] == rownames(WGCNAp2)[j] && customheatp$Response[i] == colnames(WGCNAp2)[k]){
        
        customheatp$Value[i] = WGCNAp2[j,k]
        
      }
      
    }
    
    
  }
  
}

#Clean up Response variable names
for(i in 1:nrow(customheatp)){
  
  if(customheatp$Response[i] == "Age.of.Onset"){
    customheatp$Response[i] = "Age of Onset"
  }else if(customheatp$Response[i] == "Age.of.Death"){
    customheatp$Response[i] = "Age of Death"
  }else if(customheatp$Response[i] == "Disease.Duration"){
    customheatp$Response[i] = "Disease Duration"
  }
  
}

customheatp$Module = as.character(customheatp$Module)
customheatp$Module = factor(customheatp$Module,levels=unique(customheatp$Module))


p = ggplot(customheatp,aes(Module,Response)) + geom_tile(aes(fill=Value),color="black")
p = p+scale_fill_gradient2(low="blue",mid="white",high="darkgreen",limits=c(0,4))
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "white"))
p



#Enrichment Heatmap
setwd("C:/Users/jeshima/Desktop/ALS Manuscript/Nature/Communications Submission/Peer Review")
EigGO = read.csv("EigengeneEnrichmentMatrix.csv") #Bonferroni-adjusted p-values from Table S9
rownames(EigGO) = EigGO$X; EigGO = EigGO[,-1]; EigGO = -log10(EigGO)
plotorderGO = c(rep("MEpink",4),rep("MEred",4),rep("MEtan",4),rep("MEturquoise",4),rep("MEbrown",4),rep("MEgreen",4),rep("MEpurple",4),rep("MEgrey",4),rep("MEmagenta",4),rep("MEyellow",4),rep("MEblue",4),rep("MEsalmon",4),rep("MEblack",4),rep("MEgreenyellow",4))

customheatGO = data.frame(matrix(NA,nrow(EigGO)*ncol(EigGO),3))
colnames(customheatGO) = c("Module","Response","Value")
customheatGO$Module = plotorderGO
customheatGO$Response = rep(colnames(EigGO),nrow(EigGO))

#Fill plot matrix
for(i in 1:nrow(customheatGO)){
  
  for(j in 1:nrow(EigGO)){
    
    for(k in 1:ncol(EigGO)){
      
      if(customheatGO$Module[i] == rownames(EigGO)[j] && customheatGO$Response[i] == colnames(EigGO)[k]){
        
        customheatGO$Value[i] = EigGO[j,k]
        
      }
      
    }
    
    
  }
  
}

#Clean up Response variable names
for(i in 1:nrow(customheatGO)){
  
  if(customheatGO$Response[i] == "Extracellular.Matrix"){
    customheatGO$Response[i] = "Extracellular Matrix"
  }else if(customheatGO$Response[i] == "Synaptic.Signaling"){
    customheatGO$Response[i] = "Synaptic Signaling"
  }else if(customheatGO$Response[i] == "Immune.Response"){
    customheatGO$Response[i] = "Immune Response"
  }
  
}

customheatGO$Module = as.character(customheatGO$Module)
customheatGO$Module = factor(customheatGO$Module,levels=unique(customheatGO$Module))


p = ggplot(customheatGO,aes(Module,Response)) + geom_tile(aes(fill=Value),color="black")
p = p+scale_fill_gradient2(low="blue",mid="white",high="darkgreen",limits=c(0,16))
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
p


#Regression with MEs to get subtype specificity in expression
Pheno = read.csv("Table_S13.csv")
rownames(MEs) = gsub("-","\\.",rownames(MEs))
head(MEs)

table(Pheno$Subject == rownames(MEs))

SubtypeBeta = data.frame(matrix(NA,nrow = length(names(table(Pheno$Subtype))),ncol(MEs)))
colnames(SubtypeBeta) = colnames(MEs)
rownames(SubtypeBeta) = c("GLIA","OX","TD")

for(i in 1:ncol(MEs)){
  
  for(j in 1:nrow(SubtypeBeta)){
    
    smartlevel = rownames(SubtypeBeta)[j]
    dummyreg = rep(NA,length(Pheno$Subtype))
    
    for(k in 1:length(Pheno$Subtype)){
      
      if(Pheno$Subtype[k] == smartlevel){
        dummyreg[k] = 1
      }else{
        dummyreg[k] = 0
      }
        
    }
      
    model = lm(MEs[,i] ~ dummyreg)
    tmp = summary(model)$coefficients[,1][[2]] #Beta coefficient for non-zero level
    SubtypeBeta[j,i] = tmp
  }
  
}

rownames(SubtypeBeta) = c("ALS-Glia","ALS-Ox","ALS-TD")


plotordersub = c(rep("MEpink",3),rep("MEred",3),rep("MEtan",3),rep("MEturquoise",3),rep("MEbrown",3),rep("MEgreen",3),rep("MEpurple",3),rep("MEgrey",3),rep("MEmagenta",3),rep("MEyellow",3),rep("MEblue",3),rep("MEsalmon",3),rep("MEblack",3),rep("MEgreenyellow",3))
colorder = c("ALS-Glia","ALS-Ox","ALS-TD")

customheatsub = data.frame(matrix(NA,nrow(SubtypeBeta)*ncol(SubtypeBeta),3))
colnames(customheatsub) = c("Module","Response","Value")
customheatsub$Module = plotordersub
customheatsub$Response = rep(colorder,ncol(SubtypeBeta))

for(i in 1:nrow(SubtypeBeta)){
  
  for(j in 1:ncol(SubtypeBeta)){
    
    for(k in 1:nrow(customheatsub)){
      
      if(rownames(SubtypeBeta)[i] == customheatsub$Response[k] && colnames(SubtypeBeta)[j] == customheatsub$Module[k]){
        customheatsub$Value[k] = SubtypeBeta[i,j]
      }
      
    }
    
  }
  
  
}

customheatsub$Module = as.character(customheatsub$Module)
customheatsub$Module = factor(customheatsub$Module,levels=unique(customheatsub$Module))

p = ggplot(customheatsub,aes(Module,Response)) + geom_tile(aes(fill=Value),color="black")
p = p+scale_fill_gradient2(low="blue",mid="white",high="red",limits=c(-0.08,0.08))
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
p




#Regression p-values
Subtypep = data.frame(matrix(NA,nrow = length(names(table(Pheno$Subtype))),ncol(MEs)))
colnames(Subtypep) = colnames(MEs)
rownames(Subtypep) = c("GLIA","OX","TD")

for(i in 1:ncol(MEs)){
  
  for(j in 1:nrow(Subtypep)){
    
    smartlevel = rownames(Subtypep)[j]
    dummyreg = rep(NA,length(Pheno$Subtype))
    
    for(k in 1:length(Pheno$Subtype)){
      
      if(Pheno$Subtype[k] == smartlevel){
        dummyreg[k] = 1
      }else{
        dummyreg[k] = 0
      }
      
    }
    
    model = lm(MEs[,i] ~ dummyreg)
    tmp = summary(model)$coefficients[,4][[2]] #Beta coefficient for non-zero level
    Subtypep[j,i] = tmp
  }
  
}

rownames(Subtypep) = c("ALS-Glia","ALS-Ox","ALS-TD")


plotordersub = c(rep("MEpink",3),rep("MEred",3),rep("MEtan",3),rep("MEturquoise",3),rep("MEbrown",3),rep("MEgreen",3),rep("MEpurple",3),rep("MEgrey",3),rep("MEmagenta",3),rep("MEyellow",3),rep("MEblue",3),rep("MEsalmon",3),rep("MEblack",3),rep("MEgreenyellow",3))
colorder = c("ALS-Glia","ALS-Ox","ALS-TD")

customheatsubp = data.frame(matrix(NA,nrow(Subtypep)*ncol(Subtypep),3))
colnames(customheatsubp) = c("Module","Response","Value")
customheatsubp$Module = plotordersub
customheatsubp$Response = rep(colorder,ncol(Subtypep))

for(i in 1:nrow(Subtypep)){
  
  for(j in 1:ncol(Subtypep)){
    
    for(k in 1:nrow(customheatsubp)){
      
      if(rownames(Subtypep)[i] == customheatsubp$Response[k] && colnames(Subtypep)[j] == customheatsubp$Module[k]){
        customheatsubp$Value[k] = Subtypep[i,j]
      }
      
    }
    
  }
  
  
}

customheatsubp$Module = as.character(customheatsubp$Module)
customheatsubp$Module = factor(customheatsubp$Module,levels=unique(customheatsubp$Module))

customheatsubp$Value = as.numeric(customheatsubp$Value)
customheatsubp$Value =p.adjust(customheatsubp$Value,method = "bonferroni")
refp = customheatsubp
customheatsubp$Value = -log10(customheatsubp$Value)

p = ggplot(customheatsubp,aes(Module,Response)) + geom_tile(aes(fill=Value),color="black")
p = p+scale_fill_gradient2(low="blue",mid="white",high="darkgreen",limits=c(0,55))
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
p

customheatsub #beta coeff
refp #Bonferroni adjusted p value
refp[which(refp$Value < 0.05),] #Sig boxes
#################################################################################################################


# Define variable weight containing the weight column of datTrait
DD = as.data.frame(subtypepheno$`Disease Duration`);
names(DD) = "Disease Duration"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(finaldatExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(finaldatExpr, DD, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(DD), sep="")
names(GSPvalue) = paste("p.GS.", names(DD), sep="")


module = "magenta"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for",names(DD)),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module,pch=20,cex=1.5)



#EigenGene associated with difference in age of onset and disease duration
names(finaldatExpr)[moduleColors=="purple"]



# Create the starting data frame
probes = names(finaldatExpr)
geneInfo0 = data.frame(Probes = probes,
                       #geneSymbol = annot$gene_symbol[probes2annot],
                       #LocusLinkID = annot$LocusLinkID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, DD, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Disease.Duration));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "AgeDeath_EigengeneInfo_12-2-21.csv")

########################################################################################################################

#Optional Enrichment

allLLIDs = names(finaldatExpr)

library(org.Hs.eg.db)

symbols = allLLIDs 
hs = org.Hs.eg.db
LUT = select(hs, 
       keys = symbols,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

allLLIDs = LUT$ENTREZID #This gene list is in the same order as expression data

GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "human", nBestP = 10)

#extract info
tmp = GOenr$bestPTerms
tmp2 = tmp$CC
tmp3 = tmp2$enrichment
write.csv(tmp3,"GO_CC_Eigengene.csv")

########################################################################################################################

#Build the Network based on EigenGenes

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
#dissTOM = 1-TOMsimilarityFromExpr(finaldatExpr, power = 13);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^14
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")



#Set nSelect to 1681 to get TOMplot of all genes
nSelect = 1681
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^14;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")




# Recalculate module eigengenes
MEs = moduleEigengenes(finaldatExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
DD = as.data.frame(subtypepheno$`Disease Duration`);
names(DD) = "Disease Duration"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, DD))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(10,10);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,3.5,1,3), marHeatmap = c(4,5,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)


# Plot the dendrogram
sizeGrWindow(8,8);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(5,5,3,3),
                      plotDendrograms = FALSE, xLabelsAngle = 90)



########################################################################################################################
#VisANT Input File

# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(finaldatExpr, power = 14);
# Select module
module = "yellow"
# Select module probes
probes = names(finaldatExpr)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0)
                            #probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )

