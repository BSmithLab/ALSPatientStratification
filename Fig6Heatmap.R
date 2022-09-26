#Updated figure 6 to enhance visibility and interpretability


#Figure 6A

setwd("D:/Jarrett/Research/Nat Comm Reviews/Figure Edits/Fig6 - Univariate/Heatmap")

Expr = read.csv("Fig6Heatmap.csv") #Source data: Figure 6A median-of-ratios
rownames(Expr) = Expr$Gene
Expr = Expr[,-1]

Pheno = read.csv("AllSubjects_Phenotype.csv") #Table S13 sheet 2

tmp = which(Pheno$Subject == "CGND.HRA.01732")
Pheno = Pheno[-tmp,]

GID = which(Pheno$Subtype == "GLIA")
OID = which(Pheno$Subtype == "OX")
TID = which(Pheno$Subtype == "TE")
HID = which(Pheno$Subtype == "Control")
FID = which(Pheno$Subtype == "FTLD")

Expr2 = cbind(Expr[,GID],Expr[,OID])
Expr2 = cbind(Expr2,Expr[,TID])
ALSExpr = Expr2
ALSExpr = matrix(unlist(ALSExpr),nrow=nrow(ALSExpr),ncol=ncol(ALSExpr))
ALSExpr = log2(ALSExpr)
Expr2 = cbind(Expr2,Expr[,HID])
Expr2 = cbind(Expr2,Expr[,FID])

ZExpr = log2(Expr2)
ZExpr = matrix(unlist(ZExpr),nrow=nrow(ZExpr),ncol=ncol(ZExpr))
colnames(ZExpr) = colnames(Expr); rownames(ZExpr) = rownames(Expr)


ALSmu = ALSsd = rep(NA,nrow(ALSExpr))
for(i in 1:nrow(ALSExpr)){
  ALSmu[i] = mean(ALSExpr[i,],na.rm=T)
  ALSsd[i] = sd(ALSExpr[i,],na.rm=T)
}


for(i in 1:nrow(ZExpr)){
  for(j in 1:ncol(ZExpr)){
   
    ZExpr[i,j] = (ZExpr[i,j]-ALSmu[i])/ALSsd[i]
     
  }
  
}

resplab = rep(NA,ncol(ZExpr)*nrow(ZExpr))
for(i in 1:nrow(ZExpr)){
  
  count2 = i*ncol(ZExpr)
  if(i == 1){
    resplab[1:count2] = rep(rownames(ZExpr)[i],ncol(ZExpr))
    count1 = count2+1
  }else{
    resplab[count1:count2] = rep(rownames(ZExpr)[i],ncol(ZExpr))
    count1 = count2+1
  }
  
}


heatdat = data.frame(matrix(NA,nrow=ncol(ZExpr)*nrow(ZExpr),ncol=3))
colnames(heatdat) = c("Module","Response","Value")
heatdat$Module = rep(colnames(ZExpr),nrow(ZExpr))
heatdat$Response = resplab

for(i in 1:nrow(heatdat)){
  
  for(j in 1:nrow(ZExpr)){
    
    for(k in 1:ncol(ZExpr)){
      
      if(heatdat$Module[i] == colnames(ZExpr)[k] && heatdat$Response[i] == rownames(ZExpr)[j]){
        
        heatdat$Value[i] = ZExpr[j,k]
        
      }
      
    }
    
  }
  
  if((i %% 100) == 0) cat("% Done:",i/nrow(heatdat)*100,"\n")
  
}

#write.csv(heatdat,"fig6heatmap_rawdata_log2_ALSmu.csv")
#heatdat = read.csv("fig6heatmap_rawdata_log2_ALSmu.csv")

heatdat$Response = as.character(heatdat$Response)
heatdat$Response = factor(heatdat$Response,levels=c("NANOGP4","MIRLET7BHG","MIR24-2","LINC01347","HSP90AB4P","GATA2-AS1","ENSG00000273151","ENSG00000258674","ENSG00000205041","COL3A1","CHKB-CPT1B","AGPAT4-IT1","UCP2","UBQLN2","TCIRG1","SLC17A6","SLC6A13","SERPINI1","OXR1","HTR2A","GLRA3","GAD2","GABRA1","COL18A1","TYROBP","TREM2","TNC","TMEM125","TLR7","HLA-DRA","FOLH1","CX3CR1","CHI3L2","CD44","APOC2","AIF1"))

heatdat$Module = as.character(heatdat$Module)
heatdat$Module = factor(heatdat$Module,levels=make.unique(heatdat$Module))


library(ggplot2)
p = ggplot(heatdat,aes(Module,Response)) + geom_tile(aes(fill=Value),color="black")
p = p+scale_fill_gradient2(low="blue",mid="white",high="red",limits=c(-4,4))
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
p


#Clean up the extremes (for plotting purposes, not used in calculation of DE FDR q)

cleandat = heatdat

for(i in 1:nrow(cleandat)){
  
  if(cleandat$Value[i] < -4){
    cleandat$Value[i] = -4
  }else if(cleandat$Value[i] > 4){
    cleandat$Value[i] = 4
  }
  
}

p = ggplot(cleandat,aes(Module,Response)) + geom_tile(aes(fill=Value),color="black")
p = p+scale_fill_gradient2(low="blue",mid="white",high="red",limits=c(-4,4))
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
p

#Add Subtype lines
table(Pheno$Subtype)
myline = data.frame(x = 84,y=0,xend=84,yend=36)
p = p+geom_segment(data = myline,aes(x,y,xend=xend,yend=yend),size=0.5)
myline = data.frame(x = 84+239,y=0,xend=84+239,yend=36)
p = p+geom_segment(data = myline,aes(x,y,xend=xend,yend=yend),size=0.5)
myline = data.frame(x = 84+239+127,y=0,xend=84+239+127,yend=36)
p = p+geom_segment(data = myline,aes(x,y,xend=xend,yend=yend),size=0.5)
myline = data.frame(x = 84+239+127+93,y=0,xend=84+239+127+93,yend=36)
p = p+geom_segment(data = myline,aes(x,y,xend=xend,yend=yend),size=0.5)
myline = data.frame(x = 84+239+127+93+42,y=0,xend=84+239+127+93+42,yend=36)
p = p+geom_segment(data = myline,aes(x,y,xend=xend,yend=yend),size=0.5)
p

####################################################################################################################################################################################################

#Figure 6B
#plot p-values as heatmap

load("D:/Jarrett/Research/Nat Comm Reviews/RData/ALSPatientStratification_UnivariateDatasets_RINSite_PeerReview.RData")

ncomp = 10
nfeat = 36

pdat = data.frame(matrix(NA,nrow = ncomp*nfeat,ncol=3))
colnames(pdat) = c("Module","Response","Value")

plab = rep(NA,ncomp*nfeat)
for(i in 1:nrow(ZExpr)){
  
  count2 = i*ncomp
  if(i == 1){
    plab[1:count2] = rep(rownames(ZExpr)[i],ncomp)
    count1 = count2+1
  }else{
    plab[count1:count2] = rep(rownames(ZExpr)[i],ncomp)
    count1 = count2+1
  }
  
}

pdat$Response = plab
mods = c("Glia_vs_Ox","Glia_vs_TD","Ox_vs_TD","Glia_vs_HC","Ox_vs_HC","TD_vs_HC","Glia_vs_FTLD","Ox_vs_FTLD","TD_vs_FTLD","HC_vs_FTLD")
pdat$Module = rep(mods,nfeat)

pref = data.frame(matrix(NA,nrow=length(Transcripts),ncol = ncomp))
rownames(pref) = Transcripts

#Glia vs Ox
for(i in 1:nrow(pref)){
  for(j in 1:length(rownames(filt.GO.sig))){
    if(rownames(pref)[i] == rownames(filt.GO.sig)[j]){
      pref$X1[i] = as.numeric(filt.GO.sig$padj[j])
    }
  }
}
#Glia vs TD
for(i in 1:nrow(pref)){
  for(j in 1:length(rownames(filt.GT.sig))){
    if(rownames(pref)[i] == rownames(filt.GT.sig)[j]){
      pref$X2[i] = as.numeric(filt.GT.sig$padj[j])
    }
  }
}
#Ox vs TD
for(i in 1:nrow(pref)){
  for(j in 1:length(rownames(filt.TO.sig))){
    if(rownames(pref)[i] == rownames(filt.TO.sig)[j]){
      pref$X3[i] = as.numeric(filt.TO.sig$padj[j])
    }
  }
}
#Glia vs HC
for(i in 1:nrow(pref)){
  for(j in 1:length(rownames(filt.glia.sig))){
    if(rownames(pref)[i] == rownames(filt.glia.sig)[j]){
      pref$X4[i] = as.numeric(filt.glia.sig$padj[j])
    }
  }
}
#Ox vs HC
for(i in 1:nrow(pref)){
  for(j in 1:length(rownames(filt.ox.sig))){
    if(rownames(pref)[i] == rownames(filt.ox.sig)[j]){
      pref$X5[i] = as.numeric(filt.ox.sig$padj[j])
    }
  }
}
#TD vs HC
for(i in 1:nrow(pref)){
  for(j in 1:length(rownames(filt.TE.sig))){
    if(rownames(pref)[i] == rownames(filt.TE.sig)[j]){
      pref$X6[i] = as.numeric(filt.TE.sig$padj[j])
    }
  }
}
#Glia vs FTLD
for(i in 1:nrow(pref)){
  for(j in 1:length(rownames(filt.glia.sig.ond))){
    if(rownames(pref)[i] == rownames(filt.glia.sig.ond)[j]){
      pref$X7[i] = as.numeric(filt.glia.sig.ond$padj[j])
    }
  }
}
#Ox vs FTLD
for(i in 1:nrow(pref)){
  for(j in 1:length(rownames(filt.ox.sig.ond))){
    if(rownames(pref)[i] == rownames(filt.ox.sig.ond)[j]){
      pref$X8[i] = as.numeric(filt.ox.sig.ond$padj[j])
    }
  }
}
#TD vs FTLD
for(i in 1:nrow(pref)){
  for(j in 1:length(rownames(filt.TE.sig.ond))){
    if(rownames(pref)[i] == rownames(filt.TE.sig.ond)[j]){
      pref$X9[i] = as.numeric(filt.TE.sig.ond$padj[j])
    }
  }
}
#TD vs FTLD
for(i in 1:nrow(pref)){
  for(j in 1:length(rownames(filt.COND.sig))){
    if(rownames(pref)[i] == rownames(filt.COND.sig)[j]){
      pref$X10[i] = as.numeric(filt.COND.sig$padj[j])
    }
  }
}

my36 = rownames(ZExpr)

for(i in 1:length(my36)){
  
  currentfeat = my36[i]
  
  pdat$Value[which(pdat$Response == currentfeat)] = pref[which(rownames(pref) == currentfeat),] #Columns already in the correct order
  
}

#Clean up pdat

minp = pdat
a = which(is.na(minp$Value))
b = which(minp$Value == 0)
c = c(a,b)
minp = minp[-c,]
minp = min(as.numeric(minp$Value))

cleanp = pdat
for(i in 1:nrow(cleanp)){
  if(is.na(cleanp$Value[i])){
    cleanp$Value[i] = 1
  }
  
  if(cleanp$Value[i] == 0){
    cleanp$Value[i] = minp
  }
}

cleanp$Value = as.numeric(unlist(cleanp$Value))

cleanp$Value = -log10(cleanp$Value)

cleanp$Response = as.character(cleanp$Response)
cleanp$Response = factor(cleanp$Response,levels=c("NANOGP4","MIRLET7BHG","MIR24-2","LINC01347","HSP90AB4P","GATA2-AS1","ENSG00000273151","ENSG00000258674","ENSG00000205041","COL3A1","CHKB-CPT1B","AGPAT4-IT1","UCP2","UBQLN2","TCIRG1","SLC17A6","SLC6A13","SERPINI1","OXR1","HTR2A","GLRA3","GAD2","GABRA1","COL18A1","TYROBP","TREM2","TNC","TMEM125","TLR7","HLA-DRA","FOLH1","CX3CR1","CHI3L2","CD44","APOC2","AIF1"))

cleanp$Module = as.character(cleanp$Module)
cleanp$Module = factor(cleanp$Module,levels=make.unique(cleanp$Module))

#Full unadjusted
p = ggplot(cleanp,aes(Module,Response)) + geom_tile(aes(fill=Value),color="black")
p = p+scale_fill_gradient2(low="blue",mid="white",high="darkgreen",limits=c(0,305))
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
p

#Subtype vs Subtype
a = which(cleanp$Module == "Glia_vs_Ox")
b = which(cleanp$Module == "Glia_vs_TD")
c = which(cleanp$Module == "Ox_vs_TD")
d = c(a,b,c)
subp = cleanp[d,]
p = ggplot(subp,aes(Module,Response)) + geom_tile(aes(fill=Value),color="black")
p = p+scale_fill_gradient2(low="blue",mid="white",high="darkgreen",limits=c(0,60))
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
p

#Color non-significant comparisons gray
p = ggplot(subp,aes(Module,Response)) + geom_tile(aes(fill=Value),color="black")
p = p+scale_fill_gradient2(low="blue",mid="white",high="darkgreen",limits=c(0.1,60)) #quick way - verified with table(subp$Value)
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
p


#Subtype vs HC
a = which(cleanp$Module == "Glia_vs_HC")
b = which(cleanp$Module == "Ox_vs_HC")
c = which(cleanp$Module == "TD_vs_HC")
d = c(a,b,c)
hcp = cleanp[d,]
p = ggplot(hcp,aes(Module,Response)) + geom_tile(aes(fill=Value),color="black")
p = p+scale_fill_gradient2(low="blue",mid="white",high="darkgreen",limits=c(0,305))
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
p

#Color non-significant comparisons gray
p = ggplot(hcp,aes(Module,Response)) + geom_tile(aes(fill=Value),color="black")
p = p+scale_fill_gradient2(low="blue",mid="white",high="darkgreen",limits=c(0.1,305)) #quick way - verified with table(subp$Value)
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
p

#Subtype vs FTLD
a = which(cleanp$Module == "Glia_vs_FTLD")
b = which(cleanp$Module == "Ox_vs_FTLD")
c = which(cleanp$Module == "TD_vs_FTLD")
d = c(a,b,c)
ftldp = cleanp[d,]
p = ggplot(ftldp,aes(Module,Response)) + geom_tile(aes(fill=Value),color="black")
p = p+scale_fill_gradient2(low="blue",mid="white",high="darkgreen",limits=c(0,250))
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
p

#Color non-significant comparisons gray
p = ggplot(ftldp,aes(Module,Response)) + geom_tile(aes(fill=Value),color="black")
p = p+scale_fill_gradient2(low="blue",mid="white",high="darkgreen",limits=c(0.1,250)) #quick way - verified with table(subp$Value)
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
p

# pheatdat = cleanp
# #Adjust extreme p-values for heatmap scale
# 
# for(i in 1:nrow(pheatdat)){
#   if(pheatdat$Value[i] > 40){
#     pheatdat$Value[i] = 40
#   }
# }
# 
# pheatdat$Response = as.character(pheatdat$Response)
# pheatdat$Response = factor(pheatdat$Response,levels=c("NANOGP4","MIRLET7BHG","MIR24-2","LINC01347","HSP90AB4P","GATA2-AS1","ENSG00000273151","ENSG00000258674","ENSG00000205041","COL3A1","CHKB-CPT1B","AGPAT4-IT1","UCP2","UBQLN2","TCIRG1","SLC17A6","SLC6A13","SERPINI1","OXR1","HTR2A","GLRA3","GAD2","GABRA1","COL18A1","TYROBP","TREM2","TNC","TMEM125","TLR7","HLA-DRA","FOLH1","CX3CR1","CHI3L2","CD44","APOC2","AIF1"))
# 
# pheatdat$Module = as.character(pheatdat$Module)
# pheatdat$Module = factor(pheatdat$Module,levels=make.unique(pheatdat$Module))
# 
# 
# p = ggplot(pheatdat,aes(Module,Response)) + geom_tile(aes(fill=Value),color="black")
# p = p+scale_fill_gradient2(low="blue",mid="white",high="darkgreen",limits=c(0,40))
# p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
# p
# 
