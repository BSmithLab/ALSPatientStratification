#Molecular subtypes of ALS are associated with differences in patient prognosis
#Eshima, O'Connor, Marschall, NYGC ALS Consortium, Bowser, Plaisier, Smith

#Nature Comms - Peer Review Round 1
#7/4/22

######################################################################################################################################################################################
######################################################################################################################################################################################
#################################  Reviewer 1  #######################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################

library(DESeq2)
library(stringr)
library(ggplot2)
library(beeswarm)

######################################################################################################################################################################################
#################################  QUESTION 1  #######################################################################################################################################
######################################################################################################################################################################################

#There are several non-disease factors that may partially explain the subtype classifications that I would like the authors to address:

#################################### RNA Quality Metrics

##########  RIN   #############
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping")
Meta = read.csv("GSE153960_MetaData.txt")
Clinical = read.csv("CLINICAL_DATA_PRUDENCIO.csv")

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Pub Data/Supplemental Files")
SRRs = read.csv("SRR_IDs.csv")
colnames(SRRs) = SRRs[1,]
SRRs = SRRs[-1,]

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

#Add RIN values to Supplemental Table
EshimaMetaData$RIN = rep(NA,nrow(EshimaMetaData))

for(i in 1:nrow(EshimaMetaData)){
  
  for(j in 1:nrow(Clinical)){
    
    if(EshimaMetaData$sample_id_alt[i] == Clinical$ExternalSampleId[j]){
      
      EshimaMetaData$RIN[i] = Clinical$RIN[j] 
      
    }
    
  }
  
}

ALS.RIN.Avg = mean(EshimaMetaData$RIN,na.rm=T)

which(is.na(EshimaMetaData$RIN)) #one missing value from Target ALS / NYGC


#RIN parsed by subtype
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping")
SubtypeLabels = read.csv("ALS451_coldata_SUBTYPES.csv")

convertnames = gsub("-","\\.",EshimaMetaData$sample_id_alt)
tmpMetaData = EshimaMetaData
tmpMetaData$sample_id_alt = convertnames

nsubtype = table(SubtypeLabels$Subtype)

Glia.RIN = rep(NA,nsubtype[[1]])
Ox.RIN = rep(NA,nsubtype[[2]])
TD.RIN = rep(NA,nsubtype[[3]])

count1 = count2 = count3 = 1

for(i in 1: nrow(tmpMetaData)){
  
  for(j in 1:nrow(SubtypeLabels)){
    
    if(tmpMetaData$sample_id_alt[i] == SubtypeLabels$X[j] && SubtypeLabels$Subtype[j] == "GLIA"){
      Glia.RIN[count1] = tmpMetaData$RIN[i]
      count1 = count1 + 1
    }else if(tmpMetaData$sample_id_alt[i] == SubtypeLabels$X[j] && SubtypeLabels$Subtype[j] == "OX"){
      Ox.RIN[count2] = tmpMetaData$RIN[i]
      count2 = count2 + 1
    }else if(tmpMetaData$sample_id_alt[i] == SubtypeLabels$X[j] && SubtypeLabels$Subtype[j] == "TE"){
      TD.RIN[count3] = tmpMetaData$RIN[i]
      count3 = count3 + 1
    }
    
  }
}

Glia.RIN.Avg = mean(Glia.RIN,na.rm = T)
Ox.RIN.Avg = mean(Ox.RIN,na.rm = T)
TD.RIN.Avg = mean(TD.RIN, na.rm = T)

par(mfrow=c(1,3))
#Normal Distribution isn't a terrible approximation...
hist(Glia.RIN,main="Glia RIN",breaks = 20);hist(Ox.RIN,main="Ox RIN",breaks = 20);hist(TD.RIN,main="TD RIN",breaks = 20)
par(mfrow=c(1,1))

#Stats
GO.RIN.sig = t.test(Glia.RIN,Ox.RIN) #p = 0.44
GO.RIN.sig = GO.RIN.sig$p.value
GT.RIN.sig = t.test(Glia.RIN,TD.RIN) #p = 6E-8
GT.RIN.sig = GT.RIN.sig$p.value
OT.RIN.sig = t.test(Ox.RIN,TD.RIN) #p = 7E-16
OT.RIN.sig = OT.RIN.sig$p.value

#Summary Boxplot
boxplot(Glia.RIN,Ox.RIN,TD.RIN,main="Subtype RIN",ylab="RIN",ylim = c(0,12),cex.lab = 1.5,cex.main=1.5,cex.axis=1.5,col=c("goldenrod1","navy","firebrick"),xaxt="n")
axis(at=1:3,side=1,labels=c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=1.5)
lines(c(1,2),c(10,10),type="l")
lines(c(1,1),c(10,9.5),type="l")
lines(c(2,2),c(10,9.5),type="l")
text(mean(c(1,2)),10.3,labels = paste("p-value:",formatC(GO.RIN.sig,format = "e",digits = 2)))
lines(c(2,3),c(10.75,10.75),type="l")
lines(c(2,2),c(10.75,10.25),type="l")
lines(c(3,3),c(10.75,10.25),type="l")
text(mean(c(2,3)),11.05,labels = paste("p-value:",formatC(OT.RIN.sig,format = "e",digits = 2)))
lines(c(1,3),c(11.5,11.5),type="l")
lines(c(1,1),c(11.5,11),type="l")
lines(c(3,3),c(11.5,11),type="l")
text(mean(c(1,3)),11.8,labels = paste("p-value:",formatC(GT.RIN.sig,format = "e",digits = 2)))

#####################################################################################################################

#Univariate analysis with RIN as covariate (design = ~platform + RIN + Subtype )

#Design Equation 1: design = ~ platform + Subtype (from UnivariateAnalysis script)
#Design Equation 2: design = ~ platform + RIN + Subtype

load("D:/Jarrett/Research/Fall 2021/ProgBM/SOD1/ALSPatientStratification_UnivariateDatasets_SOD1.RData")

#Filter for subtype-specific features in Fig. 6
my36 = c("AIF1","APOC2","CD44","CHI3L2","CX3CR1","FOLH1","HLA-DRA","TLR7","TMEM125","TNC","TREM2","TYROBP","COL18A1","GABRA1","GAD2","GLRA3","HTR2A","OXR1","SERPINI1","SLC6A13","SLC17A6","TCIRG1","UBQLN2","UCP2","AGPAT4-IT1","CHKB-CPT1B","COL3A1","ENSG00000205041","ENSG00000258674","ENSG00000273151","GATA2-AS1","HSP90AB4P","LINC01347","MIR24-2","MIRLET7BHG","NANOGP4")


#Clean Up Design 1
filt.glia.d1 = glia.res[rownames(glia.res) %in% my36,]
filt.ox.d1 = ox.res[rownames(ox.res) %in% my36,]
filt.TE.d1 = TE.res[rownames(TE.res) %in% my36,]
filt.GT.d1 = GT.res[rownames(GT.res) %in% my36,]
filt.GO.d1 = GO.res[rownames(GO.res) %in% my36,]
filt.TO.d1 = TO.res[rownames(TO.res) %in% my36,]
filt.glia.ond.d1 = glia.res.ond[rownames(glia.res.ond) %in% my36,]
filt.ox.ond.d1 = ox.res.ond[rownames(ox.res.ond) %in% my36,]
filt.TE.ond.d1 = TE.res.ond[rownames(TE.res.ond) %in% my36,]
filt.COND.d1 = COND.res[rownames(COND.res) %in% my36,]


#Add in RIN
FullPheno$RIN = NA
Clinical2 = Clinical
convertnames = gsub("-","\\.",Clinical2$ExternalSampleId)
Clinical2$ExternalSampleId = convertnames

for(i in 1:nrow(FullPheno)){
  
  for(j in 1:nrow(Clinical2)){
    
    if(FullPheno$Subject[i] == Clinical2$ExternalSampleId[j]){
      
      FullPheno$RIN[i] = Clinical2$RIN[j]
      
    }
    
  }
  
}

table(FullPheno$Subject == colnames(FullCount))

#Remove single sample with missing RIN
FullPheno_rin = FullPheno[-which(is.na(FullPheno$RIN)),]
FullCount_rin = FullCount[,-which(is.na(FullPheno$RIN))]
table(FullPheno_rin$Subject == colnames(FullCount_rin))

rCountData_rin = round(FullCount_rin,0)

dds_rin = DESeqDataSetFromMatrix(countData = rCountData_rin, colData = FullPheno_rin, design= ~ platform + RIN + Subtype, tidy=F) #Subtype must be last (DESeq2 vignette)
dds_rin$Subtype = relevel(dds_rin$Subtype,ref = "Control")
dseq_rin = DESeq(dds_rin,betaPrior=T)


#Pairwise "contrast()" for Design Equation 2
glia.res = results(dseq_rin,contrast = c("Subtype","GLIA","Control"))
filt.glia.d2 = glia.res[rownames(glia.res) %in% my36,]

ox.res = results(dseq_rin,contrast = c("Subtype","OX","Control"))
filt.ox.d2 = ox.res[rownames(ox.res) %in% my36,]

TE.res = results(dseq_rin,contrast = c("Subtype","TE","Control"))
filt.TE.d2 = TE.res[rownames(TE.res) %in% my36,]

GT.res = results(dseq_rin,contrast = c("Subtype","GLIA","TE"))
filt.GT.d2 = GT.res[rownames(GT.res) %in% my36,]

GO.res = results(dseq_rin,contrast = c("Subtype","GLIA","OX"))
filt.GO.d2 = GO.res[rownames(GO.res) %in% my36,]

TO.res = results(dseq_rin,contrast = c("Subtype","TE","OX"))
filt.TO.d2 = TO.res[rownames(TO.res) %in% my36,]

glia.res.ond = results(dseq_rin,contrast = c("Subtype","GLIA","OND"))
filt.glia.ond.d2 = glia.res.ond[rownames(glia.res.ond) %in% my36,]

ox.res.ond = results(dseq_rin,contrast = c("Subtype","OX","OND"))
filt.ox.ond.d2 = ox.res.ond[rownames(ox.res.ond) %in% my36,]

TE.res.ond = results(dseq_rin,contrast = c("Subtype","TE","OND"))
filt.TE.ond.d2 = TE.res.ond[rownames(TE.res.ond) %in% my36,]

COND.res = results(dseq_rin,contrast = c("Subtype","Control","OND"))
filt.COND.d2 = COND.res[rownames(COND.res) %in% my36,]


#Build design equation p-value matrix
npairs = 10
ndesign = 2
DE_RIN = data.frame(matrix(NA,length(my36),npairs*ndesign))
colnames(DE_RIN) = c("Design1_ALS-Glia_v_controls","Design2_ALS-Glia_v_controls","Design1_ALS-Ox_v_controls","Design2_ALS-Ox_v_controls","Design1_ALS-TD_v_controls","Design2_ALS-TD_v_controls","Design1_ALS-Glia_v_ALS-Ox","Design2_ALS-Glia_v_ALS-Ox","Design1_ALS-Glia_v_ALS-TD","Design2_ALS-Glia_v_ALS-TD","Design1_ALS-Ox_v_ALS-TD","Design2_ALS-Ox_v_ALS-TD","Design1_ALS-Glia_v_FTLD","Design2_ALS-Glia_v_FTLD","Design1_ALS-Ox_v_FTLD","Design2_ALS-Ox_v_FTLD","Design1_ALS-TD_v_FTLD","Design2_ALS-TD_v_FTLD","Design1_Control_v_FTLD","Design2_Control_v_FTLD")
rownames(DE_RIN) = my36


#Reference List to help with condensing dataframes (same order as DE_RIN matrix)
RefList = c("filt.glia.d1","filt.glia.d2","filt.ox.d1","filt.ox.d2","filt.TE.d1","filt.TE.d2","filt.GO.d1","filt.GO.d2","filt.GT.d1","filt.GT.d2","filt.TO.d1","filt.TO.d2","filt.glia.ond.d1","filt.glia.ond.d2","filt.ox.ond.d1","filt.ox.ond.d2","filt.TE.ond.d1","filt.TE.ond.d2","filt.COND.d1","filt.COND.d2")

for(i in 1:length(RefList)){
  
  tmp = get(RefList[i])
  
  for(j in 1:nrow(DE_RIN)){
    
    for(k in 1:nrow(tmp)){
      
      if(rownames(DE_RIN)[j] == rownames(tmp)[k]){
        
        DE_RIN[j,i] = tmp$padj[k]
        
      }
      
    }
    
  }
  
  
}


#Write out results

#write.csv(DE_RIN,"DifferentialExpression_withRINcovariate_PeerReviewRound1.csv")

#Adjusting DE for RIN variability does not seem to have much of an impact on the significance



##########  3'/5' bias and % Exon Coverage ###########
#Partial Unix script is provided below (commented) 

# for filename in *.sam
# do
# echo $filename
# java -jar /home/jeshima/miniconda3/envs/Picard/share/picard-2.27.4-0/picard.jar CollectRnaSeqMetrics -REF_FLAT refFlat.txt -I $filename -O ${filename%.*}_metrics.txt -STRAND_SPECIFICITY NONE 
# done

#############################################################################################
#############################################################################################

##########  Submitting site/hospital #############
nsubtype = table(SubtypeLabels$Subtype)

Glia.Site = rep(NA,nsubtype[[1]])
Ox.Site = rep(NA,nsubtype[[2]])
TD.Site = rep(NA,nsubtype[[3]])

count1 = count2 = count3 = 1

for(i in 1: nrow(tmpMetaData)){
  
  for(j in 1:nrow(SubtypeLabels)){
    
    if(tmpMetaData$sample_id_alt[i] == SubtypeLabels$X[j] && SubtypeLabels$Subtype[j] == "GLIA"){
      Glia.Site[count1] = tmpMetaData$project[i]
      count1 = count1 + 1
    }else if(tmpMetaData$sample_id_alt[i] == SubtypeLabels$X[j] && SubtypeLabels$Subtype[j] == "OX"){
      Ox.Site[count2] = tmpMetaData$project[i]
      count2 = count2 + 1
    }else if(tmpMetaData$sample_id_alt[i] == SubtypeLabels$X[j] && SubtypeLabels$Subtype[j] == "TE"){
      TD.Site[count3] = tmpMetaData$project[i]
      count3 = count3 + 1
    }
    
  }
}

par(mar=c(5,5,4,2))
tmp = table(Glia.Site)
barplot(c(tmp[[1]],tmp[[2]]),ylim=c(0,50),main = "ALS-Glia Sites",xlab="Site",col = c("aquamarine","coral1"),cex=2,cex.axis = 2,cex.main=2,cex.lab=2)
axis(1,c(0.65,1.9),labels = c("NYGC","Target ALS"),cex.axis=2)
title(ylab = "Frequency",mgp = c(3.5,1,0),cex.lab=2)
text(0.65,tmp[[1]]+0.05*tmp[[1]],labels = tmp[[1]],cex=2)
text(1.9,tmp[[2]]+0.05*tmp[[2]],labels = tmp[[2]],cex=2)

tmp = table(Ox.Site)
barplot(c(tmp[[1]],tmp[[2]]),ylim=c(0,250),main = "ALS-Ox Sites",xlab="Site",col = c("aquamarine","coral1"),cex=2,cex.axis = 2,cex.main=2,cex.lab=2)
axis(1,c(0.65,1.9),labels = c("NYGC","Target ALS"),cex.axis=2)
title(ylab = "Frequency",mgp = c(3.5,1,0),cex.lab=2)
text(0.65,tmp[[1]]+0.3*tmp[[1]],labels = tmp[[1]],cex=2)
text(1.9,tmp[[2]]+0.05*tmp[[2]],labels = tmp[[2]],cex=2)

tmp = table(TD.Site)
barplot(c(tmp[[1]],tmp[[2]]),ylim=c(0,100),main = "ALS-TD Sites",xlab="Site",col = c("aquamarine","coral1"),cex=2,cex.axis = 2,cex.main=2,cex.lab=2)
axis(1,c(0.65,1.9),labels = c("NYGC","Target ALS"),cex.axis=2)
title(ylab = "Frequency",mgp = c(3.5,1,0),cex.lab=2)
text(0.65,tmp[[1]]+0.1*tmp[[1]],labels = tmp[[1]],cex=2)
text(1.9,tmp[[2]]+0.05*tmp[[2]],labels = tmp[[2]],cex=2)

#Nicer plot
siteplot = data.frame(matrix(NA,nrow=6,ncol=3))
colnames(siteplot) = c("Subtype","Site","Frequency")
siteplot$Subtype = c(rep("ALS-Glia",2),rep("ALS-Ox",2),rep("ALS-TD",2))
siteplot$Site = rep(c("NYGC","Target ALS"),3)
siteplot$Frequency[1:2] = as.numeric(table(Glia.Site))
siteplot$Frequency[3:4] = as.numeric(table(Ox.Site))
siteplot$Frequency[5:6] = as.numeric(table(TD.Site))

p = ggplot(data=siteplot,aes(x=Site,y=Frequency,fill=Subtype)) + geom_bar(stat = "identity",position = position_dodge())
p = p+scale_fill_manual(values = c("goldenrod1","navy","firebrick"))
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "gray80"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "white"))
p = p+theme(axis.title.y=element_text(angle=90, vjust=4,size=24))
p = p+theme(plot.margin = unit(c(1,1,1,1), "cm"))
p = p+theme(axis.text.x = element_text(size = 24))
p = p+theme(axis.text.y = element_text(size = 24))
p = p+ xlab("")
p = p+ theme(legend.text=element_text(size=24),legend.title=element_text(size=24))
p

#################################################################################

#Univariate analysis with Site as covariate (design = ~platform + Site + Subtype )

#Design Equation 1: design = ~ platform + Subtype (from UnivariateAnalysis script)
#Design Equation 2: design = ~ platform + Site + Subtype

load("D:/Jarrett/Research/Fall 2021/ProgBM/SOD1/ALSPatientStratification_UnivariateDatasets_SOD1.RData")

#Filter for subtype-specific features in Fig. 6
my36 = c("AIF1","APOC2","CD44","CHI3L2","CX3CR1","FOLH1","HLA-DRA","TLR7","TMEM125","TNC","TREM2","TYROBP","COL18A1","GABRA1","GAD2","GLRA3","HTR2A","OXR1","SERPINI1","SLC6A13","SLC17A6","TCIRG1","UBQLN2","UCP2","AGPAT4-IT1","CHKB-CPT1B","COL3A1","ENSG00000205041","ENSG00000258674","ENSG00000273151","GATA2-AS1","HSP90AB4P","LINC01347","MIR24-2","MIRLET7BHG","NANOGP4")


#Clean Up Design 1
filt.glia.d1 = glia.res[rownames(glia.res) %in% my36,]
filt.ox.d1 = ox.res[rownames(ox.res) %in% my36,]
filt.TE.d1 = TE.res[rownames(TE.res) %in% my36,]
filt.GT.d1 = GT.res[rownames(GT.res) %in% my36,]
filt.GO.d1 = GO.res[rownames(GO.res) %in% my36,]
filt.TO.d1 = TO.res[rownames(TO.res) %in% my36,]
filt.glia.ond.d1 = glia.res.ond[rownames(glia.res.ond) %in% my36,]
filt.ox.ond.d1 = ox.res.ond[rownames(ox.res.ond) %in% my36,]
filt.TE.ond.d1 = TE.res.ond[rownames(TE.res.ond) %in% my36,]
filt.COND.d1 = COND.res[rownames(COND.res) %in% my36,]

#Add Site to Pheno
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping")
Meta = read.csv("GSE153960_MetaData.txt")
convertnames = gsub("-","\\.",Meta$sample_id_alt)
Meta$sample_id_alt = convertnames

FiltMeta = Meta[Meta$sample_id_alt %in% FullPheno$Subject,]


reftable = table(FiltMeta$sample_id_alt)

removeind = rep(NA,nrow(FiltMeta)-length(reftable))
count = 1
for(i in 1:length(reftable)) {
  
  if(reftable[[i]] > 1){
    
    tmp = names(reftable[i])
    inds = which(FiltMeta$sample_id_alt == tmp)
    
    if(length(inds) == 2){
      removeind[count] = inds[1]
      count = count+1
    }else if(length(inds) == 3){
      removeind[seq(count,count+1,1)] = inds[1:2]
      count = count+2
    }else if(length(inds) > 3){
      cat("Problem")
    }
    
    
  }
  
}

SiteMeta = FiltMeta[-removeind,]



FullPheno$Site = NA

for(i in 1:nrow(FullPheno)){
  
  for(j in 1:nrow(SiteMeta)){
    
    if(FullPheno$Subject[i] == SiteMeta$sample_id_alt[j]){
      
      FullPheno$Site[i] = SiteMeta$project[j]
      
    }
    
  }
  
}

#Clean up site names
for(i in 1:nrow(FullPheno)){
  
  if(FullPheno$Site[i] == "NYGC ALS Consortium"){
    FullPheno$Site[i] = "NYGC"
  }else{
    FullPheno$Site[i] = "TargetALS"
  }
  
}

table(FullPheno$Site)

table(FullPheno$Subject == colnames(FullCount))


rCountData_site = round(FullCount,0)

dds_site = DESeqDataSetFromMatrix(countData = rCountData_site, colData = FullPheno, design= ~ platform + Site + Subtype, tidy=F) #Subtype must be last (DESeq2 vignette)
dds_site$Subtype = relevel(dds_site$Subtype,ref = "Control")
dseq_site = DESeq(dds_site,betaPrior=T)


#Pairwise "contrast()" for Design Equation 2
glia.res = results(dseq_site,contrast = c("Subtype","GLIA","Control"))
filt.glia.d2 = glia.res[rownames(glia.res) %in% my36,]

ox.res = results(dseq_site,contrast = c("Subtype","OX","Control"))
filt.ox.d2 = ox.res[rownames(ox.res) %in% my36,]

TE.res = results(dseq_site,contrast = c("Subtype","TE","Control"))
filt.TE.d2 = TE.res[rownames(TE.res) %in% my36,]

GT.res = results(dseq_site,contrast = c("Subtype","GLIA","TE"))
filt.GT.d2 = GT.res[rownames(GT.res) %in% my36,]

GO.res = results(dseq_site,contrast = c("Subtype","GLIA","OX"))
filt.GO.d2 = GO.res[rownames(GO.res) %in% my36,]

TO.res = results(dseq_site,contrast = c("Subtype","TE","OX"))
filt.TO.d2 = TO.res[rownames(TO.res) %in% my36,]

glia.res.ond = results(dseq_site,contrast = c("Subtype","GLIA","OND"))
filt.glia.ond.d2 = glia.res.ond[rownames(glia.res.ond) %in% my36,]

ox.res.ond = results(dseq_site,contrast = c("Subtype","OX","OND"))
filt.ox.ond.d2 = ox.res.ond[rownames(ox.res.ond) %in% my36,]

TE.res.ond = results(dseq_site,contrast = c("Subtype","TE","OND"))
filt.TE.ond.d2 = TE.res.ond[rownames(TE.res.ond) %in% my36,]

COND.res = results(dseq_site,contrast = c("Subtype","Control","OND"))
filt.COND.d2 = COND.res[rownames(COND.res) %in% my36,]


#Build design equation p-value matrix
npairs = 10
ndesign = 2
DE_SITE = data.frame(matrix(NA,length(my36),npairs*ndesign))
colnames(DE_SITE) = c("Design1_ALS-Glia_v_controls","Design2_ALS-Glia_v_controls","Design1_ALS-Ox_v_controls","Design2_ALS-Ox_v_controls","Design1_ALS-TD_v_controls","Design2_ALS-TD_v_controls","Design1_ALS-Glia_v_ALS-Ox","Design2_ALS-Glia_v_ALS-Ox","Design1_ALS-Glia_v_ALS-TD","Design2_ALS-Glia_v_ALS-TD","Design1_ALS-Ox_v_ALS-TD","Design2_ALS-Ox_v_ALS-TD","Design1_ALS-Glia_v_FTLD","Design2_ALS-Glia_v_FTLD","Design1_ALS-Ox_v_FTLD","Design2_ALS-Ox_v_FTLD","Design1_ALS-TD_v_FTLD","Design2_ALS-TD_v_FTLD","Design1_Control_v_FTLD","Design2_Control_v_FTLD")
rownames(DE_SITE) = my36


#Reference List to help with condensing dataframes (same order as DE_SITE matrix)
RefList = c("filt.glia.d1","filt.glia.d2","filt.ox.d1","filt.ox.d2","filt.TE.d1","filt.TE.d2","filt.GO.d1","filt.GO.d2","filt.GT.d1","filt.GT.d2","filt.TO.d1","filt.TO.d2","filt.glia.ond.d1","filt.glia.ond.d2","filt.ox.ond.d1","filt.ox.ond.d2","filt.TE.ond.d1","filt.TE.ond.d2","filt.COND.d1","filt.COND.d2")

for(i in 1:length(RefList)){
  
  tmp = get(RefList[i])
  
  for(j in 1:nrow(DE_SITE)){
    
    for(k in 1:nrow(tmp)){
      
      if(rownames(DE_SITE)[j] == rownames(tmp)[k]){
        
        DE_SITE[j,i] = tmp$padj[k]
        
      }
      
    }
    
  }
  
  
}


#Write out results

#write.csv(DE_SITE,"DifferentialExpression_withSITEcovariate_PeerReviewRound1.csv")

#Adjusting DE for Site does not seem to have much of an impact on the significance (TE vs Ox changes)

###############################################################################
###############################################################################

########## Site & RIN as extra covariates ##########

#Design Equation 1: design = ~ platform + Subtype (from UnivariateAnalysis script)
#Design Equation 2: design = ~ platform + Site + RIN + Subtype

load("D:/Jarrett/Research/Fall 2021/ProgBM/SOD1/ALSPatientStratification_UnivariateDatasets_SOD1.RData")

#Filter for subtype-specific features in Fig. 6
my36 = c("AIF1","APOC2","CD44","CHI3L2","CX3CR1","FOLH1","HLA-DRA","TLR7","TMEM125","TNC","TREM2","TYROBP","COL18A1","GABRA1","GAD2","GLRA3","HTR2A","OXR1","SERPINI1","SLC6A13","SLC17A6","TCIRG1","UBQLN2","UCP2","AGPAT4-IT1","CHKB-CPT1B","COL3A1","ENSG00000205041","ENSG00000258674","ENSG00000273151","GATA2-AS1","HSP90AB4P","LINC01347","MIR24-2","MIRLET7BHG","NANOGP4")


#Clean Up Design 1
filt.glia.d1 = glia.res[rownames(glia.res) %in% my36,]
filt.ox.d1 = ox.res[rownames(ox.res) %in% my36,]
filt.TE.d1 = TE.res[rownames(TE.res) %in% my36,]
filt.GT.d1 = GT.res[rownames(GT.res) %in% my36,]
filt.GO.d1 = GO.res[rownames(GO.res) %in% my36,]
filt.TO.d1 = TO.res[rownames(TO.res) %in% my36,]
filt.glia.ond.d1 = glia.res.ond[rownames(glia.res.ond) %in% my36,]
filt.ox.ond.d1 = ox.res.ond[rownames(ox.res.ond) %in% my36,]
filt.TE.ond.d1 = TE.res.ond[rownames(TE.res.ond) %in% my36,]
filt.COND.d1 = COND.res[rownames(COND.res) %in% my36,]

#Add Site to Pheno
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping")
Meta = read.csv("GSE153960_MetaData.txt")
Clinical = read.csv("CLINICAL_DATA_PRUDENCIO.csv")
convertnames = gsub("-","\\.",Meta$sample_id_alt)
Meta$sample_id_alt = convertnames

#Add in RIN
FullPheno$RIN = NA
Clinical2 = Clinical
convertnames = gsub("-","\\.",Clinical2$ExternalSampleId)
Clinical2$ExternalSampleId = convertnames

for(i in 1:nrow(FullPheno)){
  
  for(j in 1:nrow(Clinical2)){
    
    if(FullPheno$Subject[i] == Clinical2$ExternalSampleId[j]){
      
      FullPheno$RIN[i] = Clinical2$RIN[j]
      
    }
    
  }
  
}

table(FullPheno$Subject == colnames(FullCount))



FiltMeta = Meta[Meta$sample_id_alt %in% FullPheno$Subject,]


reftable = table(FiltMeta$sample_id_alt)

removeind = rep(NA,nrow(FiltMeta)-length(reftable))
count = 1
for(i in 1:length(reftable)) {
  
  if(reftable[[i]] > 1){
    
    tmp = names(reftable[i])
    inds = which(FiltMeta$sample_id_alt == tmp)
    
    if(length(inds) == 2){
      removeind[count] = inds[1]
      count = count+1
    }else if(length(inds) == 3){
      removeind[seq(count,count+1,1)] = inds[1:2]
      count = count+2
    }else if(length(inds) > 3){
      cat("Problem")
    }
    
    
  }
  
}

SiteMeta = FiltMeta[-removeind,]



FullPheno$Site = NA

for(i in 1:nrow(FullPheno)){
  
  for(j in 1:nrow(SiteMeta)){
    
    if(FullPheno$Subject[i] == SiteMeta$sample_id_alt[j]){
      
      FullPheno$Site[i] = SiteMeta$project[j]
      
    }
    
  }
  
}

#Clean up site names
for(i in 1:nrow(FullPheno)){
  
  if(FullPheno$Site[i] == "NYGC ALS Consortium"){
    FullPheno$Site[i] = "NYGC"
  }else{
    FullPheno$Site[i] = "TargetALS"
  }
  
}

table(FullPheno$Site)

#Single sample with missing RIN must be removed (incomplete design equation)
table(FullPheno$Subject == colnames(FullCount))
FullPheno_sr = FullPheno[-which(is.na(FullPheno$RIN)),]
FullCount_sr = FullCount[,-which(is.na(FullPheno$RIN))]

FullPheno_sr$RIN = scale(FullPheno_sr$RIN,center = T)

rCountData_rinsite = round(FullCount_sr,0)

dds_rinsite = DESeqDataSetFromMatrix(countData = rCountData_rinsite, colData = FullPheno_sr, design= ~ platform + Site + RIN + Subtype, tidy=F) #Subtype must be last (DESeq2 vignette)
dds_rinsite$Subtype = relevel(dds_rinsite$Subtype,ref = "Control")
dseq_rinsite = DESeq(dds_rinsite,betaPrior=T)


#Pairwise "contrast()" for Design Equation 2
glia.res = results(dseq_rinsite,contrast = c("Subtype","GLIA","Control"))
filt.glia.d2 = glia.res[rownames(glia.res) %in% my36,]

ox.res = results(dseq_rinsite,contrast = c("Subtype","OX","Control"))
filt.ox.d2 = ox.res[rownames(ox.res) %in% my36,]

TE.res = results(dseq_rinsite,contrast = c("Subtype","TE","Control"))
filt.TE.d2 = TE.res[rownames(TE.res) %in% my36,]

GT.res = results(dseq_rinsite,contrast = c("Subtype","GLIA","TE"))
filt.GT.d2 = GT.res[rownames(GT.res) %in% my36,]

GO.res = results(dseq_rinsite,contrast = c("Subtype","GLIA","OX"))
filt.GO.d2 = GO.res[rownames(GO.res) %in% my36,]

TO.res = results(dseq_rinsite,contrast = c("Subtype","TE","OX"))
filt.TO.d2 = TO.res[rownames(TO.res) %in% my36,]

glia.res.ond = results(dseq_rinsite,contrast = c("Subtype","GLIA","OND"))
filt.glia.ond.d2 = glia.res.ond[rownames(glia.res.ond) %in% my36,]

ox.res.ond = results(dseq_rinsite,contrast = c("Subtype","OX","OND"))
filt.ox.ond.d2 = ox.res.ond[rownames(ox.res.ond) %in% my36,]

TE.res.ond = results(dseq_rinsite,contrast = c("Subtype","TE","OND"))
filt.TE.ond.d2 = TE.res.ond[rownames(TE.res.ond) %in% my36,]

COND.res = results(dseq_rinsite,contrast = c("Subtype","Control","OND"))
filt.COND.d2 = COND.res[rownames(COND.res) %in% my36,]


#Build design equation p-value matrix
npairs = 10
ndesign = 2
DE_RINSITE = data.frame(matrix(NA,length(my36),npairs*ndesign))
colnames(DE_RINSITE) = c("Design1_ALS-Glia_v_controls","Design2_ALS-Glia_v_controls","Design1_ALS-Ox_v_controls","Design2_ALS-Ox_v_controls","Design1_ALS-TD_v_controls","Design2_ALS-TD_v_controls","Design1_ALS-Glia_v_ALS-Ox","Design2_ALS-Glia_v_ALS-Ox","Design1_ALS-Glia_v_ALS-TD","Design2_ALS-Glia_v_ALS-TD","Design1_ALS-Ox_v_ALS-TD","Design2_ALS-Ox_v_ALS-TD","Design1_ALS-Glia_v_FTLD","Design2_ALS-Glia_v_FTLD","Design1_ALS-Ox_v_FTLD","Design2_ALS-Ox_v_FTLD","Design1_ALS-TD_v_FTLD","Design2_ALS-TD_v_FTLD","Design1_Control_v_FTLD","Design2_Control_v_FTLD")
rownames(DE_RINSITE) = my36


#Reference List to help with condensing dataframes (same order as DE_SITE matrix)
RefList = c("filt.glia.d1","filt.glia.d2","filt.ox.d1","filt.ox.d2","filt.TE.d1","filt.TE.d2","filt.GO.d1","filt.GO.d2","filt.GT.d1","filt.GT.d2","filt.TO.d1","filt.TO.d2","filt.glia.ond.d1","filt.glia.ond.d2","filt.ox.ond.d1","filt.ox.ond.d2","filt.TE.ond.d1","filt.TE.ond.d2","filt.COND.d1","filt.COND.d2")

for(i in 1:length(RefList)){
  
  tmp = get(RefList[i])
  
  for(j in 1:nrow(DE_RINSITE)){
    
    for(k in 1:nrow(tmp)){
      
      if(rownames(DE_RINSITE)[j] == rownames(tmp)[k]){
        
        DE_RINSITE[j,i] = tmp$padj[k]
        
      }
      
    }
    
  }
  
  
}


#Write out results

#write.csv(DE_RINSITE,"DE_withRINandSITEcovariates_PeerReview_ScaledRIN.csv")

#Adjusting DE for RIN and Site does not seem to have much of an impact on the significance


###############################################################################

##########  Brain Region   #############

for(i in 1:nrow(SubtypeLabels)){
  
  if(SubtypeLabels$tissue[i] == "Frontal Cortex"){
    SubtypeLabels$tissue[i] = "Cortex Frontal"
  }else if(SubtypeLabels$tissue[i] == "Lateral Motor Cortex"){
    SubtypeLabels$tissue[i] = "Cortex Motor Lateral"
  }else if(SubtypeLabels$tissue[i] == "Medial Motor Cortex"){
    SubtypeLabels$tissue[i] = "Cortex Motor Medial"
  }else if(SubtypeLabels$tissue[i] == "Other Motor Cortex"){
    SubtypeLabels$tissue[i] = "Cortex Motor Unspecified"
  }else{
    cat("Warning Flag")
  }
}

for(i in 1:nrow(SubtypeLabels)){
  if(SubtypeLabels$Subtype[i] == "TE"){
    SubtypeLabels$Subtype[i] = "TD"
  }
}

Tissuetable = table(tmpMetaData$tissue)

#[[1]] - Cortex Frontal
#[[2]] - Cortex Motor Lateral
#[[3]] - Cortex Motor Medial
#[[4]] - Cortex Motor Unspecified

#Containers
FrCor = rep(NA,Tissuetable[[1]])
LatCor = rep(NA,Tissuetable[[2]])
MedCor = rep(NA,Tissuetable[[3]])
Unsp = rep(NA,Tissuetable[[4]])

count1 = count2 = count3 = count4 = 1
for(i in 1:nrow(tmpMetaData)){
  
  for(j in 1:nrow(SubtypeLabels)){
    
    if(tmpMetaData$sample_id_alt[i] == SubtypeLabels$X[j] && SubtypeLabels$tissue[j] == "Cortex Frontal"){
      FrCor[count1] = SubtypeLabels$Subtype[j]
      count1 = count1+1
    }else if(tmpMetaData$sample_id_alt[i] == SubtypeLabels$X[j] && SubtypeLabels$tissue[j] == "Cortex Motor Lateral"){
      LatCor[count2] = SubtypeLabels$Subtype[j]
      count2 = count2+1
    }else if(tmpMetaData$sample_id_alt[i] == SubtypeLabels$X[j] && SubtypeLabels$tissue[j] == "Cortex Motor Medial"){
      MedCor[count3] = SubtypeLabels$Subtype[j]
      count3 = count3+1
    }else if(tmpMetaData$sample_id_alt[i] == SubtypeLabels$X[j] && SubtypeLabels$tissue[j] == "Cortex Motor Unspecified"){
      Unsp[count4] = SubtypeLabels$Subtype[j] 
      count4 = count4+1
    }
  }
}

#Frontal Cortex
par(mar=c(5,5,4,2))
tmp = table(FrCor)
barplot(c(tmp[[1]],tmp[[2]],tmp[[3]]),ylim=c(0,120),main = "Frontal Cortex Subtypes",xlab="Subtype",col = c("goldenrod1","navyblue","firebrick"),cex=2,cex.main=2,cex.axis=2,cex.lab=2)
axis(1,c(0.65,1.9,3.1),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=2)
title(ylab = "Frequency",mgp = c(3.5,1,0),cex.lab=2)
text(0.65,tmp[[1]]+0.15*tmp[[1]],labels = tmp[[1]],cex=2)
text(1.9,tmp[[2]]+0.05*tmp[[2]],labels = tmp[[2]],cex=2)
text(3.1,tmp[[3]]+0.1*tmp[[3]],labels = tmp[[3]],cex=2)

#Lateral Motor Cortex
tmp = table(LatCor)
barplot(c(tmp[[1]],tmp[[2]],tmp[[3]]),ylim=c(0,80),main = "Lateral Motor Cortex Subtypes",xlab="Subtype",col = c("goldenrod1","navyblue","firebrick"),cex=2,cex.main=2,cex.axis=2,cex.lab=2)
axis(1,c(0.65,1.9,3.1),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=2)
title(ylab = "Frequency",mgp = c(3.5,1,0),cex.lab=2)
text(0.65,tmp[[1]]+0.2*tmp[[1]],labels = tmp[[1]],cex=2)
text(1.9,tmp[[2]]+0.05*tmp[[2]],labels = tmp[[2]],cex=2)
text(3.1,tmp[[3]]+0.1*tmp[[3]],labels = tmp[[3]],cex=2)

#Medial Motor Cortex
tmp = table(MedCor)
barplot(c(tmp[[1]],tmp[[2]],tmp[[3]]),ylim=c(0,80),main = "Medial Motor Cortex Subtypes",xlab="Subtype",col = c("goldenrod1","navyblue","firebrick"),cex=2,cex.main=2,cex.axis=2,cex.lab=2)
axis(1,c(0.65,1.9,3.1),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=2)
title(ylab = "Frequency",mgp = c(3.5,1,0),cex.lab=2)
text(0.65,tmp[[1]]+0.2*tmp[[1]],labels = tmp[[1]],cex=2)
text(1.9,tmp[[2]]+0.05*tmp[[2]],labels = tmp[[2]],cex=2)
text(3.1,tmp[[3]]+0.1*tmp[[3]],labels = tmp[[3]],cex=2)

#Unspecified Motor Cortex
tmp = table(Unsp)
barplot(c(tmp[[1]],tmp[[2]],tmp[[3]]),ylim=c(0,40),main = "Unspecified Motor Cortex Subtypes",xlab="Subtype",col = c("goldenrod1","navyblue","firebrick"),cex=2,cex.main=2,cex.axis=2,cex.lab=2)
axis(1,c(0.65,1.9,3.1),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=2)
title(ylab = "Frequency",mgp = c(3.5,1,0),cex.lab=2)
text(0.65,tmp[[1]]+0.15*tmp[[1]],labels = tmp[[1]],cex=2)
text(1.9,tmp[[2]]+0.15*tmp[[2]],labels = tmp[[2]],cex=2)
text(3.1,tmp[[3]]+0.15*tmp[[3]],labels = tmp[[3]],cex=2)


######################################################################################################################################################################################
#################################  QUESTION 2  #######################################################################################################################################
######################################################################################################################################################################################

#Did the authors apply a cutoff for lowly expressed genes? My concern is that some of 
#the apparent subtype differences are driven by noisy gene expression. Can the authors rule this out?

load("D:/Jarrett/Research/Nat Comm Reviews/RData/ALSPatientStratification_UnivariateDatasets_RINSite_PeerReview.RData")

################## Raw Count Grid ##########################

my36 = c("AIF1","APOC2","CD44","CHI3L2","CX3CR1","FOLH1","HLA-DRA","TLR7","TMEM125","TNC","TREM2","TYROBP","COL18A1","GABRA1","GAD2","GLRA3","HTR2A","OXR1","SERPINI1","SLC6A13","SLC17A6","TCIRG1","UBQLN2","UCP2","AGPAT4-IT1","CHKB-CPT1B","COL3A1","ENSG00000205041","ENSG00000258674","ENSG00000273151","GATA2-AS1","HSP90AB4P","LINC01347","MIR24-2","MIRLET7BHG","NANOGP4")
RawCount36 = FullCount [rownames(FullCount) %in% my36,]
RawCount36 = RawCount36 [colnames(RawCount36) %in% ALSPheno$Subject]

LogCount36 = log10(RawCount36)

#Y axis is log10 raw counts
par(mfrow=c(3,3))
par(mar=c(3,4,3,3))
for(i in 1:length(my36)){
  
  plot(as.numeric(LogCount36[i,]),main=paste(rownames(LogCount36)[i],sep=""),ylab="log10 raw counts",xlab="Sample",ylim = c(0,5))
  abline(h=1,lty = "longdash",col = "red")
  
  threshpercent = round(length(which(LogCount36[i,] > 1))/length(LogCount36[i,]),3)*100
  
  if(mean(as.numeric(LogCount36[i,])) < 3.5){
    text(380,4.9,labels = paste(threshpercent,"% > 10 counts",sep=""))
  }else{
    text(75,0.5,labels = paste(threshpercent,"% > 10 counts",sep=""))
  }
  
}
par(mfrow=c(1,1))


#Y axis log scale
# par(mfrow=c(3,3))
# par(mar=c(3,4,3,3))
# for(i in 1:length(my36)){
#   
#   plot(as.numeric(RawCount36[i,]),main=paste(rownames(RawCount36)[i],sep=""),log = "y",ylab="Raw Counts (log scale)",xlab="Sample")
#   abline(h=10,lty = "longdash",col = "red")
#   
#   threshpercent = round(length(which(RawCount36[i,] > 10))/length(LogCount36[i,]),3)*100
#   
#   # if(mean(as.numeric(LogCount36[i,])) < 3.5){
#   #   text(380,4.9,labels = paste(threshpercent,"% > 10 counts",sep=""))
#   # }else{
#   #   text(75,0.5,labels = paste(threshpercent,"% > 10 counts",sep=""))
#   # }
#   
# }
# par(mfrow=c(1,1))


################## Raw Count Grouped Boxplot  ##########################

boxplotmat = LogCount36
boxplotmat = matrix(as.numeric(unlist(boxplotmat)),nrow(boxplotmat),ncol(boxplotmat))
rownames(boxplotmat) = rownames(LogCount36)
colnames(boxplotmat) = colnames(LogCount36)
typeof(boxplotmat)

#Without Legend (axis text far too small)
# decor = colorRampPalette(c("green","purple"))
# plotdecor = decor(nrow(LogCount36))
# boxplot.matrix(boxplotmat,use.cols=FALSE,xaxt = 'n',ylab="log10 raw counts",xlab="Subtype-specific feature",main="Raw Counts by feature",col=plotdecor)
# axis(1,1:36,labels = rownames(LogCount36),cex.axis=0.3)
# abline(h =1,col="red",lty="longdash")

#With Legend
decor = colorRampPalette(c("red","yellow","cyan","blue","purple","pink"))
plotdecor = decor(nrow(LogCount36))
par(mar=c(4.1, 4.1, 4.1, 2.1))
boxplot.matrix(boxplotmat,use.cols=FALSE,xaxt = 'n',ylab="log10 raw counts",main="Raw Counts by feature",col=plotdecor,cex.axis=1.5,cex.lab=1.5)
abline(h =1,col="red",lty="longdash")
legend(29.5,4.4,legend=rownames(LogCount36),col = plotdecor,pch=15,cex = 0.6,ncol=3)
title(xlab = "Subtype-specific feature",mgp=c(1,1,0),cex.lab = 1.5)

#ggplot2
library(ggplot2)

Features = rownames(boxplotmat)
tmp = stack(as.data.frame(boxplotmat))
gplotmat = data.frame(Features,tmp)

p = ggplot(gplotmat,aes(x=Features,y=boxplotmat,fill=Features)) + geom_violin(kernel = "gaussian",scale="width")
ttl = "Raw Counts by feature"
p = p +ggtitle(ttl) + xlab("Subtype-specific features") + ylab("log10 Raw Counts")
p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20),plot.title = element_text(size=24))
p = p+theme(axis.title.y=element_text(angle=90, vjust=2))
#p = p+geom_dotplot(binaxis = 'y',stackdir = 'center',dotsize = 0.1,fill="black",stackratio = 0.5)
p = p+theme(plot.title = element_text(hjust = 0.5))
p = p+geom_hline(yintercept = 1,linetype="dashed",col="red",size=1.2)
p = p+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
p


################## Mean Expression vs Log Fold Change ##########################

# rCountData_rinsite = round(FullCount_sr,0)
# dds_rinsite = DESeqDataSetFromMatrix(countData = rCountData_rinsite, colData = FullPheno_sr, design= ~ platform + Site + RIN + Subtype, tidy=F) #Subtype must be last (DESeq2 vignette)
# dds_rinsite$Subtype = relevel(dds_rinsite$Subtype,ref = "Control")
# dseq_rinsite = DESeq(dds_rinsite,betaPrior=T)


res = results(dseq_rinsite)


MAmatrix = plotMA(res,returnData = T)
rownames(MAmatrix) = rownames(FullCount)
FiltMAmatrix = MAmatrix[rownames(MAmatrix) %in% Transcripts,]

#Plot Filtered genes
plot(FiltMAmatrix$mean,FiltMAmatrix$lfc,xlim = c(0,6000),main = "Mean Expression vs LFC, Healthy Controls vs ALS",pch=1)
plot(log10(FiltMAmatrix$mean),FiltMAmatrix$lfc,main = "Mean Expression vs LFC, Healthy Controls vs ALS",pch=1)

#Give custom colors and point styles for the 36 genes selected in the paper
customcolors = rep("gray80",length(FiltMAmatrix$mean))
custompoints = rep(1,length(FiltMAmatrix$mean))

for(i in 1:length(FiltMAmatrix$mean)){
  
  for(j in 1:length(my36)){
    
    if(rownames(FiltMAmatrix)[i] == my36[j]){
      
      customcolors[i] = "red"
      custompoints[i] = 19
    }
    
    
  }
  
}

plot(log10(FiltMAmatrix$mean),FiltMAmatrix$lfc,main = "Mean Expression vs Log Fold Change",pch=custompoints,col = customcolors,xlab = "log10 mean of normalized counts", ylab = "LFC: ALS vs Control Donors")
abline(v = 1,lty = 'longdash',col="blue")
text(0.5,10,labels = "=10 normalized counts",col="blue")
threshp = round(length(which(log10(FiltMAmatrix$mean)>1))/length(FiltMAmatrix$mean),3)*100
text(2,-14,labels = paste(threshp,"% with mean expression > 10 normalized counts",sep = ""),col="blue")
legend(3.6,11.6,legend = c("Subtype-specific features (Fig. 6)","Enrichment Features"),col=c("red","gray50"),pch =c(19,1),cex = 0.9,pt.cex = 1)


meancounts = FiltMAmatrix$mean
my36mean = FiltMAmatrix$mean[rownames(FiltMAmatrix) %in% my36]

library(beeswarm)

par(mar=c(3,5,2,3))
boxplot(log10(meancounts),log10(my36mean),ylab = "log10 mean of normalized counts",xaxt="n",col=c("lightsteelblue1","lightsteelblue3"),cex.axis=1.6,cex.lab=1.5,outline=F,ylim=c(0,5))
axis(1,at=c(1,2),labels = c("Enrichment Features","Subtype-specific Features"),cex.axis=1.6)
abline(h=1,col="red",lty="longdash")
beeswarm(log10(meancounts),at=1,add = T,cex=0.75,pch=1,col="darkorange1")
beeswarm(log10(my36mean),at=2,add = T,cex=1,pch=19,col="darkorange2")


######################################################################################################################################################################################
#################################  QUESTION 3  #######################################################################################################################################
######################################################################################################################################################################################

#It is a common analysis to apply a deconvolution algorithm to estimate the proportions 
#of the major cell-types within a bulk RNA-seq sample. In human cortex samples one can estimate 
#the proportions of the major cell-types and compare proportions across groups of samples. 
#I am curious whether the cell-type proportions differ between the three subgroups, given your 
#findings of altered microglia and astrocyte marker genes in the ALS-Glia subtype.

#Deconv Stats
setwd("D:/Jarrett/Research/Nat Comm Reviews/Cell Deconvolution")
CIBERSORT_data = read.csv("CIBERSORTx_ALS5000_with_subtype.csv")

gind = which(CIBERSORT_data$Subtype == "GLIA")
oind = which(CIBERSORT_data$Subtype == "OX")
tind = which(CIBERSORT_data$Subtype == "TE")
hind = which(CIBERSORT_data$Subtype == "Control")
find = which(CIBERSORT_data$Subtype == "FTLD")

GliaDecon = CIBERSORT_data[gind,]
OxDecon = CIBERSORT_data[oind,]
TDDecon = CIBERSORT_data[tind,]
HCDecon = CIBERSORT_data[hind,]
FTDecon = CIBERSORT_data[find,]

subtypes = c("Glia","Ox","TD","HC","FTLD")
celltypes = colnames(GliaDecon)[-1:-2]
celltypes = celltypes[-11:-13]
celltypes

WilcoxRes = data.frame(matrix(NA,nrow= length(celltypes),ncol = 10))
rownames(WilcoxRes) = celltypes
colnames(WilcoxRes) = c("Glia_v_Ox","Glia_v_TD","Glia_v_HC","Glia_v_FTLD","Ox_v_TD","Ox_v_HC","Ox_v_FTLD","TD_v_HC","TD_v_FTLD","HC_v_FTLD")

quickfun = function(x,y){
  z = x[[y]]
  return(z)
}

for(i in 1:length(celltypes)){
  WilcoxRes[i,1] = wilcox.test(quickfun(GliaDecon,celltypes[i]),quickfun(OxDecon,celltypes[i]))$p.value
  WilcoxRes[i,2] = wilcox.test(quickfun(GliaDecon,celltypes[i]),quickfun(TDDecon,celltypes[i]))$p.value
  WilcoxRes[i,3] = wilcox.test(quickfun(GliaDecon,celltypes[i]),quickfun(HCDecon,celltypes[i]))$p.value
  WilcoxRes[i,4] = wilcox.test(quickfun(GliaDecon,celltypes[i]),quickfun(FTDecon,celltypes[i]))$p.value
  WilcoxRes[i,5] = wilcox.test(quickfun(OxDecon,celltypes[i]),quickfun(TDDecon,celltypes[i]))$p.value
  WilcoxRes[i,6] = wilcox.test(quickfun(OxDecon,celltypes[i]),quickfun(HCDecon,celltypes[i]))$p.value
  WilcoxRes[i,7] = wilcox.test(quickfun(OxDecon,celltypes[i]),quickfun(FTDecon,celltypes[i]))$p.value
  WilcoxRes[i,8] = wilcox.test(quickfun(TDDecon,celltypes[i]),quickfun(HCDecon,celltypes[i]))$p.value
  WilcoxRes[i,9] = wilcox.test(quickfun(TDDecon,celltypes[i]),quickfun(FTDecon,celltypes[i]))$p.value
  WilcoxRes[i,10] = wilcox.test(quickfun(HCDecon,celltypes[i]),quickfun(FTDecon,celltypes[i]))$p.value
}

AdjP = p.adjust(as.numeric(unlist(WilcoxRes)),method = "bonferroni")
AdjPmat = matrix(AdjP,nrow(WilcoxRes),ncol(WilcoxRes))

colnames(AdjPmat) = c("Glia_v_Ox","Glia_v_TD","Glia_v_HC","Glia_v_FTLD","Ox_v_TD","Ox_v_HC","Ox_v_FTLD","TD_v_HC","TD_v_FTLD","HC_v_FTLD")
rownames(AdjPmat) = celltypes

#write.csv(AdjPmat,"CIBERSORT_Stats_PeerReview.csv")

######################################################################################################################################################################################
#################################  QUESTION 4  #######################################################################################################################################
######################################################################################################################################################################################

#The increased sample size relative to Tam et al has allowed the authors to examine multiple clinical variables in association with their ALS subtype labels. 
#I am curious why the authors did not also include associations with the presence of C9orf72 repeat expansions or SOD1 mutations, the two most common genetic causes of ALS.

####################### Known genetic associations (C9orf72,SOD1) at the patient level
library(stringr)

setwd("C:/Users/jeshima/Desktop/ALS Manuscript/Nature/Communications Submission/Peer Review")
MetaData = read.csv("Table_S10.csv")

#Reference values
nsubtype = as.numeric(length(table(MetaData$FinalSubtype)))
nmutations = 2 #c9orf72 and sod1
C9_char = as.numeric(nchar("C9orf72"))
SOD_char = as.numeric(nchar("SOD1"))
MySubtypes = names(table(MetaData$FinalSubtype))
nGlia = as.numeric(length(which(MetaData$FinalSubtype == "GLIA")))
nOx = as.numeric(length(which(MetaData$FinalSubtype == "OX")))
nTD = as.numeric(length(which(MetaData$FinalSubtype == "TE")))
nDisc = as.numeric(length(which(MetaData$FinalSubtype == "Discordant")))
npatients = as.numeric(nrow(MetaData))

#use string searching method 

for(i in 1:nsubtype){
  
  currentsub = MySubtypes[i]
  
  #Build container
  tmpmat = data.frame(matrix(0,nmutations,nmutations+1))
  colnames(tmpmat) = c("C9orf72","SOD1","Unknown")
  rownames(tmpmat) = c("Positive","Negative")
  
  
  for(j in 1:npatients){
    
    if(MetaData$FinalSubtype[j] == MySubtypes[i]){
      
      muts = str_split(MetaData$Mutation[j],",") #Manually converted ; to , in excel
      nmuts = length(muts[[1]])
      
      for(k in 1:nmuts){
        
        charlen = nchar(muts[[1]][k])
        tmpstring = muts[[1]][k]
        
        if(muts[[1]][k] == "Unknown"){
          tmpmat[1,3] = tmpmat[1,3]+1
        }
        
        #C9orf72
        for(m in 1:charlen-C9_char+1){
          
          if(str_sub(tmpstring,m,m+C9_char-1) == "C9orf72"){
            
            for(p in 1:charlen){
              if(str_sub(tmpstring,p,p+2) == "Pos" || str_sub(tmpstring,p,p+2) =="pos"){ #Ignores case sensitivity
                tmpmat[1,1] = tmpmat[1,1]+1
              }else if(str_sub(tmpstring,p,p+2) == "Neg" || str_sub(tmpstring,p,p+2) =="neg"){
                tmpmat[2,1] = tmpmat[2,1]+1
              }
            }
            
          }
          
        }
        
        #SOD1
        for(n in 1:charlen-SOD_char+1){
          
          if(str_sub(tmpstring,n,n+SOD_char-1) == "SOD1"){
            
            for(q in 1:charlen){
              
              if(str_sub(tmpstring,q,q+2) == "Pos" || str_sub(tmpstring,q,q+2) == "pos"){
                tmpmat[1,2] = tmpmat[1,2]+1
              }else if(str_sub(tmpstring,q,q+2) == "Neg" || str_sub(tmpstring,q,q+2) =="neg"){
                tmpmat[2,2] = tmpmat[2,2]+1
              }
              
            }
            
          }
          
        }
        
      }
      
    }
    
    
  }
  #Return results
  tmpname = paste(MySubtypes[i],"Mutations",sep="")
  assign(tmpname,tmpmat)
}


#Plot

CleanLabels = c("Discordant","ALS-Glia","ALS-Ox","ALS-TD")

#Generate plotting data frames
plotdatpos = data.frame(matrix(NA,3,4))
colnames(plotdatpos) = MySubtypes
rownames(plotdatpos) = colnames(GLIAMutations)
plotdatpos[,1] = as.numeric(DiscordantMutations[1,])
plotdatpos[,2] = as.numeric(GLIAMutations[1,])
plotdatpos[,3] = as.numeric(OXMutations[1,])
plotdatpos[,4] = as.numeric(TEMutations[1,])

plotdatneg = data.frame(matrix(NA,2,4))
colnames(plotdatneg) = MySubtypes
rownames(plotdatneg) = colnames(GLIAMutations)[1:2]
plotdatneg[,1] = as.numeric(DiscordantMutations[2,1:2])
plotdatneg[,2] = as.numeric(GLIAMutations[2,1:2])
plotdatneg[,3] = as.numeric(OXMutations[2,1:2])
plotdatneg[,4] = as.numeric(TEMutations[2,1:2])

plotdatpos = as.matrix(plotdatpos)
colnames(plotdatpos) = CleanLabels
rownames(plotdatpos) = colnames(GLIAMutations)
plotdatneg = as.matrix(plotdatneg)
colnames(plotdatneg) = CleanLabels
rownames(plotdatneg) = colnames(GLIAMutations)[1:2]

#Positive Mutations
barplot(plotdatpos,beside=T,col = c("palegreen2","plum2","ivory3"),main = "Positive Mutations by Subtype",ylab = "Frequency",ylim=c(0,35))
legend(1,35,legend = c("C9orf72","SOD1","Unknown"),col = c("palegreen2","plum2","ivory3"),pch=15,pt.cex=1.5)
abline(h=0)

#Negative Mutations
barplot(plotdatneg,beside=T,col = c("palegreen2","plum2"),main = "Negative Mutations by Subtype",ylab = "Frequency",ylim=c(0,50))
legend(1,50,legend = c("C9orf72","SOD1"),col = c("palegreen2","plum2"),pch=15,pt.cex=1.5)
abline(h=0)

#Discordant Mutations
DiscordantMutations = as.matrix(DiscordantMutations)
customcol = c(rep(c("palegreen2","plum2"),nmutations),"gray50")
barplot(DiscordantMutations,beside=T,col = customcol,main = "Discordant Mutations",ylab = "Frequency",ylim=c(0,20),xaxt = "n")
axis(1,c(2,5,7.5),labels = c("C9orf72","SOD1","Unknown"))
legend(7,20,legend = c("Positive","Negative","Unknown"),col = c("palegreen2","plum2","gray50"),pch=15,pt.cex=1.5)
abline(h=0)


#GLIA Mutations
GLIAMutations = as.matrix(GLIAMutations)
customcol = c(rep(c("palegreen2","plum2"),nmutations),"gray50")
barplot(GLIAMutations,beside=T,col = customcol,main = "ALS-Glia Mutations",ylab = "Frequency",ylim=c(0,30),xaxt = "n")
axis(1,c(2,5,7.5),labels = c("C9orf72","SOD1","Unknown"))
legend(7,30,legend = c("Positive","Negative","Unknown"),col = c("palegreen2","plum2","gray50"),pch=15,pt.cex=1.5)
abline(h=0)


#OX Mutations
OXMutations = as.matrix(OXMutations)
customcol = c(rep(c("palegreen2","plum2"),nmutations),"gray50")
barplot(OXMutations,beside=T,col = customcol,main = "ALS-Ox Mutations",ylab = "Frequency",ylim=c(0,50),xaxt = "n")
axis(1,c(2,5,7.5),labels = c("C9orf72","SOD1","Unknown"))
legend(7,50,legend = c("Positive","Negative","Unknown"),col = c("palegreen2","plum2","gray50"),pch=15,pt.cex=1.5)
abline(h=0)


#TD Mutations
TEMutations = as.matrix(TEMutations)
customcol = c(rep(c("palegreen2","plum2"),nmutations),"gray50")
barplot(TEMutations,beside=T,col = customcol,main = "ALS-TD Mutations",ylab = "Frequency",ylim=c(0,30),xaxt = "n")
axis(1,c(2,5,7.5),labels = c("C9orf72","SOD1","Unknown"))
legend(1,30,legend = c("Positive","Negative","Unknown"),col = c("palegreen2","plum2","gray50"),pch=15,pt.cex=1.5)
abline(h=0)


#Manually clean the matrices
GLIAMutations2 = GLIAMutations
tmp = rep(NA,ncol(GLIAMutations))
GLIAMutations2 = rbind(GLIAMutations2,tmp)
rownames(GLIAMutations2) = c("Positive","Negative","Unknown")
GLIAMutations2[3,3] = GLIAMutations2[1,3]
GLIAMutations2[1,3] = NA
GLIAMutations2[2,3] = NA

DiscordantMutations2 = DiscordantMutations
tmp = rep(NA,ncol(DiscordantMutations))
DiscordantMutations2 = rbind(DiscordantMutations2,tmp)
rownames(DiscordantMutations2) = c("Positive","Negative","Unknown")
DiscordantMutations2[3,3] = DiscordantMutations2[1,3]
DiscordantMutations2[1,3] = NA
DiscordantMutations2[2,3] = NA

OXMutations2 = OXMutations
tmp = rep(NA,ncol(OXMutations))
OXMutations2 = rbind(OXMutations2,tmp)
rownames(OXMutations2) = c("Positive","Negative","Unknown")
OXMutations2[3,3] = OXMutations2[1,3]
OXMutations2[1,3] = NA
OXMutations2[2,3] = NA

TEMutations2 = TEMutations
tmp = rep(NA,ncol(TEMutations))
TEMutations2 = rbind(TEMutations2,tmp)
rownames(TEMutations2) = c("Positive","Negative","Unknown")
TEMutations2[3,3] = TEMutations2[1,3]
TEMutations2[1,3] = NA
TEMutations2[2,3] = NA



#Grouped and stacked bar chart

#Create data frame

FullBar = data.frame(matrix(NA,20,4))
colnames(FullBar) = c("Facet","Group","Mutation","Value")
FullBar$Facet = c(rep("Discordant",5),rep("GLIA",5),rep("OX",5),rep("TE",5))
FullBar$Group = rep(c(rep(c("Positive","Negative"),2),"Unknown"),nsubtype)
FullBar$Mutation = rep(c(rep("C9orf72",2),rep("SOD1",2),"Unknown"),nsubtype)

#Fill data frame
for(i in 1:nsubtype){
  
  tmpdat = get(paste(MySubtypes[i],"Mutations2",sep=""))
  
  for(j in 1:nrow(tmpdat)){
    for(k in 1:ncol(tmpdat)){
      for(l in 1:nrow(FullBar)){
        if(FullBar$Facet[l] == MySubtypes[i]){
          if(rownames(tmpdat)[j] == FullBar$Group[l] && colnames(tmpdat)[k] == FullBar$Mutation[l]){
            FullBar$Value[l] = tmpdat[j,k]
          }
          
        }
      }
      
    }
    
  }
  
}

FullBar$Facet = c(rep("Discordant",5),rep("ALS-Glia",5),rep("ALS-Ox",5),rep("ALS-TD",5))

#mycol = rep(c(rep("lightsalmon2",2),rep("palegreen3",2),"thistle3"),nsubtype)

p = ggplot(FullBar,aes(x=Group,y=Value,fill=Mutation)) + geom_bar(stat = "identity",position = "stack") + facet_grid(~Facet)
p = p+scale_fill_manual(values = c("lightsalmon2","palegreen3","thistle3"))
p = p +ggtitle("C9orf72 and SOD1 Mutations by Subtype") + xlab("") + ylab("Frequency")
p = p+theme(axis.text = element_text(size=12), axis.title = element_text(size=14),plot.title = element_text(size=22))
p = p+theme(axis.title.y=element_text(angle=90, vjust=2,size = 18))
p = p+theme(plot.title = element_text(hjust = 0.5))
p = p+theme(text = element_text(size=14))
p = p+theme(plot.margin = margin(t=10,r=10,b=10,l=10))
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "gray80"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray80"))
p = p+theme(axis.text.x = element_text(size = 13))
p = p+theme(axis.text.y = element_text(size = 16))
p = p+theme(legend.text = element_text(size = 16))
p = p+theme(legend.title = element_text(size = 16))
p


#################chi-squared stats

#Contingency Tables (mutation-wise)

CSmat = data.frame(matrix(NA,nsubtype,nmutations+2))
colnames(CSmat) = c("Subtype","Positive_C9orf72","Negative_C9orf72","Unknown")
CSmat$Subtype = c("Discordant","ALS-Glia","ALS-Ox","ALS-TD")

CSmatSOD = data.frame(matrix(NA,nsubtype,nmutations+2))
colnames(CSmatSOD) = c("Subtype","Positive_SOD1","Negative_SOD1","Unknown")
CSmatSOD$Subtype = c("Discordant","ALS-Glia","ALS-Ox","ALS-TD")

for(i in 1:nrow(CSmat)){
  
  for(j in 1:nrow(FullBar)){
    
    if(CSmat$Subtype[i] == FullBar$Facet[j] && FullBar$Group[j] == "Positive" && FullBar$Mutation[j] == "C9orf72"){
      CSmat$Positive_C9orf72[i] = FullBar$Value[j]
    }else if(CSmat$Subtype[i] == FullBar$Facet[j] && FullBar$Group[j] == "Negative" && FullBar$Mutation[j] == "C9orf72"){
      CSmat$Negative_C9orf72[i] = FullBar$Value[j]
    }else if(CSmat$Subtype[i] == FullBar$Facet[j] && FullBar$Group[j] == "Unknown" && FullBar$Mutation[j] == "Unknown"){
      CSmat$Unknown[i] = FullBar$Value[j]
    }
    
    
    if(CSmatSOD$Subtype[i] == FullBar$Facet[j] && FullBar$Group[j] == "Positive" && FullBar$Mutation[j] == "SOD1"){
      CSmatSOD$Positive_SOD1[i] = FullBar$Value[j]
    }else if(CSmatSOD$Subtype[i] == FullBar$Facet[j] && FullBar$Group[j] == "Negative" && FullBar$Mutation[j] == "SOD1"){
      CSmatSOD$Negative_SOD1[i] = FullBar$Value[j]
    }else if(CSmatSOD$Subtype[i] == FullBar$Facet[j] && FullBar$Group[j] == "Unknown" && FullBar$Mutation[j] == "Unknown"){
      CSmatSOD$Unknown[i] = FullBar$Value[j]
    }
    
  }
  
}

#Clean up
rownames(CSmat) = CSmat[,1]
CSmat = CSmat[,-1]
rownames(CSmatSOD) = CSmatSOD[,1]
CSmatSOD = CSmatSOD[,-1]

#With unknown category
chisq.test(CSmat) #p = 0.006119 (significance driven by unknown category, see results below)
chisq.test(CSmatSOD) #p = 0.483

#Without unknown category
CSmat2 = CSmat[-which(colnames(CSmat) == "Unknown")]
CSmatSOD2 = CSmatSOD[-which(colnames(CSmatSOD) == "Unknown")]

chisq.test(CSmat2) #p = 0.4726
chisq.test(CSmatSOD2) #p = 0.2148

######################################################################################################################################################################################
#################################  QUESTION 5  #######################################################################################################################################
######################################################################################################################################################################################

#Furthermore, the expression of the TDP-43 regulated cryptic splicing event in STMN2 has already been quantified in these samples (Prudencio et al, 2020). 
#Given the transcriptional dysregulation subtype, I am curious why the authors did not look at the presence of truncated STMN2 splicing in their subtypes?

load("D:/Jarrett/Research/Nat Comm Reviews/RData/ALSPatientStratification_UnivariateDatasets_RINSite_PeerReview.RData")


tmp = table(ALSPheno$Subtype)
tmp

Glia.STMN2 = rep(NA,tmp[[1]])
Ox.STMN2 = rep(NA,tmp[[2]])
TD.STMN2 = rep(NA,tmp[[3]])

Glia.tSTMN2 = rep(NA,tmp[[1]])
Ox.tSTMN2 = rep(NA,tmp[[2]])
TD.tSTMN2 = rep(NA,tmp[[3]])

Glia.tSTMN2.TPM = rep(NA,tmp[[1]])
Ox.tSTMN2.TPM = rep(NA,tmp[[2]])
TD.tSTMN2.TPM = rep(NA,tmp[[3]])

count1 = count2 = count3 = 1
for(i in 1:nrow(ALSPheno)){
  
  if(ALSPheno$Subtype[i] == "GLIA"){
    
    Glia.STMN2[count1] = ALSPheno$STMN2_TPM[i]
    Glia.tSTMN2[count1] = ALSPheno$tSTMN2_counts[i]
    Glia.tSTMN2.TPM[count1] = ALSPheno$tSTMN2_TPM[i]
    count1 = count1 + 1
    
  }else if(ALSPheno$Subtype[i] == "OX"){
    
    Ox.STMN2[count2] = ALSPheno$STMN2_TPM[i]
    Ox.tSTMN2[count2] = ALSPheno$tSTMN2_counts[i]
    Ox.tSTMN2.TPM[count2] = ALSPheno$tSTMN2_TPM[i]
    count2 = count2 + 1
    
  }else if(ALSPheno$Subtype[i] == "TE"){
    
    TD.STMN2[count3] = ALSPheno$STMN2_TPM[i]
    TD.tSTMN2[count3] = ALSPheno$tSTMN2_counts[i]
    TD.tSTMN2.TPM[count3] = ALSPheno$tSTMN2_TPM[i]
    count3 = count3 + 1
    
  }
  
}



####################### Wilcoxon Rank Sum Approach


#Individual subjects are not normally distributed... use nonparametric
hist(Glia.STMN2)
abline(v = mean(Glia.STMN2),col = "red")

hist(Ox.STMN2)
abline(v = mean(Ox.STMN2),col = "red")

hist(TD.STMN2)
abline(v = mean(TD.STMN2),col = "red")

hist(c(Glia.STMN2,Ox.STMN2,TD.STMN2),breaks = 30)


#Pairwise Wilcoxon
#Remove Glia 
tmpi = which(ALSPheno$Subtype == "GLIA")
tmpPheno = ALSPheno[-tmpi,]
OxTD_ranksum = wilcox.test(STMN2_TPM ~ Factor, data = tmpPheno)
OxTD_pval = OxTD_ranksum$p.value

#Remove Ox
tmpi = which(ALSPheno$Subtype == "OX")
tmpPheno = ALSPheno[-tmpi,]
GliaTD_ranksum = wilcox.test(STMN2_TPM ~ Factor, data = tmpPheno)
GliaTD_pval = GliaTD_ranksum$p.value

#Remove TD
tmpi = which(ALSPheno$Subtype == "TE")
tmpPheno = ALSPheno[-tmpi,]
GliaOx_ranksum = wilcox.test(STMN2_TPM ~ Factor, data = tmpPheno)
GliaOx_pval = GliaOx_ranksum$p.value

Adjustedp = p.adjust(c(OxTD_pval,GliaTD_pval,GliaOx_pval), method = "bonferroni")
OxTD_pval_adj = Adjustedp[1]
GliaTD_pval_adj = Adjustedp[2]
GliaOx_pval_adj = Adjustedp[3]

barplot(c(mean(Glia.STMN2),mean(Ox.STMN2),mean(TD.STMN2)),main = "STMN2 Counts by ALS Subtype",ylab = "Mean STMN2 Counts (TPM)",col=c("goldenrod","navy","firebrick"),ylim=c(0,140))
axis(1,c(0.65,1.9,3.1),labels = c("ALS-Glia","ALS-Ox","ALS-TD"))
lines(c(0.65,1.9),c(100,100),type = "l")
lines(c(0.65,0.65),c(95,100),type="l")
lines(c(1.9,1.9),c(95,100),type="l")
text(mean(c(0.65,1.9)),105,labels = paste("Bonferroni-adjusted p-value:",formatC(GliaOx_pval_adj,format = "e",digits = 2)))
lines(c(1.9,3.1),c(110,110),type = "l")
lines(c(1.9,1.9),c(105,110),type="l")
lines(c(3.1,3.1),c(105,110),type="l")
text(mean(c(1.9,3.1)),115,labels = paste("Bonferroni-adjusted p-value:",formatC(OxTD_pval_adj,format = "e",digits = 2)))
lines(c(0.65,3.1),c(125,125),type="l")
lines(c(0.65,0.65),c(120,125),type="l")
lines(c(3.1,3.1),c(120,125),type="l")
text(mean(c(0.65,3.1)),130,labels = paste("Bonferroni-adjusted p-value:",formatC(GliaTD_pval_adj,format = "e",digits = 2)))


par(mar=c(5.1, 5.1, 4.1, 2.1))
boxplot(Glia.STMN2,Ox.STMN2,TD.STMN2,main = "STMN2 Counts by ALS Subtype",ylab = "STMN2 Counts (TPM)",xlab = "Subtype",col=c("goldenrod","navy","firebrick"),ylim=c(0,500),cex.axis = 1.5,cex.lab=1.75,cex.main=1.75)
axis(1,c(1,2,3),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=1.75)
lines(c(1,2),c(410,410),type = "l")
lines(c(1,1),c(400,410),type="l")
lines(c(2,2),c(400,410),type="l")
text(mean(c(1,2)),420,labels = paste("Bonferroni-adjusted p-value:",formatC(GliaOx_pval_adj,format = "e",digits = 2)),cex = 1.35)
lines(c(2,3),c(440,440),type = "l")
lines(c(2,2),c(430,440),type="l")
lines(c(3,3),c(430,440),type="l")
text(mean(c(2,3)),450,labels = paste("Bonferroni-adjusted p-value:",formatC(OxTD_pval_adj,format = "e",digits = 2)),cex = 1.35)
lines(c(1,3),c(470,470),type="l")
lines(c(1,1),c(460,470),type="l")
lines(c(3,3),c(460,470),type="l")
text(mean(c(0.65,3.1)),480,labels = paste("Bonferroni-adjusted p-value:",formatC(GliaTD_pval_adj,format = "e",digits = 2)),cex = 1.35)

setwd("D:/Jarrett/Research/Nat Comm Reviews/Source")
write.csv(Glia.STMN2,"FigS15A-Glia.csv")
write.csv(Ox.STMN2,"FigS15A-Ox.csv")
write.csv(TD.STMN2,"FigS15A-TD.csv")


###################### Repeat for truncated STMN2 counts (TPM)

#Individual subjects are not normally distributed... use nonparametric
hist(Glia.tSTMN2.TPM)
abline(v = mean(Glia.tSTMN2.TPM),col = "red")

hist(Ox.tSTMN2.TPM)
abline(v = mean(Ox.tSTMN2.TPM),col = "red")

hist(TD.tSTMN2.TPM)
abline(v = mean(TD.tSTMN2.TPM),col = "red")

hist(c(Glia.tSTMN2.TPM,Ox.tSTMN2.TPM,TD.tSTMN2.TPM),breaks = 30)

par(mfrow=c(1,3))
hist(Glia.tSTMN2.TPM,main = "ALS-Glia tSTMN2 (TPM)",xlab = "tSTMN2 (TPM)",ylim=c(0,70),xlim = c(0,1),breaks = 10,cex.main=2,cex.lab=1.5,cex.axis=1.5,col="goldenrod");hist(Ox.tSTMN2.TPM,main = "ALS-Ox tSTMN2 (TPM)",xlab = "tSTMN2 (TPM)",ylim=c(0,200),xlim = c(0,1),breaks = 10,cex.main=2,cex.lab=1.5,cex.axis=1.5,col="navy");hist(TD.tSTMN2.TPM,main = "ALS-TD tSTMN2 (TPM)",xlab = "tSTMN2 (TPM)",ylim=c(0,120),xlim = c(0,1),breaks = 10,cex.main=2,cex.lab=1.5,cex.axis=1.5,col="firebrick")
hist(Glia.tSTMN2,main = "ALS-Glia tSTMN2",xlab = "tSTMN2 (raw counts)",ylim=c(0,70),xlim = c(0,25),breaks = 10,cex.main=2,cex.lab=1.5,cex.axis=1.5,col="goldenrod");hist(Ox.tSTMN2,main = "ALS-Ox tSTMN2",xlab = "tSTMN2 (raw counts)",ylim=c(0,200),xlim = c(0,25),breaks = 10,cex.main=2,cex.lab=1.5,cex.axis=1.5,col="navy");hist(TD.tSTMN2,main = "ALS-TD tSTMN2",xlab = "tSTMN2 (raw counts)",ylim=c(0,120),xlim = c(0,25),breaks = 10,cex.main=2,cex.lab=1.5,cex.axis=1.5,col="firebrick")
par(mfrow=c(1,1))

#Pairwise Wilcoxon
#Remove Glia 
tmpi = which(ALSPheno$Subtype == "GLIA")
tmpPheno = ALSPheno[-tmpi,]
OxTD_ranksum = wilcox.test(tSTMN2_TPM ~ Factor, data = tmpPheno)
OxTD_pval = OxTD_ranksum$p.value

#Remove Ox
tmpi = which(ALSPheno$Subtype == "OX")
tmpPheno = ALSPheno[-tmpi,]
GliaTD_ranksum = wilcox.test(tSTMN2_TPM ~ Factor, data = tmpPheno)
GliaTD_pval = GliaTD_ranksum$p.value

#Remove TD
tmpi = which(ALSPheno$Subtype == "TE")
tmpPheno = ALSPheno[-tmpi,]
GliaOx_ranksum = wilcox.test(tSTMN2_TPM ~ Factor, data = tmpPheno)
GliaOx_pval = GliaOx_ranksum$p.value

Adjustedp = p.adjust(c(OxTD_pval,GliaTD_pval,GliaOx_pval), method = "bonferroni")
OxTD_pval_adj = Adjustedp[1]
GliaTD_pval_adj = Adjustedp[2]
GliaOx_pval_adj = Adjustedp[3]


boxplot(Glia.tSTMN2.TPM,Ox.tSTMN2.TPM,TD.tSTMN2.TPM,main = "Truncated STMN2 by ALS Subtype",ylab = "tSTMN2 Counts (TPM)",xlab = "Subtype",col=c("goldenrod","navy","firebrick"),ylim=c(0,1.1),cex.axis = 1.5,cex.lab=1.75,cex.main=1.75)
axis(1,c(1,2,3),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=1.75)
lines(c(1,2),c(0.8,0.8),type = "l")
lines(c(1,1),c(0.75,0.8),type="l")
lines(c(2,2),c(0.75,0.8),type="l")
text(mean(c(1,2)),0.83,labels = "Bonferroni-adjusted p-value > 0.05",cex = 1.35)
lines(c(2,3),c(0.9,0.9),type = "l")
lines(c(2,2),c(0.85,0.9),type="l")
lines(c(3,3),c(0.85,0.9),type="l")
text(mean(c(2,3)),0.93,labels = "Bonferroni-adjusted p-value > 0.05",cex = 1.35)
lines(c(1,3),c(1,1),type="l")
lines(c(1,1),c(0.95,1),type="l")
lines(c(3,3),c(0.95,1),type="l")
text(mean(c(1,3)),1.03,labels = "Bonferroni-adjusted p-value > 0.05",cex = 1.35)


setwd("D:/Jarrett/Research/Nat Comm Reviews/Source")
write.csv(Glia.tSTMN2.TPM,"FigS15B-Glia.csv")
write.csv(Ox.tSTMN2.TPM,"FigS15B-Ox.csv")
write.csv(TD.tSTMN2.TPM,"FigS15B-TD.csv")

###################### Repeat for truncated STMN2 counts (raw counts)

hist(Glia.tSTMN2)
abline(v = mean(Glia.tSTMN2),col = "red")

hist(Ox.tSTMN2)
abline(v = mean(Ox.tSTMN2),col = "red")

hist(TD.tSTMN2)
abline(v = mean(TD.tSTMN2),col = "red")

hist(c(Glia.tSTMN2,Ox.tSTMN2,TD.tSTMN2),breaks = 30)

par(mfrow=c(1,3))
hist(Glia.tSTMN2,main = "ALS-Glia tSTMN2",xlab = "tSTMN2 (counts)",xlim = c(0,25),breaks = 10);hist(Ox.tSTMN2,main = "ALS-Ox tSTMN2",xlab = "tSTMN2 (counts)",xlim = c(0,25),breaks = 10);hist(TD.tSTMN2,main = "ALS-TD tSTMN2",xlab = "tSTMN2 (counts)",xlim = c(0,25),breaks = 10)
par(mfrow=c(1,1))

#Pairwise Wilcoxon
#Remove Glia 
tmpi = which(ALSPheno$Subtype == "GLIA")
tmpPheno = ALSPheno[-tmpi,]
OxTD_ranksum = wilcox.test(tSTMN2_counts ~ Factor, data = tmpPheno)
OxTD_pval = OxTD_ranksum$p.value

#Remove Ox
tmpi = which(ALSPheno$Subtype == "OX")
tmpPheno = ALSPheno[-tmpi,]
GliaTD_ranksum = wilcox.test(tSTMN2_counts ~ Factor, data = tmpPheno)
GliaTD_pval = GliaTD_ranksum$p.value

#Remove TD
tmpi = which(ALSPheno$Subtype == "TE")
tmpPheno = ALSPheno[-tmpi,]
GliaOx_ranksum = wilcox.test(tSTMN2_counts ~ Factor, data = tmpPheno)
GliaOx_pval = GliaOx_ranksum$p.value

Adjustedp = p.adjust(c(OxTD_pval,GliaTD_pval,GliaOx_pval), method = "bonferroni")
OxTD_pval_adj = Adjustedp[1]
GliaTD_pval_adj = Adjustedp[2]
GliaOx_pval_adj = Adjustedp[3]


boxplot(Glia.tSTMN2,Ox.tSTMN2,TD.tSTMN2,main = "Truncated STMN2 by ALS Subtype",ylab = "tSTMN2 Counts",xlab = "Subtype",col=c("goldenrod","navy","firebrick"),ylim=c(0,35),cex.axis = 1.5,cex.lab=1.75,cex.main=1.75)
axis(1,c(1,2,3),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=1.75)
lines(c(1,2),c(27,27),type = "l")
lines(c(1,1),c(26,27),type="l")
lines(c(2,2),c(26,27),type="l")
text(mean(c(1,2)),28,labels = "Bonferroni-adjusted p-value > 0.05",cex=1.35)
lines(c(2,3),c(30,30),type = "l")
lines(c(2,2),c(29,30),type="l")
lines(c(3,3),c(29,30),type="l")
text(mean(c(2,3)),31,labels = "Bonferroni-adjusted p-value > 0.05",cex=1.35)
lines(c(1,3),c(33,33),type="l")
lines(c(1,1),c(32,33),type="l")
lines(c(3,3),c(32,33),type="l")
text(mean(c(1,3)),34,labels = "Bonferroni-adjusted p-value > 0.05",cex=1.35)

setwd("D:/Jarrett/Research/Nat Comm Reviews/Source")
write.csv(Glia.tSTMN2,"FigS15C-Glia.csv")
write.csv(Ox.tSTMN2,"FigS15C-Ox.csv")
write.csv(TD.tSTMN2,"FigS15C-TD.csv")

################################## Differential Expression Approach

#Add STMN2 to feature list

prog_gene_positions = gene_positions[gene_positions$hgnc_symbol %in% "STMN2",]

blank = rep(NA,nrow(FullCount))
FullCount_STMN2 = cbind(blank,FullCount)
FullCount_STMN2[,1] = rownames(FullCount)

DONE = F
for(i in 1:nrow(FullCount_STMN2)){
  for(j in 1:nrow(prog_gene_positions)){
    if(DONE == F){
      if(sub("\\..*","",FullCount_STMN2[i,1]) == prog_gene_positions$ensembl_gene_id[j]){
        FullCount_STMN2[i,1] = prog_gene_positions$hgnc_symbol[j]
        DONE = T
      }
    }
  }
  DONE = F
  if((i %% 500) == 0) cat("% Done:",i/nrow(FullCount_STMN2)*100,"\n")
}

#Replace missing gene symbol with original ENSEMBL ID
for(i in 1:nrow(FullCount_STMN2)){
  if(FullCount_STMN2[i,1] == ""){
    FullCount_STMN2[i,1] = rownames(FullCount_STMN2)[i]
  }
}

which(FullCount_STMN2$blank == 'STMN2')

rownames(FullCount_STMN2) = FullCount_STMN2$blank
FullCount_STMN2 = FullCount_STMN2[,-1]


#Add C9orf72
tmpref = gene_positions[gene_positions$hgnc_symbol %in% "C9orf72",]

c9ind = which(rownames(FullCount_STMN2) == tmpref$ensembl_gene_id)
rownames(FullCount_STMN2)[c9ind] = "C9orf72"

#check the values
# stmnind = which(rownames(FullCount_STMN2) == "STMN2")
# table(colnames(FullCount_STMN2) == FullPheno$Subject)
# countcheck1 = FullCount_STMN2[stmnind,]
# countcheck2 = FullPheno$STMN2_TPM
# View(rbind(countcheck1,countcheck2))


library(DESeq2)

Transcripts = c(Transcripts,"STMN2","C9orf72")

rCountData_STMN2 = round(FullCount_STMN2,0)
setwd("C:/Users/jeshima/Desktop/ALS Manuscript/Nature/Communications Submission/Peer Review")
FullPheno_sr = read.csv("Site_RIN_FullPheno_PeerReview.csv") #This is Table S13, sheet 2, with the missing RIN sample removed

dds_STMN2 = DESeqDataSetFromMatrix(countData = rCountData_STMN2, colData = FullPheno_sr, design= ~ platform + RIN + Site + Subtype, tidy=F) #Subtype must be second (DESeq2 vignette)
dds_STMN2$SubjectSex = relevel(dds_STMN2$Subtype,ref = "Control")
dseq_STMN2 = DESeq(dds_STMN2,betaPrior=T)

#Pairwise "contrast()"
glia.res = results(dseq_STMN2,contrast = c("Subtype","GLIA","Control"))
glia.sig = glia.res[! is.na(glia.res$padj) & glia.res$padj<0.05,]
filt.glia.sig = glia.sig[rownames(glia.sig) %in% Transcripts,]

ox.res = results(dseq_STMN2,contrast = c("Subtype","OX","Control"))
ox.sig = ox.res[! is.na(ox.res$padj) & ox.res$padj<0.05,]
filt.ox.sig = ox.sig[rownames(ox.sig) %in% Transcripts,]

TE.res = results(dseq_STMN2,contrast = c("Subtype","TE","Control"))
TE.sig = TE.res[! is.na(TE.res$padj) & TE.res$padj<0.05,]
filt.TE.sig = TE.sig[rownames(TE.sig) %in% Transcripts,]

GT.res = results(dseq_STMN2,contrast = c("Subtype","GLIA","TE"))
GT.sig = GT.res[! is.na(GT.res$padj) & GT.res$padj<0.05,]
filt.GT.sig = GT.res[rownames(GT.res) %in% Transcripts,]

GO.res = results(dseq_STMN2,contrast = c("Subtype","GLIA","OX"))
GO.sig = GO.res[! is.na(GO.res$padj) & GO.res$padj<0.05,]
filt.GO.sig = GO.sig[rownames(GO.sig) %in% Transcripts,]

TO.res = results(dseq_STMN2,contrast = c("Subtype","TE","OX"))
TO.sig = TO.res[! is.na(TO.res$padj) & TO.res$padj<0.05,]
filt.TO.sig = TO.sig[rownames(TO.sig) %in% Transcripts,]

glia.res.ond = results(dseq_STMN2,contrast = c("Subtype","GLIA","OND"))
glia.sig.ond = glia.res.ond[! is.na(glia.res.ond$padj) & glia.res.ond$padj<0.05,]
filt.glia.sig.ond = glia.sig.ond[rownames(glia.sig.ond) %in% Transcripts,]

ox.res.ond = results(dseq_STMN2,contrast = c("Subtype","OX","OND"))
ox.sig.ond = ox.res.ond[! is.na(ox.res.ond$padj) & ox.res.ond$padj<0.05,]
filt.ox.sig.ond = ox.sig.ond[rownames(ox.sig.ond) %in% Transcripts,]

TE.res.ond = results(dseq_STMN2,contrast = c("Subtype","TE","OND"))
TE.sig.ond = TE.res.ond[! is.na(TE.res.ond$padj) & TE.res.ond$padj<0.05,]
filt.TE.sig.ond = TE.sig.ond[rownames(TE.sig.ond) %in% Transcripts,]

COND.res = results(dseq_STMN2,contrast = c("Subtype","Control","OND"))
COND.sig = COND.res[! is.na(COND.res$padj) & COND.res$padj<0.05,]
filt.COND.sig = COND.sig[rownames(COND.sig) %in% Transcripts,]


tmp = estimateSizeFactors(dseq_STMN2)
DESeq_NormalizedCounts_60k_STMN2 = counts(tmp,normalized = T)
DESeq_NormalizedCounts_60k_filtered_STMN2 = DESeq_NormalizedCounts_60k_STMN2[rownames(DESeq_NormalizedCounts_60k_STMN2) %in% Transcripts,]
NormCounts_STMN2 = DESeq_NormalizedCounts_60k_filtered_STMN2

#For plotting purposes (not used during statistical analysis), adjust zero count genes to 1

for(i in 1:nrow(NormCounts_STMN2)){
  for(j in 1:ncol(NormCounts_STMN2)){
    if(NormCounts_STMN2[i,j] == 0){
      NormCounts_STMN2[i,j] = 1
    }
  }
}

##Parse for plotting

#Subtype Normalized Count Matrix
GliaI = which(FullPheno$Subtype == "GLIA")
GliaSamp = colnames(NormCounts_STMN2)[GliaI]
OxI = which(FullPheno$Subtype == "OX")
OxSamp = colnames(NormCounts_STMN2)[OxI]
TEI = which(FullPheno$Subtype == "TE")
TESamp = colnames(NormCounts_STMN2)[TEI]
ONDI = which(FullPheno$Subtype == "OND")
ONDSamp = colnames(NormCounts_STMN2)[ONDI]
HCI = which(FullPheno$Subtype == "Control")
HCSamp = colnames(NormCounts_STMN2)[HCI]

Subtype = FullPheno$Subtype

for(i in 1:length(Subtype)){
  if(Subtype[i] == "GLIA"){
    Subtype[i] = "ALS-Glia"
  }else if(Subtype[i] == "OX"){
    Subtype[i] = "ALS-Ox"
  }else if(Subtype[i] == "TE"){
    Subtype[i] = "ALS-TD"
  }else if(Subtype[i] == "OND"){
    Subtype[i] = "FTLD"
  }
}

Subtype = factor(Subtype,levels = c("Control","FTLD","ALS-Glia","ALS-Ox","ALS-TD"))

#"AutoPlot" is the simple function version of lines 471-2643 in ALSPatientStratification_UnivariateAnalysis.R
AutoPlot(PlotGene = "STMN2",Focus = "Ox",NormCounts = NormCounts_STMN2,Subtype = Subtype,filt.glia.sig = filt.glia.sig,filt.ox.sig = filt.ox.sig,filt.TE.sig = filt.TE.sig,filt.GT.sig = filt.GT.sig,filt.GO.sig = filt.GO.sig,filt.TO.sig = filt.TO.sig,filt.glia.sig.ond = filt.glia.sig.ond,filt.ox.sig.ond = filt.ox.sig.ond,filt.TE.sig.ond = filt.TE.sig.ond,filt.COND.sig = filt.COND.sig)
AutoPlot(PlotGene = "STMN2",Focus = "Glia",NormCounts = NormCounts_STMN2,Subtype = Subtype,filt.glia.sig = filt.glia.sig,filt.ox.sig = filt.ox.sig,filt.TE.sig = filt.TE.sig,filt.GT.sig = filt.GT.sig,filt.GO.sig = filt.GO.sig,filt.TO.sig = filt.TO.sig,filt.glia.sig.ond = filt.glia.sig.ond,filt.ox.sig.ond = filt.ox.sig.ond,filt.TE.sig.ond = filt.TE.sig.ond,filt.COND.sig = filt.COND.sig)

######################################################################################################################################################################################
#################################  Recommendation 1  #################################################################################################################################
######################################################################################################################################################################################

#This reviewer does not find dendrograms (3C) and network visualisations (3D-E) helpful. I would strongly suggest the authors rethink how they present their WGCNA network results.

library(WGCNA)

setwd("C:/Users/jeshima/Desktop/ALS Manuscript/Nature/Communications Submission/Supplemental")

#Read in expression data
datExpr = read.csv("Table_S5.csv")
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
     xlab = "", sub = "",hang=-1)



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

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;


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
min(moduleTraitCor)
max(moduleTraitCor)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = c("Age of Onset","Age of Death","Disease Duration"),
               yLabels = names(MEs),
               ySymbols = newMEs,
               colorLabels = FALSE,
               colors = rev(blueWhiteRed(50)),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.lab = 1.1,
               cex.text = 1.3,
               zlim = c(-0.2,0.2))



corr = data.frame(moduleTraitCor)
WGCNAp =  signif(moduleTraitPvalue, 1)
WGCNAp2 = - log10(WGCNAp)

library(ggplot2)

customheat = data.frame(matrix(NA,nrow(corr)*ncol(corr),3))
colnames(customheat) = c("Module","Response","Value")
plotorder = c(rep("MEpink",3),rep("MEred",3),rep("MEtan",3),rep("MEturquoise",3),rep("MEbrown",3),rep("MEgreen",3),rep("MEpurple",3),rep("MEgrey",3),rep("MEmagenta",3),rep("MEyellow",3),rep("MEblue",3),rep("MEsalmon",3),rep("MEblack",3),rep("MEgreenyellow",3))
customheat$Module = plotorder
customheat$Response = rep(colnames(corr),nrow(corr))

#Fill plot matrix
for(i in 1:nrow(customheat)){
  
  for(j in 1:nrow(corr)){
    
    for(k in 1:ncol(corr)){
      
      if(customheat$Module[i] == rownames(corr)[j] && customheat$Response[i] == colnames(WGCNAmat)[k]){
        
        customheat$Value[i] = corr[j,k]
        
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
EigGO = read.csv("EigengeneEnrichmentMatrix.csv") #Manually cleaned Bonferroni-adjusted p-values from Table S9
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

#summary(lm(MEs[,7] ~ Pheno$Subtype))

Subtypepval = data.frame(matrix(NA,nrow = length(names(table(Pheno$Subtype))),ncol(MEs)))
colnames(Subtypepval) = colnames(MEs)
rownames(Subtypepval) = c("ALS-Glia","ALS-Ox","ALS-TD")
Subtypesign = Subtypepval
for(i in 1:ncol(MEs)){
  model = lm(MEs[,i] ~ Pheno$Subtype)
  tmp = summary(model)$coefficients
  tmp2 = matrix(unlist(tmp),ncol=4)
  
  Subtypepval[1:3,i] = tmp2[,4] #p-val
  
  for(j in 1:nrow(tmp2)){
    
    if(tmp2[j,1] < 0){
      Subtypesign[j,i] = "Down"
    }else if(tmp2[j,1]>0){
      Subtypesign[j,i] = "Up"
    }
    
  }
}

Subtypepval = t(Subtypepval)
Subtypesign = t(Subtypesign)


plotordersub = c(rep("MEpink",6),rep("MEred",6),rep("MEtan",6),rep("MEturquoise",6),rep("MEbrown",6),rep("MEgreen",6),rep("MEpurple",6),rep("MEgrey",6),rep("MEmagenta",6),rep("MEyellow",6),rep("MEblue",6),rep("MEsalmon",6),rep("MEblack",6),rep("MEgreenyellow",6))
colorder = c("ALS-Glia Upregulated","ALS-Glia Downregulated","ALS-Ox Upregulated","ALS-Ox Downregulated","ALS-TD Upregulated","ALS-TD Downregulated")

customheatsub = data.frame(matrix(NA,nrow(logSubtypepval)*ncol(logSubtypepval)*2,3))
colnames(customheatsub) = c("Module","Response","Value")
customheatsub$Module = plotordersub
customheatsub$Response = rep(colorder,nrow(logSubtypepval))

evenseq = seq(2,nrow(customheatsub),6)
oddseq = seq(1,nrow(customheatsub),6)

neworder = NA
count=1
for(i in 1:nrow(Subtypesign)){
  
  for(j in 1:ncol(Subtypesign)){
    
    if(colnames(Subtypesign)[j] == "ALS-Glia" && Subtypesign[i,j] == "Up"){
      neworder[count]=1+(6*(i-1))
      count=count+1
    }else if(colnames(Subtypesign)[j] == "ALS-Glia" && Subtypesign[i,j] == "Down"){
      neworder[count]=2+(6*(i-1))
      count=count+1
    }
    
    if(colnames(Subtypesign)[j] == "ALS-Ox" && Subtypesign[i,j] == "Up"){
      neworder[count]=3+(6*(i-1))
      count=count+1
    }else if(colnames(Subtypesign)[j] == "ALS-Ox" && Subtypesign[i,j] == "Down"){
      neworder[count]=4+(6*(i-1))
      count=count+1
    }
    
    if(colnames(Subtypesign)[j] == "ALS-TD" && Subtypesign[i,j] == "Up"){
      neworder[count]=5+(6*(i-1))
      count=count+1
    }else if(colnames(Subtypesign)[j] == "ALS-TD" && Subtypesign[i,j] == "Down"){
      neworder[count]=6+(6*(i-1))
      count=count+1
    }
    
  }
  
}

length(neworder)
length(as.list(logSubtypepval))

customheatsub$Value[neworder]= as.list(t(Subtypepval))


customheatsub$Value = as.numeric(customheatsub$Value)

customheatsub$Value =p.adjust(customheatsub$Value,method = "bonferroni")


for(i in 1:nrow(customheatsub)){
  if(is.na(customheatsub$Value[i])){
    customheatsub$Value[i] = 1
  }else if(customheatsub$Value[i] == 1){
    customheatsub$Value[i] = 1
  }
}

customheatsub$Value = -log10(customheatsub$Value)

customheatsub$Module = as.character(customheatsub$Module)
customheatsub$Module = factor(customheatsub$Module,levels=unique(customheatsub$Module))


p = ggplot(customheatsub,aes(Module,Response)) + geom_tile(aes(fill=Value),color="black")
p = p+scale_fill_gradient2(low="blue",mid="white",high="darkgreen",limits=c(0,45))
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "white"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray50"))
p


######################################################################################################################################################################################
######################################################################################################################################################################################
#################################  Reviewer 2  #######################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################

library(ggplot2)
library(DESeq2)
library(scatterplot3d)
library(plot3D)
library(rgl)
library(survival)
library(survminer)

######################################################################################################################################################################################
#################################  QUESTION 1  #######################################################################################################################################
######################################################################################################################################################################################

#While many findings were reproduced, such as the classification of the 3 ALS subtypes and the gene expression pathways enriched in each, 
#there are also different conclusions drawn. These include the observation that TARDBP is not reduced in the ALS-TD subtype 
#(wherein Tam et al. the reduced expression of TARDBP explains the loss of repression of transposons) and the shorter survival 
#of ALS-Glia subtypes (where Tam et al. observe no differences). How can the authors reconcile these different findings?

load("D:/Jarrett/Research/Nat Comm Reviews/RData/ALSPatientStratification_UnivariateDatasets_RINSite_PeerReview.RData")

setwd("C:/Users/jeshima/Desktop/ALS Manuscript/Nature/Communications Submission/Peer Review")
TamSamples = read.csv("TamSubset.csv") #Provided in Table S2


#TARDBP expression in Tam cohort (median-of-ratios scale)

TamTARDBP = NormCounts[,colnames(NormCounts) %in% TamSamples$Tam.et.al.]
TamTARDBP = TamTARDBP[which(rownames(TamTARDBP) == "TARDBP"),]
TamTARDBP

table(colnames(NormCounts)[1:451] == ALSPheno$Subject)
TamSubtype = ALSPheno$Subtype[ALSPheno$Subject %in% TamSamples$Tam.et.al.]


#Clean up labels
for(i in 1:length(TamSubtype)){
  if(TamSubtype[i] == "GLIA"){
    TamSubtype[i] = "ALS-Glia"
  }else if(TamSubtype[i] == "OX"){
    TamSubtype[i] = "ALS-Ox"
  }else if(TamSubtype[i] == "TE"){
    TamSubtype[i] = "ALS-TD"
  }
}

plotdat = data.frame(TamSubtype,TamTARDBP)
logplotdat = data.frame(TamSubtype,log2(TamTARDBP))

#log2 scale
p = ggplot(logplotdat,aes(x=TamSubtype,y=log2.TamTARDBP.,fill=TamSubtype)) + geom_violin(kernel = "gaussian",scale="width")
p = p +ggtitle("Tam Samples - TARDBP Expression") + xlab("") + ylab("log2 Median-of-Ratios Counts")
p = p+scale_fill_manual(values = c("goldenrod1","navy","firebrick"))
p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20),plot.title = element_text(size=24))
p = p+theme(panel.background = element_rect(fill = "gray90",colour = "gray80",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "gray80"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray80"))
p = p+theme(legend.position = "none")
p = p+geom_dotplot(binaxis = 'y',stackdir = 'center',dotsize = 0.35,fill="white",stackratio = 1,binwidth = 1/15)
p = p+theme(axis.title.y=element_text(angle=90, vjust=1.5))
p = p+theme(plot.title = element_text(hjust = 0.5))
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "gray80"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray80"))
p



#Sig Testing
TamCounts = FullCount[,colnames(FullCount) %in% TamSamples$Tam.et.al.]
TamCounts

TamPheno = ALSPheno[ALSPheno$Subject %in% TamSamples$Tam.et.al.,]

table(colnames(TamCounts) == TamPheno$Subject)

rCountData4 = round(TamCounts,0)

dds4 = DESeqDataSetFromMatrix(countData = rCountData4, colData = TamPheno, design= ~ platform + Subtype, tidy=F) #Subtype must be second (DESeq2 vignette)
dds4$SubjectSex = relevel(dds4$Subtype,ref = "Control")
dseq4 = DESeq(dds4,betaPrior=T)

#Glia vs TE DE... for TARDBP in the Tam cohort
tmp.res = results(dseq4,contrast = c("Subtype","GLIA","TE"))
tmp.res[which(rownames(tmp.res) == "TARDBP"),]


############################################################################################################################

#Survival in Tam cohort
setwd("C:/Users/jeshima/Desktop/ALS Manuscript/Nature/Communications Submission/Peer Review")

SurvivalDat = read.csv("Table_S10.csv")
TamSamples = read.csv("TamSubset.csv") #Provided in Table S2

TamSurvival = SurvivalDat
TamSurvival$RNAseqSample1 = gsub("-","\\.",TamSurvival$RNAseqSample1)
TamSurvival$RNAseqSample2 = gsub("-","\\.",TamSurvival$RNAseqSample2)
TamSurvival$RNAseqSample3 = gsub("-","\\.",TamSurvival$RNAseqSample3)

tampatients1 = tampatients2 = tampatients3 = rep(NA,nrow(TamSurvival))
count1 = count2 = count3 = 1
for(i in 1:nrow(TamSurvival)){
  
  for(j in 1:nrow(TamSamples)){
    if(TamSurvival$RNAseqSample1[i]==TamSamples$Tam.et.al.[j]){
      tampatients1[count1] = i
      count1 = count1+1
    }
    
    if(! is.na(TamSurvival$RNAseqSample2[i])){
      if(TamSurvival$RNAseqSample2[i] == TamSamples$Tam.et.al.[j]){
        tampatients2[count2] = i
        count2 = count2+1
      }
    }
    
    if(! is.na(TamSurvival$RNAseqSample3[i])){
      if(TamSurvival$RNAseqSample3[i] == TamSamples$Tam.et.al.[j]){
        tampatients3[count3] = i
        count3 = count3+1
      }
    }
  }
  
}

TamIndex = as.numeric(names(table(c(tampatients1,tampatients2,tampatients3))))

TamSurvival = TamSurvival[TamIndex,]


#Add status
status = rep(NA,nrow(TamSurvival))
for(i in 1:nrow(TamSurvival)){
  if(is.na(TamSurvival$AgeofOnset[i])){
    status[i] = 0
  }else if (is.na(TamSurvival$Duration[i])){
    status[i] = 0
  }else if(is.na(TamSurvival$AgeofDeath[i])){
    status[i] = 0
  }else{
    status[i]= 1
  }
}

FinalSurvival = cbind(TamSurvival,status)
FinalSurvival$Duration = as.numeric(FinalSurvival$Duration)

#With Discordant
# km = with(FinalSurvival,Surv(Duration,status))
# km_fit = survfit(Surv(Duration, status) ~ FinalSubtype, data=FinalSurvival)
# surv_pvalue(km_fit)
# print(km_fit)
# summary(km_fit, times = c(1,30,60,90*(1:10)))
# 
# par(mar = c(5, 5, 3, 2))
# plot(km_fit,col = c("darkorchid2","goldenrod","navy","firebrick"),main="Tam Patients - ALS Subtype Survival",xlab="Disease Duration (months)",ylab="Survival Probability",lwd=2.5,cex.axis = 1.5,cex.lab=1.5,cex.main=1.5,xaxt="n")
# timeseq = seq(0,156,12)
# axis(side = 1, at = timeseq,labels = T,tick = T,cex.axis = 1.5)
# legend(132,1.02,legend=c("ALS-Glia","ALS-Ox","ALS-TD","Discordant"),lty = 1,lwd=3.5,col=c("goldenrod","navy","firebrick","darkorchid2"))


#Without Discordant (Tam et al. replicated)
DI = which(FinalSurvival$FinalSubtype == "Discordant")
FinalSurvival2 = FinalSurvival[-DI,]
FinalSurvival$All = "All"

km_fit = survfit(Surv(Duration, status) ~ FinalSubtype, data=FinalSurvival2) #Without discordant
surv_pvalue(km_fit)
km_fit_All = survfit(Surv(Duration, status) ~ All, data=FinalSurvival) #All patients (including discordant)


custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5,size=16)
    )
}


fit = list(Subtype = km_fit,All = km_fit_All)
p = ggsurvplot_combine(fit,FinalSurvival2,legend = "none",palette = c("goldenrod","navy","firebrick","black"),ggtheme = custom_theme()) 
p = p +ggtitle("Tam Survival") + xlab("Time (months)")
p$plot = p$plot+ggplot2::annotate("text",x=150,y=1,label = "ALS-Glia",colour = "goldenrod",cex=5)
p$plot = p$plot+ggplot2::annotate("text",x=150,y=0.9,label = "ALS-Ox",colour = "navy",cex=5)
p$plot = p$plot+ggplot2::annotate("text",x=150,y=0.8,label = "ALS-TD",colour = "firebrick",cex=5)
p$plot = p$plot+ggplot2::annotate("text",x=150,y=0.7,label = "All Patients",colour = "black",cex=5)
p


######################################################################################################################################################################################
#################################  QUESTION 2  #######################################################################################################################################
######################################################################################################################################################################################

#Clustering: It is not clear if all or some of the 1475 TEs are included in the 10,000 most variable features used for clustering. 
#Also, what is the extent of overlap of the 10,000 features between the HiSeq and Novaseq data sets?

NovaFeatures = read.csv("NovaSeq_GeneList.csv") #Publicly available at: https://figshare.com/authors/Jarrett_Eshima/13813720
HiFeatures = read.csv("HiSeq_GeneList.csv") #Publicly available at: https://figshare.com/authors/Jarrett_Eshima/13813720

Overlap = NovaFeatures$Gene[NovaFeatures$Gene %in% HiFeatures$Gene]
OverlapPercent = length(Overlap)/length(NovaFeatures$Gene)

NovaTEs = NovaFeatures$Gene[which(nchar(NovaFeatures$Gene)>25)]
NovaTEs[1:10]

HiTEs = HiFeatures$Gene[which(nchar(HiFeatures$Gene)>25)]
HiTEs[1:10]

NovaTEPercent = length(NovaTEs)/1474 #provided in methods section
HiTEPercent = length(HiTEs)/1474

OverlapTEs = NovaTEs[NovaTEs %in% HiTEs]
OverlapTEsPercent = length(OverlapTEs)/1474


######################################################################################################################################################################################
#################################  QUESTION 3  #######################################################################################################################################
######################################################################################################################################################################################
#Feature Selection: Why were the top 1,000 features selected from the 10,000 MAD transcripts? 
#The remaining 1681 genes and TEs used for enrichment and network building are much lower than what is typically used for these analyses.

ALSCount = FullCount[,1:451]
table(colnames(ALSCount) == ALSPheno$Subject)

rCountData = round(ALSCount,0)

dds = DESeqDataSetFromMatrix(countData = rCountData, colData = ALSPheno, design= ~ sequencing_platform, tidy=F)
dds$sequencing_platform = relevel(dds$sequencing_platform,ref = "NovaSeq")
dseq = DESeq(dds,betaPrior=T)

#Pairwise "contrast()"
platform.res = results(dseq,contrast = c("sequencing_platform","NovaSeq","HiSeq"))
platform.res = platform.res[! is.na(platform.res$padj) & platform.res$padj<0.05,]
head(platform.res)

PlatformDepPercent = nrow(platform.res)/nrow(ALSCount)


######################################################################################################################################################################################
#################################  QUESTION 4  #######################################################################################################################################
######################################################################################################################################################################################
#Figure 1: Do ALS/FTLD cases segregate? 

library(DESeq2)
load("D:/Jarrett/Research/Nat Comm Reviews/RData/ALSPatientStratification_UnivariateDatasets_RINSite_PeerReview.RData")


dim(FullCount);dim(FullPheno)

#All FTLD patients were characterized using NovaSeq platform, so only the NovaSeq group needs to be run on NMF

HCind = which(FullPheno$Subtype == "Control")

table(colnames(FullCount) == FullPheno$Subject)

FullCount.ALS.FTLD = FullCount[,-HCind]
FullPheno.ALS.FTLD = FullPheno[-HCind,]

Novaind = which(FullPheno.ALS.FTLD$platform == "NovaSeq")
NovaCount.ALS.FTLD = FullCount.ALS.FTLD[,Novaind]
NovaPheno.ALS.FTLD = FullPheno.ALS.FTLD[Novaind,]
Hiind = which(FullPheno.ALS.FTLD$platform == "HiSeq")
HiCount.ALS.FTLD = FullCount.ALS.FTLD[,Hiind]
HiPheno.ALS.FTLD = FullPheno.ALS.FTLD[Hiind,]

#VST Normalization - NovaSeq
rCountData = round(NovaCount.ALS.FTLD,0)
ddsnova = DESeqDataSetFromMatrix(countData = rCountData, colData = NovaPheno.ALS.FTLD, design= ~ SubjectSex, tidy=F)
ddsnova$SubjectSex = relevel(ddsnova$SubjectSex,ref = "Male")
dseqnova = DESeq(ddsnova,betaPrior=T)
resnova = results(dseqnova)
signova = resnova[! is.na(resnova$padj) & resnova$padj<0.05,]
vsdnova = varianceStabilizingTransformation(dseqnova)
Nova.vstcounts = assay(vsdnova)


#Get feature set from Fig 1
setwd("C:/Users/jeshima/Desktop/ALS Manuscript/Nature/Communications Submission/Peer Review")

Nova_ALS_MAD10k_features = read.csv("FINAL_NOVASEQ_ALS451_MAD10k_10-4-21_SYMBOL.csv") #File generated in DifferentialExpression.R script
Nova_ALS_MAD10k_features = Nova_ALS_MAD10k_features$Gene

gene_positions = read.csv("ENSEMBL_HGNC_LUT_Dec2016Archive_PeerReview.csv") #Added to supplemental table S2
gene_positions = gene_positions[,-3]
genesymbols = gene_positions$hgnc_symbol

##Correct all feature set lists for excel autoformat driven errors (MARCH3 --> 3-Mar)
myfeaturesets = c("Nova_ALS_MAD10k_features","Hi_ALS_MAD10k_features","genesymbols")
excelgenelist = c("1-Sep","2-Sep","3-Sep","4-Sep","5-Sep","6-Sep","7-Sep","8-Sep","9-Sep","10-Sep","11-Sep","12-Sep","1-Mar","2-Mar","3-Mar","4-Mar","5-Mar","6-Mar","7-Mar","8-Mar","9-Mar","10-Mar","11-Mar","12-Mar")
truegenelist = c("SEPT1","SEPT2","SEPT3","SEPT4","SEPT5","SEPT6","SEPT7","SEPT8","SEPT9","SEPT10","SEPT11","SEPT12","MARCH1","MARCH2","MARCH3","MARCH4","MARCH5","MARCH6","MARCH7","MARCH8","MARCH9","MARCH10","MARCH11","MARCH12")

ExcelerrorLUT = data.frame(cbind(excelgenelist,truegenelist))

for(i in 1:length(myfeaturesets)){
  
  tmp = get(myfeaturesets[i])
  
  for(j in 1:length(tmp)){
    for(k in 1:nrow(ExcelerrorLUT)){
      
      if(tmp[j] == ExcelerrorLUT$excelgenelist[k]){
        
        tmp[j] = ExcelerrorLUT$truegenelist[k]
        
      }
      
    }
  }
  nam = myfeaturesets[i]
  assign(paste(nam),tmp)
  if((i %% 1) == 0) cat("Feature Set Completed:",myfeaturesets[i],"\n")
} 

gene_positions$hgnc_symbol = genesymbols

#################################### NovaSeq

#Convert ENSEMBL IDs to Gene Symbols
VSTNova_Symbol = Nova.vstcounts
blank = rep(NA,nrow(VSTNova_Symbol))
VSTNova_Symbol = cbind(blank,VSTNova_Symbol)
VSTNova_Symbol[,1] = rownames(VSTNova_Symbol)

DONE = F
for(i in 1:nrow(VSTNova_Symbol)){
  for(j in 1:nrow(gene_positions)){
    if(DONE == F){
      if(sub("\\..*","",VSTNova_Symbol[i,1]) == gene_positions$ensembl_gene_id[j]){
        VSTNova_Symbol[i,1] = gene_positions$hgnc_symbol[j]
        DONE = T
      }
    }
  }
  DONE = F
  if((i %% 1000) == 0) cat("% Done:",i/nrow(VSTNova_Symbol)*100,"\n")
}

#Replace missing gene symbol with original ENSEMBL ID
for(i in 1:nrow(VSTNova_Symbol)){
  if(VSTNova_Symbol[i,1] == ""){
    VSTNova_Symbol[i,1] = rownames(VSTNova_Symbol)[i]
  }
}


uniquern = make.names(VSTNova_Symbol[,1],unique = T)
rownames(VSTNova_Symbol) = uniquern
VSTNova_Symbol_Clean = VSTNova_Symbol[,-1]

#Clean up symbols
Clean_NovaSeq_Features = gsub("-","\\.",Nova_ALS_MAD10k_features)
Clean_NovaSeq_Features = gsub("\\|","\\.",Clean_NovaSeq_Features)
Clean_NovaSeq_Features = gsub("\\:","\\.",Clean_NovaSeq_Features)
Clean_NovaSeq_Features = gsub("\\+","\\.",Clean_NovaSeq_Features)
Clean_NovaSeq_Features = gsub("\\(","\\.",Clean_NovaSeq_Features)
Clean_NovaSeq_Features = gsub("\\)","\\.",Clean_NovaSeq_Features)
Clean_NovaSeq_Features = gsub("\\?","\\.",Clean_NovaSeq_Features)


FTLD_NovaSeq_NMF = VSTNova_Symbol_Clean[rownames(VSTNova_Symbol_Clean) %in% Clean_NovaSeq_Features,]

dim(FTLD_NovaSeq_NMF)

#write.csv(FTLD_NovaSeq_NMF,"ALS_FTLD_NovaSeq_MAD10k.csv")
#save.image("ALS_FTLD_NMF.RData")
#load("D:/Jarrett/Research/Nat Comm Reviews/RData/ALS_FTLD_NMF.RData")


#Unsupervised clustering
library(sake)
shiny::runApp(system.file("sake", package="sake"))

#Save NMF result file, unzip folder, change filename ending in _Groups


#Load NMF result file - after SAKE run is completed
# ONDnames = colnames(ONDControlCounts)
# setwd("D:/Jarrett/Research/Nat Comm Reviews/FTLD Clustering/Novaseq/k4/ALS_FTLD_NovaSeq_k4")
# NovaSeq_NMF = read.csv("Groups.csv")
# NMFres_nova = NovaSeq_NMF[NovaSeq_NMF$Sample_ID %in% ONDnames,]
# table(NMFres_nova$nmf_subtypes)
# table(NovaSeq_NMF$nmf_subtypes)


######################################################################################################################################################################################
#################################  QUESTION 5  #######################################################################################################################################
######################################################################################################################################################################################
#Figure 2: A-C: Are these enriched gene sets with adjusted p-values at a > 0.05 threshold? Is this a comprehensive list of significant gene sets? 

#Addressed using Enrichr

######################################################################################################################################################################################
#################################  QUESTION 6  #######################################################################################################################################
######################################################################################################################################################################################
#Figure 4: The value of this analysis is low. NMF was already used to define the 3 clusters. 
#From this, the modules were selected based on their associations with each subtype, and the classifiers were taken from the modules, and then reapplied to classify samples from the same data set. 
#This is self-referential. Patel et al., defined classifier genes from TGCA and applied these to classify inn new data sets of single cells in tumors. 
#Despite this re-classification, what more can be said about the hybrid states phenotypically? Do they have altered survival, age of onset, etc? 
#Visually, it is difficult to distinguish border colors from fill colors on the plot. Consider showing histograms and quantification of the classifications. 

#The code used to address this comment is an extension of ALSPatientStratification.BootstrapClassification.R (https://github.com/BSmithLab/ALSPatientStratification)


#Load in dependent libraries
library(scatterplot3d)
library(plot3D)
library(rgl)


#Load in cleaned/normalized count matrix (this file is generated in the UnivariateAnalysis.R script)
load("D:/Jarrett/Research/Nat Comm Reviews/RData/ALSPatientStratification_UnivariateDatasets_RINSite_PeerReview.RData")

#Read in phenotype data
wd2 = "C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Supervised Classification/ALS_Patel"
setwd(wd2)
Pheno = read.csv("ALS451_coldata_SUBTYPES.csv") #Table S13
SampleList = Pheno$Subject

############# Step 1: Define subtype predictors and visualize preliminary subtype scores ###############################

#Purple Eigengene - Glia; Magenta Eigengene - TD; Turquoise Eigengene - Ox
GliaPredictors = read.csv("GliaFeatures_Purple.csv") #Derived from table S8
GliaPredictors = unlist(GliaPredictors)
table(table(GliaPredictors)>1)#duplicate check
OxPredictors = read.csv("OxFeatures_Turq.csv") #Derived from table S8
OxPredictors = unlist(OxPredictors)
table(table(OxPredictors)>1)#duplicate check
TDPredictors = read.csv("TDFeatures_Mag.csv") #Derived from table S8
TDPredictors = unlist(TDPredictors)
table(table(TDPredictors)>1)#duplicate check

#Features used in more than one subtype? - No
GliaPredictors[GliaPredictors %in% OxPredictors]
GliaPredictors[GliaPredictors %in% TDPredictors]
OxPredictors[OxPredictors %in% TDPredictors]

GliaAverage = NormCounts[rownames(NormCounts) %in% GliaPredictors,]
OxAverage = NormCounts[rownames(NormCounts) %in% OxPredictors,]
TDAverage = NormCounts[rownames(NormCounts) %in% TDPredictors,]

GliaAverage = data.frame(GliaAverage);OxAverage = data.frame(OxAverage);TDAverage = data.frame(TDAverage)

ALSGlia = GliaAverage[,colnames(GliaAverage) %in% SampleList]
ALSOx = OxAverage[,colnames(OxAverage) %in% SampleList]
ALSTD = TDAverage[,colnames(TDAverage) %in% SampleList]

ALSGlia = data.frame(t(ALSGlia))
ALSOx = data.frame(t(ALSOx))
ALSTD = data.frame(t(ALSTD))

AverageGliaScore = rowMeans(ALSGlia,na.rm=T)
AverageOxScore = rowMeans(ALSOx,na.rm=T)
AverageTDScore = rowMeans(ALSTD,na.rm=T)

#Average Expression of 1690 features
CleanNormCounts = NormCounts[,colnames(NormCounts) %in% SampleList] #1690 features
CleanNormCounts = t(CleanNormCounts)
dim(CleanNormCounts)
table(rownames(ALSGlia) == rownames(CleanNormCounts))
AverageExpression = rowMeans(CleanNormCounts,na.rm = T)

#Patel et al. defined subtype scores
SubtypeScore_Glia = AverageGliaScore - AverageExpression
SubtypeScore_Ox = AverageOxScore - AverageExpression
SubtypeScore_TD = AverageTDScore - AverageExpression

SubtypeScores = data.frame(matrix(NA,length(SubtypeScore_Glia),3))
rownames(SubtypeScores) = SampleList
colnames(SubtypeScores) = c("Glia","Ox","TD")

SubtypeScores$Glia = SubtypeScore_Glia
SubtypeScores$Ox = SubtypeScore_Ox
SubtypeScores$TD = SubtypeScore_TD

#Assign colors based on unsupervised clustering results
subtypecolors = shapes = rep(NA,length(SampleList))
for(i in 1:length(SampleList)){
  
  tmp = SampleList[i]
  ind = which(Pheno$Subject == tmp)
  tmp2 = Pheno$Subtype[ind]
  
  if(tmp2 == "TE"){
    subtypecolors[i] = "firebrick"
    shapes[i] = 15
  }else if(tmp2 == "OX"){
    subtypecolors[i] = "navy"
    shapes[i] = 16
  }else{
    subtypecolors[i] = "goldenrod"
    shapes[i] = 17
  }
  
}

#Visualize subtype scoring method
plot3d(SubtypeScores[,1],SubtypeScores[,2],SubtypeScores[,3],xlab = "Glia Score",ylab = "Ox Score",zlab="TD Score",col=subtypecolors,size = 7)

#################################### Step 2: Define 5% cutoff ###################################################################################

set.seed(1234)
times = 100
npredictorsG = length(GliaPredictors)
npredictorsO = length(OxPredictors)
npredictorsT = length(TDPredictors)

RandomPredictors_Glia = data.frame(matrix(NA,times,npredictorsG))
RandomPredictors_Ox = data.frame(matrix(NA,times,npredictorsO))
RandomPredictors_TD = data.frame(matrix(NA,times,npredictorsT))
FeatureList = rownames(NormCounts)

#Using sampling with replacement from subtype predictor lists ("5% cutoff")
for(i in 1:times){
  Gindices = sample(1:length(GliaPredictors),npredictorsG,replace = T)
  RandomPredictors_Glia[i,] = GliaPredictors[Gindices]
  Oindices = sample(1:length(OxPredictors),npredictorsO,replace = T)
  RandomPredictors_Ox[i,] = OxPredictors[Oindices]
  Tindices = sample(1:length(TDPredictors),npredictorsT,replace = T)
  RandomPredictors_TD[i,] = TDPredictors[Tindices]
}


#Iterate over 100 different sets of subtype predictors (according to Patel et al.)
ScoreSamplingG = ScoreSamplingO = ScoreSamplingT = data.frame(matrix(NA,length(SampleList),times))
rownames(ScoreSamplingG) = rownames(ScoreSamplingO) = rownames(ScoreSamplingT) = SampleList
colnames(ScoreSamplingG) = colnames(ScoreSamplingO) = colnames(ScoreSamplingT) = paste("Rep",seq(1,times,1),sep="")

for(i in 1:times){
  
  tmp = RandomPredictors_Glia[i,]
  RSAverage = NormCounts[rownames(NormCounts) %in% tmp,]
  RSAverage = data.frame(RSAverage)
  RSAverage = RSAverage[,colnames(RSAverage) %in% SampleList]
  RSAverage = data.frame(t(RSAverage))
  RandomScore = rowMeans(RSAverage,na.rm=T)
  SubtypeScore_cutoff = RandomScore - AverageExpression
  values = unname(SubtypeScore_cutoff)
  ScoreSamplingG[,i] = values
  
  tmp = RandomPredictors_Ox[i,]
  RSAverage = NormCounts[rownames(NormCounts) %in% tmp,]
  RSAverage = data.frame(RSAverage)
  RSAverage = RSAverage[,colnames(RSAverage) %in% SampleList]
  RSAverage = data.frame(t(RSAverage))
  RandomScore = rowMeans(RSAverage,na.rm=T)
  SubtypeScore_cutoff = RandomScore - AverageExpression
  values = unname(SubtypeScore_cutoff)
  ScoreSamplingO[,i] = values
  
  tmp = RandomPredictors_TD[i,]
  RSAverage = NormCounts[rownames(NormCounts) %in% tmp,]
  RSAverage = data.frame(RSAverage)
  RSAverage = RSAverage[,colnames(RSAverage) %in% SampleList]
  RSAverage = data.frame(t(RSAverage))
  RandomScore = rowMeans(RSAverage,na.rm=T)
  SubtypeScore_cutoff = RandomScore - AverageExpression
  values = unname(SubtypeScore_cutoff)
  ScoreSamplingT[,i] = values
  
}

ScoreSamplingG2 = matrix(unlist(ScoreSamplingG),nrow(ScoreSamplingG),ncol(ScoreSamplingG))
rownames(ScoreSamplingG2) = rownames(ScoreSamplingG)
colnames(ScoreSamplingG2) = colnames(ScoreSamplingG)
ScoreSamplingO2 = matrix(unlist(ScoreSamplingO),nrow(ScoreSamplingO),ncol(ScoreSamplingO))
rownames(ScoreSamplingO2) = rownames(ScoreSamplingO)
colnames(ScoreSamplingO2) = colnames(ScoreSamplingO)
ScoreSamplingT2 = matrix(unlist(ScoreSamplingT),nrow(ScoreSamplingT),ncol(ScoreSamplingT))
rownames(ScoreSamplingT2) = rownames(ScoreSamplingT)
colnames(ScoreSamplingT2) = colnames(ScoreSamplingT)


#Average Across Predictor Sets
repthreshG = repthreshO = repthreshT = rep(NA,nrow(ScoreSamplingG))
for(i in 1:nrow(ScoreSamplingG2)){
  
  tmpmean = mean(ScoreSamplingG2[i,])
  tmpsigma = sd(ScoreSamplingG2[i,])
  repthreshG[i] = tmpmean + 1.645*tmpsigma
  
  tmpmean = mean(ScoreSamplingO2[i,])
  tmpsigma = sd(ScoreSamplingO2[i,])
  repthreshO[i] = tmpmean + 1.645*tmpsigma
  
  tmpmean = mean(ScoreSamplingT2[i,])
  tmpsigma = sd(ScoreSamplingT2[i,])
  repthreshT[i] = tmpmean + 1.645*tmpsigma
}

#Thresholds values are with respect to unsupervised clustering subtype proportions
table(Pheno$Subtype)
Gtmp = repthreshG[order(repthreshG,decreasing = T)][84] #84 of 451 are in the Glia cluster
Otmp = repthreshO[order(repthreshO,decreasing = T)][239] #239 of 451 are in the Ox cluster
Ttmp = repthreshT[order(repthreshT,decreasing = T)][128] #128 of 451 are in the TD cluster

#Visualize thresholds
par(mfrow=c(1,3))
hist(repthreshG); abline(v = Gtmp,col="red")
hist(repthreshO); abline(v = Otmp,col="red")
hist(repthreshT); abline(v = Ttmp,col="red")

#Initialize threshold reference container
SubtypeThresholds = data.frame(matrix(NA,nrow(ScoreSamplingG2),4))
colnames(SubtypeThresholds) = c("SampleID","GThresh","OThresh","TDThresh")
SubtypeThresholds$SampleID = rownames(ScoreSamplingG2)
SubtypeThresholds$GThresh = Gtmp
SubtypeThresholds$OThresh = Otmp
SubtypeThresholds$TDThresh = Ttmp


#Pre-calculate average Expression of 1690 features
CleanNormCounts = NormCounts[,colnames(NormCounts) %in% SampleList] #1690 features
CleanNormCounts = t(CleanNormCounts)
table(rownames(ALSGlia) == rownames(CleanNormCounts))
AverageExpression = rowMeans(CleanNormCounts,na.rm = T)

#Filter out control samples
ALSCounts = NormCounts[,colnames(NormCounts) %in% SampleList]
table(colnames(ALSCounts) == SubtypeThresholds$SampleID)

############################### Step 3: Calculate subtype scores ###################################################################################

iter = 1000 #n iterations in bootstrap
GliaScores = OxScores = TDScores = data.frame(matrix(NA,nrow(SubtypeThresholds),iter)) #containers
rownames(GliaScores) = rownames(OxScores) = rownames(TDScores) = SubtypeThresholds$SampleID

for(i in 1:iter){
  sampGP = sample(GliaPredictors,npredictorsG,replace = T)
  sampOP = sample(OxPredictors,npredictorsO,replace = T)
  sampTP = sample(TDPredictors,npredictorsT,replace = T)
  
  GliaAverage = data.frame(matrix(NA,nrow = length(sampGP),ncol = ncol(ALSCounts)))
  OxAverage = data.frame(matrix(NA,nrow = length(sampOP),ncol = ncol(ALSCounts)))
  TDAverage = data.frame(matrix(NA,nrow = length(sampTP),ncol = ncol(ALSCounts)))
  
  for(j in 1:length(sampGP)){
    index = which(rownames(ALSCounts) == sampGP[j])
    GliaAverage[j,] = ALSCounts[index,]
  }
  
  for(k in 1:length(sampOP)){
    index = which(rownames(ALSCounts) == sampOP[k])
    OxAverage[k,] = ALSCounts[index,]
  }
  
  for(l in 1:length(sampTP)){
    index = which(rownames(ALSCounts) == sampTP[l])
    TDAverage[l,] = ALSCounts[index,]
  }
  
  ALSGlia = data.frame(t(GliaAverage))
  ALSOx = data.frame(t(OxAverage))
  ALSTD = data.frame(t(TDAverage))
  
  AverageGliaScore = rowMeans(ALSGlia,na.rm=T)
  AverageOxScore = rowMeans(ALSOx,na.rm=T)
  AverageTDScore = rowMeans(ALSTD,na.rm=T)
  
  #Patel et al. defined subtype scores
  SubtypeScore_Glia = AverageGliaScore - AverageExpression
  SubtypeScore_Ox = AverageOxScore - AverageExpression
  SubtypeScore_TD = AverageTDScore - AverageExpression
  
  GliaScores[,i] = SubtypeScore_Glia
  OxScores[,i] = SubtypeScore_Ox
  TDScores[,i] = SubtypeScore_TD
  
  if((i %% 10) == 0) cat("% Done:",i/iter*100,"\n")
}

#Assess sample expression in the context of subtype thresholds
table(rownames(GliaScores) == SubtypeThresholds$SampleID) #Must be all true
BootClassification = data.frame(matrix(NA,nrow(SubtypeThresholds),3))
colnames(BootClassification) = c("GliaScore","OxScore","TDScore")
rownames(BootClassification) = SubtypeThresholds$SampleID
for(i in 1:length(AverageExpression)){
  GliaCounter = OxCounter = TDCounter = 0
  for(j in 1:iter){
    
    if(GliaScores[i,j] > SubtypeThresholds$GThresh[i]){
      GliaCounter = GliaCounter + 1
    }
    
    if(OxScores[i,j] > SubtypeThresholds$OThresh[i]){
      OxCounter = OxCounter + 1
    }
    
    if(TDScores[i,j] > SubtypeThresholds$TDThresh[i]){
      TDCounter = TDCounter + 1
    }
    
  }
  
  BootClassification$GliaScore[i] = GliaCounter/iter
  BootClassification$OxScore[i] = OxCounter/iter
  BootClassification$TDScore[i] = TDCounter/iter
  if((i %% 10) == 0) cat("% Done:",i/length(AverageExpression)*100,"\n")
}

plot3d(BootClassification$GliaScore,BootClassification$OxScore,BootClassification$TDScore,xlab = "Glia Score",ylab = "Ox Score",zlab="TD Score",col=subtypecolors,size = 7)

############################### Step 4: Establish color code ###################################################################################

#Color Code for Scoring-based Classification (similar to Patel et al.)
classcolors = rep(NA,nrow(BootClassification))
for(i in 1:length(classcolors)){
  
  if(BootClassification$GliaScore[i] > 0.5){
    if(BootClassification$OxScore[i] <= 1 && BootClassification$OxScore[i] > 0.4){
      classcolors[i] = "chartreuse2"
    }else if(BootClassification$TDScore[i] <= 1 && BootClassification$TDScore[i] > 0.4){
      classcolors[i] = "darkorange2"
    }else if(BootClassification$OxScore[i] < 0.4 && BootClassification$TDScore[i] < 0.4){  
      classcolors[i] = "goldenrod"
    }else{
      classcolors[i] = "gray50"
    }
  }else if(BootClassification$OxScore[i] > 0.5){
    if(BootClassification$GliaScore[i] <= 1 && BootClassification$GliaScore[i] > 0.4){
      classcolors[i] = "chartreuse2"
    }else if(BootClassification$TDScore[i] <= 1 && BootClassification$TDScore[i] > 0.4){
      classcolors[i] = "darkorchid2"
    }else if(BootClassification$GliaScore[i] < 0.4 && BootClassification$TDScore[i] < 0.4){
      classcolors[i] = "navy"
    }else{
      classcolors[i] = "gray50"
    }
  }else if(BootClassification$TDScore[i] > 0.5){
    if(BootClassification$GliaScore[i] <= 1 && BootClassification$GliaScore[i] > 0.4){
      classcolors[i] = "darkorange2"
    }else if(BootClassification$OxScore[i] <= 1  && BootClassification$OxScore[i] > 0.4){
      classcolors[i] = "darkorchid2"
    }else if(BootClassification$GliaScore[i] < 0.4 && BootClassification$OxScore[i] < 0.4){
      classcolors[i] = "firebrick"
    }else{
      classcolors[i] = "gray50"
    }
  }else if(BootClassification$GliaScore[i] <= 1 && BootClassification$GliaScore[i] > 0.4){
    if(BootClassification$OxScore[i] <= 1 && BootClassification$OxScore[i] > 0.4){
      classcolors[i] = "chartreuse2"
    }else if(BootClassification$TDScore[i] <= 1 && BootClassification$TDScore[i] > 0.4){
      classcolors[i] = "darkorange2"
    }else{
      classcolors[i] = "gray50"
    }
  }else if(BootClassification$OxScore[i] <= 1 && BootClassification$OxScore[i] > 0.4){
    if(BootClassification$GliaScore[i] <= 1 && BootClassification$GliaScore[i] > 0.4){
      classcolors[i] = "chartreuse2"
    }else if(BootClassification$TDScore[i] <= 1 && BootClassification$TDScore[i] > 0.4){
      classcolors[i] = "darkorchid2"
    }else{
      classcolors[i] = "gray50"
    }
  }else if(BootClassification$TDScore[i] <= 1 && BootClassification$TDScore[i] > 0.4){
    if(BootClassification$GliaScore[i] <= 1 && BootClassification$GliaScore[i] > 0.4){
      classcolors[i] = "darkorange2"
    }else if(BootClassification$OxScore[i] <= 1  && BootClassification$OxScore[i] > 0.4){
      classcolors[i] = "darkorchid2"
    }else{
      classcolors[i] = "gray50"
    }
  }else{ 
    classcolors[i] = "gray50"
  }
  
}

############################################## Step 5: Plot ###################################################################################

#Plot Subtype Scores
plot3d(BootClassification$GliaScore,BootClassification$TDScore,BootClassification$OxScore,xlab = "",zlab = "",ylab="",col=classcolors,pch=20,size=7,axes=F)
#Add bounding grid
axis3d("x",pos = c(0,0,0),tick = F,labels = F)
axis3d("y",pos = c(0,0,0),tick = F,labels = F)
axis3d("z",pos = c(0,0,0),tick = F,labels = F)
axis3d("x",pos = c(1,1,0),tick = F,labels = F)
axis3d("x",pos = c(1,0,1),tick = F,labels = F)
axis3d("y",pos = c(1,1,0),tick = F,labels = F)
axis3d("y",pos = c(0,1,1),tick = F,labels = F)
segments3d(x=c(1,1),y=c(0,0),z=c(0,1))
segments3d(x=c(0,0),y=c(1,1),z=c(0,1))
#Overlay Clustering information 
pch3d(BootClassification$GliaScore,BootClassification$TDScore,BootClassification$OxScore,pch = 21, bg = classcolors,color = subtypecolors,cex=0.1)
#Save Plot
#rgl.postscript('Patel_ALSSubtypeScore_3DScatter_Eigen2.pdf',fmt = 'pdf')


#Empty Grid for Legend
# plot3d()
# axis3d("x",pos = c(0,0,0),tick = F,labels = F)
# axis3d("y",pos = c(0,0,0),tick = F,labels = F)
# axis3d("z",pos = c(0,0,0),tick = F,labels = F)
# axis3d("x",pos = c(1,1,0),tick = F,labels = F)
# axis3d("x",pos = c(1,0,1),tick = F,labels = F)
# axis3d("y",pos = c(1,1,0),tick = F,labels = F)
# axis3d("y",pos = c(0,1,1),tick = F,labels = F)
# segments3d(x=c(1,1),y=c(0,0),z=c(0,1))
# segments3d(x=c(0,0),y=c(1,1),z=c(0,1))
# 
# 
# rgl.postscript('Patel_ALSSubtypeScore_3DScatter_Grid.pdf',fmt = 'pdf')


######################## Step 6: Clustering Visualization ###################################################################################

par(mfrow=c(1,1))

#Sample order is the same for all three
classcolors
subtypecolors
BootClassification

table(rownames(BootClassification) == Pheno$Subject)

ClustSubtype = ALSPheno$Subtype

for(i in 1:length(ClustSubtype)){
  
  if(ClustSubtype[i] == "GLIA"){
    ClustSubtype[i] = "ALS-Glia"
  }else if(ClustSubtype[i] == "OX"){
    ClustSubtype[i] = "ALS-Ox"
  }else if(ClustSubtype[i] == "TE"){
    ClustSubtype[i] = "ALS-TD"
  }
  
}

#ALS-Glia / ALS-TD Hybrids
par(mar = c(5, 6.75, 3, 2))
GliaTDind = which(classcolors == "darkorange2")
GliaTDScores = BootClassification[GliaTDind,]
GliaTDSamples = rownames(GliaTDScores)
GliaTDClust = ClustSubtype[GliaTDind]
GliaTDbar = table(factor(GliaTDClust,levels = c("ALS-Glia","ALS-Ox","ALS-TD")))
barplot(GliaTDbar,col=c("goldenrod","navy","firebrick"),ylim=c(0,16),cex.axis = 3,cex = 2.5,cex.lab=3)
title(ylab = "Frequency",mgp=c(4.5,1,0),cex.lab = 3)


#ALS-Glia / ALS-Ox Hybrids
GliaOxind = which(classcolors == "chartreuse2")
GliaOxScores = BootClassification[GliaOxind,]
GliaOxSamples = rownames(GliaOxScores)
GliaOxClust = ClustSubtype[GliaOxind]
GliaOxbar = table(factor(GliaOxClust,levels = c("ALS-Glia","ALS-Ox","ALS-TD")))
barplot(GliaOxbar,col=c("goldenrod","navy","firebrick"),ylim=c(0,16),cex.axis = 3,cex = 2.5,cex.lab=3)
title(ylab = "Frequency",mgp=c(4.5,1,0),cex.lab = 3)


#ALS-Glia
Gliaind = which(classcolors == "goldenrod")
GliaSamples = rownames(BootClassification)[Gliaind]
GliaClust = ClustSubtype[Gliaind]
Gliabar = table(factor(GliaClust,levels = c("ALS-Glia","ALS-Ox","ALS-TD")))
barplot(Gliabar,col=c("goldenrod","navy","firebrick"),main = "ALS-Glia: Clustering Subtypes",ylab="Frequency",ylim=c(0,16),cex.axis = 1.5,cex = 1.5,cex.main = 1.5,cex.lab=1.5)


#ALS-Ox
Oxind = which(classcolors == "navy")
OxSamples = rownames(BootClassification)[Oxind]
OxClust = ClustSubtype[Oxind]
Oxbar = table(factor(OxClust,levels = c("ALS-Glia","ALS-Ox","ALS-TD")))
barplot(Oxbar,col=c("goldenrod","navy","firebrick"),main = "ALS-Ox: Clustering Subtypes",ylab="Frequency",ylim=c(0,120),cex.axis = 1.5,cex = 1.5,cex.main = 1.5,cex.lab=1.5)


#ALS-TD
TDind = which(classcolors == "firebrick")
TDSamples = rownames(BootClassification)[TDind]
TDClust = ClustSubtype[TDind]
TDbar = table(factor(TDClust,levels = c("ALS-Glia","ALS-Ox","ALS-TD")))
barplot(TDbar,col=c("goldenrod","navy","firebrick"),main = "ALS-TD: Clustering Subtypes",ylab="Frequency",ylim=c(0,60),cex.axis = 1.5,cex = 1.5,cex.main = 1.5,cex.lab=1.5)


#Unclassified
Uncind = which(classcolors == "gray50")
UncSamples = rownames(BootClassification)[Uncind]
UncClust = ClustSubtype[Uncind]
Uncbar = table(factor(UncClust,levels = c("ALS-Glia","ALS-Ox","ALS-TD")))
barplot(Uncbar,col=c("goldenrod","navy","firebrick"),main = "Unclassified: Clustering Subtypes",ylab="Frequency",ylim=c(0,120),cex.axis = 2,cex = 2,cex.main = 2,cex.lab=2)


################################ Step 7: Clinical Assessment ###################################################################################

library(survival)
library(survminer)

### Survival

#Survival in Tam cohort
setwd("C:/Users/jeshima/Desktop/ALS Manuscript/Nature/Communications Submission/Peer Review")

SurvivalDat = read.csv("Table_S10.csv")

HybridSurvival = SurvivalDat
HybridSurvival$RNAseqSample1 = gsub("-","\\.",HybridSurvival$RNAseqSample1)
HybridSurvival$RNAseqSample2 = gsub("-","\\.",HybridSurvival$RNAseqSample2)
HybridSurvival$RNAseqSample3 = gsub("-","\\.",HybridSurvival$RNAseqSample3)

#The hybrid sample labels are provided so that Steps 1-6 may be skipped
GliaTDSamples = c("CGND.HRA.01403","CGND.HRA.01581","CGND.HRA.01646","CGND.HRA.02157","CGND.HRA.01736",
                  "CGND.HRA.01828","CGND.HRA.01859","CGND.HRA.00559","CGND.HRA.01263","CGND.HRA.00225",
                  "CGND.HRA.00388","CGND.HRA.00359","CGND.HRA.00587","CGND.HRA.00381","CGND.HRA.00387",
                  "CGND.HRA.00366","CGND.HRA.00563","CGND.HRA.00599","CGND.HRA.00615")
GliaOxSamples = c("CGND.HRA.00022","CGND.HRA.00092","CGND.HRA.00223","CGND.HRA.00313","CGND.HRA.00250")

#ALS-Glia / ALS-TD
HybridSamples = GliaTDSamples

hybpatients1 = hybpatients2 = hybpatients3 = rep(NA,nrow(HybridSurvival))
count1 = count2 = count3 = 1
for(i in 1:nrow(HybridSurvival)){
  
  for(j in 1:length(HybridSamples)){
    if(HybridSurvival$RNAseqSample1[i]==HybridSamples[j]){
      hybpatients1[count1] = i
      count1 = count1+1
    }
    
    if(! is.na(HybridSurvival$RNAseqSample2[i])){
      if(HybridSurvival$RNAseqSample2[i] == HybridSamples[j]){
        hybpatients2[count2] = i
        count2 = count2+1
      }
    }
    
    if(! is.na(HybridSurvival$RNAseqSample3[i])){
      if(HybridSurvival$RNAseqSample3[i] == HybridSamples[j]){
        hybpatients3[count3] = i
        count3 = count3+1
      }
    }
  }
  
}

GliaTDIndex = as.numeric(names(table(c(hybpatients1,hybpatients2,hybpatients3))))

GliaTDSurvival = HybridSurvival[GliaTDIndex,]
mean(GliaTDSurvival$Duration)
sd(GliaTDSurvival$Duration)/sqrt(length(GliaTDSurvival$Duration))


for(i in 1:length(GliaTDSurvival$FinalSubtype)){
  
  if(GliaTDSurvival$FinalSubtype[i] == "GLIA"){
    GliaTDSurvival$FinalSubtype[i] = "ALS-Glia"
  }else if(GliaTDSurvival$FinalSubtype[i] == "OX"){
    GliaTDSurvival$FinalSubtype[i] = "ALS-Ox"
  }else if(GliaTDSurvival$FinalSubtype[i] == "TE"){
    GliaTDSurvival$FinalSubtype[i] = "ALS-TD"
  }
}

#Patient-wise subtypes
PatientGliaTD = table(factor(GliaTDSurvival$FinalSubtype,levels = c("ALS-Glia","ALS-Ox","ALS-TD","Discordant")))
barplot(PatientGliaTD,col=c("goldenrod","navy","firebrick","darkorchid2"),main = "ALS-Glia + ALS-TD Hybrid: Patient Subtypes",ylab="Frequency",ylim=c(0,12),cex.axis = 1.5,cex = 1.5,cex.main = 1.5,cex.lab=1.5)

HybridSurvival$FinalSubtype[GliaTDIndex] = "ALS-Glia + ALS-TD Hybrid" #Patients with 1 or more hybrid samples are assigned the hybrid subtype

HybridSamples = GliaOxSamples

hybpatients1 = hybpatients2 = hybpatients3 = rep(NA,nrow(HybridSurvival))
count1 = count2 = count3 = 1
for(i in 1:nrow(HybridSurvival)){
  
  for(j in 1:length(HybridSamples)){
    if(HybridSurvival$RNAseqSample1[i]==HybridSamples[j]){
      hybpatients1[count1] = i
      count1 = count1+1
    }
    
    if(! is.na(HybridSurvival$RNAseqSample2[i])){
      if(HybridSurvival$RNAseqSample2[i] == HybridSamples[j]){
        hybpatients2[count2] = i
        count2 = count2+1
      }
    }
    
    if(! is.na(HybridSurvival$RNAseqSample3[i])){
      if(HybridSurvival$RNAseqSample3[i] == HybridSamples[j]){
        hybpatients3[count3] = i
        count3 = count3+1
      }
    }
  }
  
}

GliaOxIndex = as.numeric(names(table(c(hybpatients1,hybpatients2,hybpatients3))))

GliaOxSurvival = HybridSurvival[GliaOxIndex,]
mean(GliaOxSurvival$Duration)
sd(GliaOxSurvival$Duration)/sqrt(length(GliaOxSurvival$Duration))


for(i in 1:length(GliaOxSurvival$FinalSubtype)){
  
  if(GliaOxSurvival$FinalSubtype[i] == "GLIA"){
    GliaOxSurvival$FinalSubtype[i] = "ALS-Glia"
  }else if(GliaOxSurvival$FinalSubtype[i] == "OX"){
    GliaOxSurvival$FinalSubtype[i] = "ALS-Ox"
  }else if(GliaOxSurvival$FinalSubtype[i] == "TE"){
    GliaOxSurvival$FinalSubtype[i] = "ALS-TD"
  }
}

#Patient-wise subtypes
PatientGliaOx = table(factor(GliaOxSurvival$FinalSubtype,levels = c("ALS-Glia","ALS-Ox","ALS-TD","Discordant")))
barplot(PatientGliaOx,col=c("goldenrod","navy","firebrick","darkorchid2"),main = "ALS-Glia + ALS-Ox Hybrid: Patient Subtypes",ylab="Frequency",ylim=c(0,12),cex.axis = 1.5,cex = 1.5,cex.main = 1.5,cex.lab=1.5)

HybridSurvival$FinalSubtype[GliaOxIndex] = "ALS-Glia + ALS-Ox Hybrid" #Patients with 1 or more hybrid samples are assigned the hybrid subtype

for(i in 1:length(HybridSurvival$FinalSubtype)){
  
  if(HybridSurvival$FinalSubtype[i] == "GLIA"){
    HybridSurvival$FinalSubtype[i] = "ALS-Glia"
  }else if(HybridSurvival$FinalSubtype[i] == "OX"){
    HybridSurvival$FinalSubtype[i] = "ALS-Ox"
  }else if(HybridSurvival$FinalSubtype[i] == "TE"){
    HybridSurvival$FinalSubtype[i] = "ALS-TD"
  }
}

#With Discordant

#Add status
status = rep(NA,nrow(HybridSurvival))
for(i in 1:nrow(HybridSurvival)){
  if(is.na(HybridSurvival$AgeofOnset[i])){
    status[i] = 0
  }else if (is.na(HybridSurvival$Duration[i])){
    status[i] = 0
  }else if(is.na(HybridSurvival$AgeofDeath[i])){
    status[i] = 0
  }else{
    status[i]= 1
  }
}

FinalSurvival = cbind(HybridSurvival,status)
FinalSurvival$Duration = as.numeric(FinalSurvival$Duration)

km = with(FinalSurvival,Surv(Duration,status))
km_fit = survfit(Surv(Duration, status) ~ FinalSubtype, data=FinalSurvival)
surv_pvalue(km_fit)
print(km_fit)
summary(km_fit, times = c(1,30,60,90*(1:10)))

par(mar = c(5, 5, 3, 2))
plot(km_fit,col = c("goldenrod","chartreuse2","darkorange2","navy","firebrick","gray50"),main="Hybrid Subtype Survival",ylim = c(-0.02,1.03),xlab="Disease Duration (months)",ylab="Survival Probability",lwd=2.5,cex.axis = 2,cex.lab=2,cex.main=2,xaxt="n")
timeseq = seq(0,156,12)
axis(side = 1, at = timeseq,labels = T,tick = T,cex.axis = 2)
legend(131,1.05,legend=c("ALS-Glia","ALS-Glia + ALS-Ox","ALS-Glia + ALS-TD","ALS-Ox","ALS-TD","Discordant"),lty = 1,lwd=3.5,col=c("goldenrod","chartreuse2","darkorange2","navy","firebrick","gray50"))


#Without Discordant

remind = which(HybridSurvival$FinalSubtype == "Discordant")
HybridSurvivalnodisc = HybridSurvival[-remind,]

status = rep(NA,nrow(HybridSurvivalnodisc))
for(i in 1:nrow(HybridSurvivalnodisc)){
  if(is.na(HybridSurvivalnodisc$AgeofOnset[i])){
    status[i] = 0
  }else if (is.na(HybridSurvivalnodisc$Duration[i])){
    status[i] = 0
  }else if(is.na(HybridSurvivalnodisc$AgeofDeath[i])){
    status[i] = 0
  }else{
    status[i]= 1
  }
}

FinalSurvival = cbind(HybridSurvivalnodisc,status)
FinalSurvival$Duration = as.numeric(FinalSurvival$Duration)

km = with(FinalSurvival,Surv(Duration,status))
km_fit = survfit(Surv(Duration, status) ~ FinalSubtype, data=FinalSurvival)
surv_pvalue(km_fit)
print(km_fit)
summary(km_fit, times = c(1,30,60,90*(1:10)))

par(mar = c(5, 5, 3, 2))
plot(km_fit,col = c("goldenrod","chartreuse2","darkorange2","navy","firebrick"),main="Hybrid Subtype Survival",xlab="Disease Duration (months)",ylab="Survival Probability",lwd=2.5,cex.axis = 1.5,cex.lab=1.5,cex.main=1.5,xaxt="n")
timeseq = seq(0,156,12)
axis(side = 1, at = timeseq,labels = T,tick = T,cex.axis = 1.5)
legend(124,0.99,legend=c("ALS-Glia","ALS-Glia + ALS-Ox","ALS-Glia + ALS-TD","ALS-Ox","ALS-TD"),lty = 1,lwd=3.5,col=c("goldenrod","chartreuse2","darkorange2","navy","firebrick"))


#Pairwise comparisons for Glia-TD hybrid
table(HybridSurvivalnodisc$FinalSubtype)
tmp = which(HybridSurvivalnodisc$FinalSubtype == "ALS-Glia + ALS-Ox Hybrid") #Not enough Glia-Ox hybrids to be meaningful
HybridSurvivalnodisc2 = HybridSurvivalnodisc[-tmp,]
table(HybridSurvivalnodisc2$FinalSubtype)

tmp = which(HybridSurvivalnodisc2$FinalSubtype == "ALS-Ox")
HybridSurvivalnodisc2_Glia = HybridSurvivalnodisc2[-tmp,]
tmp = which(HybridSurvivalnodisc2_Glia$FinalSubtype == "ALS-TD")
HybridSurvivalnodisc2_Glia = HybridSurvivalnodisc2_Glia[-tmp,]

tmp = which(HybridSurvivalnodisc2$FinalSubtype == "ALS-Glia")
HybridSurvivalnodisc2_Ox = HybridSurvivalnodisc2[-tmp,]
tmp = which(HybridSurvivalnodisc2_Ox$FinalSubtype == "ALS-TD")
HybridSurvivalnodisc2_Ox = HybridSurvivalnodisc2_Ox[-tmp,]

tmp = which(HybridSurvivalnodisc2$FinalSubtype == "ALS-Ox")
HybridSurvivalnodisc2_TD = HybridSurvivalnodisc2[-tmp,]
tmp = which(HybridSurvivalnodisc2_TD$FinalSubtype == "ALS-Glia")
HybridSurvivalnodisc2_TD = HybridSurvivalnodisc2_TD[-tmp,]


#Add status
status = rep(NA,nrow(HybridSurvivalnodisc2_Glia))
for(i in 1:nrow(HybridSurvivalnodisc2_Glia)){
  if(is.na(HybridSurvivalnodisc2_Glia$AgeofOnset[i])){
    status[i] = 0
  }else if (is.na(HybridSurvivalnodisc2_Glia$Duration[i])){
    status[i] = 0
  }else if(is.na(HybridSurvivalnodisc2_Glia$AgeofDeath[i])){
    status[i] = 0
  }else{
    status[i]= 1
  }
}

FinalSurvival = cbind(HybridSurvivalnodisc2_Glia,status)
FinalSurvival$Duration = as.numeric(FinalSurvival$Duration)

km = with(HybridSurvivalnodisc2_Glia,Surv(Duration,status))
km_fit = survfit(Surv(Duration, status) ~ FinalSubtype, data=HybridSurvivalnodisc2_Glia)
surv_pvalue(km_fit)
print(km_fit)

status = rep(NA,nrow(HybridSurvivalnodisc2_Ox))
for(i in 1:nrow(HybridSurvivalnodisc2_Ox)){
  if(is.na(HybridSurvivalnodisc2_Ox$AgeofOnset[i])){
    status[i] = 0
  }else if (is.na(HybridSurvivalnodisc2_Ox$Duration[i])){
    status[i] = 0
  }else if(is.na(HybridSurvivalnodisc2_Ox$AgeofDeath[i])){
    status[i] = 0
  }else{
    status[i]= 1
  }
}

FinalSurvival = cbind(HybridSurvivalnodisc2_Ox,status)
FinalSurvival$Duration = as.numeric(FinalSurvival$Duration)

km = with(HybridSurvivalnodisc2_Ox,Surv(Duration,status))
km_fit = survfit(Surv(Duration, status) ~ FinalSubtype, data=HybridSurvivalnodisc2_Ox)
surv_pvalue(km_fit)
print(km_fit)


status = rep(NA,nrow(HybridSurvivalnodisc2_TD))
for(i in 1:nrow(HybridSurvivalnodisc2_TD)){
  if(is.na(HybridSurvivalnodisc2_TD$AgeofOnset[i])){
    status[i] = 0
  }else if (is.na(HybridSurvivalnodisc2_TD$Duration[i])){
    status[i] = 0
  }else if(is.na(HybridSurvivalnodisc2_TD$AgeofDeath[i])){
    status[i] = 0
  }else{
    status[i]= 1
  }
}

FinalSurvival = cbind(HybridSurvivalnodisc2_TD,status)
FinalSurvival$Duration = as.numeric(FinalSurvival$Duration)

km = with(HybridSurvivalnodisc2_TD,Surv(Duration,status))
km_fit = survfit(Surv(Duration, status) ~ FinalSubtype, data=HybridSurvivalnodisc2_TD)
surv_pvalue(km_fit)
print(km_fit)


################################################ Other Clinical Parameters

table(HybridSurvival$FinalSubtype)
TDind = which(HybridSurvival$FinalSubtype == "ALS-TD")
Oxind = which(HybridSurvival$FinalSubtype == "ALS-Ox")
Gliaind = which(HybridSurvival$FinalSubtype == "ALS-Glia")
Dind = which(HybridSurvival$FinalSubtype == "Discordant")
GTind = which(HybridSurvival$FinalSubtype == "ALS-Glia + ALS-TD Hybrid")
GOind = which(HybridSurvival$FinalSubtype == "ALS-Glia + ALS-Ox Hybrid")

######################################### Age of onset

TDaoo = as.numeric(HybridSurvival$AgeofOnset[TDind])
TDaoo = TDaoo[!is.na(TDaoo)]
OXaoo = as.numeric(HybridSurvival$AgeofOnset[Oxind])
OXaoo = OXaoo[!is.na(OXaoo)]
GLIAaoo = as.numeric(HybridSurvival$AgeofOnset[Gliaind])
GLIAaoo = GLIAaoo[!is.na(GLIAaoo)]
GTaoo = as.numeric(HybridSurvival$AgeofOnset[GTind])
GTaoo = GTaoo[!is.na(GTaoo)]
GOaoo = as.numeric(HybridSurvival$AgeofOnset[GOind])
GOaoo = GOaoo[!is.na(GOaoo)]
Discaoo = as.numeric(HybridSurvival$AgeofOnset[Dind])
Discaoo = Discaoo[!is.na(Discaoo)]

par(mar = c(5, 6, 3, 2))
boxplot(GLIAaoo,OXaoo,TDaoo,Discaoo,GTaoo,GOaoo,xaxt="n",ylim = c(25,85),main=c("Hybrid Subtype Age of Onset"),col=c("goldenrod1","navy","firebrick","gray50","darkorange2","chartreuse2"),cex.lab=2.25,pch=20,cex=1.75,cex.axis = 1.75,cex.main=2.5)
axis(at=1:6,side=1,labels=c("ALS-Glia","ALS-Ox","ALS-TD","Discordant","Glia-TD Hybrid","Glia-Ox Hybrid"),cex.axis=1.5)
title(ylab = "Age of Onset (years)",mgp=c(3.5,1,0),cex.lab = 2.5)


#Quick ANOVA to check for post-hoc tests
Age = HybridSurvival$AgeofOnset
ST = HybridSurvival$FinalSubtype
anovadat = data.frame(cbind(Age,ST))
anovadat$ST = factor(anovadat$ST,ordered = T)
levels(anovadat$ST)

anovadat = anovadat[-which(is.na(anovadat$Age)),]

oneway = aov(Age~ST,data=anovadat)
summary(oneway) #suggests no significant differences in age of onset by subtype

#Post hoc (p < 0.05)
t.test(GOaoo,OXaoo) #limited number of Glia-Ox hybrids reduces value of this analysis
t.test(GOaoo,GLIAaoo) #limited number of Glia-Ox hybrids reduces value of this analysis
t.test(GOaoo,TDaoo) #limited number of Glia-Ox hybrids reduces value of this analysis
t.test(GOaoo,GTaoo) #limited number of Glia-Ox hybrids reduces value of this analysis
t.test(GOaoo,Discaoo) #limited number of Glia-Ox hybrids reduces value of this analysis

t.test(GTaoo,OXaoo)
t.test(GTaoo,GLIAaoo)
t.test(GTaoo,TDaoo)
t.test(GTaoo,Discaoo)


########################################### Age of death 

TDaod = as.numeric(HybridSurvival$AgeofDeath[TDind])
TDaod = TDaod[!is.na(TDaod)]
OXaod = as.numeric(HybridSurvival$AgeofDeath[Oxind])
OXaod = OXaod[!is.na(OXaod)]
GLIAaod = as.numeric(HybridSurvival$AgeofDeath[Gliaind])
GLIAaod = GLIAaod[!is.na(GLIAaod)]
GTaod = as.numeric(HybridSurvival$AgeofDeath[GTind])
GTaod = GTaod[!is.na(GTaod)]
GOaod = as.numeric(HybridSurvival$AgeofDeath[GOind])
GOaod = GOaod[!is.na(GOaod)]
Discaod = as.numeric(HybridSurvival$AgeofDeath[Dind])
Discaod = Discaod[!is.na(Discaod)]

par(mar = c(5, 6, 3, 2))
boxplot(GLIAaod,OXaod,TDaod,Discaod,GTaod,GOaod,xaxt="n",ylim = c(30,90),main=c("Hybrid Subtype Age of Death"),col=c("goldenrod1","navy","firebrick","gray50","darkorange2","chartreuse2"),cex.lab=2.25,pch=20,cex=1.75,cex.axis = 1.75,cex.main=2.5)
axis(at=1:6,side=1,labels=c("ALS-Glia","ALS-Ox","ALS-TD","Discordant","Glia-TD Hybrid","Glia-Ox Hybrid"),cex.axis=1.5)
title(ylab = "Age of Death (years)",mgp=c(3.5,1,0),cex.lab = 2.5)


#Quick ANOVA to check for post-hoc tests
Age = HybridSurvival$AgeofDeath
ST = HybridSurvival$FinalSubtype
anovadat = data.frame(cbind(Age,ST))
anovadat$ST = factor(anovadat$ST,ordered = T)
levels(anovadat$ST)

anovadat = anovadat[-which(is.na(anovadat$Age)),]

oneway = aov(Age~ST,data=anovadat)
summary(oneway) #suggests no significant differences in age of onset by subtype

#Post hoc (confounded with age of onset; p < 0.05)
t.test(GOaod,GLIAaod) #limited number of Glia-Ox hybrids reduces value of this analysis
t.test(GOaod,OXaod) #limited number of Glia-Ox hybrids reduces value of this analysis
t.test(GOaod,TDaod) #limited number of Glia-Ox hybrids reduces value of this analysis
t.test(GOaod,Discaod) #limited number of Glia-Ox hybrids reduces value of this analysis
t.test(GOaod,GTaod) #limited number of Glia-Ox hybrids reduces value of this analysis


################################ Site of onset (categories: Bulbar, Limb, Other)

GLIAsoo = HybridSurvival$SiteofOnset[Gliaind]
OXsoo = HybridSurvival$SiteofOnset[Oxind]
TDsoo = HybridSurvival$SiteofOnset[TDind]
Discsoo = HybridSurvival$SiteofOnset[Dind]
GTsoo = HybridSurvival$SiteofOnset[GTind]
GOsoo = HybridSurvival$SiteofOnset[GOind]

table(HybridSurvival$SiteofOnset)

excelplot = data.frame(table(factor(GLIAsoo,levels = c("Axial","Axial and Bulbar","Axial and Limb","Bulbar","Bulbar and Limb","Generalized","Limb","Not Applicable","Unknown"))),table(factor(OXsoo,levels = c("Axial","Axial and Bulbar","Axial and Limb","Bulbar","Bulbar and Limb","Generalized","Limb","Not Applicable","Unknown"))),table(factor(TDsoo,levels = c("Axial","Axial and Bulbar","Axial and Limb","Bulbar","Bulbar and Limb","Generalized","Limb","Not Applicable","Unknown"))),table(factor(Discsoo,levels = c("Axial","Axial and Bulbar","Axial and Limb","Bulbar","Bulbar and Limb","Generalized","Limb","Not Applicable","Unknown"))),table(factor(GTsoo,levels = c("Axial","Axial and Bulbar","Axial and Limb","Bulbar","Bulbar and Limb","Generalized","Limb","Not Applicable","Unknown"))),table(factor(GOsoo,levels = c("Axial","Axial and Bulbar","Axial and Limb","Bulbar","Bulbar and Limb","Generalized","Limb","Not Applicable","Unknown"))))
colnames(excelplot) = c("ALS-Glia","Freq","ALS-Ox","Freq","ALS-TD","Freq","Discordant","Freq","Glia-TD","Freq","Glia-Ox","Freq")

#write.csv(excelplot,"Hybrid_SiteofOnset_Table.csv")

################################ FTLD Comorbidity

Glia.FTLD = HybridSurvival$Subcategory[Gliaind]
Ox.FTLD = HybridSurvival$Subcategory[Oxind]
TD.FTLD = HybridSurvival$Subcategory[TDind]
Disc.FTLD = HybridSurvival$Subcategory[Dind]
GT.FTLD = HybridSurvival$Subcategory[GTind]
GO.FTLD = HybridSurvival$Subcategory[GOind]


GCats = names(table(HybridSurvival$Subcategory[Gliaind])) #3 and 5 are FTD
blank = rep(NA,length(Gliaind))
for(j in 1:length(Glia.FTLD)){
  if(Glia.FTLD[j] == GCats[3]){
    blank[j] = "Positive"
  }else if(Glia.FTLD[j] == GCats[5]){
    blank[j] = "Positive"
  }else{
    blank[j] = "Negative"
  }
}

Glia.FTLD.binary = blank


OCats = names(table(HybridSurvival$Subcategory[Oxind])) #5 and 6 are FTD
blank = rep(NA,length(Oxind))
for(j in 1:length(Ox.FTLD)){
  if(Ox.FTLD[j] == OCats[5]){
    blank[j] = "Positive"
  }else if(Ox.FTLD[j] == OCats[6]){
    blank[j] = "Positive"
  }else{
    blank[j] = "Negative"
  }
}

Ox.FTLD.binary = blank


TCats = names(table(HybridSurvival$Subcategory[TDind])) #6 is FTD
blank = rep(NA,length(TDind))
for(j in 1:length(TD.FTLD)){
  if(TD.FTLD[j] == TCats[6]){
    blank[j] = "Positive"
  }else{
    blank[j] = "Negative"
  }
}

TD.FTLD.binary = blank


DCats = names(table(HybridSurvival$Subcategory[Dind])) #4 is FTD
blank = rep(NA,length(Dind))
for(j in 1:length(Disc.FTLD)){
  if(Disc.FTLD[j] == DCats[4]){
    blank[j] = "Positive"
  }else{
    blank[j] = "Negative"
  }
}

Disc.FTLD.binary = blank


GTCats = names(table(HybridSurvival$Subcategory[GTind])) #3 is FTD
blank = rep(NA,length(GTind))
for(j in 1:length(GT.FTLD)){
  if(GT.FTLD[j] == GTCats[3]){
    blank[j] = "Positive"
  }else{
    blank[j] = "Negative"
  }
}

GT.FTLD.binary = blank


GOCats = names(table(HybridSurvival$Subcategory[GOind])) #No FTLD patients
GO.FTLD.binary = rep("Negative",length(GOind))


#Hybrid FTLD Comorbidity

Glia.CoMorb = length(which(Glia.FTLD.binary == "Positive"))/length(Glia.FTLD.binary)
Ox.CoMorb = length(which(Ox.FTLD.binary == "Positive"))/length(Ox.FTLD.binary)
TD.CoMorb = length(which(TD.FTLD.binary == "Positive"))/length(TD.FTLD.binary)
Disc.CoMorb = length(which(Disc.FTLD.binary == "Positive"))/length(Disc.FTLD.binary)
GT.CoMorb = length(which(GT.FTLD.binary == "Positive"))/length(GT.FTLD.binary)
GO.CoMorb = length(which(GO.FTLD.binary == "Positive"))/length(GO.FTLD.binary)

par(mar = c(5, 6, 3, 2))
CoMorbData = c(Glia.CoMorb,Ox.CoMorb,TD.CoMorb,Disc.CoMorb,GT.CoMorb,GO.CoMorb)
barplot(CoMorbData*100,main="Hybrid Subtype FTLD Comorbidity",names.arg = c("ALS-Glia","ALS-Ox","ALS-TD","Discordant","Glia-TD Hybrid","Glia-Ox Hybrid"),col=c("goldenrod1","navy","firebrick","gray50","darkorange2","chartreuse2"),ylim=c(0,25),cex.main = 2.5,cex.axis = 2,cex.names = 1.5, cex.lab=2)
title(ylab = "FTLD Comorbidity (%)",mgp=c(3.5,1,0),cex.lab = 2.5)


######################################################################################################################################################################################
#################################  QUESTION 7  #######################################################################################################################################
######################################################################################################################################################################################
#Figure 6: FTLD comorbidity ALS cases should be analyzed separately if they are to be compared to FTLD cases.

library(DESeq2)
library(ggplot2)

load("D:/Jarrett/Research/Nat Comm Reviews/RData/ALSPatientStratification_UnivariateDatasets_RINSite_PeerReview.RData")

#Add RIN and Site to Pheno
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping")
Meta = read.csv("GSE153960_MetaData.txt") #Publicly available at: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA644618&o=acc_s%3Aa
Clinical = read.csv("CLINICAL_DATA_PRUDENCIO.csv") #Available by request from NYGC
convertnames = gsub("-","\\.",Meta$sample_id_alt)
Meta$sample_id_alt = convertnames

#Add in RIN
FullPheno$RIN = NA
Clinical2 = Clinical
convertnames = gsub("-","\\.",Clinical2$ExternalSampleId)
Clinical2$ExternalSampleId = convertnames

for(i in 1:nrow(FullPheno)){
  
  for(j in 1:nrow(Clinical2)){
    
    if(FullPheno$Subject[i] == Clinical2$ExternalSampleId[j]){
      
      FullPheno$RIN[i] = Clinical2$RIN[j]
      
    }
    
  }
  
}

table(FullPheno$Subject == colnames(FullCount))



FiltMeta = Meta[Meta$sample_id_alt %in% FullPheno$Subject,]


reftable = table(FiltMeta$sample_id_alt)

removeind = rep(NA,nrow(FiltMeta)-length(reftable))
count = 1
for(i in 1:length(reftable)) {
  
  if(reftable[[i]] > 1){
    
    tmp = names(reftable[i])
    inds = which(FiltMeta$sample_id_alt == tmp)
    
    if(length(inds) == 2){
      removeind[count] = inds[1]
      count = count+1
    }else if(length(inds) == 3){
      removeind[seq(count,count+1,1)] = inds[1:2]
      count = count+2
    }else if(length(inds) > 3){
      cat("Problem")
    }
    
    
  }
  
}

SiteMeta = FiltMeta[-removeind,]



FullPheno$Site = NA

for(i in 1:nrow(FullPheno)){
  
  for(j in 1:nrow(SiteMeta)){
    
    if(FullPheno$Subject[i] == SiteMeta$sample_id_alt[j]){
      
      FullPheno$Site[i] = SiteMeta$project[j]
      
    }
    
  }
  
}

#Clean up site names
for(i in 1:nrow(FullPheno)){
  
  if(FullPheno$Site[i] == "NYGC ALS Consortium"){
    FullPheno$Site[i] = "NYGC"
  }else{
    FullPheno$Site[i] = "TargetALS"
  }
  
}

table(FullPheno$Site)



setwd("C:/Users/jeshima/Desktop/ALS Manuscript/Nature/Communications Submission/Peer Review")
tmpPheno = read.csv("ALS451_coldata_SUBTYPES.csv") #Table S13

table(tmpPheno$disease_group)

#Subtype-specific features in Fig. 6
my36 = c("AIF1","APOC2","CD44","CHI3L2","CX3CR1","FOLH1","HLA-DRA","TLR7","TMEM125","TNC","TREM2","TYROBP","COL18A1","GABRA1","GAD2","GLRA3","HTR2A","OXR1","SERPINI1","SLC6A13","SLC17A6","TCIRG1","UBQLN2","UCP2","AGPAT4-IT1","CHKB-CPT1B","COL3A1","ENSG00000205041","ENSG00000258674","ENSG00000273151","GATA2-AS1","HSP90AB4P","LINC01347","MIR24-2","MIRLET7BHG","NANOGP4")
figs12 = c("AGER","AQP1","BECN2","C1D","CCDC154","ENSG00000260198","ENSG00000278434","ENSG00000279759","ENSG00000281969","HIST1H1T","HLA-DRB1","IFI30","IRF7","SELL","SERPINA1","SNX18P3","SOCS3","STH","TNRC6C-AS1","TUNAR")
#extraglia = c("ALOX5AP","APOBR","APOC1","CCR5","CR1","CD68","CLEC7A","FPR3","MSR1","NCF2","NINJ2","ST6GALNAC2","TLR8","TNFRSF25","TREM1","VRK2")
#extraox = c("B4GALT6","BECN1","COL4A6","COX4I2","CP","GABRA6","GPR22","MYH11","MYL9","NDUFA4L2","NOS3","NOTCH3","PCSK1","SOD1","TAGLN","UBQLN1")
#extraTD = c("ADAT3","COL6A3","EGLN1P1","ENSG00000263278","ENSG00000268670","ENSG00000279233","ITGBL1","KRT8P13","LINC00176","LINC00638","MIR219A2","NKX6-2","RPS20P22","SLX1B-SULT1A4","TP63","TUB-AS1")

#Get index of FTLD+ patients

ind = rep(NA,length(tmpPheno$disease_group))
for(i in 1:length(tmpPheno$disease_group)){
  
  if(tmpPheno$disease_group[i] == "ALS/FTLD"){
    ind[i] = i
  }else if(tmpPheno$disease_group[i] == "FTLD-TDP"){
    ind[i] = i
  }
  
}

ind = ind[!is.na(ind)]

#Filter
FTLDPheno = tmpPheno[ind,]
colnames(FTLDPheno)[1] = "Subject"
ALS.FTLD = FTLDPheno$Subject

Otherind = which(FullPheno$Factor == 4)
Otherind2 = which(FullPheno$Factor == 5)
Otherindf = c(Otherind,Otherind2)

HC.FTLD.Pheno = FullPheno[Otherindf,]
Controls = HC.FTLD.Pheno$Subject

FTLDSubj = c(ALS.FTLD,Controls)

FullCount.FTLD = FullCount[,colnames(FullCount) %in% FTLDSubj]
FullPheno.FTLD = FullPheno[FullPheno$Subject %in% FTLDSubj,]

dim(FullCount.FTLD)
dim(FullPheno.FTLD)

table(colnames(FullCount.FTLD) == FullPheno.FTLD$Subject)

#No missing RINs
which(is.na(FullPheno.FTLD$RIN))

FullPheno.FTLD$RIN = scale(FullPheno.FTLD$RIN,center = T)


#Run DE
rCountData.FTLD = round(FullCount.FTLD,0)
dds.FTLD = DESeqDataSetFromMatrix(countData = rCountData.FTLD, colData = FullPheno.FTLD, design= ~ platform + RIN + Site + Subtype, tidy=F) #Subtype must be second (DESeq2 vignette)
dds.FTLD$Subtype = relevel(dds.FTLD$Subtype,ref = "Control")
dseq.FTLD = DESeq(dds.FTLD,betaPrior=T)

#Pairwise "contrast()"
glia.res = results(dseq.FTLD,contrast = c("Subtype","GLIA","Control"))
glia.sig = glia.res[! is.na(glia.res$padj) & glia.res$padj<0.05,]
filt.glia.sig = glia.sig[rownames(glia.sig) %in% Transcripts,]

ox.res = results(dseq.FTLD,contrast = c("Subtype","OX","Control"))
ox.sig = ox.res[! is.na(ox.res$padj) & ox.res$padj<0.05,]
filt.ox.sig = ox.sig[rownames(ox.sig) %in% Transcripts,]

TE.res = results(dseq.FTLD,contrast = c("Subtype","TE","Control"))
TE.sig = TE.res[! is.na(TE.res$padj) & TE.res$padj<0.05,]
filt.TE.sig = TE.sig[rownames(TE.sig) %in% Transcripts,]

GT.res = results(dseq.FTLD,contrast = c("Subtype","GLIA","TE"))
GT.sig = GT.res[! is.na(GT.res$padj) & GT.res$padj<0.05,]
filt.GT.sig = GT.sig[rownames(GT.sig) %in% Transcripts,]

GO.res = results(dseq.FTLD,contrast = c("Subtype","GLIA","OX"))
GO.sig = GO.res[! is.na(GO.res$padj) & GO.res$padj<0.05,]
filt.GO.sig = GO.sig[rownames(GO.sig) %in% Transcripts,]

TO.res = results(dseq.FTLD,contrast = c("Subtype","TE","OX"))
TO.sig = TO.res[! is.na(TO.res$padj) & TO.res$padj<0.05,]
filt.TO.sig = TO.sig[rownames(TO.sig) %in% Transcripts,]

glia.res.ond = results(dseq.FTLD,contrast = c("Subtype","GLIA","OND"))
glia.sig.ond = glia.res.ond[! is.na(glia.res.ond$padj) & glia.res.ond$padj<0.05,]
filt.glia.sig.ond = glia.sig.ond[rownames(glia.sig.ond) %in% Transcripts,]

ox.res.ond = results(dseq.FTLD,contrast = c("Subtype","OX","OND"))
ox.sig.ond = ox.res.ond[! is.na(ox.res.ond$padj) & ox.res.ond$padj<0.05,]
filt.ox.sig.ond = ox.sig.ond[rownames(ox.sig.ond) %in% Transcripts,]

TE.res.ond = results(dseq.FTLD,contrast = c("Subtype","TE","OND"))
TE.sig.ond = TE.res.ond[! is.na(TE.res.ond$padj) & TE.res.ond$padj<0.05,]
filt.TE.sig.ond = TE.sig.ond[rownames(TE.sig.ond) %in% Transcripts,]

COND.res = results(dseq.FTLD,contrast = c("Subtype","Control","OND"))
COND.sig = COND.res[! is.na(COND.res$padj) & COND.res$padj<0.05,]
filt.COND.sig = COND.sig[rownames(COND.sig) %in% Transcripts,]

#Median-of-ratios counts
tmp = estimateSizeFactors(dds.FTLD)
DESeq_NormalizedCounts_60k = counts(tmp,normalized = T)
DESeq_NormalizedCounts_60k_filtered = DESeq_NormalizedCounts_60k[rownames(DESeq_NormalizedCounts_60k) %in% Transcripts,]
NormCounts.FTLD = DESeq_NormalizedCounts_60k_filtered

#For plotting purposes (not used during statistical analysis), adjust zero count genes to 1

for(i in 1:nrow(NormCounts.FTLD)){
  for(j in 1:ncol(NormCounts.FTLD)){
    if(NormCounts.FTLD[i,j] == 0){
      NormCounts.FTLD[i,j] = 1
    }
  }
}

#Clean up naming
Subtype = FullPheno.FTLD$Subtype
for(i in 1:length(Subtype)){
  if(Subtype[i] == "GLIA"){
    Subtype[i] = "ALS-Glia"
  }else if(Subtype[i] == "OX"){
    Subtype[i] = "ALS-Ox"
  }else if(Subtype[i] == "TE"){
    Subtype[i] = "ALS-TD"
  }else if(Subtype[i] == "OND"){
    Subtype[i] = "FTLD"
  }
}

Subtype = factor(Subtype,levels = c("Control","FTLD","ALS-Glia","ALS-Ox","ALS-TD"))

#save.image("D:/Jarrett/Research/Nat Comm Reviews/RData/FTLD_Univariate_RINSiteDesign.RData")
#load("D:/Jarrett/Research/Nat Comm Reviews/RData/FTLD_Univariate_RINSiteDesign.RData")

#"AutoPlot" is the simple function version of lines 471-2643 in ALSPatientStratification_UnivariateAnalysis.R
AutoPlot(PlotGene = figs12,Focus = "HC",NormCounts = NormCounts.FTLD,Subtype = Subtype,filt.glia.sig = filt.glia.sig,filt.ox.sig = filt.ox.sig,filt.TE.sig = filt.TE.sig,filt.GT.sig = filt.GT.sig,filt.GO.sig = filt.GO.sig,filt.TO.sig = filt.TO.sig,filt.glia.sig.ond = filt.glia.sig.ond,filt.ox.sig.ond = filt.ox.sig.ond,filt.TE.sig.ond = filt.TE.sig.ond,filt.COND.sig = filt.COND.sig)


######################################################################################################################################################################################
#################################  QUESTION 8  #######################################################################################################################################
######################################################################################################################################################################################

#Figure S6: B-E: Consider performing Wilcox Rank sum to assess significance of prediction accuracy.

#AUCs Obtained from Python
KNN = c(0.6612244897959185, 0.7473448563098709, 0.5906432748538011)
MLP = c(0.7510204081632653, 0.6185964181591004, 0.7901897601145722)
RF = c(0.8111801242236024, 0.728915035401916, 0.7720491705454111)
SVM = c(0.7882874889086069, 0.6775301957517701, 0.7725265544814417)

#Rank Sum Pairwise testing
wilcox.test(KNN,MLP,paired=T)
wilcox.test(KNN,RF,paired=T)
wilcox.test(KNN,SVM,paired=T)
wilcox.test(MLP,RF,paired=T)
wilcox.test(MLP,SVM,paired=T)
wilcox.test(RF,SVM,paired=T)

#############################################

#U-statistic based approach
library(Hmisc)
#Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3575184/

#Probabilities obtained from Python - Glia vs Rest
KNN_Probs_Glia = c(0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.19933788, 0.        , 0.        , 0.        ,
                   0.19508265, 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.40347443, 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.38576447,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.38807914, 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.60001987, 0.        , 0.        , 0.        , 0.        ,
                   0.19116252, 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.20388078, 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.40405038, 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.60143341, 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.39577286, 0.        , 0.4062808 , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.1978429 , 0.        , 0.        , 0.        , 0.20100644,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.        ,
                   0.        , 0.        , 0.3930101 , 0.        , 0.        ,
                   0.        , 0.        , 0.        , 0.        , 0.61478359,
                   0.19680883)

MLP_Probs_Glia = c(1.44722907e-01, 8.67503941e-03, 9.08461334e-04, 2.59842902e-03,
                   1.07578662e-02, 9.21661098e-02, 4.75085302e-03, 4.27490467e-03,
                   3.06307072e-02, 2.72340719e-03, 1.50468943e-03, 1.22972026e-03,
                   3.89033961e-04, 8.55028633e-04, 3.66154776e-04, 6.69779352e-03,
                   1.75415969e-03, 4.22460642e-04, 1.48436950e-03, 1.70226007e-03,
                   2.74643819e-04, 4.36591748e-03, 1.11216319e-02, 2.89373203e-04,
                   2.63531038e-04, 3.71461359e-03, 6.33801209e-04, 4.43772124e-02,
                   6.92693554e-03, 3.30292572e-03, 7.41502932e-03, 1.15165433e-01,
                   6.61921334e-04, 1.40039231e-03, 5.92440152e-04, 2.28662109e-04,
                   1.12474846e-01, 3.13863785e-02, 7.30714121e-03, 2.01506748e-02,
                   5.36092102e-03, 1.65973522e-01, 2.38406282e-04, 6.89727626e-02,
                   7.89662036e-03, 1.44808618e-02, 3.10161708e-03, 7.65161830e-04,
                   7.56888214e-02, 5.30984671e-03, 3.44175169e-04, 1.37504537e-03,
                   3.41125505e-04, 3.30395970e-03, 6.87458382e-04, 7.09995627e-04,
                   1.76125093e-03, 8.07579698e-04, 3.48899386e-03, 3.29944141e-04,
                   1.13639157e-03, 3.33048114e-04, 4.34666621e-04, 5.64451912e-01,
                   1.48440656e-01, 3.80476565e-03, 5.13462828e-04, 4.16174013e-04,
                   1.37298307e-03, 1.06354511e-02, 4.38190103e-04, 2.56859234e-03,
                   3.56325323e-03, 7.13065251e-04, 6.55753794e-04, 2.56392276e-04,
                   1.05422105e-03, 6.71505636e-01, 4.49868167e-02, 1.06355567e-01,
                   5.57165153e-04, 1.76179488e-04, 2.52154582e-04, 6.89371538e-03,
                   2.33476326e-02, 7.08308467e-04, 3.33186151e-04, 1.13749974e-01,
                   3.00727698e-04, 3.55317285e-01, 3.28825044e-04, 4.67150257e-02,
                   1.48060336e-03, 5.65372231e-04, 6.61805137e-02, 1.86278129e-03,
                   1.94873270e-03, 3.52772164e-03, 6.18440082e-04, 1.51095969e-03,
                   8.72577591e-03, 1.22660667e-04, 2.48116687e-04, 1.98667208e-03,
                   1.42699126e-02, 8.12176559e-01, 1.85903609e-02, 1.71989115e-01,
                   2.87393946e-03, 3.97244695e-04, 7.86786645e-02, 1.91990085e-01,
                   5.44964958e-04, 1.78152135e-03, 1.71500259e-03, 8.74902530e-04,
                   1.69983222e-03, 2.67045268e-03, 4.50282658e-04, 1.12606713e-02,
                   8.16292859e-02, 7.27001709e-04, 5.95369770e-04, 1.72430667e-03,
                   1.18766994e-03, 3.55849282e-03, 4.04145988e-04, 5.24826591e-01,
                   3.37469369e-04, 2.97368845e-04, 2.76956878e-04, 4.96696748e-04,
                   3.75867761e-02, 4.33145715e-04, 4.54537507e-04, 6.43237103e-03,
                   3.85024113e-04, 5.10234729e-04, 4.74973622e-04, 8.27481359e-04,
                   2.95941570e-03, 7.33730720e-03, 4.11317237e-03, 2.79902104e-02,
                   9.16738882e-04, 5.27278895e-01, 7.52229179e-04, 5.69727604e-04,
                   5.42457986e-04, 2.31039292e-03, 3.64600977e-04, 2.71117493e-03,
                   3.21718162e-03, 4.69706989e-02, 2.37426691e-04, 2.23689355e-02,
                   1.77837640e-02, 5.71497996e-01, 1.31535409e-02, 3.72114500e-03,
                   2.43580080e-02, 8.31463712e-03, 2.10329515e-02, 2.63849511e-02,
                   3.19673849e-03, 2.31684814e-02, 4.62531338e-03, 2.28338458e-02,
                   6.84857659e-03, 4.62984075e-04, 8.05406642e-02, 1.10226130e-02,
                   2.63107002e-03, 6.60889636e-04, 1.66873632e-02, 5.17835176e-04,
                   5.08464086e-04, 4.64109200e-02, 2.48585294e-03, 1.02052389e-03,
                   3.18030978e-02, 3.25357625e-04, 2.64397704e-03, 2.32129939e-03,
                   3.60466392e-04, 4.91762513e-04, 9.59802274e-04, 4.45126400e-02,
                   3.90329743e-04, 3.96902307e-04, 6.88101472e-04, 3.04961882e-03,
                   1.58361674e-03, 2.64605924e-04, 4.41523737e-01, 6.24756320e-04)

RF_Probs_Glia = c(0.16277056, 0.07293666, 0.13228492, 0.08610568, 0.12163417,
                  0.2527945 , 0.23992995, 0.1510989 , 0.17966696, 0.0406746 ,
                  0.15686275, 0.03819095, 0.018     , 0.07319392, 0.09794319,
                  0.15345528, 0.15018315, 0.05454545, 0.18342152, 0.22193211,
                  0.03781095, 0.19150943, 0.23871528, 0.03842459, 0.10547667,
                  0.15032081, 0.10210526, 0.33789954, 0.20481928, 0.10964083,
                  0.15447898, 0.2845462 , 0.05392157, 0.08105561, 0.00980392,
                  0.04403131, 0.30731707, 0.1798893 , 0.26644182, 0.22561493,
                  0.23492908, 0.24289643, 0.05852156, 0.22723404, 0.09541628,
                  0.1032197 , 0.09870389, 0.06565657, 0.44742063, 0.13459801,
                  0.02071006, 0.04      , 0.03434343, 0.0585975 , 0.14738806,
                  0.09496676, 0.13339731, 0.03853565, 0.11423221, 0.0805501 ,
                  0.1312559 , 0.01298701, 0.06568627, 0.41213202, 0.34869015,
                  0.09916589, 0.09803922, 0.04162331, 0.11464968, 0.24704025,
                  0.06634146, 0.05928086, 0.09961686, 0.03033268, 0.05096154,
                  0.04377432, 0.01516684, 0.36700767, 0.21020761, 0.15041021,
                  0.14951989, 0.01766438, 0.01771654, 0.2954023 , 0.21056978,
                  0.04116466, 0.06159769, 0.27288281, 0.01003009, 0.27884615,
                  0.0364532 , 0.27022654, 0.06863905, 0.03307393, 0.19748654,
                  0.12699906, 0.0489716 , 0.36052202, 0.06791339, 0.05705996,
                  0.2371741 , 0.02071563, 0.06087824, 0.13238095, 0.2026538 ,
                  0.57825752, 0.45825771, 0.302267  , 0.35239697, 0.10125362,
                  0.33661741, 0.32154341, 0.09123649, 0.2109777 , 0.14272971,
                  0.16666667, 0.07648184, 0.09439252, 0.0245821 , 0.29330254,
                  0.34528449, 0.07475728, 0.05836576, 0.03305785, 0.0742913 ,
                  0.16334283, 0.07049345, 0.39636076, 0.03195963, 0.0430622 ,
                  0.02865916, 0.08619001, 0.24017857, 0.0496592 , 0.05125628,
                  0.13660714, 0.02854331, 0.0186824 , 0.03703704, 0.04695305,
                  0.11173184, 0.16111111, 0.0610998 , 0.20990764, 0.01104418,
                  0.39984351, 0.01588878, 0.02058824, 0.12932605, 0.07574207,
                  0.05595117, 0.1025641 , 0.13830755, 0.2670112 , 0.05504587,
                  0.20054695, 0.23883162, 0.40948693, 0.4580292 , 0.21858407,
                  0.45962199, 0.29746282, 0.23657718, 0.15227934, 0.27830596,
                  0.24690339, 0.14860681, 0.21315789, 0.17186025, 0.0942029 ,
                  0.18518519, 0.26530612, 0.10598626, 0.10067764, 0.20465567,
                  0.0400782 , 0.09943715, 0.32037534, 0.16502947, 0.09425785,
                  0.28652139, 0.1107078 , 0.10684932, 0.12148594, 0.07212476,
                  0.07345739, 0.10857763, 0.20243674, 0.06230848, 0.01290963,
                  0.02237354, 0.07368421, 0.07372401, 0.03521878, 0.39583333,
                  0.05243089)

SVM_Probs_Glia = c(0.21912908, 0.09931967, 0.0496313 , 0.12969348, 0.18122786,
                   0.19308384, 0.1261687 , 0.11695307, 0.13367374, 0.09953449,
                   0.12535743, 0.0270476 , 0.0083905 , 0.06105242, 0.00824779,
                   0.10166527, 0.08479147, 0.00679712, 0.05994144, 0.06160949,
                   0.01564377, 0.10654695, 0.06265071, 0.00423097, 0.00594228,
                   0.10199616, 0.00656778, 0.23080106, 0.1206432 , 0.06140497,
                   0.16302212, 0.19891937, 0.03735073, 0.04242423, 0.01880193,
                   0.01047184, 0.22098828, 0.14322976, 0.10883977, 0.16429415,
                   0.05501473, 0.31239013, 0.00428569, 0.19128505, 0.09596011,
                   0.14021328, 0.04452702, 0.03039327, 0.42751073, 0.08974811,
                   0.01762825, 0.03989434, 0.00939162, 0.1316987 , 0.04486909,
                   0.02206657, 0.07716486, 0.03472745, 0.11692754, 0.01444405,
                   0.04871859, 0.01204273, 0.02030817, 0.48050029, 0.45400802,
                   0.09452424, 0.00349868, 0.00768133, 0.02522219, 0.09600513,
                   0.02906731, 0.08873455, 0.05985278, 0.03603877, 0.02297467,
                   0.01207672, 0.02597679, 0.38844833, 0.2755969 , 0.176225  ,
                   0.0712789 , 0.0056933 , 0.00867936, 0.18500926, 0.16982621,
                   0.03147915, 0.01848493, 0.35931185, 0.01751559, 0.31984525,
                   0.01915419, 0.16883051, 0.09489745, 0.04441085, 0.24416024,
                   0.07551758, 0.05836566, 0.2001479 , 0.06519385, 0.11149087,
                   0.1541122 , 0.0022848 , 0.00700351, 0.08584957, 0.10192888,
                   0.59611051, 0.4494023 , 0.37508522, 0.28222239, 0.04659499,
                   0.3083487 , 0.26573114, 0.01284591, 0.08630441, 0.07967092,
                   0.07602454, 0.03091961, 0.07937942, 0.02274165, 0.12657342,
                   0.30031789, 0.04617185, 0.04329545, 0.07024643, 0.04500144,
                   0.07390685, 0.01623373, 0.43639616, 0.00151226, 0.00943645,
                   0.00847237, 0.0277522 , 0.18230706, 0.02341225, 0.02251525,
                   0.17846358, 0.0253107 , 0.03296136, 0.01754913, 0.0347251 ,
                   0.16141416, 0.1253366 , 0.09725578, 0.17866191, 0.03695205,
                   0.3809077 , 0.03374162, 0.03409446, 0.03220558, 0.01935713,
                   0.01450819, 0.15033243, 0.21471869, 0.24982142, 0.00260835,
                   0.18706108, 0.29115739, 0.42202795, 0.32101539, 0.12162163,
                   0.44169825, 0.15905138, 0.13738341, 0.15170964, 0.10741649,
                   0.19274169, 0.05633034, 0.13687861, 0.10402799, 0.02293543,
                   0.27649717, 0.20369017, 0.05961949, 0.05082158, 0.1861483 ,
                   0.01491101, 0.03553347, 0.20443676, 0.066076  , 0.04979952,
                   0.24118182, 0.00387842, 0.10152046, 0.18526447, 0.01918511,
                   0.02418438, 0.07239732, 0.11668446, 0.03058396, 0.01554911,
                   0.02718844, 0.09576801, 0.11666454, 0.01138787, 0.35187521,
                   0.02731989)

TrueLabels_Glia = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                    0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0,
                    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                    0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0,
                    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0,
                    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0,
                    0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0)


#Pairwise concordance testing (U-statistic based)

#MLP vs KNN
rcorrp.cens(MLP_Probs_Glia,KNN_Probs_Glia,TrueLabels_Glia)
improveProb(MLP_Probs_Glia,KNN_Probs_Glia,TrueLabels_Glia)
#95% CI [-0.0623 0.0958]

#MLP vs RF
rcorrp.cens(MLP_Probs_Glia,RF_Probs_Glia,TrueLabels_Glia)
improveProb(MLP_Probs_Glia,RF_Probs_Glia,TrueLabels_Glia)
#95% CI [-0.0316 0.0960]

#MLP vs SVM
rcorrp.cens(MLP_Probs_Glia,SVM_Probs_Glia,TrueLabels_Glia)
improveProb(MLP_Probs_Glia,SVM_Probs_Glia,TrueLabels_Glia)
#95% CI [-0.0290 0.0836]

#RF vs KNN
rcorrp.cens(RF_Probs_Glia,KNN_Probs_Glia,TrueLabels_Glia)
improveProb(RF_Probs_Glia,KNN_Probs_Glia,TrueLabels_Glia)
#95% CI [-0.0878 0.0569]

#RF vs SVM
rcorrp.cens(RF_Probs_Glia,SVM_Probs_Glia,TrueLabels_Glia)
improveProb(RF_Probs_Glia,SVM_Probs_Glia,TrueLabels_Glia)
#95% CI [-0.02794 0.01815]

#SVM vs KNN
rcorrp.cens(SVM_Probs_Glia,KNN_Probs_Glia,TrueLabels_Glia)
improveProb(SVM_Probs_Glia,KNN_Probs_Glia,TrueLabels_Glia)
#95% CI [-0.0816 0.0605]

#All confidence intervals overlap with 0 - no model is clearly better than another


#Probabilities obtained from Python - Ox vs Rest

KNN_Probs_Ox = c(1.        , 1.        , 1.        , 1.        , 1.        ,
                 1.        , 1.        , 1.        , 1.        , 1.        ,
                 1.        , 1.        , 1.        , 1.        , 1.        ,
                 0.60452219, 1.        , 1.        , 1.        , 1.        ,
                 0.        , 1.        , 1.        , 0.        , 0.60004017,
                 1.        , 0.20074454, 1.        , 1.        , 0.        ,
                 1.        , 1.        , 1.        , 1.        , 1.        ,
                 1.        , 0.80066212, 1.        , 1.        , 1.        ,
                 0.80491735, 1.        , 1.        , 1.        , 1.        ,
                 1.        , 0.        , 1.        , 0.        , 1.        ,
                 1.        , 1.        , 1.        , 1.        , 1.        ,
                 1.        , 1.        , 1.        , 1.        , 1.        ,
                 1.        , 1.        , 0.        , 1.        , 0.        ,
                 1.        , 1.        , 1.        , 0.        , 0.        ,
                 1.        , 1.        , 1.        , 1.        , 1.        ,
                 1.        , 1.        , 0.80048596, 1.        , 1.        ,
                 0.80168254, 1.        , 1.        , 0.60446719, 1.        ,
                 1.        , 1.        , 1.        , 1.        , 1.        ,
                 1.        , 1.        , 0.61192086, 1.        , 0.80192134,
                 1.        , 1.        , 0.        , 1.        , 1.        ,
                 1.        , 0.40532451, 0.        , 1.        , 1.        ,
                 0.39998013, 0.        , 1.        , 0.        , 1.        ,
                 0.80883748, 1.        , 1.        , 1.        , 1.        ,
                 1.        , 1.        , 1.        , 1.        , 1.        ,
                 0.79611922, 1.        , 1.        , 1.        , 1.        ,
                 0.        , 1.        , 1.        , 0.61438238, 1.        ,
                 0.        , 1.        , 0.59594962, 1.        , 0.        ,
                 1.        , 1.        , 1.        , 1.        , 1.        ,
                 1.        , 1.        , 1.        , 1.        , 1.        ,
                 0.39856659, 1.        , 1.        , 1.        , 1.        ,
                 0.        , 0.        , 1.        , 1.        , 0.        ,
                 0.60422714, 1.        , 0.5937192 , 0.        , 1.        ,
                 0.        , 1.        , 1.        , 1.        , 1.        ,
                 1.        , 0.80613726, 1.        , 1.        , 0.        ,
                 0.8021571 , 1.        , 0.        , 1.        , 0.79899356,
                 0.        , 1.        , 0.19778559, 0.        , 0.        ,
                 1.        , 0.        , 1.        , 0.        , 1.        ,
                 0.19227589, 1.        , 0.6069899 , 0.        , 1.        ,
                 1.        , 1.        , 1.        , 1.        , 0.38521641,
                 0.80319117)

MLP_Probs_Ox = c(8.55165770e-01, 9.91235094e-01, 9.98996856e-01, 9.97227996e-01,
                 9.89170418e-01, 9.07758564e-01, 9.95185221e-01, 9.95647793e-01,
                 9.69301028e-01, 9.97204195e-01, 9.98429326e-01, 9.98660526e-01,
                 9.99322241e-01, 9.98998666e-01, 9.95660871e-01, 9.90802267e-01,
                 9.98169400e-01, 9.99445729e-01, 9.98432333e-01, 9.98223570e-01,
                 5.45968132e-04, 9.95546158e-01, 9.88788527e-01, 5.77520402e-04,
                 5.59949710e-01, 9.96063601e-01, 1.20364891e-01, 9.55529080e-01,
                 9.92991826e-01, 1.67362187e-04, 9.92375821e-01, 8.84697410e-01,
                 9.99127794e-01, 9.95098432e-01, 9.99251643e-01, 9.99649731e-01,
                 8.87462177e-01, 9.68535247e-01, 9.92606284e-01, 9.79778367e-01,
                 9.94208859e-01, 8.33940157e-01, 9.99233759e-01, 9.30965466e-01,
                 9.92028571e-01, 9.85429468e-01, 1.32310107e-04, 9.99117715e-01,
                 3.32095654e-03, 9.94612514e-01, 9.99495346e-01, 9.98503931e-01,
                 9.99308873e-01, 9.96608550e-01, 9.99221160e-01, 9.99164075e-01,
                 9.98172875e-01, 9.98888488e-01, 9.96438501e-01, 9.99556275e-01,
                 9.98787508e-01, 9.99470046e-01, 3.72923878e-04, 4.35216257e-01,
                 1.75476614e-04, 9.96080246e-01, 9.82933522e-01, 9.99159423e-01,
                 1.41375122e-03, 1.70965358e-04, 9.99483571e-01, 9.97351275e-01,
                 9.96328357e-01, 9.99197894e-01, 9.99067299e-01, 9.99643262e-01,
                 9.96612053e-01, 2.52385753e-01, 9.54944649e-01, 8.93538905e-01,
                 9.95298620e-01, 9.96511539e-01, 9.99335825e-01, 9.91041169e-01,
                 9.76549853e-01, 9.99165967e-01, 9.99338789e-01, 8.86187405e-01,
                 9.99560643e-01, 6.44627299e-01, 9.99593220e-01, 9.53203222e-01,
                 7.14474835e-01, 9.99314251e-01, 9.15671822e-01, 9.98062517e-01,
                 9.97744326e-01, 5.56325320e-03, 9.97903837e-01, 9.98423693e-01,
                 9.91152873e-01, 7.29874160e-01, 2.17456373e-03, 9.97930985e-01,
                 9.84817575e-01, 1.87740569e-01, 3.66028012e-03, 8.27960125e-01,
                 1.91568188e-04, 9.99301952e-01, 9.21257290e-01, 8.07955913e-01,
                 9.97522630e-01, 9.98120755e-01, 9.98210067e-01, 9.99056138e-01,
                 9.97687850e-01, 9.97245833e-01, 9.99288318e-01, 9.88646265e-01,
                 9.18210076e-01, 9.99167862e-01, 9.99294790e-01, 9.98118442e-01,
                 9.98651263e-01, 1.59793752e-04, 9.98909196e-01, 4.75068255e-01,
                 9.63510145e-01, 9.99444683e-01, 1.89474321e-03, 9.99136827e-01,
                 9.61782927e-01, 9.99441814e-01, 1.03918738e-03, 9.93467701e-01,
                 9.99517316e-01, 9.99358770e-01, 9.99418772e-01, 9.99010078e-01,
                 9.96964634e-01, 9.92596845e-01, 9.95612679e-01, 9.71901229e-01,
                 9.98988193e-01, 4.72495978e-01, 9.99114930e-01, 9.99334165e-01,
                 9.99378937e-01, 9.96471079e-01, 2.15488181e-03, 7.57218805e-05,
                 9.96707407e-01, 9.52806138e-01, 3.38645093e-03, 9.77368138e-01,
                 9.82153360e-01, 4.27223108e-01, 2.01407412e-03, 9.96129819e-01,
                 2.40252642e-04, 9.91606219e-01, 9.78894865e-01, 9.73347017e-01,
                 9.96719816e-01, 9.76770826e-01, 9.94844606e-01, 9.77089834e-01,
                 9.93067680e-01, 6.83853795e-04, 9.18680731e-01, 9.88910548e-01,
                 5.10645182e-04, 9.99237646e-01, 9.83147809e-01, 5.15452629e-03,
                 9.99371995e-01, 3.72234921e-02, 1.50076866e-04, 1.53197579e-02,
                 9.68111475e-01, 5.60051392e-02, 9.97239719e-01, 1.47565890e-03,
                 9.97735075e-01, 5.82753986e-03, 9.98641242e-01, 9.55027047e-01,
                 4.66245322e-03, 9.99508155e-01, 9.99124570e-01, 9.96871325e-01,
                 9.98313426e-01, 9.99051960e-01, 4.05715114e-01, 9.99019531e-01)

RF_Probs_Ox = c(0.81212121, 0.9059501 , 0.85661425, 0.87181996, 0.86350975,
                0.73860705, 0.75744308, 0.84249084, 0.81332165, 0.95238095,
                0.83600713, 0.95376884, 0.84      , 0.91444867, 0.70127326,
                0.53252033, 0.83882784, 0.91004785, 0.80864198, 0.76414273,
                0.0199005 , 0.8       , 0.71614583, 0.06051873, 0.42292089,
                0.81026581, 0.21684211, 0.61035008, 0.78571429, 0.0047259 ,
                0.81261426, 0.6901063 , 0.89901961, 0.70216777, 0.95980392,
                0.95009785, 0.68130081, 0.80904059, 0.71922428, 0.76929601,
                0.65868794, 0.73785518, 0.64681725, 0.76510638, 0.89522919,
                0.87689394, 0.04985045, 0.92222222, 0.08531746, 0.85546522,
                0.94477318, 0.945     , 0.93333333, 0.93756004, 0.8488806 ,
                0.86989554, 0.86276392, 0.91233141, 0.8829588 , 0.91060904,
                0.8611898 , 0.96403596, 0.01960784, 0.56110616, 0.01626016,
                0.8897127 , 0.65428277, 0.89073881, 0.11358811, 0.01657459,
                0.92487805, 0.93294461, 0.87739464, 0.95107632, 0.90192308,
                0.94552529, 0.83316481, 0.35294118, 0.78200692, 0.83865087,
                0.52812071, 0.84592738, 0.93503937, 0.33218391, 0.77291495,
                0.9246988 , 0.80173244, 0.72284003, 0.97993982, 0.6958042 ,
                0.95862069, 0.71925566, 0.45680473, 0.94552529, 0.63554758,
                0.86829727, 0.87855044, 0.06117455, 0.77362205, 0.9410058 ,
                0.72413793, 0.38606403, 0.0758483 , 0.84285714, 0.47527141,
                0.38781804, 0.03629764, 0.69101595, 0.03195963, 0.81388621,
                0.65188834, 0.67202572, 0.6362545 , 0.74785592, 0.85548617,
                0.82085561, 0.82982792, 0.89158879, 0.92232055, 0.6812933 ,
                0.52299299, 0.88932039, 0.94066148, 0.94731405, 0.9198436 ,
                0.01614435, 0.75327291, 0.58860759, 0.54331371, 0.82296651,
                0.02456499, 0.85210578, 0.72589286, 0.9406037 , 0.03517588,
                0.85267857, 0.96161417, 0.95968535, 0.90253411, 0.93706294,
                0.88175047, 0.83796296, 0.92158859, 0.76994123, 0.98493976,
                0.56103286, 0.95233366, 0.96078431, 0.86703097, 0.70522006,
                0.05493388, 0.0351377 , 0.85532302, 0.71920758, 0.07155963,
                0.72926162, 0.75515464, 0.29428848, 0.03284672, 0.71327434,
                0.05068729, 0.68853893, 0.73238255, 0.80795344, 0.68712187,
                0.7489678 , 0.63880289, 0.76666667, 0.82436261, 0.03312629,
                0.7967968 , 0.72734694, 0.04808636, 0.85866409, 0.77206596,
                0.06256109, 0.90056285, 0.13136729, 0.01375246, 0.04550379,
                0.7094431 , 0.09346642, 0.87762557, 0.03815261, 0.76413255,
                0.13124388, 0.86210641, 0.77975633, 0.04187947, 0.96921549,
                0.93093385, 0.91004785, 0.91304348, 0.90394877, 0.36354167,
                0.87511916)

SVM_Probs_Ox = c(0.76536108, 0.89542021, 0.94608501, 0.84697532, 0.8168111 ,
                 0.80466564, 0.87332856, 0.88053694, 0.86152813, 0.89757808,
                 0.87309595, 0.95956439, 0.93807905, 0.92237831, 0.82810413,
                 0.55495473, 0.91073456, 0.98729048, 0.93261055, 0.93191825,
                 0.21351413, 0.88275965, 0.92060094, 0.09174369, 0.47957103,
                 0.84048153, 0.684292  , 0.75128077, 0.87722158, 0.01205083,
                 0.819385  , 0.79649752, 0.91275279, 0.76527392, 0.95429742,
                 0.98299017, 0.77635261, 0.85469231, 0.88760771, 0.83248983,
                 0.90941305, 0.67325585, 0.91249045, 0.80711248, 0.89948481,
                 0.85511386, 0.04542395, 0.95597106, 0.08491992, 0.90120076,
                 0.94908016, 0.94908106, 0.90052436, 0.84771123, 0.94852007,
                 0.97038665, 0.9205778 , 0.91559818, 0.87824191, 0.97179462,
                 0.94522562, 0.94063129, 0.03992157, 0.40754562, 0.0366114 ,
                 0.88679321, 0.67073565, 0.89794688, 0.2806368 , 0.00203568,
                 0.96029599, 0.89337541, 0.91024013, 0.93989598, 0.91519307,
                 0.96879961, 0.84460607, 0.27631819, 0.71281577, 0.79161809,
                 0.60417215, 0.81903863, 0.9370113 , 0.65550679, 0.81305058,
                 0.95678843, 0.9411182 , 0.63952398, 0.96779412, 0.67680928,
                 0.97392125, 0.81888426, 0.55518965, 0.94353274, 0.59102914,
                 0.92195626, 0.90193814, 0.27031081, 0.79124593, 0.88718681,
                 0.83457103, 0.73066025, 0.19293161, 0.89814694, 0.86191529,
                 0.3892898 , 0.08903581, 0.62456799, 0.00926046, 0.93103594,
                 0.68983696, 0.73366079, 0.92486268, 0.89382642, 0.91845006,
                 0.92245691, 0.95803848, 0.91563937, 0.95217301, 0.86537556,
                 0.68409476, 0.94858393, 0.94034683, 0.90228801, 0.93806188,
                 0.00656368, 0.89770393, 0.53716919, 0.81419437, 0.93929544,
                 0.22808142, 0.90109108, 0.74230522, 0.9490094 , 0.14981082,
                 0.81243005, 0.96579344, 0.95091844, 0.95900995, 0.94531392,
                 0.83761519, 0.87289748, 0.86782135, 0.79982059, 0.95443171,
                 0.53320355, 0.94188378, 0.95829871, 0.96475346, 0.80716555,
                 0.10593361, 0.00283744, 0.78311151, 0.70586522, 0.19591197,
                 0.71502511, 0.70828231, 0.4814872 , 0.05029592, 0.86639479,
                 0.0052846 , 0.83352646, 0.86048559, 0.8280351 , 0.87741501,
                 0.80661053, 0.88164366, 0.85871927, 0.89120959, 0.11462291,
                 0.68792983, 0.79458579, 0.06701995, 0.94345126, 0.76228058,
                 0.18738305, 0.95156474, 0.49979658, 0.01357692, 0.27097442,
                 0.753122  , 0.36131437, 0.88972373, 0.19010904, 0.82193696,
                 0.31545466, 0.87275273, 0.82011712, 0.34112349, 0.96378917,
                 0.92792874, 0.89967008, 0.87626421, 0.91501695, 0.4300269 ,
                 0.90514835)

TrueLabels_Ox = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1,
                  1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1,
                  1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1,
                  0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1,
                  1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0,
                  1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                  0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1,
                  0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0,
                  1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1)


#Pairwise concordance testing (U-statistic based)

#MLP vs KNN
rcorrp.cens(MLP_Probs_Ox,KNN_Probs_Ox,TrueLabels_Ox)
improveProb(MLP_Probs_Ox,KNN_Probs_Ox,TrueLabels_Ox)
#95% CI [-0.0061 0.0574]

#MLP vs RF
rcorrp.cens(MLP_Probs_Ox,RF_Probs_Ox,TrueLabels_Ox)
improveProb(MLP_Probs_Ox,RF_Probs_Ox,TrueLabels_Ox)
#95% CI [-0.0822 -0.00976] * 

#MLP vs SVM
rcorrp.cens(MLP_Probs_Ox,SVM_Probs_Ox,TrueLabels_Ox)
improveProb(MLP_Probs_Ox,SVM_Probs_Ox,TrueLabels_Ox)
#95% CI [-0.127 -0.0555] *

#RF vs KNN
rcorrp.cens(RF_Probs_Ox,KNN_Probs_Ox,TrueLabels_Ox)
improveProb(RF_Probs_Ox,KNN_Probs_Ox,TrueLabels_Ox)
#95% CI [0.0376 0.106] *

#RF vs SVM
rcorrp.cens(RF_Probs_Ox,SVM_Probs_Ox,TrueLabels_Ox)
improveProb(RF_Probs_Ox,SVM_Probs_Ox,TrueLabels_Ox)
#95% CI [-0.07125 -0.01886] *

#SVM vs KNN
rcorrp.cens(SVM_Probs_Ox,KNN_Probs_Ox,TrueLabels_Ox)
improveProb(SVM_Probs_Ox,KNN_Probs_Ox,TrueLabels_Ox)
#95% CI [0.0783 0.155] *

#Probabilities obtained from Python - TD vs Rest

KNN_Probs_TD = c(0.        , 0.        , 0.        , 0.        , 0.        ,
                 0.        , 0.        , 0.        , 0.        , 0.        ,
                 0.        , 0.        , 0.        , 0.        , 0.        ,
                 0.39547781, 0.        , 0.        , 0.        , 0.        ,
                 1.        , 0.        , 0.        , 1.        , 0.39995983,
                 0.        , 0.79925546, 0.        , 0.        , 1.        ,
                 0.        , 0.        , 0.        , 0.        , 0.        ,
                 0.        , 0.        , 0.        , 0.        , 0.        ,
                 0.        , 0.        , 0.        , 0.        , 0.        ,
                 0.        , 1.        , 0.        , 0.59652557, 0.        ,
                 0.        , 0.        , 0.        , 0.        , 0.        ,
                 0.        , 0.        , 0.        , 0.        , 0.        ,
                 0.        , 0.        , 1.        , 0.        , 0.61423553,
                 0.        , 0.        , 0.        , 1.        , 1.        ,
                 0.        , 0.        , 0.        , 0.        , 0.        ,
                 0.        , 0.        , 0.19951404, 0.        , 0.        ,
                 0.19831746, 0.        , 0.        , 0.39553281, 0.        ,
                 0.        , 0.        , 0.        , 0.        , 0.        ,
                 0.        , 0.        , 0.        , 0.        , 0.19807866,
                 0.        , 0.        , 1.        , 0.        , 0.        ,
                 0.        , 0.59467549, 1.        , 0.        , 0.        ,
                 0.        , 1.        , 0.        , 1.        , 0.        ,
                 0.        , 0.        , 0.        , 0.        , 0.        ,
                 0.        , 0.        , 0.        , 0.        , 0.        ,
                 0.        , 0.        , 0.        , 0.        , 0.        ,
                 1.        , 0.        , 0.        , 0.38561762, 0.        ,
                 1.        , 0.        , 0.        , 0.        , 1.        ,
                 0.        , 0.        , 0.        , 0.        , 0.        ,
                 0.        , 0.        , 0.        , 0.        , 0.        ,
                 0.        , 0.        , 0.        , 0.        , 0.        ,
                 1.        , 1.        , 0.        , 0.        , 1.        ,
                 0.        , 0.        , 0.        , 1.        , 0.        ,
                 1.        , 0.        , 0.        , 0.        , 0.        ,
                 0.        , 0.19386274, 0.        , 0.        , 1.        ,
                 0.        , 0.        , 1.        , 0.        , 0.        ,
                 1.        , 0.        , 0.80221441, 1.        , 1.        ,
                 0.        , 1.        , 0.        , 1.        , 0.        ,
                 0.80772411, 0.        , 0.        , 1.        , 0.        ,
                 0.        , 0.        , 0.        , 0.        , 0.        ,
                 0.        )

MLP_Probs_TD = c(1.11322528e-04, 8.98664920e-05, 9.46831541e-05, 1.73574571e-04,
                 7.17159949e-05, 7.53259005e-05, 6.39259690e-05, 7.73024303e-05,
                 6.82651529e-05, 7.23979336e-05, 6.59849203e-05, 1.09753853e-04,
                 2.88725051e-04, 1.46305352e-04, 3.97297395e-03, 2.49993977e-03,
                 7.64398769e-05, 1.31810543e-04, 8.32975213e-05, 7.41695596e-05,
                 9.99179388e-01, 8.79241768e-05, 8.98414644e-05, 9.99133106e-01,
                 4.39786759e-01, 2.21785658e-04, 8.79001308e-01, 9.37080412e-05,
                 8.12388764e-05, 9.96529712e-01, 2.09149365e-04, 1.37156107e-04,
                 2.10284566e-04, 3.50117544e-03, 1.55916402e-04, 1.21607087e-04,
                 6.29765259e-05, 7.83749491e-05, 8.65745491e-05, 7.09580899e-05,
                 4.30220229e-04, 8.63208741e-05, 5.27835123e-04, 6.17717243e-05,
                 7.48085886e-05, 8.96702724e-05, 9.96766073e-01, 1.17122711e-04,
                 9.20990222e-01, 7.76397110e-05, 1.60478425e-04, 1.21023606e-04,
                 3.50001595e-04, 8.74901452e-05, 9.13811214e-05, 1.25929441e-04,
                 6.58737016e-05, 3.03931838e-04, 7.25049805e-05, 1.13780362e-04,
                 7.61002841e-05, 1.96905732e-04, 9.99192410e-01, 3.31830291e-04,
                 8.51383868e-01, 1.14988474e-04, 1.65530150e-02, 4.24403141e-04,
                 9.97213266e-01, 9.89193584e-01, 7.82389665e-05, 8.01325380e-05,
                 1.08389298e-04, 8.90412211e-05, 2.76947504e-04, 1.00346022e-04,
                 2.33372549e-03, 7.61086103e-02, 6.85337814e-05, 1.05527936e-04,
                 4.14421513e-03, 3.31228140e-03, 4.12020267e-04, 2.06511592e-03,
                 1.02514877e-04, 1.25724360e-04, 3.28025070e-04, 6.26208722e-05,
                 1.38629149e-04, 5.54161027e-05, 7.79549795e-05, 8.17524495e-05,
                 2.84044562e-01, 1.20377098e-04, 1.81476638e-02, 7.47012961e-05,
                 3.06941380e-04, 9.90909025e-01, 1.47772308e-03, 6.53473579e-05,
                 1.21351047e-04, 2.70003180e-01, 9.97577320e-01, 8.23427128e-05,
                 9.12512267e-04, 8.28720326e-05, 9.77749359e-01, 5.07601687e-05,
                 9.96934492e-01, 3.00803765e-04, 6.40455930e-05, 5.40023071e-05,
                 1.93240525e-03, 9.77237455e-05, 7.49307185e-05, 6.89592136e-05,
                 6.12318138e-04, 8.37147871e-05, 2.61398983e-04, 9.30632694e-05,
                 1.60638042e-04, 1.05136451e-04, 1.09840611e-04, 1.57250864e-04,
                 1.61067261e-04, 9.96281713e-01, 6.86657878e-04, 1.05153592e-04,
                 3.61523856e-02, 2.57948220e-04, 9.97828300e-01, 3.66476518e-04,
                 6.30296757e-04, 1.25039977e-04, 9.98506275e-01, 9.99281707e-05,
                 9.76601397e-05, 1.30995443e-04, 1.06254795e-04, 1.62440861e-04,
                 7.59503506e-05, 6.58476100e-05, 2.74148433e-04, 1.08560651e-04,
                 9.50679646e-05, 2.25127465e-04, 1.32841102e-04, 9.61076105e-05,
                 7.86053262e-05, 1.21852794e-03, 9.97480517e-01, 9.97213103e-01,
                 7.54114575e-05, 2.23162712e-04, 9.96376122e-01, 2.62926730e-04,
                 6.28755492e-05, 1.27889629e-03, 9.84832385e-01, 1.49036276e-04,
                 9.75401739e-01, 7.91441223e-05, 7.21830069e-05, 2.68032077e-04,
                 8.34452769e-05, 6.06921465e-05, 5.30080925e-04, 7.63199725e-05,
                 8.37432063e-05, 9.98853162e-01, 7.78604539e-04, 6.68385451e-05,
                 9.96858285e-01, 1.01464038e-04, 1.64828095e-04, 9.94327639e-01,
                 1.19540609e-04, 9.16365588e-01, 9.97364070e-01, 9.83659718e-01,
                 8.54268533e-05, 9.43669503e-01, 1.16304219e-04, 9.96203042e-01,
                 1.90445874e-03, 9.93680698e-01, 3.98956151e-04, 4.60312697e-04,
                 9.94947217e-01, 9.49422267e-05, 1.87328936e-04, 7.90562474e-05,
                 1.02956839e-04, 6.83434243e-04, 1.52761149e-01, 3.55712675e-04)

RF_Probs_TD = c(2.51082251e-02, 2.11132438e-02, 1.11008326e-02, 4.20743640e-02,
                1.48560817e-02, 8.59845228e-03, 2.62697023e-03, 6.41025641e-03,
                7.01139351e-03, 6.94444444e-03, 7.13012478e-03, 8.04020101e-03,
                1.42000000e-01, 1.23574144e-02, 2.00783546e-01, 3.14024390e-01,
                1.09890110e-02, 3.54066986e-02, 7.93650794e-03, 1.39251523e-02,
                9.42288557e-01, 8.49056604e-03, 4.51388889e-02, 9.01056676e-01,
                4.71602434e-01, 3.94133822e-02, 6.81052632e-01, 5.17503805e-02,
                9.46643718e-03, 8.85633270e-01, 3.29067642e-02, 2.53475061e-02,
                4.70588235e-02, 2.16776626e-01, 3.03921569e-02, 5.87084149e-03,
                1.13821138e-02, 1.10701107e-02, 1.43338954e-02, 5.08905852e-03,
                1.06382979e-01, 1.92483960e-02, 2.94661191e-01, 7.65957447e-03,
                9.35453695e-03, 1.98863636e-02, 8.51445663e-01, 1.21212121e-02,
                4.67261905e-01, 9.93676603e-03, 3.45167653e-02, 1.50000000e-02,
                3.23232323e-02, 3.84245917e-03, 3.73134328e-03, 3.51377018e-02,
                3.83877159e-03, 4.91329480e-02, 2.80898876e-03, 8.84086444e-03,
                7.55429651e-03, 2.29770230e-02, 9.14705882e-01, 2.67618198e-02,
                6.35049684e-01, 1.11214087e-02, 2.47678019e-01, 6.76378772e-02,
                7.71762208e-01, 7.36385162e-01, 8.78048780e-03, 7.77453839e-03,
                2.29885057e-02, 1.85909980e-02, 4.71153846e-02, 1.07003891e-02,
                1.51668352e-01, 2.80051151e-01, 7.78546713e-03, 1.09389243e-02,
                3.22359396e-01, 1.36408243e-01, 4.72440945e-02, 3.72413793e-01,
                1.65152766e-02, 3.41365462e-02, 1.36669875e-01, 4.27715997e-03,
                1.00300903e-02, 2.53496503e-02, 4.92610837e-03, 1.05177994e-02,
                4.74556213e-01, 2.14007782e-02, 1.66965889e-01, 4.70366886e-03,
                7.24779628e-02, 5.78303426e-01, 1.58464567e-01, 1.93423598e-03,
                3.86879731e-02, 5.93220339e-01, 8.63273453e-01, 2.47619048e-02,
                3.22074789e-01, 3.39244410e-02, 5.05444646e-01, 6.71704450e-03,
                6.15643398e-01, 8.48601736e-02, 1.14942529e-02, 6.43086817e-03,
                2.72509004e-01, 4.11663808e-02, 1.78412132e-03, 1.24777184e-02,
                9.36902486e-02, 1.40186916e-02, 5.30973451e-02, 2.54041570e-02,
                1.31722525e-01, 3.59223301e-02, 9.72762646e-04, 1.96280992e-02,
                5.86510264e-03, 8.20512821e-01, 1.76233635e-01, 1.50316456e-02,
                4.24726661e-01, 1.33971292e-01, 9.46775844e-01, 6.17042116e-02,
                3.39285714e-02, 9.73709834e-03, 9.13567839e-01, 1.07142857e-02,
                9.84251969e-03, 2.16322517e-02, 6.04288499e-02, 1.59840160e-02,
                6.51769088e-03, 9.25925926e-04, 1.73116090e-02, 2.01511335e-02,
                4.01606426e-03, 3.91236307e-02, 3.17775571e-02, 1.86274510e-02,
                3.64298725e-03, 2.19037871e-01, 8.89114954e-01, 8.62298196e-01,
                6.36942675e-03, 1.37812231e-02, 8.73394495e-01, 7.01914312e-02,
                6.01374570e-03, 2.96224589e-01, 5.09124088e-01, 6.81415929e-02,
                4.89690722e-01, 1.39982502e-02, 3.10402685e-02, 3.97672163e-02,
                3.45721694e-02, 4.12881916e-03, 2.12590299e-01, 2.01754386e-02,
                3.77714825e-03, 8.72670807e-01, 1.80180180e-02, 7.34693878e-03,
                8.45927380e-01, 4.06582769e-02, 2.32783705e-02, 8.97360704e-01,
                0.00000000e+00, 5.48257373e-01, 8.21218075e-01, 8.60238353e-01,
                4.03551251e-03, 7.95825771e-01, 1.55251142e-02, 8.40361446e-01,
                1.63742690e-01, 7.95298727e-01, 2.93159609e-02, 1.78069353e-02,
                8.95812053e-01, 1.78748759e-02, 4.66926070e-02, 1.62679426e-02,
                1.32325142e-02, 6.08324440e-02, 2.40625000e-01, 7.24499523e-02)

SVM_Probs_TD = c(1.55098384e-02, 5.26011902e-03, 4.28368495e-03, 2.33312075e-02,
                 1.96104085e-03, 2.25051462e-03, 5.02737123e-04, 2.50998479e-03,
                 4.79812943e-03, 2.88743355e-03, 1.54662623e-03, 1.33880100e-02,
                 5.35304559e-02, 1.65692733e-02, 1.63648083e-01, 3.43380005e-01,
                 4.47397357e-03, 5.91240088e-03, 7.44801216e-03, 6.47226194e-03,
                 7.70842100e-01, 1.06934015e-02, 1.67483482e-02, 9.04025340e-01,
                 5.14486694e-01, 5.75223173e-02, 3.09140222e-01, 1.79181743e-02,
                 2.13522097e-03, 9.26544199e-01, 1.75928844e-02, 4.58311144e-03,
                 4.98964781e-02, 1.92301853e-01, 2.69006511e-02, 6.53799177e-03,
                 2.65910673e-03, 2.07792715e-03, 3.55252042e-03, 3.21602445e-03,
                 3.55722216e-02, 1.43540239e-02, 8.32238613e-02, 1.60246654e-03,
                 4.55508345e-03, 4.67286910e-03, 9.10049032e-01, 1.36356704e-02,
                 4.87569344e-01, 9.05113108e-03, 3.32915886e-02, 1.10246022e-02,
                 9.00840222e-02, 2.05900677e-02, 6.61084703e-03, 7.54678096e-03,
                 2.25734083e-03, 4.96743647e-02, 4.83055644e-03, 1.37613345e-02,
                 6.05578842e-03, 4.73259820e-02, 9.39770263e-01, 1.11954087e-01,
                 5.09380583e-01, 1.86825488e-02, 3.25765669e-01, 9.43717823e-02,
                 6.94141010e-01, 9.01959184e-01, 1.06366971e-02, 1.78900384e-02,
                 2.99070894e-02, 2.40652468e-02, 6.18322561e-02, 1.91236675e-02,
                 1.29417145e-01, 3.35233474e-01, 1.15873363e-02, 3.21569106e-02,
                 3.24548951e-01, 1.75268075e-01, 5.43093401e-02, 1.59483950e-01,
                 1.71232157e-02, 1.17324238e-02, 4.03968727e-02, 1.16417023e-03,
                 1.46902854e-02, 3.34546725e-03, 6.92455594e-03, 1.22852326e-02,
                 3.49912904e-01, 1.20564023e-02, 1.64810627e-01, 2.52615968e-03,
                 3.96961943e-02, 5.29541294e-01, 1.43560216e-01, 1.32232117e-03,
                 1.13167718e-02, 2.67054952e-01, 8.00064881e-01, 1.60034907e-02,
                 3.61558328e-02, 1.45996963e-02, 4.61561891e-01, 3.46782619e-04,
                 7.08517154e-01, 2.23690687e-02, 1.81434160e-03, 6.08070227e-04,
                 6.22914132e-02, 1.98691728e-02, 1.87902016e-03, 1.51855104e-03,
                 1.10419036e-02, 4.98120695e-03, 2.50853445e-02, 8.05101846e-03,
                 1.55873472e-02, 5.24422791e-03, 1.63577128e-02, 2.74655612e-02,
                 1.69366840e-02, 9.19529461e-01, 8.60623333e-02, 2.64346520e-02,
                 1.84293378e-01, 5.12681048e-02, 7.63446219e-01, 7.11567208e-02,
                 7.53877244e-02, 2.75783495e-02, 8.27673929e-01, 9.10636462e-03,
                 8.89586163e-03, 1.61201976e-02, 2.34409247e-02, 1.99609858e-02,
                 9.70650869e-04, 1.76591821e-03, 3.49228716e-02, 2.15175009e-02,
                 8.61623908e-03, 8.58887558e-02, 2.43746078e-02, 7.60682333e-03,
                 3.04096121e-03, 1.73477328e-01, 8.79558206e-01, 8.46830124e-01,
                 2.16979612e-03, 4.43133688e-02, 8.01479683e-01, 9.79138126e-02,
                 5.60301864e-04, 9.64848521e-02, 6.28688686e-01, 1.19835808e-02,
                 5.53017144e-01, 7.42215872e-03, 2.13100123e-03, 2.02552572e-02,
                 1.51684966e-02, 6.47774827e-04, 6.20260051e-02, 4.40211530e-03,
                 4.76242133e-03, 8.62441655e-01, 3.55729922e-02, 1.72404854e-03,
                 8.73360556e-01, 5.72716858e-03, 5.15711120e-02, 7.97705934e-01,
                 1.29017938e-02, 2.95766661e-01, 9.20347085e-01, 6.79226060e-01,
                 5.69618363e-03, 6.34807212e-01, 8.75580564e-03, 6.24626485e-01,
                 1.58877930e-01, 6.60360961e-01, 5.48499439e-02, 6.31984204e-02,
                 6.28292551e-01, 2.06617189e-02, 4.48828161e-02, 4.56190986e-03,
                 7.07124813e-03, 7.35951777e-02, 2.18097890e-01, 6.75317550e-02)

TrueLabels_TD = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0,
                  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0,
                  0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0,
                  1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1,
                  0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1,
                  0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0,
                  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                  0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0)


#Pairwise concordance testing (U-statistic based)

#MLP vs KNN
rcorrp.cens(MLP_Probs_TD,KNN_Probs_TD,TrueLabels_TD)
improveProb(MLP_Probs_TD,KNN_Probs_TD,TrueLabels_TD)
#95% CI [-0.00975 0.03799]

#MLP vs RF
rcorrp.cens(MLP_Probs_TD,RF_Probs_TD,TrueLabels_TD)
improveProb(MLP_Probs_TD,RF_Probs_TD,TrueLabels_TD)
#95% CI [0.00613 0.09101] *

#MLP vs SVM
rcorrp.cens(MLP_Probs_TD,SVM_Probs_TD,TrueLabels_TD)
improveProb(MLP_Probs_TD,SVM_Probs_TD,TrueLabels_TD)
#95% CI [-0.04389 0.05242]

#RF vs KNN
rcorrp.cens(RF_Probs_TD,KNN_Probs_TD,TrueLabels_TD)
improveProb(RF_Probs_TD,KNN_Probs_TD,TrueLabels_TD)
#95% CI [-0.07089 0.00199]

#RF vs SVM
rcorrp.cens(RF_Probs_TD,SVM_Probs_TD,TrueLabels_TD)
improveProb(RF_Probs_TD,SVM_Probs_TD,TrueLabels_TD)
#95% CI [-0.07047 -0.01814] *

#SVM vs KNN
rcorrp.cens(SVM_Probs_TD,KNN_Probs_TD,TrueLabels_TD)
improveProb(SVM_Probs_TD,KNN_Probs_TD,TrueLabels_TD)
#95% CI [-0.0347 0.05440]


######################################################################################################################################################################################
#################################  QUESTION 9  #######################################################################################################################################
######################################################################################################################################################################################

#Figure S12: If these are strong candidates as markers for distinguishing FTLD and ALS/FTLD, consider supporting this idea by testing them as predictors in a validation data set.

load("D:/Jarrett/Research/Nat Comm Reviews/RData/ALSPatientStratification_UnivariateDatasets_RINSite_PeerReview.RData")

#ALS and FTLD Discriminatory features (from manuscript figures)
S12Features = c("AGER","AQP1","BECN2","C1D","CCDC154","ENSG00000260198","ENSG00000278434","ENSG00000279759","ENSG00000281969","HIST1H1T","HLA-DRB1","IFI30","IRF7","SELL","SERPINA1","SNX18P3","SOCS3","STH","TNRC6C-AS1","TUNAR")
F6Features = c("AIF1","APOC2","HLA-DRA","AGPAT4-IT1","CHKB-CPT1B","HSP90AB4P","MIR24-2")
F3Features =  c("FCGR3A")
S10Features = c("ADAT3","LINC00176","MIR219A2","SLX1B-SULT1A4","TUB-AS1")
ALS_FTLD_Features = c(S12Features,F6Features,F3Features,S10Features) 

table(colnames(FullCount) == FullPheno$Subject)

HCind = which(FullPheno$Subtype == "Control")

ALSFTLDCount = FullCount[,-HCind]
ALSFTLDPheno = FullPheno[-HCind,]

table(colnames(ALSFTLDCount) == ALSFTLDPheno$Subject)

Novaind = which(ALSFTLDPheno$platform == "NovaSeq")

ALSFTLDCount_Nova = ALSFTLDCount[,Novaind]
ALSFTLDPheno_Nova = ALSFTLDPheno[Novaind,]

table(colnames(ALSFTLDCount_Nova) == ALSFTLDPheno_Nova$Subject)

dim(ALSFTLDCount_Nova)
dim(ALSFTLDPheno_Nova)

## VST Normalization

rALSFTLD = round(ALSFTLDCount_Nova,0)
ddsALSFTLD = DESeqDataSetFromMatrix(countData = rALSFTLD, colData = ALSFTLDPheno_Nova, design= ~ SubjectSex, tidy=F)
ddsALSFTLD$SubjectSex = relevel(ddsALSFTLD$SubjectSex,ref = "Male")
dseqALSFTLD = DESeq(ddsALSFTLD,betaPrior=T)
resALSFTLD = results(dseqALSFTLD)
sigALSFTLD = resALSFTLD[! is.na(resALSFTLD$padj) & resALSFTLD$padj<0.05,]

vsdAF = varianceStabilizingTransformation(dseqALSFTLD)
vstcounts = assay(vsdAF)
ALSFTLD_VST_Nova = vstcounts[! (rownames(vstcounts) %in% rownames(sigALSFTLD)),]

## Classification Matrix

ALSFTLD_Classifier_Mat = ALSFTLD_VST_Nova[rownames(ALSFTLD_VST_Nova) %in% ALS_FTLD_Features,]

#setwd("C:/Users/jeshima/Desktop/ALS Manuscript/Nature/Communications Submission/Peer Review")
#write.csv(ALSFTLD_Classifier_Mat,"ALS_FTLD_ClassifierMatrix.csv")

#Clean up Pheno

for(i in 1:nrow(ALSFTLDPheno_Nova)){
  
  if(ALSFTLDPheno_Nova$Subtype[i] == "TE" || ALSFTLDPheno_Nova$Subtype[i] == "OX" || ALSFTLDPheno_Nova$Subtype[i] == "GLIA"){
    ALSFTLDPheno_Nova$Subtype[i] = "ALS"
  }
  
  if(ALSFTLDPheno_Nova$Subtype[i] == "OND"){
    ALSFTLDPheno_Nova$Subtype[i] = "FTLD"
  }
  
  if(ALSFTLDPheno_Nova$Factor[i] == 1 || ALSFTLDPheno_Nova$Factor[i] == 2 || ALSFTLDPheno_Nova$Factor[i] == 3){
    ALSFTLDPheno_Nova$Factor[i] = 1
  }else{
    ALSFTLDPheno_Nova$Factor[i] = 2
  }
  
}

#write.csv(ALSFTLDPheno_Nova,"ALS_FTLD_ClassifierPheno.csv")

#Classification performed in Python#

setwd("D:/Jarrett/Research/Nat Comm Reviews/FTLD Classification")

KNN_res = read.csv("Testing_ALS_FTLD_KNN_results_100CV_k5_PeerReview.csv")
MLP_res = read.csv("Testing_ALS_FTLD_MLP_results_3layer100Units_100CV_adam_tanh.csv")
RF_res = read.csv("Testing_ALS_FTLD_RF_results_ntree1000_100CV_PeerReview.csv")
SVM_res = read.csv("Testing_ALS_FTLD_SVM_results_100CV_Calibrated_PeerReview.csv")


list = c("KNN_res","MLP_res","RF_res","SVM_res")

for(j in 1:length(list)){
  
  tmp = get(list[j])
  nr = nrow(tmp)
  nc = ncol(tmp)
  
  for(i in 1:nr){
    
    if(tmp$X[i] == 0){
      tmp$X[i] = "Prediction"
    }else if(tmp$X[i] == 1){
      tmp$X[i] = "True.Subtype"
    }else if (tmp$X[i] == 2){
      tmp$X[i] = "Glia Prob"
    }else if(tmp$X[i] == 4){
      tmp$X[i] = "Ox Prob"
    }else if(tmp$X[i] == 3){
      tmp$X[i] = "TD Prob"
    }else{
      tmp$X[i] = "Subject ID"
    }
  }
  rownames(tmp) = make.unique(tmp$X)
  tmp = tmp[,-1]
  colnames(tmp) = paste("test.samp",seq(1:(nc-1)))
  assign(list[j],tmp)
  rm(tmp)
}


quickfun = function(x,colnm){
  return(x[[colnm]])
}

#Format: precision | recall | f1-score | support| Classifier
#ALS | FTLD

ngroup=2
nfold=100
CFlist = c("KNN","MLP","RF","SVM")
for(i in 1:length(list)){
  
  report = data.frame(matrix(NA,nrow=ngroup*nfold+1,ncol=6)) #no column names, header=F
  
  report[,6] = CFlist[i]
  report[,1] = c("",rep(c("ALS","FTLD"),nfold)) #may have to change the order manually depending on future reports
  report[1,] = c("","precision","recall","f1-score","support","Classifier")
  
  
  tmp = data.frame(t(get(list[i])))
  
  #Calculate Precision
  #TP/TP+FP = count of pred
  Pred = matrix(NA,nrow=nfold,ncol = nrow(tmp))
  TS = matrix(NA,nrow=nfold,ncol = nrow(tmp))
  
  for(j in 1:nfold){
    
    if(j == 1){
      handle = "Prediction"
      handle2 = "True.Subtype"
      Pred[j,] = quickfun(tmp,handle)
      TS[j,] = quickfun(tmp,handle2)
      
    }else{
      
      handle = paste("Prediction.",j-1,sep="")
      handle2 = paste("True.Subtype.",j-1,sep="")
      Pred[j,] = quickfun(tmp,handle)
      TS[j,] = quickfun(tmp,handle2)
      
    }
    alscount = FTLDcount=TPals=TPftld=0
    for(k in 1:nrow(tmp)){
      if(Pred[j,k] == "ALS"){
        alscount = alscount+1 #this is TP+FP
        
        if(Pred[j,k] == TS[j,k]){
          TPals = TPals+1
        }
        
        
      }else if(Pred[j,k] == "FTLD"){
        FTLDcount = FTLDcount+1
        
        if(Pred[j,k] == TS[j,k]){
          TPftld = TPftld+1
        }
        
      }
      
      
    }
    Precisionals = TPals/alscount
    Precisionftld = TPftld/FTLDcount
    
    start= (j*ngroup)
    stop = (j*ngroup)+1
    report$X2[start:stop] = c(Precisionals,Precisionftld)
    
  }
  
  #Calcualte Recall
  #TP/TP+FN
  #TP+FN = count of TrueSubtype (TS)
  Pred = matrix(NA,nrow=nfold,ncol = nrow(tmp))
  TS = matrix(NA,nrow=nfold,ncol = nrow(tmp))
  
  for(j in 1:nfold){
    
    if(j == 1){
      handle = "Prediction"
      handle2 = "True.Subtype"
      Pred[j,] = quickfun(tmp,handle)
      TS[j,] = quickfun(tmp,handle2)
      
    }else{
      
      handle = paste("Prediction.",j-1,sep="")
      handle2 = paste("True.Subtype.",j-1,sep="")
      Pred[j,] = quickfun(tmp,handle)
      TS[j,] = quickfun(tmp,handle2)
      
    }
    alscount = FTLDcount=TPals=TPftld=0
    for(k in 1:nrow(tmp)){
      if(TS[j,k] == "ALS"){
        alscount = alscount+1 #this is TP+FN
        
        if(Pred[j,k] == TS[j,k]){
          TPals = TPals+1
        }
        
        
      }else if(TS[j,k] == "FTLD"){
        FTLDcount = FTLDcount+1
        
        if(Pred[j,k] == TS[j,k]){
          TPftld = TPftld+1
        }
        
      }
    }
    Recallals = TPals/alscount
    Recallftld = TPftld/FTLDcount
    
    start= (j*ngroup)
    stop = (j*ngroup)+1
    report$X3[start:stop] = c(Recallals,Recallftld)
    
  }
  
  #Calculate F1-score
  for(m in 1:(nrow(report)-1)){
    
    pcn = as.numeric(report$X2[m+1])
    rcl = as.numeric(report$X3[m+1])
    
    report$X4[m+1] = (2*pcn*rcl)/(pcn+rcl)
  }
  
  
  
  #Write csv for report
  
  file = paste(list[i],"_report_",CFlist[i],".csv",sep="")
  write.table(report,file,sep = ",",col.names = F,row.names = F)
  
}

#All classifier have 100% accuracy in the testing (novaseq) cohort for all 100 rounds of CV
#No FTLD patients were characterized using HiSeq, so the holdout/validation cohort cannot be considered.


######################################################################################################################################################################################
######################################################################################################################################################################################
#################################  Reviewer 3  #######################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################


######################################################################################################################################################################################
#################################  QUESTION 1  #######################################################################################################################################
######################################################################################################################################################################################


#Region: The authors classify the regions used as "frontal" and "motor" cortices; what region(s) of the frontal lobe were used?
#Different frontal regions have different cellular compositions and are not equally affected in ALS (or in a cohort of FTLD patients used as controls).


setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping")
SubtypeLabels = read.csv("ALS451_coldata_SUBTYPES.csv")

for(i in 1:nrow(SubtypeLabels)){
  
  if(SubtypeLabels$tissue[i] == "Frontal Cortex"){
    SubtypeLabels$tissue[i] = "Cortex Frontal"
  }else if(SubtypeLabels$tissue[i] == "Lateral Motor Cortex"){
    SubtypeLabels$tissue[i] = "Cortex Motor Lateral"
  }else if(SubtypeLabels$tissue[i] == "Medial Motor Cortex"){
    SubtypeLabels$tissue[i] = "Cortex Motor Medial"
  }else if(SubtypeLabels$tissue[i] == "Other Motor Cortex"){
    SubtypeLabels$tissue[i] = "Cortex Motor Unspecified"
  }else{
    cat("Warning Flag")
  }
}

for(i in 1:nrow(SubtypeLabels)){
  if(SubtypeLabels$Subtype[i] == "TE"){
    SubtypeLabels$Subtype[i] = "TD"
  }
}

Tissuetable = table(tmpMetaData$tissue)

#[[1]] - Cortex Frontal
#[[2]] - Cortex Motor Lateral
#[[3]] - Cortex Motor Medial
#[[4]] - Cortex Motor Unspecified

#Containers
FrCor = rep(NA,Tissuetable[[1]])
LatCor = rep(NA,Tissuetable[[2]])
MedCor = rep(NA,Tissuetable[[3]])
Unsp = rep(NA,Tissuetable[[4]])

count1 = count2 = count3 = count4 = 1
for(i in 1:nrow(tmpMetaData)){
  
  for(j in 1:nrow(SubtypeLabels)){
    
    if(tmpMetaData$sample_id_alt[i] == SubtypeLabels$X[j] && SubtypeLabels$tissue[j] == "Cortex Frontal"){
      FrCor[count1] = SubtypeLabels$Subtype[j]
      count1 = count1+1
    }else if(tmpMetaData$sample_id_alt[i] == SubtypeLabels$X[j] && SubtypeLabels$tissue[j] == "Cortex Motor Lateral"){
      LatCor[count2] = SubtypeLabels$Subtype[j]
      count2 = count2+1
    }else if(tmpMetaData$sample_id_alt[i] == SubtypeLabels$X[j] && SubtypeLabels$tissue[j] == "Cortex Motor Medial"){
      MedCor[count3] = SubtypeLabels$Subtype[j]
      count3 = count3+1
    }else if(tmpMetaData$sample_id_alt[i] == SubtypeLabels$X[j] && SubtypeLabels$tissue[j] == "Cortex Motor Unspecified"){
      Unsp[count4] = SubtypeLabels$Subtype[j] 
      count4 = count4+1
    }
  }
}

#Frontal Cortex
par(mar=c(5,5,4,2))
tmp = table(FrCor)
barplot(c(tmp[[1]],tmp[[2]],tmp[[3]]),ylim=c(0,120),main = "Frontal Cortex Subtypes",xlab="Subtype",col = c("goldenrod1","navyblue","firebrick"),cex=2,cex.main=2,cex.axis=2,cex.lab=2)
axis(1,c(0.65,1.9,3.1),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=2)
title(ylab = "Frequency",mgp = c(3.5,1,0),cex.lab=2)
text(0.65,tmp[[1]]+0.15*tmp[[1]],labels = tmp[[1]],cex=2)
text(1.9,tmp[[2]]+0.05*tmp[[2]],labels = tmp[[2]],cex=2)
text(3.1,tmp[[3]]+0.1*tmp[[3]],labels = tmp[[3]],cex=2)

#Lateral Motor Cortex
tmp = table(LatCor)
barplot(c(tmp[[1]],tmp[[2]],tmp[[3]]),ylim=c(0,80),main = "Lateral Motor Cortex Subtypes",xlab="Subtype",col = c("goldenrod1","navyblue","firebrick"),cex=2,cex.main=2,cex.axis=2,cex.lab=2)
axis(1,c(0.65,1.9,3.1),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=2)
title(ylab = "Frequency",mgp = c(3.5,1,0),cex.lab=2)
text(0.65,tmp[[1]]+0.2*tmp[[1]],labels = tmp[[1]],cex=2)
text(1.9,tmp[[2]]+0.05*tmp[[2]],labels = tmp[[2]],cex=2)
text(3.1,tmp[[3]]+0.1*tmp[[3]],labels = tmp[[3]],cex=2)

#Medial Motor Cortex
tmp = table(MedCor)
barplot(c(tmp[[1]],tmp[[2]],tmp[[3]]),ylim=c(0,80),main = "Medial Motor Cortex Subtypes",xlab="Subtype",col = c("goldenrod1","navyblue","firebrick"),cex=2,cex.main=2,cex.axis=2,cex.lab=2)
axis(1,c(0.65,1.9,3.1),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=2)
title(ylab = "Frequency",mgp = c(3.5,1,0),cex.lab=2)
text(0.65,tmp[[1]]+0.2*tmp[[1]],labels = tmp[[1]],cex=2)
text(1.9,tmp[[2]]+0.05*tmp[[2]],labels = tmp[[2]],cex=2)
text(3.1,tmp[[3]]+0.1*tmp[[3]],labels = tmp[[3]],cex=2)

#Unspecified Motor Cortex
tmp = table(Unsp)
barplot(c(tmp[[1]],tmp[[2]],tmp[[3]]),ylim=c(0,40),main = "Unspecified Motor Cortex Subtypes",xlab="Subtype",col = c("goldenrod1","navyblue","firebrick"),cex=2,cex.main=2,cex.axis=2,cex.lab=2)
axis(1,c(0.65,1.9,3.1),labels = c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=2)
title(ylab = "Frequency",mgp = c(3.5,1,0),cex.lab=2)
text(0.65,tmp[[1]]+0.15*tmp[[1]],labels = tmp[[1]],cex=2)
text(1.9,tmp[[2]]+0.15*tmp[[2]],labels = tmp[[2]],cex=2)
text(3.1,tmp[[3]]+0.15*tmp[[3]],labels = tmp[[3]],cex=2)


#Cell deconvolution was performed using CIBERSORT with reference expression from Nowakowski et al. 


######################################################################################################################################################################################
#################################  QUESTION 2  #######################################################################################################################################
######################################################################################################################################################################################

#Although data from different brain regions are available from this cohort, the authors used only frontal/motor. 
#Including control regions that are not affected in ALS such as the cerebellum is helpful to address covariates 
#such as the RIN, the postmortem interval, agonal stage, and genetic factors

#Reason: Cerebellum samples are not available for all 208 patients considered in this analysis (150/208)

#Evidence
setwd("C:/Users/jeshima/Desktop/ALS Manuscript/Nature/Communications Submission/Peer Review")
SurvivalDat = read.csv("Table_S10.csv")
Subjects = SurvivalDat$SubjectID

Meta = read.csv("CLINICAL_DATA_PRUDENCIO.csv") #Available by request from NYGC
NewMeta = Meta[Meta$ExternalSubjectId %in% Subjects,]

table(NewMeta$Sample.Source)


## Genetic factors

library(stringr)

setwd("C:/Users/jeshima/Desktop/ALS Manuscript/Nature/Communications Submission/Peer Review")
MetaData = read.csv("Table_S10.csv")

#Reference values
nsubtype = as.numeric(length(table(MetaData$FinalSubtype)))
nmutations = 2 #c9orf72 and sod1
C9_char = as.numeric(nchar("C9orf72"))
SOD_char = as.numeric(nchar("SOD1"))
MySubtypes = names(table(MetaData$FinalSubtype))
nGlia = as.numeric(length(which(MetaData$FinalSubtype == "GLIA")))
nOx = as.numeric(length(which(MetaData$FinalSubtype == "OX")))
nTD = as.numeric(length(which(MetaData$FinalSubtype == "TE")))
nDisc = as.numeric(length(which(MetaData$FinalSubtype == "Discordant")))
npatients = as.numeric(nrow(MetaData))

#use string searching method 

for(i in 1:nsubtype){
  
  currentsub = MySubtypes[i]
  
  #Build container
  tmpmat = data.frame(matrix(0,nmutations,nmutations+1))
  colnames(tmpmat) = c("C9orf72","SOD1","Unknown")
  rownames(tmpmat) = c("Positive","Negative")
  
  
  for(j in 1:npatients){
    
    if(MetaData$FinalSubtype[j] == MySubtypes[i]){
      
      muts = str_split(MetaData$Mutation[j],",") #Manually converted ; to , in excel
      nmuts = length(muts[[1]])
      
      for(k in 1:nmuts){
        
        charlen = nchar(muts[[1]][k])
        tmpstring = muts[[1]][k]
        
        if(muts[[1]][k] == "Unknown"){
          tmpmat[1,3] = tmpmat[1,3]+1
        }
        
        #C9orf72
        for(m in 1:charlen-C9_char+1){
          
          if(str_sub(tmpstring,m,m+C9_char-1) == "C9orf72"){
            
            for(p in 1:charlen){
              if(str_sub(tmpstring,p,p+2) == "Pos" || str_sub(tmpstring,p,p+2) =="pos"){ #Ignores case sensitivity
                tmpmat[1,1] = tmpmat[1,1]+1
              }else if(str_sub(tmpstring,p,p+2) == "Neg" || str_sub(tmpstring,p,p+2) =="neg"){
                tmpmat[2,1] = tmpmat[2,1]+1
              }
            }
            
          }
          
        }
        
        #SOD1
        for(n in 1:charlen-SOD_char+1){
          
          if(str_sub(tmpstring,n,n+SOD_char-1) == "SOD1"){
            
            for(q in 1:charlen){
              
              if(str_sub(tmpstring,q,q+2) == "Pos" || str_sub(tmpstring,q,q+2) == "pos"){
                tmpmat[1,2] = tmpmat[1,2]+1
              }else if(str_sub(tmpstring,q,q+2) == "Neg" || str_sub(tmpstring,q,q+2) =="neg"){
                tmpmat[2,2] = tmpmat[2,2]+1
              }
              
            }
            
          }
          
        }
        
      }
      
    }
    
    
  }
  #Return results
  tmpname = paste(MySubtypes[i],"Mutations",sep="")
  assign(tmpname,tmpmat)
}


#Plot

CleanLabels = c("Discordant","ALS-Glia","ALS-Ox","ALS-TD")

#Generate plotting data frames
plotdatpos = data.frame(matrix(NA,3,4))
colnames(plotdatpos) = MySubtypes
rownames(plotdatpos) = colnames(GLIAMutations)
plotdatpos[,1] = as.numeric(DiscordantMutations[1,])
plotdatpos[,2] = as.numeric(GLIAMutations[1,])
plotdatpos[,3] = as.numeric(OXMutations[1,])
plotdatpos[,4] = as.numeric(TEMutations[1,])

plotdatneg = data.frame(matrix(NA,2,4))
colnames(plotdatneg) = MySubtypes
rownames(plotdatneg) = colnames(GLIAMutations)[1:2]
plotdatneg[,1] = as.numeric(DiscordantMutations[2,1:2])
plotdatneg[,2] = as.numeric(GLIAMutations[2,1:2])
plotdatneg[,3] = as.numeric(OXMutations[2,1:2])
plotdatneg[,4] = as.numeric(TEMutations[2,1:2])

plotdatpos = as.matrix(plotdatpos)
colnames(plotdatpos) = CleanLabels
rownames(plotdatpos) = colnames(GLIAMutations)
plotdatneg = as.matrix(plotdatneg)
colnames(plotdatneg) = CleanLabels
rownames(plotdatneg) = colnames(GLIAMutations)[1:2]

#Positive Mutations
barplot(plotdatpos,beside=T,col = c("palegreen2","plum2","ivory3"),main = "Positive Mutations by Subtype",ylab = "Frequency",ylim=c(0,35))
legend(1,35,legend = c("C9orf72","SOD1","Unknown"),col = c("palegreen2","plum2","ivory3"),pch=15,pt.cex=1.5)
abline(h=0)

#Negative Mutations
barplot(plotdatneg,beside=T,col = c("palegreen2","plum2"),main = "Negative Mutations by Subtype",ylab = "Frequency",ylim=c(0,50))
legend(1,50,legend = c("C9orf72","SOD1"),col = c("palegreen2","plum2"),pch=15,pt.cex=1.5)
abline(h=0)

#Discordant Mutations
DiscordantMutations = as.matrix(DiscordantMutations)
customcol = c(rep(c("palegreen2","plum2"),nmutations),"gray50")
barplot(DiscordantMutations,beside=T,col = customcol,main = "Discordant Mutations",ylab = "Frequency",ylim=c(0,20),xaxt = "n")
axis(1,c(2,5,7.5),labels = c("C9orf72","SOD1","Unknown"))
legend(7,20,legend = c("Positive","Negative","Unknown"),col = c("palegreen2","plum2","gray50"),pch=15,pt.cex=1.5)
abline(h=0)


#GLIA Mutations
GLIAMutations = as.matrix(GLIAMutations)
customcol = c(rep(c("palegreen2","plum2"),nmutations),"gray50")
barplot(GLIAMutations,beside=T,col = customcol,main = "ALS-Glia Mutations",ylab = "Frequency",ylim=c(0,30),xaxt = "n")
axis(1,c(2,5,7.5),labels = c("C9orf72","SOD1","Unknown"))
legend(7,30,legend = c("Positive","Negative","Unknown"),col = c("palegreen2","plum2","gray50"),pch=15,pt.cex=1.5)
abline(h=0)


#OX Mutations
OXMutations = as.matrix(OXMutations)
customcol = c(rep(c("palegreen2","plum2"),nmutations),"gray50")
barplot(OXMutations,beside=T,col = customcol,main = "ALS-Ox Mutations",ylab = "Frequency",ylim=c(0,50),xaxt = "n")
axis(1,c(2,5,7.5),labels = c("C9orf72","SOD1","Unknown"))
legend(7,50,legend = c("Positive","Negative","Unknown"),col = c("palegreen2","plum2","gray50"),pch=15,pt.cex=1.5)
abline(h=0)


#TD Mutations
TEMutations = as.matrix(TEMutations)
customcol = c(rep(c("palegreen2","plum2"),nmutations),"gray50")
barplot(TEMutations,beside=T,col = customcol,main = "ALS-TD Mutations",ylab = "Frequency",ylim=c(0,30),xaxt = "n")
axis(1,c(2,5,7.5),labels = c("C9orf72","SOD1","Unknown"))
legend(1,30,legend = c("Positive","Negative","Unknown"),col = c("palegreen2","plum2","gray50"),pch=15,pt.cex=1.5)
abline(h=0)


#Manually clean the matrices
GLIAMutations2 = GLIAMutations
tmp = rep(NA,ncol(GLIAMutations))
GLIAMutations2 = rbind(GLIAMutations2,tmp)
rownames(GLIAMutations2) = c("Positive","Negative","Unknown")
GLIAMutations2[3,3] = GLIAMutations2[1,3]
GLIAMutations2[1,3] = NA
GLIAMutations2[2,3] = NA

DiscordantMutations2 = DiscordantMutations
tmp = rep(NA,ncol(DiscordantMutations))
DiscordantMutations2 = rbind(DiscordantMutations2,tmp)
rownames(DiscordantMutations2) = c("Positive","Negative","Unknown")
DiscordantMutations2[3,3] = DiscordantMutations2[1,3]
DiscordantMutations2[1,3] = NA
DiscordantMutations2[2,3] = NA

OXMutations2 = OXMutations
tmp = rep(NA,ncol(OXMutations))
OXMutations2 = rbind(OXMutations2,tmp)
rownames(OXMutations2) = c("Positive","Negative","Unknown")
OXMutations2[3,3] = OXMutations2[1,3]
OXMutations2[1,3] = NA
OXMutations2[2,3] = NA

TEMutations2 = TEMutations
tmp = rep(NA,ncol(TEMutations))
TEMutations2 = rbind(TEMutations2,tmp)
rownames(TEMutations2) = c("Positive","Negative","Unknown")
TEMutations2[3,3] = TEMutations2[1,3]
TEMutations2[1,3] = NA
TEMutations2[2,3] = NA



#Grouped and stacked bar chart

#Create data frame

FullBar = data.frame(matrix(NA,20,4))
colnames(FullBar) = c("Facet","Group","Mutation","Value")
FullBar$Facet = c(rep("Discordant",5),rep("GLIA",5),rep("OX",5),rep("TE",5))
FullBar$Group = rep(c(rep(c("Positive","Negative"),2),"Unknown"),nsubtype)
FullBar$Mutation = rep(c(rep("C9orf72",2),rep("SOD1",2),"Unknown"),nsubtype)

#Fill data frame
for(i in 1:nsubtype){
  
  tmpdat = get(paste(MySubtypes[i],"Mutations2",sep=""))
  
  for(j in 1:nrow(tmpdat)){
    for(k in 1:ncol(tmpdat)){
      for(l in 1:nrow(FullBar)){
        if(FullBar$Facet[l] == MySubtypes[i]){
          if(rownames(tmpdat)[j] == FullBar$Group[l] && colnames(tmpdat)[k] == FullBar$Mutation[l]){
            FullBar$Value[l] = tmpdat[j,k]
          }
          
        }
      }
      
    }
    
  }
  
}

FullBar$Facet = c(rep("Discordant",5),rep("ALS-Glia",5),rep("ALS-Ox",5),rep("ALS-TD",5))

#mycol = rep(c(rep("lightsalmon2",2),rep("palegreen3",2),"thistle3"),nsubtype)

p = ggplot(FullBar,aes(x=Group,y=Value,fill=Mutation)) + geom_bar(stat = "identity",position = "stack") + facet_grid(~Facet)
p = p+scale_fill_manual(values = c("lightsalmon2","palegreen3","thistle3"))
p = p +ggtitle("C9orf72 and SOD1 Mutations by Subtype") + xlab("") + ylab("Frequency")
p = p+theme(axis.text = element_text(size=12), axis.title = element_text(size=14),plot.title = element_text(size=22))
p = p+theme(axis.title.y=element_text(angle=90, vjust=2,size = 18))
p = p+theme(plot.title = element_text(hjust = 0.5))
p = p+theme(text = element_text(size=14))
p = p+theme(plot.margin = margin(t=10,r=10,b=10,l=10))
p = p+theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5,linetype = "solid"),panel.grid.major = element_line(size = 0.35,linetype = "solid",colour = "gray80"),panel.grid.minor = element_line(size = 0.15,linetype = "solid",colour = "gray80"))
p = p+theme(axis.text.x = element_text(size = 13))
p = p+theme(axis.text.y = element_text(size = 16))
p = p+theme(legend.text = element_text(size = 16))
p = p+theme(legend.title = element_text(size = 16))
p


#################chi-squared stats

#Contingency Tables (mutation-wise)

CSmat = data.frame(matrix(NA,nsubtype,nmutations+2))
colnames(CSmat) = c("Subtype","Positive_C9orf72","Negative_C9orf72","Unknown")
CSmat$Subtype = c("Discordant","ALS-Glia","ALS-Ox","ALS-TD")

CSmatSOD = data.frame(matrix(NA,nsubtype,nmutations+2))
colnames(CSmatSOD) = c("Subtype","Positive_SOD1","Negative_SOD1","Unknown")
CSmatSOD$Subtype = c("Discordant","ALS-Glia","ALS-Ox","ALS-TD")

for(i in 1:nrow(CSmat)){
  
  for(j in 1:nrow(FullBar)){
    
    if(CSmat$Subtype[i] == FullBar$Facet[j] && FullBar$Group[j] == "Positive" && FullBar$Mutation[j] == "C9orf72"){
      CSmat$Positive_C9orf72[i] = FullBar$Value[j]
    }else if(CSmat$Subtype[i] == FullBar$Facet[j] && FullBar$Group[j] == "Negative" && FullBar$Mutation[j] == "C9orf72"){
      CSmat$Negative_C9orf72[i] = FullBar$Value[j]
    }else if(CSmat$Subtype[i] == FullBar$Facet[j] && FullBar$Group[j] == "Unknown" && FullBar$Mutation[j] == "Unknown"){
      CSmat$Unknown[i] = FullBar$Value[j]
    }
    
    
    if(CSmatSOD$Subtype[i] == FullBar$Facet[j] && FullBar$Group[j] == "Positive" && FullBar$Mutation[j] == "SOD1"){
      CSmatSOD$Positive_SOD1[i] = FullBar$Value[j]
    }else if(CSmatSOD$Subtype[i] == FullBar$Facet[j] && FullBar$Group[j] == "Negative" && FullBar$Mutation[j] == "SOD1"){
      CSmatSOD$Negative_SOD1[i] = FullBar$Value[j]
    }else if(CSmatSOD$Subtype[i] == FullBar$Facet[j] && FullBar$Group[j] == "Unknown" && FullBar$Mutation[j] == "Unknown"){
      CSmatSOD$Unknown[i] = FullBar$Value[j]
    }
    
  }
  
}

#Clean up
rownames(CSmat) = CSmat[,1]
CSmat = CSmat[,-1]
rownames(CSmatSOD) = CSmatSOD[,1]
CSmatSOD = CSmatSOD[,-1]

#With unknown category
chisq.test(CSmat) #p = 0.006119 (significance driven by unknown category, see results below)
chisq.test(CSmatSOD) #p = 0.483

#Without unknown category
CSmat2 = CSmat[-which(colnames(CSmat) == "Unknown")]
CSmatSOD2 = CSmatSOD[-which(colnames(CSmatSOD) == "Unknown")]

chisq.test(CSmat2) #p = 0.4726
chisq.test(CSmatSOD2) #p = 0.2148


######################################################################################################################################################################################
#################################  QUESTION 3  #######################################################################################################################################
######################################################################################################################################################################################

#Cellular composition is usually a driver of clustering in the analysis of bulk RNAseq, and so clusters unbiasedly identified may reflect the relative abundant of specific cell types. 
#Cell type deconvolution methods to infer cell type proportions may help address this. This is important in this paper since two of the subtypes of ALS identified correspond with glial (ALS-Glia) and neuronal (ALS-Ox) responses.

#Cell deconvolution was performed using CIBERSORT with reference expression from Nowakowski et al. 

######################################################################################################################################################################################
#################################  QUESTION 4  #######################################################################################################################################
######################################################################################################################################################################################

#Differential gene expression needs to take into account potential covariates beyond sex. 

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping")
Meta = read.csv("GSE153960_MetaData.txt") #Publicly available at: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA644618&o=acc_s%3Aa
Clinical = read.csv("CLINICAL_DATA_PRUDENCIO.csv")

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Pub Data/Supplemental Files")
SRRs = read.csv("SRR_IDs.csv") #Table S2
colnames(SRRs) = SRRs[1,]
SRRs = SRRs[-1,]

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

#Add RIN values to Supplemental Table
EshimaMetaData$RIN = rep(NA,nrow(EshimaMetaData))

for(i in 1:nrow(EshimaMetaData)){
  
  for(j in 1:nrow(Clinical)){
    
    if(EshimaMetaData$sample_id_alt[i] == Clinical$ExternalSampleId[j]){
      
      EshimaMetaData$RIN[i] = Clinical$RIN[j] 
      
    }
    
  }
  
}

ALS.RIN.Avg = mean(EshimaMetaData$RIN,na.rm=T)

which(is.na(EshimaMetaData$RIN)) #one missing value from Target ALS / NYGC


#RIN parsed by subtype
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping")
SubtypeLabels = read.csv("ALS451_coldata_SUBTYPES.csv")

convertnames = gsub("-","\\.",EshimaMetaData$sample_id_alt)
tmpMetaData = EshimaMetaData
tmpMetaData$sample_id_alt = convertnames

nsubtype = table(SubtypeLabels$Subtype)

Glia.RIN = rep(NA,nsubtype[[1]])
Ox.RIN = rep(NA,nsubtype[[2]])
TD.RIN = rep(NA,nsubtype[[3]])

count1 = count2 = count3 = 1

for(i in 1: nrow(tmpMetaData)){
  
  for(j in 1:nrow(SubtypeLabels)){
    
    if(tmpMetaData$sample_id_alt[i] == SubtypeLabels$X[j] && SubtypeLabels$Subtype[j] == "GLIA"){
      Glia.RIN[count1] = tmpMetaData$RIN[i]
      count1 = count1 + 1
    }else if(tmpMetaData$sample_id_alt[i] == SubtypeLabels$X[j] && SubtypeLabels$Subtype[j] == "OX"){
      Ox.RIN[count2] = tmpMetaData$RIN[i]
      count2 = count2 + 1
    }else if(tmpMetaData$sample_id_alt[i] == SubtypeLabels$X[j] && SubtypeLabels$Subtype[j] == "TE"){
      TD.RIN[count3] = tmpMetaData$RIN[i]
      count3 = count3 + 1
    }
    
  }
}

Glia.RIN.Avg = mean(Glia.RIN,na.rm = T)
Ox.RIN.Avg = mean(Ox.RIN,na.rm = T)
TD.RIN.Avg = mean(TD.RIN, na.rm = T)

par(mfrow=c(1,3))
#Normal Distribution isn't a terrible approximation...
hist(Glia.RIN,main="Glia RIN",breaks = 20);hist(Ox.RIN,main="Ox RIN",breaks = 20);hist(TD.RIN,main="TD RIN",breaks = 20)
par(mfrow=c(1,1))

#Stats
GO.RIN.sig = t.test(Glia.RIN,Ox.RIN) #p = 0.44
GO.RIN.sig = GO.RIN.sig$p.value
GT.RIN.sig = t.test(Glia.RIN,TD.RIN) #p = 6E-8
GT.RIN.sig = GT.RIN.sig$p.value
OT.RIN.sig = t.test(Ox.RIN,TD.RIN) #p = 7E-16
OT.RIN.sig = OT.RIN.sig$p.value

#Summary Boxplot
boxplot(Glia.RIN,Ox.RIN,TD.RIN,main="Subtype RIN",ylab="RIN",ylim = c(0,12),cex.lab = 1.5,cex.main=1.5,cex.axis=1.5,col=c("goldenrod1","navy","firebrick"),xaxt="n")
axis(at=1:3,side=1,labels=c("ALS-Glia","ALS-Ox","ALS-TD"),cex.axis=1.5)
lines(c(1,2),c(10,10),type="l")
lines(c(1,1),c(10,9.5),type="l")
lines(c(2,2),c(10,9.5),type="l")
text(mean(c(1,2)),10.3,labels = paste("p-value:",formatC(GO.RIN.sig,format = "e",digits = 2)))
lines(c(2,3),c(10.75,10.75),type="l")
lines(c(2,2),c(10.75,10.25),type="l")
lines(c(3,3),c(10.75,10.25),type="l")
text(mean(c(2,3)),11.05,labels = paste("p-value:",formatC(OT.RIN.sig,format = "e",digits = 2)))
lines(c(1,3),c(11.5,11.5),type="l")
lines(c(1,1),c(11.5,11),type="l")
lines(c(3,3),c(11.5,11),type="l")
text(mean(c(1,3)),11.8,labels = paste("p-value:",formatC(GT.RIN.sig,format = "e",digits = 2)))

#####################################################################################################################

#Univariate analysis with RIN as covariate (design = ~platform + RIN + Subtype )

#Design Equation 1: design = ~ platform + Subtype (from UnivariateAnalysis script)
#Design Equation 2: design = ~ platform + RIN + Subtype

load("D:/Jarrett/Research/Fall 2021/ProgBM/SOD1/ALSPatientStratification_UnivariateDatasets_SOD1.RData")

#Filter for subtype-specific features in Fig. 6
my36 = c("AIF1","APOC2","CD44","CHI3L2","CX3CR1","FOLH1","HLA-DRA","TLR7","TMEM125","TNC","TREM2","TYROBP","COL18A1","GABRA1","GAD2","GLRA3","HTR2A","OXR1","SERPINI1","SLC6A13","SLC17A6","TCIRG1","UBQLN2","UCP2","AGPAT4-IT1","CHKB-CPT1B","COL3A1","ENSG00000205041","ENSG00000258674","ENSG00000273151","GATA2-AS1","HSP90AB4P","LINC01347","MIR24-2","MIRLET7BHG","NANOGP4")


#Clean Up Design 1
filt.glia.d1 = glia.res[rownames(glia.res) %in% my36,]
filt.ox.d1 = ox.res[rownames(ox.res) %in% my36,]
filt.TE.d1 = TE.res[rownames(TE.res) %in% my36,]
filt.GT.d1 = GT.res[rownames(GT.res) %in% my36,]
filt.GO.d1 = GO.res[rownames(GO.res) %in% my36,]
filt.TO.d1 = TO.res[rownames(TO.res) %in% my36,]
filt.glia.ond.d1 = glia.res.ond[rownames(glia.res.ond) %in% my36,]
filt.ox.ond.d1 = ox.res.ond[rownames(ox.res.ond) %in% my36,]
filt.TE.ond.d1 = TE.res.ond[rownames(TE.res.ond) %in% my36,]
filt.COND.d1 = COND.res[rownames(COND.res) %in% my36,]


#Add in RIN
FullPheno$RIN = NA
Clinical2 = Clinical
convertnames = gsub("-","\\.",Clinical2$ExternalSampleId)
Clinical2$ExternalSampleId = convertnames

for(i in 1:nrow(FullPheno)){
  
  for(j in 1:nrow(Clinical2)){
    
    if(FullPheno$Subject[i] == Clinical2$ExternalSampleId[j]){
      
      FullPheno$RIN[i] = Clinical2$RIN[j]
      
    }
    
  }
  
}

table(FullPheno$Subject == colnames(FullCount))

#Remove single sample with missing RIN
FullPheno_rin = FullPheno[-which(is.na(FullPheno$RIN)),]
FullCount_rin = FullCount[,-which(is.na(FullPheno$RIN))]
table(FullPheno_rin$Subject == colnames(FullCount_rin))

rCountData_rin = round(FullCount_rin,0)

dds_rin = DESeqDataSetFromMatrix(countData = rCountData_rin, colData = FullPheno_rin, design= ~ platform + RIN + Subtype, tidy=F) #Subtype must be last (DESeq2 vignette)
dds_rin$Subtype = relevel(dds_rin$Subtype,ref = "Control")
dseq_rin = DESeq(dds_rin,betaPrior=T)


#Pairwise "contrast()" for Design Equation 2
glia.res = results(dseq_rin,contrast = c("Subtype","GLIA","Control"))
filt.glia.d2 = glia.res[rownames(glia.res) %in% my36,]

ox.res = results(dseq_rin,contrast = c("Subtype","OX","Control"))
filt.ox.d2 = ox.res[rownames(ox.res) %in% my36,]

TE.res = results(dseq_rin,contrast = c("Subtype","TE","Control"))
filt.TE.d2 = TE.res[rownames(TE.res) %in% my36,]

GT.res = results(dseq_rin,contrast = c("Subtype","GLIA","TE"))
filt.GT.d2 = GT.res[rownames(GT.res) %in% my36,]

GO.res = results(dseq_rin,contrast = c("Subtype","GLIA","OX"))
filt.GO.d2 = GO.res[rownames(GO.res) %in% my36,]

TO.res = results(dseq_rin,contrast = c("Subtype","TE","OX"))
filt.TO.d2 = TO.res[rownames(TO.res) %in% my36,]

glia.res.ond = results(dseq_rin,contrast = c("Subtype","GLIA","OND"))
filt.glia.ond.d2 = glia.res.ond[rownames(glia.res.ond) %in% my36,]

ox.res.ond = results(dseq_rin,contrast = c("Subtype","OX","OND"))
filt.ox.ond.d2 = ox.res.ond[rownames(ox.res.ond) %in% my36,]

TE.res.ond = results(dseq_rin,contrast = c("Subtype","TE","OND"))
filt.TE.ond.d2 = TE.res.ond[rownames(TE.res.ond) %in% my36,]

COND.res = results(dseq_rin,contrast = c("Subtype","Control","OND"))
filt.COND.d2 = COND.res[rownames(COND.res) %in% my36,]


#Build design equation p-value matrix
npairs = 10
ndesign = 2
DE_RIN = data.frame(matrix(NA,length(my36),npairs*ndesign))
colnames(DE_RIN) = c("Design1_ALS-Glia_v_controls","Design2_ALS-Glia_v_controls","Design1_ALS-Ox_v_controls","Design2_ALS-Ox_v_controls","Design1_ALS-TD_v_controls","Design2_ALS-TD_v_controls","Design1_ALS-Glia_v_ALS-Ox","Design2_ALS-Glia_v_ALS-Ox","Design1_ALS-Glia_v_ALS-TD","Design2_ALS-Glia_v_ALS-TD","Design1_ALS-Ox_v_ALS-TD","Design2_ALS-Ox_v_ALS-TD","Design1_ALS-Glia_v_FTLD","Design2_ALS-Glia_v_FTLD","Design1_ALS-Ox_v_FTLD","Design2_ALS-Ox_v_FTLD","Design1_ALS-TD_v_FTLD","Design2_ALS-TD_v_FTLD","Design1_Control_v_FTLD","Design2_Control_v_FTLD")
rownames(DE_RIN) = my36


#Reference List to help with condensing dataframes (same order as DE_RIN matrix)
RefList = c("filt.glia.d1","filt.glia.d2","filt.ox.d1","filt.ox.d2","filt.TE.d1","filt.TE.d2","filt.GO.d1","filt.GO.d2","filt.GT.d1","filt.GT.d2","filt.TO.d1","filt.TO.d2","filt.glia.ond.d1","filt.glia.ond.d2","filt.ox.ond.d1","filt.ox.ond.d2","filt.TE.ond.d1","filt.TE.ond.d2","filt.COND.d1","filt.COND.d2")

for(i in 1:length(RefList)){
  
  tmp = get(RefList[i])
  
  for(j in 1:nrow(DE_RIN)){
    
    for(k in 1:nrow(tmp)){
      
      if(rownames(DE_RIN)[j] == rownames(tmp)[k]){
        
        DE_RIN[j,i] = tmp$padj[k]
        
      }
      
    }
    
  }
  
  
}

#Univariate analysis with Site as covariate (design = ~platform + Site + Subtype )

#Design Equation 1: design = ~ platform + Subtype (from UnivariateAnalysis script)
#Design Equation 2: design = ~ platform + Site + Subtype

load("D:/Jarrett/Research/Fall 2021/ProgBM/SOD1/ALSPatientStratification_UnivariateDatasets_SOD1.RData")

#Filter for subtype-specific features in Fig. 6
my36 = c("AIF1","APOC2","CD44","CHI3L2","CX3CR1","FOLH1","HLA-DRA","TLR7","TMEM125","TNC","TREM2","TYROBP","COL18A1","GABRA1","GAD2","GLRA3","HTR2A","OXR1","SERPINI1","SLC6A13","SLC17A6","TCIRG1","UBQLN2","UCP2","AGPAT4-IT1","CHKB-CPT1B","COL3A1","ENSG00000205041","ENSG00000258674","ENSG00000273151","GATA2-AS1","HSP90AB4P","LINC01347","MIR24-2","MIRLET7BHG","NANOGP4")


#Clean Up Design 1
filt.glia.d1 = glia.res[rownames(glia.res) %in% my36,]
filt.ox.d1 = ox.res[rownames(ox.res) %in% my36,]
filt.TE.d1 = TE.res[rownames(TE.res) %in% my36,]
filt.GT.d1 = GT.res[rownames(GT.res) %in% my36,]
filt.GO.d1 = GO.res[rownames(GO.res) %in% my36,]
filt.TO.d1 = TO.res[rownames(TO.res) %in% my36,]
filt.glia.ond.d1 = glia.res.ond[rownames(glia.res.ond) %in% my36,]
filt.ox.ond.d1 = ox.res.ond[rownames(ox.res.ond) %in% my36,]
filt.TE.ond.d1 = TE.res.ond[rownames(TE.res.ond) %in% my36,]
filt.COND.d1 = COND.res[rownames(COND.res) %in% my36,]

#Add Site to Pheno
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping")
Meta = read.csv("GSE153960_MetaData.txt")
convertnames = gsub("-","\\.",Meta$sample_id_alt)
Meta$sample_id_alt = convertnames

FiltMeta = Meta[Meta$sample_id_alt %in% FullPheno$Subject,]


reftable = table(FiltMeta$sample_id_alt)

removeind = rep(NA,nrow(FiltMeta)-length(reftable))
count = 1
for(i in 1:length(reftable)) {
  
  if(reftable[[i]] > 1){
    
    tmp = names(reftable[i])
    inds = which(FiltMeta$sample_id_alt == tmp)
    
    if(length(inds) == 2){
      removeind[count] = inds[1]
      count = count+1
    }else if(length(inds) == 3){
      removeind[seq(count,count+1,1)] = inds[1:2]
      count = count+2
    }else if(length(inds) > 3){
      cat("Problem")
    }
    
    
  }
  
}

SiteMeta = FiltMeta[-removeind,]



FullPheno$Site = NA

for(i in 1:nrow(FullPheno)){
  
  for(j in 1:nrow(SiteMeta)){
    
    if(FullPheno$Subject[i] == SiteMeta$sample_id_alt[j]){
      
      FullPheno$Site[i] = SiteMeta$project[j]
      
    }
    
  }
  
}

#Clean up site names
for(i in 1:nrow(FullPheno)){
  
  if(FullPheno$Site[i] == "NYGC ALS Consortium"){
    FullPheno$Site[i] = "NYGC"
  }else{
    FullPheno$Site[i] = "TargetALS"
  }
  
}

table(FullPheno$Site)

table(FullPheno$Subject == colnames(FullCount))


rCountData_site = round(FullCount,0)

dds_site = DESeqDataSetFromMatrix(countData = rCountData_site, colData = FullPheno, design= ~ platform + Site + Subtype, tidy=F) #Subtype must be last (DESeq2 vignette)
dds_site$Subtype = relevel(dds_site$Subtype,ref = "Control")
dseq_site = DESeq(dds_site,betaPrior=T)


#Pairwise "contrast()" for Design Equation 2
glia.res = results(dseq_site,contrast = c("Subtype","GLIA","Control"))
filt.glia.d2 = glia.res[rownames(glia.res) %in% my36,]

ox.res = results(dseq_site,contrast = c("Subtype","OX","Control"))
filt.ox.d2 = ox.res[rownames(ox.res) %in% my36,]

TE.res = results(dseq_site,contrast = c("Subtype","TE","Control"))
filt.TE.d2 = TE.res[rownames(TE.res) %in% my36,]

GT.res = results(dseq_site,contrast = c("Subtype","GLIA","TE"))
filt.GT.d2 = GT.res[rownames(GT.res) %in% my36,]

GO.res = results(dseq_site,contrast = c("Subtype","GLIA","OX"))
filt.GO.d2 = GO.res[rownames(GO.res) %in% my36,]

TO.res = results(dseq_site,contrast = c("Subtype","TE","OX"))
filt.TO.d2 = TO.res[rownames(TO.res) %in% my36,]

glia.res.ond = results(dseq_site,contrast = c("Subtype","GLIA","OND"))
filt.glia.ond.d2 = glia.res.ond[rownames(glia.res.ond) %in% my36,]

ox.res.ond = results(dseq_site,contrast = c("Subtype","OX","OND"))
filt.ox.ond.d2 = ox.res.ond[rownames(ox.res.ond) %in% my36,]

TE.res.ond = results(dseq_site,contrast = c("Subtype","TE","OND"))
filt.TE.ond.d2 = TE.res.ond[rownames(TE.res.ond) %in% my36,]

COND.res = results(dseq_site,contrast = c("Subtype","Control","OND"))
filt.COND.d2 = COND.res[rownames(COND.res) %in% my36,]


#Build design equation p-value matrix
npairs = 10
ndesign = 2
DE_SITE = data.frame(matrix(NA,length(my36),npairs*ndesign))
colnames(DE_SITE) = c("Design1_ALS-Glia_v_controls","Design2_ALS-Glia_v_controls","Design1_ALS-Ox_v_controls","Design2_ALS-Ox_v_controls","Design1_ALS-TD_v_controls","Design2_ALS-TD_v_controls","Design1_ALS-Glia_v_ALS-Ox","Design2_ALS-Glia_v_ALS-Ox","Design1_ALS-Glia_v_ALS-TD","Design2_ALS-Glia_v_ALS-TD","Design1_ALS-Ox_v_ALS-TD","Design2_ALS-Ox_v_ALS-TD","Design1_ALS-Glia_v_FTLD","Design2_ALS-Glia_v_FTLD","Design1_ALS-Ox_v_FTLD","Design2_ALS-Ox_v_FTLD","Design1_ALS-TD_v_FTLD","Design2_ALS-TD_v_FTLD","Design1_Control_v_FTLD","Design2_Control_v_FTLD")
rownames(DE_SITE) = my36


#Reference List to help with condensing dataframes (same order as DE_SITE matrix)
RefList = c("filt.glia.d1","filt.glia.d2","filt.ox.d1","filt.ox.d2","filt.TE.d1","filt.TE.d2","filt.GO.d1","filt.GO.d2","filt.GT.d1","filt.GT.d2","filt.TO.d1","filt.TO.d2","filt.glia.ond.d1","filt.glia.ond.d2","filt.ox.ond.d1","filt.ox.ond.d2","filt.TE.ond.d1","filt.TE.ond.d2","filt.COND.d1","filt.COND.d2")

for(i in 1:length(RefList)){
  
  tmp = get(RefList[i])
  
  for(j in 1:nrow(DE_SITE)){
    
    for(k in 1:nrow(tmp)){
      
      if(rownames(DE_SITE)[j] == rownames(tmp)[k]){
        
        DE_SITE[j,i] = tmp$padj[k]
        
      }
      
    }
    
  }
  
  
}


#Write out results

#write.csv(DE_SITE,"DifferentialExpression_withSITEcovariate_PeerReviewRound1.csv")

#Adjusting DE for Site does not seem to have much of an impact on the significance (TE vs Ox changes)

###############################################################################
###############################################################################

########## Site & RIN as covariates ##########

#Design Equation 1: design = ~ platform + Subtype (from UnivariateAnalysis script)
#Design Equation 2: design = ~ platform + Site + RIN + Subtype

load("D:/Jarrett/Research/Fall 2021/ProgBM/SOD1/ALSPatientStratification_UnivariateDatasets_SOD1.RData")

#Filter for subtype-specific features in Fig. 6
my36 = c("AIF1","APOC2","CD44","CHI3L2","CX3CR1","FOLH1","HLA-DRA","TLR7","TMEM125","TNC","TREM2","TYROBP","COL18A1","GABRA1","GAD2","GLRA3","HTR2A","OXR1","SERPINI1","SLC6A13","SLC17A6","TCIRG1","UBQLN2","UCP2","AGPAT4-IT1","CHKB-CPT1B","COL3A1","ENSG00000205041","ENSG00000258674","ENSG00000273151","GATA2-AS1","HSP90AB4P","LINC01347","MIR24-2","MIRLET7BHG","NANOGP4")


#Clean Up Design 1
filt.glia.d1 = glia.res[rownames(glia.res) %in% my36,]
filt.ox.d1 = ox.res[rownames(ox.res) %in% my36,]
filt.TE.d1 = TE.res[rownames(TE.res) %in% my36,]
filt.GT.d1 = GT.res[rownames(GT.res) %in% my36,]
filt.GO.d1 = GO.res[rownames(GO.res) %in% my36,]
filt.TO.d1 = TO.res[rownames(TO.res) %in% my36,]
filt.glia.ond.d1 = glia.res.ond[rownames(glia.res.ond) %in% my36,]
filt.ox.ond.d1 = ox.res.ond[rownames(ox.res.ond) %in% my36,]
filt.TE.ond.d1 = TE.res.ond[rownames(TE.res.ond) %in% my36,]
filt.COND.d1 = COND.res[rownames(COND.res) %in% my36,]

#Add Site to Pheno
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping")
Meta = read.csv("GSE153960_MetaData.txt")
Clinical = read.csv("CLINICAL_DATA_PRUDENCIO.csv")
convertnames = gsub("-","\\.",Meta$sample_id_alt)
Meta$sample_id_alt = convertnames

#Add in RIN
FullPheno$RIN = NA
Clinical2 = Clinical
convertnames = gsub("-","\\.",Clinical2$ExternalSampleId)
Clinical2$ExternalSampleId = convertnames

for(i in 1:nrow(FullPheno)){
  
  for(j in 1:nrow(Clinical2)){
    
    if(FullPheno$Subject[i] == Clinical2$ExternalSampleId[j]){
      
      FullPheno$RIN[i] = Clinical2$RIN[j]
      
    }
    
  }
  
}

table(FullPheno$Subject == colnames(FullCount))



FiltMeta = Meta[Meta$sample_id_alt %in% FullPheno$Subject,]


reftable = table(FiltMeta$sample_id_alt)

removeind = rep(NA,nrow(FiltMeta)-length(reftable))
count = 1
for(i in 1:length(reftable)) {
  
  if(reftable[[i]] > 1){
    
    tmp = names(reftable[i])
    inds = which(FiltMeta$sample_id_alt == tmp)
    
    if(length(inds) == 2){
      removeind[count] = inds[1]
      count = count+1
    }else if(length(inds) == 3){
      removeind[seq(count,count+1,1)] = inds[1:2]
      count = count+2
    }else if(length(inds) > 3){
      cat("Problem")
    }
    
    
  }
  
}

SiteMeta = FiltMeta[-removeind,]



FullPheno$Site = NA

for(i in 1:nrow(FullPheno)){
  
  for(j in 1:nrow(SiteMeta)){
    
    if(FullPheno$Subject[i] == SiteMeta$sample_id_alt[j]){
      
      FullPheno$Site[i] = SiteMeta$project[j]
      
    }
    
  }
  
}

#Clean up site names
for(i in 1:nrow(FullPheno)){
  
  if(FullPheno$Site[i] == "NYGC ALS Consortium"){
    FullPheno$Site[i] = "NYGC"
  }else{
    FullPheno$Site[i] = "TargetALS"
  }
  
}

table(FullPheno$Site)

#Single sample with missing RIN must be removed (incomplete design equation)
table(FullPheno$Subject == colnames(FullCount))
FullPheno_sr = FullPheno[-which(is.na(FullPheno$RIN)),]
FullCount_sr = FullCount[,-which(is.na(FullPheno$RIN))]

FullPheno_sr$RIN = scale(FullPheno_sr$RIN,center = T)

rCountData_rinsite = round(FullCount_sr,0)

dds_rinsite = DESeqDataSetFromMatrix(countData = rCountData_rinsite, colData = FullPheno_sr, design= ~ platform + Site + RIN + Subtype, tidy=F) #Subtype must be last (DESeq2 vignette)
dds_rinsite$Subtype = relevel(dds_rinsite$Subtype,ref = "Control")
dseq_rinsite = DESeq(dds_rinsite,betaPrior=T)


#Pairwise "contrast()" for Design Equation 2
glia.res = results(dseq_rinsite,contrast = c("Subtype","GLIA","Control"))
filt.glia.d2 = glia.res[rownames(glia.res) %in% my36,]

ox.res = results(dseq_rinsite,contrast = c("Subtype","OX","Control"))
filt.ox.d2 = ox.res[rownames(ox.res) %in% my36,]

TE.res = results(dseq_rinsite,contrast = c("Subtype","TE","Control"))
filt.TE.d2 = TE.res[rownames(TE.res) %in% my36,]

GT.res = results(dseq_rinsite,contrast = c("Subtype","GLIA","TE"))
filt.GT.d2 = GT.res[rownames(GT.res) %in% my36,]

GO.res = results(dseq_rinsite,contrast = c("Subtype","GLIA","OX"))
filt.GO.d2 = GO.res[rownames(GO.res) %in% my36,]

TO.res = results(dseq_rinsite,contrast = c("Subtype","TE","OX"))
filt.TO.d2 = TO.res[rownames(TO.res) %in% my36,]

glia.res.ond = results(dseq_rinsite,contrast = c("Subtype","GLIA","OND"))
filt.glia.ond.d2 = glia.res.ond[rownames(glia.res.ond) %in% my36,]

ox.res.ond = results(dseq_rinsite,contrast = c("Subtype","OX","OND"))
filt.ox.ond.d2 = ox.res.ond[rownames(ox.res.ond) %in% my36,]

TE.res.ond = results(dseq_rinsite,contrast = c("Subtype","TE","OND"))
filt.TE.ond.d2 = TE.res.ond[rownames(TE.res.ond) %in% my36,]

COND.res = results(dseq_rinsite,contrast = c("Subtype","Control","OND"))
filt.COND.d2 = COND.res[rownames(COND.res) %in% my36,]


#Build design equation p-value matrix
npairs = 10
ndesign = 2
DE_RINSITE = data.frame(matrix(NA,length(my36),npairs*ndesign))
colnames(DE_RINSITE) = c("Design1_ALS-Glia_v_controls","Design2_ALS-Glia_v_controls","Design1_ALS-Ox_v_controls","Design2_ALS-Ox_v_controls","Design1_ALS-TD_v_controls","Design2_ALS-TD_v_controls","Design1_ALS-Glia_v_ALS-Ox","Design2_ALS-Glia_v_ALS-Ox","Design1_ALS-Glia_v_ALS-TD","Design2_ALS-Glia_v_ALS-TD","Design1_ALS-Ox_v_ALS-TD","Design2_ALS-Ox_v_ALS-TD","Design1_ALS-Glia_v_FTLD","Design2_ALS-Glia_v_FTLD","Design1_ALS-Ox_v_FTLD","Design2_ALS-Ox_v_FTLD","Design1_ALS-TD_v_FTLD","Design2_ALS-TD_v_FTLD","Design1_Control_v_FTLD","Design2_Control_v_FTLD")
rownames(DE_RINSITE) = my36


#Reference List to help with condensing dataframes (same order as DE_SITE matrix)
RefList = c("filt.glia.d1","filt.glia.d2","filt.ox.d1","filt.ox.d2","filt.TE.d1","filt.TE.d2","filt.GO.d1","filt.GO.d2","filt.GT.d1","filt.GT.d2","filt.TO.d1","filt.TO.d2","filt.glia.ond.d1","filt.glia.ond.d2","filt.ox.ond.d1","filt.ox.ond.d2","filt.TE.ond.d1","filt.TE.ond.d2","filt.COND.d1","filt.COND.d2")

for(i in 1:length(RefList)){
  
  tmp = get(RefList[i])
  
  for(j in 1:nrow(DE_RINSITE)){
    
    for(k in 1:nrow(tmp)){
      
      if(rownames(DE_RINSITE)[j] == rownames(tmp)[k]){
        
        DE_RINSITE[j,i] = tmp$padj[k]
        
      }
      
    }
    
  }
  
  
}

