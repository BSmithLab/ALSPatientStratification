####################################
### SUPERVISED CLASSIFICATION ML ###
### ***   POST-PROCESSING   ***  ###
### Dr. Barbara Smith Lab        ###
### Written by: Jarrett Eshima   ###
### Date: June, 2021             ###
####################################

#Note to users: To get this code to run properly, complete the following steps in order:
#1) New file directory (results)
#2) Create new subdirectories for each condition run (e.g. ALS, ALS+HC, ALS+OND) - also supports assessment of multiple feature sets
#3) Create four additional directories in each feature set subdirectory, corresponding to each supervised classification model (KNN, RF, SVM, MLP)
#4) Organize files accordingly, run this script

#Read in all results / report files

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Supervised Classification/451/results/")
d = "C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Supervised Classification/451/results/"

setwd("D:/Jarrett/Research/Summer 2021/FinalClassificationResults_12-13-21/ALS_CA1000/results/")
d = "D:/Jarrett/Research/Summer 2021/FinalClassificationResults_12-13-21/ALS_CA1000/results/"

myfolders = list.files() #This must be a clean directory, with only folders associated with each SET OF FULL MODELS ran (i.e. variable number of sig genes used)

models = c("SVM","RF","KNN","MLP") #List the subdirectory names, each Supervised Classification set should have all models generated
#models = "MLP"
#This loop is going to reach into each of the directory and subdirectory locations, scan for .csv files, scan csv files for "report" or "result" string to properly label object, and finally read in data and assign it to the R object
for(i in 1:length(myfolders)){
  tmp = paste(d,myfolders[i],sep = "")
  setwd(tmp)
  for(j in 1:length(models)){
    wd = paste(d,myfolders[i],"/",models[j],sep="")
    resfiles = list.files(path = wd,pattern = ".csv")
    setwd(wd)
    for(y in 1:length(resfiles)){
      for(z in 1:nchar(resfiles)[y]){
        test = substr(resfiles[y],z,z+5)

        if(test == "report"){
          tmpp = read.csv(resfiles[y]) #File generated in Python/R scripts
          nam = paste(myfolders[i],models[j],"report",sep="")
          assign(nam,tmpp)
        }else if(test == "result"){
          tmpp2 = read.csv(resfiles[y]) #File generated in Python/R scripts
          nam = paste(myfolders[i],models[j],"results",sep="")
          assign(nam,tmpp2)
        }
      }
    }
  }
}

#Report contains F1 score for each k-fold
#Results contain true and predicted labels for testing samples, for all k-folds


################################### F1 Score - Model Summary

#Generate container for average F1 scores across k-fold cross-validation
F1 = data.frame(matrix(NA,nrow=length(myfolders),ncol = length(models)))
colnames(F1) = models
rownames(F1) = myfolders
head(F1)


for(i in 1:length(myfolders)){
  for(j in 1:length(models)){
    nam = paste(myfolders[i],models[j],"report",sep="")
    temp = get(nam)

    F1[i,j] = mean(temp$f1.score,na.rm = T)
  }
}
F1


###############################################################################################################################################


############# PLOT ML METRICS FOR EACH MODEL (separately)
par(mfrow=c(1,3))
modelnames = "Combined 1000" #same order as myfolders

for(i in 1:length(myfolders)){
  for(j in 1:length(models)){
    nam = paste(myfolders[i],models[j],"report",sep="")
    tmp = get(nam)

    nTE = table(tmp$X)[[3]]
    nOX = table(tmp$X)[[2]]
    nGL = table(tmp$X)[[1]]

    TE = data.frame(matrix(NA,nTE,ncol(tmp)))
    OX = data.frame(matrix(NA,nOX,ncol(tmp)))
    Glia = data.frame(matrix(NA,nGL,ncol(tmp)))

    colnames(TE) = colnames(OX) = colnames(Glia) = colnames(tmp)
    count1 = count2 = count3 = 1
    for(k in 1:nrow(tmp)){
      if(tmp$X[k] == "TE"){
        TE[count1,] = tmp[k,]
        count1 = count1+1
      }else if(tmp$X[k] == "OX"){
        OX[count2,] = tmp[k,]
        count2 = count2+1
      }else{
        Glia[count3,] = tmp[k,]
        count3 = count3+1
      }
    }

    cleannam = paste(modelnames[i],models[j],sep=" ")
    boxplot(TE$precision,OX$precision,Glia$precision,xaxt="n",main=c(cleannam,"Precision"),cex.axis = 1.5,cex.main=1.5,col=c("firebrick","navy","goldenrod1"),xlab="Subtype",ylab="Classifier Precision",cex.lab=1.5,ylim=c(0,1))
    axis(at=1:3,side=1,labels=c("TE","OX","Glia"),cex.axis=1.25)


    boxplot(TE$recall,OX$recall,Glia$recall,xaxt="n",main=c(cleannam,"Recall"),cex.axis = 1.5,cex.main=1.5,col=c("firebrick","navy","goldenrod1"),xlab="Subtype",ylab="Classifier Recall",cex.lab=1.5,ylim=c(0,1))
    axis(at=1:3,side=1,labels=c("TE","OX","Glia"),cex.axis=1.25)


    boxplot(TE$f1.score,OX$f1.score,Glia$f1.score,xaxt="n",main=c(cleannam,"F1 Score"),cex.axis = 1.5,cex.main=1.5,col=c("firebrick","navy","goldenrod1"),xlab="Subtype",ylab="Classifier F1",cex.lab=1.5,ylim=c(0,1))
    axis(at=1:3,side=1,labels=c("TE","OX","Glia"),cex.axis=1.25)

  }

}

############# PLOT ML METRICS FOR ALL MODELS (by feature set)

par(mfrow=c(1,1))
for(i in 1:length(myfolders)){
  for(j in 1:length(models)){
    nam = paste(myfolders[i],models[j],"report",sep="")
    tmp = get(nam)

    mdl1 = paste(models[j],"F1",sep="")
    assign(mdl1,tmp$f1.score)

    mdl2 = paste(models[j],"Precision",sep = "")
    assign(mdl2,tmp$precision)

    mdl3 = paste(models[j], "Recall",sep="")
    assign(mdl3,tmp$recall)

  }

  precombdata = recombdata = f1combdata = data.frame(matrix(NA,nrow(tmp),length(models)+1))
  colnames(precombdata) = colnames(recombdata) = colnames(f1combdata) = c("Subtype",models)
  precombdata$Subtype = recombdata$Subtype = f1combdata$Subtype = tmp$X

  for(k in 1:length(models)){
    tmpf1 = get(paste(models[k],"F1",sep=""))
    tmppr = get(paste(models[k],"Precision",sep = ""))
    tmpre = get(paste(models[k],"Recall",sep = ""))

    f1combdata[,k+1] = tmpf1
    precombdata[,k+1] = tmppr
    recombdata[,k+1] = tmpre
  }

  a = f1combdata$SVM
  b = f1combdata$RF
  c = f1combdata$KNN
  d = f1combdata$MLP
  dat = c(a,b,c,d)
  vr = rep(c("SVM","RF","KNN","MLP"),each=300)
  Subtype = rep(c("Glia","OX","TE"),100)
  dat2 = data.frame(vr,Subtype,dat)
  p = ggplot(dat2,aes(x=vr,y=dat,fill=Subtype)) + geom_boxplot()
  ttl = paste(modelnames[i],"- F1 Scores - 100-fold CV",sep=" ")
  p = p +ggtitle(ttl) + xlab("Classification Model") + ylab("F1 Score")
  p = p+scale_fill_manual(values = c("goldenrod1","navy","firebrick"))
  p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20),plot.title = element_text(size=24),legend.text = element_text(size=14),legend.title = element_text(size=14))
  p = p+theme_bw()
  p = p+ylim(0,1)
  #p = p+theme(legend.position = "none")
  p = p+theme(axis.title.x=element_text(vjust=-2))
  p = p+theme(axis.title.y=element_text(angle=90, vjust=6))
  p = p+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
  p = p+theme(plot.title = element_text(hjust = 0.5))
  print(p)


  a = precombdata$SVM
  b = precombdata$RF
  c = precombdata$KNN
  d = precombdata$MLP
  dat = c(a,b,c,d)
  vr = rep(c("SVM","RF","KNN","MLP"),each=300)
  Subtype = rep(c("Glia","OX","TE"),100)
  dat2 = data.frame(vr,Subtype,dat)
  p = ggplot(dat2,aes(x=vr,y=dat,fill=Subtype)) + geom_boxplot()
  ttl = paste(modelnames[i],"- Classifier Precision - 100-fold CV",sep=" ")
  p = p +ggtitle(ttl) + xlab("Classification Model") + ylab("Precision")
  p = p+scale_fill_manual(values = c("goldenrod1","navy","firebrick"))
  p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20),plot.title = element_text(size=24),legend.text = element_text(size=14),legend.title = element_text(size=14))
  p = p+theme_bw()
  p = p+ylim(0,1)
  #p = p+theme(legend.position = "none")
  p = p+theme(axis.title.x=element_text(vjust=-2))
  p = p+theme(axis.title.y=element_text(angle=90, vjust=6))
  p = p+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
  p = p+theme(plot.title = element_text(hjust = 0.5))
  print(p)


  a = recombdata$SVM
  b = recombdata$RF
  c = recombdata$KNN
  d = recombdata$MLP
  dat = c(a,b,c,d)
  vr = rep(c("SVM","RF","KNN","MLP"),each=300)
  Subtype = rep(c("Glia","OX","TE"),100)
  dat2 = data.frame(vr,Subtype,dat)
  p = ggplot(dat2,aes(x=vr,y=dat,fill=Subtype)) + geom_boxplot()
  ttl = paste(modelnames[i],"- Classifier Recall - 100-fold CV",sep=" ")
  p = p +ggtitle(ttl) + xlab("Classification Model") + ylab("Recall")
  p = p+scale_fill_manual(values = c("goldenrod1","navy","firebrick"))
  p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20),plot.title = element_text(size=24),legend.text = element_text(size=16),legend.title = element_text(size=16))
  p = p+theme_bw()
  p = p+ylim(0,1)
  #p = p+theme(legend.position = "none")
  p = p+theme(axis.title.x=element_text(vjust=-2))
  p = p+theme(axis.title.y=element_text(angle=90, vjust=6))
  p = p+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
  p = p+theme(plot.title = element_text(hjust = 0.5))
  print(p)

}

############# PLOT ML METRICS INDIVIDUALLY (a bit messy)


#Necessary function
quickfun = function(x,colnm){
  return(x[[colnm]])
}




for(i in 1:length(models)){
  for(j in 1:length(myfolders)){

    nam = paste(myfolders[j],models[i],"report",sep="")
    tmp = get(nam)

    mdl1 = paste(myfolders[j],"F1",sep="")
    assign(mdl1,tmp$f1.score)

    mdl2 = paste(myfolders[j],"Precision",sep = "")
    assign(mdl2,tmp$precision)

    mdl3 = paste(myfolders[j], "Recall",sep="")
    assign(mdl3,tmp$recall)

  }

  modsf1 = modspre = modsrec = data.frame(matrix(NA,nrow(tmp),length(myfolders)+1))
  colnames(modsf1) = colnames(modspre) = colnames(modsrec) = c("Subtype",myfolders)
  modsf1$Subtype = modspre$Subtype = modsrec$Subtype = tmp$X

  for(k in 1:length(myfolders)){
    tmpf1 = get(paste(myfolders[k],"F1",sep=""))
    tmppr = get(paste(myfolders[k],"Precision",sep = ""))
    tmpre = get(paste(myfolders[k],"Recall",sep = ""))

    modsf1[,k+1] = tmpf1
    modspre[,k+1] = tmppr
    modsrec[,k+1] = tmpre
    
    tmpnam = paste("m",k,sep="")
    assign(tmpnam,quickfun(modsf1,myfolders[k]))
    
    tmpnam = paste("n",k,sep="")
    assign(tmpnam,quickfun(modspre,myfolders[k]))
    
    tmpnam = paste("o",k,sep="")
    assign(tmpnam,quickfun(modsrec,myfolders[k]))
  }

  dat = c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22) #This is the only line that is manual
  vr = rep(myfolders,each=300)
  Subtype = rep(c("TE","OX","Glia"),100)
  dat2 = data.frame(vr,Subtype,dat)
  p = ggplot(dat2,aes(x=vr,y=dat,fill=Subtype)) + geom_boxplot()
  ttl = paste(models[i],"- F1 Scores - 100-fold CV",sep=" ")
  p = p +ggtitle(ttl) + xlab("Feature Set") + ylab("F1 Score")
  p = p+scale_fill_manual(values = c("goldenrod1","navy","firebrick"))
  p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20),plot.title = element_text(size=24),legend.text = element_text(size=14),legend.title = element_text(size=14))
  p = p + scale_y_continuous(limits = c(0,1))
  p = p+theme_bw()
  print(p)

  
  dat = c(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16,n17,n18,n19,n20,n21,n22) #This is the only line that is manual
  vr = rep(myfolders,each=300)
  Subtype = rep(c("TE","OX","Glia"),100)
  dat2 = data.frame(vr,Subtype,dat)
  p = ggplot(dat2,aes(x=vr,y=dat,fill=Subtype)) + geom_boxplot()
  ttl = paste(models[i],"- Precision - 100-fold CV",sep=" ")
  p = p +ggtitle(ttl) + xlab("Feature Set") + ylab("Precision")
  p = p+scale_fill_manual(values = c("goldenrod1","navy","firebrick"))
  p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20),plot.title = element_text(size=24),legend.text = element_text(size=14),legend.title = element_text(size=14))
  p = p+scale_y_continuous(limits = c(0,1))
  p = p+theme_bw()
  print(p)


  dat = c(o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,o13,o14,o15,o16,o17,o18,o19,o20,o21,o22) #This is the only line that is manual
  vr = rep(myfolders,each=300)
  Subtype = rep(c("TE","OX","Glia"),100)
  dat2 = data.frame(vr,Subtype,dat)
  p = ggplot(dat2,aes(x=vr,y=dat,fill=Subtype)) + geom_boxplot()
  ttl = paste(models[i],"- Recall - 100-fold CV",sep=" ")
  p = p +ggtitle(ttl) + xlab("Feature Set") + ylab("Recall")
  p = p+scale_fill_manual(values = c("goldenrod1","navy","firebrick"))
  p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20),plot.title = element_text(size=24),legend.text = element_text(size=14),legend.title = element_text(size=14))
  p = p+scale_y_continuous(limits = c(0,1))
  p = p+theme_bw()
  print(p)


}
