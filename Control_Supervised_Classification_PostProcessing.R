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

setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Supervised Classification/Control/results/")
d = "C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Supervised Classification/Control/results/"

setwd("D:/Jarrett/Research/Summer 2021/FinalClassificationResults_12-13-21/Control/results/")
d = "D:/Jarrett/Research/Summer 2021/FinalClassificationResults_12-13-21/Control/results/"

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


#This will give the overall average of f1.scores for all 3 subtypes
for(i in 1:length(myfolders)){
  for(j in 1:length(models)){
    nam = paste(myfolders[i],models[j],"report",sep="")
    temp = get(nam)

    F1[i,j] = mean(temp$f1.score,na.rm = T)
  }
}
F1

#Get the average F1 score across supervised classifiers
# AvgF1 = rep(NA,nrow(F1))
# F1.num = matrix(as.numeric(unlist(F1)),ncol=ncol(F1),nrow=nrow(F1))
# for(i in 1:nrow(F1.num)){
#   AvgF1[i] = mean(F1.num[i,])
# }
# AvgF1 = data.frame(AvgF1)
# rownames(AvgF1) = rownames(F1)
# colnames(AvgF1) = "Average F1 score"


###############################################################################################################################################
###################################  MANUAL  ##################################################################################################
###############################################################################################################################################
#Manually reorder F1 object so that feature sets are not plotted alphabetically

F2 = F1
F2[1,] = F1[2,]
F2[2,] = F1[3,]
F2[3,] = F1[1,]

F2names = c("ALS+HC","ALS+OND","AllSubjects")
rownames(F2) = F2names
F2

###############################################################################################################################################

blank = rep(NA,nrow(F2))
F2 = cbind(blank,F2)
F2[,1] = rownames(F2)

F2$blank = factor(F2$blank,levels = F2$blank)
F2

#Bar Graphs 

library(ggplot2)
p = ggplot(F2, aes(x=blank,y=SVM)) + geom_bar(stat="identity",color = "red",fill="tomato")
p = p +ggtitle("Average F1 Scores 100-fold CV - LinearSVC") + xlab("Feature Set") + ylab("F1 Score")
p = p +scale_y_continuous(limits = c(0,1)) #+ scale_x_discrete((limits = c("Top 250 Group","Top 500 Group","Top 1000 Group","Top 500 All","Top 1000 All","Top 2000 All","Top 5000 All","MAD 10k")))
p = p+theme(axis.text = element_text(size=16), axis.title = element_text(size=18),plot.title = element_text(size=18))
p = p+theme_bw()
p


p = ggplot(F2, aes(x=blank,y=RF)) + geom_bar(stat="identity",color = "gold3",fill="gold2")
p = p +ggtitle("Average F1 Scores 100-fold CV - Random Forest, 1000 trees") + xlab("Feature Set") + ylab("F1 Score")
p = p +scale_y_continuous(limits = c(0,1)) #+ scale_x_discrete((limits = c("Top 250 Group","Top 500 Group","Top 1000 Group","Top 500 All","Top 1000 All","Top 2000 All","Top 5000 All","MAD 10k")))
p = p +theme(axis.text = element_text(size=12), axis.title = element_text(size=18),plot.title = element_text(size=18))
p = p+theme_bw()
p



p = ggplot(F2, aes(x=blank,y=KNN)) + geom_bar(stat="identity",color = "darkgreen",fill="green3")
p = p +ggtitle("Average F1 Scores 100-fold CV - K-Nearest Neighbor") + xlab("Feature Set") + ylab("F1 Score")
p = p +scale_y_continuous(limits = c(0,1)) #+ scale_x_discrete((limits = c("Top 250 Group","Top 500 Group","Top 1000 Group","Top 500 All","Top 1000 All","Top 2000 All","Top 5000 All","MAD 10k")))
p = p +theme(axis.text = element_text(size=12), axis.title = element_text(size=18),plot.title = element_text(size=18))
p = p+theme_bw()
p


p = ggplot(F2, aes(x=blank,y=MLP)) + geom_bar(stat="identity",color = "blue",fill="steelblue3")
p = p +ggtitle("Average F1 Scores 100-fold CV - Neural Net") + xlab("Feature Set") + ylab("F1 Score")
p = p +scale_y_continuous(limits = c(0,1)) #+ scale_x_discrete((limits = c("Top 250 Group","Top 500 Group","Top 1000 Group","Top 500 All","Top 1000 All","Top 2000 All","Top 5000 All","MAD 10k")))
p = p +theme(axis.text = element_text(size=12), axis.title = element_text(size=18),plot.title = element_text(size=18))
p = p+theme_bw()
p


#Point plots for F1 scores (more asthetic than descriptive)
rbPal <- colorRampPalette(c('red','blue'))
cols <- rbPal(8)[as.numeric(cut(F2$SVM,breaks = 8))]
plot.default(F2$blank,F2$SVM,pch=20,type = "p",cex = 2.5,cex.axis = 1.35,cex.lab = 1.5,cex.main = 1.75,xlab = "Features",ylab = "Avg F1 score",main="Average F1 Scores 100-fold CV - LinearSVC",xaxt="n",col=cols)
axis(at=seq(1,nrow(F2),1),labels = rownames(F2),side = 1,cex.axis = 1.15)

rbPal <- colorRampPalette(c('red','blue'))
cols <- rbPal(8)[as.numeric(cut(F2$RF,breaks = 8))]
plot.default(F2$blank,F2$RF,pch=20,type = "p",cex = 2.5,cex.axis = 1.25,cex.lab = 1.5,cex.main = 1.75,xlab = "Features",ylab = "Avg F1 score",main="Average F1 Scores 100-fold CV - Random Forest, 1000 trees",xaxt="n",col=cols)
axis(at=seq(1,nrow(F2),1),labels = rownames(F2),side = 1,cex.axis = 1.15)

rbPal <- colorRampPalette(c('red','blue'))
cols <- rbPal(8)[as.numeric(cut(F2$KNN,breaks = 8))]
plot.default(F2$blank,F2$KNN,pch=20,type = "p",cex = 2.5,cex.axis = 1.25,cex.lab = 1.5,cex.main = 1.75,xlab = "Features",ylab = "Avg F1 score",main="Average F1 Scores 100-fold CV - K-Nearest Neighbor",xaxt="n",col=cols)
axis(at=seq(1,nrow(F2),1),labels = rownames(F2),side = 1,cex.axis = 1.15)

rbPal <- colorRampPalette(c('red','blue'))
cols <- rbPal(8)[as.numeric(cut(F2$MLP,breaks = 8))]
plot.default(F2$blank,F2$MLP,pch=20,type = "p",cex = 2.5,cex.axis = 1.25,cex.lab = 1.5,cex.main = 1.75,xlab = "Features",ylab = "Avg F1 score",main="Average F1 Scores 100-fold CV - Neural Net",xaxt="n",col=cols)
axis(at=seq(1,nrow(F2),1),labels = rownames(F2),side = 1,cex.axis = 1.15)


############# PLOT WITHOUT SUBCATEGORIES
par(mfrow=c(1,3))
modelnames = c("All Subjects","ALS+HC","ALS+OND") #same order as myfolders
#1) Breakdown each model (8 feature sets X 4 models = 32) by subtype F1 score - prob a boxplot
for(i in 1:length(myfolders)){
  for(j in 1:length(models)){
    nam = paste(myfolders[i],models[j],"report",sep="")
    tmp = get(nam)

    nTE = table(tmp$X)[[4]]
    nOX = table(tmp$X)[[3]]
    nGL = table(tmp$X)[[2]]
    nCn = table(tmp$X)[[1]]

    TE = data.frame(matrix(NA,nTE,ncol(tmp)))
    OX = data.frame(matrix(NA,nOX,ncol(tmp)))
    Glia = data.frame(matrix(NA,nGL,ncol(tmp)))
    Control = data.frame(matrix(NA,nCn,ncol(tmp)))

    colnames(TE) = colnames(OX) = colnames(Glia) = colnames(Control) = colnames(tmp)
    count1 = count2 = count3 = count4 = 1
    for(k in 1:nrow(tmp)){
      if(tmp$X[k] == "TE"){
        TE[count1,] = tmp[k,]
        count1 = count1+1
      }else if(tmp$X[k] == "OX"){
        OX[count2,] = tmp[k,]
        count2 = count2+1
      }else if(tmp$X[k] == "GLIA"){
        Glia[count3,] = tmp[k,]
        count3 = count3+1
      }else{
        Control[count4,] = tmp[k,]
        count4 = count4+1
      }
    }

    cleannam = paste(modelnames[i],models[j],sep=" ")
    boxplot(TE$precision,OX$precision,Glia$precision,Control$precision,xaxt="n",main=c(cleannam,"Precision"),cex.axis = 1.5,cex.main=1.5,col=c("firebrick","navy","goldenrod1","gray50"),xlab="Group",ylab="Classifier Precision",cex.lab=1.5,ylim=c(0,1))
    axis(at=1:4,side=1,labels=c("ALS-TD","ALS-Ox","ALS-Glia","Control"),cex.axis=1.25)


    boxplot(TE$recall,OX$recall,Glia$recall,Control$recall,xaxt="n",main=c(cleannam,"Recall"),cex.axis = 1.5,cex.main=1.5,col=c("firebrick","navy","goldenrod1","gray50"),xlab="Group",ylab="Classifier Recall",cex.lab=1.5,ylim=c(0,1))
    axis(at=1:4,side=1,labels=c("ALS-TD","ALS-Ox","ALS-Glia","Control"),cex.axis=1.25)


    boxplot(TE$f1.score,OX$f1.score,Glia$f1.score,Control$f1.score,xaxt="n",main=c(cleannam,"F1 Score"),cex.axis = 1.5,cex.main=1.5,col=c("firebrick","navy","goldenrod1","gray50"),xlab="Group",ylab="Classifier F1",cex.lab=1.5,ylim=c(0,1))
    axis(at=1:4,side=1,labels=c("ALS-TD","ALS-Ox","ALS-Glia","Control"),cex.axis=1.25)

  }

}


############# PLOT WITH SUBCATEGORIES

#Maybe combine all this into a single plot? (i.e. F1 by subtype for all 4 models from a feature set)
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
  vr = rep(c("SVM","RF","KNN","MLP"),each=400)
  Subtype = rep(c("Control","Glia","Ox","TD"),100)
  dat2 = data.frame(vr,Subtype,dat)
  p = ggplot(dat2,aes(x=vr,y=dat,fill=Subtype)) + geom_boxplot()
  ttl = paste(modelnames[i],"- F1 Scores - 100-fold CV",sep=" ")
  p = p +ggtitle(ttl) + xlab("Classification Model") + ylab("F1 Score")
  p = p+scale_fill_manual(values = c("gray50","goldenrod1","navy","firebrick"))
  p = p+theme_bw()
  p = p+ylim(0,1)
  #p = p+theme(legend.position = "none")
  p = p+theme(axis.title.x=element_text(vjust=-2))
  p = p+theme(axis.title.y=element_text(angle=90, vjust=6))
  p = p+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
  p = p+theme(plot.title = element_text(hjust = 0.5))
  p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20),plot.title = element_text(size=20),legend.text = element_text(size=14),legend.title = element_text(size=14))
  print(p)


  a = precombdata$SVM
  b = precombdata$RF
  c = precombdata$KNN
  d = precombdata$MLP
  dat = c(a,b,c,d)
  vr = rep(c("SVM","RF","KNN","MLP"),each=400)
  Subtype = rep(c("Control","Glia","Ox","TD"),100)
  dat2 = data.frame(vr,Subtype,dat)
  p = ggplot(dat2,aes(x=vr,y=dat,fill=Subtype)) + geom_boxplot()
  ttl = paste(modelnames[i],"- Classifier Precision - 100-fold CV",sep=" ")
  p = p +ggtitle(ttl) + xlab("Classification Model") + ylab("Precision")
  p = p+scale_fill_manual(values = c("gray50","goldenrod1","navy","firebrick"))
  p = p+theme_bw()
  p = p+ylim(0,1)
  #p = p+theme(legend.position = "none")
  p = p+theme(axis.title.x=element_text(vjust=-2))
  p = p+theme(axis.title.y=element_text(angle=90, vjust=6))
  p = p+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
  p = p+theme(plot.title = element_text(hjust = 0.5))
  p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20),plot.title = element_text(size=20),legend.text = element_text(size=14),legend.title = element_text(size=14))
  print(p)


  a = recombdata$SVM
  b = recombdata$RF
  c = recombdata$KNN
  d = recombdata$MLP
  dat = c(a,b,c,d)
  vr = rep(c("SVM","RF","KNN","MLP"),each=400)
  Subtype = rep(c("Control","Glia","Ox","TD"),100)
  dat2 = data.frame(vr,Subtype,dat)
  p = ggplot(dat2,aes(x=vr,y=dat,fill=Subtype)) + geom_boxplot()
  ttl = paste(modelnames[i],"- Classifier Recall - 100-fold CV",sep=" ")
  p = p +ggtitle(ttl) + xlab("Classification Model") + ylab("Recall")
  p = p+scale_fill_manual(values = c("gray50","goldenrod1","navy","firebrick"))
  p = p+theme_bw()
  p = p+ylim(0,1)
  #p = p+theme(legend.position = "none")
  p = p+theme(axis.title.x=element_text(vjust=-2))
  p = p+theme(axis.title.y=element_text(angle=90, vjust=6))
  p = p+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
  p = p+theme(plot.title = element_text(hjust = 0.5))
  p = p+theme(axis.text = element_text(size=20), axis.title = element_text(size=20),plot.title = element_text(size=20),legend.text = element_text(size=14),legend.title = element_text(size=14))
  print(p)

}
