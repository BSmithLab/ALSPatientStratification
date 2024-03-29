
## Convert MLP Neural Net results to report file with performance metrics (Precision, Recall, F1 scores)

#Written By: Jarrett Eshima
#For: Dr. Barbara Smith Lab

#4 directories corresponding to the 4 classification models considered
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Supervised Classification/Holdout/results/KNN")
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Supervised Classification/Holdout/results/MLP")
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Supervised Classification/Holdout/results/RF")
setwd("C:/Users/jeshima/Documents/Smith Lab/Summer 2021/ALS Subtyping/Supervised Classification/Holdout/results/SVM")


#Necessary function
quickfun = function(x,colnm){
  return(x[[colnm]])
}

#Adjust depending on the classification model used
type = "Training"
type = "Holdout"
model = "MLP"

###########################################################################################
#Provide the list of feature sets used in MLP classifier development 

#This code supports multi-file processing, as long as the naming scheme is consistent
#list = c()

#Testing Results
list = "OA1000_adam_tanh"

for(i in 1:length(list)){
  tmp = read.csv(paste("Testing_ALS451_MLP_results_3layer100Units_100CV_",list[i],".csv",sep="")) #File generated in Python scripts
  assign(list[i],tmp)
  rm(tmp)
}

#Holdout Results
list = "OA1000_adam_tanh"

for(i in 1:length(list)){
  tmp = read.csv(paste("HiSeqHoldout_MLP_results_3layer100Units_",list[i],".csv",sep="")) #File generated in Python scripts
  assign(list[i],tmp)
  rm(tmp)
}


###########################################################################################

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


#Convert to report file (values are shifted by 1 here due to column headers in the summary report file)

#Format: precision | recall | f1-score | support| Classifier
#GLIA | OX | TD

ngroup=3
nfold=1
for(i in 1:length(list)){

  report = data.frame(matrix(NA,nrow=ngroup*nfold+1,ncol=6)) #no column names, header=F

  report[,6] = "MLP"
  report[,1] = c("",rep(c("GLIA","OX","TD"),nfold)) #may have to change the order manually depending on future reports
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
    gcount = ocount=tcount=TPg=TPo=TPt=0
    for(k in 1:nrow(tmp)){
      if(Pred[j,k] == "GLIA"){
        gcount = gcount+1 #this is TP+FP

        if(Pred[j,k] == TS[j,k]){
          TPg = TPg+1
        }


      }else if(Pred[j,k] == "OX"){
        ocount = ocount+1

        if(Pred[j,k] == TS[j,k]){
          TPo = TPo+1
        }

      }else if(Pred[j,k] == "TE"){
        tcount = tcount+1

        if(Pred[j,k] == TS[j,k]){
          TPt = TPt+1
        }
      }


    }
    Precisiong = TPg/gcount
    Precisiono = TPo/ocount
    Precisiont = TPt/tcount
    start= (j*ngroup)-1
    stop = (j*ngroup)+1
    report$X2[start:stop] = c(Precisiong,Precisiono,Precisiont)

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
    gcount = ocount=tcount=TPg=TPo=TPt=0
    for(k in 1:nrow(tmp)){
      if(TS[j,k] == "GLIA"){
        gcount = gcount+1 #this is TP+FN

        if(Pred[j,k] == TS[j,k]){
          TPg = TPg+1
        }


      }else if(TS[j,k] == "OX"){
        ocount = ocount+1

        if(Pred[j,k] == TS[j,k]){
          TPo = TPo+1
        }

      }else if(TS[j,k] == "TE"){
        tcount = tcount+1

        if(Pred[j,k] == TS[j,k]){
          TPt = TPt+1
        }
      }
    }
    Recallg = TPg/gcount
    Recallo = TPo/ocount
    Recallt = TPt/tcount
    start= (j*ngroup)-1
    stop = (j*ngroup)+1
    report$X3[start:stop] = c(Recallg,Recallo,Recallt)

  }

  #Calculate F1-score
  for(m in 1:(nrow(report)-1)){

    pcn = as.numeric(report$X2[m+1])
    rcl = as.numeric(report$X3[m+1])

    report$X4[m+1] = (2*pcn*rcl)/(pcn+rcl)
  }



  #Write csv for report

  file = paste(list[i],"_report_",model,"_",type,".csv",sep="")
  write.table(report,file,sep = ",",col.names = F,row.names = F)

}





