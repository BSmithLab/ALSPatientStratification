# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 09:17:51 2021

@author: jeshima
"""

##This is a script to perform sklearn K-nearest neighbors for supervised classification problems

#Reference: https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.KNeighborsClassifier.html

#Base
import os
from os.path import exists
import numpy as np
import pandas as pd
from copy import deepcopy
from multiprocessing import Pool
import pickle
from joblib import dump, load


#I have to append directories to my path 
import sys
sys.path.append("C:/users/jeshima/appdata/roaming/python/python38/Scripts")
sys.path.append("C:/users/jeshima/appdata/roaming/python/python38/site-packages")


#Cross validation and metrics
from sklearn.neighbors import KNeighborsClassifier
from sklearn.utils.random import sample_without_replacement
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from clusim.clustering import Clustering, print_clustering
import clusim.sim as sim

#ROC
import matplotlib.pyplot as plt
from itertools import cycle
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import roc_auc_score



###########################################################

#k-fold cross validation
Expression = pd.read_csv('pydat/PubDat/Holdout/Nova_Supervised_Classifier_OverlappingPlatform_All_1000.csv', header = 0, index_col=0, sep=',')
Expression2 = np.transpose(Expression)

Pheno = pd.read_csv('pydat/PubDat/Holdout/Nova_ALS_coldata_SUBTYPES.csv', header = 0, index_col=1, sep=',') #Use col 1 as the index so subject IDs are accessible
Pheno['Subtype']


ntot = 255
numSamples = 77 #30% - THIS IS TESTING NUMBER (30.2%)
nfolds = 100

# Initialize helper vars/indices for subsetting data (train/test)
nSub = Expression.shape[1]
allInds = np.arange(0, nSub)

#Create containers
trainInds = []
truelab = []
SID = []
testInds = []
trainLab = []

#Generate training and testing subsets according to numSamples (in testing)
for k in range(nfolds):
    samp1 = sample_without_replacement(nSub, numSamples, random_state = 1234+k) #This generates the test indices
    testInds.append([k,samp1]) #This assigns the specific samples to the test dataset
    truelab.extend([k,Pheno['Subtype'].iloc[samp1]]) #This assigns the true subtype/class for test dataset samples
    trainLab.extend([k,Pheno['Subtype'].iloc[np.setdiff1d(allInds, samp1)]])
    SID.extend([k,Pheno['Subject'].iloc[samp1]])
    trainInds.append([k,np.setdiff1d(allInds, samp1)]) #Remaining samples get assigned to training dataset

#Clean up the gene expression matrices
out = []
out2 = []
for k in range(nfolds):
    for i in trainInds[k]:
        out.extend([k,Expression2.iloc[i]])

for k in range(nfolds):
    for i in testInds[k]:
        out2.extend([k,Expression2.iloc[i]])

#Final clean up
trainExp = out
testExp = out2
testLab = truelab[1:200:2]
trainLab2 = trainLab[1:200:2]
SID2 = SID[1:200:2]

train = trainExp[3:400:4]
test = testExp[3:400:4]


#KNN Cross-Validation

DF = []
for k in range(nfolds):
    
    k2 = k+1
    print(str(nfolds)+ '-fold CV, ' +'Running KNN round '+str(k2)+'...')
    
    X = np.array(train[k])

    Y = np.array(trainLab2[k])

    clf = OneVsRestClassifier(KNeighborsClassifier(n_neighbors=5,p=1,algorithm='auto',weights="distance"))

    clf.fit(X, Y)

    #Get test set predictions
    X2 = np.array(test[k])

    #Predict new samples
    pred = clf.predict(X2)
    #Get classification probabilities
    probs = clf.predict_proba(X2)

    myseq = list(range(0,len(pred))) 

    #Obtain results from each CV iteration
    GLp = []
    TDp = []
    OXp = []
    for i in myseq:
        GLp.append(probs[i][0])
        OXp.append(probs[i][1])
        TDp.append(probs[i][2])
        
        
        
    res = []
    res.append(np.transpose(pred))
    res.append(np.transpose(testLab[k]))
    res.append(GLp)
    res.append(TDp)
    res.append(OXp)
    res.append(SID2[k])
    DF.append(pd.DataFrame(res))
    
    
    y = testLab[k]
    y2 = label_binarize(y,classes=['GLIA','OX','TE'])
    y_score = probs
    nclass = y2.shape[1]

    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(nclass):
        fpr[i], tpr[i], _ = roc_curve(y2[:, i], y_score[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])


    #Macro-averaging class-specific FPR and TPR values
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(nclass)]))
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(nclass):
        mean_tpr += np.interp(all_fpr,fpr[i],tpr[i])

    mean_tpr /= nclass
    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"],tpr["macro"])

    fig, ax = plt.subplots()
    lw =2 
    plt.plot(
        fpr["macro"],
        tpr["macro"],
        label="macro-average ROC curve (area = {0:0.2f})".format(roc_auc["macro"]),
        color="darkgreen",
        linestyle="dashed",
        linewidth=1.5,
        )

    #ALS-Glia ROC
    plt.plot(
            fpr[0],
            tpr[0],
            color="goldenrod",
            lw=lw,
            label="ALS-Glia ROC curve (area = {1:0.2f})".format(i, roc_auc[0]),
            )
    #ALS-Ox ROC
    plt.plot(
            fpr[1],
            tpr[1],
            color="navy",
            lw=lw,
            label="ALS-Ox ROC curve (area = {1:0.2f})".format(i, roc_auc[1]),
            )
    #ALS-TD ROC
    plt.plot(
            fpr[2],
            tpr[2],
            color="firebrick",
            lw=lw,
            label="ALS-TD ROC curve (area = {1:0.2f})".format(i,roc_auc[2]),
            )

    plt.plot([0, 1], [0, 1], "k--", lw=lw)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("Testing Cohort ROC")
    plt.legend(loc="lower right")
    plt.show()
    


    
DF2 = pd.concat(DF)
DF2.to_csv('pydat/results/KNN/Testing_ALS451_KNN_results_100CV_OA1000_p1_k3_dist.csv')

dump(clf,'KNN_OA1000_p1_k3_dist.joblib')

fig.savefig('pydat/ROC/NovaSeq_Test_ROC_KNN_OA1000_p1_k3_dist.svg')




#Use ALSPatientStratification_CVresults_to_F1report.R script to convert results file to classification report



