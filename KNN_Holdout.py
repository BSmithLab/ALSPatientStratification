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

#Read in Holdout dataset (HiSeq Cohort)
Holdout = pd.read_csv('pydat/PubDat/Holdout/Hi_Supervised_Classifier_OverlappingPlatform_All_1000.csv', header = 0, index_col=0, sep=',')
Holdout2 = np.transpose(Holdout)
Validation = Holdout2

#Read in Holdout labels
HoldoutPheno = pd.read_csv('pydat/PubDat/Holdout/Hi_ALS_coldata_SUBTYPES.csv', header = 0, index_col=1, sep=',') #Use col 1 as the index so subject IDs are accessible
ValidationLab = HoldoutPheno['Subtype']
SIDH = HoldoutPheno['Subject']

#Predict holdout labels using MLP classifier built on NovaSeq cohort
clf = load('KNN_OA1000_p1_k3_dist.joblib')
pred = clf.predict(Validation)
probs = clf.predict_proba(Validation)
GLp = []
TDp = []
OXp = []
myseq = list(range(0,len(pred)))
for i in myseq:
    GLp.append(probs[i][0])
    OXp.append(probs[i][1])
    TDp.append(probs[i][2])
    
res = []
res.append(np.transpose(pred))
res.append(np.transpose(ValidationLab))
res.append(GLp)
res.append(TDp)
res.append(OXp)
res.append(SIDH)
DF = []
DF.append(pd.DataFrame(res))
DF2 = pd.concat(DF)
DF2.to_csv('pydat/results/KNN/HiSeqHoldout_KNN_results_OA1000_p1_k3_dist.csv')


##Construct multiclass ROC plot
y = ValidationLab
y2 = label_binarize(y,classes=['GLIA','OX','TE'])
y_score = probs
nclass = y2.shape[1]

#Get individual class TPR and FPR values
fpr = dict()
tpr = dict()
roc_auc = dict()
roc_auc2 = dict()
for i in range(nclass):
    fpr[i], tpr[i], _ = roc_curve(y2[:, i], y_score[:, i],pos_label=1)
    roc_auc[i] = auc(fpr[i], tpr[i])

#Micro-averaging class-specific FPR and TPR values
fpr["micro"], tpr["micro"], _ = roc_curve(y2.ravel(), y_score.ravel())
roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

#Macro-averaging class-specific FPR and TPR values
all_fpr = np.unique(np.concatenate([fpr[i] for i in range(nclass)]))
mean_tpr = np.zeros_like(all_fpr)
for i in range(nclass):
    mean_tpr += np.interp(all_fpr,fpr[i],tpr[i])

mean_tpr /= nclass
fpr["macro"] = all_fpr
tpr["macro"] = mean_tpr
roc_auc["macro"] = auc(fpr["macro"],tpr["macro"])

#Plot Micro-average ROC and class-specifics (4 ROC curves total)
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


plt.plot(
        fpr[0],
        tpr[0],
        color="goldenrod",
        lw=lw,
        label="ALS-Glia ROC curve (area = {1:0.2f})".format(i, roc_auc[0]),
    )
plt.plot(
        fpr[1],
        tpr[1],
        color="navy",
        lw=lw,
        label="ALS-Ox ROC curve (area = {1:0.2f})".format(i, roc_auc[1]),
    )
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
plt.title("Holdout ROC")
plt.legend(loc="lower right")
plt.show()

fig.savefig('pydat/ROC/HiSeq_Holdout_ROC_KNN_OA1000_p1_k3_dist.svg')