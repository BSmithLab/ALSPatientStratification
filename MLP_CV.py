# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 09:17:51 2021

@author: jeshima
"""

##This is a script to perform sklearn multi-layer perceptron (neural network) for supervised classification problems

#Reference: https://scikit-learn.org/stable/modules/neural_networks_supervised.html

#Base
import os
from os.path import exists
import numpy as np
import pandas as pd
from copy import deepcopy
from multiprocessing import Pool
import pickle


#I have to append directories to my path 
import sys
sys.path.append("C:/users/jeshima/appdata/roaming/python/python38/Scripts")
sys.path.append("C:/users/jeshima/appdata/roaming/python/python38/site-packages")


#Cross validation and metrics
from sklearn.neural_network import MLPClassifier
from sklearn.utils.random import sample_without_replacement
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from clusim.clustering import Clustering, print_clustering
import clusim.sim as sim

#ML
from sklearn.neural_network import MLPClassifier



###########################################################

#k-fold cross validation

Expression = pd.read_csv('pydat/PubDat/Supervised_Classifier_CombinedPlatform_All_1000.csv', header = 0, index_col=0, sep=',')


Expression2 = np.transpose(Expression)

Pheno = pd.read_csv('pydat/PubDat/ALS451_coldata_SUBTYPES.csv', header = 0, index_col=1, sep=',') #Use col 1 as the index so subject IDs are accessible
Pheno['Subtype']


ntot = 451
numSamples = 136 #30% - THIS IS TESTING NUMBER (30.2%)
nfolds = 100
#ncores = 8 #Number of cores for multi-core processing. This code does not currently support multi-core processing

# Initialize helper vars/indices for subsetting data (train/test)
nSub = Expression.shape[1]
allInds = np.arange(0, nSub)

#Create containers
trainInds = []
truelab = []
SID = []
testInds = []
trainLab = []
mapMLP = []

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


#MLP Cross-Validation

DF = []
for k in range(nfolds):
    
    print(str(nfolds)+ '-fold CV, ' +'Running Neural Net round '+str(k)+'...')
    
    X = np.array(train[k])

    Y = np.array(trainLab2[k])

    clf = MLPClassifier(solver='adam', alpha=1e-4, hidden_layer_sizes = (100,100,100), random_state=1,learning_rate='constant',learning_rate_init=0.0001,max_iter=10000)


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
    TEp = []
    OXp = []
    for i in myseq:
        GLp.append(probs[i][0])
        TEp.append(probs[i][1])
        OXp.append(probs[i][2])
        
        
    res = []
    res.append(np.transpose(pred))
    res.append(np.transpose(testLab[k]))
    res.append(GLp)
    res.append(TEp)
    res.append(OXp)
    res.append(SID2[k])
    
        
    DF.append(pd.DataFrame(res))

    
    
DF2 = pd.concat(DF)
DF2.to_csv('pydat/results/MLP/ALS451_MLP_results_3layer100Units_100CV_GO5000.csv')

#Use ALSPatientStratification_MLP_CV_to_F1.R script to convert results file to classification report

#Hidden_layer_sizes Explanation
#There are two visible layers, the input and the output
#Input and output layers can have more than one "unit" (aka neuron)
#Units typically represent different "pieces" of data for example Gene expression after knockdown and unperturbed expression
#The hidden layers consist of int units (in our case, 100), where each input unit is mapped to each unit in the hidden layer using a weight (w_ij)
#The single unit combines the inputs using a linear combination (sum)
#The summed inputs (net input) is fed into the activation function
#The output of the activation function is fed into the next hidden layer with int units/neurons
#This layer also has transfer weights (learning) associated with each neuron in the previous hidden layer
#Process is repeated until all hidden layers are calculated
#Final hidden layer feeds into the output layer (in our case, single output is ALS subtype)

