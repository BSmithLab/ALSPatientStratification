# -*- coding: utf-8 -*-
"""
Created on Tue May 11 12:29:26 2021

@author: Samantha Oconnor, Christopher L. Plaisier, Jarrett Eshima
@institution: Arizona State University
"""


#Base
import os
from os.path import exists
import numpy as np
import pandas as pd
from copy import deepcopy
from multiprocessing import Pool
import pickle

import sys
sys.path.append("C:/users/jeshima/appdata/roaming/python/python38/Scripts")
sys.path.append("C:/users/jeshima/appdata/roaming/python/python38/site-packages")

import scanpy as sc

#Cross validation and metrics
from sklearn.utils.random import sample_without_replacement
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from clusim.clustering import Clustering, print_clustering
import clusim.sim as sim


#Custom classes for classification
import Eshima_classifiersV3 as cl

import SVM as svm

import RF as rf

import KNN as knn

# Load up loom file
print('\nLoading data...')
ALS_Subtype = sc.read_loom('pydat/Looms/ALS451_CA1000.loom')



ALS_Subtype.obs
ALS_Subtype.var_names

### Cross validation parsing
ntot = 451
numSamples = 136 #30% - THIS IS TESTING NUMBER (30.2%)
nfolds = 100
ncores = 8

# Initialize helper vars/indices for subsetting data (train/test)
nCells = ALS_Subtype.shape[0]
allInds = np.arange(0, nCells)

# Make folder to store downstream results
#for meth1 in ['SVMrej', 'RFpy', 'KNN', 'ACTINN']:
#    if not os.path.exists('pydat/results/'+meth1):
#        os.makedirs('pydat/results/'+meth1)
#        print ("Directory created")
#    else:
#        print("Directory already exists")


# Precompute testing and training data sets
trainInds = []
truelab = []
SID = []
testInds = []
mapMe_SVMrejRF = []
mapMe_KNN = []
mapMe_ACTINN = []
errorACTINN = []
for k in range(nfolds):
    samp1 = sample_without_replacement(nCells, numSamples, random_state = 1234 + k) #This generates the test indices
    testInds.append(samp1) #This assigns the specific samples to the test dataset
    truelab.extend(ALS_Subtype.obs['clust_ID'][samp1]) #This assigns the true subtype/class for test dataset samples
    SID.extend(ALS_Subtype.obs['SubID'][samp1])
    trainInds.append(np.setdiff1d(allInds, samp1)) #Remaining samples get assigned to training dataset
    mapMe_SVMrejRF.append([k, ALS_Subtype[trainInds[k],:], ALS_Subtype.obs['clust_ID'][trainInds[k]], ALS_Subtype[testInds[k],:]]) #When run in the loop, this generates 100 different training subsets
    mapMe_KNN.append([k, ALS_Subtype[trainInds[k],:], ALS_Subtype.obs['clust_ID'][trainInds[k]], ALS_Subtype[testInds[k],:]])


#Build Classifiers

#################
### SVMrej CV ###
#################

if not exists('pydat/results/SVMrej/ALS451_Subtypingreport_SVM_CV_100fold_70-30_CA1000.csv'):

    # Cross validation within dataset
    print('\nSVMrej cross-validation (k-fold = '+str(nfolds)+')...')
    with Pool(ncores) as p:
        res1 = p.map(svm.runSVMrej, mapMe_SVMrejRF)

    # Reassemble predictions
    res1 = dict(res1)
    pred = []
    for k in range(nfolds):
        pred.extend(res1[k])

    # Dataframe of true labels, predictions, probabilities for all iterations
    DF = pd.DataFrame({'Subject ID':SID,'True Labels':truelab, 'Predictions':pred})
    DF.to_csv('pydat/results/SVMrej/ALS451_Subtypingresults_SVM_CV_100fold_70-30_CA1000.csv')

    # Get classification report for each iteration
    performanceResults = []
    for k in range(nfolds):
        performanceResults.append(classification_report(truelab[slice((numSamples)*k, (numSamples)*(k+1), 1)], pred[slice((numSamples)*k, (numSamples)*(k+1), 1)], output_dict=True, zero_division=0))

    # Convert into a dataframe
    performDF = pd.concat([pd.DataFrame(i) for i in performanceResults], axis=1).T
    states1 = ['TE','OX','GLIA'] #Use this line when performing classification considering only the ALS patients
    #states1 = ['TE','OX','GLIA','Control'] #Use this line when performing classification considering ALS patients and Controls
    performDF = performDF.loc[[True if i in states1 else False for i in list(performDF.index)]]
    performDF['Classifier'] = 'SVMrej'
    performDF.to_csv('pydat/results/SVMrej/ALS451_Subtypingreport_SVM_CV_100fold_70-30_CA1000.csv')

#############
### RF CV ###
#############

if not exists('pydat/results/RFpy/ALS451_RF_CV_classification_report_CA1000_ntree1000_'+str(ALS_Subtype._n_vars)+'.csv'):

    # Cross validation within dataset
    print('\nRF cross-validation (k-fold = '+str(nfolds)+')...')
    with Pool(ncores) as p:
        res1 = p.map(rf.runRF, mapMe_SVMrejRF)

    # Reassemble predictions
    res1 = dict(res1)
    pred = []
    for k in range(nfolds):
        pred.extend(res1[k])
    
    
    # Dataframe of true labels, predictions, probabilities for all iterations
    DF = pd.DataFrame({'Subject ID':SID,'True Labels':truelab, 'Predictions':pred})
    DF.to_csv('pydat/results/RFpy/ALS451_Subtype_RF_CV_results_CA1000_ntree1000_'+str(ALS_Subtype._n_vars)+'.csv')

    # Get classification report for each iteration
    performanceResults = []
    for k in range(nfolds):
        performanceResults.append(classification_report(truelab[slice((numSamples)*k, (numSamples)*(k+1), 1)], pred[slice((numSamples)*k, (numSamples)*(k+1), 1)], output_dict=True, zero_division=0))

    # Convert into a dataframe
    performDF = pd.concat([pd.DataFrame(i) for i in performanceResults], axis=1).T
    states1 = ['TE','OX','GLIA']
    #states1 = ['TE','OX','GLIA','Control']
    performDF = performDF.loc[[True if i in states1 else False for i in list(performDF.index)]]
    performDF['Classifier'] = 'RFpy'
    performDF.to_csv('pydat/results/RFpy/ALS451_RF_CV_classification_report_CA1000_ntree1000_'+str(ALS_Subtype._n_vars)+'.csv')


##############
### KNN CV ###
##############

if not exists('pydat/results/KNN/ALS451_KNN_CV_classification_report_CA1000_k5p1_'+str(ALS_Subtype._n_vars)+'.csv'):
    # Cross validation within dataset
    print('\nKNN cross-validation (k-fold = '+str(nfolds)+')...')
    with Pool(ncores) as p:
        res1 = p.map(knn.runKNN, mapMe_KNN)
    
    
    # Reassemble predictions
    res1 = dict(res1)
    pred = []
    for k in range(nfolds):
        pred.extend(res1[k])
    

    # Dataframe of true labels, predictions, probabilities for all iterations
    DF = pd.DataFrame({'Subject ID':SID,'True Labels':truelab, 'Predictions':pred})
    DF.to_csv('pydat/results/KNN/ALS451_KNN_CV_results_CA1000_k5p1_'+str(ALS_Subtype._n_vars)+'.csv')

    # Get classification report for each iteration
    performanceResults = []
    for k in range(nfolds):
        performanceResults.append(classification_report(truelab[slice((numSamples)*k, (numSamples)*(k+1), 1)], pred[slice((numSamples)*k, (numSamples)*(k+1), 1)], output_dict=True, zero_division=0))

    # Convert into a dataframe
    performDF = pd.concat([pd.DataFrame(i) for i in performanceResults], axis=1).T
    states1 = ['TE','OX','GLIA']
    #states1 = ['TE','OX','GLIA','Control']
    performDF = performDF.loc[[True if i in states1 else False for i in list(performDF.index)]]
    performDF['Classifier'] = 'KNN'
    performDF.to_csv('pydat/results/KNN/ALS451_KNN_CV_classification_report_CA1000_k5p1_'+str(ALS_Subtype._n_vars)+'.csv')

