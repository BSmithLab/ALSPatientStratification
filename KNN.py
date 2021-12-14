# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 09:24:52 2021

@author: jeshima
"""

import classifiersV3 as cl

def runKNN(params):
        print('RF round '+str(params[0])+'...')
        KNN = cl.Classifier_KNN(params[1], params[2])
        testPredLbls = KNN.predict_labels(params[3])
        return [params[0], testPredLbls]