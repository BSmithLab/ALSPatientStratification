# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 16:20:41 2021

@author: jeshima
"""
import classifiersV3 as cl


def runSVMrej(params):
    print('SVMrej round '+str(params[0])+'...')
    svmRej = cl.Classifier_SVMrej(params[1], params[2])
    testPredLbls = svmRej.predict_labels(params[3])
    return [params[0], testPredLbls]