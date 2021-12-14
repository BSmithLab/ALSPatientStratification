# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 09:24:52 2021

@author: jeshima
"""

import classifiersV3 as cl

def runRF(params):
        print('RF round '+str(params[0])+'...')
        RF = cl.Classifier_RF(params[1], params[2])
        testPredLbls = RF.predict_labels(params[3])
        return [params[0], testPredLbls]