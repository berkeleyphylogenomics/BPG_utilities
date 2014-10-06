'''
This module performs logistic regression.

Inputs: database connection, training data, training labels, test data, test labels
Outputs: precision/recall curves

Author: Curt Hansen
Created: Aug 4, 2012
Modified:
'''

import sys, os, time
import pf_connect as db
import numpy as n
import math_functions as m
import logreg as l
import svmutil as svm


def simple_comp(X,coeffs,cutoff):
    # This function ignores the coeffs variable but includes it for consistency with other methos.
    y_pred = n.zeros((len(X),1))
    y_pred[(X>cutoff)] = 1
    return y_pred

    
def run_logreg_model(dbCur,trainY,trainX,testY,testX,time_stamp,round,num_cutoff_points):
    coeffs = l.compute_logreg_coeffs(trainX,trainY)
    performance_results,valid_ind = m.compute_pred_performance_curve(l.predict_classes,testY,testX,coeffs,num_cutoff_points)
    for idx in range(num_cutoff_points):
        if valid_ind[0,idx] == 1:
            recall = performance_results[0,idx]; precision = performance_results[1,idx]
        else:
            recall = None; precision = None
        dbCur.execute("INSERT INTO snp_results_detail (time_stamp,model_used,fold_no,seq_no,val_recall,val_precision) \
              VALUES(timestamp %s,'lr',%s,%s,%s,%s);",[time_stamp,round+1,idx,recall,precision])

    
def run_svm_model(dbCur,trainY,trainX,testY,testX,time_stamp,round,num_cutoff_points):
    trainY = (2*trainY-1).tolist(); testYSVM = (2*testY-1).tolist() # Convert labels to -1,+1 format and lists as required by libsvm.
    trainY = [i[0] for i in trainY]; testYSVM = [i[0] for i in testY] # Convert to list of floats from list of lists.
    trainX = trainX.tolist(); testX = testX.tolist() # Convert to list as required by libsvm.
    prob = svm.svm_problem(trainY,trainX)
    params = svm.svm_parameter('-b 1 -q')
    svmmodel = svm.svm_train(prob,params)
    p_label, p_acc, p_val = svm.svm_predict(testYSVM,testX,svmmodel,'-b 1')
    probs = n.array(p_val)[:,1]; probs.shape = (len(probs),1)
    performance_results,valid_ind = m.compute_pred_performance_curve(simple_comp,testY,probs,None,num_cutoff_points)
    for idx in range(num_cutoff_points):
        if valid_ind[0,idx] == 1:
            recall = performance_results[0,idx]; precision = performance_results[1,idx]
        else:
            recall = None; precision = None
        dbCur.execute("INSERT INTO snp_results_detail (time_stamp,model_used,fold_no,seq_no,val_recall,val_precision) \
              VALUES(timestamp %s,'svm',%s,%s,%s,%s);",[time_stamp,round+1,idx,recall,precision])
    
