#!/clusterfs/ohana/software/bin/python2.7
'''
This is the main SNP prediction model based on logistic regression.

Inputs: snp dataset and various feature scores (e.g., INTREPID scores) per SNP
Outputs: precision/recall curves

Author: Curt Hansen
Created: July 13, 2012
Modified:
'''

#Import modules.
import sys, os, time
import pf_connect as db
from datetime import datetime
import numpy as n
import math_functions as m
import feature_extraction as f
import stat_models as models


# Lists to record run options.
modelsUsed = []
featuresUsed = []


def askWhichModels(modelName):
    while True:
        ans = raw_input('Do you want to run the ' + modelName + ' model? (y/n): ')
        if ans == 'y':
            modelsUsed.append(modelName)
            break
        elif ans == 'n':
            break
        else:
            print "You must enter 'y' or 'n'!"


def main():
    if len(sys.argv) != 3:
        print "Usage: %s <DATASET> <NUMFOLDS>\nProgram terminated!" % sys.argv[0]
        sys.exit(1)
    dataset = sys.argv[1]
    if dataset != 'rost':
        print "Dataset '%s' not supported.\nProgram terminated!" % dataset
        sys.exit(1)
    numFolds = int(sys.argv[2])
    askWhichModels('log reg')
    askWhichModels('SVM')

    startTime = datetime.now()
    formattedStartTime = startTime.strftime("%Y%m%dd%Hh%M")
    formattedDBTimeStamp = startTime.strftime("%Y-%m-%d %H:%M:%S")
    pgmName = 'run_model'
    fout = open(pgmName+".out."+formattedStartTime,"w")
    fout.write("\nStarting program...\n")

    dbConn = db.connect_to_server('pfacts003_test','webuser')
    dbCur = db.get_cursor(dbConn)
    dbConn.set_isolation_level(0)

    #Get entire dataset with all required features.
    completeSNPList = db.get_query_list(dbCur,"SELECT id, label FROM snp_rost_detail_valid ORDER BY id ASC;")
    numSNPsPre = len(completeSNPList)
    idsWithFeatures = n.array(completeSNPList)[:,0]; y = n.array(completeSNPList)[:,1]
    idsWithFeatures.shape = (numSNPsPre,1); y.shape = (numSNPsPre,1) # Force as col vectors.
    
    X = n.ones((numSNPsPre,1)) #Prepend a col of ones for the intercept.
    featuresUsed.append('Intercept')

    validIndices,newFeatures = f.extract_intrepid_features(dbCur,idsWithFeatures)
    idsWithFeatures = idsWithFeatures[validIndices]
    y = y[validIndices]; X = n.hstack([X[validIndices,:],newFeatures])
    featuresUsed.append('INTREPID scores')
    numSNPsPost = len(y)

    validIndices,newFeatures = f.extract_blossum62_scores(dbCur,idsWithFeatures)
    idsWithFeatures = idsWithFeatures[validIndices]
    y = y[validIndices]; X = n.hstack([X[validIndices,:],newFeatures])
    featuresUsed.append('BLOSSUM62 scores')
    numSNPsPost = len(y)

    #Determine partition of dataset into folds.
    foldAssignments = m.get_shuffled_discrete_uniform_indices(0,numFolds-1,numSNPsPost)

    #Create result header record to document run.
    dbCur.execute("INSERT INTO snp_results_header (time_stamp, dataset_used) VALUES (timestamp %s, %s);",[formattedDBTimeStamp,dataset])

    #Peform k-fold cross-validation and produce results (per fold and overall).
    num_cutoff_points = 101
    results_all_folds = n.zeros((2*numFolds,num_cutoff_points))
    average_performance_results = n.zeros((2,num_cutoff_points)) #This will hold the overall results, first as a total and then as an average.
    average_performance_num_contributors = n.zeros((1,num_cutoff_points)) #This keeps track of how many rounds contributed to the total for a cutoff point.
    for round in range(0,numFolds):
        trainX = X[n.nonzero(foldAssignments != round)[0],:]
        trainY = y[n.nonzero(foldAssignments != round)[0],:]
        testX = X[n.nonzero(foldAssignments == round)[0],:]
        testY = y[n.nonzero(foldAssignments == round)[0],:]
        if 'log reg' in modelsUsed:
            models.run_logreg_model(dbCur,trainY,trainX,testY,testX,formattedDBTimeStamp,round,num_cutoff_points)
        if 'SVM' in modelsUsed:
            models.run_svm_model(dbCur,trainY,trainX,testY,testX,formattedDBTimeStamp,round,num_cutoff_points)

    l1=65
    l2=15
    fout.write("\nResults\n")
    fout.write("-"*(l1+l2)+"\n")
    fout.write("Dataset used:".ljust(l1)+str(dataset).rjust(l2)+"\n")
    fout.write("Number of folds in cross validation:".ljust(l1)+str(numFolds).rjust(l2)+"\n")
    fout.write("Number of SNPs included in dataset before feature extraction:".ljust(l1)+str(numSNPsPre).rjust(l2)+"\n")
    fout.write("Number of SNPs included in dataset post feature extraction:".ljust(l1)+str(numSNPsPost).rjust(l2)+"\n")
    fout.write("Models used:"+"\n")
    for modelDesc in modelsUsed:
        fout.write("\t"+modelDesc+"\n")
    fout.write("Features used:"+"\n")
    for featureDesc in featuresUsed:
        fout.write("\t"+featureDesc+"\n")
    fout.write("See database for storage of results.\n")
    endTime=datetime.now()
    diff=endTime-startTime
    fout.write("Elapsed time (HH:MM:SS):".ljust(l1)+str(diff).rjust(l2)+"\n\n")
    fout.close()
    dbConn.commit()
    dbConn.close()


if __name__=='__main__':
  main()
