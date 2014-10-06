#!/clusterfs/ohana/software/bin/python2.7
'''
This program produces result reports.
Author: Curt Hansen
Created: 2012.08.07
Modified:
'''

import numpy as n
import sys,os
from mlabwrap import mlab
import warnings
import pf_connect as db
warnings.simplefilter("ignore", n.ComplexWarning) #Need to include this to turn off message from inside mlabwrap.
        

def getChoice(message,validChoices,inclDesc):
    validNumbers = range(len(validChoices))
    if inclDesc:
        optionsList = ", ".join([str(each)+'-'+str(validChoices[each]) for each in validNumbers])
    else:
        optionsList = ", ".join([str(each) for each in validNumbers])
    while True:
        ans = raw_input(message+' [options: '+optionsList+']: ')
        try:
            ans = int(ans)
            if ans in validNumbers:
                return ans
            else:
                print "Error! Valid choices include: "+optionsList
        except:
            print "Error! Valid choices include: "+optionsList

def plot_pr_curve(data,startPoints,filename,legend_string):
    mlab.addpath('/home/cchansen/snp_analysis/common')
    xlabel = 'Recall'
    ylabel = 'Precision'
    xlim = n.array([0,1])
    ylim = n.array([0,1])
    mlab.create_pdf_plot(data,startPoints,'',filename,xlabel,14,ylabel,14,xlim,ylim,legend_string)

def main():
    #Check input.
    if len(sys.argv) != 2:
        print "\nUsage: %s <filename w/o pdf extension>\nProgram terminated!\n" % os.path.basename(sys.argv[0])
        mlab.close()
        sys.exit(1)
    else:
        filename = sys.argv[1]
        print "\nProgram will create a plot/figure file named "+filename+".pdf\n" 

    #Set up empty request list.
    plotsRequested = []

    #Prompt for type of report.
    options = ['None/Quit','PR Curve(s)']
    chosenReport = getChoice('Choose the number of the desired report',options,True)
    if chosenReport == 0: mlab.close(); sys.exit('Program terminated at user request!')

    dbConn = db.connect_to_server('pfacts003_test','webuser')
    dbCur = db.get_cursor(dbConn)
    dbConn.set_isolation_level(0)

    if chosenReport == 1:
        runsAvailable = db.get_query_list(dbCur,"SELECT time_stamp, dataset_used FROM snp_results_header ORDER BY time_stamp DESC;")
        if len(runsAvailable) == 0:
            print "No results available for chosen report type"
        else:
            idxFirst = 1; runsAvailable.insert(0,('None/Quit','NA'))
            print "Available runs include the following:"
            for time_stamp, dataset in runsAvailable[1:]:
                print "\t"+str(idxFirst)+"\t"+str(time_stamp)+"\t"+dataset
                idxFirst += 1
            chosenRun = getChoice("Choose the number of the desired run",range(len(runsAvailable)),False)
            if chosenRun == 0: mlab.close(); sys.exit('Program terminated at user request!')
            modelsAvailable = db.get_query_list(dbCur,"SELECT DISTINCT model_used FROM snp_results_detail WHERE time_stamp = %s ORDER BY model_used ASC;", \
                                                [runsAvailable[chosenRun][0]])
            modelsAvailable = [i[0] for i in modelsAvailable] # Flatten list of lists.
            modelsAvailable.insert(0,'None/Quit') # Add this option at start of list for user.
            print "For the run chosen, the following models have results:"
            idxSec = 1
            for each in modelsAvailable[1:]:
                print "\t"+str(idxSec)+"\t"+str(each)
                idxSec += 1
            while True:
                chosenModel = getChoice("Choose the number of a model to add/include or 'None' to continue",modelsAvailable,True)
                if chosenModel == 0:
                    break
                else:
                    chosenType = getChoice("Choose the number corresponding to level of detail or 'None' to continue",\
                                           ['None','Individual Folds','Average over Folds'],True)
                    if chosenType == 0:
                        continue
                    else:
                        plotsRequested.append((runsAvailable[chosenRun][0],modelsAvailable[chosenModel],chosenType))
                        
            # Collect data requested.
            dataToPlot = n.zeros((2,0)) # Create an empty 2-row, 0-col matrix that will used to horizontally append data.
            startPoints = n.array([]) # Initialize empty array that stores where each series starts.
            counter = 1 # Need to use base 1 because that's what Matlab uses.
            legend_string = ''
            for time_stamp,model,type in plotsRequested:
                if type == 1:
                    currFold = -1
                    data = db.get_query_list(dbCur,"SELECT fold_no, seq_no, val_recall, val_precision FROM snp_results_detail \
                           WHERE time_stamp = %s AND model_used = %s ORDER BY fold_no ASC, seq_no ASC;", [time_stamp,model])
                    for fold,seq,valrecall,valprec in data:
                        if fold != currFold:
                            startPoints = n.hstack([startPoints,n.array([counter])])
                            currFold = fold
                            legend_string = legend_string + 'Fold ' + str(fold) + ','
                        if valrecall != None and valprec != None:
                            dataToPlot = n.hstack([dataToPlot,n.array([[float(valrecall)],[float(valprec)]])])
                            counter += 1
                if type == 2:
                    legend_string = legend_string + model + ','
                    data = db.get_query_list(dbCur,"SELECT seq_no, recall, precision FROM snp_results_detail_avg \
                           WHERE time_stamp = %s AND model_used = %s ORDER BY seq_no ASC;", [time_stamp,model])
                    startPoints = n.hstack([startPoints,n.array([counter])])
                    for seq,valrecall,valprec in data:
                        if valrecall != None and valprec != None:
                            dataToPlot = n.hstack([dataToPlot,n.array([[float(valrecall)],[float(valprec)]])])
                            counter += 1
                        
    mlab.close()
    if startPoints.shape[0] <= 1: # Means that only one series is present.
        legend_string = 'NA'
    else:
        legend_string = legend_string[0:-1] # Strip off trailing comma.
    plot_pr_curve(dataToPlot,startPoints,filename,legend_string)

if __name__=='__main__':
  main()
