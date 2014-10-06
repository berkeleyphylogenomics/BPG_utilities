#!/clusterfs/ohana/software/bin/python2.6 -Wignore::DeprecationWarning
'''
This script scans through all proteins in the SNP table, finds the best GHG and all PFAM families, and runs INTREPID on
each protein/family combination.

Inputs: "snp_rost_data_unique" table in pfacts003_test database
Outputs: "snp_rost_intrepid_header" and "snp_rost_intrepid_detail" tables in pfacts003_test database 

Author: Curt Hansen
Created: 20Apr2012
Modified: 
'''


import sys, os, time
import pf_connect as db
import intrepid_functions as i
from datetime import datetime


class NullDevice:  # Can use this to capture print messages.
    def write(self,s):
        pass


def main():
    if len(sys.argv) != 4:
        print "Usage: %s <PBS_JOBID> <NUMTHREADS> <NUMPROTSPROCESSED>" % sys.argv[0]
        sys.exit(1)

    jobID=int(sys.argv[1])
    numThreads=int(sys.argv[2])
    numProtsProcessed=int(sys.argv[3])
    jobidalpha="00"+str(jobID)
#    sys.stdout=NullDevice() #Turn off print messages.
    startTime=datetime.now()
    pgmName='compute_intrepid_snp_scores'
    fout=open(pgmName+".out."+startTime.strftime("%Y%m%d.%Hh%M")+"."+jobidalpha[-3:],"w")
    fout.write("\nStarting program...\n")
    fout.write("Turned off deprecation warnings.\n")
   
    dbConn=db.connect_to_server('pfacts003_test','webuser')
    dbCur=db.get_cursor(dbConn)
    dbConn.set_isolation_level(0)

    proteinsList=db.get_query_list(dbCur,"SELECT DISTINCT u.accession FROM uniprot as u INNER JOIN snp_rost_master as s \
    ON s.accession=u.accession WHERE s.status='M' ORDER BY u.accession ASC LIMIT %s;", [numProtsProcessed])

    totCounter = 0
    totCounterThread = 0
    totResCounter = 0

    l1=65
    l2=15
    for prot in proteinsList:
        if totCounter % numThreads != jobID: #Check if this thread should be computing this protein in list.
            totCounter += 1
            continue #Skip this protein.
        
        protAccID = prot[0]
        totCounterThread += 1

        famsList=db.get_query_list(dbCur,"SELECT DISTINCT f.id, f.family_type_id, f.score \
        FROM family_sequence_taxa_2012_02_17 as fs \
        LEFT JOIN family as f ON fs.family_id=f.id WHERE fs.uniprot_accession=%s \
        ORDER BY f.family_type_id DESC, f.score DESC;",[protAccID])
        if len(famsList) == 0:  # No families found for protein.
            famid = 0
            famtype = 'N'
            famscore = 0
            msaok = False
            treeok = False
            idmapok = False
            scoresok = False
            dbCur.execute("INSERT INTO intrepid_header \
              (accession,fam_id,fam_type,fam_score,msaok,treeok,idmapok,scoresok,qty_scored,start_pos,start_time) \
               VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);",\
                   [protAccID,famid,famtype,famscore,msaok,treeok,idmapok,scoresok,0,0,startTime])
        else:
            flgGHGFound = False # Initialize to false.
            for fam in famsList: # Loop over all families for protein.
                famid = fam[0]
                famtype = fam[1]
                famscore = fam[2]
                if famtype=='C' or (famtype=='G' and flgGHGFound==False):
                    if famtype=='G':
                        flgGHGFound = True
                    results = i.run_intrepid(protAccID,famid) # Run INTREPID and get results.
                    if results[0] == False: # INTREPID run was not OK, results[1] will indicate why.
                        dbCur.execute("INSERT INTO intrepid_header \
                          (accession,fam_id,fam_type,fam_score,msaok,treeok,idmapok,scoresok,qty_scored,\
                           start_pos,start_time) \
                           VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);",\
                                 [protAccID,famid,famtype,famscore,results[1][0],results[1][1],results[1][2],\
                                  results[1][3],None,None,startTime])
                    else: # INTREPID run was OK, results[1] will be list of scores and results[2] will be start pos in seq.
                        for line in results[1]:
                            dbCur.execute("INSERT INTO intrepid_detail \
                              (accession, fam_id, position, resseq, chain, residue, icode, scoreimpmap, \
                               scorejsmap, scoreremap, scoreglobalmap, scoreglobalksmap, scoreglobalremap) \
                               VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);", \
                                 [protAccID,famid,line[0],line[1],line[2],line[3],\
                                  line[4],line[5],line[6],line[7],line[8],line[9],line[10]])
                        startpos = i.get_start_pos(dbCur,protAccID,results[2])
                        dbCur.execute("INSERT INTO intrepid_header \
                          (accession,fam_id,fam_type,fam_score,msaok,treeok,idmapok,scoresok,qty_scored,\
                           start_pos,start_time) \
                           VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);", \
                             [protAccID,famid,famtype,famscore,True,True,True,True,len(results[1]),startpos,startTime])
                        totResCounter += len(results[1]) #Increment counter of total number of detail lines added.
        totCounter+=1

    fout.write("\nResults\n")
    fout.write("-"*(l1+l2)+"\n")
    fout.write("Number of proteins scored:".ljust(l1)+str(totCounterThread).rjust(l2)+"\n")
    fout.write("Number of protein positions scored:".ljust(l1)+str(totResCounter).rjust(l2)+"\n") 
    endTime=datetime.now()
    diff=endTime-startTime
    fout.write("Elapsed time (HH:MM:SS):".ljust(l1)+str(diff).rjust(l2)+"\n\n")
    fout.close()
    dbConn.commit()
    dbConn.close()


if __name__=='__main__':
  main()
