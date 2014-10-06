#!/clusterfs/ohana/software/bin/python2.7
'''
This script deletes current versions of two pfacts003 tables for SIFT predictions.

Inputs: None.
Outputs: "snp_SIFT_pred_header" and "snp_SIFT_pred_detail" tables in pfacts003 database.

Author: Curt Hansen
Created: 2012.09.05
Modified: 
'''

import pf_connect as db

def askForConfirm():
    while True:
        ans = raw_input('Warning!! This will delete all tables/views related to SIFT results. Do you want to continue? (y/n): ')
        if ans == 'y':
            break
        elif ans == 'n':
            sys.exit('Program aborted at user request!\n')
        else:
            print "You must enter 'y' or 'n'!"

def main():
    askForConfirm()
    
    dbConn=db.connect_to_server('pfacts003','webuser')
    dbCur=db.get_cursor(dbConn)
    dbConn.set_isolation_level(0)

    try:
        dbCur.execute("DROP TABLE snp_SIFT_pred_header;")
        print "\nPrevious version of 'snp_SIFT_pred_header' deleted."
    except:
        print "\nTable 'snp_SIFT_pred_header' does not exist and will be created."
    try:
        dbCur.execute("DROP TABLE snp_SIFT_pred_detail;")
        print "Previous version of 'snp_SIFT_pred_detail' deleted."
    except:
        print "Table 'snp_SIFT_pred_detail' does not exist and will be created."

    dbCur.execute("CREATE TABLE snp_SIFT_pred_header ( \
                      id serial PRIMARY KEY, \
                      time_stamp timestamp UNIQUE, \
                      dataset_used varchar \
                      );")
    print "Created new or replacement 'snp_SIFT_pred_header' table in database."

    dbCur.execute("CREATE TABLE snp_SIFT_pred_detail ( \
                      id serial PRIMARY KEY, \
                      time_stamp timestamp REFERENCES snp_SIFT_pred_header(time_stamp), \
                      accession varchar, \
                      aa_ref varchar, \
                      aa_pos integer, \
                      aa_mut varchar, \
                      prediction integer, \
                      score1 float, \
                      score2 float, \
                      score3 float, \
                      score4 float \
                      );")
    print "Created new or replacement 'snp_SIFT_pred_detail' table in database."

    dbConn.commit()
    dbConn.close()


if __name__=='__main__':
    main()
