#!/clusterfs/ohana/software/bin/python2.7
'''
This script deletes current versions of two pfacts003 tables for precision and recall results and recreates them before a new computation in a separate program.

Inputs: None.
Outputs: "snp_results_header" and "snp_results_detail" tables in pfacts003_test database 

Author: Curt Hansen
Created: 2012.08.05
Modified: 
'''

import pf_connect as db
import sys

def askForConfirm():
    while True:
        ans = raw_input('Warning!! This will delete all tables/views related to results. Do you want to continue? (y/n): ')
        if ans == 'y':
            break
        elif ans == 'n':
            sys.exit('Program aborted at user request!\n')
        else:
            print "You must enter 'y' or 'n'!"


def main():
    askForConfirm()
    
    dbConn=db.connect_to_server('pfacts003_test','webuser')
    dbCur=db.get_cursor(dbConn)
    dbConn.set_isolation_level(0)

    try:
        dbCur.execute("DROP VIEW snp_results_detail_avg;")
        print "\nPrevious version of 'snp_results_detail_avg' deleted."
    except:
        print "\nView 'snp_results_detail_avg' does not exist and will be created."
    try:
        dbCur.execute("DROP TABLE snp_results_detail;")
        print "Previous version of 'snp_results_detail' deleted."
    except:
        print "Table 'snp_results_detail' does not exist and will be created."
    try:
        dbCur.execute("DROP TABLE snp_results_header;")
        print "Previous version of 'snp_results_header' deleted."
    except:
        print "Table 'snp_results_header' does not exist and will be created."

    dbCur.execute("CREATE TABLE snp_results_header ( \
          id serial PRIMARY KEY, \
          time_stamp timestamp UNIQUE, \
          dataset_used varchar \
          );")
    print "Created new or replacement 'snp_results_header' table in database."

    dbCur.execute("CREATE TABLE snp_results_detail ( \
      id serial PRIMARY KEY, \
      time_stamp timestamp REFERENCES snp_results_header(time_stamp), \
      model_used varchar, \
      fold_no integer, \
      seq_no integer, \
      val_recall decimal(4,3), \
      val_precision decimal(4,3) \
      );")
    print "Created new or replacement 'snp_results_detail' table in database."

    dbCur.execute("CREATE VIEW snp_results_detail_avg AS \
      SELECT time_stamp, model_used, seq_no, AVG(val_recall) as recall, AVG(val_precision) as precision \
      FROM snp_results_detail \
      GROUP BY time_stamp, model_used, seq_no \
      ORDER BY time_stamp ASC, model_used ASC, seq_no ASC;")
    print "Created new or replacement 'snp_results_detail_avg' view in database.\n"

    dbConn.commit()
    dbConn.close()


if __name__=='__main__':
  main()
