#!/clusterfs/ohana/software/bin/python2.6 -Wignore::DeprecationWarning
'''
This script deletes current versions of two pfacts003 tables for INTREPID scores and recreates them before a new computation in a separate program.

Inputs: None.
Outputs: "snp_rost_intrepid_header" and "snp_rost_intrepid_detail" tables in pfacts003_test database 

Author: Curt Hansen
Created: May 11, 2012
Modified: 
'''

import pf_connect as db

def main():
    dbConn=db.connect_to_server('pfacts003_test','webuser')
    dbCur=db.get_cursor(dbConn)
    dbConn.set_isolation_level(0)

    try:
        dbCur.execute("DROP TABLE snp_rost_intrepid_header;")
        print "\nPrevious version of 'snp_rost_intrepid_header' deleted."
    except:
        print "\nTable 'snp_rost_intrepid_header' does not exist and will be created."
    try:
        dbCur.execute("DROP TABLE snp_rost_intrepid_detail;")
        print "Previous version of 'snp_rost_intrepid_detail' deleted."
    except:
        print "Table 'snp_rost_intrepid_detail' does not exist and will be created."

    dbCur.execute("CREATE TABLE snp_rost_intrepid_header \
          (id serial PRIMARY KEY, protaccid varchar, familyid integer, \
          familytype char(1), msaok boolean, treeok boolean, idmapok boolean, scoresok boolean, familyscore integer, \
          numresscored integer, starttime varchar);")
    print "Created new or replacement 'snp_rost_intrepid_header' table in database."

    dbCur.execute("CREATE TABLE snp_rost_intrepid_detail \
      (id serial PRIMARY KEY, protaccid varchar, familyid integer, position integer,resseq integer,\
      chain varchar,residue char(1),icode char(1),scoreimpmap float,scorejsmap float,scoreremap float,\
      scoreglobalmap float,scoreglobalksmap float,scoreglobalremap float);")
    print "Created new or replacement 'snp_rost_intrepid_detail' table in database.\n"

    dbConn.commit()
    dbConn.close()


if __name__=='__main__':
  main()
