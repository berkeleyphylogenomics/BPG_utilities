#!/usr/bin/python

import os
import sys
import re
import MySQLdb
from connect_db import *



def main():
    linepatt = re.compile('Book (.+?) with (\d+) sequences')

    handle = open(sys.argv[1], 'r')

    db = MySQLdb.connect(db=phylofacts_db_name, host=phylofacts_db_host, \
                         user=phylofacts_db_user, passwd=phylofacts_db_passwd)
    cur = db.cursor(MySQLdb.cursors.DictCursor)

    for line in handle:
        match = linepatt.match(line)
        
        book_acc = match.group(1)
        seq_count = int(match.group(2))
        
        if seq_count <= 1000:
            thedir = '/home/bpg/pfacts/'
        
            if book_acc[0:3] == 'bpg':
                dir1 = book_acc[0:6]
                thedir += '%s/%s/user' % (dir1, book_acc)
        
            else:
                thedir += book_acc
        
                # get best_gathering_method_across_sw
                # for this book from the db, so we can look up the
                # directory to put the NJ tree in
        
                sql = 'select best_gathering_method_across_sw from book where scopid_or_bpgid="%s";' % book_acc
                cur.execute(sql)

                for row in cur:
                    treedir = row['best_gathering_method_across_sw']
        
                thedir += '/' + treedir

            # get alignment file name and output tree file name
            alninfilename = thedir + '/' + book_acc + '_nr100.a2m'
            treeoutfilename = thedir + '/' + book_acc + '.nj'

            print '%s [%d], building NJ tree...' % (thedir, seq_count),
            sys.stdout.flush()

            if not os.path.exists(alninfilename):
                print 'ERROR: alignment file does not exist.'

            elif os.path.exists(treeoutfilename):
                print 'OOPS: tree file already exists.'

            else:
                # print temporary alignment, build tree        
        

                print 'done.'
    
    handle.close()



if __name__ == '__main__':
    main()
    
