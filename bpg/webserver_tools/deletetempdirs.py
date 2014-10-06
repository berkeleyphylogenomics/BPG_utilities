#!/usr/bin/python

import time
import glob
import os

tempdirpath  = '/clusterfs/ohana/software/webserver/temp/'
days_to_kill = 30 # delete temporary dirs > 30 days old

current_time = time.time()

existing_dirs = glob.glob(tempdirpath + '*')

for mydir in existing_dirs:
    creation_time = os.path.getctime(mydir)

    age_days = (current_time - creation_time) / 86400 # secs in a day

    # delete temp dir if it's too old #
    if age_days > days_to_kill:
        print '[age: %d days] deleting %s' % (age_days, mydir)
        cmd = 'rm -fr %s' % mydir
        os.system(cmd)

