#!/usr/bin/python

from daemon import Daemon
import logging
import os
import re
import stat
import sys
import time
import glob


# Set up daemon logging
log_base = '/clusterfs/ohana/software/webserver/logs'
if not os.path.exists(log_base):
    os.makedirs(log_base)
DAEMON_LOG = os.path.join(log_base, 'daemon.log')
logging.basicConfig(filename=DAEMON_LOG,level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(message)s")


class RunJobsDaemon(Daemon):
    """RunJobsDaemon class for running webuser-induced jobs"""

    def run(self):

        job_prefix = '/clusterfs/ohana/software/webserver/submitted_jobs'
        work = '/clusterfs/ohana/software/webserver/temp'

        qsub_template = \
"""#!/bin/sh

#PBS -l nodes=1:ppn=1
#PBS -l walltime=60:00
#PBS -d %s

source /etc/profile.d/ohana.sh
# This is a temporary measure and needs permenantly fixed!
source /clusterfs/ohana/software/test/bin/activate
%s
"""

        allowed_job_patterns = [ re.compile(regstr) for regstr in (
            # Original legacy location
            r'^((?P<target>/clusterfs/ohana/software/webserver/incoming(_test)?/.*)\n)?(?P<exec>/clusterfs/ohana/software/webserver/(staging|production)/pfacts_django/pfacts003/queued/queued_run.py \w+( /clusterfs/ohana/software/webserver/(staging|production)(/[-.\w]+)*)*)$',
            # Experimental environment (to be phased in)
            r'^((?P<target>/clusterfs/ohana/software/webserver/incoming(_test)?/.*)\n)?(?P<exec>/clusterfs/ohana/software/(test|prod)/lib/python2.4/site-packages/pfacts003/queued/queued_run.py \w+ (/clusterfs/ohana/software/(test|prod)/lib/python2.4/site-packages/ /clusterfs/ohana/software/(test|prod)/lib/python2.4/)).*$',
        )]

        #
        # main server loop here
        #

        while True:

            newfiles = (
                filename \
                for filename in os.listdir(job_prefix) \
                if filename.endswith('.sh')
            )

            for newfile in newfiles:
                
                tmpdirname = newfile.split('.')[0]
                qsub_file_path = os.path.join(job_prefix, newfile)
                work_dir_path = os.path.join(work, tmpdirname)
                work_file_path = os.path.join(work_dir_path, 'qsub.sh')

                # make sure the script is okay
                handle = open(qsub_file_path, 'r')
                myscript = handle.read().strip()
                handle.close()
                os.unlink(handle.name)

                
                # I wish I had python 2.5 for the 'any' function
                match = None
                for pattern in allowed_job_patterns:
                    match = pattern.match(myscript)
                    if match:
                        break

                # move script to appropriate directory and run it
                # if it is okay
                if match:
                    # use an indicated path, or the default
                    target_dir_path = match.groupdict().get('target', False) or work_dir_path
                    target_file_path = os.path.join(target_dir_path, 'qsub.sh')

                    logging.info( \
                        "Moving script %s to %s" % \
                            (qsub_file_path, target_file_path))
                    try:
                        handle = open(target_file_path, 'w')
                        handle.write(qsub_template % (
                            target_dir_path, match.group('exec'),
                        ))
                        handle.close()
                        os.chmod(handle.name, stat.S_IRWXU | stat.S_IRWXG | stat.S_IROTH)
                        logging.info("Submitting script %s to queue." % qsub_file_path)
                        os.system('qsub -q web %s' % target_file_path)
                    except:
                        logging.error('lamezorz')
                        sys.exit(1)
                else:
                    # flag script as 'EVIL'
                    logging.error("Script %s did not pass validation." % qsub_file_path)
                    handle = open(qsub_file_path + '.EVIL', 'w')
                    handle.write(myscript)
                    handle.close()

            time.sleep(1)


if __name__ == "__main__":
    daemon = RunJobsDaemon('/clusterfs/ohana/software/webserver/runjobsd.pid')

    if len(sys.argv) == 2:
        if     'start' == sys.argv[1]:
            daemon.start()
        elif    'stop' == sys.argv[1]:
            daemon.stop()
        elif 'restart' == sys.argv[1]:
            daemon.restart()
        else:
            print "Unknown command"
            sys.exit(2)

        sys.exit(0)

    else:
        print "usage: %s start|stop|restart" % sys.argv[0]
        sys.exit(2)

