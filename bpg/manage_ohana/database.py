#!/bin/env python

"""Backup Ohana module

This module is used as a library (import), or a silent cron script only. For
cron script use assistance, execute this module with --help.

For interactive use, see the manage_ohana.py master script.
"""

from datetime import datetime, timedelta
import os
import sys
import subprocess
import tarfile

from optparse import OptionParser
# There is no psycopg2 on host 'db'
try:
    from psycopg2 import connect
except ImportError:
    pass

try:
    from bpg.manage_ohana.shared import current_date, stamp, logging, press_enter,\
     get_self_executable
except ImportError:
    sys.stderr.write("Activate an environment first.\n")
    sys.exit(1)

# Configurable static constants
DB_ROOT = '/clusterfs/ohana/software/db'
EXCLUDE_DB = ('template0', 'template1', 'pfacts003')
SCHEMA_ONLY_DB = ('postgres',)
FULL, INDIVIDUAL_DB, ROWS = xrange(3)
BACKUP_TYPES = (FULL, INDIVIDUAL_DB, ROWS)

# Dynamically configured global variables
dirname = "db_backup_%s" % stamp
db_root = os.path.join(DB_ROOT, dirname)


def local_logs(days=1, verbosity=False):

    """Helper function to copy postgres log files to the local $HOME directory

    This function will use the sudo permissions given to the postgres user to
    copy all log files from postgres to a subdirectory in the local $HOME
    directory. This is done because, on the host 'db', only the /home
    filesystem is writable/readable.

    This process *must* be run only on the host 'db' since: 1) that's the only
    location that the sudo permissions are granted, and 2) that's the only
    location where the postgres logs (that we're interested in) are located.

    This is an intermediate step to retrieving the logs, renaming them, storing
    a local copy, and parsing them for errors.
    """

    days = int(days)
    if days > 7:
        days = 1
    # Expands to whatever user is currently logged in
    home_directory = os.path.expanduser('~')

    for d in xrange(1, days+1):
        # Previous days logs (not today)
        day = datetime.today() - timedelta(days=d)
        source = '/var/lib/pgsql/data/pg_log/postgresql-%s.log' %\
                                                       day.strftime('%a')
        destination = os.path.normpath('%s/postgresql-%s.log' % (home_directory,
            day.strftime('%a')))

        process = subprocess.Popen(['sudo', '-u', 'postgres', 'tail', '-800000',
            source],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        stdoutdata, stderrdata = process.communicate()
        if process.returncode == 0:
            logfile = open(destination, 'w')
            logfile.write(stdoutdata)
            logfile.close()
            if verbosity:
                sys.stdout.write("Grabbed log %s.\n\n" % destination)
        else:
            sys.stderr.write(stderrdata)

        # HERE - clean-up and do wc on logs

def copy_logs(verbosity=False):

    print get_self_executable(__file__)
    press_enter()
    
    #retcode = subprocess.call(['pg_dump', '-h', 'db',
    #    '--file=%s' % filepath, db, '-o'])
    #if not retcode == 0:
    #    logging.error('%s: Dump failed!' % datetime.now())
    ##    if verbosity:
    #    sys.stderr.write("Dump failed!")

def take_db_backup(verbosity=False):

    if verbosity:
        print "Starting database backup..."

    # Create directory to put backup
    try:
        os.makedirs(db_root)
    except OSError:
        error = "ERROR: Unable to make the backup directory '%s'" % db_root
        sys.stderr.write('\n%s...\n\n' % error)

        if verbosity:
            press_enter()
        logging.error('%s: %s'% (datetime.now(), error))
        return False

    # Connect to DB
    conn = connect("dbname=postgres host=db")

    # Decide on type of backup
    backup_type = current_date.day % len(BACKUP_TYPES)

    backup_type = INDIVIDUAL_DB # HERE: hacked

    if backup_type == FULL:
        if verbosity:
            print "Full backup"
    else:
        # If INDIVIDUAL_DB or ROWS, iterate over all databases

        curs = conn.cursor()
        curs.execute("SELECT datname from pg_database")
        rows = curs.fetchall()
        databases = [db[0] for db in rows]
        databases = filter(lambda db: db not in EXCLUDE_DB, databases)
        if verbosity:
            print "Getting databases to backup..."
        logging.info('%s: Beginning backup of databases' % datetime.now())
        for db in databases:

            if backup_type == INDIVIDUAL_DB:
                filename = "%s_%s.sql" % (db, stamp)
                filepath = os.path.join(db_root, filename)


                # Dumping database
                if verbosity:
                    print "\tDumping %s database..." % db
                logging.info('%s: Dumping %s database' % (datetime.now(), db))
                begin_time = datetime.now()
                retcode = subprocess.call(['pg_dump', '-h', 'db',
                    '--file=%s' % filepath, db, '-o'])
                if not retcode == 0:
                    logging.error('%s: Dump failed!' % datetime.now())
                    if verbosity:
                        sys.stderr.write("Dump failed!")
                end_time = datetime.now()
                diff = (end_time - begin_time).seconds
                if verbosity:
                    print "\tDump complete."
                logging.info('%s: Completed in %d seconds' % (datetime.now(),
                    diff))

                # Create compressed tar file from original
                if verbosity:
                    print "\tCompressing %s database..." % db
                os.chdir(db_root)
                compressed_filename = "%s_%s.sql.tar.gz" % (db, stamp)
                tar = tarfile.open(compressed_filename, 'w:gz')
                tar.add(filename)
                tar.close()

                # At least on case, from the command line, we've seen the
                # compressed file not created (if it's too small)
                if os.path.exists(compressed_filename):
                    os.remove(filename)

                if verbosity:
                    print "\tCompression complete..."

            if backup_type == ROWS:
                pass
        del curs

        logging.info('%s: Completed backup of databases' % datetime.now())
        if verbosity:
            print "Completed backup of databases..."
            press_enter()


def main():
    """Options are passed and start_qa is called"""

    parser = OptionParser(version='%prog 0.1')
    parser.add_option('-s', '--scriptable', dest='scriptable',
        action="store_true",
        help="For cron job scripting - no user interaction is needed",
        default=False)
    parser.add_option('-l', '--locallogs', dest='locallogs',
        action="store_true",
        help="Copy logs from postgres to the local user directory (from db)",
        default=False)
    parser.add_option('-d', '--days', dest='days',
        help="Number of days of previous logs to grab.",
        default=1)

    (options, args) = parser.parse_args()

    if not options.scriptable:
        print __doc__
    else:
        if options.locallogs:
            print "Logs for ", options.days
            local_logs(days=options.days, verbosity=True)
        else:
            take_db_backup(verbosity=False)

if __name__ == '__main__':
    main()
