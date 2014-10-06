#!/bin/env python

"""Vacuum Ohana module

This module is used as a library (import), or a silent cron script only. For
cron script use assistance, execute this module with --help.

For interactive use, see the manage_ohana.py master script.
"""

from datetime import datetime
import os
import sys
from optparse import OptionParser
import subprocess

from psycopg2 import connect

from database import dirname
from bpg.manage_ohana.shared import logging, press_enter, get_databases

# Configurable static constants
DB_ROOT = '/clusterfs/ohana/software/db'
EXCLUDE_TABLES = ('sql_features','sql_implementation_info','sql_languages',
    'sql_packages','sql_sizing','sql_sizing_profiles')

VACUUM_VERBOSE_OUTPUT = 'vacuum_verbose.out'


vacuum_output = os.path.normpath("%s/%s/%s" % (DB_ROOT, dirname,
    VACUUM_VERBOSE_OUTPUT))

print vacuum_output

def vacuum_database(db_name, verbosity=False):
    if verbosity:
        sys.stdout.write("Vacuuming %s...\n" % db_name)
    conn = connect("dbname=%s host=db" % db_name)

    if verbosity:
        sys.stdout.write("\tGetting tables to vacuum %s...\n" % db_name)
    curs = conn.cursor()
    curs.execute("SELECT tablename from pg_tables")
    rows = curs.fetchall()
    tables = [tbl[0] for tbl in rows]
    tables = filter(lambda tbl: tbl not in EXCLUDE_TABLES, tables)
    del curs

    for table in tables:
        if verbosity:
            "Vacuuming %s: %s" % (db_name, table)
        logging.info('%s: Begin Vacuuming %s:%s' % (datetime.now(), db_name,
            table))
        begin_time = datetime.now()

        # Using command line/not wrapping within transaction
        sql_statement = "VACUUM VERBOSE ANALYZE %s;" % table

        # Because of ident server errrors, this may fail intermittantly
        # Retries are often successful
        retries_left = 3
 
        while retries_left > 0:
            logging.info('%s: Begin Vacuuming %s' % (datetime.now(), table))
            process = subprocess.Popen(['/usr/bin/psql', '-h', 'db',
                db_name, '-c', sql_statement],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                )

            stdoutdata, stderrdata = process.communicate()
            if process.returncode == 0:
                retries_left = 0
            else:
                logging.warning('%s: Vacuuming %s failed. Retrying.' %\
                                                       (datetime.now(), table))
                retries_left = retries_left - 1
                if retries_left == 0:
                    logging.error('%s: Unable to vacuum %s.' % (datetime.now(),
                        table))

            # Save stdout and stderr data HERE - stdout still not handled
            vacuum_verbose = open(vacuum_output, 'a')
            vacuum_verbose.write(stdoutdata)
            vacuum_verbose.write(stderrdata)
            vacuum_verbose.close()

        end_time = datetime.now()
        diff = (end_time - begin_time).seconds
        logging.info('%s: Completed in %d seconds' % (datetime.now(),
                    diff))

def vacuum_all(verbosity=False):

    if verbosity:
        sys.stdout.write("Beginning vacuuming of databases...\n")
    logging.info('%s: Beginning vacuuming of databases' % datetime.now())

    databases = get_databases()
    for db in databases:
        vacuum_database(db, verbosity=False)
    logging.info('%s: Vacuuming complete.' % datetime.now())
    if verbosity:
        sys.stdout.write("Vacuuming complete.\n")
        press_enter()


def main():
    """Options are passed and start_qa is called"""

    parser = OptionParser(version='%prog 0.1')
    parser.add_option('-s', '--scriptable', dest='scriptable',
        action="store_true",
        help="For cron job scripting - no user interaction is needed",
        default=False)

    (options, args) = parser.parse_args()

    if not options.scriptable:
        print __doc__
    else:
        vacuum_all()

if __name__ == '__main__':
    main()
