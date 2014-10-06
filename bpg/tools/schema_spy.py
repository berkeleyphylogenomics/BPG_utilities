#!/bin/env python

"""Gather database Metadata and generate Schema Spy data

    The Schema Spy utility is useful to create a browsable HTML schema
and graphical representation of relationship between tables.

    At the time of this writing, it would preferred that this would
run as a qsub job instead of running directly from Ohana. Therefore, we
submit this job and wait for the results to be printed.

    Since this requires a database connection, a password is used. It
is not stored in this script, but it retrieved and placed within the
qsub intermediate bash file. This is in the current working directory
and should, therefore, be removed immediately.

    The final Schema Spy product is in the 'schema_spy_output'
sub-directory in the current working directory. Qsub output and error
files are also placed there when completed. Monitor job status with the
'qstat' command.
"""

import os
from optparse import OptionParser
import subprocess
import sys

try:
    from pfacts003.utils.credentials import get_credentials
except ImportError:
    print """
    
    I couldn't import credentials. Are you sure you set up the
        environment? Your choices are production, staging and development.
    """
    sys.exit(1)


USER='bpg_user'
password = get_credentials(USER)

if not password:
    print "Could not get password."
    sys.exit(1)

def create_submission(working_dir, dirname='schema_spy_output'):

    contents = """#!/bin/bash
#PBS -e %(working_dir)s/schema_spy_error.log
#PBS -o %(working_dir)s/schema_spy_output.log
#PBS -N schema_spy

# WARNING!
# The bpg password had been retrieved and is included below. This
# file should therefore be removed immediately.

/usr/bin/java -jar /clusterfs/ohana/software/bin/schemaSpy_4.1.1.jar -t pgsql -host db -db pfacts003_test -s public -u %(user)s -p %(password)s -dp /clusterfs/ohana/software/bin/jdbc/postgresql-8.4-701.jdbc3.jar -o %(working_dir)s/%(dirname)s

""" % {'working_dir': working_dir,
       'dirname': dirname,
       'user': USER,
       'password': password}

    filename = os.path.join(working_dir, './schema_spy_qsub.sh')
    f = open(filename, 'w')
    f.write(contents)
    f.close()

    print "Submitting Schema Spy request to queue."
    subprocess.call(['qsub', filename, '-q', 'library'])
    print "\nUse the command 'qstat' to monitor..."


if __name__ == '__main__':
    parser = OptionParser(version='%prog 0.1')
    parser.add_option('-o', '--outdir', dest='outdir',
        help="Directory to write output",
        default=False)
    parser.add_option('-n', '--dirname', dest='dirname',
        help="Name of directory to write",
        default=False)

    (options, args) = parser.parse_args()

    if options.outdir:
        outdir = options.outdir
    else:
        outdir = os.getcwd()

    if not options.dirname:
        dirname = "schema_spy_output"
    else:
        dirname = options.dirname

    outdir = os.path.normpath(outdir)
    schema_out = os.path.join(outdir, 'schema_spy_output')
    try: 
        os.makedirs(schema_out)
    except OSError:
        print """
        
        I can't make the '%s' directory. Does it already exist?

        """ % schema_out
        sys.exit(1)

    create_submission(outdir)
