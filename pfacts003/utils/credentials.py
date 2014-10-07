#!/usr/bin/env python

"""Read and return Berkeley Phylogenomics Group credentials

    This module reads the configuration file pointed to by the variable
CREDENTIAL_FILE. This indicated file should *not* be included into the code
respository because it contains sensitive information (like what is needed for
database access).  Think of the CREDENTIAL_FILE, like a license file. The
software won't run without it but extra copies shouldn't be laying around.

    If the file needs to be created, it should be of the RFC 882 format
(http://tools.ietf.org/html/rfc822.html)(similar to a windows.ini file). Each
username listed in the file should be in its own section, indicated by
brackets. A single parameter, password, is expected in each section. An example
follows:

[flock]
password: sample_flock_pw

[webuser]
password: sample_webuser_pw

    The following two code examples (one in bash, the other in python), show
how this module may be used:

EXAMPLE
-------
Example for python, from inside of Django (see django_script.py for example of
how to access outside of Django):

from pfacts003.utils.credentials import get_credentials

password = get_credentials('the_username_you_are_using')

if not password:
    print "Could not get password."


EXAMPLE
-------
Example for command line/bash:

password=$(./credentials.py --username=the_username_you_are_using)

# Check if retrieval was successful
if [ "$?" -ne 0 ]; then
    echo "Could not get password."
else
    echo "Password: $password"
fi
"""

import ConfigParser
from optparse import OptionParser
import os


CREDENTIAL_FILE = '/clusterfs/ohana/software/db/bpg_credentials.cfg'
ERROR = 2


def get_credentials(username):
    config = ConfigParser.ConfigParser()
    config.read(CREDENTIAL_FILE)
    sections = config.sections()

    password = False

    if username in sections:
        try:
            password = config.get(username, 'password')
        except ConfigParser.NoOptionError:
            password = False

    return password


def main():
    """Options are passed and get_credentials called"""
    parser = OptionParser(version='%prog 0.1')
    parser.add_option('-u', '--username', dest='username',
        help="The username you wish to get credentials for")

    (options, args) = parser.parse_args()

    if options.username is not None:
        password = get_credentials(options.username)

    if password:
        print password
    else:
        os._exit(ERROR)


if __name__ == '__main__':
    main()
