#!/usr/bin/env python

"""Webserver

This module is used as a library (import) only. For automated use, see the
manage_ohana.py master script.
"""

import datetime
import httplib
from optparse import OptionParser
import os
import shutil
import stat
import sys

import shared
shared.setup_environment() # Jump to real environment
from bpg.manage_ohana.shared import _custom_command, config, info, error,\
    get_hostname, get_self_executable, press_enter, send_email,\
    yes_to_continue, temp_base, from_address, admin_emails


def should_we_restart():
    """Confirm if user wishes to restart webserver.

    Because this is a function that interacts directly with the user, when we
    call the "restart_webserver", we always will pass 'verbosity=True'
    """

    if yes_to_continue("Would you like to restart the webserver now?"):
        restart_webserver(verbosity=True)


def _restart_helper(verbosity=False):
    """Helper function to restart webserver

    ADVANCED!  Do not use directly

    This function will use the sudo permissions given to check syntax and
    restart the webserver.

    This process *must* be run only on the host 'makana' since: 1) that's the
    only location that the sudo permissions are granted, and 2) that's the only
    location where the webserver is installed.

    This is an intermediate step  to automating the restarting of the
    webserver.
    """
    if get_hostname() != "makana":
        error(verbosity, "Webserver can only be restarted from makana.")
        sys.exit(1)

    # First check configuration
    _custom_command("Webserver configuration",
        ['/usr/bin/sudo', '/usr/sbin/apachectl', 'configtest'],
        verbosity=verbosity)

    # Then restart server
    _custom_command("Restart Server",
        ['/usr/bin/sudo', '/sbin/service', 'httpd', 'restart'],
        verbosity=verbosity)


def restart_webserver(verbosity=False):
    """Restart the webserver via an ssh call to host(makana)
   
    This function assumes that this file that you are reading is in the
    exact same location in the second host machine makana. In the BPG
    setup, these hosts share the remoted network file system, so this is
    the case.

    Essentially, one of these two commands is executed, depending upon
    the verbosity given:

    ssh -tt makana <this_file_name>.py --restart_helper
    ssh -tt makana <this_file_name>.py --restart_helper --verbose
    """

    options = ['ssh', '-tt', 'makana',
               get_self_executable(os.path.abspath(__file__)),
               '--restart_helper']

    if verbosity:
        options.append('--verbose')

    _custom_command("Calling restart on makana",
        options,
        verbosity=verbosity)

    send_email(from_address,
               admin_emails,
               subject="Webserver Restarted",
               message="FYI: Makana webserver was just restarted",
               verbosity=False) # No need to clutter/say we're sending email

    if verbosity:
        press_enter()


def clean_webserver_temp(verbosity=False, apache=True, days=7):
    """Remove the temporary files older than 7 (or a configured number) of days
    
    Unfortunately, the webserver temp directory has files owned by both the
    webserver (Apache) and the Job Daemon (usually ran under a particular
    username).

    This function performs two different types of actions depending upon the
    Apache boolean variable.

    If apache is True, the request is made via a Django URL (because the code
    has to execute as the webserver/apache user). This code tries to delete the
    entire tree with the shutils rmtree command. Please see the associated
    Django view for comments regarding the token.

    Otherwise, each particular file is walked and, if expired, there is an
    attempt to remove the file. This is done on an individual file/directory
    basis since some files will be owned by a 'third party' userid (the userid
    running the daemon).
    """

    current_date = datetime.datetime.now()

    if apache:
        info(verbosity, "Initiating webserver temp directory cleaning request")
        data = ["Removed files via Apache."]
        token = "%dlwx3%d0K91vmL%di" % (current_date.day, current_date.month,
                                        current_date.year)

        conn = httplib.HTTPConnection("makana-test.berkeley.edu")
        conn.request("GET", "/staff/clean_webserver_temp?token=%s" % token)
        response = conn.getresponse()

        if response.status != 200:
            error_msg = ("Error %d received " % response.status, 
                         "when trying to initate cleaning request!")
            error(verbosity, error_msg)
        else:
            info(verbosity, "Cleaning request successful.")

        data.append(response.read())

    else:
        # Remove files as logged in user instead of making request to Apache
        message = "Removing files owned by user: %s." % os.getlogin()
        data = [message]
        info(verbosity, message)

        for dirpath, dirnames, filenames in os.walk(TEMP_BASE):
            for file in filenames:
                # Get the date from each file's timestamp
                f = os.path.join(dirpath, file)
                timestamp = os.stat(f)[stat.ST_CTIME]
                file_date = datetime.datetime.fromtimestamp(timestamp)

                # If the file is too old, try to remove
                if file_date < current_date - datetime.timedelta(days=days):
                    try:
                        os.remove(f)
                        data.append("Deleting: %s." % f)
                    except Exception, e:
                        data.append("Can't delete: %s; %s" % (f, e))
                else:
                    # The file is still current and should be kept
                    data.append("Skipping: %s." % f)
        info(verbosity, "Cleaning complete")

    data = "\n".join(data)

    # This will be removed once we establish a solid working procedure that
    # has been tested.
    paranoid_file = open(os.path.expanduser('~/webserver_clean.out'), 'w')
    paranoid_file.write(data)
    paranoid_file.close()

    # This will also be removed when systems are guaranteed to be running well
    send_email('bpg.shailen@gmail.com', ['bpg.shailen@gmail.com', ],
               subject="Results of Webserver Temporary Cleaning", 
               message=data, verbosity=False)

    if verbosity:
        press_enter()



def main():
    """Options are passed and possibly silent scripting called"""

    parser = OptionParser(version='%prog 0.3')

    parser.add_option('-v', '--verbose', dest='verbose',
        action="store_true",
        help="Give verbose output (default is no output or interaction)",
        default=False)
    parser.add_option('-c', '--clean_temp', dest='clean_temp',
        action="store_true",
        help="Clean the temporary files from webserver")
    parser.add_option('-a', '--apache', dest='apache',
        action="store_true",
        default=False,
        help=\
       "Use apache userid when cleaning (use in conjunction with --clean_temp)")
    parser.add_option('-r', '--restart', dest='restart_webserver',
        action="store_true",
        help="Restart webserver (to make changes take affect)",
        default=False)
    parser.add_option('-x', '--restart_helper', dest='restart_helper',
        action="store_true",
        help="Advanced feature (not intended for general use)",
        default=False)

    (options, args) = parser.parse_args()

    if options.restart_helper:
        _restart_helper(verbosity=options.verbose)
        return

    if options.restart_webserver:
        restart_webserver(verbosity=options.verbose)
        return

    if options.clean_temp:
        clean_webserver_temp(verbosity=options.verbose, apache=options.apache)

    else:
        print __doc__

if __name__ == '__main__':
    main()
