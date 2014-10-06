#!/bin/env python

"""Common functions used in the Manage Ohana package"""

import ConfigParser
from datetime import datetime
import logging
import os
import stat
import smtplib
from socket import gethostname
import subprocess
import sys


TEST = 'test'
PRODUCTION = 'prod'
environments = (TEST, PRODUCTION)


def _get_config_filename():
    """Find the correct ohana.cfg file

    Find the ohana.cfg file that is in the same directory as this
    executing shared.py file (or a default one if this function is
    imported via the python command line intereter).
    Requires global variable config to be defined an already
    initialized ConfigParser instance with ohana.cfg configuration
    """
    
    if not '__file__' in globals():
        # If ran from an python intereter, instead of a program file,
        # __file__ won't be defined. Using sensible defaulting in this
        # case.
        _default__file__ = ('/clusterfs/ohana/software/test/lib/',
            'python2.4/site-packages/bpg/manage_ohana/')
        sys.stderr.write(
            "Warning: No __file__, using '%s' as default." %
            _default__file)
        _this_file = _default__file__
    else:
        _this_file = __file__

    this_file = os.path.abspath(_this_file)
    (dir, tail) = os.path.split(this_file)
    return os.path.join(dir, 'ohana.cfg')


def setup_configuration():
    """Setup and read configuration file"""

    # Configuration_file = _get_config_filename()
    config = ConfigParser.SafeConfigParser()
    config.read(_get_config_filename())

    return config


def setup_environment(environment=TEST):
    """Setup Virtualenv environment"""

    config = setup_configuration()
    code_base = config.get('Default', 'code_base')

    activate_this = os.path.join(code_base, environment, 'bin/activate_this.py')
    execfile(activate_this, dict(__file__=activate_this))


def _set_permissions(pathname):
    """Set appropriate permissiosn on pathname

    There are certain files in the Berkeley Phylogenomics Group (BPG)
    that need to be shared with the entire group. These files should be
    owned by the 'bpg' group (not the 'basic' group) and should have
    read/write permissions for both users and groups. This function
    only augments the ownserhip and permissions of the files to add
    these conditions. No permissions are taken away (i.e., if you were
    silly enough to make a script that gave everyone execute
    permission, this function will not remove that permission for
    you).

    If errors were encountered (as we would see on a read only file
    system), then a boolean "False" is returned. True is returned
    otherwise.

    Requires global variable config to be defined an already
    initialized ConfigParser instance with ohana.cfg configuration
    """
    had_error = False

    gid = int(config.get('Default', 'bpg_group_id'))

    current_permissions = os.stat(pathname)[stat.ST_MODE]

    try:
        os.chown(pathname, -1, gid)
    except OSError:
        # Probably read-only file system like db or makana
        # So, can't change
        had_error = True

    # Augment current permissions so user/group has read/write
    new_permissions = current_permissions | stat.S_IRUSR | stat.S_IWUSR |  stat.S_IRGRP |  stat.S_IWGRP

    try:
        os.chmod(pathname, new_permissions)
    except OSError:
        # Probably read-only file system like db or makana
        # So, can't change
        had_error = True

    if os.path.isdir(pathname):
        try:
            os.chmod(pathname, new_permissions | stat.S_IXUSR | stat.S_IXGRP)
        except OSError:
            # Probably read-only file system like db or makana
            # So, can't change
            had_error = True


def _get_emergency_log_location():
    """When log file can't be written, get emergency fallback filename

    Whenever possible, we want to use the standard manage_ohana logs as
    defined in the configuration files. However, there are some
    circumstances where this file cannot be written to. For example, on
    both the 'db' and 'makana' nodes, most filepaths are mounted as read
    only (not read-write) for security reasons -- in those cases the only
    thing that you can write to is your $HOME directory.

    Requires global variable config to be defined an already initialized
    ConfigParser instance with ohana.cfg configuration
    """

    home_dir_name = config.get('Manage', 'home_dir_name')

    dir_name = os.path.join(os.path.expanduser('~'), home_dir_name)
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    # Make certain is readable by bpg group
    _set_permissions(dir_name)
    return dir_name


def get_validate_logfile():
    """Open log file to ensure we can read/write from it.

    Attempts are made to make certain the log file is in the bpg group
    and is readable by the gorup. If the file can't be opened or
    accessed (as any files not in your $HOME directory on systems like
    'makana' or 'db') a fall-back log file is created in your home
    directory.

    Requires global variable config to be defined an already initialized
    ConfigParser instance with ohana.cfg configuration
    """

    # Configurable static constants (from ohana.cfg)
    log_base = config.get('Default', 'log_base')

    # Calculate timestamp for log
    current_date = datetime.now()
    stamp = "%02d%02d%02d" % (current_date.year,
                              current_date.month,
                              current_date.day)

    log_filepath = os.path.join(log_base,
        'manage_ohana_%s.log' % stamp)


    # Several different people (and different user accounts)  may be
    # using manage_ohana. Make certain the log file  can be opened
    # (i.e., it's not on a real only file system and we have group
    # permissions to read/write to it).
    try:
        f = open(log_filepath, 'a+')
    except IOError, OSError: 
        if not _set_permissions(log_filepath):
            # Probably read only filesystem - have to bail
            # on this location regardless
            log_filepath = os.path.join(
                _get_emergency_log_location(),
                'manage_ohana_%s.log' % stamp)
            f = open(log_filepath, 'a+')
            f.close()

    return log_filepath


def press_enter():
    """Prompt and wait for Enter/Return keypress"""
    input = raw_input("Press the [Enter] key to continue> ")


def yes_to_continue(message="Do you wish to continue?"):
    """Prompt to continue. 
    
    Prompt user to continue. Only accepts a 'yes' or 'no' response and
    will not exit until one is chosen. If 'yes' is chosen, returns True.
    False is returned otherwise.
    """

    input = ""
    while input[:3].lower() not in ('yes', 'no'):
        input = raw_input("%s ('yes' or 'no')> " % message)

    if input[:3].lower() == 'yes': 
        return True
    else:
        return False


def get_hostname():
    """Get hostname (not fully qualified) for system"""

    return gethostname().split('.')[0]

def info(verbosity, msg):
    """Print info message to screen and log it

    This is a simple wrapper around logging.info(). However, if the
    verbosity variable is True, it will also report a message to the
    screen. Helpful if one is interacting with a user. These scripts
    are silent otherwise.
    """
    
    logging.info(msg)

    if verbosity:
        print "\n%s\n" % msg


def error(verbosity, msg, enter=True):
    """Write error message to stdout, log it, and prompt to continue

    This is a simple wrapper around logging.error(). However, if the
    verbosity variable is True, it will also write the error message to
    standard error.  Additionally, as errors should not pass by
    unnoticed, the user is prompted to press Enter. If verbosity is not
    True, errors will be written to logs, but will otherwise fail
    silently. This is for cron jobs and other non-interactive scripts.
    """

    logging.error(msg)
    if verbosity and msg != "":
        sys.stderr.write("\nERROR: %s\n" % msg)
        if enter:
            press_enter()


def send_email(sender, to_list, subject="", message="",
               debug=False, verbosity=False):
    info(verbosity, "Sending email to %s..." % ", ".join(to_list))

    msg = """From: %(sender)s\r
To: %(to)s\r
Subject: %(subject)s\r
\r
%(message)s
""" % {'sender': sender,
       'to': ", ".join(to_list),
       'subject': subject,
       'message': message}

    try:
        server = smtplib.SMTP('localhost')
        if debug:
            server.set_debuglevel(1)
        server.sendmail(sender, to_list, msg)
        server.quit()
    except Exception, e:
        error(verbosity, e)

    info(verbosity, "Email send sucessful...")


def get_self_executable(file_name):
    """Get the fully qualified filename of the executing file 

    This function is intended to be used by being called in a manner
    similar to the following:

    get_self_executable(os.path.abspath(__file__))

    The purpose of this is to find the fully qualified filename to
    execute on a remote machine via ssh. This assumes that the __file__
    in question is exactly the same as it is on this machine. That is a
    safe assumption in our environments were we use NFS filesystems
    mounted in the same manner across our nodes. 
    
    This "trick" (using this function to get an executable to ssh to
    another machine) is used to be able to execute commands on remote
    machines (like restarting a webserver on a different machine than
    where this code is currently running).
    """

    # To avoid permission errors, execute the python file, not the
    # compiled pyc file
    (head, tail) = os.path.split(file_name)

    (root, ext) = os.path.splitext(tail)
    if ext == '.pyc':
        ext = '.py'

    return os.path.join(head, "%s%s" % (root, ext))


def check_stale_dir(path, verbosity=False, message=None):
    """If this path exists, this should fail and possibly print error"""

    if message is None:
        message = """
   This directory (%s)
   is used during the code-push process as an intermediate step.
   It shouldn't exist during a normal push. If it still exists, then:

   a) someone else may be trying to push changes at the exact same time,
      or
   b) a previous push failed (more likely).

   Pushes are made to be 'as atomic as possible.' That is, no old
   information is removed until the new information is successfully in
   place. This way, we can:

   a) trouble-shoot what may have happened (RECOMMENDED),
   a) roll back (by manually restoring this original directory)
      (RECOMMENDED), or
   b) roll foward (by manually removing this directory)(NOT recommended
      unless you *really* know what you're doing.).

   Regardless, because this shouldn't exist, it will require manual
   intervention before we continue.
""" % path
   
    if os.path.exists(path):
        error(verbosity, "Stale directory (%s) found." % path, enter=False)
        if verbosity:
            sys.stderr.write(message)
        sys.exit(1)



def confirm_production(environment):
    """Print blazingly obvious message if this is production

    After obvious message is displayed, the user will be asked if
    they wish to continue. This routine will only return True if
    the user understands that it is production and still wishes to
    continue.
    """

    assert(environment in environments)

    production_message = """

PPPPP   RRRRR    OOOO   DDDDD   U    U   CCCC  TTTTT  IIIII   OOOO   N    N
P    P  R    R  O    O  D    D  U    U  C    C   T      I    O    O  NN   N
P    P  R    R  O    O  D    D  U    U  C        T      I    O    O  N N  N
PPPPP   RRRRR   O    O  D    D  U    U  C        T      I    O    O  N  N N
P       R  R    O    O  D    D  U    U  C        T      I    O    O  N   NN
P       R   R   O    O  D    D  U    U  C    C   T      I    O    O  N    N
P       R    R   OOOO   DDDDD    UUUU    CCCC    T    IIIII   OOOO   N    N

You are about to perform an act that will directly affect the production
systems on Ohana.

"""

    if environment == PRODUCTION:
        print production_message
        return yes_to_continue()
    else:
        return False


def _custom_command(title, command_list, verbosity=False):
    """Issue shell command, capture stdout/stderr, and log"""

    process = subprocess.Popen(command_list,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    stdoutdata, stderrdata = process.communicate()
    if process.returncode == 0:
        info(verbosity, "%s passed" % title)
        info(verbosity, stdoutdata)
    else:
        error(verbosity, "%s failed" % title)
        error(verbosity, stderrdata)


# Temporarily removed database function that follows.
#    else:
##    from psycopg2 import connect
##EXCLUDE_DB = ('template0', 'template1')
#
#def get_databases(verbosity=False):
#    """return list of databases as list"""
#    # Connect to DB
#    conn = connect("dbname=postgres host=db")
#
#    curs = conn.cursor()
#    curs.execute("SELECT datname from pg_database")
#    rows = curs.fetchall()
#    databases = [db[0] for db in rows]
#    databases = filter(lambda db: db not in EXCLUDE_DB, databases)
#
#    return databases


# Jump from this environment to the production one
setup_environment()

config = setup_configuration()
# Configure logging
logging.basicConfig(filename=get_validate_logfile(), level=logging.DEBUG)

# Get configurable parameters that can be imported from this module
code_base = config.get('Default', 'code_base')
#db_log_filepath = config.get('Default', 'db_log_filepath')
temp_base = config.get('Default', 'temp_base')
web_base = config.get('Default', 'web_base')
local_settings_file = config.get('Django', 'local_settings')
static_svn = config.get('Repository', 'static_svn')
sym_link_file = config.get('Repository', 'sym_link_file')
pfacts_svn = config.get('Repository', 'pfacts_svn')
admin_emails = config.get('Admin', 'admin_emails').split(',')
from_address = config.get('Admin', 'from_address')

