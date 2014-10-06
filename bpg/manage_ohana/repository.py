#!/usr/bin/env python

"""Repository

For automated use, see the manage_ohana.py master script.

command lines are silent by default - and use default behaviors. For any interactivity, muse explicitly set using --verbose function. This is so that cron scripts are guaranteed to run without interference. 

To import, you can set the verbose level if you wish.
"""

from optparse import OptionParser
import os
import re
import shutil
import subprocess
import sys

import shared 
shared.setup_environment() # Bootstrap environment for remaining imports
from bpg.manage_ohana.shared import TEST, PRODUCTION, _get_config_filename,\
    admin_emails, check_stale_dir, code_base, config, error, from_address,\
    get_hostname, info, confirm_production, local_settings_file, press_enter,\
    pfacts_svn, send_email, static_svn, sym_link_file, web_base, yes_to_continue 
from bpg.manage_ohana.webserver import restart_webserver, \
    should_we_restart


def _export(base_dir, repository, repo_name, verbosity=False):
    """Repository export

    It is the caller of this helper function to ensure that base_dir exists.
    """
    assert(os.path.exists(base_dir))

    if get_hostname() != "ohana":
        error(verbosity, "Export can only be done from Ohana.")
        sys.exit(1)

    info(verbosity, "Exporting repository %s as %s..." % (repository, repo_name))
    if verbosity:
        print """
    This output is being buffered so I can read the version number.
    ....  I'm not stuck, just busy exporting files....
    """

    try:
        os.chdir(base_dir)
    except OSError:
        error(verbosity, "Unable to cd to: %s" % base_dir)
        return False
    
    process = subprocess.Popen(['svn', 'export', repository, repo_name],
        stdout=subprocess.PIPE)

    stdoutdata, stderrdata = process.communicate()
    if process.returncode == 0:
        info(verbosity, "\tFiles exported successfully")
    else:
        error(verbosity,
           "Error codes %d was returned from SVN export." % process.returncode)
        error(verbosity, stdoutdata, enter=False)
        return False

    # Grab the revision exported
    revision_number = re.search(r'^Exported revision (\d*)\.$',
        stdoutdata, re.MULTILINE).group(1)

    if revision_number is None:
        error(verbosity,"Revision number not successfully retrieved")
        return False

    del stdoutdata
    # Write revision to README.txt
    readme = os.path.join(base_dir, repo_name, 'README.txt')
    readme_file = open(readme, 'w')
    readme_file.write("Exported revision %s." % revision_number)
    readme_file.close()

    message = "Exported %s to %s." % (repository,
                                      os.path.join(base_dir, repo_name))
    send_email(from_address,
               admin_emails,
               subject="New push/export",
               message=message,
               verbosity=False)

    # If we got this far without an error, return success.
    return True




def export(base_dir, repository, repo_name, verbosity=False,
           has_local_settings=False):

    _OLD = '_old'
    _NEW = '_new'
    # Step 1 of 5: Preliminary checks
    if not os.path.exists(base_dir):
        info(verbosity, "%s doesn't exist. Trying to create." % base_dir)
        try:
            os.makedirs(base_dir)
        except OSError:
            error(verbosity, "Unable to create directory: %s" % base_dir)
            return False

    exported_orig = os.path.join(base_dir, repo_name)
    exported_old = os.path.join(base_dir, "%s%s" % (repo_name, _OLD))
    exported_new = os.path.join(base_dir, "%s%s" % (repo_name, _NEW))
    check_stale_dir(exported_old, verbosity=verbosity)
    check_stale_dir(exported_new, verbosity=verbosity)

    # Step 2 of 5: Export _new copy
    if not _export(base_dir, repository, "%s%s" % (repo_name, _NEW),
                      verbosity=verbosity):
        return False

    # Step 3 of 5 (optional): Copy local_settings.py
    if has_local_settings:
       local_filename = os.path.join(exported_orig, local_settings_file)
       local_new_filename = os.path.join(exported_new, local_settings_file)

       if os.path.exists(local_filename):
           info(verbosity, "\tCopying local settings.")
           shutil.copyfile(local_filename, local_new_filename) 
       else:
           error(verbosity, "No local settings. Creating empty file")
           open(local_new_filename, 'w').close()

    # Step 4 of 5: Put new directories in place
    # This wrapped in an blanket exception because any failure places server
    # in inconsistent state
    try:
        info(verbosity, "\tPutting changes in place: Renaming directories.")
        # Move original to old
        if os.path.exists(exported_orig):
            os.rename(exported_orig, exported_old)
        # Move new to original
        os.rename(exported_new, exported_orig)
    except:
        error(verbosity, 
            "Rename of directories failed. Webserver in unstable state.")
        sys.exit(1)

    # Step 5 of 5: Put new directories in place
    # Delete original (now called _old)
    if os.path.exists(exported_old):
        info(verbosity, "\tRemoving old files.")
        shutil.rmtree(exported_old)

    info(verbosity, "Update completed")
    if verbosity:
        press_enter()


def generic_push(environment, base_dir, repository, repo_name, verbosity=False,
                has_local_settings=False, ask_restart=True):

    info(verbosity, "Beginning push of %s on %s..." % (repo_name, environment))

    # If we're interacting with the user and this is production
    # We need to make certain they know this is production
    if verbosity and environment == PRODUCTION: 
        if not confirm_production(environment):
            return False

    export(base_dir, repository, repo_name, verbosity=verbosity,
           has_local_settings=has_local_settings)
           
    info(verbosity, "Push of %s on %s completed..." % (repo_name, environment))

    if verbosity and ask_restart: 
        should_we_restart()


def push_pfacts003(environment, verbosity=False, ask_restart=True):

    base_dir = os.path.join(code_base, environment,
        'lib/python2.4/site-packages/')

    generic_push(environment, base_dir,
                 os.path.join(pfacts_svn, 'TRUNK', 'pfacts003'),
                 'pfacts003', verbosity=verbosity,
                 has_local_settings=True, ask_restart=ask_restart)


def push_bpg(environment, verbosity=False, ask_restart=True):

    base_dir = os.path.join(code_base, environment,
        'lib/python2.4/site-packages/')

    generic_push(environment, base_dir,
                 os.path.join(pfacts_svn, 'TRUNK', 'bpg'),
                 'bpg', verbosity=verbosity, ask_restart=ask_restart)


def push_bin(environment, verbosity=False, ask_restart=True):

    base_dir = os.path.join(code_base, environment)

    generic_push(environment, base_dir,
                 os.path.join(pfacts_svn, 'TRUNK', 'bin'),
                 'repo_bin', verbosity=verbosity, ask_restart=ask_restart)


def push_static(environment, verbosity=False, ask_restart=False):

    base_dir = os.path.join(web_base, environment, 'pfacts003_static_files')

    # Static files aren't cached by Apache. No need to restart server
    # after these files are pushed
    generic_push(environment, base_dir,
                 os.path.join(static_svn, 'TRUNK'),
                 'static', verbosity=verbosity, ask_restart=ask_restart)

def link_repo(environment, verbosity=False, ask_restart=False):
    """Mapping of repository command line scripts into repo_bin"""
    #ERROR: Duplicate:
        #./bin/discern.py,
        #/bpg/discern/bin/discern.py
    #ERROR: Duplicate:
        #./bin/intrepid.pl,
        #./bpg/discern/bin/intrepid.pl
    #ERROR: Duplicate:
        #./bin/intrepid.py,
        #./bpg/intrepid/intrepid.py
    #ERROR: Duplicate:
        #./bin/parse_hmmpfam_results.py,
        #./bpg/common/parsers/hmm_parsers/parse_hmmpfam_results.py
    #ERROR: Duplicate:
        #./bin/print_uniprot_orthologs_from_taxa.py,
        #./bpg/orthologs/print_uniprot_orthologs_from_taxa.py

    base = os.path.normpath(os.path.join(code_base, environment))

    src_base = os.path.join(base, 'lib%spython2.4%ssite-packages' %\
        (os.path.sep, os.path.sep))

    if not os.path.exists(os.path.join(base, 'repo_bin')):
        error(verbosity, "Warning: repo_bin doesn't exist. Pushing now.")
        push_bin(environment, verbosity=verbosity, ask_restart=ask_restart)

    # This runs out of the test environment. Thus, the symbolic
    # link file used is that one. If you don't remember to push bpg
    # to test before symbolic linking, the link file uses is the
    # one previously in the test environment bpg.
    # TODO: Ask if changed sym_links and link in.
    sym_file = _get_config_filename().replace('ohana.cfg', sym_link_file)
    sym_links = open(sym_file, 'r')

    for rel_filepath in sym_links:
        rel_filepath = rel_filepath.strip()
        (dir, filename) = os.path.split(rel_filepath)

        # Build fully qualified names
        src = os.path.normpath(os.path.join(src_base, rel_filepath))
        link_name = os.path.join(base, 'repo_bin', filename)

        if os.path.exists(link_name):
            # Something already exists where the link should be
            if os.path.islink(link_name):
                os.unlink(link_name)
            else:
                error(verbosity, "Error: Can't link because '%s' already exists" % link_name)
                continue

        # Make certain file exists before linking
        if os.path.exists(src):
            os.symlink(src, link_name)
        else:
            error(verbosity,  "Warning: %s doesn't exist. Not linking." % src)


def push_all(environment, verbosity=False):
    """Convenience function..."""

    # Start with static files - least chance of impact on production during 
    # Skip any restarting until the end - it's better for our webserver and
    # user (and it's faster).
    push_static(environment, verbosity=verbosity, ask_restart=False)

    # Next bpg
    push_bpg(environment, verbosity=verbosity, ask_restart=False)

    # Finally Django
    push_pfacts003(environment, verbosity=verbosity, ask_restart=False)

    if verbosity:
        should_we_restart()


def main():
    """Options are passed and possibly silent scripting called"""

    parser = OptionParser(version='%prog 2.0')
    parser.add_option('-p', '--production', dest='production',
        action="store_true",
        help="To choose push to production",
        default=False)
    parser.add_option('-g', '--test', dest='test',
        action="store_true",
        help="To choose push to test",
        default=False)
    parser.add_option('-a', '--all', dest='push_all',
        action="store_true",
        help="Push all (django, bpg and static) to server")
    parser.add_option('-b', '--push_bpg', dest='push_bpg',
        action="store_true",
        help="Push bpg to server")
    parser.add_option('-d', '--push_pfacts003', dest='push_pfacts003',
        action="store_true",
        help="Push django to server")
    parser.add_option('-c', '--push_static', dest='push_static',
        action="store_true",
        help="Push static to server")
    parser.add_option('-v', '--verbose', dest='verbose',
        action="store_true",
        help="Yields verbose output (scriptable has priority)",
        default=False)

    (options, args) = parser.parse_args()

    # Determine Environment: TEST or PRODUCTION
    environment = TEST
    if options.test and options.production:
        # TEST wins in a tie. You must be exact to push to production
        environment = TEST
    elif options.production:
        environment = PRODUCTION

    # Read Environment and set-up VirtualEnv
    shared.setup_environment(environment)

    # Push all and stop
    if options.push_all:
        push_all(environment, verbosity=options.verbose)
        return

    # Push any option given
    if options.push_bpg:
        push_bpg(environment, verbosity=options.verbose)
    if options.push_pfacts003:
        push_pfacts003(environment, verbosity=options.verbose)
    if options.push_static:
        push_static(environment, verbosity=options.verbose)


if __name__ == '__main__':
    main()
