#!/usr/bin/env python

"""Manage Ohana Master Script

The Ohana cluster maintenance tasks, such as webserver file updates,
database backups, etc. are encapsulated here. This serves the purpose of
automating tasks and self-documenting procedures so we can keep more
organized and the system in good working order.

Menus or menu items can be added, changed, or modified as long as the following
structure is kept:

Each menu should be a dictionary with the 

"""

import ConfigParser
from datetime import datetime
import logging
from optparse import OptionParser
import os
import sys

from psycopg2 import connect

import shared
shared.setup_environment() # Bootstrap environment for remaining imports
#from bpg.manage_ohana.database import take_db_backup, copy_logs
#from bpg.manage_ohana.diskspace import check_disk_space
#from bpg.manage_ohana.families import check_families, BOOK_LIST_FILENAME
#from bpg.manage_ohana.vacuum import vacuum_all
from bpg.manage_ohana.repository import  TEST, PRODUCTION, link_repo, \
    push_all, push_bpg, push_bin, push_pfacts003, push_static
from bpg.manage_ohana.shared import config #, db_log_filepath
from bpg.manage_ohana.webserver import clean_webserver_temp, restart_webserver


# Configurable static constants (from ohana.cfg)
screen_height = config.getint('Default', 'screen_height')


# Menu data that drives the menu function
# The key 'Q' is reserved for quitting/exiting
#family_menu = {
#    'title': 'PhyloFacts Family Management Submenu',
#    'C': ('Check Families Disk Structure', check_families,
#                                    [BOOK_LIST_FILENAME], {'verbosity': True}),
#}
#database_menu = {
#    'title': 'Database Management Submenu',
#    'B': ('Take Database Backup', take_db_backup, [], {'verbosity': True}),
#    'L': ('Get Postgres Logs', copy_logs, [], {'verbosity': True}),
#    'V': ('Vacuum Database', vacuum_all, [], {'verbosity': True}),
#}
#cleaning_menu = {
#    'title': 'Clean Temporary Fiels/Directories',
#    'A': ('Clean using Apache', clean_webserver_temp, [],
#                   {'verbosity': True, 'apache': True}),
#    'U': ('Clean using your userid', clean_webserver_temp, [],
#                   {'verbosity': True, 'apache': False}),
#}
#volume_menu = {
#    'title': 'Disk Volume/Space Submenu',
#    'C': ('Check Disk Space', check_disk_space, [], {'verbosity': True}),
#}
staging_webserver_menu = {
    'title': 'Staging Webserver Submenu',
    'A': ('Export Everything', push_all, [TEST], {'verbosity': True}),
    'B': ('Export BPG Files', push_bpg, [TEST], {'verbosity': True}),
    'D': ('Export Django Files', push_pfacts003, [TEST], {'verbosity': True}),
    'I': ('Export bin Files', push_bin, [TEST], {'verbosity': True}),
    'L': ('Symbolically Link repository files', link_repo, [TEST], {'verbosity': True}),
    'R': ('Restart Web Server', restart_webserver, [], {'verbosity': True}),
    #'C': ('Clean Temporary Files/Directories', clean_webserver_temp, [],
    #               {'verbosity': True, 'apache': False}),
    #'C': ('Clean Temporary Files/Directories', 'menu', [],
    #                                            {'menu_data': cleaning_menu}),
    'S': ('Export Static Files', push_static, [TEST], {'verbosity': True}),

}
production_webserver_menu = {
    'title': 'Production Webserver Submenu',
    'A': ('Export Everything', push_all, [PRODUCTION], {'verbosity': True}),
    'I': ('Export bin Files', push_bin, [PRODUCTION], {'verbosity': True}),
    'L': ('Symbolically Link repository files', link_repo, [PRODUCTION], {'verbosity': True}),
#    #'C': ('Clean Temporary Files/Directories', 'menu', [],
#    #                                            {'menu_data': cleaning_menu}),
    'B': ('Export BPG Files', push_bpg, [PRODUCTION], {'verbosity': True}),
    'D': ('Export Django Files', push_pfacts003, [PRODUCTION], 
                                                          {'verbosity': True}),
    'S': ('Export Static Files', push_static, [PRODUCTION],
                                                          {'verbosity': True}),
    'R': ('Restart Web Server', restart_webserver, [], {'verbosity': True}),
}
webserver_menu = {
    'title': 'Webserver Main Menu',
    #'C': ('Clean Temporary Files/Directories', 'menu', [],
    #                                            {'menu_data': cleaning_menu}),
    'P': ('Update Production Web Server', 'menu', [],
                                     {'menu_data': production_webserver_menu}),
    'R': ('Restart Web Server', restart_webserver, [], {'verbosity': True}),
    'S': ('Update Staging Web Server', 'menu', [],
                                     {'menu_data': staging_webserver_menu}),
}
main_menu = {
    'title': 'Ohana Management Main menu',
    #'B': ('Book Management', 'menu', [], {'menu_data': family_menu}),
    #'D': ('Database Management', 'menu', [], {'menu_data': database_menu}),
    #'V': ('Disk volume/space management', 'menu', [],
    #                                               {'menu_data': volume_menu}),
    'W': ('Webserver Management', 'menu', [], {'menu_data': webserver_menu}),
}

class InvalidFunctionError(Exception):
    """Invalid function (entered as string) in menu structure

    You are allowed to enter a string that represents a function, instead of an
    actual function type, in the menu structure. However, the function should
    exist in globals() when the menu is attempted to be evaludated. Otherwise,
    this exception is raised. 

    Check the function name that you entered in the menu structure.
    """

    pass

def clear_screen():
    for row in xrange(screen_height):
        print 

def menu(menu_data):
    while True:
        clear_screen()

        menu_keys = menu_data.keys()
        menu_keys.sort()

        # Print menu title
        underscore = len(menu_data['title'])
        print "\t%s" % menu_data['title']
        print "\t", "-"*underscore

        for item in menu_keys:
            if item != 'title':
                print "\t%s) %s" % (item, menu_data[item][0])
        print "\tQ) Quit (or return to previous menu)"
        print "\n"

        choice = ""
        while len(choice) < 1:
            choice = raw_input("--> ").upper()

        if choice[0] == 'Q':
            return

        if choice in menu_data.keys():
        # Call associated menu function

            # Although not really necessary, these make the code more readable
            new_menu_data = menu_data[choice]
            new_function = new_menu_data[1]
            new_args = new_menu_data[2]
            new_kwargs = new_menu_data[3]

            if isinstance(new_function, basestring):
                try:
                    new_function = globals()[new_function]
                except KeyError:
                    raise InvalidFunctionError(\
    "\n\nFunction '%s' does not exist. Check menu structure." % (new_function))
            new_function(*new_args, **new_kwargs)


if __name__ == '__main__':
    menu(main_menu) 
