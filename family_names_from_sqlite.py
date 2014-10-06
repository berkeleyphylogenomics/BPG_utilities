#!/usr/bin/python
'''
Rename families based on names stored (include different naming schemes and manually curated names) in an sqlite database.
'''
import os
import sys
from pysqlite2 import dbapi2 as sqlite3


DB_LOCATION = '/clusterfs/ohana/external/sqlite_databases/family_names.db'
#Adding naming order here, e.g. if the order preference is 2, then 1, then 3
#then do [2,1,3]
NAME_ORDER = [1] 

def add_related(name_str):
    """ If the word 'related' is not in the name, add it"""
    if not 'related' in name_str:
        return ('%s (related)' % name_str)
    
def get_db_cursor():
    """Get a cursor for the database"""
    connection = sqlite3.connect(DB_LOCATION)
    cursor = connection.cursor()
    cursor.execute('pragma foreign_keys=ON')
    return cursor, connection

def main(bpg_accession, family_type):
    """Gets a family name from the sqlite database.

    1. * (Changed to NAME_ORDER in this file) * Checks the name_order file to figure out the order (based on naming technique)
    in which names will be chosen. Converts this to a list.
    2. Connect to the database and query with the accession and name_order.
    3. Return the first value (takes care of the name order if multiple values are available)
    Add the word 'related' if needed.
    """
    name_order = NAME_ORDER
    #2.
    cursor, connection = get_db_cursor()
    if len(name_order) > 1:
        name_order = 'in (%s)' % (','.join(name_order))
    else:
        name_order = '= %s' % name_order[0]
    accession = (bpg_accession,)
    cursor.execute('''select name from family_name
    where family_bpg_accession=? and naming_method_id %s;''' % name_order, accession)
    #3.
    name  = cursor.fetchone()[0]
    if family_type == 'C':
        return name
    return add_related(name)

def get_family_name(bpg_accession, family_type):
    '''Just a convenience function that returns an empty string if there are any errors.
    '''
    try:
        return main(bpg_accession, family_type)
    except:
        return ''
    
if __name__ == "__main__":
    id_val = sys.argv[1]
    type_val = sys.argv[2]
    print get_family_name(id_val, type_val)
