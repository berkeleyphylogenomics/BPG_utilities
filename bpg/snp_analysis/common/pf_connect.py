'''
This is a module of commands used to access/query/update the pfacts003 database by the SNP analysis programs.
Author: Curt Hansen
Created: 3 May 2012 (rewrite of lost program)
Modified:
'''

import psycopg2
import psycopg2.extras
from pfacts003.utils.credentials import get_credentials


def connect_to_server(DB_NAME,USER):
    """
    Connects to postgres database and returns the cursor.
    """
    PWD = get_credentials(USER)
    conn = psycopg2.connect("dbname='%s' user='%s' host='db1' password='%s'" % (DB_NAME, USER, PWD))
    return conn


def get_cursor(conn):
    """
    Returns cursor.
    """
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    return cur


def get_query_list(cur,query,values=None):
    """
    Returns list.
    """
    if values==None:
        cur.execute(query)
    else:
        cur.execute(query,values)
    result = cur.fetchall()
    return result


def get_query_result(cur,query,values=None):
    """
    Returns single result.
    """
    if values==None:
        cur.execute(query)
    else:
        cur.execute(query,values)
    result = cur.fetchone()
    return result

