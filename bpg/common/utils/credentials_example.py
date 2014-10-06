#!/usr/bin/env python

import psycopg2
import psycopg2.extras
from pfacts003.utils.credentials import get_credentials

def main():
  # Get the password for webuser
  webuser_password = get_credentials('webuser')

  # Connect to the database as the webuser
  connection = psycopg2.connect(
    "dbname='%s' user='%s' host='db' password='%s'" %
    ('pfacts003_test', 'webuser', webuser_password))

  # Make a cursor, through which the database will give us results
  cur = connection.cursor(cursor_factory = psycopg2.extras.DictCursor)

  # Construct an SQL query string
  sql = """ SELECT tree.family_id, tree_node_name.name
              FROM family, tree, tree_node, tree_node_name
             WHERE tree.id = family.canonical_tree_id
               AND tree_node.tree_id = tree.id
               AND tree_node.left_id = 1
               AND tree_node_name.tree_node_id = tree_node.id
               AND tree_node_name.name LIKE '%Cellulase%';
        """

  # Execute the query on the database server
  cur.execute(sql)
  # Loop through the results
  for row in cur:
    # Format and print each result.
    # Newer versions of psycopg2 allow the columns to be accessed by name, but
    # these are not currently installed.
    family_accession = 'bpg%07d' % int(row[0])
    description = '"%s"' % row[1]
    print '%s: %s' % (family_accession, description)

if __name__ == '__main__':
  main()
