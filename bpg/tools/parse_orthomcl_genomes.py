#!/usr/bin/env python

import urllib2, re, MySQLdb

from connect_db import *

db = MySQLdb.connect( db=phylofacts_db_name, host=phylofacts_db_host, \
                      user=phylofacts_db_user, passwd=phylofacts_db_passwd )
cur = db.cursor( MySQLdb.cursors.DictCursor )

def get_taxon_of_abbrev():
  url = "http://www.orthomcl.org/cgi-bin/OrthoMclWeb.cgi?rm=genome"
  cell_re = re.compile('<td align="center">(.*)</td>')
  italic_re = re.compile('<i>(.*)</i>')
  f = urllib2.urlopen(url)
  line = f.readline()
  while line and line.rstrip() != '<!-- BEGIN CONTENT -->':
    line = f.readline()
  while line and line.rstrip()[0:6] != '<table':
    line = f.readline()
  while line and line.strip()[0:6] != '<table':
    line = f.readline()
  while line and line.strip() != '</tr>':
    line = f.readline()
  begin_table_row, num_cell, abbrev_cell, name_cell, end_table_row = range(5)
  def next_state(state):
    return (state + 1) % 5
  state = begin_table_row
  abbrev = ''
  name = ''
  name_of_abbrev = {}
  lines = f.readlines()
  for line in lines:
    line = line.strip()
    if state == begin_table_row:
      if line and line[0:4] == '<tr>':
        state = next_state(state)
    elif state == num_cell:
      m = cell_re.match(line)
      if m:
        state = next_state(state)
    elif state == abbrev_cell:
      m = cell_re.match(line)
      if m:
        abbrev = m.group(1)
        state = next_state(state)
    elif state == name_cell:
      m = cell_re.match(line)
      if m:
        italicized_name = m.group(1)
        m = italic_re.match(italicized_name)
        if m:
          name = m.group(1)
          # The names of these bacteria with strain specifications are not
          # found in the database, so leave off the strain specifications
          if name[0:6] == 'Vibrio' or name[0:11] == 'Escherichia':
            name = ' '.join(name.split()[0:2])
          name_of_abbrev[abbrev] = name
          state = next_state(state)
    elif state == end_table_row:
      if line and line[0:5] == '</tr>':
        state = next_state(state)
    elif line[0:8] == '</table>':
      break
  names = ['"%s"' % name_of_abbrev[abbrev] for abbrev in name_of_abbrev.keys()]
  names_clause = ','.join(names)
  sql = """SELECT id, scientific_name
            FROM  ncbi_taxonomy
            WHERE scientific_name IN (%s)""" % names_clause
  n_rows = cur.execute(sql)
  rows = cur.fetchall()
  taxon_id_of_name = {}
  for row in rows:
    taxon_id_of_name[row['scientific_name']] = row['id']
  taxon_of_abbrev = {}
  for abbrev in name_of_abbrev.keys():
    name = name_of_abbrev[abbrev]
    if name in taxon_id_of_name:
      taxon = taxon_id_of_name[name]
      taxon_of_abbrev[abbrev] = taxon
  return taxon_of_abbrev

def main():
  taxon_of_abbrev = get_taxon_of_abbrev()
  print '"Abbreviation","Taxon Id"'
  for abbrev in taxon_of_abbrev:
    taxon = taxon_of_abbrev[abbrev]
    print '%s,%s' % (abbrev, taxon)



if __name__ == '__main__':
  main()
