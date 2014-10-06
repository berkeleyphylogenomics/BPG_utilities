#!/usr/bin/env python

import os, sys, re, pg, pgdb, tempfile
from optparse import OptionParser
from Bio import SeqIO
from random import randint
from pfacts003.utils.credentials import get_credentials

def match_species(species_re, record):
  fields = record.id.split('|')
  if len(fields) >= 3:
    id_components = fields[2].split('_')
    if len(id_components) >= 2:
      m = species_re.match(id_components[1])
      if m != None:
        return True
  return False

def sql_escape(str):
  escaped_str = str.replace("\\", "")
  escaped_str = escaped_str.replace("'", "\\'")
  return escaped_str
      
def main():
  # parse command line options
  usage = "%prog [options] fasta_file_to_cluster"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-s", "--candidate_seed_species", dest="species_str",
                default="",
                help="Comma-separated mnemonics of species which may be seeds")
  (options, args) = opt_parser.parse_args()
  if len(args) != 1:
    opt_parser.error('Incorrect number of arguments')
  fasta_file = args[0]
  if (not os.path.exists(fasta_file)):
    opt_parser.error('fasta file %s not found' % args[0])
  species_list = [species for species in options.species_str.upper().split(',')
                  if species != '']
  if len(species_list) > 0:
    species_re = re.compile('|'.join(species_list))
    is_desired_species = (lambda record: match_species(species_re, record))
  else:
    is_desired_species = (lambda record: True)
  flock_password = get_credentials('flock_user')
  connection = pgdb.connect("db:flock_seeds:flock_user:%s" % flock_password)
#  connection = pgdb.connect("db:flock_seeds::")
  cursor = connection.cursor()
  handle, random_path = tempfile.mkstemp()
  need_random_name = True
  while need_random_name:
    random_name = os.path.split(random_path)[1]
    os.close(handle)
    os.unlink(random_path)
    sql = """CREATE TABLE %s (
              sequence_key VARCHAR(255) PRIMARY KEY,
              length INTEGER,
              is_seed BOOLEAN,
              has_been_clustered BOOLEAN,
              random_integer INTEGER
              )""" % random_name
    try:
      cursor.execute(sql)
      connection.commit()
      need_random_name = False
    except pg.DatabaseError:
      pass
  print random_name
  sql = """CREATE INDEX %s_is_seed
            ON %s(is_seed)""" % (random_name, random_name)
  cursor.execute(sql)
  connection.commit()
  sql = """CREATE INDEX %s_has_been_clustered
            ON %s(is_seed)""" % (random_name, random_name)
  cursor.execute(sql)
  connection.commit()
  sql = """CREATE INDEX %s_random_integer
            ON %s(random_integer)""" % (random_name, random_name)
  cursor.execute(sql)
  connection.commit()
  num_inserted_records = 0
  num_records = 0
  f = open(fasta_file, "rU")
  seq_iterator = SeqIO.parse(f, "fasta")
  for record in seq_iterator:
    num_records += 1
    if is_desired_species(record):
      sql = "SELECT sequence_key FROM %s WHERE sequence_key = '%s'" \
              % (random_name, sql_escape(record.id))
      cursor.execute(sql)
      connection.commit()
      if cursor.rowcount > 0:
        # There is a duplicate entry, which we do not expect, so report it
        print "Duplicate entry for sequence key %s" % sql_escape(record.id)
      else:
        sql = """INSERT INTO %s
                  ( sequence_key,
                    length,
                    is_seed,
                    has_been_clustered,
                    random_integer
                  )
                  VALUES
                  ( '%s',
                    %d,
                    false,
                    false,
                    %d
                  )""" % (random_name,
                          sql_escape(record.id),
                          len(record.seq),
                          randint(1,10000000))
        cursor.execute(sql)
        connection.commit()
        num_inserted_records += 1
  connection.commit()
  cursor.close()
  f.close()
  connection.close()
  print "Inserted %d of %d records" % (num_inserted_records, num_records)

if __name__ == '__main__':
  main()
