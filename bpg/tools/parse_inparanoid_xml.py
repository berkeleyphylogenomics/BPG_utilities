#!/usr/bin/env python

import os, sys, MySQLdb, glob
from connect_db import *
from xml.sax import ContentHandler, make_parser, saxutils

class InParanoidHandler(saxutils.DefaultHandler):
  def __init__(self, filename):
    self.filename = filename
    self.currentCluster = 0
    self.currentBitScore = 0
    self.db = MySQLdb.connect( db=phylofacts_db_name, host=phylofacts_db_host, 
                               user=phylofacts_db_user,
                               passwd=phylofacts_db_passwd, 
                               charset='utf8' )
    self.cur = self.db.cursor( MySQLdb.cursors.DictCursor )
  def startElement(self, name, attrs):
    if name == 'CLUSTER':
      self.currentCluster = int(attrs.get('CLUSTERNO', '0'))
      self.currentBitScore = int(attrs.get('BITSCORE', '0'))
    elif name == 'GENE':
      fields_clause = "(FILENAME, CLUSTERNO, BITSCORE"
      values_clause = "('%s', %d, %d" % (self.filename, self.currentCluster,
                                        self.currentBitScore)
      geneid = attrs.get('GENEID', '')
      if geneid != '':
        fields_clause = fields_clause + ", GENEID"
        values_clause = values_clause + ", '%s' " % attrs['GENEID']
      protid = attrs.get('PROTID', '')
      if protid != '':
        fields_clause = fields_clause + ", PROTID"
        values_clause = values_clause + ", '%s'" % attrs['PROTID']
      score = attrs.get('SCORE', '')
      if score != '':
        fields_clause = fields_clause + ", SCORE"
        values_clause = values_clause + ", %s" % attrs['SCORE']
      species = attrs.get('SPECIES', '')
      if species != '':
        fields_clause = fields_clause + ", SPECIES"
        values_clause = values_clause + ", '%s'" % attrs['SPECIES']
      boot = attrs.get('BOOT', '')
      if boot != '':
        fields_clause = fields_clause + ", BOOT"
        values_clause = values_clause + ", %s" % attrs['BOOT']
      fields_clause = fields_clause + ")"
      values_clause = values_clause + ")"
      sql = "INSERT IGNORE INTO inparanoid %s VALUES %s" % (fields_clause,
                                                      values_clause)
      self.cur.execute(sql)

def main():
  if len(sys.argv) < 3:
    print "Usage: %s <species1> <species2>" % sys.argv[0]
    sys.exit(0)
  species1 = sys.argv[1].upper()
  species2 = sys.argv[2].upper()
  xml_files = glob.glob("/usr/local/blastdb/inparanoid/InParanoid.*%s-*%s.xml"
                        % (species1, species2))
  if len(xml_files) == 0:
    print "Couldn't find InParanoid ortholog file"
    sys.exit(0)
  xml_path = xml_files[0]
  xml_filename = os.path.split(xml_path)[1]
  xml_file = open(xml_path)
  parser = make_parser()
  inparanoidHandler = InParanoidHandler(xml_filename)
  parser.setContentHandler(inparanoidHandler)
  parser.parse(xml_file)
  xml_file.close()

if __name__ == '__main__':
  main()
