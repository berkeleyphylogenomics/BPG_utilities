#!/usr/bin/env python

import os, sys, MySQLdb, string, commands, re
from connect_db import *
from xml.sax import ContentHandler, make_parser, saxutils

named_element_names = [
  'interactionDetectionMethod',
  'interactionType'
]
named_element_re = re.compile('|'.join(named_element_names))
relevant_element_names = [ 
  'bibref',
  'experimentDescription',
  'interactionDetectionMethod',
  'interactor',
  'interactionType',
  'interaction',
]
relevant_element_re = re.compile('|'.join(relevant_element_names))
expected_db = {}
expected_db['interaction'] = 'dip'
expected_db['bibref'] = 'pubmed'
expected_db['interactionDetectionMethod'] = 'psi-mi'
expected_db['interactionType'] = 'psi-mi'
db_table = {}
db_table['experimentDescription'] = 'DIP_Experiment'
db_table['interactionDetectionMethod'] = 'PSI_MI_InteractionDetectionMethod'
db_table['interactor'] = 'DIP_Node'
db_table['interactionType'] = 'PSI_MI_InteractionType'
db_table['interaction'] = 'DIP_Interaction'
nonstatechanging_element_names = [
  'xref',
  'primaryRef',
  'secondaryRef',
  'names',
  'fullName',
]
nonstatechanging_element_names_re \
  = re.compile('|'.join(nonstatechanging_element_names))


class DIP_IMEx_Handler(saxutils.DefaultHandler):
  def __init__(self, filename):
    self.seen_ids = {}
    self.current_id = {}
    for element_name in relevant_element_names:
      self.current_id[element_name] = ''
      self.seen_ids[element_name] = set()
    self.current_name = {}
    for element_name in named_element_names:
      self.current_name[element_name] = ''
    self.state_stack = {}
    self.state_stack_size = 0
    self.currentContent = ''
    self.current_interaction_imex_id = ''
    self.current_uniprot_accession = ''
    self.current_pdb = ''
    self.current_experiment_ids = set()
    self.current_node_ids = set()
    self.db = MySQLdb.connect( db=phylofacts_db_name, host=phylofacts_db_host, 
                               user=phylofacts_db_user,
                               passwd=phylofacts_db_passwd, 
                               charset='utf8' )
    self.cur = self.db.cursor( MySQLdb.cursors.DictCursor )
  def current_state(self):
    return self.state_stack[self.state_stack_size - 1]
  def characters(self, chrs):
    self.currentContent = self.currentContent + chrs
  def startElement(self, element_name, attrs):
    if nonstatechanging_element_names_re.match(element_name):
      m = relevant_element_re.match(self.current_state())
      if m and len(m.group()) == len(self.current_state()):
        if self.current_state() == 'interactor':
          if element_name == 'primaryRef' or element_name == 'secondaryRef':
            db = attrs.get('db', '')
            if db == 'dip':
              self.current_id['interactor'] = attrs.get('id', '')
            elif db == 'uniprot knowledge base':
              self.current_uniprot_accession = attrs.get('id', '')
        elif self.current_state() == 'experimentDescription':
          if element_name == 'primaryRef' or element_name == 'secondaryRef':
            db = attrs.get('db', '')
            if db == 'dip':
              self.current_id['experimentDescription'] = attrs.get('id', '')
            elif db == 'wwpdb':
              self.current_pdb = attrs.get('id', '')
        else:
          if element_name == 'primaryRef':
              db = attrs.get('db', '')
              if db == expected_db[self.current_state()]:
                self.current_id[self.current_state()] = attrs.get('id', '')
              else:
                print "Unexpected db %s for %s element" % (db, 
                                                          self.current_state())
    else:
      self.state_stack[self.state_stack_size] = element_name
      self.state_stack_size = self.state_stack_size + 1
      if element_name == 'interaction':
        self.current_interaction_imex_id = attrs.get('imexId', '')
        self.current_experiment_ids = set()
        self.current_node_ids = set()
  def endElement(self, element_name):
    if nonstatechanging_element_names_re.match(element_name):
      if element_name == 'fullName':
        if named_element_re.match(self.current_state()):
          self.current_name[self.current_state()] = self.currentContent.strip()
    else:
      self.state_stack_size = self.state_stack_size - 1
      del self.state_stack[self.state_stack_size]
      m = relevant_element_re.match(element_name)
      if m and len(m.group()) == len(element_name):
        isExperiment = element_name == 'experimentDescription'
        isNode = element_name == 'interactor'
        named_element_match = named_element_re.match(element_name)
        if self.current_id[element_name] not in self.seen_ids[element_name] \
            or element_name == 'interaction':
          if named_element_match:
            sql = """REPLACE INTO %s
                              SET id = '%s',
                                  name = '%s'
                  """ % (db_table[element_name],
                          self.current_id[element_name],
                          self.current_name[element_name])
            print sql
            n = self.cur.execute(sql)
          elif isExperiment:
            if self.current_pdb == '':
              pdb_clause = ''
            else:
              pdb_clause = "pdb_id = '%s'," % self.current_pdb
            sql = """REPLACE INTO %s
                              SET id = '%s',
                                  pubmed_id = '%s',
                                  %s
                                  psi_mi_interactiondetectionmethod_id = '%s'
                  """ % (db_table[element_name],
                        self.current_id[element_name],
                        self.current_id['bibref'],
                        pdb_clause,
                        self.current_id['interactionDetectionMethod'])
            print sql
            n = self.cur.execute(sql)
          elif isNode:
            sql = """REPLACE INTO %s
                              SET id = '%s'
                  """ % (db_table[element_name],
                          self.current_id[element_name])
            print sql
            n = self.cur.execute(sql)
            if self.current_uniprot_accession != '':
              sql = """SELECT   id AS genbank_uniprot_id
                          FROM  genbank_uniprot
                          WHERE uniprot_accession = '%s'
                    """ % self.current_uniprot_accession
              n_rows = self.cur.execute(sql)
              rows = self.cur.fetchall()
              for row in rows:
                sql = """SELECT id
                           FROM genbank_uniprot_dip
                          WHERE genbank_uniprot_id = %s
                          AND   dip_node_id = '%s'
                      """ % (row['genbank_uniprot_id'],
                            self.current_id[element_name])
                n_rows = self.cur.execute(sql)
                if n_rows == 0:
                  sql = """INSERT INTO genbank_uniprot_dip
                                        ( genbank_uniprot_id,
                                          dip_node_id )
                                  VALUES ( %s, '%s' )
                        """ % (row['genbank_uniprot_id'],
                            self.current_id[element_name])
                  print sql
                  n = self.cur.execute(sql)
          elif element_name == 'interaction':
            sql = """REPLACE INTO %s
                              SET id = '%s',
                                  imex_id = '%s',
                                  psi_mi_interactiontype_id = '%s'
                  """ % (db_table[element_name],
                          self.current_id[element_name],
                          self.current_interaction_imex_id,
                          self.current_id['interactionType'])
            print sql
            n = self.cur.execute(sql)
            for experiment_id in self.current_experiment_ids:
              sql = """SELECT   id
                          FROM  DIP_Interaction_Experiment
                          WHERE dip_interaction_id = '%s'
                          AND   dip_experiment_id = '%s'
                    """ % (self.current_id[element_name], experiment_id)
              n_rows = self.cur.execute(sql)
              if n_rows == 0:
                sql = """INSERT INTO  DIP_Interaction_Experiment
                                      ( dip_interaction_id,
                                        dip_experiment_id )
                              VALUES  ( '%s', '%s' )
                      """ % (self.current_id[element_name], experiment_id)
                print sql
                n = self.cur.execute(sql)
            for node_id in self.current_node_ids:
              sql = """SELECT   id
                          FROM  DIP_Interaction_Node
                          WHERE dip_interaction_id = '%s'
                          AND   dip_node_id = '%s'
                    """ % (self.current_id[element_name], node_id)
              n_rows = self.cur.execute(sql)
              if n_rows == 0:
                sql = """INSERT INTO  DIP_Interaction_Node
                                      ( dip_interaction_id,
                                        dip_node_id )
                              VALUES  ( '%s', '%s' )
                        """ % (self.current_id[element_name], node_id)
                print sql
                n = self.cur.execute(sql)
          self.seen_ids[element_name].add(self.current_id[element_name])
        if named_element_match:
          self.current_name[element_name] = ''
        elif isExperiment:
          self.current_experiment_ids.add(self.current_id[element_name])
          self.current_pdb = ''
        elif isNode:
          self.current_node_ids.add(self.current_id[element_name])
          self.current_uniprot_accession = ''
        elif element_name == 'interaction':
          self.current_experiment_ids = set()
          self.current_node_ids = set()
          self.current_interaction_imex_id = ''
          for element_name in relevant_element_names:
            self.current_id[element_name] = ''
    self.currentContent = ''

def main():
  filename = "/usr/local/blastdb/DIP-IMEx/dip_imex_current.xml"
  f = open(filename)
  parser = make_parser()
  handler = DIP_IMEx_Handler(filename)
  parser.setContentHandler(handler)
  parser.parse(f)
  f.close()

if __name__ == '__main__':
  main()
