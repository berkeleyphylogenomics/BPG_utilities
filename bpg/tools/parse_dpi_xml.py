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
  'interactorRef',
  'interactor',
  'interactionDetectionMethod',
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
  'attribute',
]
nonstatechanging_element_names_re \
  = re.compile('|'.join(nonstatechanging_element_names))


class DIP_Handler(saxutils.DefaultHandler):
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
    self.current_attribute = ''
    self.current_evidence_for_ac = ''
    self.current_evidence_for_name = ''
    self.current_evidence_scale_ac = ''
    self.current_evidence_scale_name = ''
    self.current_quality_status = ''
    self.current_experiment_ids = set()
    self.current_node_ids = set()
    self.interactors = {}
    self.current_interactor_id = ''
    self.db = MySQLdb.connect( db=phylofacts_db_name, host=phylofacts_db_host, 
                               user=phylofacts_db_user,
                               passwd=phylofacts_db_passwd, 
                               charset='utf8' )
    self.cur = self.db.cursor( MySQLdb.cursors.DictCursor )
    self.num_interactions = 0
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
          elif element_name == 'attribute':
            self.current_attribute = attrs.get('name', '')
            self.currentContent = ''
        elif self.current_state() == 'interaction' and \
            element_name == 'attribute':
          self.current_attribute = attrs.get('name', '')
          self.currentContent = ''
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
      elif element_name == 'interactor':
        self.current_interactor_id = attrs.get('id', '')
  def endElement(self, element_name):
    if nonstatechanging_element_names_re.match(element_name):
      if element_name == 'fullName':
        if named_element_re.match(self.current_state()):
          self.current_name[self.current_state()] = self.currentContent.strip()
      elif element_name == 'attribute':
        if self.current_attribute == 'quality-status':
          self.current_quality_status = self.currentContent.strip()
        elif self.current_attribute == 'dip:evidence-for-ac':
          self.current_evidence_for_ac = self.currentContent.strip()
        elif self.current_attribute == 'dip:evidence-for-name':
          self.current_evidence_for_name = self.currentContent.strip()
        elif self.current_attribute == 'dip:evidence-scale-ac':
          self.current_evidence_scale_ac = self.currentContent.strip()
        elif self.current_attribute == 'dip:evidence-scale-name':
          self.current_evidence_scale_name = self.currentContent.strip()
        elif self.current_attribute == 'dip:evidence-for-name':
          self.current_evidence_for_name = self.currentContent.strip()
        self.current_attribute = ''
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
            if self.current_evidence_for_ac != '':
              sql = """SELECT accession
                        FROM  DIP_Evidence_for
                        WHERE accession = '%s'
                    """ % self.current_evidence_for_ac
              n_rows = self.cur.execute(sql)
              if n_rows == 0:
                sql = """INSERT INTO  DIP_Evidence_for
                                      ( accession, name )
                              VALUES  ('%s', '%s')
                      """ % (self.current_evidence_for_ac,
                            self.current_evidence_for_name)
                print sql
                n = self.cur.execute(sql)
              sql = """SELECT id
                        FROM  DIP_Experiment_Evidence_for
                        WHERE dip_experiment_id = '%s'
                        AND   dip_evidence_for_accession = '%s'
                    """ % (self.current_id[element_name],
                          self.current_evidence_for_ac)
              n_rows = self.cur.execute(sql)
              if n_rows == 0:
                sql = """INSERT INTO  DIP_Experiment_Evidence_for
                                      ( dip_experiment_id,
                                        dip_evidence_for_accession
                                      ) VALUES (
                                        '%s', '%s'
                                      )
                      """ % (self.current_id[element_name],
                          self.current_evidence_for_ac)
                print sql
                n = self.cur.execute(sql)
            if self.current_evidence_scale_ac != '':
              sql = """SELECT accession
                        FROM  DIP_Evidence_Scale
                        WHERE accession = '%s'
                    """ % self.current_evidence_scale_ac
              n_rows = self.cur.execute(sql)
              if n_rows == 0:
                sql = """INSERT INTO  DIP_Evidence_Scale
                                      ( accession, name )
                              VALUES  ('%s', '%s')
                      """ % (self.current_evidence_scale_ac,
                            self.current_evidence_scale_name)
                print sql
                n = self.cur.execute(sql)
              sql = """SELECT id
                        FROM  DIP_Experiment_Evidence_Scale
                        WHERE dip_experiment_id = '%s'
                        AND   dip_evidence_scale_accession = '%s'
                    """ % (self.current_id[element_name],
                          self.current_evidence_scale_ac)
              n_rows = self.cur.execute(sql)
              if n_rows == 0:
                sql = """INSERT INTO  DIP_Experiment_Evidence_scale
                                      ( dip_experiment_id,
                                        dip_evidence_scale_accession
                                      ) VALUES (
                                        '%s', '%s'
                                      )
                      """ % (self.current_id[element_name],
                          self.current_evidence_scale_ac)
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
            if self.current_quality_status == '':
              quality_status_clause = ''
            else:
              quality_status_clause = "quality_status = '%s'," \
                                      % self.current_quality_status
            sql = """REPLACE INTO %s
                              SET id = '%s',
                                  imex_id = '%s',
                                  %s
                                  psi_mi_interactiontype_id = '%s'
                  """ % (db_table[element_name],
                          self.current_id[element_name],
                          self.current_interaction_imex_id,
                          quality_status_clause,
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
          self.current_evidence_for_ac = ''
          self.current_evidence_for_name= ''
          self.current_evidence_scale_ac = ''
          self.current_evidence_scale_name= ''
        elif isNode:
          self.current_node_ids.add(self.current_id[element_name])
          self.current_uniprot_accession = ''
          if self.current_interactor_id != '':
            self.interactors[self.current_interactor_id] \
                = self.current_id[element_name]
          self.current_interactor_id = ''
        elif element_name == 'interactorRef':
          interactorRef = self.currentContent.strip()
          self.current_node_ids.add(self.interactors[interactorRef])
        elif element_name == 'interaction':
          self.current_experiment_ids = set()
          self.current_node_ids = set()
          self.current_interaction_imex_id = ''
          self.current_quality_status = ''
          for element_name in relevant_element_names:
            self.current_id[element_name] = ''
          self.num_interactions = self.num_interactions + 1
    self.currentContent = ''

def main():
  filename = "/usr/local/blastdb/DIP-IMEx/dip_current.xml"
  f = open(filename)
  parser = make_parser()
  handler = DIP_Handler(filename)
  parser.setContentHandler(handler)
  parser.parse(f)
  f.close()

if __name__ == '__main__':
  main()
