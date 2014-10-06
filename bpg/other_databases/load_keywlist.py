#!/usr/bin/env python

from array import array
from optparse import OptionParser

import django.db.models.base
from django.utils.encoding import smart_unicode

from pfacts003.phylofacts.models import GO_Term, Keyword,\
    Keyword_GO_Term, KeywordKeyword, KeywordSynonym, KeywordWWW

class KeywordEntry:
  """
  Represent the flat-text file entry for a keyword  as an object.

  An example of a flat-file entry that is read and parsed follows:

  ID   2Fe-2S.
  AC   KW-0001
  DE   Protein which contains at least one 2Fe-2S iron-sulfur cluster...
  DE   atoms complexed to 2 inorganic sulfides and 4 sulfur atoms of
  DE   cysteines from the protein.
  SY   Fe2S2; [2Fe-2S] cluster; [Fe2S2] cluster; Fe2/S2 (inorganic) c...
  SY   Di-mu-sulfido-diiron; 2 iron, 2 sulfur cluster binding.
  GO   GO:0051537; 2 iron, 2 sulfur cluster binding
  HI   Ligand: Iron; Iron-sulfur; 2Fe-2S.
  HI   Ligand: Metal-binding; 2Fe-2S.
  CA   Ligand.
  //

  The two-character line codes on the left are described at the top of
  the keyword file to be read and is included here for your convenience.
  The flat-file records have
  ---------  ---------------------------     ----------------------
  Line code  Content                         Occurrence in an entry
  ---------  ---------------------------     ----------------------
  ID         Identifier (keyword)            Once; starts a kw entry
  IC         Identifier (category)           Once; starts a cat entry
  AC         Accession (KW-xxxx)             Once
  DE         Definition                      Once or more
  SY         Synonyms                        Optional; once or more
  GO         Gene ontology (GO) mapping      Optional; once or more
  HI         Hierarchy                       Optional; once or more
  WW         Relevant WWW site               Optional; once or more
  CA         Category                        Once per keyword entry;
                                             absent in category entries
  //         Terminator                      Once; ends an entry

  """
  def __init__(self, id=''):
    self.id = id
    self.accession = ''
    self.definition = array('c')
    self.synonyms = set()
    self.www = set()
    self.go = set()
    self.category_id = ''
    self.more_general_keyword_ids = set()
    self.more_general_keywords = set()
    self.more_specific_keywords = set()
    self.category = None
    self.db_object = None

  def append_definition(self, definition):
    if len(self.definition) > 0:
      self.definition.fromstring(' ')
    self.definition.fromstring(definition.strip())

  def get_definition(self):
    return self.definition.tostring()

  def update_general_and_specific_keywords(self, keyword_of_id):
    if len(self.category_id) > 0:
      self.category = keyword_of_id[self.category_id]
    self.more_general_keyords = set()
    for id in self.more_general_keyword_ids:
      more_general_keyword = keyword_of_id[id]
      self.more_general_keywords.add(more_general_keyword)
      more_general_keyword.more_specific_keywords.add(self)

def print_parsed_keywords(categories, keyword_of_id):
  id_list = keyword_of_id.keys()
  id_list.sort()
  print "Categories:"
  for category in categories:
    print "  %s" % category.id
  for id in id_list:
    keyword = keyword_of_id[id]
    if keyword.category:
      print "ID   %s." % id
    else:
      print "IC   %s." % id
    print "AC   %s" % keyword.accession
    print "DE   %s" % keyword.get_definition()
    if len(keyword.synonyms) > 0:
      print "SY   %s" % '; '.join(keyword.synonyms)
    if len(keyword.www) > 0:
      for www in keyword.www:
        print "WW   %s" % www
    if len(keyword.go) > 0:
      for go in keyword.go:
        print "GO   %s; %s" % (go.acc, go.name)
    for more_general_keyword in keyword.more_general_keywords:
      print "HI   %s." % more_general_keyword.id
    if keyword.category:
      print "CA   %s" % keyword.category.id
    print "//"

def main(options, args):

  f = open('/clusterfs/ohana/external/UniProt/to_import/keywlist.txt')
  lines = [line.rstrip() for line in f.readlines()]
  f.close()
  keyword_of_id = {}
  categories = set()
  current_keyword = None

  # Parse entry
  for line in lines:
    if line[0:2] == 'ID':
      id = line[5:].rstrip('.')
      current_keyword = KeywordEntry(id)
      keyword_of_id[id] = current_keyword
    elif line[0:2] == 'IC':
      id = line[5:].rstrip('.')
      current_keyword = KeywordEntry(id)
      keyword_of_id[id] = current_keyword
      categories.add(current_keyword)
    elif line[0:2] == 'AC':
      current_keyword.accession = line[5:].rstrip()
    elif line[0:2] == 'DE':
      current_keyword.append_definition(line[5:].rstrip())
    elif line[0:2] == 'SY':
      synonym_line = line[5:].rstrip()[:-1]
      fields = synonym_line.split(';')
      for field in fields:
        current_keyword.synonyms.add(field.strip())
    elif line[0:2] == 'WW':
      current_keyword.www.add(line[5:].rstrip())
    elif line[0:2] == 'GO':
      go_accession = line[5:].rstrip().split(';')[0]
      try:
        go_term = GO_Term.objects.get(acc__exact = go_accession)
        current_keyword.go.add(go_term)
      except GO_Term.DoesNotExist:
        print "Unrecognized GO accession %s" % go_accession
    elif line[0:2] == 'CA':
      current_keyword.category_id = line[5:].rstrip('.')
    elif line[0:2] == 'HI':
      fields = line[5:].replace(':',';').split(';')
      more_general_keyword_id = fields[-2].strip().rstrip('.')
      current_keyword.more_general_keyword_ids.add(more_general_keyword_id)
    elif line[0:2] == '//':
      current_keyword = None
  for id in keyword_of_id:
    keyword_of_id[id].update_general_and_specific_keywords(keyword_of_id)

  # Find the Django Keyword Object or create it if it doesn't exist
  for category in categories:
    category_objects = Keyword.objects.filter(accession__exact =
                                                category.accession)
    if category_objects:
      category.db_object = category_objects[0]
    else:
      category.db_object = Keyword.objects.create(id = category.id,
                                accession = category.accession,
                                definition = category.get_definition())

  for id in keyword_of_id:
    keyword = keyword_of_id[id]
    keyword_objects = Keyword.objects.filter(accession__exact =
                                                keyword.accession)
    if keyword_objects:
      keyword.db_object = keyword_objects[0]
      keyword.db_object.id = keyword.id
      keyword.db_object.definition = keyword.get_definition()
      keyword.db_object.save()
    else:
      keyword.db_object = Keyword.objects.create(identifier = keyword.id,
                                                  accession = keyword.accession,
                                          definition = keyword.get_definition(),
                                          category = keyword.category.db_object)
    keyword_go_term_objects = Keyword_GO_Term.objects.filter(
                                              keyword__exact = keyword.db_object)
    go_accessions = [go_term.acc for go_term in keyword.go]
    for object in keyword_go_term_objects:
      if object.go_term.acc not in go_accessions:
        print "Deleting obsolete GO term %s for keyword %s" \
            % (object.go_term.acc, keyword.accession)
        object.delete()
    db_go_accessions = [keyword_go_term.go_term.acc for keyword_go_term 
                          in keyword_go_term_objects]
    for go_term in keyword.go:
      if go_term.acc not in db_go_accessions:
        Keyword_GO_Term.objects.create(go_term = go_term, 
                                      keyword = keyword.db_object)
    keyword_synonym_objects = KeywordSynonym.objects.filter(
                                              keyword__exact = keyword.db_object)
    for object in keyword_synonym_objects:
      if object.synonym not in keyword.synonyms:
        print "Deleting obsolete synonym %s for keyword %s" \
            % (object.synonym, keyword.accession)
        object.delete()
    db_synonyms = set([object.synonym for object in keyword_synonym_objects])
    for synonym in keyword.synonyms:
      if synonym not in db_synonyms:
        KeywordSynonym.objects.create(synonym = synonym, 
                                      keyword = keyword.db_object)
    keyword_www_objects = KeywordWWW.objects.filter(
                                              keyword__exact = keyword.db_object)
    for object in keyword_www_objects:
      if object.www_site not in keyword.www:
        print "Deleting obsolete www %s for keyword %s" \
            % (object.www_site, keyword.accession)
        object.delete()
    db_wwws = set([object.www_site for object in keyword_www_objects])
    for www in keyword.www:
      if www not in db_wwws:
        KeywordWWW.objects.create(www_site = www, 
                                  keyword = keyword.db_object)
  # Now that all keyword objects exist in the database, update keyword_keyword
  for id in keyword_of_id:
    keyword = keyword_of_id[id]
    keyword_keyword_objects = KeywordKeyword.objects.filter(
                                    specific_keyword__exact = keyword.db_object)
    more_general_accessions = [general_keyword.accession for general_keyword
                                in keyword.more_general_keywords]
    for object in keyword_keyword_objects:
      if object.general_keyword.accession not in more_general_accessions:
        print "Deleting obsolete more general keyword %s for keyword %s" \
            % (object.general_keyword.accession, keyword.accession)
        object.delete()
    more_general_db_accessions = [object.general_keyword.accession 
                                          for object in keyword_keyword_objects]
    for more_general_keyword in keyword.more_general_keywords:
      if more_general_keyword.accession not in more_general_db_accessions:
        KeywordKeyword.objects.create(specific_keyword = keyword.db_object,
                                general_keyword = more_general_keyword.db_object)

if __name__ == '__main__':
  parser = OptionParser(usage='%prog')
  (options, args) = parser.parse_args()
  main(options, args)
