from datetime import datetime, timedelta
import random
import time
import logging
import os
import string
import StringIO
import subprocess
#import functools
import tempfile
from textwrap import fill
import time
import hashlib, base64
from lxml import etree, objectify
from itertools import groupby
import glob
from Bio import SeqIO, AlignIO
from django.contrib.auth.models import User
from django.db import models
from django.db.models import Max, Q, F, Sum
from django.db import connection
from django.db.models.query import QuerySet
from django.utils.safestring import mark_safe
from model_utils.managers import PassThroughManager
try:
    import json
except:
    import simplejson as json
from pfacts003.utils.links import make_link, family_phyloscope_url, \
                  family_phyloscope_link, ncbi_taxonomy_url, \
                  ncbi_taxonomy_link, ncbi_entrez_url, ncbi_entrez_link, \
                  uniprot_url, uniprot_link, \
                  pfam_url, pfam_link, kegg_icon, uniprot_lit_icon, \
                  swissprot_icon, uniprot_evidence_icon, sequence_orthologs_link, \
                  uniprot_ppi_link, uniprot_orthologs_url, uniprot_orthologs_link, \
                  make_kegg_map_links, \
                  uniprot_ec_link, dip_url, netscope_url, pmid_link, make_pfam_links, make_link
from pfacts003.utils.consts import max_phyloscope_tree_size
from pfacts003.utils.id_patterns import *
from pfacts003.utils.hmm import hmmbuild, hmmpress
import pfacts003.utils.hmm as hmm_utils
from pfacts003.utils.alignment import mafft
from pfacts003.utils.annotation import TreeAnnotation, GOAnnotation, Annotation, summarize, consensus
#from pfacts003.phylofacts.sqlite_naming import name as sqlite_name
from bpg.common.fatcat import fatcat_ortholog_consensus_uniprot_description, fatcat_ortholog_gene_clustering
from bpg.common.utils import get_phobius_tmhelix
from bpg.common.BPGPWID import pairwise_identity_belvu as pwid_belvu
from django.utils.datastructures import SortedDict
from bpg.common.utils.sequenceutils.fasta_utils import get_uniprot_from_header

experimental_evidence_codes = set(['EXP','IDA','IPI','IMP','IGI','IEP'])

# TODO: Hardcoding a path like this is a bad idea. This needs to go
# in some sort of config file.
# Also, according to RSD, copies of this def are floating around in the
# codebase. That should be fixed
# Note made by ST on 1/31/11
def get_dir_of_family_accession(bpg_accession):
    return '/clusterfs/ohana/bpg/pfacts/%s/%s/%s' % (bpg_accession[0:4],
                                                    bpg_accession[0:7],
                                                    bpg_accession)

def make_ranked_summary(go_annotation_dict, name_of_acc):
    summary = {}
    for term_type in go_annotation_dict.keys():
        if len(go_annotation_dict[term_type]) > 0:
            terms = [(go_annotation_dict[term_type][acc].priority,
                      name_of_acc[acc], go_annotation_dict[term_type][acc].evidence,
                      acc, go_annotation_dict[term_type][acc].name)
                      for acc in go_annotation_dict[term_type].keys()]
            terms.sort()
            summary[term_type] = [(acc, name, evidence, description)
                                  for (priority, name, evidence, acc, description)
                                  in terms]
    return summary

class BigIntegerField(models.IntegerField):
    """Custom Big Integer Type.

    Surprisingly, there is no BigInteger field currently (v1.1) in Django.
    Therefore we create our own.
    """
    empty_strings_allowed=False

    def get_internal_type(self):
        return "BigIntegerField"

    def db_type(self, connection):
        return 'bigint'


class memoized(object):
   '''Decorator. Caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned 
   (not reevaluated).
   '''
   def __init__(self, func):
      self.func = func
      self.cache = {}
   def __call__(self, *args):
      try:
         return self.cache[args]
      except KeyError:
         value = self.func(*args)
         self.cache[args] = value
         return value
      except TypeError:
         # uncachable -- for instance, passing a list as an argument.
         # Better to not cache than to blow up entirely.
         return self.func(*args)
   def __repr__(self):
      '''Return the function's docstring.'''
      return self.func.__doc__
   def __get__(self, obj, objtype):
      '''Support instance methods.'''
      return functools.partial(self.__call__, obj)


class OrthologTypes(object):
    """Ortholog Types encapsulated into an object"""

    SuperOrtholog, PHOG_T_Tight, PHOG_T_Medium, PHOG_T_Loose, \
        PHOG_T_Custom = xrange(5)

trivial_translation = string.maketrans('', '')
dotlowercase = '.' + string.lowercase

# These defs should ideally be placed in pfactoo3/utils/links.py.
# But, that will lead to circular imports since OrthologTypes is
# required by these. BUT, OrthologTypes will be going away soon
# according to RSD. We'll wait until that happens, and then
# refactor as necessary.
# Note made by ST on 1/31/11
def phog_url(phog, ortholog_type = OrthologTypes.SuperOrtholog,
                  threshold = 0.0):
    return phog.get_absolute_url(ortholog_type, threshold)

def phog_link(phog, ortholog_type = OrthologTypes.SuperOrtholog,
                    threshold = 0.0):
    return make_link(phog.get_accession(ortholog_type, threshold),
                    phog_url(phog, ortholog_type, threshold), 'phog')

class AlignedSequence(models.Model):
    id = models.AutoField(primary_key=True)
    chars = models.TextField()
    seguid = models.CharField(max_length=27, unique=True)
    class Meta:
        db_table = u'aligned_sequence'


class Sequence(models.Model):
    id = models.AutoField(primary_key=True)
    chars = models.TextField()
    seguid = models.CharField(max_length=27, unique=True)
    def wrap(self, width=80):
        return fill(self.chars, 80)
    class Meta:
        db_table = u'sequence'


preset_thresholds = {
    OrthologTypes.SuperOrtholog: 0.0,
    OrthologTypes.PHOG_T_Tight: 0.09375,
    OrthologTypes.PHOG_T_Medium: 0.296874,
    OrthologTypes.PHOG_T_Loose: 0.9375,
}
def get_ortholog_type_of_threshold(threshold):
    return dict((y,x) for x,y in \
        preset_thresholds.items()).get(threshold, OrthologTypes.PHOG_T_Custom)

def get_ortholog_type_threshold_from_phog_accession(phog_accession):
    ortholog_type = OrthologTypes.SuperOrtholog
    threshold = 0.0
    if len(phog_accession) >= 14 and phog_accession[13] == 'T':
        if phog_accession[14] == 'C':
            ortholog_type = OrthologTypes.PHOG_T_Tight
            threshold = preset_thresholds[ortholog_type]
        elif phog_accession[14] == 'M':
            ortholog_type = OrthologTypes.PHOG_T_Medium
            threshold = preset_thresholds[ortholog_type]
        elif phog_accession[14] == 'D':
            ortholog_type = OrthologTypes.PHOG_T_Loose
            threshold = preset_thresholds[ortholog_type]
        else:
            ortholog_type = OrthologTypes.PHOG_T_Custom
            threshold =  float(phog_accession[14:])
    return (ortholog_type, threshold)


class UniProtTaxonomy(models.Model):
    id = models.IntegerField(primary_key=True)
    mnemonic = models.TextField(null=True)
    scientific_name = models.TextField(null=True)
    common_name = models.TextField(null=True)
    synonym = models.TextField(null=True)
    other_names = models.TextField(null=True)
    reviewed = models.BooleanField(default=False)
    rank = models.TextField(null=True)
    parent = models.ForeignKey('self', related_name='children', null=True)
    left_id = models.IntegerField(null=True)
    right_id = models.IntegerField(null=True)

    def lineage(self):
        lineage = [self]
        s = lineage[-1].parent
        while s != None and s.scientific_name != 'root':
            lineage.append(s)
            s = s.parent
        return reversed(lineage)

    def is_leaf(self):
        return self.left_id + 1 == self.right_id

    def get_children(self):
        return UniProtTaxonomy.objects.filter(
            left_id__gte=self.left_id,
            right_id__lte=self.right_id
        )

    def __unicode__(self):
        if self.common_name == None:
            if self.scientific_name == None:
                return str(self.id)
            else:
                return self.scientific_name
        else:
            if self.scientific_name == None:
                return self.common_name
            else:
                return "%s (%s)" % (self.scientific_name, self.common_name)

    def get_absolute_url(self):
        return ncbi_taxonomy_url(self.id)

    def get_link(self):
        if self.id < 0:
            return self.scientific_name
        else:
            return ncbi_taxonomy_link(self, self.id)

    class Meta:
        db_table = u'uniprot_taxonomy'

comma_range_re = re.compile('(\d+),(\d+)$')
class SequenceHeader(models.Model):
    id = models.AutoField(primary_key=True)
    header = models.TextField(blank=True)
    sequence = models.ForeignKey(Sequence)
    uniprot = models.ForeignKey('UniProt', null=True)
    taxon = models.ForeignKey(UniProtTaxonomy, null=True)

    def get_uniprot(self):
        """ This function was made to deal with our alternate accession problem. 
        Here we see if there is a sequence header and a uniprot.  If so, we check to make sure
        the parsed header accession matches the uniprot accession.  If so, return it.  Otherwise,
        Try to find the UniProt entry with this accession.  
        """
        if (self.header):
            h = get_uniprot_from_header(self.header)
            if (h):
                return h
            else:
                return None
                #if (self.uniprot):
                #    return (self.uniprot)
        return None
            
    def as_fasta(self):
        return '>%s\n%s\n' % (self.header, self.sequence.wrap())

    def as_seqrecord(self):
        h = StringIO.StringIO(self.as_fasta())
        rec = SeqIO.read(h, 'fasta')
        h.close()
        return rec

    def parse_header(self):
        # This should be merged with other functions below, for now lets do this....
        try:
            return self._parse_header
        except:
            pass

        header_information = {'db': None, 'accession': None, 'identifier': None,
                        'species': None, 'taxon_id':None, 'gene_name': None,
                        'scientific_name': None, 'common_name':None,  'sequence_version': None,
                        'protein_existence': None, 'description': ''}

        words = self.header.split()
        if not (words[0].startswith('tr') or words[0].startswith('sp') or words[0].startswith('lcl')):
            header_information['description'] = self.header
        else:
            if '|' in words[0]:
                separator = '|'
            elif ':' in words[0]:
                separator = ':'
            else:
                header_information['description'] = self.header
                self._parse_header = header_information
                return self._parse_header
            try:
                header_information['db'] = words[0].split(separator)[0]
                header_information['accession'] = words[0].split(separator)[1]
                header_information['identifier'] = words[0].split(separator)[2]
            except:
                pass
            words.pop(0)
            key = 'description'
            while words:
                word = words.pop(0)
                if re.compile("OS=").match(word):
                    key = 'species'
                    header_information[key] = word.strip("OS=")
                elif re.compile("GN=").match(word):
                    key = 'gene_name'
                    header_information[key] = word.strip("GN=")
                elif re.compile("SV=").match(word):
                    key = 'sequence_version'
                    header_information[key] = word.strip("SV=")
                elif re.compile("PE=").match(word):
                    key = 'protein_existence'
                    header_information[key] = word.strip("PE=")
                else:
                    header_information[key] += ' ' + word
        # ridiculous :)
        header_information['description'].lstrip()
        # look for the species in our uniprot taxonomy.
        try:
            t = UniProtTaxonomy.objects.get(scientific_name = header_information['species'])
            header_information['taxon_id'] = t.id
            header_information['scientific_name'] = t.scientific_name
            if (t.common_name):
                header_information['common_name'] = t.common_name
        except:
            pass
        self._parse_header = header_information
        return self._parse_header

    def description(self):
        return getattr(self.uniprot, 'description', self.header)

    def identifier(self):
        ret = getattr(self.uniprot, 'accession', False)
        if ret:
            return ret

        ret = self.description().split()[0]

        parts = ret.split('|')
        for part in parts:
            if len(part) > 3:
                ret = part
                break
        return ret

    def __unicode__(self):
        return unicode(self.identifier())

    def get_orthologs_link(self):
        return sequence_orthologs_link(self.identifier())

    def get_identifier_url(self):
        identifier = self.identifier()
        if self.uniprot is not None and \
            (uniprot_accession_re1.match(identifier) is not None or \
              uniprot_accession_re2.match(identifier) is not None):
            return self.uniprot.get_absolute_url()
        else:
            return ncbi_entrez_url(identifier)

    @models.permalink
    def get_absolute_url(self):
        return('sequences', [str(self.id)])

    def get_link(self):
        identifier = self.identifier()
        if self.uniprot is not None and \
            (uniprot_accession_re1.match(identifier) is not None or \
              uniprot_accession_re2.match(identifier) is not None):
            return self.uniprot.get_link()
        else:
            return ncbi_entrez_link(identifier)

    def get_escaped_header_for_tree(self):
        header = self.header
        words = header.split()
        r = comma_range_re.search(header)
        if r:
            header = header[0:len(header)-len(r.group(0))].rstrip()
            range_word = '%s/%s-%s' % (words[0], r.group(1), r.group(2))
            header = '%s' % (range_word)
        else:
            header = words[0]
        escaped_header = header
        escaped_header = escaped_header.replace('(', '%28')
        escaped_header = escaped_header.replace(')', '%29')
        escaped_header = escaped_header.replace(':', '%3A')
        escaped_header = escaped_header.replace(',', '%2C')
        escaped_header = escaped_header.replace(';', '%3B')
        return(escaped_header)

    class Meta:
        db_table = u'sequence_header'

class BuildAlignmentNotes(models.Model):
    id = models.AutoField(primary_key=True)
    notes = models.TextField()

    class Meta:
        db_table = u'build_alignment_notes'

    def __unicode__(self):
        return unicode(self.notes)


class FamilyType(models.Model):
    id = models.CharField(max_length=1,
        primary_key=True)
    description = models.CharField(max_length=18)

    class Meta:
        db_table = u'family_type'

    def __unicode__(self):
        return unicode(self.description)
    __str__ = __unicode__

class PHOGRow:
    def __init__(self, phog, ortholog_type = OrthologTypes.SuperOrtholog,
                  threshold = 0.0):
        self.accession = phog.get_accession(ortholog_type, threshold)
        self.accession_link = phog_link(phog)
        self.description = phog.get_description(ortholog_type=ortholog_type,
                                                threshold=threshold)
        self.description_link = make_link(self.description,
                       phog_url(phog, ortholog_type, threshold), 'phog')

        # this takes about 5'' for bpg0135873
        taxon = phog.get_taxon(ortholog_type, threshold)
        self.taxon_id = taxon.id
        self.taxon_name = taxon.__str__()
        self.taxon_link = taxon.get_link()
        self.num_nonredundant_sequences = phog.get_num_nonredundant_sequences(
                                            ortholog_type, threshold)
        self.family = phog.tree.family.__str__()
        self.family_url = phog.tree.family.get_absolute_url()
        self.family_type = phog.tree.family.get_family_type_str()
        self.alignment_type = phog.tree.family.get_alignment_type_str()
        self.alignment_length = 0
        self.tree_node_aligned_sequences = TreeNodeAlignment.objects.filter(
              tree_node__exact = phog.tree.treenodes.get(left_id=1))
        if self.tree_node_aligned_sequences:
            aligned_seq = self.tree_node_aligned_sequences[0].aligned_sequence.chars
            self.alignment_length = len(str(aligned_seq).translate(
                                                trivial_translation, dotlowercase))
        self.tree_url = phog.get_tree_url(ortholog_type, threshold)
        self.pfam_links = make_pfam_links(phog.tree.family)

def make_go_summary_line(go_summary):
    return ', '.join(['<a href="http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=%s">%s</a> (<a href="http://www.geneontology.org/GO.evidence.shtml#%s" title="%s">%s</a>)' \
      % (acc, name, evidence.lower(), description, evidence)
      for (acc, name, evidence, description) in go_summary])


class FamilyCache(models.Model):
    id = models.AutoField(primary_key=True)
    family_id = models.IntegerField(unique=True)
    cached_data = models.TextField()

    class Meta:
        db_table = u'family_cache'

########## metacyc classes
class MetacycVersion(models.Model):
    version = models.TextField(primary_key=True)

    def __unicode__(self):
        return unicode(self.version)
        
    class Meta:
        db_table = u'metacyc\".\"metacyc_version'

class MetacycProtein(models.Model):
    id = models.ForeignKey('MetacycNode', primary_key=True, db_column='id')
    common_name = models.TextField(null=True)
    db = models.TextField(null=True)
    # Removed this field since uniprot accessions are now stored in arrays.
    #uniprot_accession = models.TextField(null=True)
    type = models.TextField()
    
    def __unicode__(self):
        return unicode(self.id)
        
    class Meta:
        db_table = u'metacyc\".\"metacyc_protein'

class MetacycReaction(models.Model):
    rxn_id = models.ForeignKey('MetacycNode', db_column='id')
    common_name = models.TextField(null=True)
    systematic_name = models.TextField(null=True)
    ec_number = models.TextField(null=True)
    db = models.TextField(null=True)    
    type = models.TextField()
    
    def __unicode__(self):
        return unicode(self.id)
        
    class Meta:
        db_table = u'metacyc\".\"metacyc_reaction'

class MetacycPathway(models.Model):
    pwy_id = models.ForeignKey('MetacycNode', db_column='id')
    common_name = models.TextField(null=True)
    db = models.TextField(null=True)    
    type = models.TextField()
    
    def __unicode__(self):
        return unicode(self.id)
        
    class Meta:
        db_table = u'metacyc\".\"metacyc_pathway'

class MetacycNode(models.Model):
    id = models.TextField(primary_key=True)
    common_name = models.TextField(null=True)
    db = models.TextField(null=True)
    type = models.TextField()

    @classmethod
    def get_membership_dict(cls, uniprot_acc_or_identifier):
        '''This function uses a raw sql query to get a membership dictionary
        that will be directly passed by the views to the template. The reason for
        bypassing the models structure is due to recursive calls that need to be
        made on the metacyc.metacyc_membership table that in this version of django
        is not supported. Please update me if a new version of django is used and the
        sql queries can be funneled through django.

        Alternatively, and this might be a better option, create views in the database
        by denormalizing the metacyc data and then create a model that calls this postgres view.
        '''
        output_dict = {}
        cursor = connection.cursor()
        #First get the uniprot accession
        uniprot_obj = UniProt.find_by_accession_or_identifier(uniprot_acc_or_identifier)
        uniprot_accession = uniprot_obj.accession
        # Get the gene product, version info and db info
        uniprot_accession_sql = """select id, db from metacyc.metacyc_protein
        where %s = ANY(uniprot_accession)"""
        cursor.execute(uniprot_accession_sql, (uniprot_accession,))
        metacyc_protein_product = cursor.fetchall()
        if not metacyc_protein_product:
            return output_dict
        metacyc_protein_product, metacyc_protein_db = cls.select_specific_db_for_protein(
            metacyc_protein_product)
        metacyc_version = MetacycVersion.objects.all()[0]
        output_dict['version'] = metacyc_version.version
        output_dict['product'] = metacyc_protein_product
        output_dict['db'] = metacyc_protein_db
        # Get recursive relations using an sql query
        all_relationships = cls.get_membership_recursively([output_dict['product']])
        # Presently we will only provide a list without the relations. Also
        # limit to results matching the same database.
        pathway_dict = cls.convert_relationships_to_dict(all_relationships, metacyc_protein_db)
        output_dict.update(pathway_dict)
        return output_dict

    @classmethod
    def get_membership_recursively(cls, list_of_metacyc_products):
        '''Uses raw sql to get a membership list of member tuples.'''
        cursor = connection.cursor()
        m_list = [str(line) for line in list_of_metacyc_products]
        m_list = ["'%s'" % line for line in m_list]
        m_list = "(%s)" % ','.join(m_list)
        sql_query = """
        WITH RECURSIVE membership( member_id, group_id)
        AS (
        SELECT
        metacyc.metacyc_membership.metacyc_member_id,
        metacyc.metacyc_membership.metacyc_group_id
        FROM
        metacyc.metacyc_membership
        WHERE
        metacyc_membership.metacyc_member_id in %s
        UNION ALL
        SELECT
        mm.metacyc_member_id,
        mm.metacyc_group_id
        FROM
        metacyc.metacyc_membership mm, membership m
        WHERE
        m.group_id = mm.metacyc_member_id
        )
        SELECT member_id, group_id
        FROM membership
        """
        cursor.execute(sql_query % m_list)
        all_relationships = cursor.fetchall()
        cursor.close()
        return all_relationships
    
    @classmethod
    def convert_to_metacyc_product(cls, uniprot_accession_list):
        '''Returns gene products that match the uniprot accession list and
        a dictionary with metacyc and uniprot info.
        If not gene products found, return None. Again raw sql
        simply because django doesnt allow the models to do this.'''
        if not uniprot_accession_list:
            return None, None
        # Generate a dictionary with uniprot accessions as keys and uniprot objects as the items
        uniprot_accession_dict = dict([(line.accession, line) for line in uniprot_accession_list])
        uniprot_accession_list_temp = [str(line) for line in uniprot_accession_dict]
        uniprot_accession_list_temp = ["'%s'" % line for line in uniprot_accession_list_temp]
        uniprot_accession_list_sql = "[%s]" % ','.join(uniprot_accession_list_temp)
        uniprot_accession_sql = """select id, uniprot_accession from metacyc.metacyc_protein
        where ARRAY%s && uniprot_accession"""
        cursor = connection.cursor()
        cursor.execute(uniprot_accession_sql % uniprot_accession_list_sql)
        metacyc_ids = cursor.fetchall()
        if metacyc_ids:
            # Get mapping of metacyc ids to uniprot accessions from original list.
            uniprot_acc_set = set(uniprot_accession_dict.keys())
            initial_dict = [(key, list(set(val) & uniprot_acc_set))
                            for key, val in metacyc_ids]
            uniprot_dict = {}
            # convert the list of uniprot accessions to uniprot objects from uniprot_accession_dict
            for metacyc_id, seq_list in initial_dict:
                uniprot_dict[metacyc_id] = [value for key, value in uniprot_accession_dict.items()
                                            if key in seq_list]
            return [m_id[0] for m_id in metacyc_ids], uniprot_dict
        else:
            return None, None
                                                                                                        
    @classmethod
    def map_metacyc_ids_to_reactions(cls, met_proteins, metacyc_uniprot_dict, relations):
        '''Takes metacyc_ids, metacyc to uniprot dictionary, relation tuples
        and returns a dictionary with reaction
        ids and pathway ids as keys and a list of uniprot ids as values.

        Within the for loop, iterate through the tuple list, add to the dictionary - values
        of uniprot accessions and keys of reactions or pathways (uses the fact that the tuples
        are listed first by metacyc proteins, then with reactions within reactions, then with
        reactions in pathways and finally pathways within pathways. So if some point, the
        recursive query within MetacycNode.get_membership_recursively changes, this will need
        to be modified).
        '''
        outdict = {}
        for member, group in relations:
            if (member in met_proteins and 'RXN' in group):
                # Add a new reaction key and uniprot
                if not group in outdict:
                    # This and the next 'set' were lists in the original incarnation, but it caused
                    # memory leaks when many uniprot accessions were present in a reaction (identical seguids)
                    # or when multiple pathways were nested within each other.
                    outdict[group] = set()
                outdict[group] = outdict[group].union(set(metacyc_uniprot_dict[member]))
            elif member in outdict:
                # Add the higher members 
                if not group in outdict:
                    outdict[group] = set()
                outdict[group] = outdict[group].union(set(outdict[member]))
        # Remove duplicate uniprot accessions.
        outdict = dict([(rxn, list(set(values))) for rxn,values in outdict.items()])
        output_dict = {}
        # Separate out reactions and pathways
        output_dict['reaction_members'] = dict([(key, val) for key,val in outdict.items()
                                                if 'RXN' in key])
        output_dict['pathway_members'] = dict([(key, val) for key,val in outdict.items()
                                                if not 'RXN' in key])
        return output_dict

    @classmethod
    def select_specific_db_for_protein(cls, uniprot_match_results):
        '''Returns the protein product matching a preferred database
        and the database.'''
        # Basically if a metacyc protein product is from one of these dbs
        # choose this over others.
        preferred_dbs = ['ecoli', 'human']
        for protein_id, db in uniprot_match_results:
            if db in preferred_dbs:
                return protein_id, db
        return protein_id, db

    @classmethod
    def obtain_pathways_and_reactions_from_relationships(cls, met_proteins,
                                                         metacyc_uniprot_dict, relationships):
        '''Compared to the function convert_relationships_to_dict, this one only looks
        at reactions and pathways and returns their id and common names. Useful for large numbers
        of relationships.'''
        mapped_dict = cls.map_metacyc_ids_to_reactions(met_proteins, metacyc_uniprot_dict,
                                                       relationships)
        reactions = mapped_dict['reaction_members'].keys()
        pathways = mapped_dict['pathway_members'].keys()
        #Get reaction and pathway objects
        reactions = MetacycReaction.objects.filter(rxn_id__in = reactions).distinct('db')
        pathways = MetacycPathway.objects.filter(pwy_id__in = pathways)
        #Removing duplicates (TODO: Ideally use the 'distinct' db function, but it didn't work at first pass)
        reaction_set = set()
        reaction_list = []
        pathway_set = set()
        pathway_list = []        
        # for reactions and identically for pathways
        for reaction in reactions:
            if not reaction.id in reaction_set:
                reaction_set.add(reaction.id)
                reaction_list.append(reaction)
        for pathway in pathways:
            if not pathway.id in pathway_set:
                pathway_set.add(pathway.id)
                pathway_list.append(pathway)
        # Now map the reaction and pathway objects to their members from mapped_dict
        output_dict = {}
        output_dict['reactions'] = {}
        output_dict['pathways'] = {}
        for reaction in reaction_list:
            output_dict['reactions'][reaction] = mapped_dict['reaction_members'][reaction.id]
        for pathway in pathway_list:
            output_dict['pathways'][pathway] = mapped_dict['pathway_members'][pathway.id]
        return output_dict
        

    @classmethod
    def convert_relationships_to_dict(cls, relationships, filter_to_db):
        '''Converts relationships from sql query to a dictionary.
        If there is no filter_to_db then return the dictionary for all
        metacyc databases'''
        output_dict = {}
        #Empty lists
        output_dict['reactions'] = []
        output_dict['gene'] = []
        output_dict['protein'] = []
        output_dict['complexes'] = []
        output_dict['pathways'] = []
        #Obtain a set of all objects
        object_set = set()
        for member, group in relationships:
            object_set.add(member)
            object_set.add(group)
        object_list = list(object_set)
        #Filter this object to get reactions, pathways, complexes and genes
        for obj in object_list:
            if 'CPLX' in obj:
                cplxes = MetacycProtein.objects.filter(id=obj, db=filter_to_db)
                if cplxes:
                    cplx = cplxes[0]
                    # Filter to the protein db 
                    output_dict['complexes'].append(
                        {'id' : cplx.id,
                         'common_name' : cplx.common_name,
                         'db' : cplx.db})
            elif 'PWY' in obj:
                pwys = MetacycPathway.objects.filter(pwy_id=obj, db=filter_to_db)
                if pwys:
                    pwy = pwys[0]
                    output_dict['pathways'].append(
                        {'id' : pwy.id,
                         'common_name' : pwy.common_name,
                         'db': pwy.db})
            elif 'RXN' in obj:
                reactions = MetacycReaction.objects.filter(rxn_id=obj, db=filter_to_db)
                if reactions:
                    reaction = reactions[0]
                    output_dict['reactions'].append(
                        {'id' : reaction.id,
                         'common_name' : reaction.common_name,
                         'systematic_name' : reaction.systematic_name,
                         'ec_number' : reaction.ec_number,
                         'db' : reaction.db})
            else:
                #For everything else, specifically for genes.
                default_obj = MetacycNode.objects.filter(id=obj, db=filter_to_db)
                if not default_obj:
                    continue
                default_obj = default_obj[0]
                if default_obj.type == 'protein':
                    continue
                if default_obj.type == 'gene':
                    output_dict['gene'].append({'id' : default_obj.id,
                                                'common_name' : default_obj.common_name,
                                                'db' : default_obj.db})
        return output_dict

    def __unicode__(self):
        return unicode(self.common_name)
        
    class Meta:
        db_table = u'metacyc\".\"metacyc_node'

class MetacycToPfacts(models.Model):
    id = models.IntegerField(primary_key=True)
    uniprot = models.ForeignKey('UniProt')
    member = models.ForeignKey('MetacycNode', related_name='mtp_member')
    group = models.ForeignKey('MetacycNode', related_name='mtp_group')

    class Meta:
        db_table = u'metacyc_to_pfacts'

########## end of metacyc classes
        
class Family(models.Model):
    id = models.AutoField(primary_key=True)
    seed_sequence_header = models.ForeignKey(SequenceHeader, null=True)
    average_chars = models.IntegerField(null=True)
    # for caching...
    #most_recent_common_taxid = models.ForeignKey(
    #    UniProtTaxonomy, db_column='most_recent_common_taxid')
    # for caching...
    #viral_most_recent_common_taxid = models.ForeignKey(
    #    UniProtTaxonomy, db_column='viral_most_recent_common_taxid',
    #    related_name="family_viral_most_recent_common_taxid")
    length_shortest = models.IntegerField(null=True)
    length_longest = models.IntegerField(null=True)
    minimum_identities = models.FloatField(null=True)
    mean_identities = models.FloatField(null=True)
    gaps = models.FloatField(null=True)
    columns_w_blosum62_lt_0 = models.IntegerField(null=True)
    longest_deletion_length = models.IntegerField(null=True)
    longest_insertion_length = models.IntegerField(null=True)
    build_database_source = models.TextField()
    build_alignment_notes = models.ForeignKey(BuildAlignmentNotes, null=True)
    family_specific_evalue_criterion = models.FloatField(null=True)
    family_specific_sw_method = models.CharField(max_length=14, null=True)
    author = models.CharField(max_length=255)
    notes = models.TextField(null=True)
    build_date = models.DateField()
    status = models.CharField(max_length=8)
    private = models.BooleanField()
    gathering_method = models.TextField()
    family_type = models.ForeignKey(FamilyType)
    partition = models.CharField(max_length=8)
    #tree_node_consensus_sequence = models.ForeignKey('TreeNodeConsensus')
    canonical_tree = models.ForeignKey('Tree',
        related_name="family_canonical_tre")
    author = models.ForeignKey(User, null=True)
    active = models.BooleanField()
    score = models.IntegerField(null=True)
    score_updated = models.DateTimeField()


    @classmethod
    def find_by_family_id(cls, id):
        try:
            family = cls.objects.get(id = id)
            if (family.status == 'bad') or not (family.active):
                raise cls.DoesNotExist
            return family
        except cls.DoesNotExist:
            raise

    @classmethod
    def get_families_containing_sequence(self, seq):
        # TODO - check to make sure hash didn't collide -> check seq also.
        # seq is ONLY an amino acid pattern, not the defline.
        hashobj = hashlib.sha1()
        hashobj.update(seq)
        seqhash = base64.b64encode(hashobj.digest()).strip('=')
        return self.objects.filter(trees__treenodes__sequence_header__uniprot__seguid=seqhash).distinct()


    def get_accession(self):
        """Return accession in correct format"""

        return u'bpg%07d' % self.id

    def get_num_informative_members(self):
        # This function returns number of TreeNode objects corresponding to leaves that
        # are in this family and have some information in their sequence header.
        # TODO make this a more efficient db call.
        return len(self.get_informative_members())

    def get_informative_members(self):
        # This function returns number of TreeNode objects corresponding to leaves that
        # are in this family and have some information in their sequence header.
        # TODO make this a more efficient db call.
        return self.canonical_root_node().get_informative_members()

    def get_alignment_length(self):
        sequences = TreeNodeAlignment.objects.filter(tree_node__exact = self.canonical_root_node())
        if not sequences:
            return 0
        sequence_string = sequences[0].aligned_sequence.chars
        return len(str(sequence_string).translate(trivial_translation, dotlowercase))

    def get_pfam_accession(self):
        if self.family_type.id == 'C':
            return self.get_pfams().pop()[0].accession
        else:
            return []

    def get_score(self):
        if self.score:
            return self.score
        elif self.family_type.id == 'C': #pfam
            # This is from Kimmen's email
            #score (PHOG) = SP_wt*(Nseqs in SwissProt) + GOEXP_wt * (Nseqs having experimentally supported GO fn/biological process or cellular location) + Taxa_wt * (Number of unique taxa)
            #
            #We could try the following
            #SP_wt = 10
            #GOEXP_wt = 2
            #Taxa_wt = 1
            self.score = 10 * self.canonical_tree.number_of_sequences_in_swissprot() + 2 * self.canonical_tree.number_of_go_evidence_supported_sequences() + 1 * self.canonical_tree.number_of_taxa()
        else: #ghg
            self.score = self.canonical_tree.number_of_taxa()

        self.score_updated = datetime.now()
        self.save()

        return self.score

    def get_superorthologous_phogs(self):
        return Phog.objects.filter(family = self.id, is_superorthologous = True)

    def get_superorthologous_phog_rows(self):
        return [ PHOGRow(phog.tree_node, OrthologTypes.SuperOrtholog , 0.0) for phog in self.get_superorthologous_phogs() ]

    # TODO Make these sourced from our local .phobius files, not a webserver request.  Will make family page load faster.

    def get_phobius_predictions(self):
        (type, start, end) = get_phobius_tmhelix.get_phobius_results_for_family(self.get_accession())

        (tmstart, tmend) = get_phobius_tmhelix.get_transmembrane_regions(type, start, end)
        (sstart, send, nstart, nend, hstart, hend, cstart, cend) = \
            get_phobius_tmhelix.get_signal_peptide_regions(type, start, end)
        return (tmstart, tmend, sstart, send, nstart, nend, hstart, hend, cstart, cend)

    def parse_data_obtained_from_cache(self, data):
        def sort_ecs_by_number(x, y):
            tokens = [x[1].split("."), y[1].split(".")]
            temps = [list(), list()]
            for i in range(len(tokens)):
                for token in  tokens[i]:
                    if token != "-":
                        if token[0].isalpha():
                            temps[i].append(int(token[1:]))
                        else:
                            temps[i].append(int(token))
                    else:
                        temps[i].append(0)
            if temps[0] > temps[1]:
                return 1
            else:
                return -1

        parsed_data = {}
        go_data, ecs, description, family_taxon = data.split(TreeNode.delimiters['field_delimiter'])
        go_data_list = [['Biological Process', []], ['Molecular Function', []], ['Cellular Component', []]]
        categories = go_data.split(TreeNode.delimiters['go_category_delimiter'])
        for i in range(len(categories)):
            records = categories[i].split(TreeNode.delimiters['record_delimiter'])
            for record in records:
                if not records[0] == '':
                    go_data_list[i][1].append(record.split(TreeNode.delimiters['component_delimiter']))
        parsed_data['go_data'] = go_data_list

        ec_links = []
        ec_url = "'http://www.ebi.ac.uk/intenz/query?cmd=SearchEC&ec="
        if ecs:
            ec_records = ecs.split(TreeNode.delimiters['record_delimiter'])
            ec_records = [record.split(TreeNode.delimiters['component_delimiter']) for record in ec_records]
            ec_records.sort(sort_ecs_by_number)
            for record in ec_records:
                ec_description, number, experimental = record
                ec_number_cleaned = number.replace(".-", "")
                if bool(int(experimental)):
                    experimental = " (In Brenda)"
                else:
                    experimental = ""
                ec_number_link = "<a href=%s%s'>%s</a>" % (ec_url,
                        ec_number_cleaned, number)
                ec_text_link = "<a href=%s%s'>%s</a> %s" % (ec_url,
                        ec_number_cleaned, ec_description, experimental)
                ec_links.append([mark_safe(ec_number_link), mark_safe(ec_text_link)])

        parsed_data['ec_links'] = ec_links
        parsed_data['description'] = description
        parsed_data['family_taxon'] = family_taxon

        return parsed_data

    def get_summary_data_from_cache(self):
        summary_data = {}
        tree_node = self.canonical_root_node()
        family_accession = self.get_accession()
        family_dir = get_dir_of_family_accession(family_accession)

        data = tree_node.get_data_from_cache()

        go_data, ecs, description, taxa = data.split(TreeNode.delimiters['field_delimiter'])
        pfams = self.get_pfams_for_caching()
        go_data_list = [['Biological Process', []], ['Molecular Function', []], ['Cellular Component', []]]
        categories = go_data.split(TreeNode.delimiters['go_category_delimiter'])
        for i in range(len(categories)):
            records = categories[i].split(TreeNode.delimiters['record_delimiter'])
            for record in records:
                go_data_list[i][1].append(record.split(TreeNode.delimiters['component_delimiter']))
        summary_data['go_data'] = go_data_list

        ec_links = []
        ec_url = "'http://www.ebi.ac.uk/intenz/query?cmd=SearchEC&ec="
        if ecs:
            ec_records = ecs.split(TreeNode.delimiters['record_delimiter'])
            for record in ec_records:
                ec_record = record.split(TreeNode.delimiters['component_delimiter'])
                ec_links.append("<a href=%s%s'>%s(%s)</a>" % (ec_url,ec_record[1],
                    ec_record[0], ec_record[1]))

        summary_data['tree_methods'] = set()
        for method in ['nj', 'ml']:
            if os.path.exists(os.path.join(family_dir, "%s.%s" % (family_accession, method))):
                summary_data['tree_methods'].add(method)


        summary_data['treeMethod'] = None
        if summary_data['tree_methods']:
            if 'nj' in summary_data['tree_methods']:
                summary_data['treeMethod'] = 'nj'
            elif 'ml' in summary_data['tree_methods']:
                summary_data['treeMethod'] = 'ml'
            else:
                pass

        summary_data['max_phyloscope_tree_size'] = max_phyloscope_tree_size
        summary_data['small_tree'] = (tree_node.num_nonredundant_sequences <= max_phyloscope_tree_size)
        summary_data['ec_links'] = mark_safe(" ".join(ec_links))
        summary_data['ebi_link'] = "http://www.ebi.ac.uk/QuickGO/GTerm"
        summary_data['pfam_links'] = mark_safe(pfams)
        summary_data['family_taxon'] = taxa
        summary_data['num_nonredundant_sequences'] = len(tree_node.get_included_leaves())
        summary_data['description'] = description
        summary_data['alignment_length'] = '' # TODO: ask RSD how this is computed
        return summary_data

    def get_family_name(self):
        """Generates the family name for a PhyloFacts family.
        1. Family should be named with trusted editbook names first.
        2. PhyloFacts family should be named after a pfam family if it is > 70%
            of the family consensus sequence length, or if it is a pfam family.
        3. If neither of these holds, name the family the way it would have been
            named before (by description voting consensus).
        """
        family_accession = self.get_accession()
        tree_node = self.canonical_root_node()
        try:
            if self.family_type.id == 'C':
                pfam = self.get_pfams().pop()[0]
                return "%s: %s" % (pfam.name, pfam.accession)
            else:
                for pfam in self.get_pfams():
                    if (pfam[2]-pfam[1]) > (0.7*self.get_alignment_length()):
                        return "%s (related)" % (pfam[0].name)
                # Just name it the way it would have been named
                return "%s (related)" % (tree_node.get_description_for_caching())
        except:
            return family_accession

    def get_family_description(self):
        # TODO merge all naming functions together.  This is using a different naming function for the sequence page (uniprot_detail.html) vs the family page.  Change this, make everything consistent.
        tree_node = self.canonical_root_node()
        if self.family_type.id == 'G':
            return "%s" % (tree_node.get_description_for_caching())
        else:
            return "%s" % (self.get_pfams().pop()[0].description)

    def get_family_descriptions(self):
        # This function gets ALL the family descriptions sorted by score and returns the top 10.
        if self.family_type.id == 'G':
            # I am a horrible person for this - JSD
            try:
                tree_node = self.canonical_root_node()
                return tree_node.get_description_for_caching(returnAllDescriptions=True)[:50]
            except:
                tree_node = self.trees.get(method='nj').treenodes.get(left_id=1)
                return tree_node.get_description_for_caching(returnAllDescriptions=True)[:50]
        else:
            return "%s" % (self.get_pfams().pop()[0].description)

    def map_member_data_to_biocyc_objects(self, biocyc_dict):
        '''Takes a dictionary containing uniprot accessions
        for reactions and pathways and maps member data onto
        them after sorting the dictionary based on number of accessions.'''
        output_dict = {}
        output_dict['reactions'] = []
        output_dict['pathways'] = []
        if not (biocyc_dict['reactions'] or biocyc_dict['pathways']):
            return output_dict
        member_data = self.get_member_data()
        #Sorting in a descending order (-1*) based on number of members
        sorted_biocyc_reaction_tuple = sorted(biocyc_dict['reactions'].items(),
                                             key=lambda e: -1*len(e[1]))
        sorted_biocyc_pathway_tuple = sorted(biocyc_dict['pathways'].items(),
                                             key=lambda e: -1*len(e[1]))
        #Mapping member data
        for biocyc_obj, uniprot_list in sorted_biocyc_reaction_tuple:
            member_dict = {'members' : [val for val in member_data if
                                        val[1] in uniprot_list]}
            output_dict['reactions'].append((biocyc_obj, member_dict))
        for biocyc_obj, uniprot_list in sorted_biocyc_pathway_tuple:
            member_dict = {'members' : [val for val in member_data if
                                        val[1] in uniprot_list]}
            output_dict['pathways'].append((biocyc_obj, member_dict))
        return output_dict

    #@memoized 
    def get_consensus_sequence(self):
        #Because the TreeNodeConsensus isn't populated fully, we have to get this from file.
        # TODO: put this information in the database when we encounter an unpopulated family
        family_dir = get_dir_of_family_accession(self.get_accession())
        if os.path.exists(family_dir):
            consensus_paths = os.path.join(family_dir, "*.con.fa")
            consensus_sequence_path = glob.glob(consensus_paths).pop()
            if not consensus_sequence_path:
                return ("","")
            file = open(consensus_sequence_path, 'r')
            cs = file.readlines()
            file.close()
            consensus_sequence = ""
            for line in cs:
                if line[0] == ">":
                    consensus_header = line.replace("\n", "")
                else:
                    consensus_sequence += line.replace("\n", "")
            return (consensus_header, consensus_sequence)
        else:
            return ("","")
            
        
    def get_summary_data(self):
        summary_data = {}
        tree_node = self.canonical_root_node()
        family_accession = self.get_accession()
        family_dir = get_dir_of_family_accession(family_accession)
        summary_data['tree_methods'] = set()
        for method in ['nj', 'ml']:
            if os.path.exists(os.path.join(family_dir, "%s.%s" % (family_accession, method))):
                if os.path.getsize(os.path.join(family_dir, "%s.%s" % (family_accession, method))) > 0:
                    summary_data['tree_methods'].add(method)

        # treeMethod refers to the method that will be used by phyloscope
        summary_data['treeMethod'] = None
        if summary_data['tree_methods']:
            if 'ml' in summary_data['tree_methods']:
                summary_data['treeMethod'] = 'ml'
            elif 'nj' in summary_data['tree_methods']:
                summary_data['treeMethod'] = 'nj'
            else:
                pass

        if len(tree_node.get_included_leaves()) <= max_phyloscope_tree_size:
            summary_data['small_tree'] = True
        else:
            summary_data['small_tree'] = False

        # for pfam families, the description is the pfam domain
        # replaced the pfam_description and other naming
        # protocols to naming from sqlite db
        summary_data['family_name'] = None
        summary_data['pfam_description'] = None
        summary_data['isPfam'] = (self.family_type.id == 'C')
        if self.family_type.id == 'C':
            pfam =  self.get_pfams().pop()[0]
            pfam_description = "%s (%s)" % (pfam.name, pfam.accession)
            summary_data['pfam_description'] = pfam_description
            summary_data['pfam_name'] = pfam.name
            summary_data['pfam_longname'] = pfam.description
        try:
            data = FamilyCache.objects.get(family_id=self.id).cached_data
            cached_data = self.parse_data_obtained_from_cache(data)
            summary_data['ec_links'] = cached_data['ec_links']
            summary_data['go_data'] = cached_data['go_data']
        #    summary_data['description'] = cached_data['description']
            summary_data['family_taxon'] = cached_data['family_taxon']
        except FamilyCache.DoesNotExist:
            summary_data['ec_links'] = ['', '']
            summary_data['go_data'] = tree_node.go_data_sorted_by_priority()
            #summary_data['description'] = self.get_family_descriptions()
            #summary_data['family_taxon'] = phog_row.taxon_name
            summary_data['family_taxon'] =  self.get_taxon()
        try:
            summary_data['description'] = self.get_family_descriptions()
        except:
            summary_data['description'] = ""

        summary_data['max_phyloscope_tree_size'] = max_phyloscope_tree_size
        summary_data['ebi_link'] = "http://www.ebi.ac.uk/QuickGO/GTerm"
        summary_data['pfam_links'] = make_pfam_links(self)
        summary_data['num_nonredundant_sequences'] = self.get_num_informative_members()
        #Obtaining biocyc reactions and pathways
        biocyc_dict = self.canonical_root_node().get_biocyc_reactions_and_pathways()
        biocyc_member_dict = self.map_member_data_to_biocyc_objects(biocyc_dict)
        summary_data['pathways'] = biocyc_member_dict['pathways']
        summary_data['reactions'] = biocyc_member_dict['reactions']
        #
        summary_data['alignment_length'] = self.get_alignment_length()
        pfamlist = list(self.get_pfams())
        sequence_hmms = []
        hm = []
        hm = self.get_hmms()
        for h in hm:
            if (h.hmm.pfam.accession in [x[0].accession for x in pfamlist]):
                sequence_hmms.append(h)

        hmms = []
        color_options = '3456789abcde'
        for shmm in sequence_hmms:
            color = ['#']
            for i in range(6):
                color.append(random.choice(color_options))

            display_styles_map = {"." : "jagged", "]" : "curved", "[" : "curved"}
            startStyle, endStyle = [display_styles_map[char] for char in shmm.match_type]

            hmms.append({
                "modelStart" : shmm.hmm_start,
                "modelEnd" : shmm.hmm_end,
                "aliStart" : shmm.sequence_start,
                "aliEnd" : shmm.sequence_end,
                "startStyle" : startStyle,
                "endStyle" : endStyle,
                "color" : ''.join(color),
                "modelLength" : shmm.hmm.length,
                "e_value" : shmm.e_value,
                "pfam_accession" : shmm.hmm.pfam.accession,
                "pfam_description" : shmm.hmm.pfam.description,
                "pfam_source" : "Pfam (version %s)" % str(shmm.hmm.pfam.overall_pfam_version),
                "pfam_name" : shmm.hmm.pfam.name
                })
        summary_data['hmms'] = hmms

        consensus_sequence_length = 0
        try:
            consensus_sequence_length = len(TreeNodeConsensus.objects.get(tree_node = self.canonical_root_node()).sequence.chars)
        except:
            pass

        summary_data["consensus_sequence_length"] = consensus_sequence_length

        (tmstart, tmend, sstart, send, nstart, nend, hstart, hend, cstart, cend) = self.get_phobius_predictions()

        summary_data['transmembrane_regions'] = zip(tmstart,tmend,["#FF0000" for i in tmstart])
        summary_data['signal_regions'] = zip(sstart,send,["#00FF00" for i in sstart])
        summary_data['n_regions'] = zip(nstart,nend,["#FFFF00" for i in nstart])
        summary_data['h_regions'] = zip(hstart,hend,["#0000FF" for i in hstart])
        summary_data['c_regions'] = zip(cstart,cend,["FFA500" for i in cstart])

        return summary_data

    def get_long_go_summary(self, go_component):
        try:
            treenode = self.canonical_root_node()
        except TreeNode.DoesNotExist:
            return []
        return filter(lambda x:x, [treenode.get_go_summary(go_component,
            experimental=ex, ortholog_type=OrthologTypes.SuperOrtholog,
            threshold=0.0) for ex in [True,False]])

    def canonical_root_node(self):
        return self.canonical_tree.treenodes.get(left_id=1)

    def get_pfams(self):
        # This function returns the pfam domains associated with a family.
        # To match the pfam website, we pick the domain with the most significant hit
        # along the protein.  We then remove any domains that overlap with
        # the chosen one.  Then we choose the next most significant hit, and remove
        # all domains overlapping with that one.  This process continues
        # until there are no overlapping domains.

        def __remove_overlapping_pfams__(pfam,pfamlist):
            retlist = []
            template = set(range(pfam[1],pfam[2]))
            for pf in pfamlist:
                template2 = set(range(pf[1], pf[2]))
                if (len(template.intersection(template2)) == 0):
                    retlist.append(pf)
            return retlist

        def __contains_overlapping_pfams__(pfamlist):
            # TODO Make this search less sucky, shouldn't be very hard.
            for pf in pfamlist:
                for pf2 in pfamlist:
                    if pf!=pf2:
                        a = set(range(pf[1], pf[2]))
                        b = set(range(pf2[1],pf2[2]))
                        if len(a.intersection(b)) != 0:
                            return True
            return False

        # I am a horrible person - JSD
        # We had a lot of NJ trees without phogs set as the canonical trees
        # This breaks using the canonical tree for things like printing the right tree
        # So I added this horribleness.
        # Moving some shit to tree node would probably help.
        pfam = self.canonical_root_node().get_seed_pfam()
        if not pfam:
            pfam = self.trees.get(method='nj').treenodes.get(left_id=1).get_seed_pfam()

        if pfam is not None:
            return set([(pfam, 1, -1)])
        pfamlist = [(sequence_hmm.hmm.pfam, sequence_hmm.sequence_start,
                        sequence_hmm.sequence_end, sequence_hmm.e_value)
                        for sequence_hmm in self.get_hmms()]
        # are there any overlap issues?
        if not __contains_overlapping_pfams__(pfamlist):
            return set([(pfam[0],pfam[1],pfam[2]) for pfam in pfamlist])
        # TODO modify this section to use python's sort rather than this
        # there a lot that can change to be more pythonic here.
        filtered = []
        while True:
            min_eval = 1.0
            for pfam in pfamlist:
                if (pfam[3]<min_eval):
                    min_eval = pfam[3]
            for pfam in pfamlist:
                if pfam[3] == min_eval:
                    pfamlist = __remove_overlapping_pfams__(pfam,pfamlist)
                    filtered.append(pfam)
            if not(__contains_overlapping_pfams__(pfamlist)):
                pflist = [(pfam[0], pfam[1], pfam[2]) for pfam in pfamlist]
                filteredlist = [(f[0], f[1], f[2]) for f in filtered]
                return set(pflist).union(set(filteredlist))

    def get_pfams_for_caching(self):
        return str(make_pfam_links(self))

    def get_member_data(self):
        """Returns member data that will be used by the views and for api calls.
        This function has been moved from views.py
        """
        cursor = connection.cursor()
        sequence_sql = """
        select uniprot.id, uniprot.accession, uniprot.uniprot_identifier,
          uniprot.in_swissprot_f, uniprot.de, sequence_header.header,
          uniprot_taxonomy.scientific_name, uniprot_taxonomy.common_name, uniprot_taxonomy.id
        from tree_node, family, tree, sequence_header, uniprot, uniprot_taxonomy
        where uniprot_taxonomy.id = uniprot.taxon_id and
          uniprot_taxonomy.id = sequence_header.taxon_id and
          sequence_header.uniprot_id = uniprot.id and
          tree_node.sequence_header_id = sequence_header.id and
          family.canonical_tree_id = tree.id and
          tree_node.tree_id = tree.id and
          tree_node.sequence_header_id is not Null and
          family.id = %s AND
          family.active IS TRUE
          """
        # Get sequence data
        cursor.execute(sequence_sql, [self.id])
        results = cursor.fetchall()
        # Modify results to include the words 'Swissprot' or 'Trembl'
        # TODO change this to accept non-uniprot sequence, this assumes its either in swissprot, or its in trembl.
        reslist = []
        for result in results:
            if result[3]:
                reslist.append((result[0], result[1], result[2], "SwissProt",
                                result[4], result[5], result[6], result[7], result[8]))
            else:
                reslist.append((result[0], result[1], result[2], "TrEMBL",
                                result[4], result[5], result[6], result[7], result[8]))
        results = reslist
        uniprot_ids = [result[0] for result in results]
        go_sql = """
        select uniprot_id, go_term.name, go_term.term_type, go_term.acc, go_evidence.code
        from uniprot_go, go_term, go_evidence
        where go_term.id = uniprot_go.go_term_id and
          go_evidence.id = uniprot_go.go_evidence_id and
          uniprot_id in %s order by uniprot_id, term_type
          """
        # Get GO results if there are uniprot ids
        go_results = []
        if uniprot_ids:
            cursor.execute(go_sql % str(tuple(uniprot_ids)))
            go_results = cursor.fetchall()
        go_map = {}

        for go_result in go_results:
            uniprot_id = go_result[0]
            go_category = go_result[2]
            go_data = (go_result[1], go_result[3], go_result[4])
            go_categories = ["molecular_function", "biological_process", "cellular_component"]
            go_categories_set = set(go_categories)
            if uniprot_id not in go_map.keys():
                go_map[uniprot_id] = [ ['Molecular Function(s)', []],
                                       ['Biological Process(es)', []],
                                       ['Cellular Component(s)', []]]
            if go_category in go_categories_set:
                category_index = go_categories.index(go_category)
                go_map[uniprot_id][category_index][1].append(go_data)

        # TODO: deal with deprecated sequences
        # TODO: deal with non-uniprot sequences
        # TODO: show fasta sequence - Identifier, Accession, Sequence (goes
        # to new page)
        # show seed sequence

        final_results = []
        for result in results:
            result = list(result)
            # if there is go data ...
            if result[0] in go_map:
                result.append(go_map[result[0]])
            else:
                result.append([])
            final_results.append(result)
        return final_results
    
    def get_hmms(self):
        """
            This gets the Pfam HMMs associated with the family consensus
            sequence.  It actually gets the SeqeunceHMM objects that
            describe the match: the e-value, the start and end indices,
            etc.
        """
        # I am a horrible person - JSD
        try:
            return SequenceHMM.objects.filter(
                sequence=self.canonical_root_node(
                    ).treenodeconsensus_set.all()[0].sequence,
                hmm__pfam__isnull=False,
            )
        except:
            return SequenceHMM.objects.filter(
                sequence=self.trees.get(method='nj').treenodes.get(left_id=1
                    ).treenodeconsensus_set.all()[0].sequence,
                hmm__pfam__isnull=False,
            )

    def get_pmid_links(self):
        pmids = self.canonical_root_node().get_pmids(
                                ortholog_type = OrthologTypes.PHOG_T_Custom,
                                threshold=10000.0)
        pmid_links = mark_safe(', '.join([pmid_link(pmid) for pmid in pmids]))
        return pmid_links

    def get_ec_links(self):
        experimental_ecs = self.canonical_root_node().get_ecs(
                              experimental=True,
                              ortholog_type = OrthologTypes.PHOG_T_Custom,
                              threshold=10000.0)
        nonexperimental_ecs = self.canonical_root_node().get_ecs(
                  experimental=False,
                  ortholog_type = OrthologTypes.PHOG_T_Custom,
                  threshold=10000.0)
        ec_links = mark_safe(','.join([uniprot_ec_link(ec, experimental=True)
                                     for ec in experimental_ecs] +
                                      [uniprot_ec_link(ec, experimental=False)
                                      for ec in nonexperimental_ecs]))
        return ec_links

    def get_taxon(self):
        return self.canonical_root_node().get_taxon(ortholog_type =
            OrthologTypes.PHOG_T_Custom, threshold = 10000.0)

    def get_permissive_phog(self):
        return self.canonical_root_node().get_phog(10000.0)

    def get_alignment_type_str(self):
        return self.family_type.id == 'G' and u'Global' or u'Local'

    def get_family_type_str(self):
        return self.family_type == 'global_homology' and 'GHG' or \
            self.family_type == 'domain' and u'Domain' or u'Other'

    def is_global_homology(self):
        return self.family_type == 'global_homology'

    def get_phyloscope_url(self):
        return family_phyloscope_url(self.id)

    def get_phyloscope_link(self):
        return family_phyloscope_link(self.id)

    @models.permalink
    def get_absolute_url(self):
        return('family', ['bpg%07d' % self.id])

    def get_link(self):
        return "<a href='%s'>bpg%07d</a>" % (
            self.get_absolute_url(),
            self.id
            )


    ''' Returns the path to the subtree hmm file
    Create the file if necessary '''
    def subtree_hmm_file(self, no_build=False):
        #TODO: Fixed the race conditions (mostly), this should be refactored, because its kind of crappy codewise
        file_location = "/clusterfs/vasudha/bpg/fat-cat/%s.hmm" % self.get_accession()
        queue = FatcatHMMBuildQueue.objects.filter(family = self)
        if queue:
            # increment the number of processes waiting for this family
            queue.count += 1
            queue.save()
            # we are waiting on this guy to finish
            while queue:
                # wait 10 seconds and poll again
                time.sleep(10)
                queue = FatcatHMMBuildQueue.objects.filter(family=self)
            return file_location
        else:        
            if os.path.exists(file_location):
                return file_location
            if no_build:
                return None
            # we must build it
            queue = FatcatHMMBuildQueue.objects.create(family = self)
            while queue:
                # wait 10 seconds and poll again
                time.sleep(10)
                queue = FatcatHMMBuildQueue.objects.filter(family=self)    
            return file_location

    def __unicode__(self):
        return unicode('bpg%07d' % self.id)

    __str__ = __unicode__

    class Meta:
        db_table = u'family'

seqhdr_re = re.compile('(SEQHDR[0-9]*)')
capturing_seqhdr_re = re.compile('(SEQHDR([0-9]*))')
class Tree(models.Model):
    id = models.AutoField(primary_key=True)
    family = models.ForeignKey(Family, related_name='trees')
    method = models.CharField(max_length=21)
    reciprocal_longest_distance = models.FloatField(null=True)
    is_rsd_rooted = models.BooleanField(default=False)

    def write_newick_to_handle(self, outf):
        family_accession = self.family.get_accession()
        family_dir = get_dir_of_family_accession(family_accession)
        if not os.path.exists(family_dir):
            raise IOException
        filename = "%s.%s" % (family_accession, self.method)
        tree_path = os.path.join(family_dir, filename)
        if not os.path.exists(tree_path):
            raise IOException
        # Now tree_path exists
        f = open(tree_path)
        tree_str = f.read()
        f.close()
        seqhdr_tree_tokens = seqhdr_re.split(tree_str)
        for i in range(len(seqhdr_tree_tokens)):
            if i % 2 == 0:
                outf.write(seqhdr_tree_tokens[i])
            else:
                m = capturing_seqhdr_re.match(seqhdr_tree_tokens[i])
                sequence_header = SequenceHeader.objects.get(id = m.group(2))
                outf.write(sequence_header.get_escaped_header_for_tree())

    def get_phogs(self, threshold=0.0):
        return TreeNode.objects.filter(tree = self,
                  duplication_distance__gte = threshold,
                  greatest_duplication_distance_of_maximal_descendant__lt
                  = threshold).order_by('left_id')

    def number_of_sequences_in_swissprot(self):
        return UniProt.objects.filter(
            sequenceheader__treenode__tree = self.id,
            in_swissprot_f = True
        ).count()

    def number_of_go_evidence_supported_sequences(self):
        return Annotation.objects.filter(
                uniprot__sequenceheader__treenode__tree = self,
                evidence__priority__lte = 5
            ).values('uniprot_id').distinct('uniprot_id').count()


    def number_of_taxa(self):
        return UniProtTaxonomy.objects.filter(
            uniprot__sequenceheader__treenode__tree = self.id
        ).distinct().count()

    def root(self):
        return TreeNode.objects.get(tree=self, left_id=1)

    def internal_nodes(self):
        return self.root().get_internal_nodes()

    def __unicode__(self):
        return unicode('tree: %s' % self.id)

    class Meta:
        db_table = u'tree'


global left_right_ids
global max_left_id
left_right_ids = dict([(key,None) for key in []])
max_left_id = None



'''
                    acc = annotation.go_term.acc
                    name_of_acc[acc] = annotation.go_term.name
                    evidence = annotation.go_evidence.evidence
                    evidence_priority = annotation.go_evidence.priority
                    if acc in gos:
                        current_priority = gos[acc].priority
                        if evidence_priority < current_priority:
                            gos[acc] = annotation.go_evidence
                    else:
                        gos[acc] = annotation.go_evidence
        rs = make_ranked_summary({go_component: gos}, name_of_acc)
        return go_component in rs and rs[go_component] or None

'''

def getTreeNodeQuerySetForTreeNodes(tree_nodes):
    if tree_nodes:
        where_clause = '(tree_node.tree_id, tree_node.left_id) IN (%s)' \
        % ','.join(['(%d,%d)' % (node.tree_id, node.left_id) for node
                                                              in tree_nodes])
    else:
        where_clause = 'FALSE'
    return TreeNode.objects.extra(where=[where_clause])

class TreeNode(models.Model):
    def __init__(self, *args, **kwargs):
        models.Model.__init__(self, *args, **kwargs)
        # these are for caching...
        self.description = ''
        self.accession = ''
        self.taxon = {}
        self.num_nonredundant_sequences = 0
        self.experimental_go_summary = {}
        self.nonexperimental_go_summary = {}
        self.have_summarized_go = {}
        # We have not summarized the GO experimental annotations
        self.have_summarized_go[True] = False
        # We have not summarized the GO nonexperimental annotations
        self.have_summarized_go[False] = False
        self.go_annotation_dict = {}
        for experimental in [True, False]:
            self.go_annotation_dict[experimental] = {}
            for go_component in ['biological_process',
                'molecular_function', 'cellular_component']:
                self.go_annotation_dict[experimental][go_component] = {}
        self.num_contained_taxa = 0
        self.included_leaves = None
        self.have_collected_ecs = False
        self.experimental_ecs = set()
        self.nonexperimental_ecs = set()
        self.have_collected_pmids = False
        self.pmids = set()
        self.has_swissprot = None
        self.have_collected_pdb_chains = False
        self.pdb_chains = set()
        # The variable below will be modified by get_pfam_domain_architecture_in_members
        self.no_pfam_domain_members = None 
        
    delimiters = { 'field_delimiter' : "<<field>>",
        'go_category_delimiter' : "<<go_category>>",
        'record_delimiter' : "<<record>>",
        'component_delimiter' : "<<component>>",
        }

    def msa(self, sequence_lambda=lambda x: x):
        return '\n'.join(
            [ ">%s\n%s" % (tna['sequence_header__header'], sequence_lambda(tna['aligned_sequence__chars']) ) for tna in
            TreeNodeAlignment.objects.filter(
                tree_node__tree__family = self.tree.family,
                sequence_header__id__in = self.get_included_leaves().values('sequence_header_id')
            ).distinct('sequence_header__header').values('sequence_header__header', 'aligned_sequence__chars')] )

    '''Only Upper Case and -'''
    def uc_msa(self):
        return self.msa(sequence_lambda=lambda x: re.sub(r'[a-z.]','',x))

    def is_in_pfam_family(self):
        return (self.tree.family.family_type_id == 'C')

    def is_in_mda_family(self):
        return (self.tree.family.family_type_id == 'G')
   
    def get_consensus_sequence_for_family(self):
        return self.tree.family.get_consensus_sequence()        

    def subtree_hmm(self):
        # Returns the subtree hmm for this node
        try:
            return self._subtree_hmm
        except:
            pass

        (hmm, err) = hmmbuild(self.uc_msa(), self.id)
        self._subtree_hmm = hmm
        return self._subtree_hmm
        
    #@memoized
    def get_treenode_names(self):
        """ Return the top maxNames ranked names for this treenode """
        try:
            return self._get_treenode_names
        except:
            pass

        all_descriptions = self.get_description(returnAll=True, score_pfam=True)
        # This is crappy...should be changed but it's an undertaking, hence this piece of code.
        if type(all_descriptions) == type(u''): 
            all_descriptions = [(all_descriptions, 1, 1)]
        sorted_descriptions = sorted(all_descriptions, key=lambda score: score[1], reverse=True)
        top_ranked_description = sorted_descriptions[0][0]
        treenode_name = "No Node Description"
        treenode_accession = "N/A"
        treenode_type = "N/A"
        if self.is_family_root():
            treenode_accession = self.tree.family.get_accession()
            if self.is_in_pfam_family():
                pfam = self.tree.family.get_pfams().pop()[0]
                treenode_name = pfam.name
                treenode_type = "Phylofacts-Pfam family"
            else:
                treenode_name = top_ranked_description + "(related)"
                treenode_type = "PhyloFacts-MDA family"
        elif self.is_leaf():
            treenode_accession = "LEAF%07d_%05d" % (self.tree.id, self.left_id)
            treenode_name = top_ranked_description
            treenode_type = "Protein"
        else:
            treenode_accession = "NODE%07d_%05d" % (self.tree.id, self.left_id)
            treenode_name = top_ranked_description + "(related)"
            treenode_type = "PhyloFacts Subtree"
        self._get_treenode_names = {'treenode_name':treenode_name, 'treenode_type':treenode_type,
                                    'treenode_accession':treenode_accession, 
                                    'all_descriptions':sorted_descriptions}
        return self._get_treenode_names
     
    #@memoized
    def get_phobius_results_for_family(self):
        accession = 'bpg%07d' % self.tree.family.id
        (type,start,end) = get_phobius_tmhelix.get_phobius_results_for_family(accession)
        (tmstart, tmend) = get_phobius_tmhelix.get_transmembrane_regions(type,start,end)
        (sstart, send, nstart, nend, hstart, hend, cstart, cend) = \
            get_phobius_tmhelix.get_signal_peptide_regions(type, start, end)
        return {'transmembrane': zip(tmstart, tmend), 'signal_peptide': zip(sstart, send)}

    def get_informative_members(self):
        # This function returns number of TreeNode objects corresponding to leaves that
        # are in this family and have some information in their sequence header.
        # TODO make this a more efficient db call.
        try:
            return self._get_informative_members
        except:
            pass
        leaves = self.get_included_leaves()
        self._get_informative_members = [leaf for leaf in leaves if leaf.sequence_header.uniprot is not None]
        return self._get_informative_members

    def get_informative_member_headers(self):
        """ Returns the sequence headers for all nodes under this one, for autocomplete"""
        leaves = self.get_informative_members()
        return [(str(leaf.sequence_header.header), leaf.sequence_header.id) for leaf in leaves]

    def get_uniprot_in_members(self):
        '''Gets the uniprot accessions for members below this node.'''
        informative_members = self.get_informative_members()
        return [leaf.sequence_header.uniprot for leaf in informative_members]
    

    def get_species_tree_json(self):
        '''Gets the json object to populate the species tree'''
        
        leaves = self.get_uniprot_in_members()
        lineages = []
        for leaf in leaves:
            this_lineage = []
            for taxa in leaf.taxon.lineage():
                if taxa.common_name is not None:
                    this_lineage.append("%s (%s)" % (taxa.scientific_name, taxa.common_name))
                else:
                    this_lineage.append(taxa.scientific_name)
            lineages.append(this_lineage)

        output_dict = {'data' : 'root',
                       'attr': {'species_count' : 0,
                                'sequence_count' : 0,
                                'ancestor': ""
                               },
                       'children' : []
                       }
        for lineage in lineages:
                current_loc = output_dict
                current_loc['attr']['sequence_count'] += 1
                current_loc['data'] = "root [%d sequence(s)]" % current_loc['attr']['sequence_count'] 
                for ancestor in lineage:
                    children = [child['attr']['ancestor'] for child in current_loc['children']]
                    if not ancestor in children:
                        current_loc['children'].append({'data' : "%s [1 sequence(s)]" % ancestor,
                                                        'attr': {'species_count' : 1, 
                                                                 'sequence_count' : 1,
                                                                 'ancestor': ancestor
                                                                },
                                                        'children' : []}) 
                        current_loc = current_loc['children'][len(current_loc['children'])-1]
                    else:
                        current_loc = [child for child in current_loc['children']
                                       if child['attr']['ancestor'] == ancestor][0]
                        current_loc['attr']['sequence_count'] += 1
                        current_loc['data'] = "%s [%d sequence(s)]" % (ancestor, current_loc['attr']['sequence_count'])
        return json.dumps(output_dict)

    def get_pfam_domain_architecture_in_members(self):
        '''Calls the pfam_architecture function and returns this for all members, grouping by
        the largest membership.'''
        try:
            return self._get_pfam_domain_architecture_in_members
        except:
            pass

        informative_members = self.get_informative_members()
        domain_architectures = [leaf.sequence_header.uniprot.pfam_architecture() for
                                leaf in informative_members]
        domain_architectures = [arch for arch in domain_architectures if arch]
        num_architectures = len(domain_architectures)
        output_dict = {}
        for arch in domain_architectures:
            domain_architecture = arch[0]
            if not domain_architecture in output_dict:
                output_dict[domain_architecture] = [arch[1]]
            else:
                output_dict[domain_architecture].append(arch[1])
        # To this dictionary add the total number of members, the number with this domain architecture
        # and the percent covered
        output_items = []
        total_num_members = len(informative_members)
        # Modify the instance variable for members with no pfam_domains
        self.no_pfam_domain_members = int(total_num_members)
        for key, value in output_dict.items():
            fraction = '%0.2f' % (len(value)*100/
                                  float(total_num_members))
            fraction_with_pfams = '%0.2f' % (len(value)*100/
                                        float(num_architectures))
            output_items.append((key, value, {'num_this_grp' : len(value),
                                              'fraction' : fraction,
                                              'fraction_with_pfams' : fraction_with_pfams,
                                              'total_num_members' : total_num_members}
                                 ))
            # Update the value for members with no pfam_domain
            self.no_pfam_domain_members -= len(value)
        self.no_pfam_domain_members = (self.no_pfam_domain_members,
                                       '%0.2f' % (self.no_pfam_domain_members * 100/float(total_num_members)))
        self._get_pfam_domain_architecture_in_members = sorted(output_items, key=lambda x: -1*len(x[1]))
        return self._get_pfam_domain_architecture_in_members

    def get_go_data_with_uniprot_members(self):
        '''This might be redundant, however, the difference here is that uniprot objects are also stored for each
        GO accession so that a datatable of members with a particular go function can be displayed.
        When another function with go_data will be used - one which retrieves the uniprot members as well, this function
        can be removed.
        There are two steps in this function. Run the function on a small node to see how this works.
        Create the output dictionary that looks like this
            {"total_go_annotations": 33, "molecular_function": [
            (GO term description, ["GO:0005524", [
            ("IEA", ["Inferred from Electronic Annotation", number_of_IEA_members])],
            [list of uniprot objects]]),
            ...
            ]
            }
        Next order this dictionary, first by evidence codes (using the reversed_order_of_importance list)
        and then by the number of uniprot members for each GO term.
        '''
        go_dict = {'molecular_function': {}, 'cellular_component': {}, 'biological_process': {}, 'total_go_annotations': 0 }
        members = self.get_uniprot_in_members()

                #UniProt.objects.filter(
                #sequenceheader__treenode__tree = self.tree,
                #sequenceheader__treenode__left_id__gte = self.left_id,
                #sequenceheader__treenode__right_id__lte = self.right_id,
            #).all()
        for key in [ 'molecular_function',  'cellular_component',  'biological_process' ]:
            for uniprot in members:
                for a in [ a.annotation for a in uniprot.annotations('go_' + key) ]:
                    if not a['accession'] in go_dict[key]:
                        go_dict['total_go_annotations'] += 1
                        go_dict[key][a['accession']] = { 'accession': a['accession'],
                                                            'name': a['description'],
                                                            'evidence': {a['evidence_code'] : { 'evidence_code': a['evidence_code'], 'name': a['evidence_description'], 'count': 1} }, 
                                                            'members': [uniprot] }
                    elif not a['evidence_code'] in go_dict[key][a['accession']]['evidence']:
                        go_dict[key][a['accession']]['evidence'][a['evidence_code']] = { 'evidence_code': a['evidence_code'], 'name': a['evidence_description'], 'count': 1}
                        go_dict[key][a['accession']]['members'].append(uniprot)
                    else:
                        go_dict[key][a['accession']]['evidence'][a['evidence_code']]['count'] += 1
                        go_dict[key][a['accession']]['members'].append(uniprot)

        output_dict = {'molecular_function': [], 'cellular_component': [], 'biological_process': [], 
                            'total_go_annotations': go_dict['total_go_annotations'] }

        def evidence_value(e):
            if e['evidence_code'] == 'IEA':
                return  e['count']
            return -1 * e['count']

        for key in [ 'molecular_function',  'cellular_component',  'biological_process' ]:
            # sort the terms by the number of sequences annotated with it
            output_dict[key] = sorted(go_dict[key].values(), key=lambda a: sum( [ e['count'] for e in a['evidence'].values() ] ) )
            for term in output_dict[key]:
                term['evidence'] = sorted(term['evidence'].values(), key=evidence_value)

        return output_dict
        
    def get_biocyc_reactions_and_pathways(self):
        '''Returns a dictionary containing pathways and reactions for members belonging to a node.
        '''
        # Get uniprot objects and map to metacyc protein ids
        uniprot_accs = self.get_uniprot_in_members()
        metacyc_ids, metacyc_uniprot_dict = MetacycNode.convert_to_metacyc_product(uniprot_accs)
        output_dict = {'pathways' : None, 'pathway_members' : None,
                       'reactions' : None, 'reaction_members' : None
                       }
        # Get membership tuples for metacyc members, and convert to a dictionary with
        # reactions and pathways
        if metacyc_ids:
            metacyc_relations = MetacycNode.get_membership_recursively(metacyc_ids)
            metacyc_dict = MetacycNode.obtain_pathways_and_reactions_from_relationships(
                metacyc_ids, metacyc_uniprot_dict,
                metacyc_relations)
            output_dict = metacyc_dict
        return output_dict

    def alignment_to_query(self, query_sequence_string):
        '''Takes a query sequence as a fasta string, the FAT-CAT node it matched to and
        returns a pairwise identity between the query and the subtree hmm.
        TODO: convert all the tempfiles to named pipes (to bypass IO altogether)
        '''
        # First get the tree node hmm
        TEMP_DIRECTORY = "/clusterfs/ohana/software/webserver/incoming/"
        family_obj = self.tree.family
        family_subtree_hmm_file = family_obj.subtree_hmm_file()
        subtree_hmm_file = tempfile.NamedTemporaryFile(mode='w+b', dir=TEMP_DIRECTORY)
        hmm_utils.hmmfetch(family_subtree_hmm_file, self.id, subtree_hmm_file.name)
        # do hmm align of the hmm with the query and the hmm consensus file.
        consensus_sequence = hmm_utils.hmmemit(subtree_hmm_file.name) + '\n'
        sequence_file = tempfile.NamedTemporaryFile(mode='w+b', dir=TEMP_DIRECTORY)
        sequence_file.write(consensus_sequence)
        sequence_file.write(query_sequence_string)
        sequence_file.seek(0)
        alignment = hmm_utils.hmmalign(consensus_sequence + query_sequence_string, subtree_hmm_file.name)
        sequence_file.close()
        subtree_hmm_file.close()
        # pretty align this alignment
        pp_alignment = tempfile.NamedTemporaryFile(mode='w+b', dir=TEMP_DIRECTORY)        
        pretty_alignment = hmm_utils.prettyalign(alignment)
        pp_alignment.write(pretty_alignment)
        pp_alignment.seek(0)
        # Return the pretty alignment and the pairwise identity
        seq0, seq1 = [str(record.seq) for record in AlignIO.read(pp_alignment.name, 'fasta')]
        pp_alignment.close()
        if self.is_leaf():
            pretty_alignment = '\n'.join([str('>%s' % self.sequence_header.header)] +
                                         pretty_alignment.split('\n')[1:])
        return [pwid_belvu(seq0, seq1)*100, pretty_alignment]

    def parent(self):
        parent_left_id = self.parent_left_id
        if parent_left_id:
            return TreeNode.objects.get(tree=self.tree, left_id=parent_left_id)
        else:
            return None

    def children(self):
        return TreeNode.objects.filter(tree=self.tree, parent_left_id=self.left_id)

    def ancestors(self):
        return TreeNode.objects.filter(tree=self.tree,
                                       left_id__lt=self.left_id,
                                       right_id__gt=self.right_id).order_by("left_id")

    def ancestors_starting_from_node(self, tree_node):
        return TreeNode.objects.filter(tree=self.tree,
                                       left_id__lt=self.left_id,
                                       left_id__gte=tree_node.left_id,
                                       right_id__gt=self.right_id).order_by("left_id")

    def is_sfld_supported(self, sfld_grouping="family"):
        '''Returns the sfld group if supported else returns None.
        Is SFLD supported if both child nodes of this node have >= 1 sfld protein
        and they belong to the same sfld group.
        Alternative sfld_groupings are '''
        try:
            child1, child2 = self.children()
            child1_sfld_value = child1.sfld_value(sfld_grouping=sfld_grouping)
            if ((not child1_sfld_value) or
                child1_sfld_value != child2.sfld_value(sfld_grouping=sfld_grouping)):
                return False
            else:
                return child1_sfld_value
        except:
            return None

    def sfld_value(self, sfld_grouping="family"):
        '''Returns the sfld group that the members belong to or None. Same function as for orthomcl.'''
        if sfld_grouping == "family":
            sfld_grouping = sfld_grouping_value = "family_id"
        else:
            sfld_grouping_value = "%s_id" % sfld_grouping
        sfld_support = self.get_sfld_families().distinct(sfld_grouping).values(sfld_grouping_value)
        if len(sfld_support) != 1:
            return None
        else:
            return sfld_support[0][sfld_grouping_value]

    def sfld_supported_top_ancestor(self, start_node=None, sfld_grouping="family"):
        '''Return the current node or topmost node object along the path that is sfld supported. If start_node
        is provided, then the ancestors are chosen such that they start from that node. Else, all the ancestors
        starting from the root node are tested'''
        if self.is_sfld_supported(sfld_grouping=sfld_grouping):
            return self
        if start_node:
            ancestors = self.ancestors_starting_from_node(start_node)
        else:
            ancestors = self.ancestors()
        for ancestor in ancestors:
            if ancestor.is_sfld_supported(sfld_grouping=sfld_grouping):
                return ancestor

    def is_orthomcl_supported(self):
        '''Returns the orthomcl group if supported else returns None.
        Is Orthomcl supported if both child nodes of this node have >= 1 orthomcl protein
        and they belong to the same orthomcl group.
        try:
            child1, child2 = self.children()
            child1_orthomcl_value = child1.orthomcl_value()
            if ((not child1_orthomcl_value) or
                child1_orthomcl_value != child2.orthomcl_value()):
                return False
            else:
                return child1_orthomcl_value
        except:
            return None'''
        try:
            return self.orthology_group.get().orthomcl
        except:
            return False
        
    def orthomcl_value(self):
        '''Returns the orthomcl group that the members belong to or None.'''
        informative_members = self.get_informative_members()
        # Get orthomcl groups for all members, flatten list and return None for no_group cases.
        orthomcl_support = [leaf.sequence_header.uniprot.get_orthomcl_group()
                            for leaf in informative_members]
        orthomcl_support = [item for sublist in orthomcl_support
                            for item in sublist]
        if 'no_group' in orthomcl_support: return None
        # remove duplicates and return none if either an empty list or if there are two groups
        orthomcl_support = set(orthomcl_support)
        if len(orthomcl_support) != 1:
            return None
        else:
            return list(orthomcl_support)[0]

    def orthomcl_supported_top_ancestor(self, start_node=None):
        '''Return the current node or topmost node object along the path that is orthomcl supported. If start_node
        is provided, then the ancestors are chosen such that they start from that node. Else, all the ancestors
        starting from the root node are tested'''
        if start_node:
            ancestors = self.ancestors_starting_from_node(start_node)
        else:
            ancestors = self.ancestors()
        for ancestor in ancestors:
            if ancestor.is_orthomcl_supported():
                return ancestor
        if self.is_orthomcl_supported():
            return self
        else:
            return None

    def is_oma_supported(self):
        '''Returns the oma group if supported else returns None.
        Is OMA supported if both child nodes of this node have >= 1 oma protein
        and they belong to the same oma group. Same function as for orthomcl
        try:
            child1, child2 = self.children()
            child1_oma_value = child1.oma_value()
            if ((not child1_oma_value) or
                child1_oma_value != child2.oma_value()):
                return False
            else:
                return child1_oma_value
        except ValueError:
            return None'''
        try:
            return self.orthology_group.get().oma
        except:
            return False       

    def oma_value(self):
        '''Returns the oma group that the members belong to or None. Same function as for orthomcl.'''
        informative_members = self.get_informative_members()
        # Get oma groups for all members, flatten list and return None for no_group cases.
        oma_support = [leaf.sequence_header.uniprot.get_oma_group()
                            for leaf in informative_members]
        oma_support = [item for sublist in oma_support
                            for item in sublist]
        if 'no_group' in oma_support: return None
        # remove duplicates and return none if either an empty list or if there are two groups
        oma_support = set(oma_support)
        if len(oma_support) != 1:
            return None
        else:
            return list(oma_support)[0]

    def oma_supported_top_ancestor(self, start_node=None):
        '''Return the current or topmost node object along the path that is oma supported. If start_node
        is provided, then the ancestors are chosen such that they start from that node. Else, all the ancestors
        starting from the root node are tested. Same function as for orthomcl'''
        if start_node:
            ancestors = self.ancestors_starting_from_node(start_node)
        else:
            ancestors = self.ancestors()
        for ancestor in ancestors:
            if ancestor.is_oma_supported():
                return ancestor
        if self.is_oma_supported():
            return self
        else:
            return None 

    def kerf_supported_top_ancestor(self, start_node=None, pwid_type="belvu", threshold=70):
        '''Return the current or topmost node object along the path that is kerf supported. If start_node
        is provided, then the ancestors are chosen such that they start from that node. Else, all the ancestors
        starting from the root node are tested. Same function as for oma and orthomcl'''
        if start_node:
            ancestors = self.ancestors_starting_from_node(start_node)
        else:
            ancestors = self.ancestors()
        for ancestor in ancestors:
            if ancestor.is_kerf_supported(pwid_type=pwid_type, threshold=threshold):
                return ancestor
        if self.is_kerf_supported(pwid_type=pwid_type, threshold=threshold):
            return self
        return None
        
    def phog_supported_top_ancestor(self, start_node=None, threshold=0):
        '''Return the current or topmost node object along the path that is phog supported. If start_node
        is provided, then the ancestors are chosen such that they start from that node. Else, all the ancestors
        starting from the root node are tested. Same function as for oma and orthomcl'''
        if start_node:
            ancestors = self.ancestors_starting_from_node(start_node)
        else:
            ancestors = self.ancestors()
        for ancestor in ancestors:
            if ancestor.is_phog(threshold=threshold):
                return ancestor
        if self.is_phog(threshold=threshold):
            return self
        else:
            return None

    def orthology_supported_ancestors(self, kerf_threshold=70):
        '''Returns the results from oma, orthomcl, kerf70 supported ancestors.'''
        try:
            return self._orthology_supported_ancestors
        except:
            pass
        self._orthology_supported_ancestors = {'oma_supported_ancestor' : self.oma_supported_top_ancestor(),
                'orthomcl_supported_ancestor' : self.orthomcl_supported_top_ancestor(),
                'kerf_supported_ancestor' : self.kerf_supported_top_ancestor(pwid_type="belvu", threshold=kerf_threshold),
                'phog_supported_ancestor' : self.get_containing_phog()}
                #'phog_supported_ancestor' : self.phog_supported_top_ancestor()}
        return self._orthology_supported_ancestors

    ''' Public function to create a phyloxml representation of this node and its descendants '''
    def phyloxml(self, **kwargs):
        try:
            phyloxml_E
        except:
            phyloxml_E = objectify.ElementMaker(annotate=False,
                namespace=u'http://www.phyloxml.org',
                nsmap={
                  u'xs': u'http://www.w3.org/2001/XMLSchema',
                  u'phy': "http://www.phyloxml.org"
                })
        return etree.tostring(phyloxml_E.phyloxml(
            phyloxml_E.phylogeny(
                self._phyloxml(**kwargs),
                rooted="true")), pretty_print=True)

    ''' Private function to create the phyloxml - doesnt make the outer phyloxml tag '''
    # Jonathan, you are a crap programmer.  Rewrite this!!!!
    def _phyloxml(self, **kwargs):
        # There is probably a better way to do this - trying not to create this many times
        # Maybe just use the BioPython PhyloXML Stuff?  It is a little janky, but might
        # Be better than reinventing the wheel
        try:
            phyloxml_E
        except:
            phyloxml_E = objectify.ElementMaker(annotate=False,
                namespace=u'http://www.phyloxml.org',
                nsmap={
                  u'xs': u'http://www.w3.org/2001/XMLSchema',
                  u'phy': "http://www.phyloxml.org"
                })

        stop = False
        if self.is_leaf():
            stop = True
        
        # attributes are the kwargs for the clade, arguments are the sub objects (clades, attributes, etc)
        attributes = {}
        annotations = []
        arguments = []
        # We use negative branch lengths as "Magic" values.  This is stupid. Turning negatives to 0 makes it so viewers work with our data
        attributes["branch_length"] = str(max(self.distance_to_parent, 0))
        if ((kwargs.has_key("annotations")) and ("percent_id" in kwargs["annotations"])):
            arguments.append(phyloxml_E.property("%.1f" % (self.minimum_pairwise_identity_belvu),
                ref="minimum_percent_id", applies_to="node", datatype = "xsd:float", unit="percent"))
        
        arguments.append(phyloxml_E.node_id(self.id))

        if not stop:
            # we are an internal node
            if (kwargs.has_key("annotations")):
                # append the mrca to this internal node if we can.  this is done in a taxonomy element
                if (("mrca" in kwargs["annotations"])):
                    try:
                        mrca = self.mrca.get().mrca
                        tax_args = []
                        tax_args.append(phyloxml_E.id(mrca.id))
                        tax_args.append(phyloxml_E.scientific_name(mrca.scientific_name))
                        if (mrca.common_name):
                            tax_args.append(phyloxml_E.common_name(mrca.common_name))
                        arguments.append(phyloxml_E.taxonomy(*tax_args))
                    except:
                        pass
                # append the consensus uniprot description for this node to the name element
                if (("consensus_descriptions" in kwargs["annotations"])):
                    try:
                        desc = self.consensus_description.get().consensus_uniprot_description
                        arguments.append(phyloxml_E.name(desc))
                    except:
                        pass
            
            for child in self.children():
                arguments.append(child._phyloxml(**kwargs)) 
                
        else:
            # we are a leaf
            sequence_arguments = []
            if (self.sequence_header.uniprot):
                # we have a uniprot entry
                arguments.append(phyloxml_E.name(self.sequence_header.uniprot.description))
                sequence_arguments.append(phyloxml_E.symbol(self.sequence_header.uniprot.uniprot_identifier))
                sequence_arguments.append(phyloxml_E.accession(self.sequence_header.uniprot.accession))
                sequence_arguments.append(phyloxml_E.name(self.sequence_header.uniprot.description))
                # annotate the phyloxml if we need to
                if (kwargs.has_key("annotations")):
                    if "GO" in kwargs["annotations"]:
                        """sequence_arguments.append(phyloxml_E.annotation(
                                phyloxml_E.desc(annotation.ontology_term.name),
                                ref=annotation.ontology_term.accession,
                                evidence=annotation.evidence.code,
                                # type="TODO"
                            ))"""
                        if Annotation.objects.filter(
                                uniprot__sequenceheader__treenode__tree = self.tree,
                                uniprot__sequenceheader__treenode__left_id__gte = self.left_id,
                                uniprot__sequenceheader__treenode__right_id__lte = self.right_id,
                                ontology_term__ontology__name = 'go').exists():
                            sequence_arguments.append(phyloxml_E.annotation(" ", type='go'))
                    #if "EC" in kwargs["annotations"]:
                    #    if self.get_ecs(experimental=False):
                    #        sequence_arguments.append(phyloxml_E.annotation(" ", type='ec'))
                    if "literature" in kwargs["annotations"]:
                        if UniProtLiterature.objects.filter(uniprot = self.sequence_header.uniprot).exists():
                            sequence_arguments.append(phyloxml_E.annotation(" ", type='literature'))
                    if "structure" in kwargs["annotations"]:
                        if UniProtPDB_Chain.objects.filter(uniprot = self.sequence_header.uniprot).exists():
                            sequence_arguments.append(phyloxml_E.annotation(" ", type='structure'))
                    if "swissProt" in kwargs["annotations"]:
                        if self.sequence_header.uniprot.in_swissprot_f:
                            sequence_arguments.append(phyloxml_E.annotation(" ", type="swissProt"))
                    if "taxonomy" in kwargs["annotations"]:
                        try:
                            taxon = self.sequence_header.uniprot.taxon
                            tax_args = []
                            tax_args.append(phyloxml_E.id(taxon.id))
                            tax_args.append(phyloxml_E.scientific_name(taxon.scientific_name))
                            if (taxon.common_name):
                                tax_args.append(phyloxml_E.common_name(taxon.common_name))
                            arguments.append(phyloxml_E.taxonomy(*tax_args))
                        except:
                            pass
            else:
                # we have no uniprot entry, we should parse the header to show what we can
                h = self.sequence_header.parse_header()
                arguments.append(phyloxml_E.name(h['description']))
                sequence_arguments.append(phyloxml_E.name(h['description']))
                if h['accession']:
                    sequence_arguments.append(phyloxml_E.accession(h['accession']))
                if h['identifier']:
                    sequence_arguments.append(phyloxml_E.symbol(h['identifier']))
                if (kwargs.has_key("annotations")):
                    if "swissProt" in kwargs["annotations"]:
                        if (h['db'] == 'sp'):
                            sequence_arguments.append(phyloxml_E.annotations(type="swissProt"))
                    if "taxonomy" in kwargs["annotations"]:
                        if (h['taxon_id']):
                            tax_args = []
                            tax_args.append(phyloxml_E.id(h['taxon_id']))
                            tax_args.append(phyloxml_E.scientific_name(h['scientific_name']))
                            if (h['common_name']):
                                tax_args.append(phyloxml_E.common_name(h['common_name']))
                            arguments.append(phyloxml_E.taxonomy(*tax_args))
                        elif (h['species']):
                            arguments.append(phyloxml_E.taxonomy(phyloxml_E.scientific_name(h['species'])))
                        else:
                            # No taxon
                            pass

            arguments.append(phyloxml_E.sequence(*sequence_arguments))

        return phyloxml_E.clade(
            *arguments,
            **attributes
        )

    ''' We need to generate 3 different types of resultsets, preferably with one query each
        The input for all of them is a treenode_id

        The outputs are:

        1: all nodes between the input and the leaves
        2: all nodes between the input and the informative nodes
        3: all nodes between the input and the leaves, and the summary from the root down.

    '''
    def get_subtree(self, type='full'):
        nodes = TreeNode.objects.filter(
            tree = self.tree,
            left_id__gte = self.left_id,
            right_id__lte = self.right_id
        ).order_by('left_id')

        root = nodes[0]
        root.children = []
        node_dict = {root.left_id: root}
        for node in nodes[1:]:
            node_dict[node.left_id] = node
            node.children = []
            node_dict[node.parent_left_id].children.append(node)

        def prune(node, pruning_criteria):
            pass

        if (type == 'full'):
            return root
        elif (type == 'summary'):
            return root
        elif (type == 'compressed'):
            return root
        else:
            raise "unknown type of subtree"

    def newick(self, type='full'):
        tree = self.get_subtree(type=type)
        def _n(node):
            # TODO make our DB less stupid
            if node.distance_to_parent == -1.0:
                distance_to_parent = 0
            else:
                distance_to_parent = node.distance_to_parent
        
            if (len(node.children)):
                return "(%s):%s" % (",".join([ _n(n) for n in node.children ]), distance_to_parent)
            else:
                return "%s:%s" % (node.name, distance_to_parent)
        return _n(tree)
  
    def phyloscope_json(self, top_scoring_node = None):
        ''' Returns the json object necessary for phyloscope - hopefully. '''
        leaves=self.get_included_leaves()
        if top_scoring_node:
            tree = self.newick_sequence_header_ids(top_scoring_node = top_scoring_node)
        else:
            tree = self.newick_sequence_header_ids()
        tree += ';'
        sequence_header_ids = [leaf.sequence_header.id for leaf in leaves]
        family = self.tree.family 
        sequence_headers = SequenceHeader.objects.filter(id__in =
                                                            sequence_header_ids)
        annotation_info = {}
        sequence_headers_of_uniprot = {}
        for sequence_header in sequence_headers:
            annotation_info[sequence_header] = {}
            annotation_info[sequence_header]['header'] = sequence_header.header
            annotation_info[sequence_header]['in_swissprot_f'] = 0
            taxon = None
            if sequence_header.uniprot is None:
                if family.seed_sequence_header \
                and sequence_header.id == family.seed_sequence_header.id:
                    annotation_info[sequence_header]['uniprot_de'] \
                    = string.upper(sequence_header.header)
                else:
                    annotation_info[sequence_header]['uniprot_de'] = sequence_header.header
                    annotation_info[sequence_header]['uniprot_id'] = 'N/A'
                    annotation_info[sequence_header]['uniprot_accession'] = 'N/A'
                    taxon = sequence_header.taxon
            else:
                if sequence_header.uniprot not in sequence_headers_of_uniprot:
                    sequence_headers_of_uniprot[sequence_header.uniprot] = set()
                sequence_headers_of_uniprot[sequence_header.uniprot].add(sequence_header)
                if family.seed_sequence_header \
                  and sequence_header.id == family.seed_sequence_header.id:
                    annotation_info[sequence_header]['uniprot_de'] \
                    = string.upper("%s %s" % (sequence_header.header,
                                              sequence_header.uniprot.description))
                else:
                    annotation_info[sequence_header]['uniprot_de'] \
                        = sequence_header.uniprot.description
                annotation_info[sequence_header]['uniprot_id'] \
                  = sequence_header.uniprot.uniprot_identifier
                annotation_info[sequence_header]['uniprot_accession'] \
                  = sequence_header.uniprot.accession
                if sequence_header.uniprot.in_swissprot_f:
                    annotation_info[sequence_header]['in_swissprot_f'] = 1
                taxon = sequence_header.uniprot.taxon
            # Remove backslashes which break JavaScript htmlEntityDecode
            annotation_info[sequence_header]['uniprot_de'] \
                = annotation_info[sequence_header]['uniprot_de'].replace('\\', '')
            annotation_info[sequence_header]['scientific_name'] = 'N/A'
            annotation_info[sequence_header]['common_name'] = 'N/A'
            annotation_info[sequence_header]['ncbi_taxid'] = 'N/A'
            if taxon is not None:
                annotation_info[sequence_header]['ncbi_taxid'] = taxon.id
                if taxon.scientific_name:
                    annotation_info[sequence_header]['scientific_name'] \
                        = taxon.scientific_name
                if taxon.common_name:
                    annotation_info[sequence_header]['common_name'] \
                        = taxon.common_name

        uniprot_ecs = UniProtEC.objects.filter(uniprot__in =
                                                  sequence_headers_of_uniprot.keys())

        for uniprot_ec in uniprot_ecs:
            for sequence_header in sequence_headers_of_uniprot[uniprot_ec.uniprot]:
                if 'ec' in annotation_info[sequence_header]:
                    annotation_info[sequence_header]['ec'] \
                    = ';'.join([annotation_info[sequence_header]['ec'],
                              uniprot_ec.ec.__unicode__()])
                else:
                    annotation_info[sequence_header]['ec'] = uniprot_ec.ec.__unicode__()

        uniprot_lits = UniProtLiterature.objects.filter(uniprot__in =
                          sequence_headers_of_uniprot.keys()).exclude(
                                                      is_large_scale_f = True)
        for uniprot_lit in uniprot_lits:
            if uniprot_lit.pmid:
                for sequence_header in sequence_headers_of_uniprot[uniprot_lit.uniprot]:
                    if 'Literature' not in annotation_info[sequence_header]:
                        annotation_info[sequence_header]['Literature'] = set()
                    annotation_info[sequence_header]['Literature'].add(uniprot_lit.pmid)

        uniprot_pdb_chains = UniProtPDB_Chain.objects.filter(uniprot__in =
                                                  sequence_headers_of_uniprot.keys())
        for uniprot_pdb_chain in uniprot_pdb_chains:
            for sequence_header in sequence_headers_of_uniprot[uniprot_pdb_chain.uniprot]:
                if 'ThreeDStructure' not in annotation_info[sequence_header]:
                    annotation_info[sequence_header]['ThreeDStructure'] = set()
                annotation_info[sequence_header]['ThreeDStructure'].add('%s%s' %
                                                    (uniprot_pdb_chain.pdb_chain.pdb.id,
                                                    uniprot_pdb_chain.pdb_chain.chain_id))

        uniprot_gos = UniProtGO.objects.filter(uniprot__in =
                                                  sequence_headers_of_uniprot.keys())

                 # PhyloScope itself knows about experimental vs. nonexperimental evidence, 
          # so just give it everything
        for uniprot_go in uniprot_gos:
            for sequence_header in sequence_headers_of_uniprot[uniprot_go.uniprot]:
                if uniprot_go.go_term.term_type not in annotation_info[sequence_header]:
                    annotation_info[sequence_header][uniprot_go.go_term.term_type] = set()
                annotation_info[sequence_header][uniprot_go.go_term.term_type].add(
                  (uniprot_go.go_term.name, uniprot_go.go_evidence.evidence))

        ret = {}
        for sequence_header in sequence_headers:
            ret[str(sequence_header.id)] = annotation_info[sequence_header]
            if 'Literature' in annotation_info[sequence_header]:
                ret[str(sequence_header.id)]['Literature'] \
                  = list(annotation_info[sequence_header]['Literature'])
            if 'ThreeDStructure' in annotation_info[sequence_header]:
                ret[str(sequence_header.id)]['ThreeDStructure'] \
                  = list(annotation_info[sequence_header]['ThreeDStructure'])
            if 'biological_process' in annotation_info[sequence_header]:
                ret[str(sequence_header.id)]['biological_process'] \
                  = list(annotation_info[sequence_header]['biological_process'])
            if 'molecular_function' in annotation_info[sequence_header]:
                ret[str(sequence_header.id)]['molecular_function'] \
                  = list(annotation_info[sequence_header]['molecular_function'])
            if 'cellular_component' in annotation_info[sequence_header]:
                ret[str(sequence_header.id)]['cellular_component'] \
                  = list(annotation_info[sequence_header]['cellular_component'])
        ret['__tree__'] = tree
        return json.dumps(ret)                              
    
    def newick_sequence_header_ids(self, type='full', top_scoring_node = None):
        tree = self.get_subtree(type=type)
        def _n(node):
            # TODO make our DB less stupid
            if node.distance_to_parent == -1.0:
                distance_to_parent = 0
            else:
                distance_to_parent = node.distance_to_parent
        
            if (len(node.children)):
                if top_scoring_node and (top_scoring_node == node.id):
                    return "(%s)TSN:%s" % (",".join([ _n(n) for n in node.children ]), distance_to_parent)
                else:
                    return "(%s):%s" % (",".join([ _n(n) for n in node.children ]), distance_to_parent)
            else:
                if top_scoring_node and (top_scoring_node == node.id):
                    return "%d:%s" % (node.sequence_header.id, distance_to_parent)
                else:
                    return "%d:%s" % (node.sequence_header.id, distance_to_parent)
        return _n(tree)

    @property
    def name(self):
        try:
            return self.sequence_header.uniprot.accession
        except:
            try:
                return self.sequence_header.header.split(' ')[0]
            except:
                return name.__unicode__()

    def get_data_from_cache(self):
        return FamilyCache.objects.get(family_id = self.__class__.cached_treenodes_and_families_map[self.id]).cached_data

    def data_for_caching(self):
        go_data = self.go_data_for_caching()
        ecs = self.get_ecs_for_caching()
        description = self.get_description_for_caching()
        taxonomic_distributrion = str(self.get_taxon())
        return self.__class__.delimiters['field_delimiter'].join([go_data, ecs,
            description, taxonomic_distributrion])

    # RSD: this is not necessary...
    def get_phog(self, threshold=0.296874):
        return PHOG(self, threshold)

    def _get_max_left_id(self):
        global max_left_id
        max_left_id = UniProtTaxonomy.objects.aggregate(m=Max('left_id'))['m']
        return max_left_id

    def _get_left_right_ids(self, key):
        global left_right_ids
        if left_right_ids.get(key, None) is None:
            n = UniProtTaxonomy.objects.filter(scientific_name=key).get()
            left_right_ids[key] = (n.left_id, n.right_id)
        return left_right_ids[key]

    id = models.AutoField(primary_key=True)
    tree = models.ForeignKey(Tree, related_name='treenodes')
    left_id = models.IntegerField()
    right_id = models.IntegerField()
    parent_left_id = models.IntegerField()
    level = models.IntegerField(null=True)
    superorthologous_node = models.ForeignKey('self', null=True,
        related_name='superorthologs')
    superorthologous_node_left_id = models.IntegerField(null=True)
    is_phogt_tight = models.NullBooleanField()
    is_phogt_medium = models.NullBooleanField()
    is_phogt_loose = models.NullBooleanField(null=True)
    phogt_tight_node = models.ForeignKey('self',
        related_name='tight_orthologs', null=True)
    phogt_medium_node = models.ForeignKey('self',
        related_name='medium_orthologs', null=True)
    phogt_loose_node = models.ForeignKey('self',
        related_name='loose_orthologs', null=True)
    duplication_distance = models.FloatField(null=True)
    greatest_duplication_distance_of_maximal_descendant \
        = models.FloatField(null=True)
    is_superorthologous = models.NullBooleanField()
    sequence_header = models.ForeignKey(SequenceHeader,null=True)
    distance_to_parent = models.FloatField(default=-1.0, null=True)
    canonical_tree_node_name = models.ForeignKey('TreeNodeName',
        related_name="tree_cannonical_tree_node_name", null=True)
    cellular_most_recent_common_taxon = models.ForeignKey(UniProtTaxonomy,
        db_column='cellular_most_recent_common_taxon', null=True)
    viral_most_recent_common_taxon = models.ForeignKey(UniProtTaxonomy,
        db_column='viral_most_recent_common_taxon',
        related_name="treenode_uniprottaxonomy_mrc_taxon",
        null=True)
    bootstrap_support = models.FloatField(null=True)
    likelihood = models.FloatField(null=True)
    minimum_pairwise_identity = models.FloatField(null=True)
    minimum_pairwise_identity_belvu = models.FloatField(null=True)

    def is_kerf_root(self, pwid_type="belvu", threshold=70):
        ''' Returns true if the node is the root of a kerf group.  False otherwise. '''
        parent = self.parent()
        if parent:
            if pwid_type == "belvu":
                return ((self.minimum_pairwise_identity_belvu >= threshold) and (parent.minimum_pairwise_identity_belvu < threshold))
            else:                
                return ((self.minimum_pairwise_identity >= threshold) and (parent.minimum_pairwise_identity < threshold))
        else:
            # This node is the root of a tree, is it a kerf group?
            if pwid_type == "belvu":
                return ((self.minimum_pairwise_identity_belvu >= threshold))
            else:                
                return ((self.minimum_pairwise_identity >= threshold))
    
    def is_kerf_supported(self, pwid_type="belvu", threshold=70):
        '''Returns true if the node is supported by kerf with the minimum pairwise identity
        specified, else returns false.'''
        if pwid_type == "default":
            pwid = self.minimum_pairwise_identity
        else: # pwid_type == "belvu":
            pwid = self.minimum_pairwise_identity_belvu
        if pwid >= float(threshold):
            return pwid
        else:
            return False

    def get_sub_phogts(self, ortholog_type, threshold):
        return get_ortholog_type_of_threshold(threshold) != \
            OrthologTypes.SuperOrtholog and self.tree.treenodes.filter(
                left_id__gt=self.left_id,
                right_id__lt=self.right_id,
                duplication_distance__gte=threshold,
            ) or TreeNode.objects.none()

    def get_included_leaves(self, ortholog_type=OrthologTypes.SuperOrtholog,
                            threshold=0.0):
        if not self.included_leaves:
            leaves = TreeNode.objects.filter(tree__exact = self.tree,
                                      left_id__gte = self.left_id,
                                      right_id__lte = self.right_id).exclude(
                                          sequence_header__exact = 0).exclude(
                                          sequence_header__isnull = True)
            self.included_leaves = leaves
        return self.included_leaves

    def get_internal_nodes(self):
        try:
            return self.included_nodes
        except:
            nodes = TreeNode.objects.filter(tree__exact = self.tree,
                                      left_id__gte = self.left_id,
                                      right_id__lte = self.right_id,
                                      sequence_header__isnull = True)
            self.included_nodes = nodes
            return self.included_nodes

    def get_containing_phog(self, threshold=0.0):
        # The tree nodes that are ancestral to the query node will have the
        # property left_id__lte = query_node.left_id and right_id__gte =
        # query_node.right_id. The tree nodes that are maximal tree nodes
        # (i.e., potential PHOGs) will have a non-null value of
        # duplication_distance. The tree nodes that meet the threshold
        # will satisfy the duplication_distance__gte = threshold.
        # However, these will only be PHOGs when none of their descendants also
        # meet the threshold, i.e., when
        # greatest_duplication_distance_of_maximal_descendant__lt = threshold.
        # When we order_by('right_id'), we are going up the tree from the
        # leaf toward the root.
        containing_phogs = TreeNode.objects.filter(tree__exact = self.tree,
                  left_id__lte = self.left_id,
                  right_id__gte = self.right_id,
                  duplication_distance__gte = threshold,
                  greatest_duplication_distance_of_maximal_descendant__lt
                  = threshold).order_by('right_id')
        if containing_phogs:
            return containing_phogs[0]
        else:
            return None

    def get_containing_kerf(self, threshold=70):
        kerfs = TreeNode.objects.filter(tree__exact = self.tree,
                  left_id__lte = self.left_id,
                  right_id__gte = self.right_id,
                  minimum_pairwise_identity_belvu__gte = threshold
                 ).order_by('left_id')
        if kerfs:
            return kerfs[0]
        else:
            return None

    def get_informative_node(self, fatcat_family):
        if ((self.minimum_pairwise_identity_belvu >= fatcat_family.fatcat_job.criteria_get_best_nodes_subtree_min_pid) and 
            (self.alignment_to_query(fatcat_family.fatcat_job.fasta)[0] >= fatcat_family.fatcat_job.criteria_get_best_nodes_query_to_subtree_hmm_pid) and
            (fatcat_family.passes_first_stage_coverage_conditions) and (fatcat_family.passes_second_stage_coverage_conditions)):
            orthology_methods = fatcat_family.fatcat_job.criteria_enclosing_clade_orthology_methods.strip().split()
            orthology_methods = [x + '_supported_ancestor' for x in orthology_methods]
            orthology_methods_required = fatcat_family.fatcat_job.criteria_enclosing_clade_orthology_methods_required
            kt = fatcat_family.fatcat_job.criteria_enclosing_clade_kerf_threshold
            o_items = [x for x in self.orthology_supported_ancestors(kerf_threshold=kt).items() if x[1] is not None and x[0] in orthology_methods]
            if o_items:
                o_items = sorted(o_items, key=lambda t: t[1].left_id)
            else:
                return None
            # this can be rewritten much more simply...and it should be.
            matched = 1
            top_ec = o_items[0]
            o_items.remove(o_items[0])
            if matched == orthology_methods_required:
                return top_ec[1]
            for (system, root) in o_items:
                if matched == orthology_methods_required:
                    return top_ec[1]
                if root == top_ec[1]:
                    matched += 1
                else:
                    top_ec = root
                    matched = 1
            return None
            #if ((o_items[0][1]) and ((o_items[0][1] == o_items[1][1]) or (o_items[0][1] == o_items[2][1]) or (o_items[0][1] == o_items[3][1]))):
            #    return o_items[0][1]
            #elif ((o_items[1][1]) and ((o_items[1][1] == o_items[2][1]) or (o_items[1][1] == o_items[3][1]))):
            #    return o_items[1][1]
            #elif ((o_items[2][1]) and ((o_items[2][1] == o_items[3][1]))):
            #    return o_items[2][1]
            #else:
            #    return None
            #largest_orthologous_clade_left_id = self.left_id + 1
            #largest_orthologous_clade = None 
            #for (system, root) in self.orthology_supported_ancestors().items():
            #    if (root and (root.left_id < largest_orthologous_clade_left_id)):
            #        largest_orthologous_clade = root 
            #        largest_orthologous_clade_left_id = root.left_id     
            #return largest_orthologous_clade        
        else:
            return None

    def get_num_included_leaves(self,
                                ortholog_type = OrthologTypes.SuperOrtholog,
                                threshold = 0.0):
        return len(self.get_included_leaves(
                                      ortholog_type = ortholog_type,
                                      threshold = threshold))

    def get_num_contained_species(self,
                                ortholog_type = OrthologTypes.SuperOrtholog,
                                threshold = 0.0):
        if self.num_contained_taxa == 0:
            leaves = self.get_included_leaves(ortholog_type, threshold)
            contained_taxa = set()
            for leaf in leaves:
                taxon = leaf.get_taxon(ortholog_type, threshold)
                if taxon:
                    contained_taxa.add(taxon.id)
            self.num_contained_taxa = len(contained_taxa)
        return self.num_contained_taxa

    def get_contained_leaves_from_taxon(self, taxon,
        ortholog_type = OrthologTypes.SuperOrtholog,
        threshold = 0.0):
        leaves = self.get_included_leaves(ortholog_type, threshold)
        leaf_ids = set()
        for leaf in leaves:
            leaf_taxon = leaf.get_taxon(ortholog_type, threshold)
            if leaf_taxon:
                if leaf_taxon.left_id >= taxon.left_id \
                    and leaf_taxon.right_id <= taxon.right_id:
                    leaf_ids.add(leaf.id)
        return leaves.filter(id__in = leaf_ids)

    def get_num_nonredundant_sequences(self,
        ortholog_type = OrthologTypes.SuperOrtholog,
        threshold = 0.0):
        if self.num_nonredundant_sequences == 0:
            leaves = self.get_included_leaves(ortholog_type, threshold)
            identifiers = set()
            for leaf in leaves:
                identifiers.add(leaf.sequence_header.identifier())
            self.num_nonredundant_sequences = len(identifiers)
        return self.num_nonredundant_sequences

    def __unicode__(self):
        return u'%s' % self.id

    def get_hyper_neighbors(self, ortholog_type = OrthologTypes.SuperOrtholog,
        threshold = 0.0):
        partner_phog_ids = set()
        partner_phog_ids2 = set()
        phog_edges = {}
        leaves = self.get_included_leaves(ortholog_type, threshold)
        for leaf in leaves:
            if leaf.sequence_header.uniprot:
                partner_uniprots =\
                    leaf.sequence_header.uniprot.get_interacting_partners()
                if len(partner_uniprots) == 0:
                    continue

                uniprot_accessions =\
                    [uniprot.uniprot_accession for uniprot in partner_uniprots]
                if ortholog_type == OrthologTypes.SuperOrtholog:
                    phogQ = Q(superorthologs__sequence_header__uniprot__uniprot_accession__in=uniprot_accessions)
                elif ortholog_type == OrthologTypes.PHOG_T_Tight:
                    phogQ = Q(tight_orthologs__sequence_header__uniprot__uniprot_accession__in=uniprot_accessions)
                elif ortholog_type == OrthologTypes.PHOG_T_Medium:
                    phogQ = Q(medium_orthologs__sequence_header__uniprot__uniprot_accession__in=uniprot_accessions)
                elif ortholog_type == OrthologTypes.PHOG_T_Loose:
                    phogQ = Q(loose_orthologs__sequence_header__uniprot__uniprot_accession__in=uniprot_accessions)
                elif ortholog_type == OrthologTypes.PHOG_T_Custom:
                    query_nodes = TreeNode.objects.filter(sequence_header__uniprot__uniprot_accession__in=uniprot_accessions)
                if ortholog_type == OrthologTypes.PHOG_T_Custom:
                    for query_node in query_nodes:
                        partner_phog_ids.add(TreeNode.objects.filter(
                                  tree__exact = query_node.tree,
                                  left_id__lte = query_node.left_id,
                                  right_id__gte = query_node.right_id,
                                  duplication_distance__gte = threshold).order_by('right_id')[0].id)
                else:
                    partner_phogs = TreeNode.objects.filter(phogQ).distinct()
                    for partner_phog in partner_phogs:
                        partner_phog_ids.add(partner_phog.id)
                ## Do the db query for each protein so that the edges
                ## connecting PHOG neighbors can be returned. (phog_edges)
                for uo in partner_uniprots:
                    ua = uo.uniprot_accession
                    if ortholog_type == OrthologTypes.SuperOrtholog:
                        phogQ = Q(superorthologs__sequence_header__uniprot__uniprot_accession=ua)
                    elif ortholog_type == OrthologTypes.PHOG_T_Tight:
                        phogQ = Q(tight_orthologs__sequence_header__uniprot__uniprot_accession=ua)
                    elif ortholog_type == OrthologTypes.PHOG_T_Medium:
                        phogQ = Q(medium_orthologs__sequence_header__uniprot__uniprot_accession=ua)
                    elif ortholog_type == OrthologTypes.PHOG_T_Loose:
                        phogQ = Q(loose_orthologs__sequence_header__uniprot__uniprot_accession=ua)
                    elif ortholog_type == OrthologTypes.PHOG_T_Custom:
                        query_nodes = TreeNode.objects.filter(sequence_header__uniprot__uniprot_accession=ua)
                    if ortholog_type == OrthologTypes.PHOG_T_Custom:
                        phog_list = []
                        for query_node in query_nodes:
                            phogs = TreeNode.objects.filter(
                                tree__exact = query_node.tree,
                                left_id__lte = query_node.left_id,
                                right_id__gte = query_node.right_id,
                                duplication_distance__gte = threshold)
                            for p in phogs:
                                phog_list.append(p)
                    else:
                        phog_list = TreeNode.objects.filter(phogQ).distinct()
                    for p in phog_list:
                        partner_phog_ids2.add(p.id)
                        key = "%s:%s" % (leaf.sequence_header.uniprot.uniprot_accession, ua)
                        if p.get_accession() not in phog_edges:
                            phog_edges[p.get_accession()] = {}
                            phog_edges[p.get_accession()][key] = (leaf.sequence_header.uniprot, uo)

        hyper_neighbors = {}
        hyper_neighbors['PPI'] = TreeNode.objects.filter(id__in = partner_phog_ids)

        return (hyper_neighbors, phog_edges)

    def get_ecs(self, experimental=True,
                ortholog_type = OrthologTypes.SuperOrtholog,
                threshold = 0.0):
        if not self.have_collected_ecs:
            experimental_ecs = set()
            nonexperimental_ecs = set()
            leaves = self.get_included_leaves(ortholog_type, threshold)
            annotations = UniProtEC.objects.filter(uniprot__in =
                            [leaf.sequence_header.uniprot
                            for leaf in leaves if leaf.sequence_header.uniprot])
            for annotation in annotations:
                if annotation.is_in_brenda_f:
                    experimental_ecs.add(annotation.ec)
                    if annotation.ec in nonexperimental_ecs:
                        nonexperimental_ecs.remove(annotation.ec)
                else:
                    if annotation.ec not in experimental_ecs:
                        nonexperimental_ecs.add(annotation.ec)
            self.experimental_ecs = experimental_ecs
            self.nonexperimental_ecs = nonexperimental_ecs
            self.have_collected_ecs = True
        if experimental:
            return self.experimental_ecs
        else:
            return self.nonexperimental_ecs

    def get_ecs_for_caching(self):
        all_ecs = [self.get_ecs(),  self.nonexperimental_ecs]
        ec_list = []
        for i in range(len(all_ecs)):
            for ec in all_ecs[i]:
                if not ec.description:
                    ec.description = ""
                else:
                    if ec.description[-1] == ".":
                        ec.description = ec.description[:-1]
                # experimental_ecs will be flagged with a 1; non-experimental ecs with a 0
                ec_list.append(self.__class__.delimiters['component_delimiter'].join([ec.description,
                    str(ec), str((i + 1) % 2)]))

        return self.__class__.delimiters['record_delimiter'].join(ec_list)

    def get_sfld_families(self):
        return SFLDFamily.objects.filter(
            efds__sfldtopfacts__uniprot__sequenceheader__treenode__tree = self.tree,
            efds__sfldtopfacts__uniprot__sequenceheader__treenode__left_id__gte = self.left_id,
            efds__sfldtopfacts__uniprot__sequenceheader__treenode__right_id__lte = self.right_id,
        ) 

    def get_pmids(self, include_large_scale = False,
                  ortholog_type = OrthologTypes.SuperOrtholog,
                  threshold = 0.0):
        if not self.have_collected_pmids:
            leaves = self.get_included_leaves(ortholog_type,
                    threshold).select_related('sequence_header__uniprot__id')
            uniprot_ids = set([leaf.sequence_header.uniprot.id for leaf in
                leaves if leaf.sequence_header.uniprot])
            if include_large_scale:
                refs = UniProtLiterature.objects.filter(uniprot__id__in = uniprot_ids,
                                                        pmid__isnull = False)
            else:
                refs = UniProtLiterature.objects.filter(uniprot__id__in = uniprot_ids,
                                                      pmid__isnull = False).exclude(
                                                      is_large_scale_f = False)
            self.pmids = set([ref.pmid for ref in refs])
            self.have_collected_pmids = True
        return self.pmids

    def get_pdb_chains(self):
        if not self.have_collected_pdb_chains:
            leaves = self.get_included_leaves().select_related(
                                                  'sequence_header__uniprot__id')
            uniprot_ids = set([leaf.sequence_header.uniprot.id for leaf in leaves
                              if leaf.sequence_header.uniprot])
            uniprot_pdb_chains = UniProtPDB_Chain.objects.filter(uniprot__id__in
                                                                = uniprot_ids)
            self.pdb_chains = set([uniprot_pdb_chain.pdb_chain for
                                    uniprot_pdb_chain in uniprot_pdb_chains])
        return self.pdb_chains

    def biocyc_data(self):
        '''Adding biocyc data collected from individual members.'''
        pass

    def go_data(self):
        go_data = {'biological_process' : {}, 'molecular_function' : {}, 'cellular_component' : {} }

        for k in go_data.keys():
            annotations = summarize(self.annotations("go_" + k))
            for a in [ a.annotation for a in annotations ]:
                go_data[k][a['description']] = a
                go_data[k][a['description']]["acc"] = a['accession']
        return go_data

    def go_data_sorted_by_priority(self):
        return self.go_data()

    def get_taxon(self, ortholog_type = OrthologTypes.SuperOrtholog,
      threshold = 0.0):
        if ortholog_type not in self.taxon:
            self.taxon[ortholog_type] = {}
        if threshold not in self.taxon[ortholog_type]:
            try:
                if self.sequence_header is None:
                    raise SequenceHeader.DoesNotExist
                if self.sequence_header.uniprot is None:
                    raise UniProt.DoesNotExist
                if self.sequence_header.uniprot.taxon:
                    self.taxon[ortholog_type][threshold] \
                       = self.sequence_header.uniprot.taxon
                elif self.sequence_header.taxon:
                    self.taxon[ortholog_type][threshold] = self.sequence_header.taxon
                # otherwise the taxon will have to remain None
            except UniProtTaxonomy.DoesNotExist:
            # the taxon will remain None
                self.taxon[ortholog_type][threshold] = None
            except UniProt.DoesNotExist:
                try:
                    self.taxon[ortholog_type][threshold] = self.sequence_header.taxon
                except UniProtTaxonomy.DoesNotExist:
                    self.taxon[ortholog_type][threshold] = None
            except SequenceHeader.DoesNotExist:
                leaves = self.get_included_leaves(ortholog_type, threshold)
                max_all_left = self._get_max_left_id()
                min_cellular_left = max_all_left
                min_bacterial_left = max_all_left
                min_archaeal_left = max_all_left
                min_eukaryotic_left = max_all_left
                min_noncellular_left = max_all_left
                max_cellular_right = 0
                max_bacterial_right = 0
                max_archaeal_right = 0
                max_eukaryotic_right = 0
                max_noncellular_right = 0
                cellular_left, cellular_right \
                    = self._get_left_right_ids('cellular organisms')
                bacterial_left, bacterial_right \
                    = self._get_left_right_ids('Bacteria')
                archaeal_left, archaeal_right \
                    = self._get_left_right_ids('Archaea')
                eukaryotic_left, eukaryotic_right \
                    = self._get_left_right_ids('Eukaryota')
                for leaf in leaves:
                    leaf_taxon = leaf.get_taxon(ortholog_type, threshold)
                    if leaf_taxon:
                        left = leaf_taxon.left_id
                        right = leaf_taxon.right_id
                        if left > cellular_left and right < cellular_right:
                # Cellular organism
                            if left < min_cellular_left:
                                min_cellular_left = left
                            if right > max_cellular_right:
                                max_cellular_right = right
                            # Bacteria?
                            if left > bacterial_left and right < bacterial_right:
                                if left < min_bacterial_left:
                                    min_bacterial_left = left
                                if right > max_bacterial_right:
                                    max_bacterial_right = right
                            # Archaea?
                            elif left > archaeal_left and right < archaeal_right:
                                if left < min_archaeal_left:
                                    min_archaeal_left = left
                                if right > max_archaeal_right:
                                    max_archaeal_right = right
                            # Eukaryota?
                            elif left > eukaryotic_left and right < eukaryotic_right:
                                if left < min_eukaryotic_left:
                                    min_eukaryotic_left = left
                                if right > max_eukaryotic_right:
                                    max_eukaryotic_right = right
                        else:
                            # Noncellular organism
                            if left < min_noncellular_left:
                                min_noncellular_left = left
                            if right > max_noncellular_right:
                                max_noncellular_right = right
                if min_cellular_left == max_all_left:
                    self.taxon[ortholog_type][threshold] = UniProtTaxonomy.objects.filter(
                      left_id__lte = min_noncellular_left,
                      right_id__gte = max_noncellular_right).order_by('right_id')[0]
                else:
                    # ignore noncellular
                    min_kingdom_lefts = [min_bacterial_left, min_archaeal_left,
                                        min_eukaryotic_left]
                    if len([left for left in min_kingdom_lefts \
                            if left <> max_all_left]) > 1:
                        # The encompassing taxon would be cellular organisms
                        # Return our fake taxon which describes which of the kingdoms
                        # is present
                        kingdom_strings = []
                        if min_bacterial_left <> max_all_left:
                            kingdom_strings = kingdom_strings + ['Bacteria']
                        if min_archaeal_left <> max_all_left:
                            kingdom_strings = kingdom_strings + ['Archaea']
                        if min_eukaryotic_left <> max_all_left:
                            kingdom_strings = kingdom_strings + ['Eukaryotes']
                        scientific_name = ', '.join(kingdom_strings)
                        self.taxon[ortholog_type][threshold] \
                            = UniProtTaxonomy.objects.get(scientific_name__exact
                                                              = scientific_name)
                    else:
                        self.taxon[ortholog_type][threshold] \
                          = UniProtTaxonomy.objects.filter(left_id__lte = min_cellular_left,
                              right_id__gte = max_cellular_right).order_by('right_id')[0]
        return self.taxon[ortholog_type][threshold]

    def is_family_root(self):
        return self.left_id == 1

    def is_leaf(self):
        return self.left_id + 1 == self.right_id

    def is_phog(self, threshold=0):
        return self.duplication_distance >= threshold and self.greatest_duplication_distance_of_maximal_descendant < threshold

    def get_seed_pfam(self):
        if self.left_id == 1 and self.tree.family.family_type.id == 'C':
            consensi = TreeNodeConsensus.objects.filter(tree_node = self, method =
                                                      'hmm')
            if consensi:
                consensus = consensi[0]
                hmm = consensus.hmm_consensus.hmm
                if hmm.pfam is not None:
                    return hmm.pfam
        return None

    def get_description(self, ortholog_type = OrthologTypes.SuperOrtholog,
        threshold = 0.0, returnAll = None, score_pfam=False):
        # For conserved region families, try to get the description from the
        # associated Pfam domain

        if returnAll is None:
            returnAll = False
        if not score_pfam:
            pfam = self.get_seed_pfam()
            if pfam is not None:
                return pfam.description
        if ((not score_pfam) and self.canonical_tree_node_name):
            return self.canonical_tree_node_name.name
        if self.description == '':
            if not self.is_leaf():
                if self.has_swissprot is None:
                    self.has_swissprot = False
                votes_for_description = {}
                votes_for_uninformative_description = {}
                raw_votes = {}
                for ortholog in self.get_included_leaves(ortholog_type, threshold):
                    description = ortholog.sequence_header.description()
                    if not raw_votes.has_key(description):
                        raw_votes[description] = 1
                    else:
                        raw_votes[description] += 1
                    if not description or uncharacterized_re.search(description) \
                        or predicted_re.search(description) or \
                        hypothetical_re.search(description):
                        if not votes_for_uninformative_description.has_key(description):
                            votes_for_uninformative_description[description] = 0
                        votes_for_uninformative_description[description] += 1
                    else:
                        if description not in votes_for_description:
                            votes_for_description[description] = 0
                        if ortholog.sequence_header.uniprot and \
                            ortholog.sequence_header.uniprot.in_swissprot_f:
                            self.has_swissprot = True
                            votes_for_description[description] += 10
                        # in some other 3rd party db that is not swissprot
                        else:
                            votes_for_description[description] += 1
                if len(votes_for_description) > 0:
                    used_votes_for_description = votes_for_description
                else:
                    used_votes_for_description = votes_for_uninformative_description
                winning_votes = 0
                winning_description = ''
                if returnAll:
                    returnlist = []
                    for key, value in used_votes_for_description.items():
                        returnlist.append((key, value, raw_votes[key]))
                    return returnlist
                for description, votes in used_votes_for_description.items():
                    # take the descriptions with the max votes; if there is a
                    # tie, take the shortest one(since it is likely to be the
                    # most general); if there is still a tie, take
                    # the lexicographically least one (this may not be
                    # happening...)
                    if votes > winning_votes or votes == winning_votes and \
                      (len(description) < len(winning_description) or
                        description < winning_description):
                        winning_votes = votes
                        winning_description = description
                self.description = winning_description
            # RSD: we should never enter the elif and not find a
            # sequence_header
            elif self.sequence_header:
                self.description = self.sequence_header.description()
                return self.description
            else:
                # this is problematic, since this will simply assign the id
                # to the description; we don't want to display the id
                self.description =  str(self)
        return self.description


    def get_description_for_caching(self, returnAllDescriptions=None):
        '''
        Returns either a string (if returnAllDescriptions is none) or a list
        containing tuples (description, score, and the number of members with
        this description)
        '''
        if returnAllDescriptions is None:
            returnAllDescriptions = False
        if returnAllDescriptions:
            descriptionlist = self.get_description(returnAll=True)
            sortedlist = sorted(descriptionlist, key=lambda score: score[1], reverse=True)
            retlist =[]
            for desc, score, raw in sortedlist:
                r = desc.strip()
                if r[-1] == ";":
                    retlist.append((str(r[0:-1]), score, raw))
                else:
                    retlist.append((str(r), score, raw))
            return retlist
        else:
            description = self.get_description().strip()
            if description[-1] == ";":
                description = description[0:-1]
            return description

    def get_has_swissprot(self):
        if self.has_swissprot is None:
            self.has_swissprot = False
            for ortholog in self.get_included_leaves():
                if ortholog.sequence_header.uniprot and \
                    ortholog.sequence_header.uniprot.in_swissprot_f:
                    self.has_swissprot = True
        return self.has_swissprot

    def get_accession(self, ortholog_type = OrthologTypes.SuperOrtholog,
                      threshold=0.0):
        '''
        Called in __unicode__ function when print is called on this node.
        '''
        if self.accession == '':
            if self.greatest_duplication_distance_of_maximal_descendant is not None\
                and (self.duplication_distance is None or
                    self.greatest_duplication_distance_of_maximal_descendant
                      < self.duplication_distance):
                self.accession = 'PHOG%07d_%05d' % (self.tree.id, self.left_id)
            elif self.right_id == self.left_id + 1:
                self.accession = 'LEAF%07d_%05d' % (self.tree.id, self.left_id)
            else:
                self.accession = 'NODE%07d_%05d' % (self.tree_id, self.left_id)
        return self.accession

    def get_tree_url(self, ortholog_type = OrthologTypes.SuperOrtholog,
                            threshold=0.0):
        return "/phog/tree/%s/" % self.get_accession(ortholog_type, threshold)
    def get_absolute_url(self, ortholog_type = OrthologTypes.SuperOrtholog,
                          threshold=0.0):
        return "/phog/%s/" % self.get_accession(ortholog_type, threshold)

    def annotations(self, type, uniprot_list=None):
        def go_annotations(term_type):
            return [ GOAnnotation({
                        'accession': a['ontology_term__accession'],
                        'description': a['ontology_term__name'],
                        'evidence_code' : a['evidence__code'],
                        'evidence_description' : a['evidence__name'],
                        'evidence_priority': a['evidence__priority'],
                    }) for a in Annotation.objects.filter(
                uniprot__sequenceheader__treenode__tree = self.tree,
                uniprot__sequenceheader__treenode__left_id__gte = self.left_id,
                uniprot__sequenceheader__treenode__right_id__lte = self.right_id,
                ontology_term__ontologysubsetmembership__subset = term_type
            ).values('ontology_term__accession', 'ontology_term__name', 'evidence__code', 'evidence__name', 'evidence__priority') ]


        if type == 'ec':
            return [ ec.as_annotation() for ec in  EC.objects.filter(
                uniprotec__uniprot__sequenceheader__treenode__tree = self.tree,
                uniprotec__uniprot__sequenceheader__treenode__left_id__gte = self.left_id,
                uniprotec__uniprot__sequenceheader__treenode__right_id__lte = self.right_id,
            ).all() ]
        if type == 'go_biological_process':
            return go_annotations("biological_process")
        if type == 'go_molecular_function':
            return go_annotations("molecular_function")
        if type == 'go_cellular_component':
            return go_annotations("cellular_component")
        if type == 'sfld':
            return [ sfld.as_annotation() for sfld in SFLDFamily.objects.filter(
                efds__sfldtopfacts__uniprot__sequenceheader__treenode__tree = self.tree,
                efds__sfldtopfacts__uniprot__sequenceheader__treenode__left_id__gte = self.left_id,
                efds__sfldtopfacts__uniprot__sequenceheader__treenode__right_id__lte = self.right_id,
            ).all() ]
        if type == 'uniprot_descriptions':
            return [ Annotation(u.de) for u in UniProt.objects.filter(
                sequenceheader__treenode__tree = self.tree,
                sequenceheader__treenode__left_id__gte = self.left_id,
                sequenceheader__treenode__right_id__lte = self.right_id,
            ).all() ]

    #@memoized
    def average_distance_down(self):
        if hasattr(self, '_average_distance_down') and True==False:
            return self._average_distance_down
        distances = []
        for leaf in self.get_included_leaves():
            distance = TreeNode.objects.filter(
                    tree=self.tree, 
                    left_id__gt=self.left_id, 
                    left_id__lte=leaf.left_id, 
                    right_id__lt=self.right_id, 
                    right_id__gte=leaf.right_id, 
            ).aggregate(Sum('distance_to_parent'))['distance_to_parent__sum']
            if distance != None:
                distances.append(distance)

        if len(distances):
            self._average_distance_down = sum(distances) / len(distances)
        else:
            self._average_distance_down = 0
        return self._average_distance_down

    def get_member_fasta(self):
        leaves = self.get_uniprot_in_members()
        sequences = []
        for leaf in leaves:
            if (leaf.fasta_from_db):
                sequences.append(leaf.fasta_from_db)
        return '\n'.join(sequences)

    def __unicode__(self):
        return unicode(self.get_accession() or self.id)

    class Meta:
        db_table = u'tree_node'

class TreeNodeConsensus(models.Model):
    id = models.AutoField(primary_key=True)
    tree_node = models.ForeignKey(TreeNode)
    sequence = models.ForeignKey(Sequence)
    method = models.TextField()
    uppercase_threshold = models.FloatField(null=True)
    hmm_consensus = models.ForeignKey('HMM_Consensus', null=True)
    class Meta:
        db_table = u'tree_node_consensus'

class UniProtQuerySet(QuerySet):
    def orthologs(self, threshold):
        return UniProt.objects.filter(
              phogmembership__ancestor_tree_node__phog_members__uniprot__in = self,
              phogmembership__ancestor_tree_node__phog_members__uniprot__accession__isnull=False, #This is to force a join with the uniprot table back again because the Django ORM is retarded and they insist it is perfect
              phogmembership__phog_duplication_distance__gte = threshold,
              phogmembership__phog_greatest_duplication_distance_of_maximal_descendant__lt = threshold
        )

class UniProtGene(models.Model):
    id = models.AutoField(primary_key=True)
    uniprot = models.ForeignKey('UniProt')
    name = models.TextField()
    name_type = models.TextField()

    class Meta:
        db_table = u'uniprot_gene'

class UniProtAlternateAccessions(models.Model):
    id = models.AutoField(primary_key=True)
    primary_accession = models.TextField();
    alternate_accession = models.TextField();

    class Meta:
        db_table = u'uniprot_alternate_accession'

class UniProt(models.Model):
    id = models.AutoField(primary_key=True)
    uniprot_identifier = models.CharField(max_length=12)
    accession = models.CharField(unique=True, max_length=6)
    taxon = models.ForeignKey('UniProtTaxonomy')
    de = models.CharField(max_length=255)
    seguid = models.TextField() # This field type is a guess.
    in_swissprot_f = models.BooleanField()
    description = models.TextField()
    is_fragment = models.BooleanField()
    is_precursor = models.BooleanField()
    pfam = models.ManyToManyField('Pfam', db_table='uniprot_pfam')
    
    objects = PassThroughManager.for_queryset_class(UniProtQuerySet)()
    sequence_chars = models.TextField()
    sequence_length = models.IntegerField()
    existence_proof = models.TextField()
    is_deleted = models.BooleanField()

    def annotations(self, type):
        def go_annotations(term_type):
            return [ GOAnnotation({
                        'accession': a['ontology_term__accession'],
                        'description': a['ontology_term__name'],
                        'evidence_code' : a['evidence__code'],
                        'evidence_description' : a['evidence__name'],
                        'evidence_priority': a['evidence__priority'],
                    }) for a in Annotation.objects.filter(
                uniprot = self,
                ontology_term__ontologysubsetmembership__subset = term_type
            ).distinct('ontology_term__accession').values('ontology_term__accession', 'ontology_term__name', 'evidence__code', 'evidence__name', 'evidence__priority') ]


        if type == 'ec':
            return [ ec.as_annotation() for ec in  EC.objects.filter(
                uniprotec__uniprot = self,
            ).all() ]
        if type == 'go_biological_process':
            return go_annotations("biological_process")
        if type == 'go_molecular_function':
            return go_annotations("molecular_function")
        if type == 'go_cellular_component':
            return go_annotations("cellular_component")
        if type == 'sfld':
            return [ sfld.as_annotation() for sfld in SFLDFamily.objects.filter(
                efds__sfldtopfacts__uniprot = self,
            ).all() ]
        if type == 'uniprot_descriptions':
            return [ Annotation(u.de) for u in UniProt.objects.filter(
                sequenceheader = self,
            ).all() ]

    @classmethod
    def find_by_accession_or_identifier(cls, acc_or_ident):
        try:
            return cls.objects.get(accession = acc_or_ident)
        except cls.DoesNotExist:
            return  cls.objects.get(uniprot_identifier = acc_or_ident)
        except cls.DoesNotExist:
            raise

    @classmethod
    def find_by_sequence(cls, seq):
        hashobj = hashlib.sha1()
        hashobj.update(seq)
        seqhash = base64.b64encode(hashobj.digest()).strip('=')
        retobj = cls.objects.filter(seguid=seqhash)
        if len(retobj) != 0:
            return cls.objects.filter(seguid=seqhash)[0]
        else:
            return []

    @property
    def has_medium_phogs(self):
        if not hasattr(self, '_has_medium_phogs'):
            self._has_medium_phogs = bool(self.sequenceheader_set.all().filter(
                treenode__phogt_medium_node__isnull=False)[:1])
        return self._has_medium_phogs

    @property
    def fasta(self):
        '''This is pretty terrible'''
        if not hasattr(self, '_fasta'):
            self._fasta = subprocess.Popen(['/clusterfs/ohana/software/bin/fastacmd',
                '-d', '/clusterfs/ohana/external/UniProt/current/protein',
                '-s', self.uniprot_identifier,
            ], stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0] \
            or None
        return self._fasta

    @property
    def fasta_from_db(self):
        ''' Source fasta from the database instead '''
        try:
            return self._fasta_from_db
        except:
            pass

        if (self.sequence_chars):
            aa_seq = self.sequence_chars
        elif (self.chars):
            aa_seq = self.chars
        else:
            self._fasta_from_db = None
            return self._fasta_from_db
        if self.in_swissprot_f:
            db = "sp"
        else:
            db = "tr"
        sequence = ">" + db + "|" + self.accession + "|" + self.uniprot_identifier + " " + self.description + " OS=%s" % self.taxon.scientific_name
        genes = UniProtGene.objects.filter(uniprot = self)
        if genes:
            sequence += " GN=%s\n" % genes[0].name
        else:
            sequence += "\n"
        for (idx, char) in enumerate(aa_seq):
            if ((idx % 60) == 59):
                sequence += '\n'
            sequence += char 
        self._fasta_from_db = sequence
        return self._fasta_from_db

    @property
    def sequence_len(self):
        if self.sequence_length:
            return self.sequence_length
        elif self.chars:
            return len(self.chars)
        return None
                    
    @property
    def wrapped_chars(self):
        return self.fasta and self.fasta.split('\n', 1)[1] or None

    @property
    def chars(self):
        return self.wrapped_chars and self.wrapped_chars.replace('\n', '')

    def get_pfacts_url(self):
        '''kludgy, but get_absolute_url is already taken by an outside link'''
        return '/phylofacts/sequence/UniProt/%s/' % self.uniprot_identifier

    class Meta:
        db_table = u'uniprot'

    def __init__(self, *args, **kwargs):
        models.Model.__init__(self, *args, **kwargs)
        self.partners = UniProt.objects.none()
        self.have_sought_partners = False
        self.have_summarized_go=False
        self.experimental_go_summary = {}
        self.nonexperimental_go_summary = {}

    def get_treenodes(self, method = 'ml'):
        ''' Returns the treenodes that are associated with this UniProt object. '''
        families = self.ghg_families()
        for domain, domain_object in self.pfam_families().items():
            families += domain_object['families']
        return list(TreeNode.objects.filter(tree__family__in = families, 
                    sequence_header__uniprot = self, tree__family__status = "draft", 
                    tree__family__active = True, tree__method = method))

    def get_orthology_groups(self):
        try:
            return self._get_orthology_groups
        except:
            pass
        treenodes = self.get_treenodes()
        orthology_dict = {}
        for node in treenodes:
            family_accession = node.tree.family.get_accession()
            if family_accession not in orthology_dict:
                orthology_dict[family_accession] = node.orthology_supported_ancestors()
            else:
                for key in orthology_dict[family_accession].keys():
                    # get the largest orthology group containing this sequence in this family, 
                    # if there are multiple copies of a sequence in the family
                    # this may not make the most sense....
                    if (node.orthology_supported_ancestors()[key] and 
                        (not orthology_dict[family_accession][key] or 
                        (orthology_dict[family_accession][key].left_id > 
                        node.orthology_supported_ancestors()[key].left_id))):
                        orthology_dict[family_accession][key] = node.orthology_supported_ancestors()[key]
        self._get_orthology_groups = orthology_dict
        return self._get_orthology_groups

    def get_pdb_structures(self):
        return PDB.objects.filter( pdb_chain__uniprotpdb_chain__uniprot = self ).distinct()


    def pfam_families_for_current_object(self):
        '''
        Get Pfam families containing the current uniprot object
        '''
        families = Family.objects.filter(
          canonical_tree__treenodes__sequence_header__uniprot = self,
          family_type__id = 'C',
          active = True
          ).exclude(
          status = "bad"
          )
        return families

    def get_families_containing_me(self):
        ''' Another function to get all active not bad families containing the uniprot object in question. '''
        return Family.objects.filter(id__in=TreeNode.objects.filter(sequence_header__uniprot = self,
                    tree__family__active=True, tree__family__status = 'draft',
                    tree__family__canonical_tree=F('tree')).values_list('tree__family__id', 
                    flat=True).distinct('tree__family__id'))

    # The next two methods return the families in the score order, and make sure
    # that all of them have a score (this is why we can't just sort them with sql)

    # returns a hash of lists, the keys being the pfam domains
    def pfam_families(self):
        '''
        Get PFAM family map for the current uniprot object
        '''
        families = self.pfam_families_for_current_object()
        pfam_family_map = {}
        for family in families:
            pfam = list(family.get_pfams())[0][0] # extract the pfam object
            pfam_accession = pfam.accession # there are
            # different versions of pfam, so the same
            # accession can be found in differently
            # versioned objects
            if pfam_accession not in pfam_family_map:
                pfam_family_map[pfam_accession] = {}
                pfam_family_map[pfam_accession]['pfam_object'] = pfam
                pfam_family_map[pfam_accession]['families'] = set()
            pfam_family_map[pfam_accession]['families'].add(family)

        for pfam_accession in pfam_family_map:
            pfam_family_map[pfam_accession]['families'] = list(pfam_family_map[pfam_accession]['families'])
            pfam_family_map[pfam_accession]['families'].sort(key = lambda f: -1 * f.get_score()  )

        return pfam_family_map

    # returns a list
    def ghg_families(self):
        families = Family.objects.filter(
                canonical_tree__treenodes__sequence_header__uniprot = self,
                family_type__id = 'G',
                active = True
            ).exclude(
                status = "bad"
            )
        return sorted(families, key=lambda family: -1 * family.get_score())

    def pfam_architecture(self):
        '''Reads and returns the pfam architecture for this uniprot protein.'''
        uniprot_pfam_objs = UniProtPfam.objects.filter(uniprot=self)
        if not uniprot_pfam_objs: return None
        # Get a list of pfam names and coordinates as [(pfam_name, (start, end)), ...]
        pfam_domain_list = []
        for uniprot_pfam_obj in uniprot_pfam_objs:
            for coords in uniprot_pfam_obj.coordinates.all():
                pfam_domain_list.append((uniprot_pfam_obj.pfam.name,
                                         (coords.sequence_start, coords.sequence_end)))
        if not pfam_domain_list: return None
        # Get the architecture by first sorting the pfam_domain_list by the start values,
        # obtaining the first value in the set and grouping into an architecture
        pfam_domain_list = sorted(pfam_domain_list, key=lambda x: x[1][0])
        architecture = zip(*pfam_domain_list)[0]
        architecture = ', '.join([x[0] for x in groupby(architecture)])
        return (architecture, {'uniprot_accession' : self.accession,
                            'pfam_coordinates' : pfam_domain_list})

    def get_orthomcl_group(self):
        '''Returns the orthomcl group that this protein belongs to, by accessing the
        associated OrthoMclGroupUniProt object.'''
        return self.orthomcl_group.all().distinct().values_list(
            'orthomcl_group_id', flat=True)

    def get_oma_group(self):
        '''Returns the oma group that this protein belongs to, by accessing the
        associated OMAGroupUniProt object.'''
        return self.oma_group.all().distinct().values_list(
            'oma_group_id', flat=True)

    def get_sfld_group(self):
        '''Returns the sfld group that this protein belongs to..
        '''
        return SFLDFamily.objects.filter(
            efds__sfldtopfacts__uniprot = self
            )

    def get_usual_summaries(self):
        '''This needs work, but first get_go_summary needs to return UniProt_GO
        objects, instead of a tuple of strings.'''
        summaries = {}
        for x in ('biological_process', 'molecular_function',
            'cellular_component'):
            summary = (self.get_go_summary(x, experimental=True) \
                      or []) \
                    + (self.get_go_summary(x, experimental=False) \
                      or [])
            if summary:
                summaries[x] = []
                for (acc, name, evidence, long_evidence) in summary:
                    summaries[x] = summaries[x] + [\
                      { 'go_term': {'acc': acc, 'name': name}, \
                        'go_evidence': {'evidence': evidence, \
                                      'name': long_evidence}}]
        return summaries

    def get_interacting_partners(self):
        if not self.have_sought_partners:
            irefindices1 = IRefIndexSimple.objects.filter(uniprot_a = self)
            partner_set = set([irefindex.uniprot_b.id for irefindex in irefindices1])
            irefindices2 = IRefIndexSimple.objects.filter(uniprot_b = self)
            partner_set = partner_set | \
                          set([irefindex.uniprot_a.id for irefindex in irefindices2])
            self.partners = UniProt.objects.filter(id__in=list(partner_set))
            self.have_sought_partners = True
        return self.partners

    def get_pfam_hits(self):
        return Pfam25PfamA_RegFullSignificant.objects.filter(pfamseq__in =
                    Pfam25PfamSeq.objects.filter(pfamseq_acc = self.accession))

    def has_experimental_evidence(self):
        evidence_codes = set([annotation.go_evidence.evidence for annotation in
                              self.go_annotations.all()])
        return (len(evidence_codes & experimental_evidence_codes) > 0)

    def has_literature(self):
        return self.literature.exclude(is_large_scale_f = True).count() > 0

    def kegg_sequence_id(self): #not yet implemented
        for kegg in self.keggs.all():
            return kegg.kegg_sequence_id
        return ''

    def biocyc_id(self): #not yet implemented
        for biocyc in self.biocycs.all():
            return biocyc.biocyc_id
        return ''

    def get_description(self):
        return self.de or self.description or 'N/A'

    def __unicode__(self):
        return self.in_swissprot_f and self.uniprot_identifier or self.accession or self.de or 'N/A'

    def get_evidence_icon(self):
        return self.has_experimental_evidence() and uniprot_evidence_icon(self.accession) or ''

    def get_literature_icon(self):
        return self.has_literature() and uniprot_lit_icon(self.accession) or ''

    def get_swissprot_icon(self):
        return self.in_swissprot_f and swissprot_icon(self.accession) or ''

    def get_ec(self):
        uniprot_ec_objects = UniProtEC.objects.filter(uniprot__exact = self)
        if uniprot_ec_objects:
            return (uniprot_ec_objects[0].ec, uniprot_ec_objects[0].is_in_brenda_f)
        else:
            return (None, False)

    def get_ec_link(self):
        ec, is_in_brenda = self.get_ec()
        if ec:
            return uniprot_ec_link(ec, experimental=is_in_brenda)
        else:
            return ''

    def get_kegg_map_ids(self):
        ec, is_in_brenda = self.get_ec()
        if ec:
            kegg_map_ec_objects = KEGG_Map_EC.objects.filter(ec__exact = ec)
            if kegg_map_ec_objects:
                return [kegg_map_ec.kegg_map for kegg_map_ec in kegg_map_ec_objects]
            else:
                return []
        else:
            return []


    def get_kegg_links(self):
        kegg_maps = self.get_kegg_map_ids()
        if len(kegg_maps) > 0:
            return make_kegg_map_links(kegg_maps)
        else:
            return ''

    def get_ppi_link(self): #not yet implemented
        return bool(self.get_interacting_partners()) and uniprot_ppi_link(self.uniprot_identifier) or ''

    def get_netscope_url(self, ortholog_type = 0, threshold = 0.0): #not yet implemented
        return bool(self.get_interacting_partners()) and netscope_url(self.uniprot_identification, ortholog_type, threshold) or ''

    def get_orthologs_url(self):
        return uniprot_orthologs_url(self.accession)

    def get_orthologs_link(self):
        return uniprot_orthologs_link(self.accession)

    def get_absolute_url(self):
        return uniprot_url(self.accession)

    def get_link(self):
        return uniprot_link(self.accession)

    def __unicode__(self):
        return (self.in_swissprot_f and self.uniprot_identifier or self.accession) or self.de or 'N/A'

class UniProtIDMapping(models.Model):
    uniprot_accession = models.TextField(primary_key=True)
    database = models.TextField()
    foreign_id = models.TextField()
    class Meta:
        db_table = u'uniprot_id_mapping'

    @classmethod
    def find_uniprot(self, foreign_id):
        '''Tries to find a uniprot object in our DB based on a foreign key
            returns None if None is found
        '''

        hits = self.objects.filter(foreign_id=foreign_id)

        for hit in hits:
            uniprot_object = UniProt.objects.get(accession=hit.uniprot_accession)
            if uniprot_object:
                return uniprot_object
        return None


class PPI_DetectionMethod(models.Model):
    id = models.AutoField(primary_key=True)
    mi_ontology_id = models.CharField(max_length=40, null=True)
    short_label = models.CharField(max_length=255)
    class Meta:
        db_table = u'ppi_detection_method'

class PPI_SourceDatabase(models.Model):
    id = models.AutoField(primary_key=True)
    mi_ontology_id = models.CharField(max_length=40, null=True)
    source_name = models.CharField(max_length=255)
    class Meta:
        db_table = u'ppi_source_database'

class IRefIndexSimple(models.Model):
    id = models.AutoField(primary_key=True)
    uniprot_a = models.ForeignKey(UniProt, related_name = 'a_interactions')
    uniprot_b = models.ForeignKey(UniProt, related_name = 'b_interactions')
    ppi_detection_method = models.ForeignKey(PPI_DetectionMethod, null = True)
    ppi_source_database = models.ForeignKey(PPI_SourceDatabase)
    interaction_identifier_in_source_db = models.CharField(max_length=255)
    interaction_type = models.CharField(max_length=1)
    num_participants = models.IntegerField()
    class Meta:
        db_table = u'irefindex_simple'

class CSA_Residue(models.Model):
    id = models.AutoField(primary_key=True)
    res_num = models.IntegerField()
    site_num = models.IntegerField()
    res_type = models.CharField(max_length=3)
    res_inc = models.CharField(max_length=1)
    catalytic = models.BooleanField()
    hetero = models.BooleanField()
    chem_func = models.CharField(max_length=3)
    uniprot = models.ForeignKey(UniProt)
    notes = models.TextField(null=True)
    pdb_chain = models.ForeignKey('PDB_Chain')
    class Meta:
        db_table = u'csa_residue'

class CSA_ResidueFunction(models.Model):
    id = models.AutoField(primary_key=True)
    residue = models.ForeignKey(CSA_Residue)
    function = models.TextField()
    target_type = models.TextField()
    class Meta:
        db_table = u'csa_residue_function'

class CSA_ResidueTarget(models.Model):
    id = models.AutoField(primary_key=True)
    residue = models.ForeignKey(CSA_Residue)
    target_chain_id = models.CharField(max_length=1)
    target_res_num = models.IntegerField()
    target_hetero = models.BooleanField()
    target_res_type = models.CharField(max_length=3)
    class Meta:
        db_table = u'csa_residue_target'

class CSA_EvidenceType(models.Model):
    id = models.AutoField(primary_key=True)
    type = models.CharField(max_length=80)
    class Meta:
        db_table = u'csa_evidence_type'

class CSA_Evidence(models.Model):
    id = models.AutoField(primary_key=True)
    type = models.ForeignKey(CSA_EvidenceType)
    residue = models.ForeignKey(CSA_Residue)
    source = models.CharField(max_length=10)
    score = models.CharField(max_length=12)
    created = models.DateTimeField(default=datetime.now)
    updated = models.DateTimeField(default=datetime.now)

    class Meta:
        db_table = u'csa_evidence'


class Draft(models.Model):
    id = models.AutoField(primary_key=True)
    family = models.ForeignKey(Family)
    ready_for_review = models.BooleanField(default=False)
    reviewed = models.BooleanField(default=False)
    checked_out = models.BooleanField()
    checked_out_in_by = models.ForeignKey(User)
    check_out_in_time = models.DateTimeField(default=datetime.now)
    class Meta:
        db_table = u'draft'

class FamilyChecked(models.Model):
    id = models.AutoField(primary_key=True)
    family = models.ForeignKey(Family)
    check_type = models.TextField()
    passed = models.BooleanField(default=False)
    updated = models.DateTimeField(default=datetime.now)
    notes = models.TextField(null=True)
    class Meta:
        db_table = u'family_checked'

class FamilyEdit(models.Model):
    id = models.AutoField(primary_key=True)
    family = models.ForeignKey(Family)
    editor = models.ForeignKey(User)
    updated = models.DateTimeField(default=datetime.now)

    class Meta:
        db_table = u'family_edit'

class Village(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(unique=True, max_length=35)
    display_name = models.CharField(max_length=255)
    password = models.CharField(max_length=255, null=True)
    user_specific = models.BooleanField(default=False)
    class Meta:
        db_table = u'village'

class FamilyVillage(models.Model):
    family = models.ForeignKey(Family)
    village = models.ForeignKey(Village)
    class Meta:
        db_table = u'family_village'

class Genbank(models.Model):
    id = models.CharField(max_length=255, primary_key=True)
    accession = models.CharField(max_length=255)
    version = models.IntegerField()
    taxon = models.ForeignKey(UniProtTaxonomy)
    description = models.TextField(null=True)
    sequence = models.ForeignKey(Sequence, null=True)
    class Meta:
        db_table = u'genbank'

class GenbankPDB(models.Model):
    id = models.AutoField(primary_key=True)
    genbank = models.ForeignKey(Genbank)
    pdb_chain = models.ForeignKey('PDB_Chain')
    class Meta:
        db_table = u'genbank_pdb'

class ThirdParty(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.TextField()
    urlformat = models.TextField()
    taxon = models.ForeignKey('UniProtTaxonomy', null=True)
    class Meta:
        db_table = u'third_party'

class ThirdPartySequenceHeader(models.Model):
    id = models.AutoField(primary_key=True)
    thirdparty = models.ForeignKey(ThirdParty)
    sequence_header = models.ForeignKey(SequenceHeader)
    class Meta:
        db_table = u'third_party_sequence_header'

class GenbankThirdPartySequenceHeader(models.Model):
    id = models.AutoField(primary_key=True)
    genbank = models.ForeignKey(Genbank)
    thirdparty_sequence_header = models.ForeignKey(ThirdPartySequenceHeader)
    class Meta:
        db_table = u'genbank_thirdparty_sequence_header'

class GO_AssocRel(models.Model):
    id = models.IntegerField(primary_key=True)
    from_id = models.IntegerField()
    to_id = models.IntegerField()
    relationship_type_id = models.IntegerField()
    class Meta:
        db_table = u'go_assoc_rel'

class GO_Association(models.Model):
    id = models.IntegerField(primary_key=True)
    term = models.ForeignKey('GO_Term')
    gene_product = models.ForeignKey('GO_GeneProduct')
    is_not = models.IntegerField(null=True)
    role_group = models.IntegerField(null=True)
    assocdate = models.IntegerField(null=True)
    source_db = models.ForeignKey('GO_DB', null=True)
    class Meta:
        db_table = u'go_association'

class GO_AssociationProperty(models.Model):
    id = models.AutoField(primary_key=True)
    association = models.ForeignKey(GO_Association)
    relationship_type = models.ForeignKey('GO_Term',
       related_name='go_association_relationship_type')
    term = models.ForeignKey('GO_Term',
       related_name='go_association_term')

    class Meta:
        db_table = u'go_association_property'

class GO_AssociationQualifier(models.Model):
    id = models.IntegerField(primary_key=True)
    association_id = models.IntegerField()
    term_id = models.IntegerField()
    value = models.CharField(max_length=255)
    class Meta:
        db_table = u'go_association_qualifier'

class GO_AssociationSpeciesQualifier(models.Model):
    id = models.IntegerField(primary_key=True)
    association_id = models.IntegerField()
    species_id = models.IntegerField()
    class Meta:
        db_table = u'go_association_species_qualifier'

class GO_DB(models.Model):
    id = models.IntegerField(primary_key=True)
    name = models.CharField(max_length=55)
    fullname = models.CharField(max_length=255)
    datatype = models.CharField(max_length=255)
    generic_url = models.CharField(max_length=255)
    url_syntax = models.CharField(max_length=255)
    url_example = models.CharField(max_length=255)
    uri_prefix = models.CharField(max_length=255)
    class Meta:
        db_table = u'go_db'

class GO_DBxref(models.Model):
    id = models.AutoField(primary_key=True)
    xref_dbname = models.CharField(max_length=55)
    xref_key = models.CharField(max_length=255)
    xref_keytype = models.CharField(max_length=32)
    xref_desc = models.CharField(max_length=255)
    class Meta:
        db_table = u'go_dbxref'

class GO_Evidence(models.Model):
    id = models.IntegerField(primary_key=True)
    code = models.CharField(max_length=8)
    association_id = models.IntegerField()
    dbxref_id = models.IntegerField()
    seq_acc = models.CharField(max_length=255)

    def __unicode__(self):
        return "code = %s; association_id = %s; dbxref_id = %s; seq_acc = %s" % (self.code, self.association_id, self.dbxref_id, self.seq_acc)

    __str__ = __unicode__

    def get_absolute_url(self):
        return "http://www.geneontology.org/GO.evidence.shtml#%s" % self.code.lower()

    class Meta:
        db_table = u'go_evidence'

class GO_EvidenceDbxref(models.Model):
    evidence_id = models.IntegerField()
    dbxref_id = models.IntegerField()
    class Meta:
        db_table = u'go_evidence_dbxref'

class GO_EvidencePriority(models.Model):

    # if GO evidence codes ever change, this constant should change too.
    # (ideally this would be a staticmethod property, but these functions
    # don't play well together)
    experimental_priority = 5

    evidence = models.CharField(unique=True, max_length=3)
    name = models.CharField(max_length=100)
    priority = models.IntegerField(unique=True)
    id = models.IntegerField(primary_key=True)

    @classmethod
    def lookup_by_priority(cls):
        lookup = {}
        objects = cls.objects.all()
        for object in objects:
            lookup[priority] = [object.evidence, object.name]
        return lookup

    def __unicode__(self):
        return "evidence = %s; name = %s; priority = %s" % (self.evidence, self.name, self.priority)

    __str__ = __unicode__


    class Meta:
        db_table = u'go_evidence_priority'


class GO_GeneProduct(models.Model):
    id = models.IntegerField(primary_key=True)
    symbol = models.CharField(max_length=128)
    dbxref_id = models.IntegerField()
    species_id = models.IntegerField()
    type_id = models.IntegerField()
    full_name = models.TextField()
    class Meta:
        db_table = u'go_gene_product'

class GO_GeneProductAncestor(models.Model):
    gene_product_id = models.IntegerField()
    ancestor_id = models.IntegerField()
    class Meta:
        db_table = u'go_gene_product_ancestor'

class GO_GeneProductCount(models.Model):
    term_id = models.IntegerField()
    code = models.CharField(max_length=8)
    speciesdbname = models.CharField(max_length=55)
    species_id = models.IntegerField()
    product_count = models.IntegerField()
    class Meta:
        db_table = u'go_gene_product_count'

class GO_GeneProductHomology(models.Model):
    gene_product1_id = models.IntegerField()
    gene_product2_id = models.IntegerField()
    relationship_type_id = models.IntegerField()
    class Meta:
        db_table = u'go_gene_product_homology'

class GO_GeneProductHomolset(models.Model):
    id = models.IntegerField(primary_key=True)
    gene_product_id = models.IntegerField()
    homolset_id = models.IntegerField()
    class Meta:
        db_table = u'go_gene_product_homolset'

class GO_GeneProductProperty(models.Model):
    gene_product_id = models.IntegerField()
    property_key = models.CharField(max_length=64)
    property_val = models.CharField(max_length=255)
    class Meta:
        db_table = u'go_gene_product_property'

class GO_GeneProductSeq(models.Model):
    gene_product_id = models.IntegerField()
    seq_id = models.IntegerField()
    is_primary_seq = models.IntegerField()
    class Meta:
        db_table = u'go_gene_product_seq'

class GO_GeneProductSubset(models.Model):
    gene_product_id = models.IntegerField()
    subset_id = models.IntegerField()
    class Meta:
        db_table = u'go_gene_product_subset'

class GO_GeneProductSynonym(models.Model):
    gene_product_id = models.IntegerField()
    product_synonym = models.CharField(max_length=255)
    class Meta:
        db_table = u'go_gene_product_synonym'

class GO_GraphPath(models.Model):
    id = models.IntegerField(primary_key=True)
    term1_id = models.IntegerField()
    term2_id = models.IntegerField()
    relationship_type_id = models.IntegerField()
    distance = models.IntegerField()
    relation_distance = models.IntegerField()
    class Meta:
        db_table = u'go_graph_path'

class GO_GraphPath2Term(models.Model):
    graph_path_id = models.IntegerField()
    term_id = models.IntegerField()
    rank = models.IntegerField()
    class Meta:
        db_table = u'go_graph_path2term'

class GO_Homolset(models.Model):
    id = models.IntegerField(primary_key=True)
    symbol = models.CharField(max_length=128)
    dbxref_id = models.IntegerField()
    target_gene_product_id = models.IntegerField()
    taxon_id = models.IntegerField()
    type_id = models.IntegerField()
    description = models.TextField()
    class Meta:
        db_table = u'go_homolset'

class GO_InstanceData(models.Model):
    release_name = models.CharField(max_length=255)
    release_type = models.CharField(max_length=255)
    release_notes = models.TextField()
    class Meta:
        db_table = u'go_instance_data'

class GO_RelationComposition(models.Model):
    id = models.IntegerField(primary_key=True)
    relation1_id = models.IntegerField()
    relation2_id = models.IntegerField()
    inferred_relation_id = models.IntegerField()
    class Meta:
        db_table = u'go_relation_composition'

class GO_RelationProperties(models.Model):
    relationship_type_id = models.IntegerField()
    is_transitive = models.IntegerField()
    is_symmetric = models.IntegerField()
    is_anti_symmetric = models.IntegerField()
    is_cyclic = models.IntegerField()
    is_reflexive = models.IntegerField()
    is_metadata_tag = models.IntegerField()
    class Meta:
        db_table = u'go_relation_properties'

class GO_Seq(models.Model):
    id = models.IntegerField(primary_key=True)
    display_id = models.CharField(max_length=64)
    description = models.CharField(max_length=255)
    seq = models.TextField()
    seq_len = models.IntegerField()
    md5checksum = models.CharField(max_length=32)
    moltype = models.CharField(max_length=25)
    timestamp = models.IntegerField()
    class Meta:
        db_table = u'go_seq'

class GO_SeqDbxref(models.Model):
    seq_id = models.IntegerField()
    dbxref_id = models.IntegerField()
    class Meta:
        db_table = u'go_seq_dbxref'

class GO_SeqProperty(models.Model):
    id = models.IntegerField(primary_key=True)
    seq_id = models.IntegerField()
    property_key = models.CharField(max_length=64)
    property_val = models.CharField(max_length=255)
    class Meta:
        db_table = u'go_seq_property'

class GO_SourceAudit(models.Model):
    source_id = models.CharField(max_length=255)
    source_fullpath = models.CharField(max_length=255)
    source_path = models.CharField(max_length=255)
    source_type = models.CharField(max_length=255)
    source_md5 = models.TextField() # This field type is a guess.
    source_parsetime = models.IntegerField()
    source_mtime = models.IntegerField()
    class Meta:
        db_table = u'go_source_audit'

class GO_Species(models.Model):
    id = models.IntegerField(primary_key=True)
    ncbi_taxa_id = models.IntegerField()
    common_name = models.CharField(max_length=255)
    lineage_string = models.TextField()
    genus = models.CharField(max_length=55)
    species = models.CharField(max_length=255)
    parent_id = models.IntegerField()
    left_value = models.IntegerField()
    right_value = models.IntegerField()
    taxonomic_rank = models.CharField(max_length=255)
    class Meta:
        db_table = u'go_species'

class GO_Term(models.Model):
    id = models.IntegerField(primary_key=True)
    name = models.CharField(max_length=255)
    term_type = models.CharField(max_length=55)
    acc = models.CharField(max_length=255)
    is_obsolete = models.IntegerField()
    is_root = models.IntegerField()
    is_relation = models.IntegerField()

    def __unicode__(self):
        return ('''name = %s; term_type = %s; acc = %s; is_obsolete = %s;
        is_root = %s, is_relation = %s''' % (self.name, self.term_type, self.acc,
        self.is_obsolete, self.is_root, self.is_relation))

    __str__ = __unicode__

    class Meta:
        db_table = u'go_term'

class GO_Term2Term(models.Model):
    id = models.IntegerField(primary_key=True)
    relationship_type_id = models.IntegerField()
    term1_id = models.IntegerField()
    term2_id = models.IntegerField()
    complete = models.IntegerField()
    class Meta:
        db_table = u'go_term2term'

class GO_Term2TermMetadata(models.Model):
    id = models.IntegerField(primary_key=True)
    relationship_type_id = models.IntegerField()
    term1_id = models.IntegerField()
    term2_id = models.IntegerField()
    class Meta:
        db_table = u'go_term2term_metadata'

class GO_TermAudit(models.Model):
    term_id = models.IntegerField()
    term_loadtime = models.IntegerField()
    class Meta:
        db_table = u'go_term_audit'

class GO_TermDbxref(models.Model):
    term_id = models.AutoField(primary_key=True)
    dbxref_id = models.IntegerField()
    is_for_definition = models.IntegerField()
    class Meta:
        db_table = u'go_term_dbxref'

class GO_TermDefinition(models.Model):
    term_id = models.IntegerField()
    term_definition = models.TextField()
    dbxref_id = models.IntegerField()
    term_comment = models.TextField()
    reference = models.CharField(max_length=255)
    class Meta:
        db_table = u'go_term_definition'

class GO_TermProperty(models.Model):
    term_id = models.IntegerField()
    property_key = models.CharField(max_length=64)
    property_val = models.CharField(max_length=255)
    class Meta:
        db_table = u'go_term_property'

class GO_TermSubset(models.Model):
    term_id = models.IntegerField()
    subset_id = models.IntegerField()
    class Meta:
        db_table = u'go_term_subset'

class GO_TermSynonym(models.Model):
    term_id = models.IntegerField()
    term_synonym = models.CharField(max_length=996)
    acc_synonym = models.CharField(max_length=255)
    synonym_type_id = models.IntegerField()
    synonym_category_id = models.IntegerField()
    class Meta:
        db_table = u'go_term_synonym'

class Pfam(models.Model):
    id = models.AutoField(primary_key=True)
    accession = models.CharField(max_length=255)
    name = models.CharField(max_length=255)
    description = models.TextField()
    hmmlength = models.IntegerField()
    hmm_type = models.CharField(max_length=8, null=True)
    version = models.IntegerField(null=True)
    overall_pfam_version = models.FloatField()
    pfam_clan_id = models.IntegerField()
    def __unicode__(self):
        return u'%s.%d: %s' % (self.accession, self.version, self.name)
    def __str__(self):
        return self.__unicode__()
    def get_absolute_url(self):
        return pfam_url(self.name)
    def get_link(self):
        if self.name in ['TM', 'SP']:
            return self.name
        else:
            return pfam_link(self.name, self.description)
    def get_clan(self):
        return PfamClan.objects.filter(id=self.pfam_clan_id).distinct()
    class Meta:
        db_table = u'pfam'

class PfamClan(models.Model):
    id = models.IntegerField(primary_key=True)
    description = models.TextField()
    accession = models.TextField()
    class Meta:
        db_table = u'pfam_clan'

class Pfam25PfamA(models.Model):
    auto_pfama = models.IntegerField(primary_key=True)
    pfama_acc = models.CharField(unique=True, max_length=7)
    pfama_id = models.CharField(unique=True, max_length=16)
    previous_id = models.TextField()
    description = models.CharField(max_length=100)
    author = models.TextField()
    deposited_by = models.CharField(max_length=100)
    seed_source = models.TextField()
    type = models.TextField() # This field type is a guess.
    comment = models.TextField()
    sequence_ga = models.FloatField()
    domain_ga = models.FloatField()
    sequence_tc = models.FloatField()
    domain_tc = models.FloatField()
    sequence_nc = models.FloatField()
    domain_nc = models.FloatField()
    buildmethod = models.TextField()
    model_length = models.IntegerField()
    searchmethod = models.TextField()
    msv_lambda = models.FloatField()
    msv_mu = models.FloatField()
    viterbi_lambda = models.FloatField()
    viterbi_mu = models.FloatField()
    forward_lambda = models.FloatField()
    forward_tau = models.FloatField()
    num_seed = models.IntegerField()
    num_full = models.IntegerField()
    updated = models.DateTimeField()
    created = models.DateTimeField()
    version = models.SmallIntegerField()
    number_archs = models.IntegerField()
    number_species = models.IntegerField()
    number_structures = models.IntegerField()
    number_ncbi = models.IntegerField()
    number_meta = models.IntegerField()
    average_length = models.FloatField()
    percentage_id = models.SmallIntegerField()
    average_coverage = models.FloatField()
    change_status = models.TextField()
    seed_consensus = models.TextField()
    full_consensus = models.TextField()
    number_shuffled_hits = models.IntegerField()
    class Meta:
        db_table = u'pfam25_pfama'

class Pfam25PfamSeq(models.Model):
    auto_pfamseq = models.IntegerField(primary_key=True)
    pfamseq_id = models.CharField(max_length=12)
    pfamseq_acc = models.CharField(unique=True, max_length=6)
    seq_version = models.SmallIntegerField()
    crc64 = models.CharField(max_length=16)
    md5 = models.CharField(max_length=32)
    description = models.TextField()
    evidence = models.SmallIntegerField()
    length = models.IntegerField()
    species = models.TextField()
    taxonomy = models.TextField()
    is_fragment = models.SmallIntegerField()
    sequence = models.TextField()
    updated = models.DateTimeField()
    created = models.DateTimeField()
    ncbi_taxid = models.IntegerField()
    genome_seq = models.SmallIntegerField()
    auto_architecture = models.IntegerField()
    treefam_acc = models.CharField(max_length=8)
    class Meta:
        db_table = u'pfam25_pfamseq'

class Pfam25PfamA_RegFullSignificant(models.Model):
    auto_pfama_reg_full = models.IntegerField(primary_key=True)
    pfamA = models.ForeignKey(Pfam25PfamA, db_column='auto_pfama')
    pfamseq = models.ForeignKey(Pfam25PfamSeq, db_column='auto_pfamseq')
    seq_start = models.IntegerField()
    seq_end = models.IntegerField()
    ali_start = models.IntegerField()
    ali_end = models.IntegerField()
    model_start = models.IntegerField()
    model_end = models.IntegerField()
    domain_bits_score = models.FloatField()
    domain_evalue_score = models.CharField(max_length=15)
    sequence_bits_score = models.FloatField()
    sequence_evalue_score = models.CharField(max_length=15)
    cigar = models.TextField()
    in_full = models.SmallIntegerField()
    tree_order = models.IntegerField()
    domain_order = models.SmallIntegerField()
    domain_oder = models.SmallIntegerField()
    class Meta:
        db_table = u'pfam25_pfama_reg_full_significant'

class Treefam(models.Model):
    accession = models.CharField(max_length=40, primary_key=True)
    hmmlength = models.IntegerField()
    family_type = models.CharField(max_length=10)
    class Meta:
        db_table = u'treefam'

class TreeCut(models.Model):
    id = models.IntegerField(primary_key=True)
    tree = models.ForeignKey(Tree)
    method = models.TextField()
    parameter = models.FloatField()
    class Meta:
        db_table = u'tree_cut'

class HMM(models.Model):
    id = models.AutoField(primary_key=True)
    length = models.IntegerField()
    hmm_type = models.CharField(max_length=10)
    method = models.CharField(max_length=255)
    pfam = models.ForeignKey(Pfam, null=True)
    treefam_accession = models.ForeignKey(Treefam,
        db_column='treefam_accession', null=True)
    tigrfam_accession = models.ForeignKey('Tigrfam',
        db_column='tigrfam_accession', null=True)
    tree_node = models.ForeignKey(TreeNode, null=True)
    tree_cut = models.ForeignKey(TreeCut, null=True)
    source_hmm = models.ForeignKey('HMM', null=True)
    class Meta:
        db_table = u'hmm'

class HMM_Consensus(models.Model):
    id = models.AutoField(primary_key=True)
    hmm = models.ForeignKey(HMM)
    sequence = models.ForeignKey(Sequence)
    class Meta:
        db_table = u'hmm_consensus'

class SpeciesNearestMaximalTreeNode(models.Model):
    id = models.AutoField(primary_key=True)
    tree_node = models.ForeignKey(TreeNode)
    species = models.ForeignKey(UniProtTaxonomy, db_column='species')
    sequence_header = models.ForeignKey(SequenceHeader)
    next_nearest_distance = models.FloatField(null=True)
    next_nearest_sequence_header = models.ForeignKey(SequenceHeader,
                          null=True,
                          db_column='next_nearest_sequence_header_id',
                          related_name='sibling_proximal_subtrees')
    distance = models.FloatField(null=True, default='-1')
    is_descendant = models.BooleanField()
    is_unique = models.BooleanField()
    class Meta:
        db_table = u'species_nearest_maximal_tree_node'

class Metropolis(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(unique=True, max_length=255)
    display_name = models.CharField(max_length=255)
    password = models.CharField(max_length=255, null=True)
    class Meta:
        db_table = u'metropolis'

class PDB(models.Model):
    id = models.CharField(primary_key=True, max_length=4)
    resolution = models.FloatField(null=True)
    date_original = models.DateField(null=True)
    class Meta:
        db_table = u'pdb'

    def file_path(self):
        return '/clusterfs/ohana/external/pdb_structures/%s/pdb%s.ent.gz' % (self.id[1:3], self.id)


class PDB_Chain(models.Model):
    id = models.AutoField(primary_key=True)
    pdb = models.ForeignKey(PDB)
    chain_id = models.CharField(max_length=1)
    full_sequence = models.ForeignKey(Sequence, null=True)
    all_residues_have_atom_records_f = models.NullBooleanField()
    description = models.TextField(null=True)
    class Meta:
        db_table = u'pdb_chain'

class SCOP(models.Model):
    id = models.AutoField(primary_key=True)
    class_letter = models.CharField(max_length=4)
    fold_number = models.IntegerField(null=True)
    superfamily_number = models.IntegerField(null=True)
    family_number = models.IntegerField(null=True)
    family = models.ForeignKey(Family, null=True)
    description = models.TextField(null=True)
    class Meta:
        db_table = u'scop'

class PDB_SCOP(models.Model):
    id = models.AutoField(primary_key=True)
    pdb_chain = models.ForeignKey(PDB_Chain)
    scop = models.ForeignKey(SCOP)
    start_residue = models.IntegerField(null=True)
    end_residue = models.IntegerField(null=True)
    scop_version = models.FloatField()
    scop_identifier = models.CharField(max_length=7, unique=True)
    class Meta:
        db_table = u'pdb_scop'

class Replacement(models.Model):
    id = models.AutoField(primary_key=True)
    old_family = models.ForeignKey(Family)
    new_family = models.ForeignKey(Family,
        related_name="replacement_new_family")
    class Meta:
        db_table = u'replacement'

class ResidueScore(models.Model):
    id = models.AutoField(primary_key=True)
    family = models.ForeignKey(Family, null=True)
    residue_number = models.IntegerField()
    score = models.FloatField()
    method = models.CharField(max_length=40)
    class Meta:
        db_table = u'residue_score'

class SequenceHMM(models.Model):
    id = models.AutoField(primary_key=True)
    hmm = models.ForeignKey(HMM)
    sequence = models.ForeignKey(Sequence)
    aligned_sequence = models.ForeignKey(AlignedSequence, null=True)
    bit_score = models.FloatField(null=True)
    e_value = models.FloatField()
    n_aligned_chars = models.IntegerField(null=True)
    n_core_aligned_chars = models.IntegerField(null=True)
    n_windows_below_bl62_cutoff = models.IntegerField(null=True)
    updated = models.DateTimeField(default=datetime.now)
    sequence_type = models.CharField(max_length=255)
    accessed_time = models.DateTimeField(default=datetime.now, null=True)
    hmm_start = models.IntegerField()
    hmm_end = models.IntegerField()
    sequence_start = models.IntegerField()
    sequence_end = models.IntegerField()

    # a dot denotes that the match is not at the domain boundary (this
    # has a jagged edge) and a square bracket indicates that the full
    # HMM was matched at that boundary (depicted as a rounded edge).
    match_type = models.CharField(null=True, max_length=2)

    class Meta:
        db_table = u'sequence_hmm'

class Tigrfam(models.Model):
    id = models.CharField(max_length=255)
    accession = models.CharField(max_length=255, primary_key=True)
    name = models.CharField(max_length=255)
    description = models.TextField()
    global_trusted_cutoff = models.IntegerField()
    frag_trusted_cutoff = models.IntegerField()
    global_noise_cutoff = models.IntegerField()
    frag_noise_cutoff = models.IntegerField()
    isology = models.CharField(max_length=255)
    gene_symbol = models.CharField(max_length=10)
    expanded_name = models.CharField(max_length=255)
    notes = models.TextField()
    ec_id = models.IntegerField()
    class Meta:
        db_table = u'tigrfam'

class TigrfamThirdParty(models.Model):
    tigrfam_acc = models.ForeignKey(Tigrfam, db_column='tigrfam_acc')
    thirdparty = models.ForeignKey(ThirdParty)
    class Meta:
        db_table = u'tigrfam_thirdparty'

class TreeCutTreeNode(models.Model):
    id = models.IntegerField(primary_key=True)
    tree_cut = models.ForeignKey(TreeCut)
    tree_node = models.ForeignKey(TreeNode)
    subfamiy_node_id = models.CharField(max_length=20)
    class Meta:
        db_table = u'tree_cut_tree_node'

class TreeML_Model(models.Model):
    id = models.IntegerField(primary_key=True)
    model = models.CharField(max_length=20)
    tree = models.ForeignKey(Tree)
    class Meta:
        db_table = u'tree_mlmodel'

class TreeNodeLiterature(models.Model):
    id = models.IntegerField(primary_key=True)
    tree_node = models.ForeignKey(TreeNode)
    reference = models.TextField()
    url = models.CharField(max_length=255)
    auth_user = models.ForeignKey(User)
    class Meta:
        db_table = u'tree_node_literature'

class TreeNodeName(models.Model):
    id = models.AutoField(primary_key=True)
    tree_node = models.ForeignKey(TreeNode)
    name = models.TextField()
    source = models.CharField(max_length=10)
    method = models.CharField(max_length=255)
    auth_user_id = models.IntegerField()
    created_at = models.DateTimeField(default=datetime.now)
    class Meta:
        db_table = u'tree_node_name'

class TreeNodeAlignment(models.Model):
    id = models.AutoField(primary_key=True)
    tree_node = models.ForeignKey(TreeNode)
    sequence_header = models.ForeignKey(SequenceHeader)
    aligned_sequence = models.ForeignKey(AlignedSequence)
    class Meta:
        db_table = u'tree_node_alignment'

class TreeNodeAlignmentConservation(models.Model):
    id = models.AutoField(primary_key=True)
    tree_node = models.ForeignKey(TreeNode)
    column_index = models.IntegerField()
    blosum62_conservation_score = models.FloatField()
    conserved_residue = models.CharField(max_length=1)
    class Meta:
        db_table = u'tree_node_alignment_conservation'


class TreeNodePDB(models.Model):
    id = models.AutoField(primary_key=True)
    sequence_hmm = models.ForeignKey(SequenceHMM)
    tree_node = models.ForeignKey(TreeNode)
    pdb_chain = models.ForeignKey(PDB_Chain)
    class Meta:
        db_table = u'tree_node_pdb'

    def __unicode__(self):
        return u'%s' % self.id
    __str__ = __unicode__

class TreeNodeAnnotation(models.Model):
    id = models.IntegerField(primary_key=True)
    tree_node = models.ForeignKey(TreeNode)
    annotation = models.TextField()
    auth_user = models.ForeignKey(User)
    class Meta:
        db_table = u'tree_node_annotation'

class TreeNodeSignalPeptide(models.Model):
    id = models.AutoField(primary_key=True)
    seq_start = models.IntegerField()
    seq_end = models.IntegerField()
    method = models.TextField()
    tree_node = models.ForeignKey(TreeNode)
    class Meta:
        db_table = u'tree_node_signal_peptide'

class TreeNodeTransmembrane(models.Model):
    id = models.AutoField(primary_key=True)
    seq_start = models.IntegerField()
    seq_end = models.IntegerField()
    method = models.TextField()
    tree_node = models.ForeignKey(TreeNode)
    class Meta:
        db_table = u'tree_node_transmembrane'


class TreeNodeSpecies(models.Model):
    id = models.IntegerField(primary_key=True)
    tree_node = models.ForeignKey(TreeNode)
    species = models.ForeignKey(UniProtTaxonomy,
                              db_column = 'species')
    is_maximal = models.IntegerField()
    class Meta:
        db_table = u'tree_node_species'

class EC(models.Model):
    id = models.AutoField(primary_key=True)
    class_number = models.IntegerField()
    subclass_number = models.IntegerField(null=True)
    subsubclass_number = models.IntegerField(null=True)
    enzyme_number = models.IntegerField(null=True)
    description = models.CharField(max_length=255, null=True)
    family = models.ForeignKey(Family, null=True)
    is_preliminary_f = models.NullBooleanField()
    def __unicode__(self):
        if self.enzyme_number:
            if self.is_preliminary_f:
                return u'%d.%d.%d.n%d' % (self.class_number, self.subclass_number,
                                        self.subsubclass_number, self.enzyme_number)
            else:
                return u'%d.%d.%d.%d' % (self.class_number, self.subclass_number,
                                        self.subsubclass_number, self.enzyme_number)
        elif self.subsubclass_number:
            return u'%d.%d.%d.-' % (self.class_number, self.subclass_number,
                                    self.subsubclass_number)
        elif self.subclass_number:
            return u'%d.%d.-.-' % (self.class_number, self.subclass_number)
        else:
            return u'%d.-.-.-' % self.class_number

    # This is slow, I'm sorry
    def as_annotation(self):
        output = []
        if self.class_number:
            ec = EC.objects.get(class_number=self.class_number,subclass_number__isnull=True,subsubclass_number__isnull=True,enzyme_number__isnull=True)
            output.append({ 'value': ec.class_number, 'description': ec.description})
        if self.subclass_number:
            ec = EC.objects.get(class_number=self.class_number,subclass_number=self.subclass_number,subsubclass_number__isnull=True,enzyme_number__isnull=True)
            output.append({ 'value': ec.subclass_number, 'description': ec.description})
        if self.subsubclass_number:
            ec = EC.objects.get(class_number=self.class_number,subclass_number=self.subclass_number,subsubclass_number=self.subsubclass_number,enzyme_number__isnull=True)
            output.append({ 'value': ec.subsubclass_number, 'description': ec.description})
        if self.enzyme_number:
            output.append({ 'value': self.enzyme_number, 'description': self.description})
        return TreeAnnotation(output)

    def __str__(self):
        return self.__unicode__()
    class Meta:
        db_table = u'ec'

class KEGG_Map(models.Model):
    # Note this is *not* an AutoField, it is set by KEGG
    id = models.IntegerField(primary_key=True)
    title = models.TextField()
    def __unicode__(self):
        return u'map%05d: %s' % (self.id, self.title)
    def __str__(self):
        return self.__unicode__()
    class Meta:
        db_table = u'kegg_map'

class KEGG_Map_EC(models.Model):
    id = models.AutoField(primary_key=True)
    kegg_map = models.ForeignKey(KEGG_Map)
    ec = models.ForeignKey(EC)
    class Meta:
        db_table = u'kegg_map_ec'

class UniProtEC(models.Model):
    id = models.AutoField(primary_key=True)
    uniprot = models.ForeignKey(UniProt, related_name='ec_annotations')
    ec = models.ForeignKey(EC)
    is_in_brenda_f = models.BooleanField()
    class Meta:
        db_table = u'uniprot_ec'

class UniProtDatIndex(models.Model):
    id = models.AutoField(primary_key=True)
    file_char = BigIntegerField()
    has_accession_split = models.NullBooleanField(default=False)
    uniprot_accession = models.CharField(max_length=6)
    uniprot = models.ForeignKey('UniProt')
    last_update = models.DateTimeField(auto_now=True)

    def __unicode__(self):
        return u'%s: %s' % (self.uniprot_accession,
            self.file_char)
    __str__ = __unicode__

    class Meta:
        db_table = u'uniprot_dat_index'

class UniProtGenbank(models.Model):
    id = models.AutoField(primary_key=True)
    uniprot_accession = models.CharField(max_length=6, null=True)
    genbank_accession = models.CharField(max_length=255, null=True)
    class Meta:
        db_table = u'uniprot_genbank'

class UniProtGO(models.Model):
    id = models.AutoField(primary_key=True)
    go_term = models.ForeignKey(GO_Term)
    go_evidence = models.ForeignKey(GO_EvidencePriority)
    uniprot = models.ForeignKey(UniProt, related_name='go_annotations')

    def __unicode__(self):
        return (unicode(self.go_term) + "\n" + unicode(self.go_evidence))

    __str__ = __unicode__

    class Meta:
        db_table = u'uniprot_go'

class UniProtLiterature(models.Model):
    id = models.AutoField(primary_key=True)
    uniprot = models.ForeignKey(UniProt, related_name='literature')
    title = models.TextField()
    authors = models.TextField(null=True)
    is_large_scale_f = models.NullBooleanField()
    medline_ui = models.CharField(max_length=9, null=True)
    pmid = models.IntegerField(null=True)
    doi = models.TextField(null=True)
    agricola = models.CharField(max_length=11, null=True)
    class Meta:
        db_table = u'uniprot_literature'

class UniProtThirdPartySequenceHeader(models.Model):
    id = models.AutoField(primary_key=True)
    uniprot = models.ForeignKey(UniProt)
    thirdparty_sequence_header = models.ForeignKey(ThirdPartySequenceHeader)
    class Meta:
        db_table = u'uniprot_thirdparty_sequence_header'

class VillageMetropolis(models.Model):
    metropolis = models.ForeignKey(Metropolis)
    village = models.ForeignKey(Village)
    class Meta:
        db_table = u'village_metropolis'


#############
# New stuff #
#############

class FeatureKey(models.Model):
    id = models.AutoField(primary_key=True)
    key_name = models.CharField(unique=True, max_length=8)
    description = models.TextField(null=True)

    def __unicode__(self):
        return u'%s' % self.key_name
    __str__ = __unicode__

    class Meta:
        db_table = u'feature_key'



class Keyword(models.Model):
    accession = models.CharField(primary_key=True, max_length=7)
    identifier = models.CharField(unique=True, max_length=255)
    definition = models.TextField()
    category = models.ForeignKey('self',
        db_column='category_accession', null=True)

    def __unicode__(self):
        return u'Keyword: %s' % self.accession
    __str__ = __unicode__

    class Meta:
        db_table = u'keyword'

class Keyword_GO_Term(models.Model):
    id = models.AutoField(primary_key=True)
    keyword = models.ForeignKey(Keyword, db_column='keyword_accession')
    go_term = models.ForeignKey(GO_Term)
    class Meta:
        db_table = u'keyword_go_term'

class KeywordSynonym(models.Model):
    id = models.AutoField(primary_key=True)
    keyword = models.ForeignKey(Keyword, db_column='keyword_accession')
    synonym = models.TextField()
    class Meta:
        db_table = u'keyword_synonym'

class KeywordWWW(models.Model):
    id = models.AutoField(primary_key=True)
    keyword = models.ForeignKey(Keyword, db_column='keyword_accession')
    www_site = models.URLField()
    class Meta:
        db_table = u'keyword_www'

class KeywordKeyword(models.Model):
    id = models.AutoField(primary_key=True)
    specific_keyword = models.ForeignKey(Keyword,
        db_column='specific_keyword_accession',
        related_name='keyword_keyword_specific')
    general_keyword = models.ForeignKey(Keyword,
        db_column='general_keyword_accession',
        related_name='keyword_keyword_general')
    class Meta:
        db_table = u'keyword_keyword'

class NonExperimentalQualifier(models.Model):
    id = models.AutoField(primary_key=True)
    description = models.CharField(max_length=40)

    def __unicode__(self):
        return u'%s' % self.description
    __str__ = __unicode__

    class Meta:
        db_table = u'nonexperimental_qualifier'

class Organelle(models.Model):
    id = models.AutoField(primary_key=True)
    description = models.CharField(max_length=255)
    plastid_type = models.CharField(max_length=255, null=True)
    plasmid_name = models.CharField(max_length=255, null=True)
    class Meta:
        db_table = u'organelle'

class PostTranslationalModificationType(models.Model):
    id = models.AutoField(primary_key=True)
    modification = models.CharField(unique=True, max_length=40)
    description = models.TextField()

    def __unicode__(self):
        return u'%s' % self.modification
    __str__ = __unicode__

    class Meta:
        db_table = u'posttranslational_modification_type'

class UniProtFeature(models.Model):
    id = models.AutoField(primary_key=True)
    uniprot = models.ForeignKey(UniProt)
    feature_key = models.ForeignKey(FeatureKey)
    from_residue = models.IntegerField(null=True)
    to_residue = models.IntegerField(null=True)
    from_residue_is_uncertain = models.NullBooleanField()
    to_residue_is_uncertain = models.NullBooleanField()
    extends_n_terminally = models.NullBooleanField()
    extends_c_terminally = models.NullBooleanField()
    description = models.TextField(null=True)
    nonexperimental_qualifier = models.ForeignKey(
        NonExperimentalQualifier, null=True)
    feature_identifier = models.CharField(max_length=14, null=True)
    posttranslational_modification_type = models.ForeignKey(
        PostTranslationalModificationType, null=True)
    dbsnp_rs_number = models.IntegerField(null=True)
    class Meta:
        db_table = u'uniprot_feature'

class UniProtOrganelle(models.Model):
    id = models.AutoField(primary_key=True)
    uniprot = models.ForeignKey('UniProt')
    organelle = models.ForeignKey(Organelle)
    class Meta:
        db_table = u'uniprot_organelle'

class UniProtHostOrganism(models.Model):
    id = models.AutoField(primary_key=True)
    uniprot = models.ForeignKey(UniProt)
    host_organism = models.ForeignKey(UniProtTaxonomy)
    class Meta:
        db_table = u'uniprot_host_organism'

class UniProtKeyword(models.Model):
    id = models.AutoField(primary_key=True)
    uniprot = models.ForeignKey(UniProt)
    keyword = models.ForeignKey(Keyword, db_column='keyword_accession')

    def __unicode__(self):
        return u'(%s, %s)' % (self.uniprot, self.keyword)

    class Meta:
        db_table = u'uniprot_keyword'

class UniProtPDB_Chain(models.Model):
    id = models.AutoField(primary_key=True)
    uniprot = models.ForeignKey(UniProt)
    pdb_chain = models.ForeignKey(PDB_Chain)
    from_residue = models.IntegerField(null=True)
    to_residue = models.IntegerField(null=True)

    class Meta:
        db_table = u'uniprot_pdb_chain'
    

class UniProtPfam(models.Model):
    id = models.AutoField(primary_key=True)
    uniprot = models.ForeignKey(UniProt)
    pfam = models.ForeignKey(Pfam)
    class Meta:
        db_table = u'uniprot_pfam'

class UniProtPfamCoordinates(models.Model):
    id = models.AutoField(primary_key=True)
    uniprot_pfam = models.ForeignKey(UniProtPfam, related_name='coordinates')
    sequence_start = models.IntegerField()
    sequence_end = models.IntegerField()
    def coordinates(self):
        return (self.sequence_start, self.sequence_end)
    
    def __unicode__(self):
        return u'(%s, %s)' % (self.sequence_start, self.sequence_end)

    class Meta:
        db_table = u'uniprot_pfam_coordinates'

class UniProtGeneID(models.Model):
    id = models.AutoField(primary_key=True)
    uniprot = models.ForeignKey(UniProt)
    geneid = models.IntegerField()
    class Meta:
        db_table = u'uniprot_gene_id'

#""" This is a materialized view - it doesn't use the fancy matview stuff that someone wrote though.  We're not really sure where
#  we are going with phogs/edges/trees.  Don't ever do any updates to this table."""
class Phog(models.Model):
    phog_accession = models.TextField(primary_key=True)
    tree_node = models.ForeignKey(TreeNode)
    tree = models.ForeignKey(Tree)
    left_id = models.IntegerField()
    right_id = models.IntegerField()
    duplication_distance = models.FloatField()
    greatest_duplication_distance_of_maximal_descendant = models.FloatField()
    family = models.ForeignKey(Family)
    is_superorthologous = models.BooleanField()

    class Meta:
        db_table = u'phog'
#
#'''The tables below are inherited in the order listed.  They provide the ability to add arbitrary edges between nodes in our tree node table
#One possible end goal would be a 'closure tree', where all possible connections are represented.  Kimmen isn't sold on that idea yes, so right
#now it might be a bit of overkill.  There are some SPs to help with their creation.  Do \df on the psql prompt to see what functions are installed
#
#Also, the related names are poorly selected.  Plese fix them.
#'''
#class TreeEdge(models.Model):
#    id = models.AutoField(primary_key=True)
#    ancestor_tree_node = models.ForeignKey(TreeNode, related_name='edge_descendants')
#    descendant_tree_node = models.ForeignKey(TreeNode, related_name='edge_ancestors')
#    family = models.ForeignKey(Family)
#    tree = models.ForeignKey(Tree)
#
#    class Meta:
#        db_table = u'tree_edge'
#
#
#''' Subclass of TreeEdge where the descendant is a leaf '''
#''' It's very possible that this can be done with model inheritence, but django has burnt me enough times that I'll look into that later'''
#class TreeLeafEdge(models.Model):
#    id = models.AutoField(primary_key=True)
#    ancestor_tree_node = models.ForeignKey(TreeNode, related_name='leaves')
#    descendant_tree_node = models.ForeignKey(TreeNode, related_name='leaf_ancestors')
#    family = models.ForeignKey(Family)
#    tree = models.ForeignKey(Tree)
#    sequence_header = models.ForeignKey(SequenceHeader)
#    uniprot = models.ForeignKey(UniProt)
#    taxon= models.ForeignKey(UniProtTaxonomy)
#
#    class Meta:
#        db_table = u'tree_leaf_edge'
#
#''' Special subclass of TreeEdge for members of PHOGs '''
#class PHOGMembership(models.Model):
#    id = models.AutoField(primary_key=True)
#    ancestor_tree_node = models.ForeignKey(TreeNode, related_name='phog_members' )
#    descendant_tree_node = models.ForeignKey(TreeNode, related_name='phog_membership')
#    family = models.ForeignKey(Family)
#    tree = models.ForeignKey(Tree)
#    sequence_header = models.ForeignKey(SequenceHeader)
#    uniprot = models.ForeignKey(UniProt)
#    taxon= models.ForeignKey(UniProtTaxonomy)
#    phog_accession = models.TextField()
#    left_id = models.IntegerField() # This is actually part of a FKey.  Django can't handle complex Fkeys
#    phog_duplication_distance = models.FloatField()
#    phog_greatest_duplication_distance_of_maximal_descendant = models.FloatField()
#
#    class Meta:
#        db_table = u'phog_membership'




''' FatCat Queue Tables - If you need to, you can drop referential integrity on these, we will probably never use it
    These are here to store fatcat jobs, and maintain the information.  The scripts dont talk to them directly,
    but instead POST to the API.  If you change the API, you can change these.
'''

class OhanaJobs(models.Model):
    id = models.AutoField(primary_key = True)
    submitting_user = models.ForeignKey(User)
    submitting_ip_address = models.TextField()
    job_create_time = models.DateTimeField(default=datetime.now)
    job_type = models.TextField()

    class Meta:
        db_table = u'ohana_job_table'

class FatcatJobStatus(models.Model):
    id = models.AutoField(primary_key=True)
    status = models.TextField()

    class Meta:
        db_table = u'fatcat_job_status'

class FatcatJob(models.Model):
    id = models.AutoField(primary_key=True)
    fasta_header = models.TextField()
    fasta_sequence = models.TextField()
    created_at = models.DateTimeField(default=datetime.now)
    updated_at = models.DateTimeField(default=datetime.now)
    status = models.ForeignKey(FatcatJobStatus, related_name="jobs")
    get_families_pbs_job_id = models.TextField()
    get_best_nodes_pbs_job_id = models.TextField()
    user_email = models.TextField()
    criteria_get_families_e_value = models.FloatField()
    criteria_get_families_mda_query_coverage = models.FloatField()
    criteria_get_families_pfam_query_coverage = models.FloatField()
    criteria_get_families_mda_hmm_coverage = models.FloatField()
    criteria_get_families_pfam_hmm_coverage = models.FloatField()
    criteria_get_best_nodes_e_value = models.FloatField()
    criteria_get_best_nodes_mda_query_coverage = models.FloatField()
    criteria_get_best_nodes_pfam_query_coverage = models.FloatField()
    criteria_get_best_nodes_mda_hmm_coverage = models.FloatField()
    criteria_get_best_nodes_pfam_hmm_coverage = models.FloatField()
    criteria_get_best_nodes_query_to_subtree_hmm_pid = models.FloatField()
    criteria_get_best_nodes_subtree_min_pid = models.FloatField()
    criteria_enclosing_clade_orthology_methods = models.TextField()
    criteria_enclosing_clade_kerf_threshold = models.IntegerField()
    criteria_enclosing_clade_orthology_methods_required = models.IntegerField()
    criteria_cluster_similarity = models.IntegerField()
    criteria_ortholog_coverage = models.IntegerField()
    criteria_query_coverage = models.IntegerField()
    criteria_minimum_pid_to_query_for_orthology = models.IntegerField()
    criteria_consensus_uniprot_base_parameter = models.FloatField()
    criteria_consensus_uniprot_threshold_for_high_confidence = models.FloatField()
    criteria_consensus_uniprot_threshold_for_medium_confidence = models.FloatField()
    submitter_ip_address = models.TextField()
    submitter_user_account = models.ForeignKey(User)
    fatcat_job_name = models.TextField()
    email_subject = models.TextField()
    ohana_job = models.ForeignKey(OhanaJobs, related_name='fatcat_job')
 
    class Meta:
        db_table = u'fatcat_job'

    @property
    def is_deleted(self):
        # returns true if the response file is not where it should be
        return ((self.status_id == 9) and (not os.path.exists('/clusterfs/ohana/software/webserver/temp/fatcat/%d/job%d.httpResponse' % (self.id, self.id))))
        
    @property
    def is_running(self):
        qs = eval(OhanaQueueStatus.objects.get(id=1).status)
        if 'web' in qs:
            if (self.get_best_nodes_pbs_job_id):
                if self.get_best_nodes_pbs_job_id.split('.')[0] in qs['web']:
                    return True
                else:
                    return False
            elif (self.get_families_pbs_job_id):
                if self.get_families_pbs_job_id.split('.')[0] in qs['web']:
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False
    
    @property
    def is_error(self):
        return ((self.status_id < 9) and (not self.is_running))
 
    @property
    def fasta(self):
        return '>' + self.fasta_header + "\n" + self.fasta_sequence

    @property
    def stage_one_mda_query_coverage_criteria(self):
        return round(float(100*self.criteria_get_families_mda_query_coverage), 2)
    
    @property
    def stage_one_pfam_query_coverage_criteria(self):
        return round(float(100*self.criteria_get_families_pfam_query_coverage), 2)

    @property
    def stage_one_mda_hmm_coverage_criteria(self):
        return round(float(100*self.criteria_get_families_mda_hmm_coverage), 2)

    @property
    def stage_one_pfam_hmm_coverage_criteria(self):
        return round(float(100*self.criteria_get_families_pfam_hmm_coverage), 2)
    
    @property
    def stage_two_mda_query_coverage_criteria(self):
        return round(float(100*self.criteria_get_best_nodes_mda_query_coverage), 2)
    
    @property
    def stage_two_pfam_query_coverage_criteria(self):
        return round(float(100*self.criteria_get_best_nodes_pfam_query_coverage), 2)

    @property
    def stage_two_mda_hmm_coverage_criteria(self):
        return round(float(100*self.criteria_get_best_nodes_mda_hmm_coverage), 2)

    @property
    def stage_two_pfam_hmm_coverage_criteria(self):
        return round(float(100*self.criteria_get_best_nodes_pfam_hmm_coverage), 2)
    
    def has_informative_families(self):
        return len(self.informative_families()) > 0

    def has_ambiguous_families(self):
        return len(self.ambiguous_families()) > 0

    def informative_families(self):
        return [ f for f in FatcatJobFamily.objects.filter(fatcat_job=self, best_tree_node__isnull=False).order_by('best_node_e_value').all() if f.informative ]

    def ambiguous_families(self):
        return [ f for f in FatcatJobFamily.objects.filter(fatcat_job=self, best_tree_node__isnull=False).order_by('best_node_e_value').all() if not f.informative ]

    def all_families(self):
        return [ f for f in FatcatJobFamily.objects.filter(fatcat_job=self, family_e_value__isnull=False).order_by('family_e_value').all() ]

    def paralogs(self):
        try:
            return self._paralogs
        except:
            pass
        (o, p) = self.get_orthologs_and_paralogs()
        self._paralogs = p
        self._orthologs = o
        return self._paralogs

    def orthologs(self):
        try:
            return self._orthologs
        except:
            pass
        (o,p) = self.get_orthologs_and_paralogs()
        self._paralogs = p
        self._orthologs = o
        return self._orthologs

    def get_orthologs_and_paralogs(self):
        orthologs = {}
        orthology_supported_ancestors = {}

        def _traverse_tree(node, node_orthology_support):
            """ Recursively traverses this subtree, builds orthologs list, and does the 
            orthology support stuff... """
            for (system, root) in orthology_supported_ancestors.items():
                if root == node:
                    node_orthology_support[system] = root
            if not (node.is_leaf()):
                for child in node.children():
                    _traverse_tree(child, node_orthology_support)
            else:
                # build the list of orthologs
                if ((node.sequence_header.uniprot) and (node.sequence_header.uniprot.fasta_from_db) and (self.fasta)):
                    u = node.sequence_header.uniprot
                    u.orthology_support = {
                        'kerf70': node_orthology_support['kerf_supported_ancestor'], 
                        'oma': node_orthology_support['oma_supported_ancestor'], 
                        'orthomcl': node_orthology_support['orthomcl_supported_ancestor'], 
                        'phogt0': node_orthology_support['phog_supported_ancestor']
                    }
                    # do alignment of candidate ortholog to query
                    ali = AlignIO.read(StringIO.StringIO(mafft(u.fasta_from_db+'\n'+self.fasta)), 'fasta')
                    ortholog_alignment = str(ali[0].seq)
                    query_alignment = str(ali[1].seq)
                    u.pwid_to_query = 100*pwid_belvu(ortholog_alignment, query_alignment)
                    u.alignment = repr('>' + str(ali[0].id) + '\n' + ortholog_alignment + '\n>' + str(ali[1].id) + '\n' + query_alignment)
                    # do overlap of query and ortholog
                    query_length = len(self.fasta_sequence)
                    ortholog_length = u.sequence_len
                    aligned_chars = 0
                    for (index, char) in enumerate(ortholog_alignment):
                        if (char.isupper() and query_alignment[index].isupper()):
                            aligned_chars += 1
                    u.query_coverage = float(100.0*aligned_chars/query_length)
                    u.ortholog_coverage = float(100.0*aligned_chars/ortholog_length)
                    u.family = family.family
                    orthologs.setdefault(u.taxon.id,[]).append(u)
            return 

        for family in self.informative_families():
            informative_node = family.best_tree_node.get_informative_node(fatcat_family=family)
            orthology_supported_ancestors = family.best_tree_node.orthology_supported_ancestors(kerf_threshold=self.criteria_enclosing_clade_kerf_threshold)
            node_orthology_support_default = {
                        'oma_supported_ancestor':False, 
                        'orthomcl_supported_ancestor':False, 
                        'kerf_supported_ancestor':False, 
                        'phog_supported_ancestor':False
            }
            _traverse_tree(informative_node, node_orthology_support_default)    
        (o,p) = fatcat_ortholog_gene_clustering(orthologs, self.criteria_cluster_similarity, self.criteria_query_coverage, self.criteria_ortholog_coverage, self.criteria_minimum_pid_to_query_for_orthology)
        # write these to a file.
        file = open('/clusterfs/ohana/software/webserver/temp/fatcat/%d/Job_%d_orthologs.csv' % (self.id, self.id),'w')
        file.write('# UniProt accession, Taxon name, Taxon ID, UniProt description, % ID to query, Query coverage, Ortholog coverage\n')
        for ortholog_cluster in o:
            for uniprot in ortholog_cluster:
                file.write('"%s","%s",%d,"%s",%.2f,%.2f,%.2f\n' % (uniprot[0].accession, uniprot[0].taxon.scientific_name, uniprot[0].taxon.id, uniprot[0].description, uniprot[0].pwid_to_query, uniprot[0].query_coverage, uniprot[0].ortholog_coverage))
        file.close()                
        return (o,p)

#    def orthologs_with_distances(self):
    def _old_orthologs(self):
        # This will become orthologs once it is done
        # This code is crap.
        perspective_orthologs = dict()
        for family in self.informative_families():
            for leaf in family.best_tree_node.get_informative_node(fatcat_family=family).get_informative_members():
                if not perspective_orthologs.has_key(leaf.sequence_header.taxon_id):
                    perspective_orthologs[leaf.sequence_header.taxon_id] = []
                perspective_orthologs[leaf.sequence_header.taxon_id].append({
                    'leaf': leaf,
                    'fatcat_hit': family.best_tree_node
                })

        orthologs = []
        for taxon_id, po_list in perspective_orthologs.items():
            distance = float("inf")
            o = None
            for po in po_list:
                if po['leaf'].left_id > po['fatcat_hit'].left_id and po['leaf'].right_id < po['fatcat_hit'].right_id:
                    po_distance = po['fatcat_hit'].average_distance_down() + \
                        ( TreeNode.objects.filter(
                                tree=po['fatcat_hit'].tree, 
                                left_id__gt=po['fatcat_hit'].left_id, 
                                left_id__lte=po['leaf'].left_id, 
                                right_id__lt=po['fatcat_hit'].right_id, 
                                right_id__gte=po['leaf'].right_id, 
                        ).aggregate(Sum('distance_to_parent'))['distance_to_parent__sum'] or 0)
                else:
                    po_distance = po['fatcat_hit'].average_distance_down()
                    mrca = TreeNode.objects.filter(
                                tree=po['fatcat_hit'].tree, 
                                left_id__lte=min(po['fatcat_hit'].left_id, po['leaf'].left_id),
                                right_id__gte=max(po['fatcat_hit'].right_id, po['leaf'].right_id),
                            ).order_by('-left_id')[0]
                    #distance between fatcat hit and most recent common ancestor
                    po_distance += \
                        ( TreeNode.objects.filter(
                                tree=po['fatcat_hit'].tree, 
                                left_id__lte=po['fatcat_hit'].left_id, 
                                left_id__gt=mrca.left_id, 
                                right_id__gte=po['fatcat_hit'].right_id, 
                                right_id__lt=mrca.right_id, 
                        ).aggregate(Sum('distance_to_parent'))['distance_to_parent__sum'] or 0)
                    # distance betweeen candidate ortholog and mcra
                    po_distance += \
                        ( TreeNode.objects.filter(
                                tree=po['leaf'].tree, 
                                left_id__lte=po['leaf'].left_id, 
                                left_id__gt=mrca.left_id, 
                                right_id__gte=po['leaf'].right_id, 
                                right_id__lt=mrca.right_id, 
                        ).aggregate(Sum('distance_to_parent'))['distance_to_parent__sum'] or 0)

                if po_distance < distance:
                    distance = po_distance
                    o = po['leaf']
            if o:
                u = o.sequence_header.uniprot
                u.distance = po_distance
                u.pwid_to_query, u.alignment = o.alignment_to_query(self.fasta)
                u.pwid_to_query = str(u.pwid_to_query)
                u.alignment = repr(u.alignment)
                orthologs.append(u)
                
        return orthologs

    def summary(self):
        try:
            return self._summary
        except:
            pass
        summary = { 
            "ec": [], 
            "go_molecular_function": [], 
            "go_biological_process": [], 
            "go_cellular_component": [],
            "uniprot_descriptions": [] 
        }
        o = self.orthologs()
        uniprots = set()
        rep_uniprots = set()
        # get list of uniprots of orthologs
        for uniprot_tuple in o:
            for uni in uniprot_tuple:
                rep_uniprots.add(uni[0])
        rep_uniprots = list(rep_uniprots)
        
        go_molecular = {}
        go_biological = {}
        go_cellular = {}

        for uni in rep_uniprots:
            # get molecular function
            for annotation in Annotation.objects.filter(uniprot__accession = uni.accession,
                ontology_term__ontologysubsetmembership__subset = "molecular_function").values(
                'ontology_term__accession','ontology_term__name','evidence__code','evidence__name',
                'evidence__priority'):
                t = ({'accession': annotation['ontology_term__accession'], 
                     'description': annotation['ontology_term__name'],
                     'evidence_code': annotation['evidence__code'],
                     'evidence_priority': annotation['evidence__priority'],
                     'evidence_description': annotation['evidence__name']}, uni)
                if t[0]['accession'] in go_molecular:
                    if t in go_molecular[t[0]['accession']]:
                        pass
                    else:
                        go_molecular[t[0]['accession']].append(t)
                else:
                    go_molecular[t[0]['accession']] = [t]
            # get biological process
            for annotation in Annotation.objects.filter(uniprot__accession = uni.accession,
                ontology_term__ontologysubsetmembership__subset = "biological_process").values(
                'ontology_term__accession','ontology_term__name','evidence__code','evidence__name',
                'evidence__priority'):
                t = ({'accession': annotation['ontology_term__accession'], 
                     'description': annotation['ontology_term__name'],
                     'evidence_code': annotation['evidence__code'],
                     'evidence_priority': annotation['evidence__priority'],
                     'evidence_description': annotation['evidence__name']}, uni)
                if t[0]['accession'] in go_biological:
                    if t in go_biological[t[0]['accession']]:
                        pass
                    else:
                        go_biological[t[0]['accession']].append(t)
                else:
                    go_biological[t[0]['accession']] = [t]
            # go cellular component
            for annotation in Annotation.objects.filter(uniprot__accession = uni.accession,
                ontology_term__ontologysubsetmembership__subset = "cellular_component").values(
                'ontology_term__accession','ontology_term__name','evidence__code','evidence__name',
                'evidence__priority'):
                t = ({'accession': annotation['ontology_term__accession'], 
                     'description': annotation['ontology_term__name'],
                     'evidence_code': annotation['evidence__code'],
                     'evidence_priority': annotation['evidence__priority'],
                     'evidence_description': annotation['evidence__name']}, uni)
                if t[0]['accession'] in go_cellular:
                    if t in go_cellular[t[0]['accession']]:
                        pass
                    else:
                        go_cellular[t[0]['accession']].append(t)
                else:
                    go_cellular[t[0]['accession']] = [t]
        for (go_acc, annotation_tuple_list) in go_molecular.items():
            go_molecular[go_acc] = sorted(annotation_tuple_list, key=lambda t: t[0]['evidence_priority'])
        for (go_acc, annotation_tuple_list) in go_biological.items():
            go_biological[go_acc] = sorted(annotation_tuple_list, key=lambda t: t[0]['evidence_priority'])
        for (go_acc, annotation_tuple_list) in go_cellular.items():
            go_cellular[go_acc] = sorted(annotation_tuple_list, key=lambda t: t[0]['evidence_priority'])
        if go_molecular:
            summary["go_molecular_function"] = sorted(go_molecular.items(), key=lambda t: t[1][0][0]['evidence_priority'])
        if go_biological:
            summary["go_biological_process"] = sorted(go_biological.items(), key=lambda t: t[1][0][0]['evidence_priority'])
        if go_cellular:
            summary["go_cellular_component"] = sorted(go_cellular.items(), key=lambda t: t[1][0][0]['evidence_priority'])
        # do ec
        ecd = {}
        for ec in EC.objects.filter(uniprotec__uniprot__in = uniprots):
            if not ec.is_preliminary_f:
                if ec in ecd:
                    ecd[ec] += 1
                else:
                    ecd[ec] = 1
        if ecd:
            summary["ec"] = sorted(ecd.items(), key=lambda t: t[1], reverse=True)
        # do consensus description
        summary['uniprot_descriptions'] = fatcat_ortholog_consensus_uniprot_description(rep_uniprots, fatcat_query=self.fasta, base=self.criteria_consensus_uniprot_base_parameter)
        self._summary = summary
        return summary

    def results(self):
        return [ f.as_result() for f in FatcatJobFamily.objects.filter(
                fatcat_job=self, best_tree_node__isnull=False).order_by('best_node_e_value') ]

class FatcatJobPfams(models.Model):
    id = models.AutoField(primary_key=True)
    fatcat_job = models.ForeignKey(FatcatJob, related_name="pfams")
    pfam_description = models.TextField()
    pfam_accession = models.TextField()
    pfam_shortname = models.TextField()
    ali_from = models.IntegerField()
    ali_to = models.IntegerField()
    hmm_from = models.IntegerField()
    hmm_to = models.IntegerField()
    c_evalue = models.FloatField()
    i_evalue = models.FloatField()

    class Meta:
        db_table = u'fatcat_job_pfams'

class FatcatJobFamily(models.Model):
    id = models.AutoField(primary_key=True)
    fatcat_job = models.ForeignKey(FatcatJob, related_name="families")
    family = models.ForeignKey(Family)
    family_type = models.TextField()
    best_tree_node = models.ForeignKey(TreeNode)
    family_e_value = models.FloatField()
    best_node_e_value = models.FloatField()

    @property
    def informative(self):
        try:
            return self._informative
        except:
            pass

        if self.incomplete:
            return None
        self._informative = self.best_tree_node.get_informative_node(fatcat_family=self) 
        return self._informative

    @property
    def incomplete(self):
        try:
            j = self.best_tree_node
        except:
            return True
        return False

    @property
    def hmm_coverage(self):
        fatcat_match = FatcatJobFamilyMatchCoord.objects.get(fatcat_job_family = self)
        return float(fatcat_match.best_hmm_to - fatcat_match.best_hmm_from)/fatcat_match.best_hmm_length

    @property
    def query_coverage(self):
        fatcat_match = FatcatJobFamilyMatchCoord.objects.get(fatcat_job_family = self)
        return float(fatcat_match.best_ali_to - fatcat_match.best_ali_from)/len(self.fatcat_job.fasta_sequence)        

    @property
    def family_hmm_coverage(self):
        fatcat_match = FatcatJobFamilyMatchCoord.objects.get(fatcat_job_family = self)
        return float(fatcat_match.hmm_to - fatcat_match.hmm_from)/fatcat_match.hmm_length

    @property
    def readable_family_hmm_coverage(self):
        return 100*self.family_hmm_coverage

    @property
    def family_query_coverage(self):
        fatcat_match = FatcatJobFamilyMatchCoord.objects.get(fatcat_job_family = self)
        return float(fatcat_match.ali_to - fatcat_match.ali_from)/len(self.fatcat_job.fasta_sequence)

    @property
    def readable_family_query_coverage(self):
        return 100 * self.family_query_coverage

    @property
    def passes_first_stage_coverage_conditions(self):
        if self.family_type == "Pfam":
            return ((self.family_hmm_coverage >= self.fatcat_job.criteria_get_families_pfam_hmm_coverage))
        else:
            return ((self.family_hmm_coverage >= self.fatcat_job.criteria_get_families_mda_hmm_coverage) and
                    (self.family_query_coverage >= self.fatcat_job.criteria_get_families_mda_hmm_coverage))

    @property
    def passes_second_stage_coverage_conditions(self):
        if self.family_type == "Pfam":
            return ((self.hmm_coverage >= self.fatcat_job.criteria_get_best_nodes_pfam_hmm_coverage)) 
        else:
            return ((self.hmm_coverage >= self.fatcat_job.criteria_get_best_nodes_mda_hmm_coverage) and
                    (self.query_coverage >= self.fatcat_job.criteria_get_best_nodes_mda_query_coverage))
    
    def get_family_consensus_pfam_mda(self):
        # for a pfam family, return the pfam domain name.  for a ghg return the consensus architecture
        # if there is one, else return an empty string.
        if self.family.family_type.id == 'C':
            return self.family.get_pfams().pop()[0].name
        else:
            try:
                return self.family.canonical_root_node().get_pfam_domain_architecture_in_members()[0][0]
            except:
                return ""

    def get_family_taxon(self):
        try:
            return self.family.get_taxon()
        except:
            return {'id':0, 'scientific_name': 'N/A', 'common_name': None}

    def get_family_match_positions(self):
        return [ { "alignment_match_from": c.ali_from, "alignment_match_to": c.ali_to,
                                "hmm_match_from": c.hmm_from, "hmm_match_to": c.hmm_to } 
                                    for c in self.coords.all()  ]

    # TODO memoize wrapper
    def as_result(self):
        try:
            return self._as_result
        except:
            pass

        if self.informative:
            tree_node = self.informative
            informative = True
        else:
            tree_node = self.best_tree_node
            informative = False
        query_to_hmm_pid, alignment = self.best_tree_node.alignment_to_query(self.fatcat_job.fasta)
        alignment = repr(alignment)
        if self.family.canonical_root_node().get_taxon().common_name:
            family_taxon_common_name = self.family.canonical_root_node().get_taxon().common_name
        else:
            family_taxon_common_name = ""
        try:
            pfam_mda = self.best_tree_node.get_pfam_domain_architecture_in_members()[0][0]
        except:
            pfam_mda = ""
        try:
            family_pfam_mda = self.family.canonical_root_node().get_pfam_domain_architecture_in_members()[0][0]
        except:
            family_pfam_mda = ""
        if self.family_type == "Pfam":
            family_pfam_mda = self.family.get_pfams().pop()[0].name
        try:
            family_taxonomic_distribution = self.family.canonical_root_node().get_taxon().scientific_name
            family_taxonomic_distribution_id = self.family.canonical_root_node().get_taxon().id
        except:
            family_taxonomic_distribution = "N/A"
            family_taxonomic_distribution_id = 0
        try:
            node_taxonomic_distribution = self.best_tree_node.get_taxon().scientific_name
            node_taxonomic_distribution_id = self.best_tree_node.get_taxon().id
        except:
            node_taxonomic_distribution = "N/A"
            node_taxonomic_distribution_id = 0
        self._as_result = {
            "informative": informative,
            "family_id": self.family_id,
            "family_name": self.family.canonical_root_node().get_treenode_names()['treenode_name'],
            "family_type": self.family_type, 
            "node_id": self.best_tree_node_id,
            "node": tree_node,
            "query_to_hmm_pid": query_to_hmm_pid,
            "query_coverage": 100*self.query_coverage,
            "hmm_coverage": 100*self.hmm_coverage,
            "family_query_coverage": 100*self.family_query_coverage,
            "family_taxonomic_distribution": family_taxonomic_distribution,
            "family_taxonomic_distribution_common_name": family_taxon_common_name,
            "family_taxonomic_distribution_taxon_id": family_taxonomic_distribution_id,
            "family_hmm_coverage": 100*self.family_hmm_coverage,
            "family_evalue": self.family_e_value,
            "family_pfam_mda": family_pfam_mda,
            "subtree_min_pairwise_id": self.best_tree_node.minimum_pairwise_identity_belvu,
            "alignment":  alignment,            
            "node_description": tree_node.get_description(),
            "node_e_value": self.best_node_e_value,
            "most_common_mda_string": pfam_mda,
            "number_of_taxa": tree_node.get_num_contained_species(),
            "node_taxonomic_distribution": node_taxonomic_distribution,
            "node_taxon_id": node_taxonomic_distribution_id,
            "ec": [],#[ a.annotation for a in summarize(tree_node.annotations('ec')) ],
            "go_molecular_function": [],#[ a.annotation for a in summarize(tree_node.annotations('go_molecular_function') )],
            "go_biological_process": [],#[ a.annotation for a in summarize(tree_node.annotations('go_biological_process') )], 
            "go_cellular_component": [],#[ a.annotation for a in summarize(tree_node.annotations('go_cellular_component') )],
            "biocyc_pathways": [],
            "biocyc_reactions": [],
            "sfld" : [],#[ a.annotation for a in summarize(tree_node.annotations('sfld')) ],
            "match_positions": [ { "alignment_match_from": c.ali_from, "alignment_match_to": c.ali_to,
                                "best_alignment_from": c.best_ali_from, "best_alignment_to": c.best_ali_to,
                                "hmm_match_from": c.hmm_from, "hmm_match_to": c.hmm_to } 
                                    for c in self.coords.all()  ]
        }
        return self._as_result

    class Meta:
        db_table = u'fatcat_job_family'

class FatcatJobFamilyMatchCoord(models.Model):
    id = models.AutoField(primary_key=True)
    fatcat_job_family = models.ForeignKey(FatcatJobFamily, related_name="coords")
    ali_from = models.IntegerField()
    ali_to = models.IntegerField()
    hmm_from = models.IntegerField()
    hmm_to = models.IntegerField()
    best_hmm_from = models.IntegerField()
    best_hmm_to = models.IntegerField()
    best_ali_from = models.IntegerField()
    best_ali_to = models.IntegerField()
    hmm_length = models.IntegerField()
    best_hmm_length = models.IntegerField()

    class Meta:
        db_table = u'fatcat_job_family_match_coord'

# If we want this to really scale, use rabbit or something
# only one daemon is popping stuff from this
class FatcatHMMBuildQueueStatus(models.Model):
    id = models.AutoField(primary_key=True)
    status = models.TextField()

    class Meta:
        db_table = u'fatcat_hmm_build_queue_status'

class FatcatHMMBuildQueue(models.Model):
    id = models.AutoField(primary_key=True)
    family = models.ForeignKey(Family)
    created_at = models.DateTimeField(default=datetime.now)
    count = models.IntegerField(null=False, default=1)
    status = models.ForeignKey(FatcatHMMBuildQueueStatus, related_name="jobs", default=1)
    pbs_job_id = models.TextField()

    class Meta:
        db_table = u'fatcat_hmm_build_queue'

    @classmethod
    def add_family(cls, family):
        try:
            # The unique may cause this to error
            cls.objects.create(family=family)
        except:
            cls.objects.filter(family=family).update(count=F('count') + 1)

    @classmethod
    def pop(cls, status, **update_parameters):
        try:
            row = cls.objects.filter(status=status).order_by('-count', 'created_at')[0]
        except:
            return None

        rows_updated = cls.objects.filter(id=row.id, status=row.status).update(status=F('status') + 1, **update_parameters)
        if not rows_updated:
            return None
        return row

    @classmethod
    def estimated_workers(cls):
        return cls.objects.filter(status__in=[2,3]).count()


class FxnSitePredictionJobStatus(models.Model):
    id = models.AutoField(primary_key=True)
    status = models.TextField()

    class Meta:
        db_table = u'fxn_site_prediction_job_status'

class FxnSitePredictionJob(models.Model):
    id = models.AutoField(primary_key=True)
    fasta_header = models.TextField()
    fasta_sequence = models.TextField()
    created_at = models.DateTimeField(default=datetime.now)
    status = models.ForeignKey(FxnSitePredictionJobStatus, related_name="jobs")
    pbs_job_id = models.TextField()
    jackhmmer_iterations = models.IntegerField()
    user_email = models.TextField()
    jackhmmer_evalue = models.FloatField()
    treecut_pid = models.FloatField()
    num_homologs = models.IntegerField()

    @property
    def fasta(self):
        return '>' + self.fasta_header + "\n" + self.fasta_sequence

    class Meta:
        db_table = u'fxn_site_prediction_job'

class SFLDSuperfamily(models.Model):
    superfamily_id            = models.AutoField(primary_key=True)
    superfamily_name          = models.TextField()
    superfamily_evidence_code = models.TextField()

    class Meta:
        db_table = u'sfld\".\"sfld_superfamily'

class SFLDSubgroup(models.Model):
    subgroup_id    = models.AutoField(primary_key=True)
    subgroup_name  = models.TextField()
    superfamily = models.ForeignKey(SFLDSuperfamily, related_name="subgroups")

    class Meta:
        db_table = u'sfld\".\"sfld_subgroup'

class SFLDReaction(models.Model):
    reaction_id   = models.AutoField(primary_key=True)
    reaction_name = models.TextField()
    ec_number     = models.TextField()

    class Meta:
        db_table = u'sfld\".\"sfld_reaction'

class SFLDFamily(models.Model):
    family_id            = models.AutoField(primary_key=True)
    family_name          = models.TextField()
    family_evidence_code = models.TextField()
    subgroup          = models.ForeignKey(SFLDSubgroup, related_name="families")
    superfamily       = models.ForeignKey(SFLDSuperfamily, related_name="families")
    reaction          = models.ForeignKey(SFLDReaction, related_name="families")

    def as_annotation(self):
        return TreeAnnotation([
            { 'value': self.superfamily.superfamily_id, 'name': self.superfamily.superfamily_name },
            { 'value': self.subgroup.subgroup_id, 'name': self.subgroup.subgroup_name },
            { 'value': self.family_id, 'name': self.family_name }, #TODO reaction
        ])

    class Meta:
        db_table = u'sfld\".\"sfld_family'

class SFLDEFD(models.Model):
    efd_id         = models.AutoField(primary_key=True)
    efd_name       = models.TextField()
    family      = models.ForeignKey(SFLDFamily, related_name="efds")
    subgroup    = models.ForeignKey(SFLDSubgroup, related_name="efds")
    superfamily = models.ForeignKey(SFLDSuperfamily, related_name="efds")
    family_evidence_code = models.TextField()
    superfamily_evidence_code = models.TextField()
    
    def in_phylofacts_family(self):
        return TreeNode.objects.filter(
            sequence_header__uniprot__sfldtopfacts__efd = self,
            tree__family__status = 'draft',
            tree__family__active = True
        ).count() > 0

    class Meta:
        db_table = u'sfld\".\"sfld_efd'

class SFLDFasta(models.Model):
    efd         = models.ForeignKey(SFLDEFD, primary_key=True)
    seguid         = models.TextField()
    fasta_sequence = models.TextField()
    
    class Meta:
        db_table =u'sfld\".\"sfld_fasta'


#I'm just a view, don't write to me
class SFLDtoPfacts(models.Model):
    id         = models.AutoField(primary_key=True)
    efd     = models.ForeignKey(SFLDEFD)
    uniprot = models.ForeignKey(UniProt)

    class Meta:
        db_table = u'sfld_to_pfacts'

class SFLDTaxonomy(models.Model):
    id = models.IntegerField(primary_key = True)
    efd_id = models.IntegerField()
    taxonomy_id = models.TextField()

    class Meta:
        db_table = u'sfld\".\"sfld_efd_taxonomy_id'

class OrthoMclGroupUniProt(models.Model):
    id = models.AutoField(primary_key = True)
    orthomcl_group_id = models.TextField()
    uniprot = models.ForeignKey(UniProt, related_name="orthomcl_group")

    def __unicode__(self):
        return unicode(self.orthomcl_group_id)

    class Meta:
        db_table = u'orthomcl\".\"orthomcl_group_uniprot'

class OMAGroupUniProt(models.Model):
    id = models.AutoField(primary_key = True)
    oma_group_id = models.TextField()
    uniprot = models.ForeignKey(UniProt, related_name="oma_group")

    def __unicode__(self):
        return unicode(self.oma_group_id)

    class Meta:
        db_table = u'oma\".\"oma_group_uniprot'

class Ontology(models.Model): 
    name = models.TextField(primary_key = True)
    date_updated = models.DateField(default=datetime.now())

    class Meta:
        db_table = u'ontology'

class OntologyTerm(models.Model):
    accession = models.TextField(primary_key = True)
    ontology = models.ForeignKey(Ontology, related_name='terms', db_column='ontology')
    name = models.TextField()

    class Meta:
        db_table = u'ontology_term'

class OntologySubset(models.Model): 
    name = models.TextField(primary_key = True)   
    ontology = models.ForeignKey(Ontology, related_name='subsets', db_column='ontology')

    class Meta:
        db_table = u'ontology_subset'

class OntologySubsetMembership(models.Model):
    id = models.AutoField(primary_key = True)
    term = models.ForeignKey(OntologyTerm, db_column='term')
    subset = models.ForeignKey(OntologySubset, db_column='subset')

    class Meta:
        db_table = u'ontology_subset_membership'

class OntologyDBXRef(models.Model):
    id = models.AutoField(primary_key = True)           
    from_term = models.ForeignKey(OntologyTerm, related_name='dbxref', db_column='from_accession')
    to_accession = models.TextField()
    type = models.TextField()

    class Meta:
        db_table = u'ontology_dbxref'

class OntologyClosure(models.Model):
    id  = models.AutoField(primary_key = True)
    ancestor = models.ForeignKey(OntologyTerm, db_column='ancestor')
    descendant = models.ForeignKey(OntologyTerm, related_name='closure', db_column='descendant')
    type = models.TextField()

    class Meta:
        db_table = u'ontology_closure'

class Evidence(models.Model):
    code = models.TextField(primary_key = True)
    name = models.TextField()
    priority = models.IntegerField()

    class Meta:
        db_table = u'evidence'

class Annotation(models.Model):
    id  = models.AutoField(primary_key = True)
    uniprot = models.ForeignKey(UniProt)
    ontology_term = models.ForeignKey(OntologyTerm, db_column='ontology_accession')
    source = models.TextField()
    evidence = models.ForeignKey(Evidence, db_column='evidence_code')
    assigned_by = models.TextField()
    date_assigned = models.DateField()

    class Meta:
        db_table = u'annotation'

    def as_annotation(self):
        #TODO annotation types!
        return GOAnnotation({
                    'accession': self.ontology_term.accession,
                    'description': self.ontology_term.name,
                    'evidence_code' : self.evidence.code,
                    'evidence_description' : self.evidence.name,
                    'evidence_priority': self.evidence.priority,
                }) 

class Kerf70(models.Model):
    id = models.AutoField(primary_key = True)
    kerf_treenode_id = models.ForeignKey(TreeNode, db_column='kerf_tree_node_id')
    uniprot_accession = models.TextField()

    class Meta:
        db_table = u'kerf2_70_table'

class ConsensusUniProtTreeNodeNames(models.Model):
    id = models.AutoField(primary_key = True)
    tree_node = models.ForeignKey(TreeNode, related_name='consensus_description')
    consensus_uniprot_description = models.TextField()

    class Meta:
        db_table = u'consensus_uniprot_names'

class TreeNodeMRCA(models.Model):
    id = models.AutoField(primary_key = True)
    tree_node = models.ForeignKey(TreeNode, related_name='mrca')
    mrca = models.ForeignKey(UniProtTaxonomy)

    class Meta:
        db_table = u'tree_node_mrca'

class OhanaQueueStatus(models.Model):
    id = models.AutoField(primary_key = True)
    qstat_required = models.BooleanField()
    status = models.TextField()

    class Meta:
        db_table = u'ohana_queue_status'

class PhyloFactsUserProfile(models.Model):
    id = models.AutoField(primary_key = True)
    user = models.ForeignKey(User, related_name='profile')
    institution = models.TextField()
    institution_address = models.TextField()
    institution_country = models.TextField()
    user_position = models.TextField()
    has_ever_logged_in = models.BooleanField()
    maximum_fatcat_jobs = models.IntegerField()   
    maximum_phylobuilder_jobs = models.IntegerField()
 
    class Meta:
        db_table = u'phylofacts_user_profile'

class TreeNodeOrthology(models.Model):
    id = models.AutoField(primary_key = True)
    root = models.ForeignKey(TreeNode, related_name='orthology_group')
    oma = models.TextField()
    orthomcl = models.TextField()

    class Meta:
        db_table = u'tree_node_orthology'

class PhyloFactsEmail(models.Model):
    id = models.AutoField(primary_key = True)
    recipient = models.TextField()
    subject = models.TextField()
    text_body = models.TextField()
    html_body = models.TextField()
    
    class Meta:
        db_table = u'phylofacts_email'

class SequenceFamily(models.Model):
    id = models.AutoField(primary_key = True)
    sequence = models.ForeignKey(SequenceHeader, related_name = 'seq_hdr')
    family = models.ForeignKey(Family, related_name = 'sequence_family')

    class Meta:
        db_table = u'sequence_family'

# For the DIP interactions...shouldn't really use these except for timing interolog stuff
# If we want to use this for more, we need to extend the model

class DIP(models.Model):
    id = models.AutoField(primary_key = True)
    uniprot_interactor_a = models.ForeignKey(UniProt, related_name = 'int_a')
    uniprot_interactor_b = models.ForeignKey(UniProt, related_name = 'int_b')
    detection_methods = models.TextField()
    first_author_string = models.TextField()
    publication_string = models.TextField()
    interaction_types = models.TextField()
    source_db = models.TextField()
    interaction_identifier = models.TextField()
    confidence = models.TextField()

    class Meta:
        db_table = u'dip'

class UniProtXList(models.Model):
    id = models.AutoField(primary_key = True)
    uniprot_accession = models.TextField()
    cross_reference_accession = models.TextField()
    cross_reference_type = models.TextField()

    class Meta:
        db_table = u'uniprot_x_list'
