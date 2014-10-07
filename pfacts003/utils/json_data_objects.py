#!/usr/bin/env python

from pfacts003.phylofacts.models import TreeNode, SequenceHeader, UniProt, \
    UniProtEC, EC, UniProtLiterature, UniProtPDB_Chain, UniProtGO, Family, \
    OrthologTypes, TreeNodeAlignment, TreeNodeAlignmentConservation, Tree
import os, re, string, glob
import json
from utils import remove_duplicates
from Bio import SubsMat
from Bio.SubsMat import MatrixInfo
from Bio.pairwise2 import dictionary_match
import cStringIO
from bpg.common.utils.phog_summary_tree import phog_summary_tree

blosum62_of_residues = dictionary_match(SubsMat.SeqMat(MatrixInfo.blosum62))

def _uniq(seq, idfun=None):
    # order preserving
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        # in old Python versions:
        # if seen.has_key(marker)
        # but in new ones:
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result

def get_alignment_rows(tree_node, pageSize=20, pageNum=1):
  end_index = pageNum * pageSize
  start_index = end_index - pageSize
  alignment_objects = TreeNodeAlignment.objects.filter(tree_node = tree_node)
  alignment_object_of_sequence_header = {}
  for alignment_object in alignment_objects:
    alignment_object_of_sequence_header[alignment_object.sequence_header] \
        = alignment_object
  if len(alignment_object_of_sequence_header) != len(alignment_objects):
    print "Discrepancy: multiple alignment objects w/ same seq header"
  leaves = tree_node.get_included_leaves(
                ortholog_type = OrthologTypes.PHOG_T_Custom,
                threshold=10000.0).order_by('left_id')
  ordered_alignment_objects = [alignment_object_of_sequence_header[
                                leaf.sequence_header] for leaf in leaves]
  alignment_rows = {}
  row_index = 0
  info_of_column = {}
  conservation_infos \
      = TreeNodeAlignmentConservation.objects.filter(tree_node = tree_node)
  for info in conservation_infos:
    info_of_column[info.column_index] = info
  row_index = 0
  for i in xrange(start_index, end_index):
    alignment_object = ordered_alignment_objects[i]
    row_index += 1
    alignment_row = {}
    alignment_row['identifier'] \
        = alignment_object.sequence_header.identifier()
    alignment_row['identifier'] = string.replace(alignment_row['identifier'],
                                                  '"', '\\"')
    alignment_row['description'] \
        = alignment_object.sequence_header.description()
    if len(alignment_row['description']) > 25:
      alignment_row['description'] = alignment_row['description'][0:25] + '...'
    alignment_row['residueNum'] = {}
    alignment_row['class'] = {}
    seq = alignment_object.aligned_sequence.chars
    k = 1
    for j in xrange(len(seq)):
      residue = seq[j]
      spec={}
      alignment_row['col%d' % j] = residue
      alignment_row['residueNum']['col%d' % j] = -1
      alignment_row['class']['col%d' % j] = 'none'
      if residue.isupper():
        if blosum62_of_residues(residue, info_of_column[j].conserved_residue) \
            >= info_of_column[j].blosum62_conservation_score:
          if info_of_column[j].blosum62_conservation_score >= 3:
            alignment_row['class']['col%d' % j] = 'align_high'
          elif info_of_column[j].blosum62_conservation_score >= 1.5:
            alignment_row['class']['col%d' % j] = 'align_moderate'
          elif info_of_column[j].blosum62_conservation_score >= 0.5:
            alignment_row['class']['col%d' % j] = 'align_low'
        alignment_row['residueNum']['col%d' % j] = k
        k += 1
      elif residue.islower():
        alignment_row['residueNum']['col%d' % j] = k
        k += 1
    alignment_rows[row_index] = alignment_row
  retf = cStringIO.StringIO()
  retf.write('{')
  retf.write('data:[')
  firstRow=True
  for row_index in alignment_rows:
    alignment_row = alignment_rows[row_index]
    if firstRow:
      firstRow=False
    else:
      retf.write(',')
    retf.write('{')
    retf.write('"seqid":')
    retf.write('"%s"' % alignment_row['identifier'])
    retf.write(',')
    retf.write('"desc":')
    retf.write('"%s"' % alignment_row['description'])
    retf.write(",")
    for j in xrange(len(seq)):
      retf.write('"c%d":' % j)
      retf.write('"%s"' % seq[j])
      retf.write(",")
    retf.write('"resnos":')
    retf.write('[')
    for j in xrange(len(seq)):
      if j > 0:
        retf.write(',')
      retf.write('%d' % alignment_row['residueNum']['col%d' % j])
    retf.write(']')
    retf.write(',')
    retf.write('"classes":')
    retf.write('[')
    for j in xrange(len(seq)):
      if j > 0:
        retf.write(',')
      retf.write('"%s"' % alignment_row['class']['col%d' % j])
    retf.write(']')
    retf.write('}')
  retf.write(']')
  retf.write(',')
  retf.write('pageInfo:{totalRowNum:%d},' % len(alignment_objects))
  retf.write("recordType:'object'")
  retf.write('}')
  ret_string = retf.getvalue()
  ret_string = ret_string.replace(unichr(1), ' ')
  return ret_string

def _annotate_nhx_element(bpg_accession, sequence_header_ids):
  """ Helper function for annotateTree. """
  family = Family.objects.get(id = int(bpg_accession[4:]))
  sequence_headers = SequenceHeader.objects.filter(id__in =
                                                    sequence_header_ids)
  annotation_info = {}
  sequence_headers_of_uniprot = {}
  for sequence_header in sequence_headers:
    annotation_info[sequence_header] = {}
    annotation_info[sequence_header]['header'] = sequence_header.header
    annotation_info[sequence_header]['in_swissprot_f'] = 0
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
  return ret

def fileExtension(urlTreeMethod):
  """Return the correct file extent for the tree method passed in the URL; default is
    'nj' for Neighbor-Joining."""
  if (urlTreeMethod =='ml'):
    return 'phyml'
  else:
    return 'nj'

trivial_translation = string.maketrans('', '')
seqhdr_re = re.compile('SEQHDR([0-9]+)')

def _annotate_phog_element(family_accession, sequence_header_ids,
                            tree_method='ml', threshold=0.0):
  """Create a dictionary of annotations of the PHOGs in a family.

  Here sequence_header_ids are the ids of the sequence_headers in the family
  that are *not* contained in any PHOGs at this threshold.
  """
  
  tree = Tree.objects.get(family__id = int(family_accession[4:]), method =
                          tree_method)
  phogs = tree.get_phogs(threshold)
  annotation_info = {}
  leaves_in_phogs = TreeNode.objects.filter(tree__exact = tree,
              sequence_header__isnull = False).exclude(sequence_header__in =
              sequence_header_ids).order_by('left_id')
  leaf_index = 0
  for phog in phogs:
    annotation_info[phog] = {}
    annotation_info[phog]['uniprot_de'] = phog.get_description()
    annotation_info[phog]['in_swissprot_f'] = 0
    if phog.get_has_swissprot():
      annotation_info[phog]['in_swissprot_f'] = 1
    annotation_info[phog]['scientific_name'] = 'N/A'
    annotation_info[phog]['common_name'] = 'N/A'
    annotation_info[phog]['ncbi_taxid'] = 'N/A'
    taxon = phog.get_taxon()
    if taxon is not None:
      annotation_info[phog]['ncbi_taxid'] = taxon.id
      if taxon.scientific_name:
        annotation_info[phog]['scientific_name'] = taxon.scientific_name
      if taxon.common_name:
        annotation_info[phog]['common_name'] = taxon.common_name
    ecs = phog.get_ecs()
    if ecs and len(ecs) > 0:
      annotation_info[phog]['ec'] = ';'.join([ec.__unicode__() for ec in ecs])
    pmids = phog.get_pmids()
    if pmids and len(pmids) > 0:
      annotation_info[phog]['Literature'] = pmids
    pdb_chains = phog.get_pdb_chains()
    if pdb_chains and len(pdb_chains) > 0:
      annotation_info[phog]['ThreeDStructure'] = pdb_chains
    go_data = phog.go_data()
    for term_type in go_data:
      annotation_info[phog][term_type] = set()
      for name in go_data[term_type]:
        annotation_info[phog][term_type].add((name,
                                    go_data[term_type][name]['evidence_code']))

  ret = {}
  for phog in phogs:
    ret[phog.get_accession()] = annotation_info[phog]
    if 'Literature' in annotation_info[phog]:
      ret[phog.get_accession()]['Literature'] \
          = list(annotation_info[phog]['Literature'])
    if 'ThreeDStructure' in annotation_info[phog]:
      ret[phog.get_accession()]['ThreeDStructure'] \
          = ['%s%s' % (pdb_chain.pdb_id, pdb_chain.chain_id) for pdb_chain in
            annotation_info[phog]['ThreeDStructure']]
    if 'biological_process' in annotation_info[phog]:
      ret[phog.get_accession()]['biological_process'] \
          = list(annotation_info[phog]['biological_process'])
    if 'molecular_function' in annotation_info[phog]:
      ret[phog.get_accession()]['molecular_function'] \
          = list(annotation_info[phog]['molecular_function'])
    if 'cellular_component' in annotation_info[phog]:
      ret[phog.get_accession()]['cellular_component'] \
          = list(annotation_info[phog]['cellular_component'])
  return ret

def annotatedSummaryTree(family_accession, tree_method = 'ml', threshold=0.0):
  _writer = json.JsonWriter()
  newick_tree = phog_summary_tree(family_accession, method=tree_method,
                      threshold=threshold)
  ret = {}
  sequence_header_ids = [int(seqhdr_str) for seqhdr_str 
                          in seqhdr_re.findall(newick_tree)]
  ret['SequenceHeaderAnnotations'] = _annotate_nhx_element(family_accession,
                                                          sequence_header_ids)
  ret['PHOG_Annotations'] = _annotate_phog_element(family_accession, 
                                                sequence_header_ids,
                                                tree_method,
                                                threshold)
  ret['__tree__'] = newick_tree
  ret_string = _writer.write(ret)
  return ret_string.replace(unichr(1), ' ')

                          
def annotatedTree(accession, tree_method = 'nj'):
    """ Retrieve tree annotation-information from a family, based on bpg accession number,
        and return JSON-formatted representation of the needed fields. """
    _writer = json.JsonWriter()
    family_root_dir = '/clusterfs/ohana/bpg/pfacts/'
    family_dir = os.path.join(family_root_dir, accession[0:4], accession[0:7],
                            accession)
    tree_files = glob.glob('%s/%s.%s' % (family_dir, accession, tree_method))
    if len(tree_files) == 0:
      print "No tree file found: %s, %s, %s" \
            % (family_dir, accession, tree_method)
      return ''
    f = open(tree_files[0])
    nhxTree = f.read().translate(trivial_translation, string.whitespace)
    f.close()

    sequence_header_ids = [int(seqhdr_str) for seqhdr_str 
                            in seqhdr_re.findall(nhxTree)]

    ret = {}
    ret = _annotate_nhx_element(accession, sequence_header_ids)
    nhxTree = nhxTree.replace('SEQHDR','')
    ret['__tree__'] = nhxTree
    ret_string = _writer.write(ret)
    return ret_string.replace(unichr(1), ' ')

def superorthologousNodes(accession, tree_method = 'ml'):
    _writer = json.JsonWriter()
    querySet = TreeNode.objects.filter(duplication_distance__gte = 0.0,
                  greatest_duplication_distance_of_maximal_descendant = -1,
                  tree__family__id__exact = int(accession[4:]),
                  tree__method = tree_method)
    left_ids = [tree_node.left_id for tree_node in querySet]
    ret = {}
    ret['superorthologous_left_ids'] = left_ids
    return(_writer.write(ret))

def hyperNeighborhood(phog, ortholog_type, threshold, taxon_id):
  edges = set()
  hyper_neighbors = phog.get_hyper_neighbors(ortholog_type, threshold)
  if taxon_id:
    taxon = NCBITaxonomy.objects.get(id__exact = taxon_id)
  nodes = {}
  nodes[phog.get_accession(ortholog_type, threshold)] = {}
  nodes[phog.get_accession(ortholog_type, threshold)]['description'] \
      = phog.get_description(ortholog_type, threshold)
  if taxon_id:
    nodes[phog.get_accession(ortholog_type, threshold)]['genes_from_taxon'] \
        = remove_duplicates([leaf.sequence_header.identifier() for leaf in
            phog.get_contained_leaves_from_taxon(taxon, ortholog_type, threshold)])
  for type in hyper_neighbors.keys():
    for neighbor in hyper_neighbors[type]:
      nodes[neighbor.get_accession(ortholog_type, threshold)] = {}
      nodes[neighbor.get_accession(ortholog_type, threshold)]['description'] \
          = neighbor.get_description(ortholog_type, threshold)
      if taxon_id:
        nodes[neighbor.get_accession(ortholog_type, 
                                      threshold)]['genes_from_taxon']\
            = remove_duplicates([leaf.sequence_header.identifier() for leaf in
                neighbor.get_contained_leaves_from_taxon(taxon, ortholog_type, 
                                                          threshold)])
  for type in hyper_neighbors.keys():
    for hyper_neighbor in hyper_neighbors[type]:
      edges.add(((phog.get_accession(ortholog_type, threshold), 
                  hyper_neighbor.get_accession(ortholog_type, threshold)), type))
      neighbors_of_hyper_neighbors \
        = hyper_neighbor.get_hyper_neighbors(ortholog_type, threshold)
      if type in neighbors_of_hyper_neighbors:
        relevant_neighbors \
          = neighbors_of_hyper_neighbors[type] & hyper_neighbors[type]
        for neighbor in relevant_neighbors:
          edges.add(((hyper_neighbor.get_accession(ortholog_type, threshold), 
                    neighbor.get_accession(ortholog_type, threshold)), type))
  ret = {}
  ret['query'] = phog.get_accession(ortholog_type, threshold)
  ret['edges'] = list(edges)
  ret['nodes'] = nodes
  _writer = json.JsonWriter()
  # print ret
  return(_writer.write(ret))

# This is called by run_interactome_viewer to get the list of edges in non-JSON form.
# It probably should not be here since it is not returning a JSON object.
# It might be replaced by using hyperNeighborhood.
def simpleHyperNeighborhood(phog, ortholog_type, threshold, taxon_id):
  edges = set()
  hyper_neighbors = phog.get_hyper_neighbors(ortholog_type, threshold)
  if taxon_id:
    taxon = NCBITaxonomy.objects.get(id__exact = taxon_id)
  nodes = {}
  nodes[phog.get_accession(ortholog_type, threshold)] = {}
  nodes[phog.get_accession(ortholog_type, threshold)]['description'] \
      = phog.get_description(ortholog_type, threshold)
  if taxon_id:
    nodes[phog.get_accession(ortholog_type, threshold)]['genes_from_taxon'] \
        = [leaf.sequence_header.identifier() for leaf in
            phog.get_contained_leaves_from_taxon(taxon, ortholog_type, threshold)]
  for type in hyper_neighbors.keys():
    for neighbor in hyper_neighbors[type]:
      nodes[neighbor.get_accession(ortholog_type, threshold)] = {}
      nodes[neighbor.get_accession(ortholog_type, threshold)]['description'] \
          = neighbor.get_description(ortholog_type, threshold)
      if taxon_id:
        nodes[neighbor.get_accession(ortholog_type, 
                                      threshold)]['genes_from_taxon']\
            = [leaf.sequence_header.identifier() for leaf in
                neighbor.get_contained_leaves_from_taxon(taxon, ortholog_type, 
                                                          threshold)]
  for type in hyper_neighbors.keys():
    for hyper_neighbor in hyper_neighbors[type]:
      edges.add(((phog.get_accession(ortholog_type, threshold), 
                  hyper_neighbor.get_accession(ortholog_type, threshold)), type))
      neighbors_of_hyper_neighbors \
        = hyper_neighbor.get_hyper_neighbors(ortholog_type, threshold)
      if type in neighbors_of_hyper_neighbors:
        relevant_neighbors \
          = neighbors_of_hyper_neighbors[type] & hyper_neighbors[type]
        for neighbor in relevant_neighbors:
          edges.add(((hyper_neighbor.get_accession(ortholog_type, threshold), 
                    neighbor.get_accession(ortholog_type, threshold)), type))
  ret = {}
  ret['query'] = phog.get_accession(ortholog_type, threshold)
  ret['edges'] = edges
  ret['nodes'] = nodes 
  return ret

def _getParents(taxa, _cur):
    sql = "SELECT id, parent_id FROM ncbi_taxonomy WHERE id IN (%s)" % ",".join(taxa)
    _cur.execute(sql)
    rows = _cur.fetchall()
    return rows

def _recursiveGetParent(tree, _cur):
    # A tree is a list of lists, where each list in the tree starts with some taxon and follows its parents
    # up toward the root. We grow the lists in the tree by fetching the parent of the last element of each
    # list. If at any time the same taxid appears in all lists in the tree, that's the parent.
    current = [x[-1] for x in tree]
    data = _getParents(current, _cur)
    # Put the parents we just fetched into the tree.
    for taxon in tree:
        row = [x for x in data if str(x['id']) == taxon[-1]]
        if 0 == len(row): # this taxon wasn't found in the database, skip it
            tree.remove(taxon)
            return _recursiveGetParent(tree, _cur)
        if 0 == row[0]['parent_id']:
            taxon.append('1')
        else:
            taxon.append(str(row[0]['parent_id']))
    
    # Does the last element in any taxon's list appear in all branches of the tree?
    for taxon in tree:
        count = 0
        for x in tree:
            if taxon[-1] in x: count += 1
        # If so, that's the parent; return it.
        if count == len(tree): 
            return taxon[-1]
    # Otherwise, recurse and try the next level up.
    return _recursiveGetParent(tree, _cur)

def _getlin(taxid):
    _db = MySQLdb.connect(db=phylofacts_db_name, host=phylofacts_db_host, user=phylofacts_db_user, passwd=phylofacts_db_passwd)
    _cur = _db.cursor(MySQLdb.cursors.DictCursor)
    if 0 == len(taxid): 
        return _writer.write(None)
    sql = "SELECT id, scientific_name, common_name, parent_id FROM ncbi_taxonomy WHERE id = %s"
    num_rows = _cur.execute(sql % taxid)
    if 0 == num_rows: 
        return None
    row = _cur.fetchone()
    ids = [row]
    while ids[-1]['parent_id'] != 1:
        _cur.execute(sql % ids[-1]['parent_id'])
        row = _cur.fetchone()
        ids.append(row)
    return ids

def getMRCA(taxid=[], withNonCellular=False):
    _writer = json.JsonWriter()
    _db = MySQLdb.connect(db=phylofacts_db_name, host=phylofacts_db_host, \
                      user=phylofacts_db_user, passwd=phylofacts_db_passwd)
    _cur = _db.cursor(MySQLdb.cursors.DictCursor)
    if 0 == len(taxid): 
        return _writer.write(None)
    if isinstance(taxid, list):
        taxid = _uniq(taxid)
    if 1 == len(taxid): 
        taxid = taxid[0]
        
    if isinstance(taxid, basestring):
        if _getlin(taxid)[-1]['id'] != 131567 and not withNonCellular:
            return _writer.write(None)
        else:
            sql = "SELECT id, scientific_name, common_name FROM ncbi_taxonomy WHERE id = %s" % taxid
    else:
        if withNonCellular:
            tree = [[taxon] for taxon in taxid]
        else:
            tree = []
            for taxon in taxid:
                lineage = _getlin(taxon)
                if lineage is not None and lineage[-1]['id'] == 131567:
                    tree.append([taxon])
        if len(tree) == 0:
            return _writer.write(None)
        elif len(tree) == 1:
            sql = "SELECT id, scientific_name, common_name FROM ncbi_taxonomy WHERE id = %s" % tree[0]
        else:
            parent = _recursiveGetParent(tree, _cur)
            sql = "SELECT id, scientific_name, common_name FROM ncbi_taxonomy WHERE id = %s" % parent
    numRows = _cur.execute(sql)
    if numRows > 0:
        row = _cur.fetchone()
        return _writer.write(row)
    else:
       return _writer.write(None)

def getLineage(taxid):
    _writer = json.JsonWriter()
    return _writer.write(_getlin(taxid))
