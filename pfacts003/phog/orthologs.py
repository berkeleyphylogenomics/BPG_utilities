#from pfacts003.phylofacts.models import TreeNode, UniProt, UniProtGI, OrthologTypes, \

import time

from django.db.models import Q

from pfacts003.utils.id_patterns import *
from pfacts003.utils.messages import unrecognized_identifier, not_in_phog
from pfacts003.phylofacts.models import TreeNode, UniProt, OrthologTypes, \
                              get_ortholog_type_of_threshold, \
                              getTreeNodeQuerySetForTreeNodes

def getPHOGQuerySet(phog_accession, orderByTaxonomicDistribution = False,
                    ortholog_type = OrthologTypes.SuperOrtholog,
                    threshold = 0.0):
  m = phog_re.match(phog_accession)
  if m is None:
    return (None, None, 'Unrecognized PHOG accession')
  tree_id = int(phog_accession[4:11])
  left_id = int(phog_accession[12:17])
  phog = TreeNode.objects.get(tree__id = tree_id, left_id = left_id)
  if phog == None:
    return (None, None, 'Nonexistent PHOG accession')
  if phog.tree.family.status == 'bad':
    return (None, None, 'Bad PHOG accession')
  all_orthologs = phog.get_included_leaves(ortholog_type, threshold)
  ortholog_dict = {}
  ortholog_set = set()
  for ortholog in all_orthologs:
    identifier = ortholog.sequence_header.identifier()
    if identifier not in ortholog_dict:
      ortholog_dict[identifier] = ortholog
      ortholog_set.add(ortholog)
  if orderByTaxonomicDistribution:
    orthologs = getTreeNodeQuerySetForTreeNodes(ortholog_set).order_by(
                          'sequence_header__uniprot__taxon__left_id')
  else:
    orthologs = getTreeNodeQuerySetForTreeNodes(ortholog_set).order_by(
                      'sequence_header__uniprot__taxon__scientific_name')
  return (phog, orthologs, None)
  

def getOrthologQuerySet(sequence_id, ortholog_type, threshold = 0.0,
                        orderByTaxonomicDistribution = False):

  """Get orthologs and PHOGs from database

  Although this function is titled "Get Ortholog Query Set", what is
  actually returned is a 6-tuple representing the results of that
  query. The ortholog query set is contained within that 6-tuple.

  The two mandatory arguments expected are:
      1) The identifier to query against (sequence_id as a string), and
      2) The ortholog type (as defined in the
         pfacts003.phylofacts.models.OrthologTypes object

  The optional arguments are:
      1) A boolean to determine if the results should be ordered by
         taxonomic distribution, and 
      2) A threshold value. Although the threshold value is pre-
         determined by pfacts003.phylofacts.models.preset_thresholds for
         all but custom ortholog types, it is currently expected that one
         looks up and passes the correct threshold value even when using
         these predetermined types.

  A 6-tuple is always returned.
  
  The first item in the 6-tuple returned is an instance the Django
  UniProt model that corresponds to the 'sequence_id' parameter given.
  If no object was found, this position in the 6-tuple is None.

  The second item in the 6-tuple is a list of PHOGs (Django TreeNode
  objects) that represent the PhyloFacts Orthology Groups that represent
  this sequence. This is a redundant item as the same results can be
  derived from the values of the fifth item in this 6-tuple. If none
  were found, this value will be None.

  The third object returned is a set of PHOGs (Django TreeNode objects)
  that represents the best PHOGs for the query. It is a subset of the
  second item in the 6-tuple. If none were found, this value will be
  None

  The fourth object returned is a django.db.models.query.QuerySet
  represented as a list of all orthologs (TreeNode objects) found for
  the query. This is a redundant item as the same results can be derived
  from the keys of the fifth item in this 6-tuple.

  The fifth object returned is a PHOG-Ortholog dictionary for looking up
  which PHOGs an Ortholog corresponds to.

  And, the final object, the sixth object, in the 6-tuple represents any
  error messages to display. This was written to directly interface with
  the web-server and thus contains HTML mark-up embedded directly within
  the message.

  An example use of this function follows:

  >>> from pfacts003.phylofacts.models import OrthologTypes, preset_thresholds
  >>> from pfacts003.phog.orthologs import getOrthologQuerySet
  >>> from pprint import pprint
  >>> query = getOrthologQuerySet('PE2R2_HUMAN',
  ...  OrthologTypes.PHOG_T_Tight,
  ...  preset_thresholds[OrthologTypes.PHOG_T_Tight])
  >>> pprint(query)
  (<UniProt: PE2R2_HUMAN>,
   [<TreeNode: PHOG069746640>],
   set([<TreeNode: PHOG069746640>]),
   [<TreeNode: 69746871>, <TreeNode: 69746819>, <TreeNode: 69746856>,
    <TreeNode: 69746881>, <TreeNode: 69746831>, <TreeNode: 69746745>,
    <TreeNode: 69746758>, <TreeNode: 69746762>, <TreeNode: 69746776>,
    <TreeNode: 69746844>, <TreeNode: 69746701>, <TreeNode: 69746890>],
   {<TreeNode: 69746844>: <TreeNode: PHOG069746640>,
    <TreeNode: 69746856>: <TreeNode: PHOG069746640>,
    <TreeNode: 69746871>: <TreeNode: PHOG069746640>,
    <TreeNode: 69746881>: <TreeNode: PHOG069746640>,
    <TreeNode: 69746890>: <TreeNode: PHOG069746640>,
    <TreeNode: 69746701>: <TreeNode: PHOG069746640>,
    <TreeNode: 69746745>: <TreeNode: PHOG069746640>,
    <TreeNode: 69746758>: <TreeNode: PHOG069746640>,
    <TreeNode: 69746762>: <TreeNode: PHOG069746640>,
    <TreeNode: 69746776>: <TreeNode: PHOG069746640>,
    <TreeNode: 69746819>: <TreeNode: PHOG069746640>,
    <TreeNode: 69746831>: <TreeNode: PHOG069746640>},
   None)
  """

  # If the custom value is the same as one of the default types, use the
  # default instead -- as they have been pre-computed. Using the
  # default types is faster.
  if ortholog_type == OrthologTypes.PHOG_T_Custom:
    ortholog_type = get_ortholog_type_of_threshold(threshold)

  # Construct a Q object that has all of the PHOGs containing this
  # sequence.
  if uniprot_accession_re1.match(sequence_id) != None or \
      uniprot_accession_re2.match(sequence_id) != None:
    # This sequence_id is a UniProt accession
    query_nodes = TreeNode.objects.filter(
            sequence_header__uniprot__accession__exact=sequence_id)

    # Find all of the UniProt objects matching this sequence_id
    # if none are found, there is no need to continue
    ownUniProts = list(UniProt.objects.filter(
                          accession__exact=sequence_id))
    if len(ownUniProts) == 0:
      return (None, None, None, None, {}, unrecognized_identifier)

  elif uniprot_identifier_re1.match(sequence_id) != None or \
      uniprot_identifier_re2.match(sequence_id) != None or \
      uniprot_identifier_re3.match(sequence_id) != None:
    # This sequence_id is a UniProt identifier
    query_nodes = TreeNode.objects.filter(
            sequence_header__uniprot__uniprot_identifier__exact=sequence_id)

    # Find all of the UniProt objects matching this sequence_id
    ownUniProts = list(UniProt.objects.filter(
                          uniprot_identifier__exact=sequence_id))
    
    if len(ownUniProts) == 0:
      return (None, None, None, None, {}, unrecognized_identifier)

  elif gi_re.match(sequence_id) != None:
    # This sequence id is a Genbank ID
    gis = list(UniProtGI.objects.filter(genbank_id__exact=sequence_id))

    # Find all of the UniProt objects matching this GenBank ID
    ownUniProts = [gi.uniprot for gi in gis]
    if len(ownUniProts) == 0:
      return (None, None, None, None, {}, unrecognized_identifier)

    uniprot_accessions = [uniprot.accession for uniprot in ownUniProts]
    query_nodes = TreeNode.objects.filter(
        sequence_header__uniprot__accession__in=uniprot_accessions).exclude(tree__family__status__exact = 'bad')
  else:
    # This point in code should never be reached due to the input
    # validation above. As a security measure, if it is reached, return
    # the unrecognized identifier message
    return (None, None, None, None, {}, unrecognized_identifier)

  # Theoretically, there should only be one UniProt object. If there is more than one,
  # just take the first
  ownUniProt = ownUniProts[0]

  # If the ortholog type was custom, at this point we have the tree_node
  # objects corresponding to the query sequence, but not a Q object to
  # query for the PHOGs containing them.  Now we need to find the
  # containing PHOGs.
  phog_set = set()
  query_nodes = query_nodes.filter(tree__family__active = True)
  for query_node in query_nodes:
    containing_phog = query_node.get_containing_phog(threshold=threshold)
    if containing_phog:
      phog_set.add(containing_phog)
  phogs = getTreeNodeQuerySetForTreeNodes(phog_set)
  # Find all tree_node objects belonging to the PHOGs.
  all_ortholog_set = set()
  phog_of_tree = {}
  for phog in phogs:
    # Store which PHOG came from which tree so that we can use it to
    # find which PHOG each custom ortholog came from.
    phog_of_tree[phog.tree] = phog
    for ortholog in phog.get_included_leaves(ortholog_type, threshold):
      all_ortholog_set.add(ortholog)
  all_orthologs = getTreeNodeQuerySetForTreeNodes(all_ortholog_set)
  if not bool(all_orthologs):
    return (None, None, None, None, {}, not_in_phog(sequence_id))

  # Make a set of orthologs that is nonredundant in terms of identifiers
  # For each identifier, there may be several PHOGs displaying the
  # orthology between this identifier and the query.  We will consider
  # the one of these containing the most sequences to be the best.  So,
  # collect all the different pairs of (number of sequences contained in
  # that PHOG, tree node id of that ortholog) for each identifier.
  ortholog_dict = {}
  checked_ortholog_ids = set()
  phog_of_ortholog = {}
  for ortholog in all_orthologs:
    if ortholog.id in checked_ortholog_ids:
      continue
    identifier = ortholog.sequence_header.identifier()
    if identifier not in ortholog_dict:
      ortholog_dict[identifier] = set()
    phog_of_ortholog[ortholog] = phog_of_tree[ortholog.tree]
    num_contained_leaves = phog_of_ortholog[ortholog].get_included_leaves(
                                            ortholog_type, threshold).count()
    ortholog_dict[identifier].add((num_contained_leaves, ortholog))
    checked_ortholog_ids.add(ortholog.id)
  best_ortholog_set = set()
  for identifier in ortholog_dict.keys():
    candidate_ids = list(ortholog_dict[identifier])
    if len(candidate_ids) > 0:
      # Sort so that the one with the most sequences in the containing
      # PHOG comes first.
      candidate_ids.sort(reverse=True)
      num_contained_leaves, ortholog = candidate_ids[0]
      best_ortholog_set.add(ortholog)
  # Now we can get one tree_node object per identifier, using the
  # best_ortholog_ids.
  if orderByTaxonomicDistribution:
    orthologs = getTreeNodeQuerySetForTreeNodes(
                                  best_ortholog_set).distinct().order_by(
                                  'sequence_header__uniprot__taxon__left_id')
  else:
    orthologs = getTreeNodeQuerySetForTreeNodes(
                      best_ortholog_set).distinct().order_by(
                      'sequence_header__uniprot__taxon__scientific_name')
  # The best PHOGs are the ones that contained the best orthologs,
  # namely the ones containing the most sequences.
  best_phogs = set([phog_of_ortholog[ortholog] for ortholog in orthologs])

  return (ownUniProt, phogs, best_phogs, orthologs, phog_of_ortholog, None)
