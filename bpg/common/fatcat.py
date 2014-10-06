#!/usr/bin/python

import re
from pfacts003.utils.id_patterns import *
from bpg.common.BPGPWID import pairwise_identity_belvu as pwid_belvu

def consensus_uniprot_description(treenodes, base=0.5, fatcat_query=None):
    """ This function will take a list of treenodes and rank their descendents' uniprot descriptions according 
    to the same scoring function as before.  Swissprot sequences get 10 times the weight,
    hypothetical, predicted, etc are ranked as uninformative. 

    This function returns a list of tuples, each one having 3 values, the header as a string,
    the score as an integer, and the raw number of votes for that description.
    """

    leaves = []

    # get all the leaves below these nodes.

    for treenode in treenodes:
        for leaf in treenode.get_informative_members():
            if leaf.sequence_header.description():
                leaves.append(leaf)

    descriptions = {}
    votes = {}
    return_list = []
    total_weighted_votes = 0

    for leaf in leaves:
        description = leaf.sequence_header.description()
        if fatcat_query:
            ortholog_percentID = leaf.alignment_to_query(fatcat_query)[0]
            weight = pow(base, (10-(ortholog_percentID//10)))
        else:
            weight = 1
        if leaf.sequence_header.uniprot and leaf.sequence_header.uniprot.in_swissprot_f:
            descriptions[description] = descriptions.get(description, 0) + 10*weight
            total_weighted_votes += 10*weight
        else:
            descriptions[description] = descriptions.get(description, 0) + 1*weight
            total_weighted_votes += 1*weight
        # get the raw votes
        votes[leaf.sequence_header.description()] = votes.get(leaf.sequence_header.description(), 0) + 1
    
    for d, v in descriptions.items():
        return_list.append((d, v, votes[d], float(1.0*v/total_weighted_votes)))
    
    return sorted(return_list, key=lambda desc: desc[3], reverse=True)


def fatcat_ortholog_consensus_uniprot_description(uniprot_list, fatcat_query, base=0.5):
    """ This function will take a list of treenodes and rank their descendents' uniprot descriptions according 
    to the same scoring function as before.  Swissprot sequences get 10 times the weight,
    hypothetical, predicted, etc are ranked as uninformative. 

    This function returns a list of tuples, each one having 3 values, the header as a string,
    the score as an integer, and the raw number of votes for that description.
    """

    descriptions = {}
    votes = {}
    return_list = []
    total_weighted_votes = 0

    for leaf in uniprot_list:
        description = leaf.description
        ortholog_percentID = leaf.pwid_to_query
        weight = pow(base, (10-(ortholog_percentID//10)))
        if leaf.in_swissprot_f:
            descriptions[description] = descriptions.get(description, 0) + 10*weight
            total_weighted_votes += 10*weight
        else:
            descriptions[description] = descriptions.get(description, 0) + 1*weight
            total_weighted_votes += 1*weight
        # get the raw votes
        votes[description] = votes.get(description, 0) + 1
    
    for d, v in descriptions.items():
        return_list.append((d, v, votes[d], float(1.0*v/total_weighted_votes)))
    
    return sorted(return_list, key=lambda desc: desc[3], reverse=True)

def fatcat_ortholog_gene_clustering(orthologs,h=5, qcov=70, ocov=70, pwid_thresh=50):
    """ This piece of code implements the ortholog clustering for the fatcat orthologs tab.
    The input is a dictionary from fatcats orthologs function which contains keys for each unique
    taxon in the candidate orthologs.  The value of each key is a list of uniprot objects for that taxon.

    The output of the function will be a list of orthologs with a cluster number attached to them.  
    Within each cluster, there will be exactly one representative.  
    This is a boolean attached to the uniprot object. """

    def _count_orthology_supported_ancestors(orthology_support):
        count = 0
        for (system, contains) in orthology_support.items():
            if contains:
                count+=1
        return count

    def _make_cluster_row_tuple(uniprots):
        cluster_row = []
        # orders all of this data
        for uni_list in uniprots:
            cluster_row.append(tuple(sorted(uni_list, key=lambda u: _count_orthology_supported_ancestors(u.orthology_support), reverse=True)))
        # sort by pwid to query, if there is a swissprot accession, move it to the front of the list
        cluster_row = sorted(cluster_row, key=lambda t: t[0].pwid_to_query, reverse=True)
        for t in cluster_row:
            if t[0].in_swissprot_f:
                cluster_row.remove(t)
                cluster_row.insert(0, t)
                break
        return tuple(cluster_row)
    
    uniprots_encountered = {}
    clusters = []
    paralogs = []
    orthologs_passed = {}

    # remove orthologs which don't pass bi-directional coverage
    for taxon_id in orthologs.keys():
        for uni in orthologs[taxon_id]:
            if ((uni.query_coverage < qcov) or (uni.ortholog_coverage < ocov) or (uni.pwid_to_query < pwid_thresh)):
                uniprots_encountered.setdefault(uni.accession, []).append(uni)
            else:
                orthologs_passed.setdefault(taxon_id, []).append(uni)

    for (accession, uni_list) in uniprots_encountered.items():
        paralogs.append(_make_cluster_row_tuple([uni_list]))

    uniprots_encountered = {}

    for taxon_id in orthologs_passed.keys():
        if uniprots_encountered:
            clusters.append(_make_cluster_row_tuple(uniprots_encountered.values()))
        orthologs_from_taxon = sorted(orthologs_passed[taxon_id], key=lambda o: o.pwid_to_query, reverse=True)
        cluster_pid = orthologs_from_taxon[0].pwid_to_query
        uniprots_encountered = {}
        for uni in orthologs_from_taxon:
            if (cluster_pid - uni.pwid_to_query) < h:
                uniprots_encountered.setdefault(uni.accession, []).append(uni)
            else:
                clusters.append(_make_cluster_row_tuple(uniprots_encountered.values()))
                cluster_pid = uni.pwid_to_query
                uniprots_encountered = {}
                uniprots_encountered.setdefault(uni.accession, []).append(uni)
    if uniprots_encountered:
        clusters.append(_make_cluster_row_tuple(uniprots_encountered.values()))

    # separate clusters into orthologs and paralogs.    
 
    ortholog_taxa = {}

    for cluster in clusters:
        if cluster[0][0].taxon.id in ortholog_taxa:
            if cluster[0][0].pwid_to_query > ortholog_taxa[cluster[0][0].taxon.id][0][0].pwid_to_query:
                ortholog_taxa[cluster[0][0].taxon.id] = cluster
            else:
                paralogs.append(cluster)
        else:
            ortholog_taxa[cluster[0][0].taxon.id] = cluster

    return (ortholog_taxa.values(), paralogs)   
