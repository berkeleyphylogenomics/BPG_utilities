# FAT-CAT pipeline utils:
#
# - trim a list of fasta sequences by getting the closest k to the first one
# - mask a gappy alignment (either chop columns or convert to lowercase)
# - wrapper for FastTree
#
# Dave Dineen

# get_top_k
from Bio import SeqIO
import subprocess

# mask
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
import copy
from pfacts003.phylofacts.models import Family, SequenceFamily, SequenceHeader, UniProt
from pfacts003.fatcat.consts import *
# uclust
import os, tempfile

# Given a bunch of sequences
# gets the closest k to the first one
# which is the query sequence

# basically a wrapper for getclosestquick
# which uses the k-tuple scoring code from Clustal Omega

def get_top_k_to_query(input_seqs, k=500):
	args = ["/clusterfs/ohana/software/gettopn/src/gettopn", "-i", "-", "--threads=1", "-n", str(k)]
	
	process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	return process.communicate(input_seqs)[0]

def test_get_top_k_to_query():
	test_seqs = open("bpg0164151.unaligned", "r").read()
	print get_top_k_to_query(test_seqs, 2)
    
# Given an alignment and assuming first sequence is the 'query'
# masks columns with more than % gaps but only if the query has a gap

# returns the alignment with the masked columns chopped out
# *and* modifies the alignment with the masked columns set to lowercase/dots

# to lowercase
def mask_msa(alignment, thresh_percent):
    threshold = thresh_percent * len(alignment)
    
    # count the gaps
    mask = [0] * len(alignment[0])
    for record in alignment:
        for i in xrange(len(record)):
            # query must have gap
            if record[i] is "-":# and alignment[0, i] is "-":
                mask[i] += 1
    
    # get the masked alignment
    for record in alignment:
        newseq = ""
        for i in xrange(len(mask)):
            if mask[i] <= threshold:
                newseq += record[i]
            elif record[i] is "-":
                newseq += "."
            else:
                newseq += record[i].lower()
        record.seq = Seq(newseq)
    
    return masked

# cuts columns
def chop_msa(alignment, thresh_percent):
    threshold = thresh_percent * len(alignment)
    
    # count the gaps
    mask = [0] * len(alignment[0])
    for record in alignment:
        for i in xrange(len(record)):
            # query must have gap
            if record[i] is "-":# and alignment[0, i] is "-":
                mask[i] += 1

    # get the chopped alignment
    new_alignment = MultipleSeqAlignment([])
    for record in alignment:
        newseq = ''
        for i in xrange(len(mask)):
            if mask[i] <= threshold:
                newseq += record[i]
        new_alignment.append(SeqRecord(Seq(newseq, IUPAC.protein), id=record.id, name=record.name, description=record.description))
    return new_alignment

def test_mask():
    #test_seqs = open("bpg0164151.unaligned", "r").read()
    test_seqs = AlignIO.read("top5.fa", "fasta")
    chopped_seqs = chop_msa(test_seqs, 0.5)
    print chopped_seqs.format("fasta")
    
# give it an alignment and it comes back with a tree
# pipes all the way, no temp files necessary

def fast_tree(alignment):
    args = ["/clusterfs/ohana/software/bin/FastTree", "-quiet", "-nopr"]
    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.communicate(alignment.format("fasta"))[0]

# multithreaded version of above
def fast_tree_mp(alignment):
    args = ['/clusterfs/ohana/software/bin/FastTreeMP', '-quiet', '-nopr']
    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.communicate(alignment.format('fasta'))[0]

def test_fasttree():
	test_seqs = AlignIO.read("chopped.fa", "fasta")
	print fast_tree(test_seqs)

# fast clustering of fasta files using bob edgar's usearch
# takes a string containing the seqs in fasta as input
# returns list of (cluster #, uniprot accession)
# makes a bunch of temp files unfortunately
def uclust(input_fasta, temp_directory, identity=0.95):
    # put input fasta into a temp file
    f_taxon = open(tempfile.mkstemp(prefix='uclust', dir=temp_directory)[1], 'w')
    f_taxon.write(input_fasta)
    inputfile = f_taxon.name
    f_taxon.close()

    args = ['/clusterfs/ohana/software/bin/usearch', '-cluster_fast', inputfile, '-id', str(identity), '-uc', inputfile + '.uc']

    # future work: do something more intelligent than just sending the output to null?
    process = subprocess.call(args, stdout=open(os.devnull, "w"), stderr=open(os.devnull, "w"), shell=False)
    
    clusters = []
    f_clusters = open(inputfile + '.uc', 'r')
    for line in f_clusters:
        if line[0] != 'C':
            clusters.append([int(line.split('\t')[1]), line.split('\t')[8].split('|')[1]])
    f_clusters.close()
    return clusters
    
def chunkify(list, num_sublists):
    # this function will take a list and make sublists returning every nth item
    return [list[list_num::num_sublists] for list_num in range(0, num_sublists)]

def ellis_island_get_best_family(uniprot_list, family_list):
    ''' This function takes a list of uniprot accessions and returns a single family that contains the
    greatest number of accessions in the uniprot_list '''
    # This query does everything pretty much
    raw_query = '''WITH seq_fam AS (SELECT sequence_family.family_id AS fam_id, uniprot.accession AS acc FROM uniprot, sequence_family, sequence_header WHERE  uniprot.accession = any (%s) AND sequence_family.family_id = any (%s) AND sequence_header.uniprot_id = uniprot.id AND sequence_family.sequence_id = sequence_header.id GROUP  BY acc, fam_id),  fam_fam AS (SELECT fam_id, Count(fam_id) AS num_members FROM seq_fam GROUP BY fam_id ORDER BY num_members DESC) SELECT fam_fam.fam_id as id, family.seed_sequence_header_id, family.average_chars, family.length_shortest, family.length_longest, family.minimum_identities, family.mean_identities, family.gaps, family.columns_w_blosum62_lt_0, family.longest_deletion_length, family.longest_insertion_length, family.build_database_source, family.family_specific_evalue_criterion, family.family_specific_sw_method, family.notes, family.build_date, family.status, family.private, family.gathering_method, family.family_type_id, family.partition, family.canonical_tree_id, family.build_alignment_notes_id, family.author_id, family.active, family.score, family.score_updated FROM fam_fam, family where family.id = fam_fam.fam_id limit 1'''
    print 'uniprot list'
    print uniprot_list
    print 'family list'
    print family_list
    return Family.objects.raw(raw_query, (uniprot_list, family_list))[0]

# version of ellis_island_get_best_family where we want most taxon matches rather than most accessions
def ellis_island_get_most_taxa(uniprot_list, family_list):
    ''' This function takes a list of uniprot accessions and returns a single family that contains the
    greatest number of accessions in the uniprot_list '''
    # This query does everything pretty much
    raw_query = '''WITH seq_fam AS (SELECT DISTINCT sequence_family.family_id AS fam_id, uniprot.taxon_id AS acc FROM uniprot, sequence_family, sequence_header WHERE  uniprot.accession = any (%s) AND sequence_family.family_id = any (%s) AND sequence_header.uniprot_id = uniprot.id AND sequence_family.sequence_id = sequence_header.id GROUP  BY acc, fam_id),  fam_fam AS (SELECT fam_id, Count(fam_id) AS num_members FROM seq_fam GROUP BY fam_id ORDER BY num_members DESC) SELECT fam_fam.fam_id as id, family.seed_sequence_header_id, family.average_chars, family.length_shortest, family.length_longest, family.minimum_identities, family.mean_identities, family.gaps, family.columns_w_blosum62_lt_0, family.longest_deletion_length, family.longest_insertion_length, family.build_database_source, family.family_specific_evalue_criterion, family.family_specific_sw_method, family.notes, family.build_date, family.status, family.private, family.gathering_method, family.family_type_id, family.partition, family.canonical_tree_id, family.build_alignment_notes_id, family.author_id, family.active, family.score, family.score_updated FROM fam_fam, family where family.id = fam_fam.fam_id limit 1'''

    return Family.objects.raw(raw_query, (uniprot_list, family_list))[0]
    
def get_members_in_accession_list_in_family(uniprot_list, family):
    ''' This function takes a uniprot list, and a family, and returns a uniprot accession list with
    the members in that family.'''
    print 'family id ' + str(family.id)
    print 'uniprot list ' + str(uniprot_list)
    print uniprot_list
    return [UniProt.objects.get(id=u['sequence__uniprot_id']).accession for u in SequenceFamily.objects.filter(sequence__uniprot__accession__in=uniprot_list, family=family).values('sequence__uniprot_id').distinct('sequence__uniprot_id')]

# DD: taxon version of the above
# ie get members on list in family plus members in that family with same taxon
# simple!
def get_members_for_elimination(uniprot_list, family):
    #raw_query = '''SELECT DISTINCT uniprot.* FROM sequence_family, sequence_header, uniprot WHERE sequence_family.family_id=%s AND uniprot.accession=any (%s) AND sequence_header.taxon_id=any (SELECT taxon_id FROM uniprot WHERE accession = any (%s)) AND sequence_family.sequence_id = sequence_header.id AND sequence_header.uniprot_id = uniprot.id'''
    raw_query = '''SELECT uniprot.* FROM sequence_header, uniprot WHERE uniprot.accession=any (%s) AND sequence_header.taxon_id=any (SELECT taxon_id FROM sequence_family, sequence_header WHERE sequence_family.family_id=%s AND sequence_header.id = sequence_family.sequence_id) AND sequence_header.uniprot_id = uniprot.id'''

    return list(set([u.accession for u in UniProt.objects.raw(raw_query, [uniprot_list, family.id])]))

    
def find_families_for_uniprot_accession_list(uniprot_list, family_list=None):
    ''' This function takes a list of uniprot accessions and returns a list of unique families that contain
        one or more sequences from the list. '''
    return [Family.objects.get(id=fam['family']) for fam in SequenceFamily.objects.filter(sequence__uniprot__accession__in=uniprot_list, family__id__in=family_list).order_by('family').values('family').distinct('family')]

def find_uniprots_for_family_list(family_list):
    ''' This function takes a list of family ids as input and returns the unique set of family members 
    that have uniprot records, and also where the fasta sequence is stored in the database. ''' 
    return [UniProt.objects.get(id=u['sequence__uniprot_id']) for u in SequenceFamily.objects.filter(family__id__in=family_list, sequence__uniprot__isnull=False, sequence__uniprot__sequence_chars__isnull=False).order_by('sequence__uniprot__accession').values('sequence__uniprot_id').distinct('sequence__uniprot__accession')]

def cluster_member_list(list, cluster_number_start=1):
    ''' This function implements the stage 4 clustering.  The idea is exactly the same as uclust does really,
    except it uses member to query pid to cluster the members.  Maybe we should switch to uclust....
    probably should. '''
    taxonomies = {}
    cluster_num = cluster_number_start

    for member in list:
        taxon = member.uniprot.taxon.id
        if member.uniprot.taxon.id in taxonomies:
            for cluster in taxonomies[taxon]:
                if ((cluster['members'][0].member_to_query_pid - member.member_to_query_pid) <= STAGE4_CLUSTERING_THRESHOLD):
                    # add member to this cluster
                    cluster['members'].append(member)
                    # is this better than the current representative?
                    if member.uniprot.in_swissprot_f and not cluster['representative'].uniprot.in_swissprot_f:
                        cluster['representative'] = member
                    break    
 
            else:     
                taxonomies[taxon].append({'members': [member],
                                          'representative': member,
                                          'cluster_num': cluster_num})
                cluster_num += 1
        else:
            taxonomies[taxon] = [{'members': [member],
                                  'representative': member,
                                  'cluster_num': cluster_num}]
            cluster_num += 1
    return cluster_num, taxonomies
