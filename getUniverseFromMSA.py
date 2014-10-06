#!/usr/bin/python

import os
import sys
import cPickle
from optparse import OptionParser
from Bio import SeqIO
from Bio.Blast import NCBIStandalone
from Bio.Blast import NCBIXML



def getUniverseFromMSA(fname, propid, blastdb, ecutoff, use_subfamily_seeds):

    # executables
    uniqueseq = '/usr/bin/uniqueseq'
    blastall = '/usr/bin/blastall'
    fastacmd = '/usr/bin/fastacmd'

    if use_subfamily_seeds:
      f = open("/home/ruchira/pfam_subdir_dict.pkl")
      pfam_subdir_dict = cPickle.load(f)
      f.close()
      fam = os.path.split(os.getcwd())[1]
      if fam not in pfam_subdir_dict:
        print "Unknown PFAM family %s, exiting..." % fam
        sys.exit(1)
      seq_path_str = ' '.join(list(pfam_subdir_dict[fam]))
      os.system('cat %s > gufmsa.seeds.fa' % seq_path_str)
    else:
      # remove redundant sequences

      print 'making %s unique at %.4f proportion identical...' % (fname, propid),
      sys.stdout.flush()
      
      cmd = '%s gufmsa -alignfile %s -percentid %f &> /dev/null' % (uniqueseq, fname, propid)
      os.system(cmd)    

      handle = open('gufmsa.a2m', 'r')
      outf = open('gufmsa.seeds.fa', 'w')
      
      seq_count = 0
      
      for x in SeqIO.parse(handle, 'fasta'):
          seq_count += 1
      
          header = x.description
          sequence = x.seq.tostring()
          
          print >>outf, '>%s' % header
          
          sequence = sequence.replace('-', '')
          sequence = sequence.replace('.', '')
          sequence = sequence.upper()
          
          print >>outf, sequence
      
      outf.close()
      handle.close()

      print 'done; %d sequences remaining.' % seq_count
    
    # blast non-redundant sequences
    
    print 'blasting remaining sequences against %s...' % blastdb,
    sys.stdout.flush()
    
    result_handle, error_handle = NCBIStandalone.blastall(blastall, 'blastp', blastdb, 'gufmsa.seeds.fa', expectation=ecutoff, descriptions=10000, alignments=10000, filter=False)
    
    blasthit_ids = set([])
    
    for blast_record in NCBIXML.parse(result_handle):
        for alignment in blast_record.alignments:
            blasthit_ids.add(alignment.hit_id)
                
    outf = open('gufmsa.blasthits.ids', 'w')
    
    for x in blasthit_ids:
        print >> outf, x
        
    outf.close()
                
    print 'done; %d potential homologs found.' % len(blasthit_ids)
    
    # retrieve full-length sequences from the db
    
    print 'Retrieving full-length potential homolog sequences...',
    sys.stdout.flush()
    
    cmd = '%s -i gufmsa.blasthits.ids -d %s > TheUniverse.fa' % (fastacmd, blastdb)
    os.system(cmd)
    
    print 'done; universe is in universe.fa'
    
    

def main():
    
    parser = OptionParser(usage='Usage: [-p prop_id] [-d blastdb] [-e ecutoff] alnfile')
    
    parser.add_option('-p', '--propid', dest='propid', default=0.3, help='proportion-identical sequences to remove from alnfile before blast, default=0.3')
    parser.add_option('-d', '--blastdb', dest='blastdb', default='/clusterfs/ohana/external/UniProt/current/protein', help='blast database to search against, default=UniProt')
    parser.add_option('-e', '--ecutoff', dest='ecutoff', default=100.0, help='e-value cutoff for blast results, default=100.0')
    parser.add_option('-s', '--use_subfamily_seeds', dest='use_subfamily_seeds',
                action="store_true", default=False,
                help="Whether to use seeds from previously defined subfamilies")
    
    (options, args) = parser.parse_args()
    
    if len(args) != 1:
        print 'ERROR: Input alignment filename required.'
        sys.exit(1)
        
    fname = args[0]
    
    try:
        propid = float(options.propid)
    except ValueError:
        parser.error('nice try; -p must be a number')
        
    if propid < 0:
        parser.error('come on; -p must be at least zero')

    try:
        ecutoff = float(options.ecutoff)
    except ValueError:
        parser.error('nice try; -e must be a number')
        
    if ecutoff < 0:
        parser.error('come on; -e must be at least zero')

        
    blastdb = options.blastdb

    use_subfamily_seeds = options.use_subfamily_seeds
    
    getUniverseFromMSA(fname, propid, blastdb, ecutoff, use_subfamily_seeds)
    return 0



if __name__ == '__main__':
    main()
    
