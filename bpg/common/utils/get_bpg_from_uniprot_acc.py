import re,sys
from optparse import OptionParser
from pfacts003.phylofacts.models import * 

if len(sys.argv) < 3:
    sys.stderr.write('Usage: %s -i file with uniprot accesions -o output file' % sys.argv[0])
    sys.exit(1)
Parser = OptionParser('usage: sys.argv[0] -i file with uniprot accesions -o output file')
Parser.add_option('-o', '--output',
                dest = 'output',
                default = True,
                help = "bpg accessions of books containing uniprot accession in the input file")
Parser.add_option('-i', '--input',
                dest = 'input',
                default = True,
                help = "file containing uniprot accessions from a given species")
(CmdLineOps, Args) = Parser.parse_args()

f = open (CmdLineOps.input,'r')
w = open (CmdLineOps.output,'w')
#print 'working for %s' % name	
for line in f:
	m = re.match ('(\w+)',line)
	accession = m.group(0)
	#accession = line.rstrip('\n')
	tree_nodes = TreeNode.objects.filter(sequence_header__uniprot__accession = accession).exclude(tree__family__status__exact = 'bad')
	families_all = [tree_node.tree.family for tree_node in tree_nodes]
	families = set (families_all)
	if len (families) ==0:
	   try:
		uniprot_dat_index = UniProtDatIndex.objects.get(uniprot_accession__exact
                                                            = accession)
        	uniprot = uniprot_dat_index.uniprot
        #	sequence_header = _sequence_header(description, sequence,
        #                                   uniprot=uniprot, taxon=uniprot.taxon)
		tree_nodes = TreeNode.objects.filter(sequence_header__uniprot__accession = uniprot.accession).exclude(tree__family__status__exact = 'bad')
		new_families_all = [tree_node.tree.family for tree_node in tree_nodes]
		new_families = set (new_families_all)
		if  len (new_families) ==0:
			w.write ('%s no family found\n' % accession)
		#	print '%s no family found' % accession
		else: 
			w.write ('%s (changed)\n' % accession)
		#	print '%s (changed)' % accession
			for new_family in new_families:
				w.write ('%s %s\n' % (new_family.__str__(),new_family.family_type_id))
				#print '%s %s' % (new_family.__str__(),new_family.family_type_id) 
	   except UniProtDatIndex.DoesNotExist:
		w.write ('%s no family found\n' % accession)
               # print '%s no family found' % accession

	else:
		w.write ('%s\n' % accession)
		#print '%s' % accession
		for family in families:
			w.write ('%s %s\n' % (family.__str__(),family.family_type_id))
                 #       print '%s %s' % (family.__str__(),family.family_type_id)
			#w.write ('%s' % family.__str__())
			#print '%s' % family.__str__()
f.close()
w.close()
