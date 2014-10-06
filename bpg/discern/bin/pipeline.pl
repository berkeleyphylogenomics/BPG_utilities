#! /usr/bin/perl -w
# 
# Input: PDB ID p:q (p=pdb id, q=chain id).
# # This directory must contain the input files: p.pdb, aln.txt (the MSA) and tree.txt (the Tree), config.txt (these are used in INTREPID)
# For the contents of config.txt, refer the input to run_intrepid.pl




my $bindir="/clusterfs/ohana/software/discern/bin";
my $naccessdata="/clusterfs/ohana/software/naccess2.1.1";
my $debug = 1;

if (@ARGV < 4 ){
	print<<ENDL;
	$0 <PDB Id> <Weight file> <Mean parameter file> <Stddev parameter file>
ENDL
	die;
}

my $pdbid = $ARGV[0];
my $weightfile = $ARGV[1];
my $meanfile = $ARGV[2];
my $stddevfile = $ARGV[3];
my $bias = -6.339;

my ($j,$k) = split (/:/,$pdbid);
print "$j-$k\n";

run ("$bindir/dsspcmbi $j.pdb $j.dssp");
run ("$bindir/naccess -r $naccessdata/vdw.radii  -s $naccessdata/standard.data  $j.pdb");
run ("$bindir/run_ligsite.pl $pdbid");

# run intrepid
run ("$bindir/intrepid.pl config.txt > intrepid.log; mv output.aux intrepid.aux");

# Generates a graph file from a PDB file
run ("$bindir/pdb2graph.pl $j:$k > $j.graph1");

# Generates centrality features in centrality.txt
run ("$bindir/analyze-graph $j.graph1 &> analyze.log; mv centrality.txt $j.central");

# Generates a per-residue features file 
run ("$bindir/integrate.pl $j:$k > $j.features");

# Combines features from neighboring residues
# Also outpus neighbors.txt and distance.txt that contains lists of neighbors and their distances
run ("$bindir/enlargefeatures.pl $j:$k > $j.final-features");


run ("$bindir/rank.pl $j.final-features $weightfile $meanfile $stddevfile $bias > $j.rank");


sub run {
	my $cmd = $_[0];
	print "$cmd\n" if ($debug>0);
	system "$cmd";
	print "returned = $?\n" if ($debug>0);
}
