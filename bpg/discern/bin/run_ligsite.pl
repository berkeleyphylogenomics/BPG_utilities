#! /usr/bin/perl

my $pdb = $ARGV[0];
my ($j, $chainId) = split (/:/,$pdb);
my $bindir="/clusterfs/ohana/software/discern/bin";
print "lcs  -i $j.pdb -s 1 -n 3 >log \n";
system "lcs  -i $j.pdb -s 1 -n 3 >log";
print "mv pocket_r.pdb $j.lcs\n";
system "mv pocket_r.pdb $j.lcs";
print "$bindir/ligsitepocket.pl $pdb $j.lcs > $j.pocket\n";
system "$bindir/ligsitepocket.pl $pdb $j.lcs > $j.pocket";

