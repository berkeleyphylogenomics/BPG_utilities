#! /usr/bin/perl -w


my $pdbid = $ARGV[0];
my $bindir="/clusterfs/ohana/software/discern/bin";
my ($j,$k) = split (/:/,$pdbid);
print "$j-$k\n";


system ("$bindir/pdb2fasta.pl $pdbid  > $j.fa");
system ("echo \">$j\" > er; cat  $j.fa >> er ; mv er $j.fa");
system ("mkdir fpdir");
chdir "fpdir";
system ("$bindir/run_fp.sh ../$j.fa");
chdir "..";
system ("$bindir/msa-nospace.pl fpdir/makebook_input_alignment $j >aln.txt");
system ("ln -s fpdir/makebook_input_alignment.nj nj.txt");

open (F,">config.txt");
print F<<ENDL;
msa_file 	aln.txt
tree_file nj.txt
sequence_id $j
structure_id $pdbid
ENDL
close F;

system ("$bindir/pipeline.pl $pdbid  $bindir/parameters.train-all/weights.txt $bindir/parameters.train-all/mean.txt $bindir/parameters.train-all/stddev.txt ");

