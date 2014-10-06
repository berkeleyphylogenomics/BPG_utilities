#! /usr/bin/perl -w
#use lib '/home/sriram_s/perl/lib/perl5/site_perl';
#use lib '/home/sriram_s/perl/lib/perl5/5.8.7';
#use lib '/home/sriram_s/perl/lib/perl5/site_perl/5.8.7';

use strict;
use Getopt::Long;
$|=1;

use Bio::AlignIO;
my $file = $ARGV[0];
my $id = $ARGV[1];
my $verbose = 0;

GetOptions ( "v"	=> \$verbose);

if (@ARGV < 2){
	print "Usage:$0 <MSA file> <Reference sequence>\n";
	print "Prints out an MSA in which only columns that have no gaps in the reference sequence are retained\n";
	die;
}
my $sequence;

	
my $sequences	= 	Bio::AlignIO->new (-file 	=> $file);
my $aln=$sequences->next_aln();
my @ids=();
foreach my $seq ($aln->each_seq){
	if ($seq->display_id()=~/$id/){
		$sequence= $seq;
		last;
	}
}
if ($verbose){
	print "Sequence:";
	print $sequence->seq()."\n";
}

foreach  my $seq ($aln->each_seq){
	push (@ids, $seq->display_id());
}

my @seqs=();
for (my $i=1; $i<=$aln->length(); $i++){
	if ($sequence->subseq($i,$i) !~ /[A-Z]/
		&& $sequence->subseq($i,$i) !~ /[a-z]/){
		next;
	}
	my $slice	= slice ($aln, $i, $i);
	my $j=0;
	foreach my $seq ($slice->each_seq){
		my $res	= $seq->seq();
		$res =~ tr/a-z/A-Z/;
		$res =~ s/\./-/;
		if ( !defined $seqs[$j]){
			$seqs[$j] = "";
		}
		$seqs[$j]= $seqs[$j].$res;
		$j++;
	}
}

for (my $i=0; $i<@seqs; $i++){	
	print ">".$ids[$i]."\n".$seqs[$i]."\n";
}

sub slice{
	my $aln				= shift;
	my ($start, $end)	= @_;
	my $newAln				= new Bio::SimpleAlign();
	foreach my $seq ($aln->each_seq()){
		my $slice=$seq->subseq($start, $end);
		my $newSeq=Bio::LocatableSeq->new(	-seq=>$slice,
											-id=>$seq->id,
											-start=>$seq->start,
											-end=>$seq->end,
											-alphabet=>"protein"
											);
		$newAln->add_seq($newSeq);
	}
	return $newAln;
}


