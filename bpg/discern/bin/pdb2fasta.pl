#! /usr/bin/perl -w

use lib '/usr/lib/perl5/site_perl';
use lib '/usr/lib/perl5/5.8.3';
use lib '/usr/lib/perl5/site_perl/5.8.3';

#use lib '/home/sriram_s/perl/lib/perl5/site_perl';
#use lib '/home/sriram_s/perl/lib/perl5/5.8.7';
#use lib '/home/sriram_s/perl/lib/perl5/site_perl/5.8.7';


use strict;
use Getopt::Long;
use Bio::Structure::IO;
$| = 1;
my $pdbId;
my $chainId;
my $pdb;
my $aaFile		= "aa.txt";
my $dataDir		= "/clusterfs/ohana/software/discern/data";
$dataDir		= $ENV{'DATA_DIR'} if defined $ENV{'DATA_DIR'};
my $verbose 	= 0;

my %charmap  = (
				'A'	=> 0,
				'B'	=> 1,
				'C'	=> 2,
				'D'	=> 3,
				'E'	=> 4,
				'F'	=> 5,
				'G'	=> 6,
				'H' => 7,
				'I' => 8,
				'J' => 9,
				'K'	=>10,
				'L'	=>11,
				'M'	=>12,
				'N'	=>13,
				'O'	=>14,
				'P'	=>15,
				'Q'	=>16,
				'R' =>17,
				'S' =>18,
				'T' =>19,
				'U' =>20,
				'V' =>21,
				'W' =>22,
				'X' =>23,
				'Y' =>24,
				'Z' =>25 );



GetOptions 	( 	"verbose"	=> \$verbose);
if (@ARGV < 1){
	printusage ();
}

$pdb = $ARGV[0];

die if (!defined $pdb);
($pdbId, $chainId) = split (/:/, $pdb);
# Need this to handle PDB structures which have no chains.
$chainId = "default" if (!defined $chainId);
$chainId = "default" if ($chainId eq "");
my $pdbFile = $pdbId."\.pdb";
$aaFile 		= $dataDir."/".$aaFile;
my %aaMap		= readMap($aaFile);
my $index		= 0;

my %resMap		= ();
my $minPos;
my $maxPos;
my $structio 	= Bio::Structure::IO -> new (-file 	=> $pdbFile, -format=>'pdb');
my $struct 		= $structio -> next_structure;
my @tmp = $struct->get_chains;
print "Chains=".scalar @tmp."\n" if ($verbose);

for my $chain ($struct -> get_chains ){
	my $chainid = $chain->id;
	print "Chainid=".$chainid."\n" if ($verbose);
	if($chainid eq $chainId){
		for my $res ($struct->get_residues($chain)){
			my $resid  = $res->id;
			my ($residue,$pos)=split (/-/,$resid);
			# Checking for strange residues or HETATM records.
			next if !defined $aaMap{$residue};
			if (!( $pos =~ /^[0-9]+$/)){
				my @tmp  = split (/\./,$pos);
				my $dec = ($charmap{$tmp[1]}+1)/30;
				$pos = $tmp[0] + $dec;
			}	
#			print $pos."\t".$residue."\n";
			$resMap{$pos}=$aaMap{$residue};
		}
	}
}
#print "minpos=$minPos\tmaxpos=$maxPos\n";

foreach my $key (sort {$a<=>$b} keys %resMap){
	print $resMap{$key};
}



sub readMap{
	my $mapFile = shift (@_);
	my %map = ();
	open (FILE,$mapFile) or die "$0: Error opening $mapFile";	
	my @lines = <FILE>;
	foreach my $line (@lines){
		chomp ($line);
		my @tmp = split (/\s+/,$line);
		$map {$tmp[0]} = $tmp[1];
	}
	close FILE;
	return %map;
}

sub printusage {
	print <<ENDL;
$0 [Options] <pdbid:chain>
Options:
-v	: verbose
ENDL
}
