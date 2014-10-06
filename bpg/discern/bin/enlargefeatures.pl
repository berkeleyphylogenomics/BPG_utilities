#! /usr/bin/perl -w
#
use lib '/usr/lib/perl5/site_perl';
use lib '/usr/lib/perl5/5.8.3';
use lib '/usr/lib/perl5/site_perl/5.8.3';

#use lib '/home/sriram_s/perl/lib/perl5/site_perl';
#use lib '/home/sriram_s/perl/lib/perl5/5.8.7';
#use lib '/home/sriram_s/perl/lib/perl5/site_perl/5.8.7';

use Bio::Structure::IO;

my %features = ();
my %labelmap = ();
my @coords=();
# Positions in the 3d structure
my @positions=();
my $index		= 0;
my $cutoff 		= 5;
my $neighbors  = 10;

my $data_dir		= "/clusterfs/ohana/software/discern/data";
$data_dir		= $ENV{'DATA_DIR'} if defined $ENV{'DATA_DIR'};
my $aa_file		= "aa.txt";
$aa_file		= $data_dir."/".$aa_file;

my %aaMap	= readMap($aa_file);

if (@ARGV < 1){
	printUsage ();
}
	
my $pdb =  $ARGV[0];
my ($structId, $chainId) = split (/:/,$pdb);
if ( !defined $chainId || $chainId eq "" ){
	$chainId = "default" ;
}
else
{
    $chainId = uc($chainId);
}

my $structio 	= Bio::Structure::IO -> new (-file 	=> "$structId.pdb", -format=>'pdb');
my $struct 		= $structio -> next_structure;
my $features =  "$structId.features";
readFeatures ($features);
open (NEIGBHORS, ">neighbors.txt");
open (DIST, ">distance.txt");


for my $chain ($struct -> get_chains ){
	my $chainid = $chain->id;
#	print "chainid = " .$chainId."\n";
#	print "Examining ".$chainid."\n";

	if($chainid eq $chainId){
		for my $res ($struct->get_residues($chain)){
			my $resid  = $res->id;
			my ($residue,$pos)=split (/-/,$resid);
#			print "Examining residues $residue at $pos\n";
	# Checking for strange residues or HETATM records.
			next if !defined $aaMap{$residue};
			next if (!isnumeric ($pos));
			my @atoms  = $struct->get_atoms($res);
			$coords [$index] = \@atoms;

#			print "Adding residues $residue at $pos\n";
			$positions[$index]	= $pos;
			$index++;
		}
	}
}


#for (my $i = 0; $i < @coords; $i++){
#	print $coords[$i]."\n";
#}

my @edge =  ();
my $edges = 0;
for (my $i = 0; $i < @coords; $i++){
	my @atoms1 = @{$coords[$i]};
	my @index = ();
	my @neighbors = ();
	my @distances = ();
	my @neighborpos = ();
	for (my $j = 0; $j < @coords; $j++){
		next if ($i==$j);
		my @atoms2 = @{$coords[$j]};

		my $mindist ;
		for (my $i1 = 0; $i1 < @atoms1; $i1 ++){
			for (my $i2 = 0; $i2 < @atoms2; $i2++){
				my $distance = dist ($atoms1[$i1], $atoms2[$i2]);
				if ( $i1 ==0 && $i2 == 0 ) {
					$mindist = $distance;
				}
				if ($distance < $mindist) {
					$mindist = $distance ;
				}
			}
		}
		push (@distances, $mindist);
		push (@neighbors, $positions[$j]);
	}
	mysort(\@distances, \@index);

	
	my @f = ();
	my @dist = ();
	$f[0] = $features {$positions[$i]};
	push (@neighborpos, $positions[$i]);
	for (my $j = 0 ; $j < @index && $j < $neighbors ; $j++){
		my $tmp = $neighbors[$index[$j]];
		push (@f, $features{$tmp});
		push (@neighborpos, $neighbors[$index[$j]]);
		push (@dist, $distances[$j]);
	}
	for (my $j = @index; $j < $neighbors; $j++){
		push (@f, $features{$positions[$i]});
		push (@neighborpos, $positions[$i]);
		push (@dist, $distances[$j]);
	}

	my $posi = $positions[$i];
	my $newf =  join(",",@f);
	my $label = $labelmap {$posi};
	print $posi."\t".$newf."\t".$label."\n";

	my $neighbors = join (",",@neighborpos);
	print NEIGBHORS "$neighbors\n";

	print DIST join(",",@dist)."\n";

#	print "coord=".$positions[$i]."\t";
#	print join (",",@distances)."\n";
#	print join (",",@index)."\n";
#	print join (",",@neighbors)."\n";

#	die if ($i > 2);
}

close NEIGBHORS;
close DIST;

sub dist{
	my ($atom1,$atom2)=@_;
	my ($x1,$y1,$z1)=($atom1->x(),$atom1->y(), $atom1->z());
	my ($x2,$y2,$z2)=($atom2->x(),$atom2->y(), $atom2->z());
	return sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
}



sub printUsage{
	print "Usage:$0 <pdb:chain> \n";
    die;
}

sub readMap{
	my $mapFile = shift (@_);
	my %map = ();
	open (FILE,$mapFile);	
	my @lines = <FILE>;
	foreach my $line (@lines){
		chomp ($line);
		my @tmp = split (/\s+/,$line);
		$map {$tmp[0]} = $tmp[1];
	}
	close FILE;
	return %map;
}

sub swap {
	my $a = $_[0];
	my $b = $_[1];

	return ($b,$a);
}
	

sub mysort {
	my @a = @{$_[0]};
	my @index =  ();
	for (my $i = 0; $i < @a; $i++){
		$index [$i] = $i;
	}
	for (my $i = 0; $i < @a; $i++){
		for (my $j = $i + 1; $j < @a; $j++){
			if ($a[$i] > $a[$j]){
				($a[$i],$a[$j]) = swap ($a[$i],$a[$j]);

				($index[$i],$index[$j]) = swap ($index[$i],$index[$j]);
			}
		}
	}

	@{$_[1]} = @index;
	@{$_[0]} = @a;
}


sub readFeatures  {
	my $file = $_[0];

	open (FILE, $file) or die "Cannot open $file";
	while (<FILE>){
		chomp;
		next if (/^#/);
		my @tmp = split (/\s+/);
		my $pos  = $tmp[0];
		my $f = $tmp[1];
		my $label = $tmp[2];

		$features {$pos} = $f;
		$labelmap{$pos} = $label;
	}
	close FILE;
}

sub isnumeric{
	my $x = $_[0];
	
	return 1 if ($x =~ m/^-?\d+$/);
	return 0;
}
