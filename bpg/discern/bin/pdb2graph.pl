#! /usr/bin/perl -w
#
use lib '/usr/lib/perl5/site_perl';
use lib '/usr/lib/perl5/5.8.3';
use lib '/usr/lib/perl5/site_perl/5.8.3';

#use lib '/home/sriram_s/perl/lib/perl5/site_perl';
#use lib '/home/sriram_s/perl/lib/perl5/5.8.7';
#use lib '/home/sriram_s/perl/lib/perl5/site_perl/5.8.7';

use Bio::Structure::IO;

my @coords=();
# Positions in the 3d structure
my @positions=();
my $index		= 0;
my $cutoff 		= 5;

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

$chainId = uc($chainId);

my $structio 	= Bio::Structure::IO -> new (-file 	=> "$structId.pdb", -format=>'pdb');
my $struct 		= $structio -> next_structure;
open (DIST, ">graph-distance.txt");

for my $chain ($struct -> get_chains ){
	my $chainid = $chain->id;
#	print "chainid = " .$chainId."\n";
#	print "Examining ".$chainid."\n";

	if($chainid eq $chainId){
		for my $res ($struct->get_residues($chain)){
			my $resid  = $res->id;
			my ($residue,$pos)=split (/-/,$resid);
	# Checking for strange residues or HETATM records.
			next if !defined $aaMap{$residue};
			next if $pos eq "";
			my @atoms  = $struct->get_atoms($res);
			$coords [$index] = \@atoms;

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

if ( $cutoff > 0 ) {
	for (my $i = 0; $i < @coords; $i++){
		my @atoms1 = @{$coords[$i]};
		for (my $j = $i + 1; $j < @coords; $j++){
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
			if ($mindist < $cutoff){
				$edge [$edges] = "$i\t$j";
				print DIST "$mindist\n";
				$edges ++;
			}
		}
	}
} else {
	for (my $i = 0; $i < @coords; $i++){
		my @atoms1 = @{$coords[$i]};
		for (my $j = $i + 1; $j < @coords; $j++){
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
			if ($mindist < $cutoff){
				$edge [$edges] = "$i\t$j";
				$edges ++;
			}
		}
	}
}
close DIST;


print @coords."\t".$edges."\n";
for (my $i = 0; $i < @coords; $i++){
	print $positions[$i]."\n";
}
for (my $i = 0; $i < $edges; $i++){
	print $edge[$i]."\n";
}
print "\n";

sub dist{
	my ($atom1,$atom2)=@_;
	my ($x1,$y1,$z1)=($atom1->x(),$atom1->y(), $atom1->z());
	my ($x2,$y2,$z2)=($atom2->x(),$atom2->y(), $atom2->z());
	return sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
}



sub printUsage{
	print "Usage:$0 <pdb:chain>\n";
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
