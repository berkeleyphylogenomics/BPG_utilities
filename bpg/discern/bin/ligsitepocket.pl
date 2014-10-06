#! /usr/bin/perl -w
#
use lib '/usr/lib/perl5/site_perl';
use lib '/usr/lib/perl5/5.8.3';
use lib '/usr/lib/perl5/site_perl/5.8.3';

#use lib '/home/sriram_s/perl/lib/perl5/site_perl';
#use lib '/home/sriram_s/perl/lib/perl5/5.8.7';
#use lib '/home/sriram_s/perl/lib/perl5/site_perl/5.8.7';
#use lib '/sw/lib/perl5/5.8.6';

use Bio::Structure::IO;
use Bio::Structure::Atom;

$|=1;
my %features = ();
my %labelmap = ();
my @coords=();
# Positions in the 3d structure
my @positions=();
my $cutoff 		= 10;
my $pockets = 3;

my $data_dir		= "/clusterfs/ohana/software/discern/data";
$data_dir		= $ENV{'DATA_DIR'} if defined $ENV{'DATA_DIR'};
my $aa_file		= "aa.txt";
$aa_file		= $data_dir."/".$aa_file;

my %aaMap	= readMap($aa_file);

if (@ARGV < 1){
	printUsage ();
}
	
my $pdb =  $ARGV[0];
my $lcs =  $ARGV[1];
my ($structId, $chainId) = split (/:/,$pdb);
if ( !defined $chainId || $chainId eq "" ){
	$chainId = "default" ;
}

$chainId = uc($chainId);

my @pocketcoords = ();
my %pocketmap = ();

readlcs ( $structId, $pockets, \@pocketcoords);
#print "pockets\n";
#foreach my $pocket (@pocketcoords){ 
#	print $pocket->x().",".$pocket->y().",".$pocket->z()."\n";
#}
readpdb ($structId, $chainId, \@pocketcoords, $cutoff, \%pocketmap);
foreach my $key (sort {$a<=>$b} keys %pocketmap){
	print $key."\t".$pocketmap{$key}."\n";
}

sub readlcs {
	my $structId = $_[0];
	my $pockets  = $_[1];

	open (F,"$structId.lcs");

	my $count = 0;
	while (<F>){ 
		
		chomp;
		my @tmp = split(/\s+/);
		$atom = Bio::Structure::Atom->new( -id  => 'junk');
		$atom->x($tmp[$#tmp-2]);
		$atom->y($tmp[$#tmp-1]);
		$atom->z($tmp[$#tmp]);
		$coords[$count]  = $atom;
		$count ++;
		if ( $count >= $pockets) {
			last;
		}
	}
	

	@{$_[2]} = @coords;
}

sub readpdb {
	my $structId = $_[0];
	my $chainId = $_[1];
	my @pockets = @{$_[2]};
	my $cutoff = $_[3];
	my $structio 	= Bio::Structure::IO -> new (-file 	=> "$structId.pdb", -format=>'pdb');
	my $struct 		= $structio -> next_structure;
	my %pocketmap = ();

	for my $chain ($struct -> get_chains ){
		my $chainid = $chain->id;

		if($chainid eq $chainId){
				for my $res ($struct->get_residues($chain)){
						my $resid  = $res->id;
						my ($residue,$pos)=split (/-/,$resid);
#			print "Examining residues $residue at $pos\n";
# Checking for strange residues or HETATM records.
						next if !defined $aaMap{$residue};
						next if (!isnumeric ($pos));
						my @atoms  = $struct->get_atoms($res);
						foreach my $atom (@atoms){
								if ($atom->pdb_atomname =~ /CA/){
										$coords		= $atom;
								}
						}

						my $pocket = -1;
						my $dist = $cutoff;
#						print "pos=$pos\t";
						for (my $i = 0; $i < @pockets; $i++) {
								my $tmpdist = dist ($coords, $pockets[$i]);
#								print "dist$i=$tmpdist\t";
								if ($tmpdist < $dist) { 
										$dist = $tmpdist;
										$pocket = $i;
								}
						}
#						print "\n";

						$pocketmap {$pos} = $pocket;
				}
		}
	}
	%{$_[4]} = %pocketmap;
}


sub dist{
	my ($atom1,$atom2)=@_;
	my ($x1,$y1,$z1)=($atom1->x(),$atom1->y(), $atom1->z());
	my ($x2,$y2,$z2)=($atom2->x(),$atom2->y(), $atom2->z());
	return sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
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

sub isnumeric{
	my $x = $_[0];
	
	return 1 if ($x =~ m/^-?\d+$/);
	return 0;
}

sub printUsage { 
	print<<ENDL;
	$0 <PDB id:Chain id> <LCS file>
ENDL
	die;
}
