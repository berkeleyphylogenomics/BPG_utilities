#! /usr/bin/perl -w
use lib '/usr/lib/perl5/site_perl';
use lib '/usr/lib/perl5/5.8.3';
use lib '/usr/lib/perl5/site_perl/5.8.3';

#use lib '/home/sriram_s/perl/lib/perl5/site_perl';
#use lib '/home/sriram_s/perl/lib/perl5/5.8.7';
#use lib '/home/sriram_s/perl/lib/perl5/site_perl/5.8.7';

# Put all the feature files together
use Bio::Structure::IO;
$| = 1;
my $verbose = 0;
my $data_dir		= "/clusterfs/ohana/software/discern/data";
$data_dir		= $ENV{'DATA_DIR'} if defined $ENV{'DATA_DIR'};
my $aa_file		= "aa.txt";
$aa_file		= $data_dir."/".$aa_file;
my $metadata = 0;

my %aaMap	= readMap($aa_file);
my %atomMap 	= (
					'N'=>1,
					'CA'=>1,
					'C'=>1,
					'O'=>1);

my %serialMap = ();
my %ssMap = (	'H' => "0,0,0,0,0,0,1",
				'B'	=> "0,0,0,0,0,1,0",
				'E'	=> "0,0,0,0,1,0,0",
				'G'	=> "0,0,0,1,0,0,0",
				'I'	=> "0,0,1,0,0,0,0",
				'T'	=> "0,1,0,0,0,0,0",
				'S'	=> "1,0,0,0,0,0,0"
				 );


my %chargedMap  = (
				'H'	=> 1,
				'R' => 1,
				'K' => 1,
				'E'	=> 1,
				'D' => 1);

my %polarMap = (
				'Q'	=> 1,
				'T' => 1,
				'S'	=> 1,
				'N' => 1,
				'C' => 1,
				'Y' => 1
				);

my %hydrophobicMap = (
				'A'	=> 1,
				'F'	=> 1,
				'G' => 1,
				'I'	=> 1,
				'L' => 1,
				'M' => 1,
				'P' => 1,
				'V' => 1,
				'W' => 1 );
				
my $count =  0;
foreach my $key (keys %aaMap){
	my @tmp =  ();
	for (my $i = 0 ; $i < 20; $i++) { 
		$tmp[$i] = 0;
	}
	$tmp [$count] = 1;
	if ($metadata) {
		print "#\t".$key."\t".join(",",@tmp)."\n";
	}
	$serialMap {$key} = join(",",@tmp);
	$count ++;
}

if (@ARGV < 1){
	printusage ();
}

my $dir = $ARGV[0];
my $givencat = 0;
my $catfile;
if (@ARGV > 1 ) {
	$givencat = 1;
	$catfile = $ARGV[1];
}

my $scripthome = "/home/sriram_s/bpgtools/scripts";

chdir "$dir";

my ($structId, $chainId) = split (/:/,$dir);
my $pdbchainId;
if ($verbose){
	print $structId.":".$chainId."\n";
}
if (!defined $chainId || $chainId eq "" ){
	$pdbchainId = "default";
} else {
	$pdbchainId = uc($chainId);
}


#my $seqfile = "$structId-$pdbchainId.aux";
my $seqfile = "intrepid.aux";
my $sdpfile = "intrepid-sdp.aux";
my $pdbfile = "$structId.pdb";
my $dsspfile = "$structId.dssp";
my $rsafile = "$structId.rsa";
my $pocketfile = "$structId.pocket";
my $sprefile = "$structId-$pdbchainId.spre";
my $dmrefile = "$structId-$pdbchainId.dmre";
my $centralityfile = "$structId.central";
my %catmap =  ();
my %ncatmap = ();
my %dsspmap = ();
my %rsamap = ();
my %seqmap = ();
my %stringmap = ();
my %spremap = ();
my %dmremap = ();
my %nummap = ();

$stringmap {"structval"} = "Catalytic neighborhood, Centrality";
$nummap{"structval"} = "28 - centrality";

readpdb ($pdbfile);
readdssp ($dsspfile);
readrsa ($rsafile);
readseq ($seqfile);
#readsdp ($sdpfile);
readpocket ($pocketfile);
%centralitymap = readcat ($centralityfile);
if ($givencat) {
	%catmap = readcat ($catfile);
}
#%ncatmap = readcat ($ncfile);
#%spremap = readjmol ($sprefile);
#%dmremap = readjmol ($dmrefile);

if ($verbose){ 
print "pdbmap\n";
printmap (\%pdbmap);
print "dsspmap\n";
printmap (\%dsspmap);
print "rsamap\n";
printmap (\%rsamap);
print "catmap\n";
printmap (\%catmap);
}
$count = 0 ;

foreach my $key (sort {$a<=>$b}  keys %pdbmap){
	if ( $count == 0 ) {
		my $string = "";
		if ($metadata ){
			my $string1 = "#".$nummap{"pdb"}.",".$nummap{"seq"}.",".$nummap{"structval"}.",".$nummap{"dssp"}.",".$nummap{"rsa"}.",".$nummap{"pocket"}."\n";
			print $string1;
		}
	}
	my $cat = -1;
	$count ++;

	$cat = 1 if (defined $catmap{$key});
	print $key."\n" if ($verbose);

	print $key."\t".$pdbmap{$key}.",".$seqmap{$key}.",".$centralitymap{$key}.",".$dsspmap{$key}.",".$rsamap{$key}.",".$pocketmap{$key}."\t".$cat."\n";
}


chdir "..";

sub readcat {
	my $file = $_[0];
	my %map = ();

	open (FILE, $file) or die "Error in $0: Cannot open $file";

	print "readcat \n" if $verbose;
	while (<FILE>){
		chomp;
		my @tmp = split (/\s+/);
		my $key = shift (@tmp);
		my $value = join (",",@tmp);
		$map {$key} = $value;
		print "$key\t$value\n" if $verbose;
	}
	return %map;
}


sub readdssp {
	my $file = $_[0];

	open (FILE, $file) or die "Error in $0: Cannot open $file";
	
	print "readdssp\n" if $verbose;
	my $flag = 0;
	while (<FILE>){
		chomp;
		$_ =~s/^\s+//g;
		if ($verbose){
			print ;
			print "\n";
		}
		if (/^#/) {
			$flag = 1;
			next;
		}
		next if (!$flag);
		my @tmp = split (/\s+/);
		my $key = $tmp[1];
		my $value = "";
		my $chain = "";
		if ($chainId eq "") {
			$value = $tmp[3];
		} else {
			$value = $tmp[4];
			$chain = $tmp[2];
			next if ($chain ne $chainId);
		}
		
		if (!defined $ssMap {$value} ) {
			$value = "0,0,0,0,0,0,0";
		} else {
			$value = $ssMap {$value};
		}
		if ($verbose ){
			print "chain=$chain\tvalue=$value\n";
		}
		if (isnumeric ($key)){
			$dsspmap {$key} = $value;
			print "key $key is numeric\n" if ($verbose);
		}
	}

	
	foreach my $key (sort {$a<=>$b}  keys %pdbmap){
		$dsspmap {$key} = "0,0,0,0,0,0,0" if (!defined $dsspmap{$key});
	}

	$stringmap {"dssp"} = " Secondary Structure";
	$nummap {"dssp"} = "29-35:SS";
}

sub readrsa  {
	my $file =  $_[0];

	open (FILE, $file) or die "Error in $0: Cannot open $file";
	
	my $flag = 0;
	print "readrsa\n" if $verbose;
	while (<FILE>){
		chomp;
		next if ( ! /^RES/);
		my @tmp = split (/\s+/);
		shift @tmp;shift @tmp;

		my $chain ;
		if ( $chainId  ne "" ) {
			$chain = shift @tmp;
			next if ($chain ne $chainId);
		}
		my $key = shift @tmp;
		my $value = join (",",@tmp);
		if (isnumeric ($key)){
			$rsamap {$key} = $value;
			print "key $key is numeric\n" if ($verbose);
		}
	}
	my @null = ();
	for (my $i = 0 ; $i < 10; $i++){
		$null[$i] = "-1";
	}
	my $nullstr = join (",", @null);
	
	foreach my $key (sort {$a<=>$b}  keys %pdbmap){
		$rsamap {$key} = $nullstr if (!defined $rsamap{$key});
	}
	$stringmap {"rsa"} = "RSA features";
	$nummap {"rsa"} = "36-45: RSA features";
}


sub readpdb {
	my $file = $_[0];
	my $structio 	= Bio::Structure::IO -> new (-file 	=> $pdbfile, -format=>'pdb');
	my $struct 		= $structio -> next_structure;
	for my $chain ($struct -> get_chains ){
		my $chainid = $chain->id;

		if($chainid eq $pdbchainId){
			print "Found chain:$chainid\n" if $verbose;
			for my $res ($struct->get_residues($chain)){
				my $resid  = $res->id;
				my ($residue,$pos)=split (/-/,$resid);
				next if !($pos=~ m/^-?\d+$/ );
				print "Residue $residue at $pos \n" if ($verbose);
# Checking for strange residues or HETATM records.
				next if !defined $aaMap{$residue};
				my @atoms  = $struct->get_atoms($res);
				my $bfactor = 0;
				foreach my $atom (@atoms){
					$bfactor += $atom->tempfactor();
				}
				$bfactor /= @atoms;
				my $tmpres  = $aaMap {$residue};
				my $charged = 0;
				my $polar  = 0;
				my $hydrophobic = 0;

				$charged = 1 if (defined $chargedMap {$tmpres} ) ;
				$polar = 1 if (defined $polarMap{$tmpres});
				$hydrophobic = 1 if (defined $hydrophobicMap{$tmpres});

				
				
				$pdbmap {$pos} = $serialMap{$residue}.",".$charged.",".$polar.",".$hydrophobic.",".$bfactor;
				
			}
		}
	}
	$stringmap{"pdb"} = "Residue type, B-factor";
	$nummap{"pdb"} = "1-20: Residue type, 21-23: C/P/H, 24:B-factor";
}


sub readpocket {
	my $file = $_[0];
	print "readpocket\n" if ($verbose);
	open (FILE, $file) or die "Cannot open $file";

	while (<FILE>){
		if (/^#/){
			next;
		}

		chomp;
		my @tmp = split (/\s+/);
		my $key = $tmp[0];
		my $value = $tmp[1];
		my @val = (0,0,0);
	#	print join (",",@val)."\n";
		if ($value >= 0 && $value < @val){
			$val[$value] = 1;
		}
		$pocketmap {$key} = join(",",@val);
	}
	$stringmap {"pocket"} = "Pocket";
	$nummap {"pocket"} = "46-48:pocket";
}


sub readseq {
	my $file = $_[0];
	print "readseq\n" if ($verbose);

	open (FILE, $file) or die "Cannot open $file";

	while (<FILE>){
		if (/^#/){
			next;
		}

		chomp;
		my @tmp = split (/\|/);
		my $key =  $tmp[1];
		my $value = $tmp[5].",".$tmp[6].",".$tmp[9];
	#	print "$key\t $value\n";
		$seqmap {$key} = $value;
	}
	foreach my $key (sort {$a<=>$b}  keys %pdbmap){
		$seqmap {$key} = "0,0,0" if (!defined $seqmap{$key});
	}
	
	$stringmap {"seq"} = "Sequence conservation:1,2,3";
	$nummap{"seq"} = "25-27: int-lo,int-js,global-js";
}


sub readsdp {
	my $file = $_[0];
	print "readseq\n" if ($verbose);

	open (FILE, $file) or die "Cannot open $file";

	while (<FILE>){
		if (/^#/){
			next;
		}

		chomp;
		my @tmp = split (/\|/);
		my $key =  $tmp[1];
		my $value = $tmp[8];
	#	print "$key\t $value\n";
		$sdpmap {$key} = $value;
	}
	foreach my $key (sort {$a<=>$b}  keys %pdbmap){
		$sdpmap {$key} = "0" if (!defined $sdpmap{$key});
	}
	
	$stringmap {"sdp"} = "SDP";
	$nummap {"sdp"} = "43:intrepid-spec";
}



sub readjmol {
	my $file = $_[0];
	print "readseq\n" if ($verbose);

	open (FILE, $file) or die "Cannot open $file";

	while (<FILE>){
		if (/^#/){
			next;
		}

		chomp;
		my @tmp = split (/\|/);
		my $key =  $tmp[1];
		my $value = $tmp[5];
		$map {$key} = $value;
	}
	foreach my $key (sort {$a<=>$b}  keys %pdbmap){
		$map {$key} = "0" if (!defined $map{$key});
	}
	
	$stringmap {"re"} = "RE";
	return %map;
}




# readMap	- Loads a map of key-value pairs from a file.
#				Each key-value pair must be separated by newlines.
#				The key must be separated from the value by whitespaces.
# Param		- The file containing the key-value pairs.
# Return	- The map
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



sub printusage {
	print <<ENDL;
$0 <pdb:chain id>
ENDL
	die;
}

sub isnumeric{
	my $x = $_[0];
	
	return 1 if ($x =~ m/^-?\d+$/);
	return 0;
}

sub	printmap {
	my %map = %{$_[0]};
	foreach my $key (sort {$a<=>$b} keys %map){
		print $key.":".$map{$key}."\n";
	}
}
