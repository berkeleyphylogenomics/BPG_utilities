#! /usr/bin/perl 

my $x = $ARGV[0];
my $y = $ARGV[1];


my $ylfile = $x."_".$y.".yl";
my $ssfile =  $x."_".$y.".ss2";

my %globalmap = ();
my %ylmap  =();
my $ssmap = ();
my $ylfill = "";
my $ssfill = "";

$ylfill = createmap ($ylfile, \%ylmap);
$ssfill = createmap ($ssfile, \%ssmap);


my @t = keys %globalmap;

foreach my $key (keys %globalmap ) { 
	$ssmap {$key} = $ssfill if (!defined $ssmap{$key});
	$ylmap {$key} = $ylfill if (!defined $ylmap{$key});
	print $key."\t".$ylmap{$key}."\t".$ssmap{$key}."\n";
}


sub createmap {
	my $f = $_[0];

	open (F,$f);
	my %map  =();
	my @average = ();
	my $count = 0;
	while (<F>) {
		chomp;
		my @tmp = split (/\s+/);
		my $key = shift @tmp;
		$globalmap {$key} = 1;
		my $value = join ("\t",@tmp);
		$map{$key} = $value;
		if (@average == 0 ) {
			@average = @tmp;
			$count = 1;
		} else {
			for (my $i = 0; $i < @average; $i++){ 
				$average[$i] = $average[$i] + $tmp[$i];
			}
			$count ++;
		}
	}
	close F;
	for (my $i = 0; $i < @average; $i++){ 
		$average[$i] /= $count;
	}
	my $averagestr = join ("\t",@average);
	%{$_[1]} = %map;
	return $averagestr;

}
