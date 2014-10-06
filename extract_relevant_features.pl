#! /usr/bin/perl

my %map = (
	1	=> 1,
	4	=> 1
	);

my $f = $ARGV[0];

open (F,$f);

while (<F> ) {
	chomp;
	my @tmp  = split (/\s+/);
	my $out = $tmp[0];
	for (my $i=1; $i<@tmp;$i++){
		if (defined $map{$i}){
			$out=$out."\t".$tmp[$i];	
		}	
	}
	print $out."\n";		
}

close F;
