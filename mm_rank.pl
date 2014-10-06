#! /usr/bin/perl

my @w =qw(0.67 6.79);
my $bias = -9.21;
my $f = $ARGV[0];
open (F,$f);
while (<F> ){
	chomp ;
	my @tmp = split (/\s+/);
	my $key = shift @tmp;
	my $score = $bias;
	for (my $i=  0; $i < @tmp; $i++){
		$score += $tmp[$i];	
	}

	print $key."\t".$score."\n";
}

close F;
