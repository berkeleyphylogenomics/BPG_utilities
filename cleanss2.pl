#! /usr/bin/perl -w

my $f = $ARGV[0];
open (F,$f);

while (<F>) {
	chomp;
	next if ($.==1);
	my @tmp = split (/\_/);
	my @tmp1 = split(/,/,$tmp[$#tmp]);	
	print $tmp[2]."_".$tmp[3]."\t".$tmp1[1]."\t".$tmp1[2]."\t".$tmp1[3]."\n";	
}

close F;
