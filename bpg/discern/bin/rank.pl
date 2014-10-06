#! /usr/bin/perl -w

if (@ARGV < 4){
	print<<ENDL;
$0 <feature file> <weight file> <mean file> <std file>
ENDL
	die;
}

$|=1;
my $featurefile = $ARGV[0];
my $debug = 0;
my $wfile = $ARGV[1];
my $avefile = $ARGV[2];
my $stdfile = $ARGV[3];
my $bias = 0;
my $dfile;

$bias =  $ARGV[4] if (@ARGV > 4);
$dfile = $ARGV[5] if (@ARGV > 5);

my $singlesitedim = 48;

my @w = read_weights ($wfile);
my @ave = read_weights ($avefile);
my @stddev = read_weights ($stdfile);
my @d = read_weights ($dfile) if (@ARGV > 5);
classifySimple ($featurefile, \@w, \@ave, \@stddev);



sub read_weights {
	my $f = $_[0];

	my @w = ();
	open (F, $f)  or die "Cannot open $f\n";
	while (<F>) {
		chomp;
		next if (/^#/);
		push (@w, $_);
	}
	close F;
	return @w;
}


sub classifySimple {
	my $f = $_[0];
	my @w = @{$_[1]};
	my @ave = @{$_[2]};
	my @stddev = @{$_[3]};

	open (F, $f)  or die "Cannot open $f\n";
	my $i = 0;
	while (<F>) {
		chomp;
		next if (/^#/);
		my ($x,$y,$z) = split (/\s+/);
		my @features = split (",",$y);
		
		

		for (my $k = 0; $k < @features; $k++) {
			$features[$k] -= $ave[$k];
			$features[$k] /= $stddev[$k];
		}

	
#		die "Number of features+1 != Number of weights" if (@features+1 != @w);
		my $score = 0;
		for (my $j = 0; $j < @features; $j++){
			$score += $w[$j]*$features[$j];
		}
  		$score += $bias;

		if ($debug > 0) {
			for (my $j = 0; $j < @features; $j++){
				print $features[$j]."\t";
			}
			for (my $j = 0; $j < @features; $j++){
				print $w[$j]*$features[$j]."\t";
			}
			print "\n";
		}
		print "$x\t$score\n";
		$i ++;
	}
	close F;
}


sub classifyWithDistanceWeighting {
	my $f = $_[0];
	my @w = @{$_[1]};
	my @d = @{$_[2]};
	my @ave = @{$_[3]};
	my @stddev = @{$_[4]};

	open (F, $f)  or die "Cannot open $f\n";
	my $i = 0;
	while (<F>) {
		chomp;
		next if (/^#/);
		my ($x,$y) = split (/\s+/);
		my @features = split (",",$y);
		my $tmp = $d[$i];
		my @dist = split (/,/,$tmp);
		my @nfeatures;
		for (my $k = 0; $k < $singlesitedim; $k++){
			$nfeatures[$k] = $features[$k];
		}
		$numNN = @features/$singlesitedim - 1;
		$numNN --;
		my @dw = ();
		my $denom = 0;
		for (my $j = 0; $j < $numNN; $j++){
			$dw[$j] = 1/(1+0.01*$dist[$j]);	
#			$dw[$j] = 1;
			$denom += $dw[$j];
		}
		for (my $j = 0; $j < $numNN; $j++){
			$dw[$j] /= $denom;
		}

		for (my $j = 0; $j < $numNN; $j++){
			for (my $k = 0; $k < $singlesitedim; $k++){
				$nfeatures[$singlesitedim+$k] += $dw[$j]*$features[($j+1)*$singlesitedim+$k];
			}
		}
		if ($debug > 0) {
			for (my $j = $singlesitedim; $j <  @nfeatures; $j++) {
				print $nfeatures[$j].",";
			}	
			print "\n";
		}

		for (my $k = 0; $k < 2 * $singlesitedim; $k++) {
			$nfeatures[$k] -= $ave[$k];
			$nfeatures[$k] /= $stddev[$k];
		}
		if ($debug > 0) {
			for (my $j = $singlesitedim; $j <  @nfeatures; $j++) {
				print $nfeatures[$j].",";
			}
			print "\n";
		}

	
#		die "Number of features+1 != Number of weights" if (@features+1 != @w);
		my $score = 0;
		my $score1 = 0;
		for (my $j = 0; $j < @nfeatures; $j++){
			if ($j == $singlesitedim) {
				$score1 = $score;
			}
			$score += $w[$j]*$nfeatures[$j];
		}
		my $score2 = $score - $score1;
  		$score += $bias;

		if ($debug > 0) {
			print "$x\t$score\t$score1\t$score2\n";
		}
		print "$x\t$score\n";
		$i ++;
	}
	close F;
}
