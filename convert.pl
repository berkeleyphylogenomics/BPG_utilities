#!/usr/bin/perl
#
my $usage = "Usage: convert.pl <input file>\n";

my @states;        # state transition probabilities
my @matches;       # match emission probabilities
my @inserts;       # insert emission probabilities
my $M;             # number of positions (match nodes) in the model



##
## read in a SAM 3.0 ASCII file
## 
sub readBlock
{
    local *HANDLE = shift;  # file handle to read from
    my $filename = shift;   # name of the file
    my $cols = shift;       # number of columns to read
    my $rows = shift;       # number of rows to read
    my $array = [];         # array to return
    my $i, @a;

    for $i (1 .. $rows)
    {
	$_ = <HANDLE>;
	@a = split /\s+/;
	die "Parsing error on line $. of $filename!\n$usage\n..."
	    unless @a == $cols;

	push @{$array}, @a;
    };

    return $array;
};

sub readSAM
{
    my $filename = shift;      # name of the file to parse
    my $position = -2;         # current position in the model
    my $last = 0;              # is this the last position?

    my $array=[], @a, $i;

    @states=();
    @matches=();
    @inserts=();

    open( SAM, "$filename" ) 
	or die "Cannot open $filename: $!\n$usage\n...";
    
    while(<SAM>)
    {
	next unless $position >=0 or /^Begin/;

	if( /^Begin/ ) 	{
	    $position=0;
	} elsif ( /^End/ ) {
	    $position++;
	    $last=1;
	} elsif( /^\s*(\S+)/ ) {
	    $position++;
	    die "Parsing error on line $. of $filename!\n..."
		unless $1 == $position;	    
	} else {
	    die "Parsing error on line $. of $filename!\n...";
	};

	push @states,  readBlock( *SAM, $filename, 3, 3 );
	push @matches, readBlock( *SAM, $filename, 5, 4 );
	push @inserts, readBlock( *SAM, $filename, 5, 4 );

	last if $last;
    };
    close SAM;

    $M=$position-1;
};


##
## write a SAM 3.0 ASCII file
## 
sub pad
{
    my $len  = shift;    # pad $text so that the final length is $len
    my $text = shift;    # the text to pad
    return " "x($len-length $text) . $text;
};

sub printBlock
{
    local *HANDLE = shift;  # file handle to write to
    my $array = shift;      # the array to write
    my $cols  = shift;      # ... number of columns
    my $rows  = shift;      # ... number of rows
    my $i, $j;

    for $i (0 .. $rows-1) {
	for $j (0 .. $cols-1) {
	    printf HANDLE "%1.6f ", $array->[$i*$cols+$j];
	};
	print HANDLE "\n";
    };
};

sub writeSAM
{
    my $filename = shift;   # file name to write to
    my $position;           # current position in the model

    open( SAM, ">$filename" ) or die "Cannot open $filename: $!\n...";
    print SAM "MODEL\nalphabet protein\n";
    for $position (0 .. $M+1)
    {
	if( $position == 0 ) { 
	    print SAM "Begin\n";
	} elsif ( $position == $M+1 ) {
	    print SAM "End\n";
	} else {
	    print SAM pad(4,$position)."\n";
	};

	printBlock(*SAM, $states[$position],  3, 3);
	printBlock(*SAM, $matches[$position], 5, 4);
	printBlock(*SAM, $inserts[$position], 5, 4);
    };

    print SAM "\nENDMODEL\n";
    close SAM;
};


##
## read a HMMER save file
##

# the NULE line 
my @NULE = ( 595, -1558, 85, 338, -294, 453, -1158, 197, 249, 902, 
	     -1085, -142, -21, -313, 45, 531, 201, 384, -1998, -644 );

sub score
{
    my $prob = shift;       # convert this probability to a score
    my $nullProb = shift;   # given this null probability
    
    return $prob != 0 ? int( 0.5+1000*log($prob/$nullProb)/log(2.0) ) : "*";
};

sub prob
{
    my $score = shift;      # convert this score to a probability
    my $nullProb = shift;   # given this null probability

    return $score=~ /\*/ ? 0 : $nullProb*exp(log(2)*$score/1000);
};

sub readLine
{
    local *HANDLE = shift;  # file handle to read from
    my $filename = shift;   # name of the file
    my $number = shift;     # number of data expected on the line
    my @array;

    $_ = <HANDLE>;
    @array = split /\s+/;
    shift @array;

    die "Parsing error on line $. of $filename\n..."
	unless @array == $number;

    return @array;
};

sub readHMMER
{
    my $filename = shift;   # name of the file to read from

                            # transition probabilities: $xy = P(X->Y)
    my $dd, $dm, $md, $mm, $mi, $im, $ii, $bm, $me, $bd;
    my $position;           # current position in the model
    my $i, @a, $total;


    @states=();
    @matches=();
    @inserts=();

    open( HMMER, "$filename" ) 
	or die "Cannot open $filename: $!\n$usage\n...";

    ##
    ## header
    ##
    while(<HMMER>) 
    {
	if( /^LENG\s\s(\d+)/ )
	{
	    $M=$1;

	    # "allocate" @states [is this necessary?]
	    for $i (0 .. $M+1)
	    {
		push @states,  [0,0,0, 0,0,0, 0,0,0];
	    };

	    # the "Begin" position values
	    push @matches, [ 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0 ];
	    push @inserts, [ 0.077850, 0.016640, 0.057325, 0.065451, 0.039844,
			     0.062177, 0.024314, 0.060479, 0.064418, 0.082109,
			     0.024547, 0.047302, 0.040276, 0.041960, 0.048771,
			     0.067948, 0.059895, 0.072042, 0.011740, 0.034912 ];
	}
	elsif( /^HMM\s/ )
	{
	    <HMMER>;
	    last;
	};
    };

    ##
    ## "3-line"
    ##
    @a = readLine( *HMMER, $filename, 3 );
    $bd = prob( $a[2], 1 );


    ##
    ## data
    ##
    for $position (1 .. $M)
    {
	##
	## match emissions
	##
	@a = readLine( *HMMER, $filename, 22 );
	shift @a;  # what's this last number ??!
	pop @a;

	for $i (0 .. 19) { $a[$i] = prob( $a[$i], prob($NULE[$i],1/20) ); };
	push @matches, [@a];

	##
	## insert emissions
	##
	@a = readLine( *HMMER, $filename, 21 );
	shift @a;
	for $i (0 .. 19) { $a[$i] = prob( $a[$i], prob($NULE[$i],1/20) ); };
	push @inserts, [@a];

	##
	## state transitions
	##
	@a = readLine( *HMMER, $filename, 10 );
	shift @a;

	for $i (0 .. 8) { $a[$i] = prob( $a[$i], 1 ); };
	($mm, $mi, $md, $im, $ii, $dm, $dd, $bm, $me) = @a;

	if( $position == $M ) 
	{
	    $dd=0; $md=0; 
	    $dm=1; $mm=1; $im=1;
	           $mi=0; $ii=0;
	} 
	elsif ( $position == 1 ) 
	{
	    $total = $bm + $bd;
	    $bm /= $total;
	    $bd /= $total;

	    $states[1]->[4] = $bm;
	    $states[1]->[1] = $bd;
	    $states[1]->[5] = 1;
	};

	# renormalisation
	$total = $mi + $md + $mm;
	if( $total ) {
	    $mi /= $total;
	    $md /= $total;
	    $mm /= $total;
	};
	$total = $dd + $dm;
	if( $total ) {
	    $dd /= $total;
	    $dm /= $total;
	};
	$total = $ii + $im;
	if( $total ) {
	    $ii /= $total;
	    $im /= $total;
	};

	# save the probabilities
	$states[$position+1]->[4] = $mm;
	$states[$position]->[7]   = $mi;
	$states[$position+1]->[1] = $md;
	$states[$position+1]->[5] = $im;
	$states[$position]->[8]   = $ii;
	$states[$position+1]->[3] = $dm;
	$states[$position+1]->[0] = $dd;
    };

    close HMMER;

    ##
    ## the "End" + last position values
    ##
    push @matches, [ 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0 ];
    pop  @inserts;
    push @inserts, [ 0.077850, 0.016640, 0.057325, 0.065451, 0.039844,
		     0.062177, 0.024314, 0.060479, 0.064418, 0.082109,
		     0.024547, 0.047302, 0.040276, 0.041960, 0.048771,
		     0.067948, 0.059895, 0.072042, 0.011740, 0.034912 ];
    push @inserts, [ 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0 ];
};


##
## write a HMMER save file
##
my $header =
    "ALPH  Amino\n".
    "XT      -8455     -4  -1000  -1000  -8455     -4  -8455     -4 \n".
    "NULT      -4  -8455\n".
    "NULE     595  -1558     85    338   -294    453  -1158    197    249    902  -1085   -142    -21   -313     45    531    201    384  -1998   -644 \n".
    "HMM        A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y    \n".
    "         m->m   m->i   m->d   i->m   i->i   d->m   d->d   b->m   m->e\n";

sub printLine 
{
    local *HANDLE = shift;   # file handle to print to

    print HANDLE pad(6, shift);
    while( @_ ) { print HANDLE pad(7, shift); };
    print HANDLE "\n";
};

sub writeHMMER
{
    my $filename = shift;    # file name to read from
                             # transition probabilities: $xy = P(X->Y)
    my $dd, $dm, $md, $mm, $mi, $im, $ii, $bm, $me, $bd, $denom, $total;
    my $i, $j, @a;

    open( HMMER, ">$filename" ) or die "Cannot open $filename: $!\n...";

    ##
    ## header
    ##
    print HMMER "HMMER2.0\n";
    print HMMER "NAME  $filename\n";
    print HMMER "LENG  $M\n";
    print HMMER $header;

    ##
    ## B->D1 line
    ##
    $bd    = $states[1]->[1];
    $denom = $bd + 2*$states[1]->[4];
    printLine( *HMMER, "", score(1-$bd/$denom,1), "*", score($bd/$denom,1) );

    for $i (1 ... $M)
    {
	##
	## match emissions
	##
	@a=();       
	for $j (0 .. 19) {
	    push @a, score( $matches[$i]->[$j], prob($NULE[$j],1/20) );
	};
	printLine( *HMMER, $i, @a, $i );

	##
	## insert emissions
	##
	@a=();
	for $j (0 .. 19) {
	    push @a,  $i==$M ? "*" : score( $inserts[$i]->[$j], prob($NULE[$j],1/20) );
	};
	printLine( *HMMER, "-", @a );


	##
	## state transitions
	##
	if( $i == $M )
	{
	    $mm=$mi=$md=$im=$ii=$dm=$dd=0;
	    $me=1;
	}
	else 
	{
	    $mm = $states[$i+1]->[4];
	    $mi = $states[$i]->[7];
	    $md = $states[$i+1]->[1];
	    $im = $states[$i+1]->[5];
	    $ii = $states[$i]->[8];
	    $dm = $states[$i+1]->[3];
	    $dd = $states[$i+1]->[0];
	    $me = 1/(2*($M-1)-$i);

	    # renormalization
	    $total = $dd + $dm;
	    $dd /= $total;
	    $dm /= $total;
	    $total = $ii + $im;
	    $id /= $total;
	    $im /= $total;
	    $total = $me + $mi + $md + $mm;
	    $me /= $total;
	    $mi /= $total;
	    $md /= $total;
	    $mm /= $total;
	};
	$bm = ($i==1)  ?  $states[1]->[4]  :  $states[1]->[4]/($M-1);
	$bm = $bm / $denom;
	
	printLine( *HMMER, "-", 
		   score($mm, 1), score($mi, 1), score($md, 1),
		   score($im, 1), score($ii, 1), score($dm, 1),
		   score($dd, 1), score($bm, 1), score($me, 1) );
    };

    print HMMER "//\n";
    close HMMER;
};


##
## main program
## 
my $file = $ARGV[0];
die "$usage\n" unless @ARGV>=1;

if( $file =~ s/.mod$// )
{
    # file name ends in .mod => assuming a SAM file

    unless( $file =~ s/\.asc$// )  
    {
	# does not end in .asc.mod => assuming a binary file
	# => create the ascii .asc.mod file as well
	`hmmconvert $file.asc -model_file $file.mod -binary_output 0`;
	wait;
    };
    readSAM("$file.asc.mod");
    writeHMMER("$file.con.hmm");
}
elsif( $file =~ s/\.hmm$// )
{
    # file name ends in .hmm => assuming a HMMER file
    
    readHMMER("$file.hmm");
    writeSAM("$file.con.asc.mod");

    # create the binary .mod file as well
    `hmmconvert $file.con -model_file $file.con.asc.mod -binary_output 1`;
    wait;
}
else
{
    die "The input file must be either *.mod, *.asc.mod or *.hmm!\n$usage\n...";
};


