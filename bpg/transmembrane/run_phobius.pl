#!/usr/local/bin/perl -w


sub consensus {

        # Given the name of the book, return the consensus sequence in FASTA
        # format.

	($book) = @_;
	if (! ($book =~ /(bpg\d{6})/)) {
	    print STDERR "$book is not bpg file.\n";

            # If file with the name other than bpg given, do nothing.
	    next;
        } 
        $group = substr( $book, 0, 6 );
	$alignment_file = "/usr/local/pfacts/$group/$book/user/$book.a2m";
	use File::Temp qw/ tempfile tempdir /;
	($fh_1, $temp_1) = tempfile (DIR => "/home/bpg/tmp" );
	my $cmd = "alignment_consensus -c 0 $alignment_file > $temp_1";

	system( $cmd ) == 0 or die "alignment_consensus didn't work";
	
        # Eliminate gaps in the sequence.
        # Add header.
	
	open ($fh_1, "$temp_1")
		|| die "Cannot open the file \"$temp_1\": $! .\n";
	@data = <$fh_1>;
	close $fh_1;
	foreach (@data) {s/\-//g}
	unshift @data, "\>$book\n";
	my $cons_seq = join("", @data);

        # Run phobius and parse the output.
	my @result = `curl -F "protseq=$cons_seq" -F "format=nog" http://phobius.binf.ku.dk/cgi-bin/predict.pl 2> phobius_result.err`;
	
	if (! @result) {
		die "No Phobius result returned.\n";	
	}

	my @output = ();
	my $fname = "";	
	foreach $line (@result) {
                	if ($line =~ /ID\s+(\S+)/) {
                	        $fname = $1;
                	} elsif ($line =~ /\bTRANSMEM\s+(\d+)\s+(\d+)/) {
                	        push @output, [$fname, "TM", $1, $2];
                	} elsif ($line =~ /\bSIGNAL\s+(\d+)\s+(\d+)/) {
                	        push @output, [$fname, "SP", $1, $2];
        	        }
	}

	foreach $region (@output) {
        	$region = join "\t", @$region;
        	print "$region\n";
	}

    }


##Main part of the program.
##When only the program name is typed in, return the explanation of the program.
if (! @ARGV) {
print <<EOF;

Given files containing lists of bpg file names, return their TM SP regions
to the standard output.
If consensus sequences of those bpg sequences contain any gaps, gaps would
be eliminated when they are run on phobius.

./tmsp.pl filename1 filename2 ...


EOF
}

##After given the appropriate output file name, run main program for each file
##and print the progress.

foreach $file (@ARGV) {consensus($file);}
