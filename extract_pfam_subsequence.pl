# This script will extract the subsequence of pfam domains, and only chose one
# subsequence if multiple domains match the exact same region
use Getopt::Std;
getopts ("i:o:a:");
die "Usage: -a original fasta file -i input file -o output file" unless (defined $opt_i && defined $opt_o);
open (R, "$opt_i") || die "Cannot open $opt_i";
open (O, ">$opt_o");
@id_check=();
while (<R>)
{
	next if (/^Query/);
  $mark=0;
	$count=0;
	@line=split /,/;
	$head=`grep ">$line[0] " $opt_a`;
	#print $head;
	chomp $head;
	@desc=split/\s+/,$head;
	shift @desc;
	$desc=join (" ",@desc);	
	$id=$line[0]."_".$line[3]."_".$line[4].",".$line[1].",".$desc;
	#`mkdir New_pfam_seed/$id` unless (-e "New_pfam_seed/$id");
	#$file_name="$id.fasta";
  $id_check=$line[0]."_".$line[3]."_".$line[4];
  foreach $id_prev (@id_prev)
  {if ($id_check eq $id_prev) {$mark=1;last;}}
  next if ($mark);
	print O ">$id\n$line[5]";
	push @id_prev,$id_check;
	#print "Finish $id\n";
	
}
close (O);
close (R);

