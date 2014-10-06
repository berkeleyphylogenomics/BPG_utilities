#!/usr/bin/perl
#
# $Author: mm238 $ $Revision: 1.2 $ $Date: 2005/02/15 14:24:03 $
#
# Copyright (C) 2002 MRC and Martin Madera
#
# This program may be freely redistributed and/or modified 
# under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 2 of the 
# License, or (at your option) any later version.
#
use lib "/clusterfs/ohana/software/lib/perl/lib/perl5/site_perl/5.8.8"; # insert your lib directory
use HiddenMarkovModel;
use strict;

my $usage = "Usage: model_convert.pl <*input file> <*output file> <fasta file>\n";
my ($infile, $outfile, $fastafile) = @ARGV;
my $model;

die "$usage\n" unless @ARGV>=2;
$model = HiddenMarkovModel->clever_read($infile);
$model->clever_write($outfile,$fastafile);

