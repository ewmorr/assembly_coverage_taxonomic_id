#!/usr/bin/perl
#Eric Morisson
#040117
#get_seq_lens.pl

use strict;
use warnings;

my $fasta = $ARGV[0];
open(FAS, "$fasta") || die "Can't open fasta.\n";
chomp(my@fasta = <FAS>);

my $fas = join(":&:&:&:", @fasta);
my@fas = split(">", $fas);
shift@fas;

foreach my $seq(@fas)
	{
	my @seq = split(":&:&:&:", $seq);
	my $id = shift@seq;
	print $id, "\t", length(join("", @seq)), "\n";
	}