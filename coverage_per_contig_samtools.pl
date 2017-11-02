#!/usr/bin/perl
#Eric Morrison
#11-17-16
#coverage_per_contig.pl [input coverage.txt from samtools with -a]
#This script calculates a coverage per contig/node from assembly by totaling depth of coverage at ea. position and dividing by number of positions (i.e. contig length). The appropriate input file is generated by samtools depth -a [bamfile]
# Output is node name\tcontig length (bp)\tcoverage. Genome wide values are calculated by totaling all depths and all positions.

use strict;
use warnings;

my $in = $ARGV[0];
open (IN, "$in") || die "Can't open input.\n";

chomp(my @in = <IN>);

if(scalar@in == 1)
	{
	$in[0] =~ s/\r|\r\n|\n/\n/g;
	@in = split("\n", $in[0]);
	}

#Hashes for node values
my %nodesCov;
my %nodesLen;

#Scalars for genome-wide values
my $genCov;
my $genLen;

foreach my $line (@in)
	{
	my @line = split("\t", $line);
	
	$nodesCov{$line[0]} += $line[2];
	$nodesLen{$line[0]}++;
	
	$genCov += $line[2];
	$genLen++;
	}

foreach my $node (sort {$nodesLen{$b} <=> $nodesLen{$a} } (keys %nodesLen) )
	{
	print $node, "\t", $nodesLen{$node}, "\t", $nodesCov{$node}/$nodesLen{$node}, "\n";
	}
	
print "genome\t", $genLen, "\t", $genCov/$genLen, "\n";

