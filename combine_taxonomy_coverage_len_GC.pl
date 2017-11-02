#!/usr/bin/perl
#Eric Morrison
#03-30-16
#combine_taxonomy_coverage_len_GC.pl
#This script combines coverage estimates and sequence length from bed tools and coverage_per_contig.pl, 
# taxonomy information from taxdump_plast_parse.pl, and calculates GC content for each sequence
# in the source fasta file, then returns this information in a single output file.

use strict;
use warnings;

sub usage(){
	print STDERR q(
Usage: combine_taxonomy_coverage_len_GC.pl [taxonomy file from taxdum_blast_parse.pl] [coverage_by_sequence] [contigs fasta] [output.txt]

This script combines coverage estimates from bed tools and coverage_per_contig.pl, taxonomy information from taxdump_plast_parse.pl, and calculates GC content and length for each sequence in the source fasta file, then returns this information in a single output file.

The output is equivalent to a parsed JSON file from J. Sevigny's blobplot pipeline.

);
exit;
}
if(scalar(@ARGV) == 0 || $ARGV[0] eq "-h")
	{
	&usage();
	}

my $tax = $ARGV[0];
my $cov = $ARGV[1];
my $fas = $ARGV[2];
my $out = $ARGV[3];

open(TAX, "$tax") || die "Can't open blob.\n";
open(COV, "$cov") || die "Can't open coverage.\n";
open(FAS, "$fas") || die "Can't open fasta.\n";
open(OUT, ">$out") || die "Can't open output.\n";

chomp(my @tax = <TAX>);
chomp(my @cov = <COV>);
chomp(my @fas = <FAS>);

if(scalar@tax == 1)
	{
	$tax[0] =~ s/\r|\r\n|\n/\n/g;
	@tax = split("\n", $tax[0]);
	}
if(scalar@cov == 1)
	{
	print "linefeed\n";
	$cov[0] =~ s/\r|\r\n|\n/\n/g;
	@cov = split("\n", $cov[0]);
	}

#process fasta
my $fasta = join(":&:&:&:", @fas);
@fas = split(">", $fasta);
shift@fas;
my %gc;
my %len;
foreach my $seq(@fas)
	{
	my @seq = split(":&:&:&:", $seq);
	my $id = shift@seq;
	$seq = join("", @seq);
	$len{$id} = length($seq);
	my $gcs = () = $seq =~ /G|C/gi; #get number occurences G or C
	$gc{$id} = $gcs/$len{$id};
	#print $len{$id}, "\n";
	}

#Process taxa info to hash
my %tax;
foreach my $taxln (@tax)
	{
	my @taxln = split("\t", $taxln);
	my $taxID = shift@taxln;
	my $spp;
	foreach my $tax (@taxln)
		{
		$tax =~ s/^\s|\s$//g; #Remove leading and trailing whitespaces
		$tax =~ s/\s/./g; #Replace remaining whitespaces with dots
		$tax{$taxID} .= $tax."\t";
		$spp = $tax;
		}
	if(@taxln < 7)	
		{
		for(my $i = @taxln; $i < 7; $i++)
			{
			$tax{$taxID} .= $spp."\t";
			}
		}
	}

#process coverage
my %cov;
foreach my $covln (@cov)
	{
	my @covln = split("\t", $covln);
	$cov{$covln[0]} = $covln[2];
	}

print OUT "NodeName\tDomain\tKingdom\tPhylum\tOrder\tFamily\tGenus\tSpecies\tCoverage\tLength\tGC\n";
foreach my $node (sort{$len{$b} <=> $len{$a}} (keys %len) )
	{
	#print $node, "\n";
	if(defined($cov{$node}) == 0)
		{
		$cov{$node} = "NA";
		}
	if(defined($tax{$node}) == 0)
		{
		$tax{$node} = "no-hit\t" x 7;
		}
	print OUT $node, "\t", $tax{$node};
	if($cov{$node} eq "NA")
		{
		print OUT $cov{$node};
		}else{
		printf OUT "%.3f", $cov{$node};
		}
	print OUT "\t", $len{$node}, "\t"; 
	printf OUT "%.4f", $gc{$node}; print OUT "\n";
	}



	
	
	
	
	
