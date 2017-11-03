#!/usr/bin/perl
#Eric Morrison
#03/29/17
#taxdump_blast_parse.pl [blast file] [nodes.dmp] [names.dmp] [output.txt]
#Blast file should have -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' or at least -outfmt '6 qseqid staxids'
#Taxonomy assignments are made based on staxids and comparison to names.dmp and nodes.dmp

use strict;
use warnings;

sub usage(){
	print STDERR q(
Usage: taxdump_blast_parse.pl [blast file] [nodes.dmp] [names.dmp] [output.txt]

Blast file should have -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' or at least -outfmt '6 qseqid staxids'

Taxonomy assignments are made based on staxids and comparison to names.dmp and nodes.dmp from NCBI taxdump.

);
exit;
}
if(scalar(@ARGV) == 0 || $ARGV[0] eq "-h")
	{
	&usage();
	}

my $blast = $ARGV[0];
my $nodes = $ARGV[1];
my $names = $ARGV[2];
my $out = $ARGV[3];

open(BLS, "$blast") || die "Can't open blast file.\n";
open(NOD, "$nodes") || die "Can't open nodes files.\n";
open(NAM, "$names") || die "Can't open names file.\n";
open(OUT, ">$out") || die "Can't open output.\n";


chomp(my @blast = <BLS>);
$/ = "\t\|\n"; #taxdump line feeds
chomp(my @nodes = <NOD>);
chomp(my @names = <NAM>);

my %nodes;
foreach my $nod (@nodes) #process nodes.dmp to hash
		{
		my @nod = split('\t\|\t', $nod);
		$nodes{$nod[0]} = [@nod]
		}

my %names;
foreach my $nam (@names) #process scientific names from names.dmp to a hash referenced by unique id
	{
	my@nam = split('\t\|\t', $nam);
	if($nam[3] eq 'scientific name')
		{
		#if($nam[0] =~ /100272/){print "found $nam[0]\n"};
		$names{$nam[0]} = [@nam];
		}
	}

print OUT "NodeName\tDomain\tKingdom\tPhylum\tOrder\tFamily\tGenus\tSpecies\n";

my%doneAlready;
foreach my $bls (@blast)
	{
	my @bls = split("\t", $bls);
	if(defined($bls[0]) == 0)
		{
		next;
		}
	if($bls =~ /.*(vector|plasmid).*/i)#vector sequences including cloning, shuttle, integrative, expression are problematic
		{	#plasmid sequences also seem to cause trouble. These technical sequences do not have proper taxonomy.
		print "vector\t", $bls[0], "\n";
		next;
		}
	if($bls[1] =~ /(\d+);\d+/) #There are some non-conforming taxids that have two ids split by a semicolon
		{						#It appears that in some cases the first matches the correct taxonomy but in others the second does
		#print $bls[0], "\n";
		$bls[1] = $1;			#However, this is a difference at the species level so either ID should work for higher taxonomy 
		}
		if(defined($nodes{$bls[1]}) == 0 || defined($names{$bls[1]}) == 0 )
		{
		print "no_dmp_entry\t", $bls[0], "\t$bls[1]\t", "\n";
		next;
		}
if($bls[1] =~ /\D+/) #If species ID is not digits skip to next
		{
		print "not_digits\t", $bls[0], "\n";
		next;
		}
	if(defined($bls[1]) != 1)
		{
		print "no_sp_id\t", $bls[0], "\n";
		next;
		}
	if($bls[1] == 155900)
		{
		print "uncultured organism.\n";
		next;
		}

	if(defined($doneAlready{$bls[0]}) == 1) #only consider top blast hit
		{
		next;
		}
	$doneAlready{$bls[0]} = 1;



	
	#initialize by species tax id for searching names.dmp and nodes.dmp
	my @taxID;
	$taxID[0] = $bls[1];
	my @nodNam;
	
	#iteratively search %nodes and %names to find taxonomy for hits 
	my $killWhile = 0; #This variable prevents an infinite while loop in the case of unexpected non-conforming taxonomic info that is not already parsed in this script
	while(defined($nodNam[0]) == 0 || $nodNam[0] ne "superkingdom")
		{
		if(defined($nodes{$taxID[0]}) == 0 || defined($names{$taxID[0]}) == 0 )
			{
			print "There is a problem with the taxonomy information for $bls[0]. Species ID is undefined. Please review blast output. Exiting without completing parsing.\n";
			exit;
			}
		$killWhile++;
		if($killWhile >= 50) # Change this value to allow more or less while loop iterations
			{
			print "\nThere is a problem with the taxonomy information for $bls[0]. Taxonomic lineage has greater than 50 entries.\nPlease review blast output or increase the allowed number of iterations in this script (L126). Exiting without completing parsing.\n\n";
			exit;
			}

		unshift(@nodNam, $nodes{$taxID[0]}[2]);#This contains the taxonomic level (e.g. "species") associated with each ID number
		unshift(@taxID, $nodes{$taxID[0]}[1]);#This contains ID number. Add to tax ID array second so this doesn't index taxID array too soon
		}
	shift@taxID;#remove leading parent ID for superkingdom so that values are aligned between arrays


#Add an extra bacteria reference
#At the kingdom slot
	if($names{$taxID[0]}[1] eq "Bacteria")
		{
		splice(@taxID, 1, 0, $taxID[0]);
		splice(@nodNam, 1, 0, "kingdom");
		}

#Hash for parsing of taxonomy
	my %taxPrint;
	for(my $i = 0; $i < @nodNam; $i++)
		{
		if($nodNam[$i] eq "superkingdom" ||
			$nodNam[$i] eq "kingdom" ||
			$nodNam[$i] eq "phylum" ||
			$nodNam[$i] eq "order" ||
			$nodNam[$i] eq "family" ||
			$nodNam[$i] eq "genus" ||
			$nodNam[$i] eq "species")
			{
			$taxPrint{$nodNam[$i]} = $names{$taxID[$i]}[1];	
			}
		}

#Printing
	print OUT $bls[0], "\t";
	my $lastDefinedLevel;
	if(defined($taxPrint{"superkingdom"}) == 1)
		{
		$lastDefinedLevel = $taxPrint{"superkingdom"};
		print OUT $lastDefinedLevel, "\t";
		}else{
		$lastDefinedLevel = "undef";
		print OUT $lastDefinedLevel, "\t";
		}
	if(defined($taxPrint{"kingdom"}) == 1)
		{
		$lastDefinedLevel = $taxPrint{"kingdom"};
		print OUT $lastDefinedLevel, "\t";
		}elsif($lastDefinedLevel eq "undef")
		{
		print OUT $lastDefinedLevel, "\t";
		}else{
		print OUT $lastDefinedLevel, "-undef\t";
		}
	if(defined($taxPrint{"phylum"}) == 1)
		{
		$lastDefinedLevel = $taxPrint{"phylum"};
		print OUT $lastDefinedLevel, "\t";
		}elsif($lastDefinedLevel eq "undef")
		{
		print OUT $lastDefinedLevel, "\t";
		}else{
		print OUT $lastDefinedLevel, "-undef\t";
		}
	if(defined($taxPrint{"order"}) == 1)
		{
		$lastDefinedLevel = $taxPrint{"order"};
		print OUT $lastDefinedLevel, "\t";
		}elsif($lastDefinedLevel eq "undef")
		{
		print OUT $lastDefinedLevel, "\t";
		}else{
		print OUT $lastDefinedLevel, "-undef\t";
		}
	if(defined($taxPrint{"family"}) == 1)
		{
		$lastDefinedLevel = $taxPrint{"family"};
		print OUT $lastDefinedLevel, "\t";
		}elsif($lastDefinedLevel eq "undef")
		{
		print OUT $lastDefinedLevel, "\t";
		}else{
		print OUT $lastDefinedLevel, "-undef\t";
		}
	if(defined($taxPrint{"genus"}) == 1)
		{
		$lastDefinedLevel = $taxPrint{"genus"};
		print OUT $lastDefinedLevel, "\t";
		}elsif($lastDefinedLevel eq "undef")
		{
		print OUT $lastDefinedLevel, "\t";
		}else{
		print OUT $lastDefinedLevel, "-undef\t";
		}
	if(defined($taxPrint{"species"}) == 1)
		{
		$lastDefinedLevel = $taxPrint{"species"};
		print OUT $lastDefinedLevel, "\t";
		}elsif($lastDefinedLevel eq "undef")
		{
		print OUT $lastDefinedLevel, "\t";
		}else{
		print OUT $lastDefinedLevel, "-undef\t";
		}
	
	print OUT "\n";
	}
