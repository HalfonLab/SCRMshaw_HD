#!/usr/bin/perl

## This script uses a homology file to translate the genes in a gene file of a reference
## species to a species of interest.

use strict;
use File::Basename;
use Getopt::Long;

my $help		= undef;

GetOptions("h|help=s"  => \$help);

my $hfile 	= $ARGV[0];		## Homology file : group-id spc gene.
my $sfile 	= $ARGV[1];		## Set file : set of dmel genes seperated by line spaces.
my $odir 		= $ARGV[2];		## Output directory.
my $argc		= $#ARGV + 1;	## Number of arguments.


## Check exit conditions
if (defined($help)) {
	Usage();
	exit(-1);
} elsif ($argc != 3) {
	print "Incorrect # of args given: $argc\n";
	Usage();
	exit(-1);
}

my %dmel2id = ();					## Hash : dmel symbol -> (group-id, group-id, group-id,..)
my %id2gene = ();					## Hash : group-id		-> (gene, gene, gene, gene, ...)

my ( $name, $path, $suffix ) = fileparse($sfile);
$name =~ s/\.genes$//;
my $oname = "$odir/$name";

## Make output directory if it doesn't exist
if ( !-d $odir ) { system("mkdir -p $odir"); }


open(HFILE,"<$hfile") or die ("Could not open $hfile for reading\n");
while (my $line = <HFILE>) {
    chomp($line);
		my ($id, $spc, $gene) = split(/\s+/,$line);
		if (uc($spc) =~ /DMEL/) {
			push(@{$dmel2id{$gene}}, $id);
		} else {
			push(@{$id2gene{$id}}, $gene);
		}
}
close(HFILE);

open(SFILE, "<$sfile") or die ("Could not open $sfile for reading\n");
open(OFILE, ">$oname.hom") or die("Could not open $oname.hom for reading\n");
while(my $line = <SFILE>) {
	chomp($line);
  my $dmel = $line;
	if (!defined($dmel2id{$dmel})) { next; }
	my @ids = @{$dmel2id{$dmel}};	
	foreach my $id (@ids) {
		if (!defined($id2gene{$id})) { next; }
		my @genes = @{$id2gene{$id}};
		foreach my $gene(@genes) {
			print OFILE "$dmel\t$gene\n";
		}
	}
}
close(OFILE);
close(SFILE);
#warn("sort -uo $oname.hom $oname.hom");
system("sort -uo $oname.hom $oname.hom");
#warn("cut -f 2 $oname.hom > $oname.genes; sort -uo $oname.genes $oname.genes");
system("cut -f 2 $oname.hom > $oname.genes; sort -uo $oname.genes $oname.genes");

sub Usage {
	my ( $name, $path, $suffix ) = fileparse($0);
	print "Usage: $name [OPTIONS] HFILE SFILE ODIR\n";
	print " HFILE    Ortholog group file.\n"; 
	print " SFILE    Gene set file.\n";
	print " ODIR     Output directory.\n";
	print " Options\n";
	print " -h       Print help and exit\n";
}
