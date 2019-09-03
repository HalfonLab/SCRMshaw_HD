#!/usr/bin/perl

=head1 Description

    This file inputs a k-mer file and sums the individual
    k-mer lines to produce a single k-mer count of all sequences represented.

=head1 Usage

    Usage: kmers_sum.pl [OPTIONS] KFILE K
    
    KFILE File containing k-mers.
    K     k-mer size.
    
    --help    Print help and exit.

=cut

use strict;
use File::Basename;
use Getopt::Long;

my @count	=	();
my $help	= undef;

## Get Options
GetOptions( "h|help" => \$help);

my $kfile	= $ARGV[0];
my $k		= $ARGV[1];
my $max		= 4**$k;

## Check exit conditions
die `pod2text $0` if (@ARGV != 2 || $help);

open(IFILE, "<$kfile") or die("Couldn't open $kfile for reading\n");
while (my $line = <IFILE>) {
	chomp($line);
	my @vals = split("\t",$line);
	my $len = $#vals+1 - 3;
	if ($len != $max) { die("Incorrect length for kmer of lengths $k. $len != $max\n"); }
	for (my $i = 0; $i < $max + 1; $i++) {
		$count[$i] += $vals[$i + 2];
	}
}
close(IFILE);

print  "0\t-";
for (my $i = 0; $i < $max + 1; $i++) {
	print  "\t".$count[$i];
}
print  "\n";



