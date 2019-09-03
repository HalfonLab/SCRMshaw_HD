#!/usr/bin/perl

=head1 Description

    This file concatenates k-mers from different files such that sequences have distinct ids.

=head1 Usage

    perl $name [OPTIONS] KDIR K
    
    KDIR        Directory containing k-mers.
    K           k-mer size.
    
    --help      Print help and exit.
    --rc        Use reverse complement.

=cut

use strict;
use File::Basename;
use Getopt::Long;

my @count	=	();
my $rc		=	undef;;
my $help	= undef;

## Get Options
GetOptions( 
    "h|help" => \$help, 
    "r|rc" => \$rc
);

my $kdir 	= $ARGV[0];
my $k			= $ARGV[1];
my $suf		=	defined($rc) ? "fasta\.rc\.$k" : "fasta\.$k";
my $max		= 4**$k;

## Check exit conditions
die `pod2text $0` if (@ARGV != 2 || $help);

my $count = 0;
while(my $ifile = <$kdir/*>) {
	chomp($ifile);
	if ($ifile !~ /$suf$/) { next; }
	#print STDERR "$ifile\n";
	open(IFILE, "<$ifile") or die("Couldn't open $ifile for reading\n");
	while (my $line = <IFILE>) {
		chomp($line);
		my @vals = split("\t",$line);
		my $len = $#vals+1 - 3;
		if ($len != $max) { die("Incorrect length for kmer of lengths $k. $len != $max\n"); }
		print$count."\t".join("\t", @vals[1..$#vals])."\n";
		$count++;
	}
	close(IFILE);
}



