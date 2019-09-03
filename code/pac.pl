#!/usr/bin/perl

=head1 Description

This script runs PAC by accepting 
two types of k-mer sum files (typically 6 and 1) and a ymf file.
It needs the sequence, window, and k-mer directories of the species.
and an output directory to place scores.

Original code recommends running PACrc with no reverse complement for dmel sets.
Original code recommends running PACRC with reverse complement for dros sets.

=head1 Usage

    perl $name [OPTIONS] PFILE NFILE YFILE ODIR SDIR KDIR
    
    PFILE    File containing crm k-mers
    NFILE    File containing crm n-mers
    YFILE    File containing YMF k-mers
    ODIR     Output directory
    SDIR     Sequence directory
    KDIR     K-mers directory
    
    --help              Print help and exit
    --kmer <int>        k-mer length, default 6
    --nmer <int>        n-mer length, default 1
    --rc                Use reverse complements
    --shift <int>       Shift length. Default 250
    --win <int>         Window length. Default 500

=cut

use strict;
use FindBin qw($Bin $Script);
use lib '$Bin/../lib/';
use Bio::SeqIO;
use File::Basename;
use Getopt::Long;

my $help	= undef;
my $rc		= undef;
my $wlen    = 500;
my $slen    = 250;
my $k	    = 6;
my $n		= 1;

GetOptions(
	"h|help=s"  => \$help,
	"k|kmer=i"  => \$k,
	"n|nmer=i" => \$n,
	"r|rc" 			=> \$rc,
	"s|shift=i" => \$slen,
	"w|win=i"   => \$wlen
);

## input:
my $kfile   = $ARGV[0];			## CRM K mer file.
my $nfile   = $ARGV[1];			## CRM l mer file.
my $yfile	= $ARGV[2];			## YMF k-mer file.
my $odir 	= $ARGV[3];			## Output directory.
my $sdir    = $ARGV[4];			## Sequence directory.
my $kdir    = $ARGV[5];			## Window k-mer directory.
my $ksuf    =	defined($rc) ? "rc.$k" : "$k";
my $nsuf	=	defined($rc) ? "rc.$n" : "$n";

## Check exit conditions
die `pod2text $0` if (@ARGV != 6 || $help);

## Make output directory if it doesn't exist
if ( !-d $odir ) { system("mkdir -p $odir"); }

while(my $sfile = <$sdir/*fasta>){
	chomp($sfile);
	my ($name, $path, $suffix) = fileparse($sfile);
    my $fasta = Bio::SeqIO->new(-file=>$sfile, 		## Open sequence file to read sequence
								-format=>'Fasta');
    my $fh = $fasta->next_seq();		## Get handle to next sequence. Only one in chromosome files.
	my $num = length($fh->seq());	
	if ($num >= $wlen) {
	    my $count =	int(($num-$wlen)/$slen) +1;
		print "$name\n";
	    system("$Bin/../bin/pac $kfile $kdir/$name.$ksuf $nfile $kdir/$name.$nsuf $count $k $n $slen $yfile > $odir/$name");
	} else {
		open(OFILE, ">$odir/$name") or die("Couldn't open $odir/$name for writing\n");	
		print OFILE "0 0\n";
		close OFILE;
	}
}


