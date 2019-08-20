#!/usr/bin/perl
=head1 Description

This script run HexMCD by accepting positive and negative k-mer sum files (k = 6 and rc)
It needs sequence and window directories as well as an output directory.

=head1 Usage

    perl hexmcd.pl [OPTIONS] PFILE NFILE ODIR SDIR KDIR
    
    PFILE    File containing crm k-mers
    NFILE    File containing neg k-mers
    ODIR     Output directory
    SDIR     Sequence directory
    KDIR     K-mers directory
    
    --help              Print help and exit
    --kmer <int>        k-mer length, default 6
    --rc                Use reverse complements
    --win <int>         Window length. Default 500
    --shift <int>       Shift length. Default 250

=cut
use strict;
use FindBin qw($Bin $Script);
use lib '$Bin/../lib';
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;

my $help	= undef;
my $rc		= undef;
my $wlen    = 500;
my $slen    = 250;
my $k		= 6;

GetOptions(
	"h|help=s"  => \$help,
	"k|kmer=i"  => \$k,
	"r|rc" 		=> \$rc,
	"s|shift=i" => \$slen,
	"w|win=i"   => \$wlen
);

my $pfile = $ARGV[0];			## Positive CRM file.
my $nfile = $ARGV[1];			## Negative CRM file.
my $odir  = $ARGV[2];			## Output directory.
my $sdir  = $ARGV[3];			## Sequence directory.
my $kdir  = $ARGV[4];			## Window directory.
my $suf	  =	defined($rc) ? "rc.$k" : "$k";


## Check exit conditions
die `pod2text $0` if (@ARGV != 5 || $help);

## Make output directory if it doesn't exist
if ( !-d $odir ) { system("mkdir -p $odir"); }

while(my $sfile = <$sdir/*fasta>){
	chomp($sfile);
	my ($name, $path, $suffix) = fileparse($sfile);
    my $fasta = Bio::SeqIO->new(-file=>$sfile, 	## Open sequence file to read sequence
								-format=>'Fasta');
    my $fh = $fasta->next_seq();	## Get handle to next sequence. Only one in chromosome files.
	my $num = length($fh->seq());	
	if ($num >= $wlen) {
	   my $count =	int(($num-$wlen)/$slen) + 1;
	   system("$Bin/../bin/compare_hex $pfile $nfile $kdir/$name.$suf $count $slen > $odir/$name");
	} else {
		open(OFILE, ">$odir/$name") or die("Couldn't open $odir/$name for writing\n");	
		print OFILE "0 0\n";
		close OFILE;
	}
}

