#!/usr/bin/perl

=head1 Description

msIMM socring scheme

=head1 Usage

    perl imm.pl [OPTIONS] PFILE NFILE ODIR SDIR WDIR

    PFILE    Positive CRM file.
    NFILE    Negative CRM file.
    ODIR     Output directory.
    SDIR     Sequence FASTA directory.
    WDIR     Window directory.
    
    --help             Print help and exit
    --kmer <int>       k-mer size to use. Default 6
    --win <int>        Window length. Default 500

=cut

use strict;
use FindBin qw($Bin $Script);
use lib '$Bin/../lib/';
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;

my $k = 6;
my $wlen = 500;
my $help = undef;

GetOptions(
	"h|help=s"  => \$help,
	"k|kmer=i"  => \$k,
	"w|win=i"   => \$wlen
);

## input:
my $pfile = $ARGV[0];			## Positive CRM file.
my $nfile = $ARGV[1];			## Negative CRM file.
my $odir 	= $ARGV[2];			## Output directory.
my $sdir  = $ARGV[3];			## Sequence directory.
my $wdir  = $ARGV[4];			## Window directory.

## Check exit conditions
die `pod2text $0` if (@ARGV != 5 || $help);

## Make output directory if it doesn't exist
if ( !-d $odir ) { system("mkdir -p $odir"); }

system("$Bin/../bin/imm_build -r -k 6 < $pfile > $odir/crms.model");
system("$Bin/../bin/imm_build -r -k 6 < $nfile > $odir/neg.model");

while(my $sfile = <$sdir/*fasta>){
    chomp($sfile);
	my ($name, $path, $suffix) = fileparse($sfile);
#	print STDERR "$name\n";
    my $fasta = Bio::SeqIO->new(-file=>$sfile, 				## Open sequence file to read sequence
							-format=>'Fasta');
    my $fh = $fasta->next_seq();											## Get handle to next sequence. Only one in chromosome files.
	my $num = length($fh->seq());		#
	if ($num >= $wlen) {
	    system("$Bin/../bin/imm_score -f -n $odir/neg.model $odir/crms.model $wdir/$name > $odir/$name");
	} else {
		open(OFILE, ">$odir/$name") or die("Couldn't open $odir/$name for writing\n");	
		print OFILE "0 0\n";
		close OFILE;
	}
}
