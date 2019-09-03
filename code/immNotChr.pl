#!/usr/bin/perl

=head1 Description

msIMM socring scheme

=head1 Usage

    perl imm.pl [OPTIONS] PFILE NFILE ODIR WDIR

    PFILE    Positive CRM file.
    NFILE    Negative CRM file.
    ODIR     Output directory.
    WDIR     Window directory.
    
    --help             Print help and exit
    --kmer <int>       k-mer size to use. Default 6
    --win <int>        Window length. Default 500

=cut

use strict;
use FindBin qw($Bin $Script);
use lib '$Bin/../lib/';
use Getopt::Long;
use File::Basename;

my $k = 6;
my $wlen = 500;
my $help = undef;

GetOptions(
	"h|help=s"  => \$help,
	"k|kmer=i"  => \$k,     ## kmer length
	"w|win=i"   => \$wlen   ## window length
);

## input:
my $pfile = $ARGV[0];			## Positive CRM file. e.g. crms.fasta
my $nfile = $ARGV[1];			## Negative CRM file. e.g. neg.fasta
my $odir  = $ARGV[2];			## Output directory. e.g. scores/imm/mapping0.ap
my $wdir  = $ARGV[3];			## Window directory. e.g. fasta/windows

## Check exit conditions
die `pod2text $0` if (@ARGV != 4 || $help);

## Make output directory if it doesn't exist
if ( !-d $odir ) { system("mkdir -p $odir"); }

##======= train crm model =========##
system("$Bin/../bin/imm_build -r -k 6 < $pfile > $odir/crms.model");
##======= train neg model ==========##
system("$Bin/../bin/imm_build -r -k 6 < $nfile > $odir/neg.model");
##======== predict crm on each window =========##
while(my $wfile = <$wdir/*fasta>){
    my $name = (split /\//,$wfile)[-1];     ## e.g. 2L.fasta
	system("$Bin/../bin/imm_score -f -n $odir/neg.model $odir/crms.model $wfile > $odir/$name");
}
