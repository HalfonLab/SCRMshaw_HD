#!/usr/bin/perl
#modified 01-13-2020 by Hasiba Asma(HA)
#following changes were made to fix PAC issue
#instead of moving files back and forth, i have changed the directory to the respective training set's directory
#then ran preproc and statsvar from Bin and finally at the end changed the directory back to where it was
=head1 Description

    This script runs ymf, which is needed by PAC.

=head1 Usage

    Usage: ymf.pl [OPTIONS] PFILE NFILE ODIR
    
    PFILE FASTA file containing positive crm sequences.
    NFILE FASTA file containing negative sequences.
    ODIR  Output directory.
    
    --help          Print help and exit.
    --kmer <int>    Set k-mer size, default=6
    --size <int>    Set the number of ymf words to return. Default 100.

=cut

use strict;
use FindBin qw($Bin $Script);
#use lib '/home/sinhas/sinhas/lib/';
use File::Basename;
use Getopt::Long;
use Bio::SeqIO;
use File::Spec;
#my $root 	= "/shared-mounts/sinhas-storage1/kruppel/pipeline";
my $k = 6;					## k-mer length
my $size = 100;				## Top words
my $tfile = "table.3";	## Intermediary files
my $wfile = "powers.3";
my $bfile = "powersGeneralized.3.bin";
my $rfile = "results";
my $help = undef;

GetOptions (
	"h|help"	=> \$help,
	"k|kmer=i"  => \$k,
	"s|size=i"	=> \$size
);

my $pfile = $ARGV[0];		## Positive FASTA file
my $nfile = $ARGV[1];		## Negative FASTA file
my $odir = $ARGV[2];		## Output directory.

#absolute path
$pfile = File::Spec->rel2abs($pfile);
$nfile = File::Spec->rel2abs($nfile);
$odir = File::Spec->rel2abs($odir);
## Check exit conditions
die `pod2text $0` if (@ARGV != 3 || $help);
#HA--orig dir
my $origDir=`pwd`;

## Make output directory if necessary
if (! -d $odir) {	system("mkdir -p $odir");}

#HA--change dir
chomp($odir);
chdir($odir);
my $p=`pwd`;
chomp($p);
#warn "inside tset dir that is my pwd: $p \n";

my ($nneg) = qx/grep \">\" $nfile | wc -l/;		## Number of negative sequences.
chomp($nneg);

## Preprocess negative file for YMF
#warn "$Bin/../bin/preproc $nneg $nfile\n";
system("$Bin/../bin/preproc $nneg $nfile");

#HA--comment the following
## Move intermediary files
#system("mv $tfile $odir/_$tfile; mv $wfile $odir/_$wfile; mv $bfile $odir/_$bfile");
#HA--to preserve naming
system("cp $tfile _$tfile; cp $wfile _$wfile; cp $bfile _$bfile");

## use the above tables in YMF runs on the training crms
#warn("$Bin/../bin/statsvar $Bin/../bin/stats.config 30000 6 $odir/ -sort $pfile");
#system("$Bin/../bin/statsvar $Bin/../src/ymf/config/stats.config 30000 6 $odir/ -sort $pfile");
system("$Bin/../bin/statsvar $Bin/../src/ymf/config/stats.config 30000 6 $p/ -sort $pfile");
#HA--comment the following
## Move intermediary files
#system("mv $rfile $odir/$rfile");

#HA--fix path
## read in the words file
#open(IFILE, "<$odir/$rfile") or die("Couldn't open $odir/results for reading\n");
open(IFILE, "<$rfile") or die("Couldn't open $odir/results for reading\n");
my $header = <IFILE>; #toss out the header
my $words = "";
for (my $i = $size; $i; $i--){
    my $line = <IFILE>;
    chomp($line);
    my @wordline = split(" ",$line);
    $words = $words.$wordline[0].'N';
}
close(IFILE);

#HA--fix path
#open( OFILE,">$odir/ymf.$size.fasta") or die("Couldn't open $odir/ymf.$size.fasta writing\n");
open( OFILE,">ymf.$size.fasta") or die("Couldn't open $odir/ymf.$size.fasta writing\n");
print OFILE ">ymfwords\tsize=$size,k=$k\n$words\n";
close(OFILE);

#HA--- undo remove
## Remove YMF intermediary output
#system("rm $odir/_$tfile $odir/_$wfile $odir/_$bfile $odir/$rfile");

#HA--change back to orig directory
chomp($origDir);
chdir($origDir);
my $cd=`pwd`;
#warn "back to orig directory so my pwd is: $cd \n";
