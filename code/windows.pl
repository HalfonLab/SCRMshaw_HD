#!/usr/bin/perl

#modified 08-09-2019

=head1 Description

    This script iterates over each FASTA file in a directory and segments the input 
    into overlapping windows, each of which is a FASTA sequence in the output file.
    Additionally, it computers k-mer files for each of the window sequences.

=head1 Usage

    perl windows.pl [OPTIONS] CDIR WDIR KDIR

    CDIR     Chromosome directory with masked FASTA files
    WDIR     Window directory for output window FASTA files
    KDIR     k-mer directory for output kmer window files
    
    --help              Print help and exit
    --kmer  <str>       CSV list of k-mer lengths to get k-mer files
    --kmerc <str>       CSV list of k-mer length to get k-mer with rc files
    --win   <int>       Window length. Default 500
    --shift <int>       Shift length. Default 250
    --offset <int>	Left boundary i.e. starting points to take windows on each chromosome

=cut

use strict;
use FindBin qw($Bin $Script);
use lib '$Bin/../lib/';
use File::Basename;
use Bio::SeqIO;
#use Bio::DB::Fasta;
use Getopt::Long;

my $help	=	undef;	## 
my $kstr	= "1,6";		## CSV list of kmers to get windows for.
my $nstr 	= "1,6";		## CSV list of kmers with reverse complement to get distributions.
my $wlen	= 500;		## Window length.
my $slen	= 250;		## Shift length.
my $offset	= 0;		##start position on chromosome
GetOptions (
	"h|help"	=> \$help,
	"k|kmer=s"	=> \$kstr,
	"n|kmerc=s"	=> \$nstr,
	"s|shift=i" => \$slen, 
	"w|win=i"	=> \$wlen,
	"o|offset=i" => \$offset,
);

my $cdir	= $ARGV[0];		##	Chrom  directory.	
my $wdir	= $ARGV[1];		## 	Window directory.
my $kdir	= $ARGV[2];		## 	kmer directory.
#my ($name, $path, $suffix) = fileparse($0);

## Check exit conditions
die `pod2text $0` if (@ARGV != 3 || $help); 

## Make output directory if necessary
if (! -d $wdir) {	system("mkdir -p $wdir");}
if (! -d $kdir) {	system("mkdir -p $kdir");}

while(<$cdir/*.fasta>){
	chomp(my $line=$_);
#	my ($name, $path, $suffix) = fileparse($_);
#	my $fasta = Bio::DB::Fasta->newFh($line, -format=>'Fasta');
    my $name = basename($line);
#	my $chr = $name;
#	$chr =~ s/\.fasta//g;
	my $fasta = Bio::SeqIO->new(-file=>"$cdir/$name", -format=>'Fasta');
#    my $fasta = Bio::DB::Fasta->newFh("$cdir/$name", -format=>'Fasta');
	while (my $fh = $fasta->next_seq()) {
		my $id 		= $fh->id();				## 	Get ID
        my $seq 	= uc($fh->seq());				## 	Uppercase sequence.
        # $seq =~ tr/atgcn/ATGCN/;
        chomp $seq;
        #print "$seq\n";
		my $len 	= length($seq) - $offset;			## 	Length of sequence.
		my $des		= $fh->desc();			## 	Get description.
		my $wfile	=	"$wdir/$name";		## 	Windowed output file.
		my $lb		= $offset;				## 	Left boundary of current window: need to subtract 1?
		my $rb		= $len - $wlen + 1;	## 	Right side barrier of all windows.
		my $num		= 1;								## 	Number of windows.
		## Slide a window
		open(WFILE, ">$wfile") or die("Couldn't open $wfile for writing\n");
		while($lb < $rb) {
            my $subseq	=	substr($seq, $lb, $wlen);
            #  my $subseq = $fh->subseq($lb+1=>$lb+$slen+1);
			print WFILE ">$lb\n$subseq\n";
			$lb += $slen;
			$num++;
		}
		close(WFILE);
		foreach my $k (split(/\,/,$kstr)) {
			system("perl $Bin/kmers.pl $wfile $k > $kdir/$name\.$k");
		}

		foreach my $n (split(/\,/,$nstr)) {
			system("perl $Bin/kmers.pl --rc $wfile $n > $kdir/$name\.rc\.$n");
		}
	}
}

