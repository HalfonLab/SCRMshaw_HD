#!/usr/bin/perl

=head1 Description

This file accepts a FASTA file of a genome, 
and outputs FASTA files for each chr separately.

=head1 Usage

    perl genome2chr.pl [Options] FASTA OUTDIR

    FASTA       Fasta file.
    OUTDIR      output directory

    --exon <str>    Takes exon file and will mask exons on the chr, file format(tab delimited): chr start stop strand symbol

=cut

use strict;
use FindBin qw($Bin $Script);
use lib "$Bin/../lib";
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;

my $Exon;

GetOptions(
    "exon:s"=>\$Exon        ## option to take exon file, but the code can run without exon
);


my $ffile   = $ARGV[0];		## Input FASTA file.
my $odir	= $ARGV[1];		## Output directory for each chr.
my %seen	= ();			## Hash of chrs seen.
my %exons	= ();			## Hash of chrs to exon start and stop.


## Check exit conditions
die `pod2text $0` if (@ARGV != 2);

##========= Read gene file to get which crhs have genes ==========##
my %keep    = ();       ## Hash of chrs to keep.
my %exons   = ();       ## Hash of chrs to exon start and stop.
if ($Exon){             ## if has exon as input
    open(IFILE,$Exon) or die("Couldn't open exon file: $Exon for reading\n");
    while (my $line = <IFILE>) {
        chomp($line);
		my ($chr, $start, $stop, $strand, $gene) = split(/\s+/,$line); ## 2R 50001 50235 + eve:1
		$keep{$chr} = 1;        ## keep{2R} = 1
        ($start,$stop)=($stop,$start) if ($start > $stop);  ## swap start and stop if coordinates are reversed
        push @{$exons{$chr}}, "$start:$stop";       ## @{$exons{2R}}, "50001:50235"
    }
}
close(IFILE);

##============ Mask exons in each chr ==============##
my $fo = Bio::SeqIO->new(-file=>$ffile,-format=>'Fasta');
while (my $f = $fo->next_seq()){
    my $title= $f->id();
    my $chr=(split /\s+/,$title)[0];
    my $seq = $f->seq();
    chomp $seq;
    $seq =~ s/\s+//g;
    $seq =~ tr/atgc/ATGC/;
    my $len = length($seq);
   # if ($Exon && !defined($keep{$chr}))  { next;}   ## next if there is exon file but not keep{chr}
    if (defined($seen{$chr}))   { next;}        ## next if $seen{$chr}
	$seen{$chr} = 1;            ## seen{chr}
	my $ofile_dir	=	"$odir/$chr.fasta";     ## output to file
    # my $ofile   =   "$chr.fasta";
    open(OFILE,">$ofile_dir") or die("Couldn't open $ofile_dir for writing\n");
    print OFILE "\>$chr\n";
    if ($Exon){     ## yes if has exon input
	    foreach my $exon (@{$exons{$chr}}) {    ## genes in chr
	        my ($start, $stop) = split(":",$exon);
            substr($seq,$start-1,$stop-$start) =~ tr/ATGC/N/;   ## mask exon
	    }
    }
    for (my $i = 0; $i < $len; $i += 60) {      ## output seq, 60bp per row
        print OFILE substr($seq, $i, 60)."\n";
    }
    close(OFILE);
}

