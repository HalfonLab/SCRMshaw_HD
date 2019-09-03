#!/usr/bin/perl 

=head1 Description

    This script is for randomly extracting background non-coding sequences
    that have the same GC content as training CRM sequences. And these sequences 
    can be used as background training set for SCRM pipeline.

=head1 Usage

    perl randomWithSameGC.pl --crm CRM.fa --output OUTPUT.NAME --size SIZE --genomedir GENOMEDIR --gene GENE

    --crm <str>         CRM fasta file
    --output <str>      name of output background fasta file
    --size <int>        each output sequence is this many times the corresponding input sequences length
    --genomedir <str>   a directory that has all of the chromosomes of a genome
    --gene <str>        a gene file with format (delimited by "\t"): chrID, geneID, start, stop
    --prefix <str>      prefix of the chromosome fasta files(e.g. if filename is "chr1.fa" and chrID is "1", then the prefix should be "chr"), default=""
    --suffix <str>      suffix of the chromosome fasta files(e.g. if filename is "chr1.fa" and chrID is "1", then the suffix should be ".fa"), default=""
    --help              display help information to the screen

=cut

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::DB::Fasta;

#my $hg17dir = "/mounts/veda/disks/1/shared-data/GENOME/H.sapiens/hg17/";
#my $hg17suffix = ".fa";

#my $genefile = "/mounts/veda/disks/1/shared-data/GENOME/H.sapiens/data_hg17.chr";
my $proximity = 1000; ## extracted sequences at least this much away from genes

my ($infile, $outfile, $scalelength, $hg17dir, $genefile, $Help);
my ($Prefix, $Suffix) = ("", "");

GetOptions(
    'crm:s'=>\$infile,
    'output:s'=>\$outfile,
    'size:i'=>\$scalelength,
    'genomedir:s'=>\$hg17dir,
    'gene:s'=>\$genefile,
    'prefix:s'=>\$Prefix,
    'suffix:s'=>\$Suffix,
    'help'=>\$Help
);

#my $infile = $ARGV[0];   ## CRMs fasta file
#my $outfile = $ARGV[1];  ## output fasta file
#my $scalelength = $ARGV[2]; ## each output sequence is this many times the corresponding input sequences length
#my $hg17dir = $ARGV[3];  ## directory that has all of the chromosomes of genome sequence
#my $genefile = $ARGV[4]; ## gene file with format (delimited by "\t"): chrID, geneID, start, stop

die `pod2text $0` if (!$infile || !$outfile || !$scalelength || !$hg17dir || !$genefile || $Help); 

## check the exists of the modules
my $require_package = eval
{
    require Bio::SeqIO;
    require Bio::DB::Fasta;
    Bio::SeqIO->import();
    Bio::DB::Fasta->import();
    1;
};

if ($require_package)
{
    warn "Bio::SeqIO loaded and imported successfully\n";
    warn "Bio::DB::Fasta loaded and imported successfully\n";
}
else
{
    die("Bio::SeqIO or Bio::DB::Fasta cannot be imported. Please check to see if the modules have already been installed.\n");
}

open(OUT,">$outfile");

## read and store the gene annotations
open(G,"<$genefile");
my %genes = ();
my %genestop = ();
while (<G>) {
    chomp;
    my ($chr, $gene, $start, $stop) = split(/\s+/);
    if ($start > $stop) { ($start, $stop) = ($stop, $start); }
    push @{$genes{$chr}}, $start;
    $genestop{$chr}{$start} = $stop;
    ## warn "Seen gene $chr $start $stop\n";
}
close(G);
foreach my $chr (keys %genes) {
    my @sortedgenes = sort {$a <=> $b} @{$genes{$chr}};
    @{$genes{$chr}} = @sortedgenes;
}

## record approximate lengths of each chromosome
my %lengths;
my %index2chr;
my $count = 0;
for my $chr (keys %genes) {
    my @genelist = @{$genes{$chr}};
    my $maxstartpos = 0;
    foreach my $start (@genelist) {
	if ($start > $maxstartpos) { $maxstartpos = $start; }
    }
    $lengths{$chr} = $maxstartpos;
    $index2chr{$count} = $chr;
    $count++;
}
my @chrs = keys %lengths;
my $numchrs = $#chrs+1;

## read the promoter database
my %numwindows;
my $count = 0;
my $fo = Bio::SeqIO->new(-file=>$infile,
			 -format=>'Fasta');
while (my $f = $fo->next_seq()) {
    my $origid = $f->id();
    my $origseq = $f->seq();
    my $origlength = length($origseq);
    my $gcpct = GCPCT($origseq);

    ## now get a random sequence of this length time scale factor
    my $numAttempts = 0;
    my $found = 0;
    while ($found==0 && $numAttempts < 10000) {
	## choose a random non-coding window of length $origlength*$scalelength
	my $r = int(rand($numchrs));
	my $chr = $index2chr{$r};
	$r = int(rand($lengths{$chr}));
	my ($start,$stop) = ($r,$r+$origlength*$scalelength-1);
	my $name = "$chr\_$start\_$stop";
	## warn "Trying random window $chr $start $stop\n";

	## check if it overlapping with or close to a gene
	## if so next
	if (OverlapsGene($chr, $start, $stop)==1) {
	    ## warn "$chr $start $stop Overlaps Gene\n";
        $numAttempts++; # added by Wei
	    next;
	}

	$numAttempts++;

	## get the sequence 
	my $hdb = Bio::DB::Fasta->new("$hg17dir/$Prefix$chr$Suffix");
	my $hseq = $hdb->seq("$Prefix$chr",$start => $stop);

	my @Ns = ($hseq =~ m/N/g);
	if ($#Ns>=0) { next; }

	## compute gcpct
	my $thisgcpct = GCPCT($hseq);
	if (abs(int($gcpct) - int($thisgcpct)) > 1) {
	    next;
	}

	## and print it
	my $title = "$origid\_$name";
	print OUT ">$title\n$hseq\n";

	$found = 1;
    }
    if ($found == 0) {
	warn "Couldnt find length-matching random sequence for $origid\n";
    }
}    
close(OUT);


sub OverlapsGene {
    my ($chr,$start,$stop) = @_;
    my $length = $stop - $start + 1;

    ## warn "checking overlap for $chr $start $stop\n";
    ## find the nearby entry in the array    
    my @genelist = @{$genes{$chr}};
    my $nindex = bsearch($start, \@genelist);
    ## nindex^th entry (startpoint) in genelist is close to "$start"
    ## check all genes within 10 genes to left and right of nindex
    for (my $ind = max(0,$nindex-10); $ind <= ($nindex+10); $ind++) {
	if ($#genelist < $ind) { next; }
	my ($istart, $istop) = ($genes{$chr}[$ind], $genestop{$chr}{$genes{$chr}[$ind]});
	if ((OverlapsInterval($start,$stop,$istart,$istop) == 1) or (NearbyIntervals($start,$stop,$istart,$istop) == 1)) {
	    return 1;
	}
    }
    return 0;
}

sub max {
    my ($a, $b) = @_;
    if ($a < $b) { return $b; } 
    return $a;
}

sub OverlapsInterval {
    my ($b1,$e1, $b2,$e2) = @_;
    #warn "Checking overlap of intervals $b1 $e1 with $b2 $e2\n";
    if ($b2 >= $b1 and $b2 <= $e1) { return 1; }
    if ($b1 >= $b2 and $b1 <= $e2) { return 1; }
    return 0;
}

sub NearbyIntervals {
    my ($b1,$e1, $b2,$e2) = @_;
    #warn "Checking closeness of intervals $b1 $e1 with $b2 $e2\n";
    my $max = -1;
    if (max(abs($e1-$b2),abs($e2-$b1)) < $proximity) { return 1; }
    return 0;
}

# Binary search, array passed by reference

# search array of integers a for given integer x
# return index where found or -1 if not found
sub bsearch {
    my ($x, $a) = @_;            # search for x in array a
    my ($l, $u) = (0, @$a - 1);  # lower, upper end of search interval
    my $i;                       # index of probe
    while ($l <= $u) {
	if ($u <= ($l+2)) { return $l; }
	$i = int(($l + $u)/2);
	if ($a->[$i] < $x) {
	    $l = $i+1;
	}
	elsif ($a->[$i] > $x) {
	    $u = $i-1;
	} 
	else {
	    return $i; # found
	}
    }    
    return -1;         # not found
}

sub GCPCT {
    my $seq = shift;
    my @gc = ($seq =~ m/[GCgc]/g);
    my $numGC = $#gc+1;
    my @at = ($seq =~ m/[ATat]/g);
    my $numAT = $#at+1;
    my $sum = $numGC+$numAT;
    my $gcpct = -1;
    if ($sum > 0) {
	$gcpct = int($numGC*100/($numGC+$numAT));	
    }
    return $gcpct;
}
