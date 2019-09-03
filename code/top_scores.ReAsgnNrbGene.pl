#!/usr/bin/perl

=head1 Description

    This script takes scores and sequences of windows
    from given method (e.g. msIMM) and crm train set,
    along with a gene file as input, and will 
    output top-scored crms with their nearest genes. 
    Once assigned to a crm, the same nearest gene can 
    be re-assigned to other crm. There are two options
    to determine the number of hits to output: 1) constraint
    number of windows uisng --topw; 2) constraint
    number of genes using --topg. Please note that both
    options cannot be used in the same time.
    This version of top_score can take Universal Gene Set.

=head1 Usage

    perl top_scores.pl [OPTIONS] ODIR WDIR GENE OUTFILE
    
    ODIR     Directory containing scores for FASTA files
    SDIR     Directory containing window files
    GENE      GENE file with gene locations in species
    
    --help              Print help and exit.
    --topw <int>        Maximum number of windows. e.g. 2000. Default=0
    --topg <int>        Maximum number of genes. e.g. 200. Default=0
    --win <int>         Window length. Default 500.
    --universe          Universal Gene Set

=cut

use strict;
use FindBin qw($Bin $Script);
use lib '$Bin/../lib/';
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;

my $help 	= undef;
my $wlen	= 500;      ## window length
my $hitw    = 0;        ## number of windows
my $hitg    = 0;        ## number of genes
my $Universe = undef;

GetOptions(
    "h|help" => \$help, 
    "w|win=i" => \$wlen, 
    "topw=i" => \$hitw,     ## number of top windows
    "topg=i" => \$hitg,      ## number of top genes
    "universe=s" => \$Universe
);

my $odir 	= $ARGV[0];  	    ## Directory of window scores
my $sdir 	= $ARGV[1]; 	    ## Directory of window sequences
my $gff  	= $ARGV[2];		    ## genes files.
my $outFile = $ARGV[3];         ## output hits file
#my $universe = "/home/n-z/weiyang4/scrm/GBE/data.fromMajid/GeneSet/Universe/dmel_universe_cg.genes";

my %g2pos	= ();				## g2pos{chr} array of geneIDs sorted by 5'+pos
my %pos2g	= ();				## pos2g{chr}{5'+pos} = (3'+ pos, strand, geneID)
my %score   = ();				## score{chr:w} = score.
my %uhit	= ();				## Universe genes has: 1/0 if reported already/not.
my $nhitw 	= 0;				## Total windows reported.
my $nhitg   = 0;                ## Total genes reported.
my %sequence;                   ## sequence{chr:w} = CTCCAAAATGT...

## Check exit conditions
die `pod2text $0` if (@ARGV != 4 || $help || ($hitw == 0 && $hitg == 0) || ($hitw > 0 && $hitg > 0));

if (defined $Universe){ ## return 1 if defined universe
    open (UNI,$Universe) or die ("Couldn't open $Universe\n");
    while (<UNI>){
        chomp(my $line = $_);
        $uhit{$line} = 1;       ## mark genes in the universe set
    }
    close UNI;
}


open OUT,">$outFile";
##============ Step 1: Parse gff file ============##
open(IFILE, "<$gff ") or die ("Couldn't open $gff for reading\n");
while (<IFILE>) {
	chomp($_);
	my ($chr, $u, $d, $str, $id) = split(/\s+/, $_);    ## e.g. 2L, 7529, 9484, +, FBgn0031208
    ($u, $d) = (min($u, $d), max($u, $d));
	@{$pos2g{$chr}{$u}} = ($d, $str, $id);      ## e.g. @{$pos2g{2L}{7529}} = (9484, +, FBgn0031208)
	push (@{$g2pos{$chr}}, $u);                 ## e.g. @{$g2pos{2L}}, 7529
    if (!defined $Universe){
        $uhit{$id} = 1;         ## no universe gene set, mark genes in the gff file
    }
}
close(IFILE);

##============ Step 2: Sort g2pos by upstream/5'+ pos on chr ============##
foreach my $chr (keys %g2pos) { 
	@{$g2pos{$chr}} = sort { $a <=> $b } @{ $g2pos{$chr} };
}

##============ Step 3: Read score for each chr ===========##
while (<$odir/*.fasta>) {
	chomp(my $file = $_);
	my $name = (split /\//,$file)[-1];						## e.g. 2L.fasta
	my $chr = (split /\.fa/,$name)[0];					## e.g. 2L
	open(IFILE, "<$file") or die("Couldn't open $file for reading\n");
	while (<IFILE>) {
		chomp;													
		my ($w, $score) = split(/\s+/);                 ## e.g. 5500, 2.22802
		$score{"$chr:$w"} = $score;			            ## e.g. score{2L:5500} = 2.22802
	}
	close(IFILE);
}

##============ Step 4: Read window seq for each chr ===========##
while (<$sdir/*.fasta>) {
    chomp(my $file = $_);
    my $name = (split /\//,$file)[-1];      ## e.g. 2L.fasta
    my $chr = (split /\.fa/,$name)[0];        ## e.g. 2L
    my $fo = Bio::SeqIO->new(-file=>$file,-format=>'Fasta');
    while (my $f = $fo->next_seq()) {    
        my $title = $f->id();
        my $seq = $f->seq();
        my $seq_name = $1 if($title =~ /^(\S+)/);   ## e.g. 750

	$seq_name = 0 if($title == 0);  ## fixes mysterious issue where above regex fails if $f->id() returns 0

        chomp $seq;
        $seq =~ s/\s+//g;
        $sequence{"$chr:$seq_name"} = $seq;     ## e.g. sequence{2L:750} = CTCCAAAATGT...

    }
}

##============ Step 5: Sort all windows by score in descending order ============##
##	Until nhit = hitw or no more windows, test windows for hits.
##	A window w should be reported as a hit if : 
##	 	1. w's nearest gene is in the Universe Set  AND
##		2. w's nearest gene hasn't been reported yet with any window.
##	To report window, print the following info in FASTA format:
##		w's position, score, closest gene inf, and sequence
foreach my $key(sort{ $score{$b} <=> $score{$a}} keys %score) {
    if ($hitw > 0){
    	# warn $hitw." > 0\n";
        if ($nhitw >= $hitw) {last;}        ## don't output more than 2000
        my ($chr, $w)   = split(/\:/, $key );       ## chr, w
        my $score   = $score{$key};                 ## window score, e.g. 2.22802
        my $csv = closestGenes($chr,$w);            ## get 2 closest genes to window.
        my @colsep  = split(/\:/, $csv);            ## split by :
        my $gene = $colsep[0];                      ## closest gene
        my $wseq = $sequence{$key};                 ## window seq

        $nhitw = $nhitw + 1;                        ## incr hits and print hit to STDOUT
        print OUT ">$key\t$score\t$csv\n";      ## ex:>chr:w score neighbor-info
        print OUT "$wseq\n";                    ## ex:seq[w ... w + wlen - 1]        
    }
   
    if ($hitg > 0){
    	# warn $hitg." > 0\n";
        if ($nhitg >= $hitg) {last;}
	    my ($chr, $w) = split(/\:/, $key );			## chr, w
	    # warn "chr=".$chr." w=".$w."\n";
        my $score = $score{$key};					## window score, e.g. 2.22802
	    my $csv	= closestGenes($chr,$w);		    ## get 2 closest genes to window.
	    # warn "csv=".$csv."\n";
	    if ($csv =~ /Unknown/) { next; }			## skip if gene unknown.
	    my @colsep = split(/\:/, $csv);				## split by :
	    my $gene = $colsep[0];                      ## closest gene
        my $wseq = $sequence{$key};                 ## window seq
        if ($uhit{$gene} == 1){             		## return 1 if this gene mark in universe/gff
	        $nhitg = $nhitg + 1;					## incr hits and print hit to STDOUT
            $uhit{$gene} = 2;
        }
	    print OUT ">$key\t$score\t$csv\n";      ## ex:>chr:w score neighbor-info
	    print OUT "$wseq\n";				    ## ex:seq[w ... w + wlen - 1]
    }
}
close OUT;

##=============== subroutine ===============##

## Input-Output: string of path to FASTA file - Array(id, desc, seq) from file
sub ReadFASTA {
	my $file = shift;
	my $fo   = Bio::SeqIO->new(-file =>"$file",-format =>'Fasta');
	my $f    = $fo->next_seq();
	return ($f->id(), $f->desc(), $f->seq());
}

## Input-Output: string - bool if has '-' symbol of minus strand.
sub minus {
	my $str = shift;
	return $str =~ /\-/;
}

## Input-Output: string - bool if has '+' symbol of plus strand.
sub plus {
	my $str = shift;
	return $str =~ /\+/;
}

## Input:	 	Array	(w, u, d)
##					w		Upstream window boundary on + strand
##					u		Upstream gene boundary on + strand
##					d		Downstream gene boundary on + strand
## Output:	True if $w in gene=[$u, $d]
sub overlap {
	my ($w, $u, $d) = @_;
	($w-$u)*($w-$d) <= 0; 	## Return true if $w in [$u,$d]
}

## Input:	 	Array	(w, u, d)
##					w		Upstream window boundary on + strand
##					u		Upstream gene boundary on + strand
##					d		Downstream gene boundary on + strand
## Output:	Minimum distance of $p to gene [$u, $d], or 0 if within.
sub distToGene {
	my ($w, $u, $d) = @_;
	return overlap($w,$u,$d) ? 0 		 : ## if overlap, dist is 0
				min(abs($u-$w), abs($d-$w)); ## else, min dist of w2d and w2u
}

## Input-Output:  	scalars a, b - max of a,b
sub max {
	my ( $a, $b ) = @_;
	return ($a < $b) ? $b : $a;
}

## Input-Output:	scalars a,b - minimum of a,b
sub min {
	my ( $a, $b ) = @_;
	return ($a < $b) ? $a : $b;
}

## Input-Output:		scalar a - absolute value of a
sub abs {
	my $a = shift;
	return ($a < 0) ? -$a : $a;
}

## Input:	Scalars chr and w which are the chromosome and pos of window w.
## Output:	String with info on closest gene or gene pair (order by proximity)
## Example: 500:1000:+:AAEL0021,2000:2500:+:AAEL0883
## Notes:
##			The function may return only one gene.
##			The closest gene is start=500,  stop=100,  sense = +, id = AAEL0021
##			2nd closest gene is start=2500, stop=2500, sense = +, id = AAEL0883
##			":" breaks fields within a gene.
##			"," breaks between genes.
sub closestGenes {
	my ($chr, $w) = @_;						        ## chr, w
	unless (defined($g2pos{$chr})) {
		return "Unknown";                           ## return Unknown if chr is not in g2pos
	}
	my %dist 		= ();							## distance of genes to w.
	my $csv	 		= undef;					    ## csv of nearby gene info.
	my @genes	 	= @{$g2pos{$chr}};				## genes on chr sorted by 5'+.
	my $i 			= bsearch($w, \@genes );		## index of closest gene to w
	my @radius	= ($i - 20, $i + 20);				## consider +- 20 genes around i.
	my $gcount  = 0;								## gene count.
	if ($i < 0) {	return "Unknown"; }				## if $w not found, return Unknown

	##======= Step G1: compute min dist of each gene to window =========##
	foreach my $u (@genes[max(0,$radius[0]) .. min($#genes,$radius[1])]) {
		my ($d, $str, $id)	= @{$pos2g{$chr}{$u}};
		$dist{"$u:$d:$str:$id"}	= distToGene($w, $u, $d); 
	};

	##======= Step G2: sort genes asc by upstream,+ pos and get closest 1 to 2 genes =========##
	foreach my $key (sort {$dist{$a} <=> $dist{$b}} keys %dist) {
		if ($gcount == 2) { last; } 					# 2 genes max
		my ($u, $d, $str, $id) = split(":",$key);		# get gene info
		my ($uL, $dL) = ("upstream","downstream");	    # set up/down stream labels
		if (minus($str)){($uL,$dL) = ($dL,$uL); }		# if - stran, swap labels.
		my $dw = $dist{$key};							# min dist from gene to window.
		my ($mind, $loc) = !$dw  ? (0,"inside"):		# if dist 0, window inside gene
											 $w>$u ? ($w - $d, $dL):	        # if [u,d,w], label dL & w - d
															 ($u - $w, $uL);	# if [w,u,d], label uL & u - w
		$csv = defined($csv) ? "$csv," : $csv;			# add ',' before 2nd gene
		$csv = $csv."$id:$id:$mind:$loc";				# add gene info
		$gcount++;										# incr gene count
	}
	return $csv;
}

## Binary Search: Search an array a for a value x.
## Input: 	scalars : x (query number) and A (array reference to search).
## Output: 	scalar	: index of x in A or -1.
sub bsearch {
	my ($x, $A) = @_;            			  	 	# scan for x in ref to array, A.
	my ($l, $r) = (0, @$A - 1);  			 		# index of array A, assume A[0...|A|-1] is sorted asc..
	while ($l <= $r) {								# while our interval is valid
		if (($r - $l + 1) < 4) { return $l;} 	    # if 3 or less elem, return l.
		my $i = int(($l + $r)/2);					# let i be the midpt of [l,r]
		my $v = $A->[$i];							# let v be the value at midpt i
		if ($v == $x) { return $i;}					# if x is midpt, return its index i.
		($l, $r) = ($x > $v) ?						# if x > midpt val, scan A[i+1..r]
					($i + 1,$r) : ($l,$i - 1);      # if x < midpt val, scan A[l..i-1]
	}
	return -1;                      				# return -1 if x not found
}

