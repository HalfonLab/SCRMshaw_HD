#!/usr/bin/perl

=head1 Description

    This script accepts directories of scores, sequences file
    and a universe file of genes along with a GFF of gene locations
    and finds the top number of highest scoring sequences along with nearby genes
    and outputs them in descending order in FASTA format.
    This number is typically 2000, but can be set by the user.
    This version of top_score doesn't need universal gene set.

=head1 Usage

    perl top_scores.pl [OPTIONS] ODIR SDIR GFF
    
    ODIR     Directory containing scores for FASTA files
    SDIR     Directory containing sequence files
    GFF      GFF file with gene locations in species
    
    --help              Print help and exit.
    --top <int>         Maximum number of hits. Default 2000.
    --win <int>         Window length. Default 500.

=cut

use strict;
use FindBin qw($Bin $Script);
use lib '$Bin/../lib/';
use Bio::SeqIO;
use Bio::DB::Fasta;
use Getopt::Long;
use File::Basename;

my $help 	= undef;
my $fb   	= undef;
my $verb 	= undef;
my $wlen	= 500;
my $mhit    = 2000;
#my ($start_run, $end_run, $run_time);
#$start_run=time();

GetOptions(
    "h|help" => \$help, 
    "w|win=i" => \$wlen, 
    "t|top=i" => \$mhit
);

my $odir 	= $ARGV[0];  	    ## Directory of window scores.
my $sdir 	= $ARGV[1]; 	    ## Directory of sequences files.
#my $ufile   = $ARGV[2];		    ## Universe file.
my $gff  	= $ARGV[2];		    ## GFF3 genes files.
my %g2pos	= ();				## g2pos{chr} array of geneids sorted by 5'+pos
my %pos2g	= ();				## pos2g{chr}{pos} = (3'+pos, strand, geneid)
my %score   = ();					## score{chr:x} = score.
my %uhit	= ();				## Universe genes has: 1/0 if reported already/not.
my $nhit 	= 0;				## Total hits reported.

## Check exit conditions
die `pod2text $0` if (@ARGV != 3 || $help);

## Step 1 : Read universe geneset.
#open(IFILE, "<$ufile" ) or die("Can not open $ufile for reading. Exit.\n");
#while (<IFILE> ) {
#	chomp($_);
#	$uhit{$_} = 1;	## 0 indicates gene hasn't been reported
#}
#close(IFILE);

## Step 2: Parse gff file.
# aaegcont1.4(chr) VectorBase gene 1(b) 20(e) . -(sense) . ID=AAEL021(id)
open(IFILE, "<$gff ") or die ("Couldn't open $gff  for reading\n");
while (<IFILE>) {
	chomp($_);
	my ($chr, $u, $d, $str, $id) = split(/\s+/, $_);
    ($u, $d) = (min($u, $d), max($u, $d));
    # ($u, $d) = ($d, $u) if ($u > $d);
	@{$pos2g{$chr}{$u}} 	= ($d, $str, $id);
	push (@{$g2pos{$chr}}, $u);
}
close(IFILE);

## Step 3: Sort g2pos by upstream,+ pos on chr.
foreach my $chr (keys %g2pos) { 
	@{$g2pos{$chr}} = sort { $a <=> $b } @{ $g2pos{$chr} };
}

## Step 4: Read score for each chr.
while (<$odir/*.fasta>) {
	my $file = $_;
	chomp($file);
	$file =~ /([^\/]*)$/;								## match /....$
	my $chr = $1;												## set to match
	$chr 		=~ s/\.fasta//g;						## remove .fasta suffix
	open(IFILE, "<$file") or die("Couldn't open $file for reading\n");
	while (<IFILE>) {
		chomp;														## window score
		my ($w, $score) 	 = split(/\s+/);## 0   -3.20086 -0.00640 500
		$score{"$chr:$w"} = $score;			## score{chr:w} = score
	}
	close(IFILE);
}

## Step 5: Sort all windows by score in descending order.
##	Until nhit = mhit or no more windows, test windows for hits.
##	A window w should be reported as a hit if : 
##	 	1. w's nearest gene is in the universe  AND
##		2. w's nearest gene hasn't been reported yet with any window.
##	To report window, print the following info in FASTA format:
##		w's position, score, closest gene inf, and sequence
foreach my $key(sort{ $score{$b} <=> $score{$a}} keys %score) {
	if ($nhit >= $mhit) {last;}						## don't output more than 200
	my ($chr, $w)	= split(/\:/, $key );		## chr and window pos
	my $ifile = "$sdir/$chr.fasta";				## exon/repeat masked chr file name.
	unless(-e "$ifile") { 								## if chr doesn't exist
		print STDERR "$ifile doesn't exist. Skipping.\n";
		next; 
	}
	my ($id, $desc, $seq) = ReadFASTA("$ifile"); # read FASTA file
	my $wseq		= substr($seq, $w, $wlen);	## seq[w ... w + wlen - 1]  # use Bio::database?
	my $score 	= $score{$key};							## window score
	my $csv			= closestGenes($chr,$w);		## get 2 closest genes to window.
	if ($csv 		=~ /Unknown/) { next; }			## skip if gene unknown.
	my @colsep 	= split(/\:/, $csv);				## split by :
	my $gene		= $colsep[0];

##	if ($uhit{$gene} == 1) {							## if gene unreported & in universe
	if (!$uhit{$gene})
    {
        $uhit{$gene} = 2;									## mark gene as reported (1)
		$nhit = $nhit + 1;							## incr hits and print hit to STDOUT
		print STDOUT ">$key\t$score\t$csv\n"; ## ex:>chr:w score neighbor-info
		print STDOUT "$wseq\n";								## ex:seq[w ... w + wlen - 1]
	}
}

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

## Input:		Scalars chr and w which are the chromosome and pos of window w.
## Output:	String with info on closest gene or gene pair (order by proximity)
## Example: 500:1000:+:AAEL0021,2000:2500:+:AAEL0883
## Notes:
##			The function may return only one gene.
##			The closest gene is start=500,  stop=100,  sense = +, id = AAEL0021
##			2nd closest gene is start=2500, stop=2500, sense = +, id = AAEL0883
##			":" breaks fields within a gene.
##			"," breaks between genes.
sub closestGenes {
	my ($chr, $w) = @_;										## Window chr and start pos.
	unless (defined($g2pos{$chr})) {
		return "Unknown";
	}
	my %dist 		= ();												## dist of genes to w.
	my $csv	 		= undef;										## csv of nearby gene info.
	my @genes	 	= @{$g2pos{$chr}};				## genes on chr sorted by u.
	my $i 			= bsearch($w, \@genes );		## index of $w or closest in genes.
	my @radius	= ($i - 20, $i + 20);				## consider +- 20 genes around i.
	my $gcount  = 0;												## gene count.
	if ($i < 0) {	return "Unknown"; }				## if $w not found, return Unknown

	## Step G1: compute min dist of each gene to window
	foreach my $u (@genes[max(0,$radius[0]) .. min($#genes,$radius[1])]) {
		my ($d, $str, $id)	= @{$pos2g{$chr}{$u}};
		$dist{"$u:$d:$str:$id"}	= distToGene($w, $u, $d); 
	};

	## Step G2: sort genes asc by upstream,+ pos and get closest 1 to 2 genes.
	foreach my $key (sort {$dist{$a} <=> $dist{$b}} keys %dist) {
		if ($gcount == 2) { last; } 								# 2 genes max
		my ($u, $d, $str, $id) = split(":",$key);		# get gene info
		my ($uL, $dL) = ("upstream","downstream");	# set up/down stream labels
		if (minus($str)){($uL,$dL) = ($dL,$uL); }		# if - stran, swap labels.
		my $dw = $dist{$key};												# min dist from gene to window.
		my ($mind, $loc) = !$dw  ? (0,"inside"):		# if dist 0, window inside gene
											 $w>$u ? ($w - $d, $dL):	# if [u,d,w], label dL & w - d
															 ($u - $w, $uL);	# if [w,u,d], label uL & u - w
		$csv = defined($csv) ? "$csv," : $csv;			# add ',' before 2nd gene
		$csv = $csv."$id:$id:$mind:$loc";						# add gene info
		$gcount++;																	# incr gene count
	}
	return $csv;
}

## Binary Search: Search an array a for a value x.
## Input: 	scalars : x (query number) and A (array reference to search).
## Output: 	scalar	: index of x in A or -1.
sub bsearch {
	my ($x, $A) = @_;            			  	 	# scan for x in ref to array, A.
	my ($l, $r) = (0, @$A - 1);  			 			# assume A[0...|A|-1] is sorted asc..
	while ($l <= $r) {											# while our interval is valid
		if (($r - $l + 1) < 4) { return $l;} 	# if 3 or less elem, return l.
		my $i = int(($l + $r)/2);							# let i be the midpt of [l,r]
		my $v = $A->[$i];											# let v be the value at midpt i
		if ($v == $x) { return $i;}						# if x is midpt, return its index i.
		($l, $r) = ($x > $v) ?								# if x > midpt val, scan A[i+1..r]
							 ($i + 1,$r) : ($l,$i - 1); # if x < midpt val, scan A[l..i-1]
	}
	return -1;                      				# return -1 if x not found
}

#$end_run=time();
#$run_time=$end_run - $start_run;
#warn "top_window finished in $run_time seconds\n";
