#!/usr/bin/perl

=head1 Description

This file parses a gff file and selects a particular element.

Output format (TAB): chr  start  end  strand  id

=head1 Usage

perl gff3.pl <gff3> <exon|gene> > <output>

=cut

use strict;

## Arguments
my $ifile = $ARGV[0];		## Input file.
my $elem =	$ARGV[1];		## Element to filter.

## Check exit conditions
die `pod2text $0` if (@ARGV != 2);


# Open gff file
open(IFILE, $ifile) or die("Couldn't open gff file: $ifile for reading\n");
while(my $line = <IFILE>){
	chomp($line);
    next if $line =~ /#/;
	## chr  source  element  start  stop  strand . - . ID=######
	my ($seqid,  $source, $type, 
			$start,  $end,    $score,	
			$strand, $phase,	$attr)	= split('\s+', $line);
	if ($type ne $elem) { next; }		## if not indicated type, next
	my $chr = $1 if $seqid 	=~ /(\S+)/;

	my $id = "";
	if($attr =~ /;/){
	    $id  = $1 if $attr =~ /[ID|Name]=(\S+?);/;
	} else {
	    $id  = $1 if $attr =~ /[ID|Name]=(\S+)$/;
	}

#    my $id  = $1 if $attr =~ /[ID|Name|Parent]=(\S+?)/; #old line
	print "$chr\t$start\t$end\t$strand\t$id\n";
}
close(IFILE);

