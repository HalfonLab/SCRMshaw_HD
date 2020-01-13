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

my %feat_hash = (); 		## hash of features to accept

#check if there is a feature_list and make hash of features
if ($elem eq 'gene'){
	$feat_hash{'gene'}++;
	
} elsif ($elem eq 'exon') {
	$feat_hash{'exon'}++;
	
} else {
	open(FEATURES, "<$elem") or die "Unable to open features file: $!";
	while(<FEATURES>){
		chomp;
		$feat_hash{$_}++;  #store the desired features
	}	
	
	close(FEATURES);
}	


# Open gff file
open(IFILE, $ifile) or die("Couldn't open gff file: $ifile for reading\n");

while(my $line = <IFILE>){
	chomp($line);
    next if $line =~ /^#/;
	## chr  source  element  start  stop  strand . - . ID=######
	my ($seqid,  $source, $type, 
			$start,  $end,    $score,	
			$strand, $phase,	$attr)	= split('\s+', $line);
	
	#if an accepted feature type, parse the line
	if (exists $feat_hash{$type}) {   	
		my $chr = $1 if $seqid 	=~ /(\S+)/;
		my $id = "";
		if($attr =~ /;/){
			$id  = $1 if $attr =~ /[ID|Name]=(\S+?);/;
		} else {
			$id  = $1 if $attr =~ /[ID|Name]=(\S+)$/;
		}

		print "$chr\t$start\t$end\t$strand\t$id\n";
		
	} else {
		#print "$type not recognized\n";
		next; 	## if not indicated type, next
	}		
}

close(IFILE);

