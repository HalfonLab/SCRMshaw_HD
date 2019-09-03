#!/usr/bin/perl

## This script accepts an in-paranoid file and outputs a homology file
## containing ortho groups in the following format:
## GROUP-NUMBER SPC-ID GENE-ID

## Set script conditions
use strict;
use lib '/home/sinhas/sinhas/lib/';
use File::Basename;
use Getopt::Long;

my $help	= undef;		##	Help option
my $spc		=	"ALT";		## 	Species Name

## Get Options
GetOptions(
	"h|help" 		=> \$help,
	"s|spc=s"   => \$spc
);

my $hfile	= $ARGV[0];		## Homology file
my $argc  = $#ARGV + 1;	## Number of arguments given to script.
$spc = uc($spc);				## Conver to uppercase

## Check exit conditions
if (defined($help)) {
	Usage();
	exit(-1);
} elsif ($argc != 1) {
	print "Incorrect # of args given: $argc\n";
	Usage();
	exit(-1);
}

open(HFILE, "<$hfile") or die "Can not open $hfile.\n";
my $gnum = 0;
while(my $line = <HFILE>){
 chomp($line);
 if ($line =~ /Group of orthologs \#(\d+)\. Best.*/){
    $gnum = $1;				## 	group number
		next;
	} elsif ($line =~ /Bootstrap/ || $line =~ /Score/ || $line !~ /%/) {
		next;
	} 
	my @sep = split('\s+', $line);
	my $len	=	$#sep + 1;
	my ($sgene, $sperc, $dgene, $dperc) = @sep;
	if ($len == 4){	
		if ($sperc =~ /100/) {
			print "$gnum\t$spc\t$sgene\n";
		}
		if ($dperc =~ /100/) {
			print "$gnum\tDMEL\t$dgene\n";
		}
	} elsif ($len == 2) {
		if ($sperc =~ /100/) {
			print "$gnum\t$spc\t$sgene\n";
		}
	} 
}
close(HFILE);

sub Usage {
	my ($name, $path, $suffix) = fileparse($0);
	print "Usage: $name [OPTIONS] FILE \n";
	print " FILE     Input Fasta File.\n";
	print " -h       Print help and exit\n";
	print " -s=spc   Species name. Automatically upperased. Default ALT\n";
}
