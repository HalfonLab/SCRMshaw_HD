#!/usr/bin/perl -w

=head1 Description

    This script accepts directories of scores, sequences file
    and outputs them in descending order in FASTA format.
    This number is typically 2000, but can be set by the user.
    This version of top_score doesn't need universal gene set.

=head1 Usage

    perl top_scores.pl [OPTIONS] ODIR SDIR GFF

    ODIR     Directory containing scores for FASTA files
    SDIR     Directory containing sequence files

    --help              Print help and exit.
    --top <int>         Maximum number of hits. Default 2000.
    --win <int>         Window length. Default 500.

=cut


use strict;
use File::Basename;
use Getopt::Long;

my %score;
my $Win = 500;
my $Top = 2000;
my $help;


GetOptions(
    "win:i"=>\$Win,
    "top:i"=>\$Top,
    "help"=>\$help
);


my $score_dir = $ARGV[0];
my $chr_dir = $ARGV[1];

die `pod2text $0` if (@ARGV != 2 || $help);

while (<$score_dir/*.fasta>)
{
    chomp(my $file=$_);
    my $chr=basename($file);
    $chr =~ s/\.fasta//g;
    open SFILE,$file or die("cannot open score file: $file\n");
    while (<SFILE>)
    {
        chomp(my $line=$_);
        my ($w, $s)=split(" ",$line);
        $score{"$chr:$w"}=$s;
    }
    close SFILE;
}

my $count=0;
for my $key (sort {$score{$b} <=> $score{$a}} keys %score)
{
    last unless ($count < $Top);
    my ($chr,$w) = split(":",$key);
    if (-e "$chr_dir/$chr.fasta")
    {
        my $cfile="$chr_dir/$chr.fasta";
        open FA,$cfile or die("cannot open fasta file: $cfile\n");
        $/=">";<FA>;$/="\n";
        while (<FA>)
        {
            my $title=$_;
            my $id=(split(" ",$title))[0];
            $/=">";
            chomp(my $seq=<FA>);
            $/="\n";
            $seq =~ s/\s+//g;
            chomp $seq;
            my $crm = substr($seq, $w, $Win);
            print ">$key\t$score{$key}\n$crm\n";
        }
        close FA;
    }
    $count++;
}

