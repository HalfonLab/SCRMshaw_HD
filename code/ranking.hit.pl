#!/usr/bin/perl -w

=head1 Description

    This script is used to generate the local rank list for hits file.
    Output format: GeneID, CRM, Score, Rank 

=head1 Usage
    
    perl ranking.hit.pl [Options] HitsFile

    HitsFile    Output file from scoring schemes, this file can be found in hits/ directory.

    --distance <int>    Set the distance from CRMs to gene locus, this scripts will look at both of up- and downstream regions, default=50000

=cut


use strict;
use File::Basename;
use Getopt::Long;
#use Bio::SeqIO;

my $Distance=50000;

GetOptions(
    "distance:i"=>\$Distance
);

die `pod2text $0` if ( @ARGV != 1);

my $hits=$ARGV[0];
my %gene;
my %cand_crm;
my $file_name=basename($hits);

#my $fo = Bio::SeqIO->new(-file=>$hits,-format=>'Fasta');
open FA,$hits;
$/ = ">";<FA>;$/ = "\n";
#while (my $f =  $fo->next_seq()) {
while (<FA>) { 
    my $title = $_;
    #my $title = $f->desc();
    #print "$title\n";
    # my $seq = $f->seq();
    $/ = ">";
    my $seq = <FA>;
    chomp $seq;
    $/ = "\n";
    chomp $seq;
    my @array=split /\s+/,$title;
    my $score=$array[1];
    my $crm_id=$array[0];
    my @gene_unit=split /,/,$array[2];
    foreach (@gene_unit)
    {
        my $unit=$_;
        my @gene_loc=split /\:/,$unit;
        my $gene_id=$gene_loc[0];
        my $gene_dist=$gene_loc[2];
        if ($gene_dist <= $Distance)
        {
            push @{$gene{$gene_id}{$score}},$crm_id;
        }
    }
}
close FA;
#open OUT1, ">$file_name.globalRank.lst";

for my $gene_id (sort keys %gene)
{
    my $count=1;
    for my $score (sort {$b <=> $a} keys %{$gene{$gene_id}})
    {
        my $gene_num=0;
        foreach (@{$gene{$gene_id}{$score}})
        {
            my $crm_id=$_;
            print "$gene_id\t$crm_id\t$score\t$count\n";
            $gene_num++;
        }
        $count+=$gene_num;
    }
}
