#!/usr/bin/perl -w
#modified 08-09-2019
=head1 Description

    This is a genome-wide CRM prediction pipeline. It includes scoring schemes such as IMM, HexMCD and PAC.

=head1 Usage

    perl scrm.pl [options] --genome DNAFastaFile --traindirlst TrainingDirectoryList

    --genome <str>          a fasta file with DNA sequence from target species for CRM prediction, mandatory argument
    --traindirlst <str>     a directory list, each directory in the list stores fasta files as training set (crms.fasta and neg.fasta), mandatory argument
    --gff <str>             gff3 file, need to input both of exon and gene files if there is no gff3 file
    --exon <str>            exon file (mandatory argument if not using --gff), format: scaffold_id, start, end, string, gene_id
    --gene <str>            gene file (mandatory argument if not using --gff), format: scaffold_id, start, end, string, gene_id
    --universeMap <str>     mapping between CRM training set and Universal Gene Set
    --wlen <int>            length of windows seperated from the genome, default=500
    --wshift <int>          length of window shifts, default=250
    --wkmer <str>           kmer length to count kmer frequency in each window, default=1,6
    --imm                   use IMM scoring scheme
    --hexmcd                use HexMCD scoring scheme
    --pac                   use PAC scoring scheme
    --thitw <int>           maximum number of hits for top-scored windows, default=0 (e.g. 2000)
    --thitg <int>           maximum number of genes nearest to top-scored windows, default=0 (e.g. 300)
    --distance <int>        distance between CRMs and its nearby genes (both of up and downstream), default=50000
    --lb <int>              left boundary i.e. starting points to take windows on each chromosome, default = 0
    --step <str>            choose the steps to be ran in the pipeline, default=123. Step1 is to parse gff file and to process target sequence, step2 is to process the training data sets, and step3 is to run the scoring schemes.
    --outdir <str>          output directory, default=./
    --features <str>		a file containing a list of features (e.g. 'gene', 'ncRNA') to parse from the GFF file of genome annotations for the "gene" file; default='gene'
    --help                  output help information to screen

=head1 Example

    perl ./bin/scrm.pl --gff eve.gff3 --imm --hexmcd --pac --genome eve.location.fa --traindirlst trainingSet.lst

=cut

use strict;
use FindBin qw($Bin $Script);
use Getopt::Long;
use File::Basename;
use File::Spec;
use Data::Dumper;

my $t_start=time();

my ($Help, $Gff, $Gene, $Exon, $Imm, $Hexmcd, $Pac, $Genome, $Universe_orth, $Train_dir_lst, $feature_list);
my $Outdir = ".";
my ($Wlen, $Wshift) = (500, 250);
my ($Wkmer) = ("1,6");
my ($Thitw, $Thitg) = (0, 0);       ## Thitw = 0, Thitg = 0
my ($Plmer) = (1);      ## hard coded lmer for PAC
my ($Ysize) = (100);    ## ymf default parameters: kmer=6, size=100
my $Step = "123";       ## default is to run all steps
my ($Distance) = (50000);
my $UniverseMap = undef;
my $offset	= 	0;

GetOptions(
    "genome:s"=>\$Genome,
    "traindirlst:s"=>\$Train_dir_lst,
    "gff:s"=>\$Gff,
    "gene:s"=>\$Gene,
    "exon:s"=>\$Exon,
    "universeMap:s"=>\$UniverseMap,
    "imm"=>\$Imm,
    "hexmcd"=>\$Hexmcd,
    "pac"=>\$Pac,
    "outdir:s"=>\$Outdir,
    "wlen:i"=>\$Wlen,
    "wshift:i"=>\$Wshift,
    "wkmer:s"=>\$Wkmer,
    #  "plmer:i"=>\$Plmer,      ##  not allow change of lmer to prevernt ambiguity
    "thitw:i"=>\$Thitw,
    "thitg:i"=>\$Thitg,
    "step:s"=>\$Step,
    "distance:i"=>\$Distance,
    "lb:i"=>\$offset,	##left boundary, i.e., starting points to take windows on chr
    "features:s"=>\$feature_list,
    "help"=>\$Help
);


my $Kmer = (split(",",$Wkmer))[1];      ## Kmer = 6

die `pod2text $0` if (!$Genome && !$Train_dir_lst || $Help || ($Thitw == 0 && $Thitg == 0) || ($Thitw > 0 && $Thitg > 0));
die ("Aborting because the required argument --genome wasn't specified! Use --help more for help information.") if (!$Genome);
die ("Aborting because the required argument --traindirlst wasn't specified! Use --help for more help information.") if (!$Train_dir_lst);


##============ check the existence of the Bio Perl modules ===========##
my $require_package = eval
{
    require Bio::SeqIO;
    # require Bio::DB::Fasta;
    Bio::SeqIO->import();
    #  Bio::DB::Fasta->import();
    1;
};

if ($require_package){
    warn "Bio::SeqIO loaded and imported successfully\n";
    # warn "Bio::DB::Fasta loaded and imported successfully\n";
}
else{
    die("Bio::SeqIO cannot be imported. Please check to see if the modules have already been installed.\n");
}

my $hexmcd_param = "-rc";
my $pac_param = "-rc";

`mkdir $Outdir` unless (-d "$Outdir");
`mkdir $Outdir/fasta` unless (-d "$Outdir/fasta");
`mkdir $Outdir/fasta/chr` unless (-d "$Outdir/fasta/chr");
`mkdir $Outdir/fasta/kmers` unless (-d "$Outdir/fasta/kmers");
`mkdir $Outdir/fasta/windows` unless (-d "$Outdir/fasta/windows");
`mkdir $Outdir/gff` unless (-d "$Outdir/gff");
`mkdir $Outdir/scores` unless (-d "$Outdir/scores");
`mkdir $Outdir/hits` unless (-d "$Outdir/hits");
`mkdir $Outdir/training` unless (-d "$Outdir/training");
`mkdir $Outdir/scores/imm` unless (-d "$Outdir/scores/imm");
`mkdir $Outdir/hits/imm` unless (-d "$Outdir/hits/imm");
`mkdir $Outdir/scores/hexmcd` unless (-d "$Outdir/scores/hexmcd");
`mkdir $Outdir/hits/hexmcd` unless (-d "$Outdir/hits/hexmcd");
`mkdir $Outdir/scores/pac` unless (-d "$Outdir/scores/pac");
`mkdir $Outdir/hits/pac` unless (-d "$Outdir/hits/pac");

##============ convert to absolute file paths (added Oct 2019) ===========##

if (defined $Gff){
	$Gff = File::Spec->rel2abs($Gff);
}	

if (defined $Exon){
	$Exon = File::Spec->rel2abs($Exon);
}

if (defined $Gene){
	$Gene = File::Spec->rel2abs($Gene);	
}

$Genome = File::Spec->rel2abs($Genome);	

if (defined $feature_list){
	$feature_list = File::Spec->rel2abs($feature_list);	
} else {
	$feature_list = "gene";
}		

##================== Step 1 ===================##
if ($Step =~ /1/){
    warn "Step1: preprocessing gff file and target sequence...\n";
    ##========== parse gff file ===============##
    if (defined $Gff){      ## if has gff file as input
        warn "Parsing gff3 file...\n";
        ##======== generate exon file ===========##
        `perl $Bin/../code/gff3.pl $Gff exon > $Outdir/gff/exons`;
        $Exon = "$Outdir/gff/exons";
        ##======== generate gene file ===========##
        `perl $Bin/../code/gff3.pl $Gff $feature_list > $Outdir/gff/genes`;
        $Gene = "$Outdir/gff/genes";
    }
    else{   ## if no gff file as input
        ##======= parse exon file ============##
        if (defined $Exon){     ## if exon is defined
            `cp $Exon $Outdir/gff/exons`;
            warn "Found exon file, it will be used to mask exons in chr.\n";
        }
        ##======== parse gene file =============##
        if (defined $Gene){     ## if gene is defined
            `cp $Gene $Outdir/gff/genes`;
            warn "Found gene file, it will be used to generate hits file.\n";
        }
    }

    ##========== mask exon in genome seq, and separate into chr ============##
    warn "Seperating target sequence...\n";
    if (defined $Exon){         ## if exon file is defined
        `perl $Bin/../code/genome2chr.NoRepMsk.pl --exon $Exon $Genome $Outdir/fasta/chr`;
    }
    else{       ## if exon file is undef
        `perl $Bin/../code/genome2chr.NoRepMsk.pl $Genome $Outdir/fasta/chr`;
    }
    ##========== create windows and kmers for chr ============##
    warn "Creating window and kmer files for target sequence...\n";
    `perl $Bin/../code/windows.pl -win $Wlen -shift $Wshift -kmer $Wkmer -kmerc $Wkmer -o $offset $Outdir/fasta/chr $Outdir/fasta/windows $Outdir/fasta/kmers`;
}

##=========== record training set, and this is done outside of step condition ===========##
my %train = ();     ## hash training set
open (TSD,$Train_dir_lst) || die ("Couldn't open training directory list: $Train_dir_lst\n");
while (<TSD>){
    chomp;
    my $line=$_;
    my $dir_name=basename($line);   ## e.g. dir_name = mapping0.ap
    $train{$dir_name}=1;        ## e.g. train{mapping0.ap} = 1
}

##========== mapping between training set and universal gene set ===========##
my %crm2uni = ();           ## hash CRM train set and its corresponding universal gene set
if (defined $UniverseMap){
    open (UNI,$UniverseMap) || die ("Couldn't open $UniverseMap\n");
    while (<UNI>){
        chomp(my $line = $_);
        my ($crm, $univerFile) = split /\s+/,$line;
        $crm2uni{$crm} = $univerFile;       ## crm2uni{mapping0.ap} = dmel_universe_FB_cg.genes
        #warn "crm: $crm, universal gene set: $univerFile\n";
    }
    close UNI;
}

##===================== Step 2 ====================##
if ($Step =~ /2/){
    warn "Step2: processing training data...\n";
    open (LST,$Train_dir_lst) || die ("Couldn't open training directory list: $Train_dir_lst\n");
    while (<LST>){
        chomp;
        my $line=$_;
        my $dir_name=basename($line);       ## e.g. dir_name = mapping0.ap
        warn "Creating training set $dir_name...\n";
        my $k1=(split(",",$Wkmer))[0];    ## assume only two kmers are alllowed, default=1 and 6
        my $k6=(split(",",$Wkmer))[1];    ## assume only two kmers are alllowed, default=1 and 6
        `mkdir $Outdir/training/$dir_name` unless (-d "$Outdir/training/$dir_name");
        `cp $line/crms.fasta $Outdir/training/$dir_name`;
        `cp $line/neg.fasta $Outdir/training/$dir_name`;
        `perl $Bin/../code/ymf.pl --kmer $Kmer --size $Ysize $line/crms.fasta $line/neg.fasta $Outdir/training/$dir_name`;   ## ymf default parameters: kmer=6, size=100
        `perl $Bin/../code/kmers.pl -r $line/crms.fasta $k1 > $Outdir/training/$dir_name/crms.fasta.rc.$k1`;
        `perl $Bin/../code/kmers.pl -r $line/crms.fasta $k6 > $Outdir/training/$dir_name/crms.fasta.rc.$k6`;
        `perl $Bin/../code/kmers.pl -r $line/neg.fasta $k6 > $Outdir/training/$dir_name/neg.fasta.rc.$k6`;
        `perl $Bin/../code/kmers_sum.pl $Outdir/training/$dir_name/crms.fasta.rc.$k1 $k1 > $Outdir/training/$dir_name/crms.fasta.sum.rc.$k1`;
        `perl $Bin/../code/kmers_sum.pl $Outdir/training/$dir_name/crms.fasta.rc.$k6 $k6 > $Outdir/training/$dir_name/crms.fasta.sum.rc.$k6`;
        `perl $Bin/../code/kmers_sum.pl $Outdir/training/$dir_name/neg.fasta.rc.$k6 $k6 > $Outdir/training/$dir_name/neg.fasta.sum.rc.$k6`;
        `perl $Bin/../code/kmers.pl -r $Outdir/training/$dir_name/ymf.$Ysize.fasta $Kmer > $Outdir/training/$dir_name/ymf.$Ysize.rc.$Kmer`;
        `perl $Bin/../code/kmers_sum.pl $Outdir/training/$dir_name/ymf.$Ysize.rc.$Kmer $Kmer > $Outdir/training/$dir_name/ymf.$Ysize.sum.rc.$Kmer`;
    }
}

##===================== Step 3 ========================##
if ($Step =~ /3/){
    warn "Step3: Scoring CRMs...\n";
    $Gene = "$Outdir/gff/genes" if $Gff;        ## find gene file if only run step 3
    foreach my $train_set (sort keys %train){   ## e.g. train_set = mapping0.ap
        my $k1=(split(",",$Wkmer))[0];      ## k1 = 1
        my $k6=(split(",",$Wkmer))[1];      ## k2 = 6

        ##================= IMM ===================##
        if (defined $Imm){
            warn "IMM scoring $train_set...\n";
            `mkdir $Outdir/scores/imm/$train_set` unless (-d "$Outdir/scores/imm/$train_set");
            `mkdir $Outdir/hits/imm/$train_set` unless (-d "$Outdir/hits/imm/$train_set");
            `perl $Bin/../code/immNotChr.pl -kmer $Kmer -win $Wlen $Outdir/training/$train_set/crms.fasta $Outdir/training/$train_set/neg.fasta $Outdir/scores/imm/$train_set $Outdir/fasta/windows`;
            warn "Getting top scored CRMs from IMM $train_set results...\n";
            if (defined $Gene) {
                warn "Gene is defined: ".$Gene."\n";
                if (defined $crm2uni{$train_set}) {
                    warn "Universe gene set is defined\n";
                    `perl $Bin/../code/top_scores.ReAsgnNrbGene.pl -universe $crm2uni{$train_set} -win $Wlen -topw $Thitw -topg $Thitg $Outdir/scores/imm/$train_set $Outdir/fasta/windows $Outdir/gff/genes $Outdir/hits/imm/$train_set/$train_set.hits`;
                } else {
                    warn "Universe gene set is not defined\n";
                    `perl $Bin/../code/top_scores.ReAsgnNrbGene.pl -win $Wlen -topw $Thitw -topg $Thitg $Outdir/scores/imm/$train_set $Outdir/fasta/windows $Outdir/gff/genes $Outdir/hits/imm/$train_set/$train_set.hits`;  
                }

                `perl $Bin/../code/ranking.hit.pl --distance $Distance $Outdir/hits/imm/$train_set/$train_set.hits > $Outdir/hits/imm/$train_set/$train_set.hits.rankLst`;
                `perl $Bin/../code/add.rank2hits.pl $Outdir/hits/imm/$train_set/$train_set.hits.rankLst $Outdir/hits/imm/$train_set/$train_set.hits > $Outdir/hits/imm/$train_set/$train_set.hits.ranked`;
            } else {
                $Thitw = 2000 if $Thitw == 0;     ## set Thitw = 2000 if Thitw is undef
                `perl $Bin/../code/top_windows.pl -win $Wlen -top $Thitw $Outdir/scores/imm/$train_set $Outdir/fasta/chr > $Outdir/hits/imm/$train_set/$train_set.hits`;
            }
        }

        ##==================== HexMCD ===================##
        if (defined $Hexmcd){
            warn "HexMCD scoring $train_set...\n";
            `mkdir $Outdir/scores/hexmcd/$train_set` unless (-d "$Outdir/scores/hexmcd/$train_set");
            `mkdir $Outdir/hits/hexmcd/$train_set` unless (-d "$Outdir/hits/hexmcd/$train_set");
            `perl $Bin/../code/hexmcd.pl $hexmcd_param -kmer $Kmer -win $Wlen -shift $Wshift $Outdir/training/$train_set/crms.fasta.sum.rc.$k6 $Outdir/training/$train_set/neg.fasta.sum.rc.$k6 $Outdir/scores/hexmcd/$train_set $Outdir/fasta/chr $Outdir/fasta/kmers`;
            warn "getting top scored CRMs from HexMCD $train_set results...\n";
            if (defined $Gene){
                warn "Gene is defined: ".$Gene."\n";
                if (defined $crm2uni{$train_set}){
                    warn "Universe gene set is defined\n";
                    `perl $Bin/../code/top_scores.ReAsgnNrbGene.pl -universe $crm2uni{$train_set} -win $Wlen -topw $Thitw -topg $Thitg $Outdir/scores/hexmcd/$train_set $Outdir/fasta/windows $Outdir/gff/genes $Outdir/hits/hexmcd/$train_set/$train_set.hits`;
                }
                else{
                    warn "Universe gene set is not defined\n";
                    `perl $Bin/../code/top_scores.ReAsgnNrbGene.pl -win $Wlen -topw $Thitw -topg $Thitg $Outdir/scores/hexmcd/$train_set $Outdir/fasta/windows $Outdir/gff/genes $Outdir/hits/hexmcd/$train_set/$train_set.hits`;    
                }
                `perl $Bin/../code/ranking.hit.pl --distance $Distance $Outdir/hits/hexmcd/$train_set/$train_set.hits > $Outdir/hits/hexmcd/$train_set/$train_set.hits.rankLst`;
                `perl $Bin/../code/add.rank2hits.pl $Outdir/hits/hexmcd/$train_set/$train_set.hits.rankLst $Outdir/hits/hexmcd/$train_set/$train_set.hits > $Outdir/hits/hexmcd/$train_set/$train_set.hits.ranked`;
            }
            else{
                $Thitw = 2000 if $Thitw == 0;     ## set Thitw = 2000 if Thitw is undef
                `perl $Bin/../code/top_windows.pl -win $Wlen -top $Thitw $Outdir/scores/hexmcd/$train_set $Outdir/fasta/chr > $Outdir/hits/hexmcd/$train_set/$train_set.hits`;
            }
        }
        
        ##====================== PAC ==================##
        if (defined $Pac){
            warn "PAC scoring $train_set...\n";
            `mkdir $Outdir/scores/pac/$train_set` unless (-d "$Outdir/scores/pac/$train_set");
            `mkdir $Outdir/hits/pac/$train_set` unless (-d "$Outdir/hits/pac/$train_set");
            `perl $Bin/../code/pac.pl $pac_param -kmer $Kmer -nmer $Plmer -win $Wlen -shift $Wshift $Outdir/training/$train_set/crms.fasta.sum.rc.$k6 $Outdir/training/$train_set/crms.fasta.sum.rc.$k1 $Outdir/training/$train_set/ymf.$Ysize.sum.rc.$Kmer $Outdir/scores/pac/$train_set $Outdir/fasta/chr $Outdir/fasta/kmers`;
            warn "getting top scored CRMs from PAC $train_set results...\n";
            if (defined $Gene){
                warn "Gene is defined: ".$Gene."\n";
                if (defined $crm2uni{$train_set}){
                    warn "Universe gene set is defined\n";
                    `perl $Bin/../code/top_scores.ReAsgnNrbGene.pl -universe $crm2uni{$train_set} -win $Wlen -topw $Thitw -topg $Thitg $Outdir/scores/pac/$train_set $Outdir/fasta/windows $Outdir/gff/genes $Outdir/hits/pac/$train_set/$train_set.hits`;
                } else {
                    warn "Universe gene set is not defined\n";
                    `perl $Bin/../code/top_scores.ReAsgnNrbGene.pl -win $Wlen -topw $Thitw -topg $Thitg $Outdir/scores/pac/$train_set $Outdir/fasta/windows $Outdir/gff/genes $Outdir/hits/pac/$train_set/$train_set.hits`;
                }
                `perl $Bin/../code/ranking.hit.pl --distance $Distance $Outdir/hits/pac/$train_set/$train_set.hits > $Outdir/hits/pac/$train_set/$train_set.hits.rankLst`;
                `perl $Bin/../code/add.rank2hits.pl $Outdir/hits/pac/$train_set/$train_set.hits.rankLst $Outdir/hits/pac/$train_set/$train_set.hits > $Outdir/hits/pac/$train_set/$train_set.hits.ranked`;
            }
            else{
                $Thitw = 2000 if $Thitw == 0;     ## set Thitw = 2000 if Thitw is undef 
                `perl $Bin/../code/top_windows.pl -win $Wlen -top $Thitw $Outdir/scores/pac/$train_set $Outdir/fasta/chr > $Outdir/hits/pac/$train_set/$train_set.hits`;
            }
        }
    }
}
warn "Finished!\n";

my $t_end=time();
my $duration=$t_end-$t_start;
print "total time: $duration seconds\n";
