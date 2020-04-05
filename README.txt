###############################################################################################
#                                                                                             #
#                       SCRMshaw: Genome-wide CRM Prediction Program                          #
#                                                                                             #
#                 this version is the "HD" version (see Asma and Halfon, 2019)                #
#                        previous versions are not actively maintained                        #
#                 this version will work identical to earlier versions if run                 #
#                 with --lb option set to 0 (default)                                         #
#                                                                                             #
#              maintained by the Halfon lab (http://halfonlab.ccr.buffalo.edu)                #
#                                                                                             #
#-------------------------------Original Version credits:-------------------------------------#
#                                                                                             #
#                           Authors: Casey Hanson and Wei Yang                                #
#                                                                                             #
#                                                                                             #
#                           Email: weiyang4@illinois.edu                                      #
#                           Date: Sept 2014                                                   #
#                                                                                             #
###############################################################################################


    1. INTRODUCTION
    2. INSTALLATION
    3. RUNNING SCRMshaw
    4. REFERENCES


            1. INTRODUCTION
            ---------------

SCRMshaw is a genome-wide CRM prediction program (SCRMshaw.tar.gz), described and used in the following papers: Majid Kazemian et al 2011 and Miriam R. Kantorovitz and Majid Kazemian et al 2009. For details about the methods and procedures of this program please see the METHOD part of these two papers. The SCRMshaw program learns parameters from a list of CRMs that regulates the same spatial/temporal gene expression pattern, and predicts CRMs with similar functionality genome-wide. This verson of SCRMshaw program is designed for predicting CRMs in the genome sequence of model/non-model species given the known CRMs provided by users. Users can also choose three different scoring schemes as described in the papers to predict CRMs: IMM, PAC and HexMCD.

For a detailed protocol see Kazemian and Halfon (2019). Methods Mol. Biol. 1858:117-139. For details on the "HD" variant protocol see Asma and Halfon (2019). BMC Bioinformatics 20:174.



            2. INSTALLATION
            ---------------

1. UNPACK

The tar-archive contains one directory 'SCRMshaw' with the following sub-directories:

bin/
src/
code/
include/
example/

2. COMPILE (with GNU C++ compiler)

> cd SCRMshaw
> make

Now, the main SCRMshaw pipeline /SCRMshaw/code/scrm.pl should be able to run and print out the help information.

>perl /SCRMshaw/code/scrm.pl

            3. RUNNING SCRMshaw
            ---------------

1. INPUT:

SCRMshaw is easy to run and only has two required inputs: a target sequence FASTA file and a list of training set directories. 

1) The target sequence file. This file must be in uncompressed FASTA format (names of sequences should match sequence IDs in gff, exon and gene files if have them as input):
    
    >name_of_sequence_1
    GCAATCGGAGAGCGGGAGAGACGGAGAGAGTGCGCAGCCGCTGACTTTGA
    >name_of_sequence_2
    TATTAACTTTGGAATTGTGACAGCTGTTTAGCTTCTTAGCTCCGATGTGA

We highly recommend you to mask Tandem Repeats in your target sequence using Tandem Repeat Finder (Please keep reading for detail information, we provided an instruction from running the Repeat Finder below).

2) A list of Training set directories. provides the path to the training set directories, and each of your training set directory should contain two files: training CRM -> 'crms.fasta' and background sequence -> 'neg.fasta'. Your can add more than one training set into the list, but please make sure to have identical names for the two files within each training set directory: 'crms.fasta' and 'neg.fasta'(the current pipeline can only recognize files with the same names). Also, if you do not have background sequence 'neg.fasta', we provided you a script to generate it: /SCRMshaw/code/randomWithSameGC.pl. This script will take known CRM FASTA file and its source genome as input, and will randomly extract sequences from the genome with identitical GC content as CRM sequences. And we highly recommend you to mask Tandem Repeats in both of your CRM and background sequence. Here is a simply way to mask Tandem Repeats:

Step1. Download Tandem Repeat Finder from here: http://tandem.bu.edu/trf/trf.download.html
Step2. Run it with: TandemRepeatFinder InputFasta 2 7 7 80 10 50 500 -m -h

2. USAGE:

Highly Recommended OS: We have had issues with running it on MacOS, so we highly recommend you to run it on Linux/Ubuntu.
Before trying to run your own data, we highly recommend you to do a test run using the toy data in the 'example' directory. What you need to do is to go into 'example' directory, copy the command line in README file, paste the command line and run it. This would help to familiarize youself to SCRMshaw as well as the input and parameter settings.

>perl /SCRMshaw/code/scrm.pl [options] --genome TargetSequence --traindirlst TrainingSetList

Also, note that BioPerl modules: Bio::SeqIO is needed in SCRMshaw. Please make sure these modules have already been properly installed and are recognizable. The pipeline itself will also check the status of BioPerl modules.

Here is a useful instruction for installing BioPerl: http://www.bioperl.org/wiki/Installing_BioPerl_on_Unix

3. PARAMETERS:

        --genome <str>          target sequence file in FASTA format with one or multiple sequences, this file is required to perform the CRM prediction. We recommend you to mask Tandem Repeats in these sequences.
        --traindirlst <str>     a list of Training set directories, each training set directory should have two files: CRM sequence -> 'crms.fasta' and background sequence -> 'neg.fasta'. This is another required input for running SCRMshaw.
        --gff <str>             gff3 file, this file will be used to generate 'exons' and 'genes' files. 
        --exon <str>            exon file, format: scaffold_id, start, end, string, gene_id. This file will be used to mask exon regions in target sequence
        --gene <str>            gene file, format: scaffold_id, start, end, string, gene_id. This file will be used to report top-scored windows nearby gene locus as well as local rank of CRMs per gene
        --universeMap <str>     mapping between CRM training set and Universal Gene Set
        --wlen <int>            size of windows, default=500. The target sequence will be cut into window size sequences, and each window will be scored by scoring schemes.
        --wshift <int>          size of window shifts, default=250. This is the length of overlap between two windows and hence should be smaller than the value in '--wlen'.
        --lb <int>              left boundary i.e. starting point to take windows on each chromosome/scaffold, default = 0
        --wkmer <str>           size of 'k' for counting k-mer frequency in each window, default=1,6. Here we generate k-mer frequency files for PAC and HexMCD scoring schemes with k = 1 and 6. We do not recommed you to change this parameter unless you are an experienced user.
        --imm                   use IMM scoring scheme to score predicted CRMs.
        --hexmcd                use HexMCD scoring scheme to score predicted CRMs.
        --pac                   use PAC scoring scheme to score predicted CRMs.
        --thitw <int>           maximum number of hits for top-scored windows, default=0 (e.g. 2000)
        --thitg <int>           maximum number of genes nearest to top-scored windows, default=0 (e.g. 300)
        --distance <int>        distance between CRMs and its nearby genes (both of up and downstream), default=50000. For example, if the length of a gene is 1000, the program will look at +50000,1000,-50000, so it is a total of 101000 length region. This parameter is used to define the region for local rank of CRMs per gene. 
        --step <str>            choose the steps to be ran in the pipeline, default=123. Step1 is to process gff3 file (if you use gff as input) and target sequence file, step2 is to process the training data sets, and step3 is to run the scoring schemes. We recommend you to run all steps in the first run so that all directories needed in other steps can be created properly. After the first run, however, you can choose to run step 2 and 3 only to perform CRM prediction with added training sets and do not need to re-run the gff and genome preprocess steps.
        --outdir <str>          output directory, default=./ . You can redirect the output directory using this option. 
        --features <str>		a file containing a list of features (e.g. 'gene', 'ncRNA') to parse from the GFF file of genome annotations for the "gene" file; default='gene'.

4. OUTPUT:

The output will have sub-directories:

gff/
fasta/
scores/
hits/
training/

1) 'gff' stores exons and genes files.
2) 'fasta' stores target sequence, k-mer frequency and windows files of DNA sequence.
3) 'scores' stores the scores of CRMs from each scoring schemes.
4) 'hits' stores the top-scored windows
5) 'training' stores the training sequences and k-mer files

5. EXAMPLES:

More details about the meaning and format of output directories and files can be found in the 'example' directory, in which there is a toy dataset, output of the toy dataset and the README file. We highly recommend you to do the test run by using the toy dataset in this directory and please do read the README file.


         4. REFERENCES
         ---------------

Asma H and Halfon MS. Computational Enhancer Prediction: Evaluation and Improvements. BMC Bioinformatics 2019;20(1):174.

Kazemian M and Halfon MS. CRM Discovery Beyond Model Insects. Methods in Molecular Biology 2019;1858:117–39. 

Kazemian M, Zhu Q, Halfon MS, Sinha S. Improved accuracy of supervised CRM discovery with interpolated Markov models and cross-species comparison. Nucleic Acids Res. 2011;39(22):9463-72.

Kantorovitz MR, Kazemian M, Kinston S, et al. Motif-blind, genome-wide discovery of cis-regulatory modules in Drosophila and mouse. Dev Cell. 2009;17(4):568-79.
