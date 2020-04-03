Toy example readme page -- Sept 2014
    Modified -- Mar 2015 

This directory includes a toy dataset to run SCRMshaw.
===========================================================================

Input dataset:

___ ./test_data/            
                |___ 'eve.location.fa'
                |
                |___ 'trainingSet.lst'
                |
                |___ 'mapping0.ap'
                |
                |___ 'eve.gff3'

## data description:
1) 'eve.location.fa': sequences in this FASTA file are from Dmel's chromosome 2R between 5816824 to 5918300. This region harbors gene 'eve' as well as its 50kb flanking sequences. This sequence will be used as a query sequence to predict CRMs. Tandem Repeats in this sequence has already been masked (detailed info to mask repeat can be found in README.txt under '/SCRMshaw/' directory).

2) 'trainingSet.lst': a list of CRM training sets. User can include multiple training sets in the list.

3) 'mapping0.ap': a directory contains training data: CRM sequences -> 'crms.fasta', background sequences -> 'neg.fasta'. Both of them are in FASTA format. Please note that both of crm and neg sequences should be from the same genome. You can use script '/SCRMshaw/code/randomWithSameGC.pl' to generate background sequences if you don't have one. To mask Tandem Repeats in both of crm and neg sequences is highly recommended (detailed info to mask repeat can be found in README.txt under '/SCRMshaw/' directory).

4) 'eve.gff3': a gff3 file of gene eve from Dmel. Ideally this file should have all of the genes within the target sequences. A gff3 file is an optional input in the pipeline, and it will be used to generate 'exons' and 'genes' files. The exon file will is needed to mask exon regions in the query sequences and a gene file will be used to report top-scored windows. Instead of using a gff3 file, you can also choose to input 'exons' and 'genes' files directly. Please note that you will still be able to get top-scored windows without gff or gene file.
=====================================================================================

Command line to run SCRMshaw on the toy data (you can run it in your current directory (assuming you are in this example directory) by copy-and-paste this command line):

## Test run with gff file as input:
>perl ../code/scrm.pl --thitg 300 --gff ./test_data/eve.gff3 --imm --hexmcd --pac --genome ./test_data/eve.location.fa --traindirlst ./test_data/trainingSet.lst --outdir ./test_output_1

## Test run without gff file as input:
>perl ../code/scrm.pl --thitw 2000 --imm --hexmcd --pac --genome ./test_data/eve.location.fa --traindirlst ./test_data/trainingSet.lst --outdir ./test_output_2
=====================================================================================

Results after the first test run with gff file as input (the second example run will give you similar results except that in the 'hits' directory there are only top-scored windows (*.hits) without its nearest gene locus info):

___ ./test_output_1/
                    |
                    |___ fasta/
                    |          | 
                    |          |___ chr/
                    |          |               
                    |          |___ kmers/
                    |          |               
                    |          |___ windows/
                    |              
                    |___ gff/
                    |
                    |___ hits/
                    |         |              
                    |         |___ hexmcd/
                    |         |
                    |         |___ imm/
                    |         |              
                    |         |___ pac/
                    |
                    |___ scores/
                    |           |                
                    |           |___ hexmcd/
                    |           |                
                    |           |___ imm/        
                    |           |                
                    |           |___ pac/
                    |
                    |___ training/

## directory description:
1) chr/: This directory stores target sequence FASTA files, and each file will only have one sequence. Exon regions will be masked on the sequences if you provide the gff or exon file.

2) kmers/: This directory stores k-mer frequency files for each target sequence. For example, in this case './kmers/2R.fasta.1' is 1-mer and './kmers/2R.fasta.6' is 6-mer; and '/kmers/2R.fasta.rc.1' is 1-mer and '/kmers/2R.fasta.rc.6' is 6-mer of reverse complement sequences.

3) windows/: This directory stores target sequence FASTA files and each sequence in the files are in window-length (in this case: 500bp). Basically the pipeline will break target sequence down to window size for scoring.

4) gff/: This directory stores the 'exons' and 'genes' files. Both of them have the format (delimited by "\t"): chr_ID, start, end, string, gene_ID. This directory will be empty if you don't have gff, exons or genes file as input.

5) hits/: This directory stores top-scored windows from each scoring scheme.

6) hexmcd/: This directory stores top-scored windows from HexMCD scoring scheme. You can find sub-directories that have the same names as your training sets. Each directory has the top-scored windows (only when you privode gff or gene file. If you don't have gff or gene file, you will only get 'mapping0.ap.hits' in this case):

    ## 'mapping0.ap.hits' -> top scored windows in FASTA format:
    >chromosome_ID:CRM_locus    score   gene_ID:gene_symbol:distance_from_CRM:[within/upstream/downstream]
    AATGCCCAGAGAATGGGCAACAAGTAGCGGCGAATTAGCAATCCTATCATGCTTTTATGGCCGGCCAACTCTTGCC

    ## 'mapping0.ap.hits.ranked' -> top scored windows with local rank info in FASTA format:
    >chromosome_ID:CRM_locus    score   gene_ID:gene_symbol:distance_from_CRM:[within/upstream/downstream]:lr=rank
    AATGCCCAGAGAATGGGCAACAAGTAGCGGCGAATTAGCAATCCTATCATGCTTTTATGGCCGGCCAACTCTTGCC

    ## 'mapping0.ap.hits.rankLst' -> rank info for top scored windows in format:
    Gene_ID, CRM_ID, Score, Rank

7) imm/: This directory stores the top-scored windows from IMM scoring scheme. You will find the same files as described above.

8) pac/: This directory stores the top-scored windows from PAC scoring scheme. You will find the same files as described above.

9) scores/: This directory stores scored CRMs from each scoring scheme.

10) hexmcd/: This directory stores scored CRMs from HexMCD scoring scheme. You can find sub-directories that have the same names as your training sets. Each directory has files named after target sequence with format (delimited by " "): CRM_locus, Score

11) imm/: This directory stores the scores scored CRMs from IMM scoring scheme. You will find the same files as described above, except that there will be two extra files: 'crms.model' and 'neg.model', but you can ignore them.

12) pac/: This directory stores scored CRMs from PAC scoring scheme. You will find the same files as described above.

13) training/: This directory stores k-mer frequency files for each training set.
