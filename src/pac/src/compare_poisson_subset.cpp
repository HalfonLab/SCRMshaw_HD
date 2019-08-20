#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "listop.h"
#include "genome.h"
#include "stats.h"
#include <string.h> 

/************************************
A count file has one line for each sequence.
The first two columns (tab separated) are names, the first one is important.
The third column is length of the sequence.
Fourth column onwards are counts of all 4^k words for some fixed k.

File1: All sequences are from first set (or species). Here we have only ONE sequence
File2: All sequences are from second set(or species). Here These are windows from a sequence in spc 2.
Bkgfile1: One for each sequence in File1; provides background
Bkgfile2: One; provides global background for the sequence of spc2.
*************************************/


void ReadFile(char *filename, int **counts, int *lens, int &gcount, int &nvals);
double Compare(int *counts1, int len1, Genome *bkg1, int *counts2, int len2, Genome *bkg2, int k, int mo, int *counts3);
double CompareIID(int *counts1, int len1, double *bkg1, int *counts2, int len2, double *bkg2, int k, int *counts3);

int main(int argc, char **argv)
{
  if (argc < 9) { 
    printf("usage: %s <file1> <file2> <bkgfile1> <bkgfile2> <maxnumgenes> <wordlength> <bkgwordlength> <windowshift> <file3>[-diag]\n",argv[0]);
    exit(1);
  }
  char *file1 = argv[1];
  char *file2 = argv[2];
  char *bkgfile1 = argv[3];
  char *bkgfile2 = argv[4];
  int  maxnumgenes = atoi(argv[5]);
  int  wordlen = atoi(argv[6]); // this is k
  int  markovlen = atoi(argv[7]); // this is markov order + 1
  int windowshift = atoi(argv[8]); 
  char *file3 = argv[9]; // subset of words
    bool diagonal_only = false;
  if (argc > 10 && !strcmp(argv[10],"-diag")) {
     diagonal_only = true;
  }

  int num = 1;
  int four2k = 4096; 

  int **counts1 = new int*[num];
  int *lens1 = new int[num];
  int genec1, nvals1;
  ReadFile(file1, counts1, lens1, genec1, nvals1);

  int **counts3 = new int*[num];
  int *lens3 = new int[num];
  int genec3, nvals3;
  ReadFile(file3, counts3, lens3, genec3, nvals3);

  if (genec1 != 1) {
    fprintf(stderr,"Error: genec1 not equal 1\n");
    exit(1);
  }
  if (nvals1 != pow(4,wordlen)) {
    fprintf(stderr,"Error: incorrect array dimensions for counts\n");
    exit(1);
  }

  if (genec3 != 1 || nvals3 != nvals1) {
    fprintf(stderr,"Error: incorrect array dimensions for counts in file3\n");
    exit(1);
  }
  
  int **bkgcounts1 = new int*[maxnumgenes];
  int *bkglens1 = new int[maxnumgenes]; // these wont be needed
  int bkggenec1, bkgnvals1;
  ReadFile(bkgfile1, bkgcounts1, bkglens1, bkggenec1, bkgnvals1);

  if (bkgnvals1 != pow(4,markovlen)) {
    fprintf(stderr,"Error: incorrect array dimensions for background counts\n");
    exit(1);
  }

  // compute the background model tables from bkgcounts*
  // Here genec1 == 1 so it is ok !
  Genome **bkg1 = new Genome*[genec1];  
  double **bkgarray1 = new double*[genec1];
  if (markovlen > 1) 
    for (int i=0; i<genec1; i++) {
      bkg1[i] = new Genome(markovlen-1);
      bkg1[i]->ReadCounts(bkgcounts1[i]);
    }
  else {
    for (int g=0; g<genec1; g++) {
      bkgarray1[g] = new double[4];
      double sum = 0;
      for (int i=0; i<4; i++) sum += bkgcounts1[g][i];
      //      if (sum <= 0) { fprintf(stderr,"Error: no background counts\n"); exit(1); }
      if (sum > 0){
	for (int i=0; i<4; i++) bkgarray1[g][i] = bkgcounts1[g][i]/sum;
      }
    }
  }

   //-----------------
   int **counts2 = new int*[num];
   int *lens2 = new int[num];
   int genec2, nvals2;

   int **bkgcounts2 = new int* [num];
   int *bkglens2 = new int[num]; // these wont be needed
   int bkggenec2, bkgnvals2;
 
   FILE *fp2 = fopen(file2, "r");
   FILE *bkgfp2 = fopen(bkgfile2, "r");

  // fprintf(stderr,"Reading %s\n",file2); // changed by Wei
 //  fprintf(stderr, "Reading %s\n", bkgfile2); changed by Wei

   char line2[500000];
   int genecount2 = 0;
   int numvals2 = -1;
   char bkgline2[500000];
   int bkggenecount2 = 0;
   int bkgnumvals2 = -1;

   counts2 = new int * [num];
   counts2[0] = new int [four2k];
   bkgcounts2 = new int * [num];
   bkgcounts2[0] = new int [four2k];

   Genome **bkg2 = new Genome*[num]; // 
   double **bkgarray2 = new double*[num];
   bkg2[0] = new Genome(markovlen-1);
   bkgarray2[0] = new double[4];

   int j=0; 
   while (fgets(line2, 500000, fp2)) {

        std::list <char *> *toks2 = Split("\t",line2);
        for (int i=0; i<2; i++) Shift(toks2);    
        lens2[0] = atoi(Shift(toks2));
        if (numvals2 < 0) numvals2 = toks2->size();    
        else {
        if (numvals2 != toks2->size()) {
	   fprintf(stderr,"Error: rows of unequal sizes\n");
	   exit(1);
        }
      }
      //counts2[0] = new int[numvals2];    
      for (int i=0; i<numvals2; i++) {
          counts2[0][i] = atoi(Shift(toks2));
      }
      genecount2++;
      //fprintf(stderr,"Read gene %d with %d vals\n",genecount2,numvals2);
      delete toks2;

      // --------------------------------
      fgets(bkgline2, 500000, bkgfp2);
      std::list<char *> *bkgtoks2 = Split("\t",bkgline2);
      for (int i=0; i<2; i++) Shift(bkgtoks2);    
      bkglens2[0] = atoi(Shift(bkgtoks2));
      if (bkgnumvals2 < 0) bkgnumvals2 = bkgtoks2->size();    
      else {
        if (bkgnumvals2 != bkgtoks2->size()) {
	   fprintf(stderr,"Error: rows of unequal sizes\n");
	   exit(1);
        }
      }
      //bkgcounts2[0] = new int[bkgnumvals2];    
      for (int i=0; i<bkgnumvals2; i++) {
          bkgcounts2[0][i] = atoi(Shift(bkgtoks2));
      }
      bkggenecount2++;
      //fprintf(stderr,"Read gene %d with %d vals\n",bkggenecount2,bkgnumvals2);
      delete bkgtoks2;
    
      // -----------------------------------   
      double sumbkgd; 
      Genome **bkg2 = new Genome*[num]; // 
      double **bkgarray2 = new double*[num];
      if (markovlen > 1) { 
        bkg2[0] = new Genome(markovlen-1);
        bkg2[0]->ReadCounts(bkgcounts2[0]);
      }
      else {
        bkgarray2[0] = new double[4];
        double sum = 0;
        for (int i=0; i<4; i++) sum += bkgcounts2[0][i];
        //  if (sum <= 0) { fprintf(stderr,"Error: no background counts\n"); exit(1); }
        sumbkgd = sum;
        if (sum > 0){
           for (int i=0; i<4; i++) bkgarray2[0][i] = bkgcounts2[0][i]/sum;
        }
      }
   
      int windowindex = j * windowshift;
      if (markovlen > 1) {
				if (bkg1[0]->valid && bkg2[0]->valid && lens1[0] >= 200 && lens2[0] >= 200) {
					double val = Compare(counts1[0], lens1[0], bkg1[0], counts2[0], lens2[0], bkg2[0], wordlen, markovlen-1, counts3[0]);
					printf("%d %.4f\n",windowindex, val);
				} else {
					printf("%d 0.0000\n",windowindex);
				}
    	} else {
  	 		if (lens1[0] >= 200 && lens2[0] >= 200 && sumbkgd > 0) {
  	     	double val = CompareIID(counts1[0], lens1[0], bkgarray1[0], counts2[0], lens2[0], bkgarray2[0], wordlen, counts3[0]);
  	     	printf("%d %.4f\n",windowindex, val);	  
        	} else {
  	    	printf("%d 0.0000\n",windowindex);
					}	
      }
      
      j++; 
  }
  genec2 = genecount2;
  nvals2  = numvals2;
  
  bkggenec2 = bkggenecount2;
  bkgnvals2  = bkgnumvals2;

  fclose(fp2);
  fclose(bkgfp2);

//  fprintf(stderr,"Done with reading %s\n",file2);  // changed by Wei
//  fprintf(stderr, "Done with reading %s\n", bkgfile2);  // changed by Wei


  delete [] counts1[0];
  delete [] counts1; 
  delete [] lens1;  
  delete [] counts2[0];
  delete [] counts2;
  delete [] lens2; 
  delete [] bkgcounts1[0];
  delete [] bkgcounts1;
  delete [] bkglens1; 
  delete [] bkgcounts2[0];
  delete [] bkgcounts2;
  delete [] bkglens2; 
  delete [] bkg2[0];
  delete [] bkg2;
  delete [] bkgarray2[0];
  delete [] bkgarray2; 
  delete [] counts3[0];
  delete [] counts3;
  delete [] lens3; 

}

double Compare(int *counts1, int len1, Genome *bkg1, int *counts2, int len2, Genome *bkg2, int k, int mo, int *counts3)
  // counts1[] and counts2 have the counts of all 4^k k-mers
  // counts3[l] is 1 if l is a word in our subset and 0 otherwise.
  // mo is the markov order for which bkg1 and bkg2 have the Markov chain and associated matrices
  // compare using Poisson statistic: van Helden's ADDITIVE model if _ADDITIVE is defined else product model
{
  int nvals = (int)pow(4, k);

  // Compute score.
  // score = (sum__{i=0}^{nvals - 1}(1-Pr( X >= min(counts1[i], counts2[i]))))/nvals,
  // where X = min(X1, X2), with X1 is known and X2 ~ POisson(mean2).

  int nwords = 0; // number of words in our subset.
  double sum = 0;
  for (int l=0; l<nvals; l++) {
    nwords += counts3[l];
    //   double mean1 = (len1-k+1) * ComputeWordProb(l, bkg1, k, mo);
    double mean2 = (len2-k+1) * ComputeWordProb(l, bkg2, k, mo);
    int cmin =  std::min(counts1[l] , counts2[l]);

    //    double pm = PoissonMin(cmin, mean1, mean2);
    //    double s = (1 - pm);

    double s =  PoissonDistribution(cmin - 1, mean2);
    double sprod = 1-PoissonDistribution(cmin, mean2);

if (s < 0) { 
      fprintf(stderr,"Warning: (1-poissonmin)/nvals for cmin=%d mean2=%g is ",cmin,mean2);
      fprintf(stderr,"%g\n",s );
      s = 0;
    }
    //fprintf(stderr,"HERE HERE\n");
    //if (counts3[l] > 0) {
    // fprintf(stderr,"%.3f (%d, %.2f, %.2f)\t",s,cmin,mean1,mean2);
    //}
#ifdef _ADDITIVE
    sum += s*counts3[l]; 
#else
    sum += log(sprod)*counts3[l];
#endif
  } 
  // fprintf(stderr,"\n");

#ifdef _ADDITIVE
  double score = sum/nwords;
#else
  double a = double(sum)/double(nwords);
  double score = 1-exp(a);
  // fprintf(stderr,"a is %g\n",a);
#endif
//  fprintf(stderr,"score is "); // not output scores, changed by Wei
//  fprintf(stderr,"%g\n",score);
  
  return score;
}

double CompareIID(int *counts1, int len1, double *bkg1, int *counts2, int len2, double *bkg2, int k, int *counts3)
{
  // counts*[] has the counts of all 4^k k-mers
  // bkg1 and bkg2 have the bkg counts
  // compare using Poisson statistic: van Helden's ADDITIVE model if _ADDITIVE is defined else product model
  int nvals = (int)pow(4, k);

  // Compute score.
  // score = (sum__{i=0}^{nvals - 1}(1-Pr( X >= min(counts1[i], counts2[i]))))/nvals,
  // where X = min(X1, X2), with X1 known and X2 ~ POisson(mean2).

  int nwords =0;
  double sum = 0;
  for (int l=0; l<nvals; l++) {
    nwords += counts3[l];
    double mean2 = (len2-k+1) * ComputeWordProbIID(l, bkg2, k);
    int cmin =  std::min(counts1[l] , counts2[l]);
    double s = PoissonDistribution(cmin - 1, mean2);
    double sprod = 1-PoissonDistribution(cmin, mean2);
    if (s < 0) { 
      fprintf(stderr,"Warning: (1-poissonmin)/nvals for cmin=%d mean2=%g is ",cmin,mean2);
      fprintf(stderr,"%g\n",s );
      s = 0;
    }
#ifdef _ADDITIVE
    sum += s*counts3[l]; 
#else
    sum += log(sprod)*counts3[l];
#endif
  } 
//  fprintf(stderr,"\n"); // not output scores, changed by Wei
#ifdef _ADDITIVE
  double score = sum/nwords;
#else
  double a = double(sum)/double(nwords);
  double score = 1-exp(a);
  // fprintf(stderr,"a is %g\n",a);
#endif
//  fprintf(stderr,"score is "); // not output scores, changed by Wei
//  fprintf(stderr,"%g\n",score);
  
  return score;
}


void ReadFile(char *filename, int **counts, int *lens, int &gcount, int &nvals)
{
  FILE *fp = fopen(filename,"r");
//  fprintf(stderr,"Reading %s\n",filename); // changed by Wei
  char line[500000];
  int genecount = 0;
  int numvals = -1;
  while (fgets(line, 500000, fp)) {
      std::list<char *> *toks = Split("\t",line);
    for (int i=0; i<2; i++) {
			Shift(toks);    
		}
    lens[genecount] = atoi(Shift(toks));
    if (numvals < 0) { 
			numvals = toks->size();    
		} else if (numvals != toks->size()) {
			fprintf(stderr,"Error: rows of unequal sizes\n");
			exit(1);
    }
    counts[genecount] = new int[numvals];    
    for (int i=0; i<numvals; i++) {
      counts[genecount][i] = atoi(Shift(toks));
    }
    genecount++;
    delete toks;
  }
  gcount = genecount;
  nvals  = numvals;
//  fprintf(stderr, "Done with reading %s\n", filename); //changed by Wei
  return;
}
