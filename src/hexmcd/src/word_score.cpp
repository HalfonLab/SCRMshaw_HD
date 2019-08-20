#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "listop.h"

#define four2k 4096 //k is 6

void ReadFile(char *filename, int **counts, int *lens, int &gcount, int &nvals);

 //computes HexDiff word scores, given positive and negative training sets.
int main(int argc, char **argv) {
 if (argc < 2) { 
    printf("usage: %s <file1_pos> <file2_neg>  \n",argv[0]);
    exit(1);
  }
	char *file1 		= argv[1]; // has the k-word count (with rc) for the positive training set (one seq).
	char *file2 		= argv[2]; // has k-word count (with rc) for the negative training set (one seq).
//int wordlen 		= atoi(argv[3]); // this is k  . for hexdiff, k=6.
	int maxnumgenes = 1;
  int **counts1 	= new int*[maxnumgenes];
  int *lens1 			= new int[maxnumgenes];
  int genec1, nvals1;  //genec1=1, nvals = 4^k

  
  ReadFile(file1, counts1, lens1, genec1, nvals1);

  int **counts2 = new int*[maxnumgenes];
  int *lens2 = new int[maxnumgenes];
  int genec2, nvals2; // genec2 =1
  ReadFile(file2, counts2, lens2, genec2, nvals2);

  if (genec1 != 1 || genec2 !=1 || nvals1 != four2k || nvals2 != four2k ) {
    fprintf(stderr,"Error: inconsistent array dimensions \n");
    fprintf(stderr,"%d %d %d %d\n",genec1, genec2, nvals1, nvals2);
    exit(1);
  }
 
  //compute frequencies of k-words
  float f1[four2k]; // has freq of k-mers for the positive set
  float sum1 = 0;
  for (int l=0; l<four2k; l++) {
    f1[l] = counts1[0][l]+0.1; // pseudocounts added by SS
    sum1 += f1[l];
  }
  for (int l=0; l<four2k; l++) {
    f1[l] /= sum1;
  }
  float f2[four2k];// has k-mers freq for the neg set
  float sum2 = 0;
  for (int l=0; l<four2k ; l++) {
    f2[l] = counts2[0][l]+0.1; // pseudocounts added by SS
    sum2 += f2[l];
  }
  for (int l=0; l< four2k ; l++) {
    f2[l] /= sum2;
  }
  
  // compute k-word scores
  for (int l=0; l<four2k; l++) {
    double val = (f1[l])/(f2[l]) ; 
    //      int windowindex = j * windowshift;
    printf("%g %d\n",val, l);
  }
}


void ReadFile(char *filename, int **counts, int *lens, int &gcount, int &nvals)
{
  FILE *fp = fopen(filename,"r");
  fprintf(stderr,"Reading %s\n",filename);
  char line[500000];
  int genecount = 0;
  int numvals = -1;
  while (fgets(line, 500000, fp)) {
      std::list<char *> *toks = Split("\t",line);
    for (int i=0; i<2; i++) Shift(toks);
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
    fprintf(stderr,"Read gene %d with %d vals\n",genecount,numvals);
    delete toks;
  }
  gcount = genecount;
  nvals  = numvals;
  return;
}
