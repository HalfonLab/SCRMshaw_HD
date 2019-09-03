#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "listop.h"

#define four2k 4096 //k is 6
void readFile(char *filename, int **counts, int *lens, int &gcount, int &nvals);
void countMismatch(int *counts1, float *counts2, int nvals, int k);

//computes HexDiff scores (for all windows in one sequence).
int main(int argc, char **argv) {
	
	// Check arguments
	if (argc < 5) { 
    printf("usage: %s <pos_file> <neg_file> <window_file> <maxnumwindows> <winshift>\n",argv[0]);
    exit(1);
  }

	// Process arguments
 	char *file1 	= argv[1]; 				// File with k-mer counts with rc for the positive training set (one seq).
 	char *file2 	= argv[2]; 				// File with k-mer counts with rc for the negative training set (one seq).
 	char *file4 	= argv[3];  			// File with k-mer counts with rc for all windows of one seq.
 	int maxgenes 	= atoi(argv[4]); 	// this is (max) number of windows
	int winshift 	= atoi(argv[5]);

	// Read pos-kword count
	int numwin 		= 2;
  int **counts1 = new int*[numwin];
  int *lens1 		= new int[numwin];
  int genec1, nvals1; //genec1=1, nvals = 4^6
  readFile(file1, counts1, lens1, genec1, nvals1);
  
	// Count pos k-word mismatch
  float *counts11 = new float[four2k];
  countMismatch(counts1[0], counts11, four2k, 6);

	// Read neg k-word count
  int **counts2 = new int*[numwin];
  int *lens2 		= new int[numwin];
  int genec2, nvals2; // genec2 =1
  readFile(file2, counts2, lens2, genec2, nvals2);
	
	// Count neg k-word mismatch
  float *counts22=new float[four2k];
  countMismatch(counts2[0], counts22, four2k, 6);
 
  if (genec1 != 1 || genec2 !=1 || nvals1 != four2k || nvals2 != four2k ) {
    fprintf(stderr,"Error: inconsistent array dimensions \n");
    fprintf(stderr,"%d %d %d %d\n",genec1, genec2, nvals1, nvals2);
    exit(1);
  }

  // compute p1=pr(words |MC of positives) 
  float lambda=0;
  float f1[four2k];
  for (int l=0; l<four2k; l++) {
		float sum=0;
    long pre = l >> 2; // all but last char of l
    for (int i=0; i<4; i++) {
      long l2 = (pre << 2) + i; // pre followed by i
      sum += counts11[l2];
    }
		if(sum==0){
			f1[l]=0;
		} else{
			f1[l]=(counts11[l]+lambda)/(sum+4*lambda);
		}
  }
	// compute p2=pr(words |MC of negatives)
  float f2[four2k];
  for (int l=0; l<four2k ; l++) {
		float sum=0;
    long pre = l >> 2; // all but last char of l
    for (int i=0; i<4; i++) {
    	long l2 = (pre << 2) + i; // pre followed by i
      sum += counts22[l2];
    }
		if(sum==0){
			f2[l]=0;
		} else{
     	f2[l]=(counts22[l]+lambda)/(sum+4*lambda);
		}
  }

  // compute the word weight as log(p1/p2)
  float wordscore[four2k];
  for (int l=0; l<four2k; l++) {
    if(f1[l]==0 || f2[l]==0){
			wordscore[l]=0;
    } else{
	wordscore[l] = log((f1[l])/(f2[l]));
    }
  }

  // Load windows one by one, compute and output window scores:
  int *counts4 = new int[four2k]; // make it one dimension array
  int lens4;                       // make it a number 
  int genec4, nvals4; // genec4 has number of windows
  
  FILE *fp = fopen(file4,"r");
//  fprintf(stderr,"Reading %s\n",file4); // changed by Wei
  
	char line[500000];
  int genecount = 0;
  int numvals = -1;
  while (fgets(line, 500000, fp)) {
      std::list<char *> *toks = Split("\t",line);
    for (int i=0; i<2; i++) {
			Shift(toks);
		}
    lens4 = atoi(Shift(toks));    
    if (numvals < 0) {
			numvals = toks->size();
		} else if (numvals != toks->size()) {
			fprintf(stderr,"Error: rows of unequal sizes\n");
			exit(1);
    }
         
		// Score
    float val=0;
    for (int i=0; i < four2k; i++) {
      counts4[i] = atoi(Shift(toks));
      val += counts4[i] * wordscore[i];
    }
    
    int windowindex = genecount * winshift;
    double nsc=0;
    if(lens4) {
			nsc=val/lens4;
		}
 //   printf("%d %g %.5f %d\n",windowindex, val, nsc, lens4);  
    printf("%d %g\n",windowindex, val);
    genecount++;
    delete toks;
  }

  genec4 = genecount; // the number of windows .. 
  nvals4  = numvals;
//  fprintf(stderr, "Done with reading.\n"); // changed by Wei
  fclose(fp);
    
  delete[] lens1;
  delete[] lens2;
  delete[] counts11;
  delete[] counts22;
  
  //delete counts1,2,4
  for(int i=0; i<genec1; i++){
		delete[] counts1[i];
  }
  for(int i=0; i<genec2; i++){
  	delete[] counts2[i];
  }

  delete[] counts1;
  delete[] counts2;
  delete[] counts4;

}

void readFile(char *filename, int **counts, int *lens, int &gcount, int &nvals)
{
	FILE *fp = fopen(filename,"r");
//  fprintf(stderr,"Reading %s\n",filename); //changed by Wei

  char line[500000];
  int genecount = 0;
  int numvals 	= -1;
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
    //fprintf(stderr,"Read gene %d with %d vals\n",genecount,numvals);
    delete toks;
  }
  gcount = genecount;
  nvals  = numvals;
//  fprintf(stderr, "Done with reading.\n"); // changed by Wei
  fclose(fp);
  return;
}

// gets the original count array (counts1) and set mismatch count array(counts2) 
void countMismatch(int *counts1, float *counts2, int nvals, int k)
{
	double wt=0.25;
  for (int l=0; l<nvals; l++) {
    counts2[l]=counts1[l]; // the no mismatch term
    for (int j=1; j<=k; j++) { // enumerate all words mismatching at jth pos
      long pre, suf;
      if (j == 1) {
        suf = l & ((1 << (2*(k-1)))-1); // last k-1 chars of l
        long lfirst = l>>(2*(k-1)); // first char of l
        for (int i=0; i<4; i++) {
          if (i==lfirst) continue;
          long l2 = (i << (2*(k-1))) + suf; // i followed by suf
          counts2[l] += wt*counts1[l2];
        }
      }
      if (j == k) {
        pre = l >> 2; // all but last char of l
        long llast = l & ((1 << 2)-1); // last char of l
        for (int i=0; i<4; i++) {
          if (i==llast) continue;
          long l2 = (pre << 2) + i; // pre followed by i
          counts2[l] += wt*counts1[l2];
        }
      }
      if (j != 1 && j != k) { // the general case
        pre = l >> (2*(k-j+1)); // first j-1 chars of l
        suf = l & ((1 << (2*(k-j)))-1); // last (k-j) chars of l
        long lj = (l >> (2*(k-j))) & ((1 << 2)-1);
        for (int i=0; i<4; i++) {
          if (i==lj) continue;
          long l2 = (pre << (2*(k-j+1))) + suf; // pre followed by i
          counts2[l] += wt*counts1[l2];
        }
      }
    }
  }
}
