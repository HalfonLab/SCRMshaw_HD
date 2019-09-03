/* Copyright (c) 2000 by Martin Tompa and Saurabh Sinha.
 * All rights reserved.  Redistribution is not permitted without the 
 * express written permission of the authors.
 * The program YMF implements an algorithm for identifying
 * likely transcription factor binding sites in yeast, described in the
 * following paper:
 "A Statistical Method for Finding Transcription Factor Binding Sites"
 by Saurabh Sinha and Martin Tompa,
 Eighth International Conference on Intelligent Systems for
 Molecular Biology, San Diego, USA, August 2000, 344-354.
 */

/********************************
 * YMF Version 2.0 update
 * ANSI scoping problem with for loops fixed
 * *******************************/

/********************************
 * YMF Version 2.1 updated
 * Up to 27 spacers allowed now
 * *******************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <list>
#include "preproc-genome.h"
#include "matrix.hpp"


Genome::Genome()
{
  numChrom = 0;
  P = new double *[64];
  Q = new double *[64];
  QPPP = new double *[64];
  QPPQPP = new double *[64];
  PTmp = new long *[64];
  p = new double[64];

  int i;
  for (i=0; i < 64; i++) {
    P[i] = new double[64];
    Q[i] = new double[64];
    QPPP[i] = new double[64];
    QPPQPP[i] = new double[64];
    PTmp[i] = new long[64];
  }


  for (i=0; i < 64; i++) {
    for (int j = 0; j < 64; j++) {
      P[i][j] = 0;
      Q[i][j] = 0;
      PTmp[i][j] = 0;
      QPPP[i][j] = 0;
      QPPQPP[i][j] = 0;
    }
    p[i] = 0;
  }
}	

void Genome::ReadChromAndUpdateP(char *fname)
{
  printf("Reading DNA file %s ... \n", fname);
  char *chromTmp = new char[CHROM_SIZE];
  char *filename = new char[128];
	
  strcpy(filename,GENOME_DIR);
  strcat(filename,fname);
  FILE *fpG = fopen(filename,"r");
  if (fpG == NULL) {
    printf("Error opening genome data file %s ... Aborting\n",filename);
    exit(0);
  }

  int ch = 0;
  while ((ch = fgetc(fpG)) != EOF && ch != '>');	// skip to the first '>'
  if (ch == EOF) return;
  while ((ch = fgetc(fpG)) != EOF && ch != '\n'); // skip the line
  if (ch == EOF) return;

  int i;	
  int ptr = 0;
  while (1) {
    ch = fgetc(fpG);	
    if (ch == '>' || ch == EOF) {	// end of a region ... process it 
      if (ch == '>')	
	while ((ch = fgetc(fpG)) != EOF && ch != '\n');

      chromTmp[ptr] = 0;
      int numNuc = ptr;
	
      if (ch=='N') ptr = -1;
      if (numNuc < 4) continue;
		
      char *curStr = new char[4];
      memcpy(curStr, chromTmp, 3);
      curStr[3] = 0;
      char *dummy;
		
      long curNum = strtol(curStr, &dummy, 4);
      delete [] curStr;
				
      long modulus = 16;
		
      for (i=3; i<numNuc; i++) {
	long nextNum = chromTmp[i] - '0';
	nextNum = (curNum%modulus)*4 + nextNum;
	PTmp[curNum][nextNum]++;
	curNum = nextNum;
      }
      ptr = 0;
      if (ch == EOF) break;
      else continue;
    }
				
    if (ch != '\n' && ch != '\r') {
      switch(ch) {
      case 'a':
      case 'A': ch = '0'; break;
      case 'c':
      case 'C': ch = '1'; break;
      case 'g':
      case 'G': ch = '2'; break;
      case 't':
      case 'T': ch = '3'; break;
      case 'y':
      case 'Y': 
      case 'r':
      case 'R':
      case 's':
      case 'S': 
      case 'w': 
      case 'W':
      case 'm':
      case 'M': 
      case 'k':
      case 'K': 
      case 'h':
      case 'H':
      case 'b':
      case 'B':
      case 'v':
      case 'V':
      case 'd':
      case 'D':
      case '-':
      case 'N': 
      case 'n': ch='N'; break;
      default : printf("Invalid character %c encountered at position %d... aborting\n", ch,ptr);
	printf("%s\n",chromTmp);
	exit(0);
      }
      if (ch!='N')
	chromTmp[ptr++] = ch;
    }
  }

  delete [] chromTmp;
  printf("Done reading DNA file %s ... \n", fname);
  fclose(fpG);

}

void Genome::NormalizeP()
{
  int i;
  for (i=0; i < 64; i++) {
    int numDN = 0;
    int j;
    for (j=0; j < 64; j++)
      numDN += PTmp[i][j];
    if (numDN == 0) {
      printf("Some row in P has sum = 0 ... Aborting\n");
      exit(0);
    }
    for (j=0; j < 64; j++)
      P[i][j] = (double) ((double)PTmp[i][j])/((double)numDN);
  }
}

FILE *fpT;

void Genome::ReadAllChromosomesAndComputeP(int argc, char **argv)
{
  char *dummy;
  if (argc < 3) {
    printf("usage: %s <numRegions> <file1> ...\n",argv[0]);
    printf("numRegions: maximum number of regions in the input files\n");
    printf("file1, file2, ... : the input files, containing the background sequences\n");
    exit(0);
  }
  numChrom = strtol(argv[1],&dummy, 10);
  int numFiles = argc - 2;

  for (int i=0; i<numFiles; i++) {
    ReadChromAndUpdateP(argv[2+i]);
  }
	
  NormalizeP();

  fpT = fopen("table.3","w");	
  for (int j=0; j < 64; j++) {
    for (int k=0; k < 64; k++)
      fprintf(fpT,"%f ",P[j][k]);
    fprintf(fpT,"\n");
  }

}

void Genome::ComputeStationaryDistribution()
{
  matrix In(64,64), I(64,64), Out(64,64), P_I(64,64), pp(64,64), P_I_pp(64,64);	
  double p[64];
  int i,j,k;

  for(i=0; i<64; i++)
    for (j=0; j<64; j++) {
      In(i+1,j+1) = P[i][j];
      Out(i+1,j+1) = P[i][j];
      if (i==j) 
	I(i+1,j+1) = 1;
      else 
	I(i+1,j+1) = 0;
    }

  for (i=0; i<6; i++)
    Out = Out*Out;

  for (i=0; i<64; i++) {
    p[i] = Out(1,i+1);
    fprintf(fpT,"%f ",p[i]);
  }
  fprintf(fpT,"\n");

  // also compute the Q


  P_I = In - I;
  for (i=0; i<64; i++)
    for (j=0; j<64; j++)
      pp(i+1,j+1) = p[j];

  P_I_pp = P_I + pp;
  Out = P_I_pp.inv();
  for (i=0; i<64; i++)
    for (j=0; j<64; j++)
      Q[i][j] = Out(i+1,j+1);

  for (j=0; j < 64; j++) {
    for (k=0; k < 64; k++)
      fprintf(fpT,"%f ",Q[j][k]);
    fprintf(fpT,"\n");
  }

  matrix QPPPt(64,64), QPPQPPt(64,64);
  QPPPt = ((Out*In)*In)*In;
  QPPQPPt = (((Out*In)*In)*Out)*(In*In);

  for (i=0; i<64; i++)
    for (j=0; j<64; j++)
      QPPP[i][j] = QPPPt(i+1,j+1);

  for (j=0; j < 64; j++) {
    for (k=0; k < 64; k++)
      fprintf(fpT,"%f ",QPPPt(j+1,k+1));
    fprintf(fpT,"\n");
  }
	
  for (i=0; i<64; i++)
    for (j=0; j<64; j++)
      QPPQPP[i][j] = QPPQPPt(i+1,j+1);

  for (j=0; j < 64; j++) {
    for (k=0; k < 64; k++)
      fprintf(fpT,"%f ",QPPQPPt(j+1,k+1));
    fprintf(fpT,"\n");
  }

  FILE * fpP = fopen("powers.3","w"); // record powers from 2 to 30
  FILE * fpGB = fopen("powersGeneralized.3.bin","wb"); // record generalized powers from 1 to 16
  matrix PT(64,64);
  matrix GPT(729,729);
  PT = In;

  for (int powers = 2; powers<31; powers++) {
    // also compute the generalized  powers table
    for (j=0; j<729; j++) {
      for (k=0; k<729; k++) {
	GPT(j+1,k+1) = 0;
      }
    }
    for (j=0; j<64; j++) {
      for (k=0; k<64; k++) {
	int indexj = j%4 + 9*((j/4)%4) + 81*((j/16)%4);	
	int indexk = k%4 + 9*((k/4)%4) + 81*((k/16)%4);	
	GPT(indexj+1,indexk+1) = PT(j+1,k+1);
      }
    }
		
    // create all columns of size 64 first
    for (k=0; k < 729; k++) {
      int base[3] = {0,0,0};	
      base[0] =  k%9;
      base[1] = (k/9)%9;
      base[2] = (k/81)%9;
      std::list<int> L[3];
      int l;
      for (l=0; l<3; l++) {
	switch (base[l]) {
	case 0:	L[l].push_back(0);break; // A
	case 1: L[l].push_back(1);break; // C		
	case 2: L[l].push_back(2);break; // G		
	case 3: L[l].push_back(3);break; // T		
	case 4: L[l].push_back(1);	L[l].push_back(2); break;	// S
	case 5: L[l].push_back(0);	L[l].push_back(3); break;	// W
	case 6: L[l].push_back(0);	L[l].push_back(2); break;	// R
	case 7: L[l].push_back(1);	L[l].push_back(3); break;	// Y
	case 8: L[l].push_back(1);	L[l].push_back(0); L[l].push_back(2); L[l].push_back(3); break;	// N
	}
      }
      for (std::list<int>::iterator iter1=L[0].begin(); iter1 != L[0].end(); iter1++) {
	int base0 = *iter1;	
	for (std::list<int>::iterator iter2=L[1].begin(); iter2 != L[1].end(); iter2++) {
	  int base1 = *iter2;
	  for (std::list<int>::iterator iter3=L[2].begin(); iter3 != L[2].end(); iter3++) {
	    int base2 = *iter3;
	    int columnNo = base0 + 9*base1 + 81*base2;
	    // copy (add) that column from PT onto GPT
	    for (j=0; j<729; j++) {
	      if (k!=columnNo) 	
		GPT(j+1,k+1) += GPT(j+1,columnNo+1);
	    }
	  }
	}
      }
    }

    // now create all the rows
    for (j=0; j < 729; j++) {
      int base[3] = {0,0,0};	
      base[0] =  j%9;
      base[1] = (j/9)%9;
      base[2] = (j/81)%9;
      std::list<int> L[3];
      int l;
      for (l=0; l<3; l++) {
	switch (base[l]) {
	case 0:	L[l].push_back(0);break; // A
	case 1: L[l].push_back(1);break; // C		
	case 2: L[l].push_back(2);break; // G		
	case 3: L[l].push_back(3);break; // T		
	case 4: L[l].push_back(1);	L[l].push_back(2); break;	// S
	case 5: L[l].push_back(0);	L[l].push_back(3); break;	// W
	case 6: L[l].push_back(0);	L[l].push_back(2); break;	// R
	case 7: L[l].push_back(1);	L[l].push_back(3); break;	// Y
	case 8: L[l].push_back(1);	L[l].push_back(0); L[l].push_back(2); L[l].push_back(3); break;	// N
	}
      }
      float sumOfWeights = 0;
      for (std::list<int>::iterator iter1=L[0].begin(); iter1 != L[0].end(); iter1++) {
	int base0 = *iter1;	
	for (std::list<int>::iterator iter2=L[1].begin(); iter2 != L[1].end(); iter2++) {
	  int base1 = *iter2;
	  for (std::list<int>::iterator iter3=L[2].begin(); iter3 != L[2].end(); iter3++) {
	    int base2 = *iter3;
	    int rowNo = base0 + 9*base1 + 81*base2;
	    int rowNo4 = base0 + 4*base1 + 16*base2;
	    sumOfWeights += p[rowNo4];
	    // copy (add) that row from GPT onto GPT
	    for (k=0; k<729; k++) {
	      if (j!=rowNo) 	
		GPT(j+1,k+1) += GPT(rowNo+1,k+1)*p[rowNo4];
	    }
	  }
	}
      }
      if (base[0] < 4 && base[1] < 4 && base[2] < 4) continue;
      for (k=0; k<729; k++)
	GPT(j+1,k+1) = GPT(j+1,k+1)/sumOfWeights;
    }

    for (j=0; j < 729; j++) {
      for (k=0; k < 729; k++) {
	float val = GPT(j+1,k+1);	
	fwrite(&val,sizeof(float),1,fpGB);	
      }
    }
		
    // and compute the powers of course
    PT = In*PT;
    for (j=0; j < 64; j++) {
      for (k=0; k < 64; k++) {
	fprintf(fpP,"%f ",PT(j+1,k+1));
      }
      fprintf(fpP,"\n");
    }
  }       
   
  fclose(fpP);
  fclose(fpGB);
}










