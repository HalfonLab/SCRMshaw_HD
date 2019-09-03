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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <list>
#include "genome.h"
#include "matrix.hpp"
#include <math.h>

int Genome::Order()
{
  return morder;
}

Genome::Genome(int mord)
{
  valid = false;
  morder = mord;
  numChrom = 0;

  int expmorder = (int)pow(4,mord);
  
  P = new double*[expmorder];
  Q = new double*[expmorder];
  QP = new double*[expmorder];
  QPPQ = new double*[expmorder];
  PTmp = new long*[expmorder];
  p = new double[expmorder];

  int i;
  for (i=0; i < expmorder; i++) {
    P[i] = new double[expmorder];
    Q[i] = new double[expmorder];
    QP[i] = new double[expmorder];
    QPPQ[i] = new double[expmorder];
    PTmp[i] = new long[expmorder];
  }


  for (i=0; i < expmorder; i++) {
    for (int j = 0; j < expmorder; j++) {
      P[i][j] = 0;
      Q[i][j] = 0;
      PTmp[i][j] = 0;
      QP[i][j] = 0;
      QPPQ[i][j] = 0;
    }
    p[i] = 0;
  }
}	

void Genome::ReadCounts(int *bkgcounts)
{
  int modulus = (int)pow(4,morder);  
  int expmorder = (int)pow(4,morder+1);
  for (int i=0; i<expmorder; i++) {
    int count = bkgcounts[i];
    int pre = int(i/4);
    int suf = i%modulus;
    PTmp[pre][suf] = count;
  }

  NormalizeP();
  ComputeStationaryDistribution();
}

void Genome::NormalizeP()
{
  int expmorder = (int)pow(4,morder);
  for (int i=0; i < expmorder; i++) {
    long numDN = 0;
    int j;
    for (j=0; j < expmorder; j++)
      numDN += PTmp[i][j];
    if (numDN == 0) {
      fprintf(stderr,"Warning: Some row in P has sum = 0\n");
      valid=false;
      return;
    }
    for (j=0; j < expmorder; j++)
      P[i][j] = (double) ((double)PTmp[i][j])/((double)numDN);
  }
  valid = true;
}

void Genome::ComputeStationaryDistribution()
{
  int emo = (int)pow(4, morder);
  matrix In(emo,emo), I(emo,emo), Out(emo,emo), P_I(emo,emo), pp(emo,emo), P_I_pp(emo,emo);
  int i,j,k;

  for(i=0; i<emo; i++)
    for (j=0; j<emo; j++) {
      In(i+1,j+1) = P[i][j];
      Out(i+1,j+1) = P[i][j];
      if (i==j) 
	I(i+1,j+1) = 1;
      else 
	I(i+1,j+1) = 0;
    }
  
  // just raise it to some  power (e.g., 2^10) to get the "stationary" distribution
  for (i=0; i<10; i++)
    Out = Out*Out;

  for (i=0; i<emo; i++) {
    p[i] = Out(1,i+1);
  }

  // also compute the Q, QP, QPPQ
  P_I = In - I;
  for (i=0; i<emo; i++)
    for (j=0; j<emo; j++)
      pp(i+1,j+1) = p[j];

  P_I_pp = P_I + pp;
  Out = P_I_pp.inv();
  for (i=0; i<emo; i++)
    for (j=0; j<emo; j++)
      Q[i][j] = Out(i+1,j+1);

  matrix QPt(emo,emo), QPPQt(emo,emo);
  QPt = Out*In; for (int i=0; i<morder-1; i++) QPt = QPt*In;
  QPPQt = (((Out*In)*In)*Out); for (int i=0; i<morder-1; i++) QPPQt = QPPQt*In;
  
  for (i=0; i<emo; i++)
    for (j=0; j<emo; j++)
      QP[i][j] = QPt(i+1,j+1);
	
  for (i=0; i<emo; i++)
    for (j=0; j<emo; j++)
      QPPQ[i][j] = QPPQt(i+1,j+1);

  return;
}










