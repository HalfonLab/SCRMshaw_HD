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


#include "genome.h"
#include "stats.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>





double ComputeWordProb(long word, Genome *bkg, int k, int mo)
  // need to know mo, since that decides the prefix
  // need to know total length of word
{
  long prefix = word>>(2*(k-mo)); // equivalent to word/pow(4,k-mo); gets first mo chars of word
  long suffix = word & ((1<<(2*(k-mo)))-1) ; // bit-and with a number that is all ones in the last 2*(k-mo) bits, gets last (k-mo) chars of word
  return bkg->p[prefix]*ComputeWordProbGivenPrefix(prefix,suffix,bkg,k,mo);  
}

double ComputeWordProbGivenPrefix(long prefix, long suffix, Genome *bkg, int k, int mo)
  // calculate only if k-mo>=1; otherwise return 1
  // computes p_*(prefix suffix)
{
  if (k-mo<=0) return 1;
  
  double prob = 1;
  long pre = prefix; // assumed of length mo
  long jcopy = suffix; // assumed of length k-mo

  // need to iterate through successive positions of suffix
  for (int l=0; l<k-mo; l++) {
    long jcopyfirst = jcopy>>(2*(k-mo-1)); // get first char of jcopy, which is of length k-mo; equivalent to jcopy/pow(4,k-mo-1)
    long suf = ((pre & ((1<<(2*(mo-1)))-1))<<2) + jcopyfirst; // last (mo-1) chars of pre, cat jcopyfirst; equivalent to (pre%pow(4,mo-1))*4+jcopyfirst;

    prob *= bkg->P[pre][suf];
    
    // erase first char of jcopy, and preserve its length at k-mo by shifting one char to left; equivalent to (jcopy%pow(4,k-mo-1))*4;
    jcopy = (jcopy & ((1<<(2*(k-mo-1)))-1))<<2; 
    pre = suf;
  }
  return prob;
}

double ComputeWordProbIID(long word, double *bkg, int k)
{
  double prob = 1;
  for (int l=0; l<k; l++) { 
    int ch = word%4; // the last char in word
    prob *= bkg[ch];
    word >> 2; // remove the last char (2 bits) of word
  }
  return prob;
}

double PoissonDistribution(int x, double mean)
{
  double sum = 0;
  for(int j=0; j<=x; j++) {
    double term = 1;
    for (int k=1; k<=j; k++) {
      term *= mean/k;
    }
    sum += term;
  }
  // fprintf(stderr,"sum and poisson distrib are ");
  // fprintf(stderr,"%g\t%g\n", sum, sum/exp(mean));
  return sum/exp(mean); // this is sum_{j=0}^x (m^j/j!)*e^{-m} 
}

double PoissonMin(int c, double mean1, double mean2)
  // computes Pr(X >= c), where X = min(X1, X2) and X1 ~ Poisson(mean1), X2 ~ Poisson(mean2).
{
  if (c > 0) {
    return (1- PoissonDistribution(c-1, mean1))*(1- PoissonDistribution(c-1, mean2));
  }
  else {
    return 1;
  }
}

double Factorial(int n){
  if (n < 0)  return 0;
  double f = 1;
  for (int i=1; i <= n; i++){
    f *= i;
  }
  
  // fprintf(stderr,"j and j! are ");
  // fprintf(stderr,"%d\t%g\n",n, f);
  return f; 
}



