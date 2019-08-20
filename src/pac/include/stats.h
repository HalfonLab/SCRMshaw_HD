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

#ifndef _STAT_H_
#define _STAT_H_ 1

#include "genome.h"

double ComputeWordExpectation(long word, int slen, int k, Genome *bkg);
double ComputeWordProb(long word, Genome *bkg, int k, int mo);
double ComputeWordProbGivenPrefix(long prefix, long suffix, Genome *bkg, int k, int mo);
double ComputeWordProbIID(long word, double *bkg, int k);
double PoissonDistribution(int x, double mean);
double PoissonMin(int c, double mean1, double mean2);
double Factorial(int n);


#endif
