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

#ifndef _GENOME_H_
#define _GENOME_H_ 1

#define CHROM_SIZE 25600000
#define GENOME_DIR ""


class Genome {  
 public:
  bool valid;   // has this been properly trained
  int morder;   // markov order
  int numChrom;	// number of chromosomes 
  char **chrom;	// the chromosomes ... may not be read
  double **P;	// the transition matrix
  double **Q;	// Q = inverse(P - I + 1p')
  double **QP;  // QPP^(morder-1)
  double **QPPQ;// QPPQP^(morder-1)
  double *p;	// the left eigenvector with unit eigenvalue (or the stationary distribution)
  long **PTmp;
  
  Genome(int morder=3); // markov order is 3 by default
  void ReadCounts(int *bkgcounts); // create markov chain 
  int Order();
  void NormalizeP(); // normalize markov chain to probabilities (from counts)
  void ComputeStationaryDistribution();
};



#endif
