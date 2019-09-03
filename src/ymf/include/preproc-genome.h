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
	int numChrom;	// number of chromosomes 
	char **chrom;	// the chromosomes ... may not be read
	double **P;	// the transition matrix
	double **Q;	// Q = inverse(P - I + 1p')
	double **QPPP;
	double **QPPQPP;
	double *p;	// the left eigenvector with unit eigenvalue (or the stationary distribution)
	long **PTmp;
	
public:
	Genome();
	void ReadChromAndUpdateP(char *filename);
	void ReadAllChromosomesAndComputeP(int argc, char **argv);
	void NormalizeP();
	void ComputeStationaryDistribution();
	void WriteToFile();
};



#endif
