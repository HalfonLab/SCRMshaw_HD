#ifndef _GENOME_H_
#define _GENOME_H_ 1

#define CHROM_SIZE 1600000
#define GENOME_DIR ""

class Genome {
public:
	Genome();
	int numChrom;	// number of chromosomes 
	char **chrom;	// the chromosomes ... may not be read
	double P[64][64];	// the transition matrix
	double P_numNplus1[16][64][64]; // P^(numN+1)
	double Q[64][64];	// Q = inverse(P - I + 1p')
	double QPPP[64][64];
	double QPPQPP[64][64];
	double p[64];	// the left eigenvector with unit eigenvalue (or the stationary distribution)
	long PTmp[64][64];
	
};

#endif
