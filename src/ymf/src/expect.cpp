#include "genome.h"
#include "stats.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

double ComputeSimpleProbability(char *W, Genome &genome,int numN)
{
	// ignore numN

	int len = strlen(W);
	if (W[0] == 'N' || W[1] == 'N' || W[2] == 'N') {
		printf("Too many N's in the oligo ... Please give a larger lenOligo or smaller numN\n");
		exit(0);
	}

	char *dummy;
	char tmp[4];
	memcpy(tmp,W,3); tmp[3] = 0;
	long p0 = strtol(tmp,&dummy,4);
	double prob = genome.p[p0];

	int p = p0;
	int ns = 0;
	int lastTriplet = 0;
	int spacersSeen = 0;

	for (int j=3; j<len; j++) {
		if (W[j] == 'N') {
			if (W[j-1] != 'N') {
				lastTriplet = p;
				ns = 1;
				spacersSeen = 1;
			}
			else ns++;
		}
		long nT = 0;
		if (W[j] != 'N') {
			if (spacersSeen==1) {
				if (W[j-1]!='N' && W[j-2]!='N') {
					spacersSeen = 0;
					memcpy(tmp,&(W[j-2]),3);
					nT = strtol(tmp,&dummy,4);
					prob *= genome.P_numNplus1[ns][lastTriplet][nT];
				}
				else {
					ns++;
				}
			}
			else {
				nT = ((p%16)*4) + (W[j] - '0');
				prob *= genome.P[p][nT];
			}
		}
		p = nT;
	}
	return prob;
}

double ComputeProbability(int numWords, struct AllPatterns *W, Genome &genome,int numN)
{
	// n is the sequence length
	// takes an array of numWords strings (of same length) and computes 
	// E n(W)
	
	if (numWords < 1) {
		printf("Too few members of set W ... Aborting\n");
		exit(0);
	}
	int len = strlen(W[0].pattern);

	double probability = 0;

	for (int i=0; i<numWords; i++) {	// for each word do
		double prob = ComputeSimpleProbability(W[i].pattern,genome,numN);
		probability += prob;
	}

	return probability;
}
