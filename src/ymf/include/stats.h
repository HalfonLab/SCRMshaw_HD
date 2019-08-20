#ifndef _STAT_H_
#define _STAT_H_ 1

#include "genome.h"

struct AllPatterns {
	char *pattern;
	int count;	// if the pre is less than the suf
	int count2;	// if the pre is greater than the suf
	char preOrSuf; // = 0 if pre=suf, = -1 if pre < suf, = 1 if pre > suf
	long value;
};

struct results {
	char *pattern;
	int cnt;
	double exp;
	double var;
	double zsc;
};


double ComputePseudoPalindromeCorrection(char *curPattern,int sizeClosure,struct AllPatterns *closure,int *LenRegion,Genome &genome);
double ComputeOverlapExpectation(int numWords, struct AllPatterns *W, int *n, Genome &genome);
double ComputeAlmostVariance(int numWords, struct AllPatterns *W, int *n, Genome &genome, double expectation,int numN, int numRegions);
double ComputeProbability(int numWords,struct AllPatterns *W,Genome &genome,int);
double ComputeSimpleProbability(char *W, Genome &genome,int numN);
double ComputePalindromeCorrection(int numWords, struct AllPatterns *W, int *n, Genome &genome,int numN);
int IsPalindrome(char *W,char *&retPal);

#endif
