#include "genome.h"
#include "stats.h"
#include "math.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

extern int *totalLenRegion;
extern int *totalSquareLenRegion;

double ComputeOverlapExpectation(int numWords, struct AllPatterns *W, int *n, Genome &genome)
{
	int len = strlen(W[0].pattern);
	int numN = 0;
	for (int i=0; i<len; i++) {
		if (W[0].pattern[i]=='N')
			numN++;
	}
	
	char *tmp = new char[2*len];

	double overlapExpectation = 0;
	for (int i=0; i<numWords; i++) {
		for (int j=0; j<numWords; j++) {
			for (int k=1; k<len; k++) {	// WORRY !!!!
				// is i,j an overlap at shift k ? 
				int ptr = 0;
				int l;
				for (l=0; l<k; l++,ptr++)
					tmp[ptr] = W[i].pattern[l];
				for (l=0; l<len; l++)
					tmp[ptr++] = W[j].pattern[l];
				tmp[ptr] = 0;					
				int isO = 1;
				for (int m=k; m<len; m++) {
					if (W[i].pattern[m] == 'N' && tmp[m] == 'N') {
						continue;
					}
					if (W[i].pattern[m] != 'N' && tmp[m] == 'N') {
						tmp[m] = W[i].pattern[m];
						continue;
					}
					if (W[i].pattern[m] == 'N' && tmp[m] != 'N') continue;
					if (W[i].pattern[m] == tmp[m]) continue;
					isO = 0;
					break;
				}

				if (isO)
					overlapExpectation += ComputeSimpleProbability(tmp, genome, numN); // numN is vestigial
				
			}
		}
	}
	return overlapExpectation*totalLenRegion[numN];	
}

char flip(char c) 
{
	if (c=='0') return '3';
	if (c=='1') return '2';
	if (c=='2') return '1';
	if (c=='3') return '0';
	if (c=='N') return 'N';
	return c;
}


int IsPalindrome(char *W,char *&retPal)
// checks if this is a palindrome or a semi-palindrome. if its the latter, 
// it returns the actual palindromic instantiation
{

	int len = strlen(W);
	char *pal = new char[len+1];
	for (int c=0; c<len+1; c++)
		pal[c] = 0;
	int isP = 1;
	char oChar = 0;
	for (int i=0; i<len; i++) {
		switch (W[i]) {
		case '0':
			oChar = W[len-1-i];
			if (oChar != '3' && oChar != 'N') 
				isP = 0;
			else 
				pal[i] = W[i];
			break;
		case '3':
			oChar = W[len-1-i];
			if (oChar != '0' && oChar != 'N') 
				isP = 0;
			else 
				pal[i] = W[i];
			break;
		case '1': 
			oChar = W[len-1-i];
			if (oChar != '2' && oChar != 'N') 
				isP = 0;
			else 
				pal[i] = W[i];
			break;
		case '2': 
			oChar = W[len-1-i];
			if (oChar != '1' && oChar != 'N') 
				isP = 0;
			else 
				pal[i] = W[i];
			break;
		case 'N': 
			oChar = W[len-1-i];
			pal[i] = flip(oChar);
			break;
		default : printf("Invalid character encountered .... Aborting\n");
			exit(0);
		}
		if (isP == 0) break;
	}
	retPal = pal;
	return isP;
}

int countN(char *pal)
{
	int len = strlen(pal);
	int count = 0;
	for (int i=0; i<len; i++)
		if (pal[i]=='N') 
			count++;
	return count;
}

double ComputePalindromeCorrection(int numWords, struct AllPatterns *W, int *LenRegion, Genome &genome,int numN)
{
	double correction = 0;
	for (int i=0; i<numWords; i++) {
		char *pal;
		if (IsPalindrome(W[i].pattern,pal)) {
			int newNumN = countN(pal);
			correction += totalLenRegion[newNumN]*ComputeSimpleProbability(pal,genome,newNumN);
		}
		delete [] pal;
	}
	return correction;
}


double ComputeAlmostVariance(int numWords, struct AllPatterns *W, int *LenRegion, Genome &genome, double probability,int numN, int numRegions)
// the expectation has already been computed and is going to be used again, hence
// is passed as an argument
{

	if (numWords < 1) {
		printf("Too few members of set W ... Aborting\n");
		exit(0);
	}
	int len = strlen(W[0].pattern);
	

	double *pStar = new double[numWords];			// should pre-compute these 

	for (int i=0; i<numWords; i++) {
		if (W[i].pattern[0] == 'N') {
			printf("Too many Ns in the oligo ... change numN\n");
			exit(0);
		}
	
		char *dummy;
		char tmp[4];
		memcpy(tmp,W[i].pattern,3); tmp[3] = 0;
		long p0 = strtol(tmp,&dummy,4);
		double prob = 1;

		int p = p0;
		int ns = 0;
		int lastTriplet = 0;
		int spacersSeen = 0;
		
		for (int j=3; j<len; j++) {
           	if (W[i].pattern[j]== 'N') {
                if (W[i].pattern[j-1] != 'N') {
                    lastTriplet = p;
                    ns = 1;
                    spacersSeen = 1;
                }
                else ns++;
            }
            long nT = 0;
            if (W[i].pattern[j] != 'N') {
                if (spacersSeen==1) {
                    if (W[i].pattern[j-1]!='N' && W[i].pattern[j-2]!='N') {
                        spacersSeen = 0;
                        memcpy(tmp,&(W[i].pattern[j-2]),3);
                        nT = strtol(tmp,&dummy,4);
                        prob *= genome.P_numNplus1[ns][lastTriplet][nT];
                    }
                    else {
                        ns++;
                    }
                }
                else {
                    nT = ((p%16)*4) + (W[i].pattern[j] - '0');
                    prob *= genome.P[p][nT];
                }
            }
            p = nT;
		}
		pStar[i] = prob;
	}

	
	double variance = probability*totalLenRegion[numN];
	variance += (ComputePalindromeCorrection(numWords, W, LenRegion, genome,numN));
 
	// skipping overlap term computation, since this is the estimate

	double A = 0;

	for (int i=0; i<numRegions; i++) {
	double q = LenRegion[i] - 2*len + 2;

	char *dummy; char tmp[4];

	for (int i1=0; i1<numWords; i1++) {
		for(int i2=0; i2<numWords; i2++) {

            memcpy(tmp,&(W[i1].pattern[len-3]),3); tmp[3] = 0;
            long i1last = strtol(tmp, &dummy, 4);
            memcpy(tmp,W[i2].pattern,3); tmp[3] = 0;
            long i2first = strtol(tmp, &dummy, 4);
            memcpy(tmp,W[i1].pattern,3); tmp[3] = 0;
            long i1first = strtol(tmp, &dummy, 4);

			double tmpA = pStar[i1]*pStar[i2];
			double eQPPPe = genome.QPPP[i1last][i2first];
			double eQPPQPPe = genome.QPPQPP[i1last][i2first];
			tmpA *= genome.p[i1first] * ( (q*(q+1)/2)*genome.p[i2first] - (q-1)*eQPPPe - eQPPQPPe );
			A += tmpA;
		}
	}
	}

	variance += 2*A;
	variance -= (totalSquareLenRegion[numN]*probability*probability);
	return variance;
}

double ComputePseudoPalindromeCorrection(char *curPattern,int sizeClosure,struct AllPatterns *closure,int *LenRegion,Genome &genome)
{
	int len = strlen(curPattern);
	int i;
	int pre = (len+1)/2;

	int size = 1;
	int countWS = 0;
	


	for (i=0; i<pre; i++) {
		char l = curPattern[i];
		char r = curPattern[len-1-i];

		switch(l) {
		case '0': if (!(r=='N' || r=='3' || r=='Y' || r=='W')) return 0;
			break;
		case '1': if (!(r=='N' || r=='2' || r=='R' || r=='S')) return 0;
			break;
		case '2': if (!(r=='N' || r=='1' || r=='Y' || r=='S')) return 0;
			break;
		case '3': if (!(r=='N' || r=='0' || r=='R' || r=='W')) return 0;
			break;
		case 'R': if (!(r=='N' || r=='3' || r=='1' || r=='Y' || r=='S' || r=='W')) return 0;
			if (r=='Y' || r=='N') size *= 2;
			break;
		case 'Y': if (!(r=='N' || r=='2' || r=='0' || r=='R' || r=='S' || r=='W')) return 0;
			if (r=='R' || r=='N') size *=2;
			break;
		case 'W': if (r=='1' || r=='2') return 0; 
			if (r=='N' || r== 'W') {
				size *=2;
				countWS++;
			}
			break;
		case 'S': if (r=='0' || r=='3') return 0;
			if (r=='N' || r=='S') {
				size *=2;
				countWS++;
			}
			break;
		case 'N': if (r=='W' || r=='S' || r=='R' || r=='Y') size *=2;
			if (r=='W' || r=='S') 
				countWS++;
			break;
		default: printf("Invalid character encountered ... aborting\n");
			exit(0);
		}
	}

	if (size < 2 || countWS==0) return 0;

	char **setDups = new char*[size];
	for (i=0; i<size; i++) 
		setDups[i] = new char[len+1];

	int curSize = 1;
	for (i=0; i<pre; i++) {
		char l = curPattern[i];
		char r = curPattern[len-1-i];
		
		switch(l) {
		case '0': 
			for (int j=0; j<curSize; j++) {
				setDups[j][i] = '0';
				setDups[j][len-1-i] = '3';
			}
			break;
		case '1': 
			for (int j=0; j<curSize; j++) {
				setDups[j][i] = '1';
				setDups[j][len-1-i] = '2';
			}
			break;
		case '2': 
			for (int j=0; j<curSize; j++) {
				setDups[j][i] = '2';
				setDups[j][len-1-i] = '1';
			}
			break;
		case '3': 
			for (int j=0; j<curSize; j++) {
				setDups[j][i] = '3';
				setDups[j][len-1-i] = '0';
			}
			break;
		case 'R': 
			if (r=='3' || r=='W') {
				for (int j=0; j<curSize; j++) {
					setDups[j][i] = '0';
					setDups[j][len-1-i] = '3';
				}
				break;
			}
			if (r=='1' || r=='S') {
				for (int j=0; j<curSize; j++) {
					setDups[j][i] = '2';
					setDups[j][len-1-i] = '1';
				}
				break;
			}
			if (r=='Y' || r=='N') {
				for (int j=0; j<curSize; j++) {
					if (curSize>=size)
						printf("Error: size underestimated\n");
					memcpy(setDups[j+curSize],setDups[j],len+1);
					setDups[j][i] = '0'; setDups[j][len-1-i] = '3';
					setDups[j+curSize][i] = '2'; setDups[j+curSize][len-1-i] = '1';
				}
				curSize *= 2;
				break;
			}
			break;
		case 'Y': 
			if (r=='2' || r=='S') {
				for (int j=0; j<curSize; j++) {
					setDups[j][i] = '1';
					setDups[j][len-1-i] = '2';
				}
				break;
			}
			if (r=='0' || r=='W') {
				for (int j=0; j<curSize; j++) {
					setDups[j][i] = '3';
					setDups[j][len-1-i] = '0';
				}
				break;
			}
			if (r=='R' || r=='N') {
				for (int j=0; j<curSize; j++) {
					if (curSize>=size)
						printf("Error: size underestimated\n");
					memcpy(setDups[j+curSize],setDups[j],len+1);
					setDups[j][i] = '3'; setDups[j][len-1-i] = '0';
					setDups[j+curSize][i] = '1'; setDups[j+curSize][len-1-i] = '2';
				}
				curSize *= 2;
				break;
			}
			break;
		case 'W': 
			if (r=='0' || r=='R') {
				for (int j=0; j<curSize; j++) {
					setDups[j][i] = '3';
					setDups[j][len-1-i] = '0';
				}
				break;
			}
			if (r=='3' || r=='Y') {
				for (int j=0; j<curSize; j++) {
					setDups[j][i] = '0';
					setDups[j][len-1-i] = '3';
				}
				break;
			}
			if (r=='N' || r=='W') {
				for (int j=0; j<curSize; j++) {
					if (curSize>=size)
						printf("Error: size underestimated\n");
					memcpy(setDups[j+curSize],setDups[j],len+1);
					setDups[j][i] = '0'; setDups[j][len-1-i] = '0';
					setDups[j+curSize][i] = '3'; setDups[j+curSize][len-1-i] = '3';
				}
				curSize *= 2;
				break;
			}
			break;
		case 'S': 
			if (r=='2' || r=='R') {
				for (int j=0; j<curSize; j++) {
					setDups[j][i] = '1';
					setDups[j][len-1-i] = '2';
				}
				break;
			}
			if (r=='1' || r=='Y') {
				for (int j=0; j<curSize; j++) {
					setDups[j][i] = '2';
					setDups[j][len-1-i] = '1';
				}
				break;
			}
			if (r=='N' || r=='S') {
				for (int j=0; j<curSize; j++) {
					if (curSize>=size)
						printf("Error: size underestimated\n");
					memcpy(setDups[j+curSize],setDups[j],len+1);
					setDups[j][i] = '1'; setDups[j][len-1-i] = '1';
					setDups[j+curSize][i] = '2'; setDups[j+curSize][len-1-i] = '2';
				}
				curSize *= 2;
				break;
			}
			break;
		case 'N': 
			if (r=='0') {
				for (int j=0; j<curSize; j++) {
					setDups[j][i] = '3';
					setDups[j][len-1-i] = '0';
				}
				break;
			}
			if (r=='1') {
				for (int j=0; j<curSize; j++) {
					setDups[j][i] = '2';
					setDups[j][len-1-i] = '1';
				}
				break;
			}
			if (r=='2') {
				for (int j=0; j<curSize; j++) {
					setDups[j][i] = '1';
					setDups[j][len-1-i] = '2';
				}
				break;
			}
			if (r=='3') {
				for (int j=0; j<curSize; j++) {
					setDups[j][i] = '0';
					setDups[j][len-1-i] = '3';
				}
				break;
			}
			if (r=='N') {
				for (int j=0; j<curSize; j++) {
					setDups[j][i] = 'N';
					setDups[j][len-1-i] = 'N';
				}
				break;
			}
			if (r=='R') {
				for (int j=0; j<curSize; j++) {
					memcpy(setDups[j+curSize],setDups[j],len+1);
					setDups[j][i] = '3'; setDups[j][len-1-i] = '0';
					setDups[j+curSize][i] = '1'; setDups[j+curSize][len-1-i] = '2';
				}
				curSize *= 2;
				break;
			}
			if (r=='Y') {
				for (int j=0; j<curSize; j++) {
					memcpy(setDups[j+curSize],setDups[j],len+1);
					setDups[j][i] = '0'; setDups[j][len-1-i] = '3';
					setDups[j+curSize][i] = '2'; setDups[j+curSize][len-1-i] = '1';
				}
				curSize *= 2;
				break;
			}
			if (r=='S') {
				for (int j=0; j<curSize; j++) {
					memcpy(setDups[j+curSize],setDups[j],len+1);
					setDups[j][i] = '1'; setDups[j][len-1-i] = '1';
					setDups[j+curSize][i] = '2'; setDups[j+curSize][len-1-i] = '2';
				}
				curSize *= 2;
				break;
			}
			if (r=='W') {
				for (int j=0; j<curSize; j++) {
					memcpy(setDups[j+curSize],setDups[j],len+1);
					setDups[j][i] = '0'; setDups[j][len-1-i] = '0';
					setDups[j+curSize][i] = '3'; setDups[j+curSize][len-1-i] = '3';
				}
				curSize *= 2;
				break;
			}
			break;
		default: printf("Invalid character encountered ... aborting\n");
			exit(0);
		}
	}


	// the setDups is ready
	if (curSize != size) {
		printf("Internal error in ComputePseudoPalindromeCorrection ... aborting\n");
		printf("cursize = %d size = %d on %s\n",curSize,size,curPattern);
		exit(0);
	}

	double retVal = 0;
	for (i=0; i<size; i++) {
		setDups[i][len-1] = 0;
		int newNumN = countN(setDups[i]);
		retVal += totalLenRegion[newNumN]*ComputeSimpleProbability(setDups[i],genome,newNumN);
	}

	return retVal;
}

	
