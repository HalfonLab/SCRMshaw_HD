#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "stats.h"
#include "genome.h"

double minBest;
int minBestIndex;
struct results **best;

int NUM_CATEGORY;
int NUM_RYWS;
int MIN_TOTAL_COUNT;
int MAX_NUM_REGIONS;

int categ;
int **category;
int perCategoryFound;
int nextIndex;
FILE *fp;

int *totalLenRegion;
int *totalSquareLenRegion;

int oddLength;

void WritePattern(int level, int k, int &curIndex, char *curPattern, struct AllPatterns *allP);
void GenerateMotifAndTest(int level,int lenOligo,int numRYWS,int numN,char *curPattern,struct AllPatterns *allP,int *LenRegion,int numRegions, Genome &genome,int,int);
void ProcessMotif(char *curPattern,struct AllPatterns *allP,int *LenRegion,int numRegions,Genome &genome,int, int);

int isFilterable(char *pattern) {
  int lenTotal = strlen(pattern);

  int lenActual = 0;
  for (int i=0; i<lenTotal; i++)
    if (pattern[i] != 'N') lenActual++;

  int numRY = 0;
  for (int i=0; i<lenTotal; i++) 
    if (pattern[i] == 'R' || pattern[i] == 'Y') numRY++;

  int allT = 0;
  for (int i=0; i<lenTotal; i++) 
    if (pattern[i] != '3' && pattern[i] != 'N' && pattern[i] != 'Y') allT++;
  if (allT <= 1) return 1;

  int allA = 0;
  for (int i=0; i<lenTotal; i++) 
    if (pattern[i] != '0' && pattern[i] != 'N' && pattern[i] != 'R') allA++;
  if (allA <= 1) return 1;

  int polyAT0 = 0;
  int expecting = 0;
  for (int i=0; i<lenTotal; i++) {
    if (pattern[i] != (expecting+'0') && pattern[i] != 'N'
	&& pattern[i] != (expecting==0?'R':'Y')) polyAT0++;
    expecting = -(expecting - 3); // ha ha
  }
  if (polyAT0 <= 1) return 1;

  int polyAT1 = 0;
  expecting = 3;
  for (int i=0; i<lenTotal; i++) {
    if (pattern[i] != (expecting+'0') && pattern[i] != 'N'
	&& pattern[i] != (expecting==0?'R':'Y')) polyAT1++;
    expecting = -(expecting - 3); // ha ha
  }
  if (polyAT1 <= 1) return 1;

  // mathieu
  // remove poly(AT)
  int AT=0;
  for (int i=0; i<lenTotal; i++) {

    if (pattern[i] != 'N' && pattern[i] != '0' && pattern[i] != '3' && pattern[i]!='W') AT++;


  }
  if (!AT) return 1; 
  return 0;
}


void WritePattern(int level, int k, int &curIndex, char *curPattern, struct AllPatterns *allP)
{
  char *dummy;
  if (level == k) {
    allP[curIndex].pattern = new char[k+1];
    allP[curIndex].count = 0;
    allP[curIndex].count2 = 0;
    strcpy(allP[curIndex].pattern,curPattern);
    allP[curIndex++].value = strtol(curPattern,&dummy,4);
    return;
  }
	
  char *pattern = new char[level+2];
  strcpy(pattern,curPattern);
  pattern[level] = '0';
  pattern[level+1] = 0;
  WritePattern(level+1,k,curIndex,pattern,allP);
  pattern[level] = '1';
  WritePattern(level+1,k,curIndex,pattern,allP);
  pattern[level] = '2';
  WritePattern(level+1,k,curIndex,pattern,allP);
  pattern[level] = '3';
  WritePattern(level+1,k,curIndex,pattern,allP);
  delete [] pattern;
	
  return;
}

Genome::Genome()
{
  numChrom = 0;
  for (int i=0; i < 64; i++) {
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

int main(int argc, char **argv) {
		
  if (argc < 6) {
    printf("Usage : %s configfile lenRegion lenOligo pathoforganismtables [-sort] [-UB=<maxoccurrencesperregion>:not supported] RegionFile1 RegionFile2  ...\n",argv[0]);
    exit(0);
  }
	
  // then initialize and read in the Markov chain data.
  Genome genome;
	
  char organism[1024];
  strcpy(organism,argv[4]);
  strcat(organism,"_table.3");
  printf("Reading tables from %s ...\n",organism);
  std::fstream fs(organism,std::ios::in);
  if (fs.is_open()==0) {
    printf("Error opening file %s\n",organism);
    printf("Error opening file %s\nThis file is needed for correct execution of the software\nand should have come with the software package\nIf you cannot locate this file, and do not have the tools\nto create it (preproc), please write to saurabh@cs.washington.edu\nAborting now...\n",organism);
	
    exit(0);
  }


  int i;
  for (i=0; i<64; i++) {
    for (int j=0; j<64; j++) {
      fs >> genome.P[i][j];
    }
  }
	
  for (i=0; i<64; i++) {
    fs >> genome.p[i];
  }
	
  for (i=0; i<64; i++) {
    for (int j=0; j<64; j++) {
      fs >> genome.Q[i][j];
    }
  }
  for (i=0; i<64; i++) {
    for (int j=0; j<64; j++) {
      fs >> genome.QPPP[i][j];
    }
  }
	
  for (i=0; i<64; i++) {
    for (int j=0; j<64; j++) {
      fs >> genome.QPPQPP[i][j];
    }
  }
	
	
  char *dummy;
	
  int lenRegion = strtol(argv[2], &dummy, 10);
  if (lenRegion < 10) {
    printf("Short upstream region ... Aborting\n");
    exit(0);
  }
	
  int lenOligo = strtol(argv[3], &dummy, 10);
  if (lenOligo < 6) {
    printf("Short oligo ... Aborting\n");
    exit(0);
  }
	
  oddLength = 0;
  if ((lenOligo%2)!=0) oddLength =1;
	
  // read in the config file also 
  FILE *config = fopen(argv[1],"r");
  int MAX_LINE = 128;
  char line[MAX_LINE];
	
  if (!(fgets(line, MAX_LINE, config))) {
    printf("Error reading config file line 1\n");
    exit(0);
  }
  sscanf(line, "%d", &NUM_CATEGORY);
  int begin, end, per;
  int maxNumN = 0;
	
  category = new int *[NUM_CATEGORY];
  for (i=0; i<NUM_CATEGORY; i++) {
    category[i] = new int[3];	
    if (!(fgets(line, MAX_LINE, config))) {
      printf("Error reading config file line %d\n",1+i);
      exit(0);
    } 
    sscanf(line, "%d..%d %d", &begin, &end, &per);
    printf("begin = %d end = %d per = %d\n",begin,end,per);
    if (begin > end) {
      printf("Category begin is > Category end on line %d in config file\n",1+i);
      exit(0);
    }
    if (maxNumN < end)
      maxNumN = end;
    category[i][0] = begin;
    category[i][1] = end;
    category[i][2] = per;
  }
	
  if (!(fgets(line, MAX_LINE, config))) {
    printf("Error reading config file line %d\n",1+i+1);
    exit(0);
  }
  sscanf(line, "%d", &NUM_RYWS);
  printf("NUM_RYWS = %d\n",NUM_RYWS);
  printf("This means that the space of motifs examined contains motifs with at most %d characters from the set R,Y,W,S\n",NUM_RYWS);
	
  if (!(fgets(line, MAX_LINE, config))) { 
    printf("Error reading config file line %d\n",1+i+2);
    exit(0);
  }
  sscanf(line, "%d", &MIN_TOTAL_COUNT);
  printf("MIN_TOTAL_COUNT = %d\n",MIN_TOTAL_COUNT);
  printf("This means that the space of motifs examined contains motifs that occur at least %d times in the input set of regions\n",MIN_TOTAL_COUNT);
	
  if (!(fgets(line, MAX_LINE, config))) { 
    printf("Error reading config file line %d\n",1+i+3);
    exit(0);
  }
  sscanf(line, "%d", &MAX_NUM_REGIONS);
  printf("NUM_REG = %d\n",MAX_NUM_REGIONS);
  printf("This means that at most %d regions have been input for analysis. This number is meant only to be an upper bound on the number of regions. Also, it has nothing to do with the number of files in which the regions have been included\n",MAX_NUM_REGIONS);
	
		
  int numRYWS = NUM_RYWS;
	
  int noN;
	
  // the powers of P ... read in the preprocessed data
	
	// P^1
  for (int i1=0; i1<64; i1++) {
    for (int i2=0; i2<64; i2++) {
      genome.P_numNplus1[0][i1][i2] = genome.P[i1][i2];
    }
  }
	
	
  strcpy(organism,argv[4]);
  strcat(organism,"_powers.3");
  printf("Reading powers of P from %s ...\n",organism);
  std::fstream fsP(organism,std::ios::in);
  if (fsP.is_open()==0) {
    printf("Error opening file %s\n",organism);
    printf("Error opening file %s. This file is needed for correct execution of the software and should have come with the software package. If you cannot locate this file, and do not have the tools. to create it (preproc), please write to saurabh@cs.washington.edu. Aborting now...\n",organism);
    exit(0);
  }
	
  for (noN=1; noN<=maxNumN+3; noN++) {
    for (int i=0; i<64; i++) {
      for (int j=0; j<64; j++) {
	fsP >> genome.P_numNplus1[noN][i][j];
      }
    }
  }
	
  // the genome data read in. Now read in the input region set
	
  char **Upstream = new char *[MAX_NUM_REGIONS];
  printf("Will read at most %d regions of length at most %d each for an oligo length of %d\n",MAX_NUM_REGIONS,lenRegion,lenOligo);

  bool sortresults = false;
  int basearg = 5;
  if (strstr(argv[basearg],"sort")) {
    basearg++;
    sortresults = true;
  }

  int *LenRegion = new int[MAX_NUM_REGIONS];

  int numRegions=-1;
  int ptr = 0;

  for (i=0; i<(argc-basearg); i++) {
    FILE *fpG = fopen(argv[basearg+i], "r");
		
    if (fpG == NULL) {
      printf("Error opening sequence data file %s ... Aborting\n",argv[basearg+i]);
      exit(0);
    }
	  	    
    int ch = 0;
    while ((ch = fgetc(fpG)) != EOF){
      if (ch == '>') {	
				// skip to end of line
	while ((ch = fgetc(fpG)) != EOF && ch != '\n' && ch != '\r');
	if (ch == EOF) break;
				// get ready for a new region
	if (numRegions > -1) {
	  LenRegion[numRegions] = ptr;
	  Upstream[numRegions++][ptr] = 0;
	  ptr = 0;
	  Upstream[numRegions] = new char[lenRegion+1];
	}
	else {
	  numRegions = 0;
	  ptr = 0;
	  Upstream[numRegions] = new char[lenRegion+1];
	}
      }
      if (ch != '\n' && ch != '\r') {
	if (ptr==lenRegion) continue;
	switch (ch) {
	case 'a':
	case 'A': Upstream[numRegions][ptr++] = '0';break;
	case 'c':
	case 'C': Upstream[numRegions][ptr++] = '1';break;
	case 'g':
	case 'G': Upstream[numRegions][ptr++] = '2';break;
	case 't':
	case 'T': Upstream[numRegions][ptr++] = '3';break;
	case 'n':
	case 'N': break;
	default: printf("Invalid character %c encountered in data file ... Aborting\n",ch);
	  exit(0);
	}
      }
    }
  }
  if (numRegions > -1) {
    LenRegion[numRegions] = ptr;
    Upstream[numRegions++][ptr] = 0;
  }
  printf("%d valid input regions have been read in\n",numRegions);
	

	// initialize the results data structure

  best = new struct results *[NUM_CATEGORY];	
  for (int ind=0; ind<NUM_CATEGORY; ind++) {
    best[ind] = new struct results[category[ind][2]];
    for (int pc=0; pc<category[ind][2]; pc++) {
      best[ind][pc].pattern = new char[maxNumN+lenOligo+1];
      best[ind][pc].cnt = 0;
      best[ind][pc].exp = 0;
      best[ind][pc].var = 0;
      best[ind][pc].zsc = -1000;
    }
  }
  minBest = -1000;
  minBestIndex = 0;
	
  // n-l+1 and sq(n-l+1) for each spacer length 0..25

  totalLenRegion = new int[26];
  totalSquareLenRegion = new int[26];
  for (int i=0; i<26; i++) {
    totalLenRegion[i] = 0;
    totalSquareLenRegion[i] = 0;
    for (int j= 0; j<numRegions; j++) {
      totalLenRegion[i] += (LenRegion[j]-lenOligo-i+1);
      totalSquareLenRegion[i] += (LenRegion[j]-lenOligo-i+1)*(LenRegion[j]-lenOligo-i+1);
    }
  }
	
  // generate all patterns, and store them in the array allP of records.
	
  int numAllP = (int)(pow(4,lenOligo)); 
  struct AllPatterns *allP = new struct AllPatterns[numAllP];
  int curIndex = 0;
  WritePattern(0,lenOligo,curIndex,"",allP);
	

	// all set ! 
	// now do the stats for each value of numN
	
  int numN;
  int pre,suf; 
  char *pref, *suff; 

  // fpCategory = NULL;
	
  for (categ=0; categ < NUM_CATEGORY; categ++) {
    perCategoryFound = 0;
    nextIndex = 0;
    minBest = -1000;
    for (numN=category[categ][0]; numN<=category[categ][1]; numN++) {
      printf("Processing patterns of length %d\n",numN+lenOligo);

    // first do a count of all the patterns
	  
      int i;
      for (i=0; i<numAllP; i++) {
	allP[i].count = 0;
	allP[i].count2 = 0;
      }
	  
      int length = lenOligo + numN;	
      char *tmp = new char[length+1];
	  
	  
      // first do the count for pre <= suf.
      pre = lenOligo/2;
      suf = lenOligo - pre;
      pref = new char[pre+1];
      suff = new char[suf+1];
	  
      long suf2 = 2*suf;
      long modulusPre  = (long)pow(4,pre-1);
      long modulusSuf  = (long)pow(4,suf-1);
      int Nsuf = numN + suf;
      for (i=0; i<numRegions; i++) {
	    
	// get ready for processing this region
	    
	memcpy(tmp,&(Upstream[i][0]),length);
	tmp[length] = 0;
	memcpy(pref,tmp,pre);
	pref[pre] = 0;
	long prefix = strtol(pref,&dummy,4);
	memcpy(suff,&(tmp[length-suf]),suf);		
	suff[suf] = 0;		
	long suffix = strtol(suff,&dummy,4);
	long value = (prefix<<suf2) + suffix;
	allP[value].count++;
	    
	// now process the region 
	    
	for (int j=1; j<(LenRegion[i]-length+1); j++) {  
	  int j1 = j+pre-1;
	  prefix = ((prefix%modulusPre)<<2) + (Upstream[i][j1]-'0');
	  suffix = ((suffix%modulusSuf)<<2) + (Upstream[i][j1+Nsuf]-'0');
	  long value = (prefix<<suf2) + suffix;				
	  allP[value].count++;
	}
      }
	  
      delete [] pref;
      delete [] suff;
	  
      if (oddLength) {	// do the above again, but with reversed pre and suf
	// and update count2 now.
	// this handles the case pre > suf
	suf = lenOligo/2;
	pre = lenOligo - suf;
	pref = new char[pre+1];
	suff = new char[suf+1];
	    
	long suf2 = 2*suf;
	long modulusPre  = (long)pow(4,pre-1);
	long modulusSuf  = (long)pow(4,suf-1);
	int Nsuf = numN + suf;
	for (i=0; i<numRegions; i++) {
	      
				// get ready for processing this region
				
	  memcpy(tmp,&(Upstream[i][0]),length);
	  tmp[length] = 0;
	  memcpy(pref,tmp,pre);
	  pref[pre] = 0;
	  long prefix = strtol(pref,&dummy,4);
	  memcpy(suff,&(tmp[length-suf]),suf);		
	  suff[suf] = 0;		
	  long suffix = strtol(suff,&dummy,4);
	  long value = (prefix<<suf2) + suffix;
	  allP[value].count2++;
				
				// now process the region 
				
	  for (int j=1; j<(LenRegion[i]-length+1); j++) {  
	    int j1 = j+pre-1;
	    prefix = ((prefix%modulusPre)<<2) + (Upstream[i][j1]-'0');
	    suffix = ((suffix%modulusSuf)<<2) + (Upstream[i][j1+Nsuf]-'0');
	    long value = (prefix<<suf2) + suffix;  
	    allP[value].count2++;
	  }
	}
	delete [] pref;
	delete [] suff;
      }
      delete [] tmp;
	  
      printf("Done counting for patterns of length %d\n",length);
      // generate all motifs with numN Ns and test them ....
      GenerateMotifAndTest(0,length,numRYWS,numN,"",allP,LenRegion,numRegions,genome,numRYWS,numN);		
      printf("Done stats for patterns of length %d\n",length);
    }
  }
	
  fp = fopen("results","w");
  for (categ = 0; categ<NUM_CATEGORY; categ++) {
    fprintf(fp,"The best %d candidates in category %d are :\n",category[categ][2],categ+1);
    for (int perCat = 0; perCat < category[categ][2]; perCat++) {
      char * curPattern = best[categ][perCat].pattern;
      int notEnded = 1;
      for (int j=0; (j<lenOligo+maxNumN+1) && (notEnded==1); j++) {
	switch(curPattern[j]) {
	case  0 : notEnded = 0; break;
	case '0': fprintf(fp,"A");break;
	case '1': fprintf(fp,"C");break;
	case '2': fprintf(fp,"G");break;
	case '3': fprintf(fp,"T");break;
	case 'R': fprintf(fp,"R");break;
	case 'Y': fprintf(fp,"Y");break;
	case 'W': fprintf(fp,"W");break;
	case 'S': fprintf(fp,"S");break;
	case 'N': fprintf(fp,"N");break;
	default: printf("Internal error : Unexpected character \n");
	  exit(0);
	}
      }
      fprintf(fp," %d %.2f %.2f %.2f\n", best[categ][perCat].cnt, best[categ][perCat].zsc, best[categ][perCat].exp, best[categ][perCat].var );
    }
  }
	
  fclose(fp);

  if (sortresults) {
    system("sort -g -k 3 -r results >__tmp__\n");
    system("mv __tmp__ results\n");
  }

  return 0;
}

void GenerateMotifAndTest(int level,int lenOligo,int numRYWS,int numN,char *curPattern,struct AllPatterns *allP,int *LenRegion,int numRegions, Genome &genome,int orignumRYWS, int orignumN)
  // generates all motifs with pre <= suf
{
  if (level == lenOligo) {
    ProcessMotif(curPattern,allP,LenRegion,numRegions,genome,orignumRYWS,orignumN);
    return;
  }
	
  char *newPattern = new char[lenOligo+1];
  strcpy(newPattern,curPattern);
	
  if (level == ((lenOligo-orignumN)/2)  && orignumN != 0) {
    for (int i=0; i<orignumN; i++)
      strcat(newPattern,"N");
    GenerateMotifAndTest(level+numN,lenOligo,numRYWS,0,newPattern,allP,LenRegion,numRegions,genome,orignumRYWS,orignumN);
  }
  else {
    if (numRYWS > 0) {
      strcat(newPattern,"R");
      GenerateMotifAndTest(level+1,lenOligo,numRYWS-1,numN,newPattern,allP,LenRegion,numRegions,genome,orignumRYWS,orignumN);
      newPattern[level] = 0;
      strcat(newPattern,"Y");
      GenerateMotifAndTest(level+1,lenOligo,numRYWS-1,numN,newPattern,allP,LenRegion,numRegions,genome,orignumRYWS,orignumN);
      newPattern[level] = 0;
      strcat(newPattern,"W");
      GenerateMotifAndTest(level+1,lenOligo,numRYWS-1,numN,newPattern,allP,LenRegion,numRegions,genome,orignumRYWS,orignumN);
      newPattern[level] = 0;
      strcat(newPattern,"S");
      GenerateMotifAndTest(level+1,lenOligo,numRYWS-1,numN,newPattern,allP,LenRegion,numRegions,genome,orignumRYWS,orignumN);

    }
    newPattern[level] = 0;	// restore newPattern to curPattern
			
    strcat(newPattern,"0");
    GenerateMotifAndTest(level+1,lenOligo,numRYWS,numN,newPattern,allP,LenRegion,numRegions,genome,orignumRYWS,orignumN);
    newPattern[level] = 0;
		
		
    strcat(newPattern,"1");
    GenerateMotifAndTest(level+1,lenOligo,numRYWS,numN,newPattern,allP,LenRegion,numRegions,genome,orignumRYWS,orignumN);
    newPattern[level] = 0;
		
		
    strcat(newPattern,"2");
    GenerateMotifAndTest(level+1,lenOligo,numRYWS,numN,newPattern,allP,LenRegion,numRegions,genome,orignumRYWS,orignumN);
    newPattern[level] = 0;
		
		
    strcat(newPattern,"3");
    GenerateMotifAndTest(level+1,lenOligo,numRYWS,numN,newPattern,allP,LenRegion,numRegions,genome,orignumRYWS,orignumN);
    newPattern[level] = 0;
  }
	
  delete [] newPattern;
  return;
}

int shouldBeProcessed(char *curPattern,int lenOligo)
{
  char *revComp = new char[lenOligo+1];
  for (int i=0; i<lenOligo; i++) {
    switch(curPattern[i]) {
    case '0': revComp[lenOligo-1-i] = '3';break;
    case '1': revComp[lenOligo-1-i] = '2';break;
    case '2': revComp[lenOligo-1-i] = '1';break;
    case '3': revComp[lenOligo-1-i] = '0';break;
    case 'R': revComp[lenOligo-1-i] = 'Y';break;
    case 'Y': revComp[lenOligo-1-i] = 'R';break;
    case 'W': revComp[lenOligo-1-i] = 'W';break;
    case 'S': revComp[lenOligo-1-i] = 'S';break;
    case 'N': revComp[lenOligo-1-i] = 'N';break;
    default: printf("Invalid character encountered : aborting\n");exit(0);
    }
  }

  revComp[lenOligo] = 0;
  if (strcmp(curPattern,revComp)<=0) {
    delete [] revComp;
    return 1;
  }

  delete [] revComp;
  return 0;
}

void ProcessMotif(char *curPattern,struct AllPatterns *allP,int *LenRegion,int numRegions,Genome &genome,int numRYWS, int numN)
{
  // if (isFilterable(curPattern)) {
  //   return;
  // }

  int lenOligo = strlen(curPattern);

	// do not process for motifs if their reverse will be processed

  if (!shouldBeProcessed(curPattern,lenOligo))
    return;

  long maxSizeClosure = (long)pow(2,numRYWS);
	
  maxSizeClosure *= 2;
  struct AllPatterns *closure = new struct AllPatterns[maxSizeClosure];
	
  for (int c=0; c<maxSizeClosure; c++) {
    closure[c].count = 0;
    closure[c].count2 = 0;
    if (!oddLength)
      closure[c].preOrSuf = 0;
    closure[c].value = 0;
  }
	
  int curSize = 0;
	
  if (curPattern[0] == '0') {
    closure[0].pattern = new char[lenOligo+1];
    closure[0].pattern[0] = '0';
    curSize = 1;
  }
	
  if (curPattern[0] == '1') {
    closure[0].pattern = new char[lenOligo+1];
    closure[0].pattern[0] = '1';
    curSize = 1;
  }
  if (curPattern[0] == '2') {
    closure[0].pattern = new char[lenOligo+1];
    closure[0].pattern[0] = '2';
    curSize = 1;
  }
  if (curPattern[0] == '3') {
    closure[0].pattern = new char[lenOligo+1];
    closure[0].pattern[0] = '3';
    curSize = 1;
  }
  if (curPattern[0] == 'N') {
    closure[0].pattern = new char[lenOligo+1];
    closure[0].pattern[0] = 'N';
    curSize = 1;
  }
	
  if (curPattern[0] == 'R') {
    closure[0].pattern = new char[lenOligo+1];
    closure[0].pattern[0] = '0';
    closure[1].pattern = new char[lenOligo+1];
    closure[1].pattern[0] = '2';
    curSize = 2;
  }
	
  if (curPattern[0] == 'Y') {
    closure[0].pattern = new char[lenOligo+1];
    closure[0].pattern[0] = '1';
    closure[1].pattern = new char[lenOligo+1];
    closure[1].pattern[0] = '3';
    curSize = 2;
  }

  if (curPattern[0] == 'W') {
    closure[0].pattern = new char[lenOligo+1];
    closure[0].pattern[0] = '0';
    closure[1].pattern = new char[lenOligo+1];
    closure[1].pattern[0] = '3';
    curSize = 2;
  }
	
  if (curPattern[0] == 'S') {
    closure[0].pattern = new char[lenOligo+1];
    closure[0].pattern[0] = '1';
    closure[1].pattern = new char[lenOligo+1];
    closure[1].pattern[0] = '2';
    curSize = 2;
  }
	
  int i;
  for (i=1; i<lenOligo; i++) {
    if (curPattern[i] == '0') {
      for (int j=0; j<curSize; j++) 
	closure[j].pattern[i] = '0';
    }
    if (curPattern[i] == '1') {
      for (int j=0; j<curSize; j++) 
	closure[j].pattern[i] = '1';
    }
    if (curPattern[i] == '2') {
      for (int j=0; j<curSize; j++) 
	closure[j].pattern[i] = '2';
    }
    if (curPattern[i] == '3') {
      for (int j=0; j<curSize; j++) 
	closure[j].pattern[i] = '3';
    }
    if (curPattern[i] == 'N') {
      for (int j=0; j<curSize; j++) 
	closure[j].pattern[i] = 'N';
    }
    if (curPattern[i] == 'R') {
      for (int j=0; j<curSize; j++) {
	closure[j].pattern[i] = '0';
	closure[curSize+j].pattern = new char[lenOligo+1];
	memcpy(closure[curSize+j].pattern,closure[j].pattern,i);
	closure[curSize+j].pattern[i] = '2';
      }
      curSize *= 2;
    }
    if (curPattern[i] == 'Y') {
      for (int j=0; j<curSize; j++) {
	closure[j].pattern[i] = '1';
	closure[curSize+j].pattern = new char[lenOligo+1];
	memcpy(closure[curSize+j].pattern,closure[j].pattern,i);
	closure[curSize+j].pattern[i] = '3';
      }
      curSize *= 2;
    }
    if (curPattern[i] == 'W') {
      for (int j=0; j<curSize; j++) {
	closure[j].pattern[i] = '0';
	closure[curSize+j].pattern = new char[lenOligo+1];
	memcpy(closure[curSize+j].pattern,closure[j].pattern,i);
	closure[curSize+j].pattern[i] = '3';
      }
      curSize *= 2;
    }
    if (curPattern[i] == 'S') {
      for (int j=0; j<curSize; j++) {
	closure[j].pattern[i] = '1';
	closure[curSize+j].pattern = new char[lenOligo+1];
	memcpy(closure[curSize+j].pattern,closure[j].pattern,i);
	closure[curSize+j].pattern[i] = '2';
      }
      curSize *= 2;
    }

  }
	
  int sizeClosure = curSize;
	
  for (i=0; i<sizeClosure; i++) {
    closure[i].pattern[lenOligo] = 0;
    if (oddLength)		// till now, since we have followed the original motif, pre < suf
      closure[i].preOrSuf = -1;
  }
	
  for (i=0; i<sizeClosure; i++) {
    closure[sizeClosure+i].pattern = new char[lenOligo+1];
    if (oddLength)
      closure[sizeClosure+i].preOrSuf = 1;
    for (int j=0; j<lenOligo; j++) {
      switch (closure[i].pattern[lenOligo-1-j]) {
      case '0': closure[i+sizeClosure].pattern[j] = '3'; break;
      case '1': closure[i+sizeClosure].pattern[j] = '2'; break;
      case '2': closure[i+sizeClosure].pattern[j] = '1'; break;
      case '3': closure[i+sizeClosure].pattern[j] = '0'; break;
      case 'N': closure[i+sizeClosure].pattern[j] = 'N'; break;
      default : printf("Invalid character %c internally .. aborting\n",closure[sizeClosure+i].pattern[j]);
	exit(0);
      }
    }
    closure[sizeClosure+i].pattern[lenOligo] = 0;
  }
  sizeClosure *= 2;
	
	
	// the entire closure has been computed
	
	// now computing the count of each member of the closure
	

  long totalCount = 0;
  char *dummy;
  for (i=0; i<sizeClosure; i++) {
    int pre = 0;
    if ((numN > 0) && oddLength) {
      while (closure[i].pattern[pre]!='N')
	pre++;
    }
    else 
      pre = (lenOligo-numN)/2;
    int suf = lenOligo - numN - pre;
		
    char *pref = new char[pre+1];
    memcpy(pref,closure[i].pattern,pre);
    pref[pre] = 0;
		
    long prefix = strtol(pref,&dummy,4);
		
    char *suff = new char[suf+1];
		
    memcpy(suff,&(closure[i].pattern[lenOligo-suf]),suf);		
    suff[suf] = 0;
		
    long suffix = strtol(suff,&dummy,4);
    long sufPower = (long)pow(4,suf);
    long value = prefix*sufPower + suffix;
	
    if (closure[i].preOrSuf <=0) {
      totalCount += allP[value].count;
    }
    else {
      totalCount += allP[value].count2;
    }
    delete [] pref;
    delete [] suff;
  }
	
  if (totalCount < MIN_TOTAL_COUNT) {		// THIS is where the pruning is done
    for (int iprime=0; iprime<sizeClosure;iprime++) 
      delete [] closure[iprime].pattern;
    delete [] closure;
    return;
  }
	
  double nW = totalCount;
	

  double PnW = ComputeProbability(sizeClosure,closure,genome,numN);
  double EnW = PnW*totalLenRegion[numN];
  double evW = EnW - totalSquareLenRegion[numN]*PnW*PnW;			// make a first estimate of vW
  double ese;
  if (evW > 0) 
    ese = sqrt(evW);
  else 
    ese = -1;
  double estimatedZScore = (nW-EnW)/(ese);
  if (nW <= EnW || ((evW > 0) && estimatedZScore < minBest)) {// negative or poor z-score, so dont compute variance
    for (i=0; i<sizeClosure;i++) 
      delete [] closure[i].pattern;
    delete [] closure;
    return;
  }
	
  evW = ComputeAlmostVariance(sizeClosure,closure,LenRegion, genome, PnW, numN, numRegions);
  if (evW > 0) 
    ese = sqrt(evW);
  else 
    ese = -1;
  estimatedZScore = (nW-EnW)/(ese);

  if ((evW > 0) && (estimatedZScore < minBest)) {// negative or poor z-score, so dont compute variance
    for (i=0; i<sizeClosure;i++) 
      delete [] closure[i].pattern;
    delete [] closure;
    return;
  }
	
  double vW  = evW + 2*ComputeOverlapExpectation(sizeClosure,closure,LenRegion, genome);
  vW += 2*ComputePseudoPalindromeCorrection(curPattern,sizeClosure,closure,LenRegion,genome);
				// this last term was introduced by the W,S characters
  if (vW <= 0) {
    printf("Internal error on %s : negative variance \n",curPattern);
    printf("%s %f %f %f\n", curPattern, nW, EnW, vW);
    for (int j=0; j<sizeClosure;j++) 
      delete [] closure[j].pattern;
    delete [] closure;
    return;
  }
  double se = sqrt(vW);
  double z_score = (nW-EnW)/(se);
	
  int index;
  if (z_score > minBest) { // if this z_score is greater than the minimum of the z-scores found so far

    if (perCategoryFound)
      index = minBestIndex;
    else {
      index = nextIndex++;
      if (nextIndex == category[categ][2])
	perCategoryFound=1;
    }

    strcpy(best[categ][index].pattern,curPattern);
    best[categ][index].exp = EnW;
    best[categ][index].var = vW;
    best[categ][index].zsc = z_score;
    best[categ][index].cnt = (int)nW;
		
    minBest = z_score;
    minBestIndex = index;

    if (perCategoryFound) {
      for (int ind = 0; ind<category[categ][2]; ind++) {
	if ((best[categ][ind].zsc > - 1000) && (best[categ][ind].zsc < minBest)) {
	  minBest = best[categ][ind].zsc;
	  minBestIndex = ind;
	}
      }
    }
    else {
      minBest = -1000;
    }
  }


  for (i=0; i<sizeClosure;i++) 
    delete [] closure[i].pattern;
  delete [] closure;

	

	
  return;
}





