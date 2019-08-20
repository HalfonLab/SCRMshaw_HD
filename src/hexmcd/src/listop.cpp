#include "listop.h"
#include <string>
#include <cstring>
#include <stdlib.h> 
#include <cstdlib>
#include <algorithm>

std::list<char *> *Split(char *toks, char * str)
{
    std::list<char *> *L = new std::list<char *>;
  if (str == NULL) return L;
  
  char *head = strtok(str,toks);
  L->push_back(head);
  while (char *tok = strtok(NULL, toks)) {
    L->push_back(tok);
  }
  
  return L;
}


char *Join(char *tok, std::list<char *> *l) 
{
  int *lengths = new int[l->size()];
  int totalLength = 0;
  // first compute the lengths of each token in list
  int i = 0;
  for (std::list<char *>::iterator it = l->begin(); it != l->end(); it++) {
    char *cur_tok = (char *)*it;
    lengths[i] = strlen(cur_tok);
    totalLength += lengths[i];
    i++;    
  }
  if (i==0) {
    char *r = new char[1];
    r[0] = 0;
    return r;
  }

  // allocate memory for string
  char *retstr = new char[totalLength+l->size()*strlen(tok)+1];

  int seplen = strlen(tok);
  int sptr = 0;
  // now concatenate
  i = 0;
  for (std::list<char *>::iterator it = l->begin(); it != l->end(); it++) {
    char *cur_tok = (char *)*it;
    if (i>0) {
      for (int j=0; j<seplen; ) retstr[sptr++] = tok[j++];
    }
    for (int j=0; j<lengths[i]; ) retstr[sptr++] = cur_tok[j++];
    i++;    
  }  
  retstr[sptr] = 0;

  delete [] lengths;

  return retstr;
}

char *Join(char *tok, std::list<string> *l) 
{
  int *lengths = new int[l->size()];
  int totalLength = 0;
  // first compute the lengths of each token in list
  int i = 0;
  for (std::list<string>::iterator it = l->begin(); it != l->end(); it++) {
    string cur_tok = (string)*it;
    lengths[i] = cur_tok.length();
    totalLength += lengths[i];
    i++;    
  }
  if (i==0) {
    char *r = new char[1];
    r[0] = 0;
    return r;
  }

  // allocate memory for string
  char *retstr = new char[totalLength+l->size()*strlen(tok)+1];

  int seplen = strlen(tok);
  int sptr = 0;
  // now concatenate
  i = 0;
  for (std::list<string>::iterator it = l->begin(); it != l->end(); it++) {
    char *cur_tok = (char *)(it->c_str());
    if (i>0) {
      for (int j=0; j<seplen; ) retstr[sptr++] = tok[j++];
    }
    for (int j=0; j<lengths[i]; ) retstr[sptr++] = cur_tok[j++];
    i++;    
  }  
  retstr[sptr] = 0;

  delete [] lengths;

  return retstr;
}

char *Shift(std::list<char *> *l)
{
  char *r = l->front();
  l->pop_front();
  return r;
}

std::list<int> *StringListToIntList(std::list<char *> *l)
{
    std::list<int> *rlist = new std::list<int>;
  for (std::list<char *>::iterator it = l->begin(); it != l->end(); it++) {
    char *tok = (char *)(*it);
    rlist->push_back(atoi(tok));
  }
  return rlist;
}

std::vector<int> *StringListToIntVector(std::list<char *> *l) {
    std::vector<int> *rvec = new std::vector<int>;
  for (std::list<char *>::iterator it = l->begin(); it != l->end(); it++) {
    char *tok = (char *)(*it);
    rvec->push_back(atoi(tok));
  }
  return rvec;
}

void Chomp(char *line)
{
  int linelen = strlen(line);
  if (line[linelen-1] == '\n') line[linelen-1] = 0;
}

char* itoa( int value, char* result, int base ) {
  
  // check that the base if valid  
  if (base < 2 || base > 16) { *result = 0; return result; }
  
  char* out = result;  
  int quotient = value;
  
  do {    
    *out = "0123456789abcdef"[ std::abs( quotient % base ) ];    
    ++out;    
    quotient /= base;    
  } while ( quotient );
    
  // Only apply negative sign for base 10  
  if ( value < 0 && base == 10) *out++ = '-';     
  std::reverse( result, out );  
  *out = 0;
  
  return result;  
}
