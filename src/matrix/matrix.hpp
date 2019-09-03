/*******************************************************
Copyright (c) 2000 by Martin Tompa and Saurabh Sinha.
All rights reserved.  Redistribution is not permitted without the 
express written permission of the authors.
The program YMF implements an algorithm for identifying
likely transcription factor binding sites in yeast, described in the
following paper:
"Statistical Method for Finding Transcription Factor Binding Sites"
by Saurabh Sinha and Martin Tompa,
Eighth International Conference on Intelligent Systems for
Molecular Biology, San Diego, USA, August 2000, 344-354.
*********************************************************/
#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_

class matrix {
public:
  matrix(int i, int j);
  double& operator () (int i, int j);
  friend matrix operator * (const matrix& l, const matrix& r);  
  friend matrix operator + (const matrix& l, const matrix& r);  
  friend matrix operator - (const matrix& l, const matrix& r);
  matrix inv();
  void print();

  matrix(const matrix& other);
  ~matrix();
  matrix& operator = (const matrix& other);

private:
  int m;
  int n;
  double **mat;
};

#endif
