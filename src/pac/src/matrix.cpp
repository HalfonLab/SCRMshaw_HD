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
#include "matrix.hpp"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void lubksb(double **a, int n, int *indx, double b[]);
static void ludcmp(double **a, int n, int *indx, double *d);

matrix::matrix(int i, int j)
{
  m = i; n = j;
  mat = new double *[m];
  for (int k=0; k<m; k++) {
    mat[k] = new double[n];
  }
}

double& matrix::operator() (int i, int j)
{
  return mat[i-1][j-1];
}

matrix operator* (const matrix& l, const matrix& r)
{
  assert(l.n == r.m);
  matrix ret(l.m,r.n);
  for (int i=0; i<ret.m; i++) {
    for (int j=0; j<ret.n; j++) {
      double sum = 0;
      for (int k=0; k<l.n; k++) {
	sum += (l.mat[i][k]*r.mat[k][j]);
      }
      ret.mat[i][j] = sum;
    }
  }
  return ret;
}

matrix operator+ (const matrix& l, const matrix& r)
{
  assert(l.m == r.m);
  assert(l.n == r.n);
  matrix ret(l.m,l.n);
  for (int i=0; i<ret.m; i++) {
    for (int j=0; j<ret.n; j++) {
      ret.mat[i][j] = l.mat[i][j]+r.mat[i][j];
    }
  }
  return ret;
}

matrix operator- (const matrix& l, const matrix& r)
{
  assert(l.m == r.m);
  assert(l.n == r.n);
  matrix ret(l.m,l.n);
  for (int i=0; i<ret.m; i++) {
    for (int j=0; j<ret.n; j++) {
      ret.mat[i][j] = l.mat[i][j]-r.mat[i][j];
    }
  }
  return ret;
}

void matrix::print()
{
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      printf("%f ",mat[i][j]);
    }
    printf("\n");
  }
  return;
}

matrix matrix::inv()
{
  assert(m==n);
  matrix ret(n,n);

  // first copy the matrix into temporary space, since the procedure is in-place
  double **a = new double*[n];
  int j;
  for (j=0; j<n; j++) {
    a[j] = new double[n];
    for (int i=0; i<n; i++) 
      a[j][i] = mat[j][i];
  }

  // now the actual inversion
  double **y = ret.mat;
  double *col = new double[n];
  int *indx = new int[n];;
  double d;

  // lu decomposition of a in-place
  ludcmp(a,n,indx,&d);

  // inverse by columns
  for (j=0; j<n; j++) {
    int i;
    for (i=0; i<n; i++) col[i] = 0;
    col[j] = 1;
    lubksb(a,n,indx,col);
    for (i=0; i<n; i++) y[i][j] = col[i];
  }

  // delete the temporary spaces
  delete [] col;
  delete [] indx;
  for (j=0; j<n; j++) 
    delete [] a[j];
  delete [] a;

  return ret;
}

void lubksb(double **a, int n, int *indx, double b[])
{
  int i, ii=0,ip,j;
  double sum;

  for (i=1; i<=n; i++) {
    ip = indx[i-1];
    sum = b[ip-1];
    b[ip-1] = b[i-1];
    if (ii) {
      for (j=ii;j<=i-1;j++) sum -= a[i-1][j-1]*b[j-1];
    }
    else 
      if (sum) ii=i;
    b[i-1] = sum;
  }
  for (i=n; i>=1; i--) {
    sum = b[i-1];
    for (j=i+1; j<=n; j++) sum -= a[i-1][j-1]*b[j-1];
    b[i-1] = sum/a[i-1][i-1];
  }
}

#define TINY 1.0e-20

void ludcmp(double **a, int n, int *indx, double *d)
{
  int i, imax, j, k;
  double big,dum,sum,temp;
  double *vv = new double[n]; 
  for (i=1; i<=n; i++)
    vv[i-1] = 1.0;

  *d = 1.0;
  for (i=1; i<=n; i++) {
    big = 0.0;
    for (j=1; j<=n; j++) 
      if ((temp=fabs(a[i-1][j-1]))>big) big = temp;
    if (big==0.0) {
      printf("Singular matrix in routine ludcmp\n");
      exit(1);
    }
    vv[i-1] = 1.0/big;
  }
  for (j=1; j<=n; j++) {
    for (i=1; i<j; i++) {
      sum = a[i-1][j-1];
      for (k=1; k<i; k++) sum -= a[i-1][k-1]*a[k-1][j-1];
      a[i-1][j-1] = sum;
    }
    big = 0.0;
    for (i=j; i<=n; i++) {
      sum = a[i-1][j-1];
      for (k=1; k<j; k++)
	sum -= a[i-1][k-1]*a[k-1][j-1];
      a[i-1][j-1] = sum;
      if ((dum = vv[i-1]*fabs(sum))>=big) {
	big = dum;
	imax = i;
      }
    }
    if (j != imax) {
      for (k=1; k<=n; k++) {
	dum = a[imax-1][k-1];
	a[imax-1][k-1] = a[j-1][k-1];
	a[j-1][k-1] = dum;
      }
      *d = -(*d);
      vv[imax-1]=vv[j-1];
    }
    indx[j-1] = imax;
    if (a[j-1][j-1] == 0.0) a[j-1][j-1] = TINY;
    if (j!=n) {
      dum = 1.0/(a[j-1][j-1]);
      for (i=j+1; i<=n; i++) a[i-1][j-1] *= dum;
    }
  }
  delete [] vv;

  return;
}

matrix::matrix(const matrix& other)
{
  m = other.m;
  n = other.n;
  mat = new double *[m];
  for (int i=0; i<m; i++) {
    mat[i] = new double[n];
    for (int j=0; j<n; j++)
      mat[i][j] = other.mat[i][j];
  }
}

matrix::~matrix()
{
  for (int i=0; i<m; i++) 
    delete [] mat[i];
  delete [] mat;
}

matrix& matrix::operator = (const matrix& other)
{
  if (&other != this) {
    m = other.m;
    n = other.n;
    mat = new double *[m];
    for (int i=0; i<m; i++) {
      mat[i] = new double[n];
      for (int j=0; j<n; j++)
	mat[i][j] = other.mat[i][j];
    }
  }
}

