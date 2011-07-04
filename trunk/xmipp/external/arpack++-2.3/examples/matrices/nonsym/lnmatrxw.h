/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LNMatrxW.h
   Class template for the m by n matrix A with entries

          A(i,j) = k*(s(i))*(t(j)-1), if j-nl<= i <= j;
                   k*(t(j))*(s(i)-1), if j < i <= j+nu,

   where s(i) = i/(m+1), t(j) = j/(n+1) and k = 1/(n+1).       

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LNMATRXW_H
#define LNMATRXW_H


template<class ARFLOAT, class ARINT>
void MatrixW(int m, int n, int nl, int nu, int& nnz, 
             ARFLOAT* &A, ARINT* &irow, ARINT* &pcol)
{

  int   i, j, l, u, mn, p;
  ARFLOAT h, k;

  // Defining constants.

  l = nl-m+n;
  l = (l<0)?0:((l>nl)?nl:l);
  u = nu-n+m;
  u = (u<0)?0:((u>nu)?nu:u);
  h = 1.0/ARFLOAT(m+1);
  k = 1.0/ARFLOAT(n+1);

  // Defining the number of nonzero elements and lines of the matrix.

  mn  = (m<n)?m:n;
  nnz = mn*(nl+nu+1)-(l*l+l)/2-(u*u+u)/2;

  // Creating output vectors.

  A    = new ARFLOAT[nnz];
  irow = new ARINT[nnz];
  pcol = new ARINT[n+1];

  // Filling A, irow and pcol.

  pcol[0] = 0;
  p       = 0;

  for (j=0; j!=n; j++) {

    for(i=((j-nu)>0)?(j-nu):0; ((i<=j)&&(i<m)); i++) {
      irow[p] = i;
      A[p++]  = k*(ARFLOAT(i+1)*h-1.0)*(ARFLOAT(j+1)*k);
    }

    for(i=j+1; ((i<=(j+nl))&&(i<m)); i++) {
      irow[p] = i;
      A[p++]  = k*(ARFLOAT(j+1)*k-1.0)*(ARFLOAT(i+1)*h);
    }

    pcol[j+1] = p;

  }

} // MatrixW.

#endif // LNMATRXW_H

