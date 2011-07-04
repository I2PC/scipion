/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LNMatrxC.h
   Function template for the stiffness matrix obtained from
   the finite element discretization of the 1-dimensional
   convection diffusion operator d^2u/dx^2 + rho*(du/dx) on
   the interval [0,1] with zero Dirichlet boundary condition
   using linear elements.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LNMATRXC_H
#define LNMATRXC_H

template<class ARFLOAT, class ARINT>
void StiffnessMatrix(ARINT n, ARFLOAT rho, ARINT& nnz,
                     ARFLOAT* &A, ARINT* &irow, ARINT* &pcol)
{

  ARINT   i, j;
  ARFLOAT dd, dl, du, s, h;

  // Defining constants.

  const ARFLOAT one = 1.0;
  const ARFLOAT two = 2.0;

  h  = one/ARFLOAT(n+1);
  s  = rho/two;
  dd = two/h;
  dl = -one/h - s;
  du = -one/h + s;

  // Defining the number of nonzero matrix elements.

  nnz = 3*n-2;

  // Creating output vectors.

  A    = new ARFLOAT[nnz];
  irow = new ARINT[nnz];
  pcol = new ARINT[n+1];

  // Filling A, irow and pcol.

  pcol[0] = 0;
  j = 0;
  for (i=0; i!=n; i++) {
    if (i != 0) {
      irow[j] = i-1;
      A[j++]  = du;
    }
    irow[j] = i;
    A[j++]  = dd;
    if (i != (n-1)) {
      irow[j] = i+1;
      A[j++]  = dl;
    }
    pcol[i+1] = j;
  }

} //  StiffnessMatrix.

#endif // LNMATRXC_H

