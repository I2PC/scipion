/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE BNMatrxB.h
   Function template for the stiffness matrix obtained from
   the finite element discretization of the 1-dimensional
   convection diffusion operator d^2u/dx^2 + rho*(du/dx) on
   the interval [0,1] with zero Dirichlet boundary conditions
   using linear elements.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef BNMATRXB_H
#define BNMATRXB_H

template<class ARFLOAT, class ARINT>
void StiffnessMatrix(ARINT n, ARFLOAT rho, ARINT& nL, ARINT& nU, ARFLOAT* &A)
{

  ARINT   i;
  ARFLOAT dd, dl, du, h;

  // Defining constants.

  h  = 1.0/ARFLOAT(n+1);
  dd = 2.0/h;
  dl = -1.0/h - 0.5*rho;
  du = -1.0/h + 0.5*rho;

  // Defining the lower and the upper bandwidth.

  nL = 1;
  nU = 1;

  // Creating output vector A.

  A  = new ARFLOAT[3*n];
  for (i=1; i<(3*n); i+=3) {
    if (i-1)   A[i-1] = du;
               A[i]   = dd;
    if (n-i-1) A[i+1] = dl;
  }

} //  StiffnessMatrix.

#endif // BNMATRXB_H

