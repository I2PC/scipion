/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE BNSymGSC.cc.
   Example program that illustrates how to solve a real nonsymmetric
   generalized eigenvalue problem in complex shift and invert mode 
   using the ARluNonSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in shift and
      inverse mode, where A is a block tridiagonal matrix and B is
      a tridiagonal matrix. Each diagonal block of A is a tridiagonal
      matrix with 4 on the main diagonal, -1-rho*h/2 on the subdiagonal
      and -1+rho*h/2 on the superdiagonal. Each subdiagonal block of A
      is an identity matrix. B has 4 on the main diagonal and 1 on
      the other two diagonals.
      The shift sigma is a complex number.

   2) Data structure used to represent matrices A and B:

      {ndiagL, ndiagU, A}: matrix A data in band format. The columns
      of A are stored sequentially in vector A. ndiagL and ndiagU
      supply the lower and upper bandwidth of A, respectively. 
      {ndiagL, ndiagU, B}: matrix B in band format.

   3) Library called by this example:

      The LAPACK package is called by ARluNonSymGenEig to solve
      some linear systems involving (A-sigma*B). This is needed to
      implement the shift and invert strategy.

   4) Included header files:

      File             Contents
      -----------      -------------------------------------------
      bnmatrxd.h       BandMatrixD, a function that generates
                       matrix A in band format.
      bnmatrxe.h       BandMatrixE, a function that generates 
                       matrix B in band format.
      arbnsmat.h       The ARbdNonSymMatrix class definition.
      arbgnsym.h       The ARluNonSymGenEig class definition.
      lnsymsol.h       The Solution function.

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "bnmatrxd.h"
#include "bnmatrxe.h"
#include "arbnsmat.h"
#include "arbgnsym.h"
#include "lnsymsol.h"


int main()
{

  // Defining variables;

  int     nx;
  int     n;       // Dimension of the problem.
  int     ndiagL;  // Lower bandwidth of A and B.
  int     ndiagU;  // Upper bandwidth of A and B.
  double  rho;     // Parameter used to define A.
  double* valA;    // pointer to an array that stores the elements of A.
  double* valB;    // pointer to an array that stores the elements of B.

  // Creating matrices A and B.

  nx  = 10;
  rho = 100.0;
  BandMatrixD(nx, rho, n, ndiagL, ndiagU, valA);
  ARbdNonSymMatrix<double, double> A(n, ndiagL, ndiagU, valA);

  BandMatrixE(n, ndiagL, ndiagU, valB);
  ARbdNonSymMatrix<double, double> B(n, ndiagL, ndiagU, valB);

  // Defining what we need: the four eigenvectors nearest to 0.4+0.6I.

  ARluNonSymGenEig<double> dprob(4L, A, B, 'I', 0.4, 0.6);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, B, dprob);

} // main.

