/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE DNSymGSC.cc.
   Example program that illustrates how to solve a real nonsymmetric
   dense generalized eigenvalue problem in complex shift and invert 
   mode using the ARluNonSymGenEig class.

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

      Although A and B are very sparse in this example, they are 
      stored here as dense matrices (by columns).

   3) Library called by this example:

      The LAPACK package is called by ARluNonSymGenEig to solve
      some linear systems involving (A-sigma*B). This is needed to
      implement the shift and invert strategy.

   4) Included header files:

      File             Contents
      -----------      -------------------------------------------
      dnmatrxd.h       DenseMatrixD, a function that generates
                       matrix A.
      dnmatrxe.h       DenseMatrixE, a function that generates 
                       matrix B.
      ardnsmat.h       The ARdsNonSymMatrix class definition.
      ardgnsym.h       The ARluNonSymGenEig class definition.
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

#include "dnmatrxd.h"
#include "dnmatrxe.h"
#include "ardnsmat.h"
#include "ardgnsym.h"
#include "lnsymsol.h"


int main()
{

  // Defining variables;

  int     nx;
  int     n;       // Dimension of the problem.
  double  rho;     // Parameter used to define A.
  double* valA;    // pointer to an array that stores the elements of A.
  double* valB;    // pointer to an array that stores the elements of B.

  // Creating matrices A and B.

  nx  = 10;
  rho = 100.0;
  DenseMatrixD(nx, rho, n, valA);
  ARdsNonSymMatrix<double, double> A(n, valA);

  DenseMatrixE(n, valB);
  ARdsNonSymMatrix<double, double> B(n, valB);

  // Defining what we need: the four eigenvectors nearest to 0.4+0.6I.

  ARluNonSymGenEig<double> dprob(4L, A, B, 'I', 0.4, 0.6);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, B, dprob);

} // main.

