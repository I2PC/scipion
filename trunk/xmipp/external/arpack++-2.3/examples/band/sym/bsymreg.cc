/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE BSymReg.cc.
   Example program that illustrates how to solve a real
   symmetric banded standard eigenvalue problem in regular 
   mode using the ARluSymStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in regular 
      mode, where A is derived from the standard central difference
      discretization of the 2-dimensional Laplacian on the unit 
      square with zero Dirichlet boundary conditions.

   2) Data structure used to represent matrix A:

      {nsdiag, A}: matrix A data in symmetric band format. The elements 
      of the main diagonal and the first nsdiag subdiagonals of A are 
      stored sequentially, by columns, in vector A.  

   3) Included header files:

      File             Contents
      -----------      -------------------------------------------
      bsmatrxa.h       BandMatrixA, a function that generates
                       matrix A in band format.
      arbsmat.h        The ARbdSymMatrix class definition.
      arbssym.h        The ARluSymStdEig class definition.
      lsymsol.h        The Solution function.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "bsmatrxa.h"
#include "arbsmat.h"
#include "arbssym.h"
#include "lsymsol.h"


int main()
{

  // Defining variables;

  int     nx;
  int     n;       // Dimension of the problem.
  int     nsdiag;  // Lower (and upper) bandwidth of A.
  double* A;       // Pointer to an array that stores the elements of A.

  // Creating a 100x100 matrix.

  nx  = 10;
  BandMatrixA(nx, n, nsdiag, A);
  ARbdSymMatrix<double> matrix(n, nsdiag, A);

  // Defining what we need: the four eigenvectors of A with largest magnitude.

  ARluSymStdEig<double> dprob(4, matrix);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(matrix, dprob);

} // main.

