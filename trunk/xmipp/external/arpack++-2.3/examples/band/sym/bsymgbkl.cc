/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE BSymGBkl.cc.
   Example program that illustrates how to solve a real symmetric
   generalized eigenvalue problem in buckling mode using the
   ARluSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in buckling
      mode, where A is the one dimensional discrete Laplacian on
      the interval [0, 1], with zero Dirichlet boundary conditions, 
      and B is the mass matrix formed by using piecewise linear 
      elements on [0, 1].

   2) Data structure used to represent matrices A and B:

      {nsdiag, A}: matrix A data in symmetric band format. The elements 
      of the main diagonal and the first nsdiag subdiagonals of A are 
      stored sequentially, by columns, in vector A.  
      {nsdiag, B}: matrix B in band format.

   3) Library called by this example:

      The LAPACK package is called by ARluSymGenEig to solve
      some linear systems involving (A-sigma*B).

   4) Included header files:

      File             Contents
      -----------      -------------------------------------------
      bsmatrxb.h       BandMatrixB, a function that generates
                       matrix A in band format.
      bsmatrxc.h       BandMatrixC, a function that generates 
                       matrix B in band format.
      arbsmat.h        The ARbdSymMatrix class definition.
      arbgsym.h        The ARluSymGenEig class definition.
      lsymsol.h        The Solution function.

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "bsmatrxb.h"
#include "bsmatrxc.h"
#include "arbsmat.h"
#include "arbgsym.h"
#include "lsymsol.h"


int main()
{

  // Defining variables;

  int     n;       // Dimension of the problem.
  int     nsdiag;  // Lower (and upper) bandwidth of A (and B).
  double* valA;    // pointer to an array that stores the elements of A.
  double* valB;    // pointer to an array that stores the elements of B.

  // Creating matrices A and B.

  n   = 100;
  BandMatrixB(n, nsdiag, valA);
  ARbdSymMatrix<double> A(n, nsdiag, valA);

  BandMatrixC(n, nsdiag, valB);
  ARbdSymMatrix<double> B(n, nsdiag, valB);

  // Defining what we need: the four eigenvectors nearest to 1.0.

  ARluSymGenEig<double> dprob('B', 4L, A, B, 1.0);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, B, dprob);

} // main.

