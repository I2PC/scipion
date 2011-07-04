/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE DCompGRe.cc.
   Example program that illustrates how to solve a complex dense
   generalized eigenvalue problem in regular mode using the 
   ARluCompGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in regular mode,
      where A and B are derived from the finite element discretization
      of the 1-dimensional convection-diffusion operator
                       (d^2u/dx^2) + rho*(du/dx)
      on the interval [0,1], with zero boundary conditions, using
      piecewise linear elements.

   2) Data structure used to represent matrices A and B:

      Although A and B are very sparse in this example, they are 
      stored here as dense matrices (by columns).

   3) Library called by this example:

      The LAPACK package is called by ARluCompGenEig to solve
      some linear systems involving B.

   4) Included header files:

      File             Contents
      -----------      ---------------------------------------------
      dcmatrxb.h       CompMatrixE, a function that generates 
                       matrix A.
      dcmatrxc.h       CompMatrixF, a function that generates 
                       matrix B.
      ardnsmat.h       The ARdsNonSymMatrix class definition.
      ardgcomp.h       The ARluCompGenEig class definition.
      lcompsol.h       The Solution function.
      arcomp.h         The "arcomplex" (complex) type definition.

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "arcomp.h"
#include "dcmatrxb.h"
#include "dcmatrxc.h"
#include "ardnsmat.h"
#include "ardgcomp.h"
#include "lcompsol.h"


int main()
{

  // Defining variables;

  int               n;            // Dimension of the problem.
  arcomplex<double> rho;          // Parameter used to define A.
  arcomplex<double> *valA, *valB; // pointers to arrays that store
                                  // the elements of A and B.

  // Creating complex matrices A and B.

  n   =  100;
  rho = arcomplex<double>(10.0, 0.0);
  CompMatrixB(n, rho, valA);
  ARdsNonSymMatrix<arcomplex<double>, double> A(n, valA);

  CompMatrixC(n, valB);
  ARdsNonSymMatrix<arcomplex<double>, double> B(n, valB);

 // Defining what we need: the four eigenvectors with largest magnitude.

  ARluCompGenEig<double> dprob(4L, A, B);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, B, dprob);

} // main.
