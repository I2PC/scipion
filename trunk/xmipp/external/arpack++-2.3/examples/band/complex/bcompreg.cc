/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE BCompReg.cc.
   Example program that illustrates how to solve a complex standard
   eigenvalue problem in regular mode using the ARluCompStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in regular mode,
      where A is obtained from the standard central difference
      discretization of the convection-diffusion operator
                    (Laplacian u) + rho*(du / dx)
      on the unit square [0,1]x[0,1] with zero Dirichlet boundary
      conditions.

   2) Data structure used to represent matrix A:

      {ndiagL, ndiagU, A}: matrix A data in band format. The columns
      of A are stored sequentially in vector A. ndiagL and ndiagU
      supply the lower and upper bandwidth of A, respectively. 

   3) Included header files:

      File             Contents
      -----------      ---------------------------------------------
      bcmatrxa.h       CompMatrixA, a function that generates matrix
                       A in band format.
      arbnsmat.h       The ARbdNonSymMatrix class definition.
      arbscomp.h       The ARluCompStdEig class definition.
      lcompsol.h       The Solution function.
      arcomp.h         The "arcomplex" (complex) type definition.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "arcomp.h"
#include "arbscomp.h"
#include "arbnsmat.h"
#include "bcmatrxa.h"
#include "lcompsol.h"


int main()
{

  // Defining variables;

  int                nx;
  int                n;      // Dimension of the problem.
  int                ndiagL; // Lower bandwidth of A and B.
  int                ndiagU; // Upper bandwidth of A and B.
  arcomplex<double>* valA;   // pointer to an array that stores
                             // the elements of A.

  // Creating a complex matrix.

  nx = 10;
  CompMatrixA(nx, n, ndiagL, ndiagU, valA);
  ARbdNonSymMatrix<arcomplex<double>, double> A(n, ndiagL, ndiagU, valA);

  // Defining what we need: the four eigenvectors of A with largest magnitude.

  ARluCompStdEig<double> dprob(4L, A);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, dprob);

} // main

