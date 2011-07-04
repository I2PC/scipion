/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE DCompShf.cc.
   Example program that illustrates how to solve a complex dense
   standard eigenvalue problem in shift and invert mode using the
   ARluCompStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in shift and invert
      mode, where A is derived from the central difference discretization
      of the convection-diffusion operator
                    (Laplacian u) + rho*(du / dx)
      on the unit square [0,1]x[0,1] with zero Dirichlet boundary
      conditions.

   2) Data structure used to represent matrix A:

      Although A is very sparse in this example, it is stored
      here columnwise as a dense matrix. 

   3) Library called by this example:

      The LAPACK package is called by ARluCompStdEig to solve
      some linear systems involving (A-sigma*I). This is needed to
      implement the shift and invert strategy.

   4) Included header files:

      File             Contents
      -----------      ---------------------------------------------
      dcmatrxa.h       CompMatrixB, a function that generates 
                       matrix A.
      ardnsmat.h       The ARdsNonSymMatrix class definition.
      ardscomp.h       The ARluCompStdEig class definition.
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
#include "ardnsmat.h"
#include "ardscomp.h"
#include "dcmatrxa.h"
#include "lcompsol.h"


int main()
{

  // Defining variables;

  int                nx;
  int                n;      // Dimension of the problem.
  arcomplex<double>* valA;   // pointer to an array that stores
                             // the elements of A.

  // Creating a complex matrix.

  nx = 10;
  CompMatrixA(nx, n, valA);
  ARdsNonSymMatrix<arcomplex<double>, double> A(n, valA);

  // Defining what we need: the four eigenvectors of F nearest to 0.0.

  ARluCompStdEig<double> dprob(4L, A, arcomplex<double>(0.0, 0.0));

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, dprob);

} // main.

