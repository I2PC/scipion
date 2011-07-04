/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE DNSymShf.cc.
   Example program that illustrates how to solve a real nonsymmetric
   dense standard eigenvalue problem in shift and invert mode using \
   the ARluNonSymStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in shift and invert
      mode, where A is derived from the standard central difference
      discretization of the 2-dimensional convection-diffusion operator
                       (Laplacian u) + rho*(du/dx)
      on a unit square with zero Dirichlet boundary conditions.
      The shift is a real number.

   2) Data structure used to represent matrix A:

      Although A is very sparse in this example, it is stored
      here columnwise as a dense matrix. 

   3) Library called by this example:

      The LAPACK package is called by ARluNonSymStdEig to solve
      some linear systems involving (A-sigma*I). This is needed to
      implement the shift and invert strategy.

   4) Included header files:

      File             Contents
      -----------      --------------------------------------------
      dnmatrxa.h       DenseMatrixA, a function that generates
                       matrix A.
      ardnsmat.h       The ARdsNonSymMatrix class definition.
      ardsnsym.h       The ARluNonSymStdEig class definition.
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

#include "dnmatrxa.h"
#include "ardnsmat.h"
#include "ardsnsym.h"
#include "lnsymsol.h"


int main()
{

  // Defining variables;

  int     nx;
  int     n;       // Dimension of the problem.
  double  rho;     // Parameter used to define A.
  double* A;       // pointer to an array that stores the elements of A.

  // Creating a 100x100 matrix.

  nx  = 10;
  rho = 10.0;
  DenseMatrixA(nx, rho, n, A);
  ARdsNonSymMatrix<double, double> matrix(n, A);

  // Defining what we need: the four eigenvectors of A nearest to 1000.0.

  ARluNonSymStdEig<double> dprob(4L, matrix, 1000.0, "LM", 20);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(matrix, dprob);

} // main
