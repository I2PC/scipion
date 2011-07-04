/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE DSymReg.cc.
   Example program that illustrates how to solve a real
   symmetric dense standard eigenvalue problem in regular 
   mode using the ARluSymStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in regular 
      mode, where A is derived from the standard central difference
      discretization of the 2-dimensional Laplacian on the unit 
      square with zero Dirichlet boundary conditions.

   2) Data structure used to represent matrix A:

      Although A is very sparse in this example, it is stored
      here as a dense symmetric matrix. The lower triangular part
      of A is stored, by columns, in the vector A.

   3) Included header files:

      File             Contents
      -----------      -------------------------------------------
      dsmatrxa.h       DenseMatrixA, a function that generates
                       matrix A.
      ardsmat.h        The ARdsSymMatrix class definition.
      ardssym.h        The ARluSymStdEig class definition.
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

#include "dsmatrxa.h"
#include "ardsmat.h"
#include "ardssym.h"
#include "lsymsol.h"


int main()
{

  // Defining variables;

  int     nx;
  int     n;   // Dimension of the problem.
  double* A;   // Pointer to an array that stores the lower triangular
               // elements of A.

  // Creating a 100x100 matrix.

  nx  = 10;
  DenseMatrixA(nx, n, A);
  ARdsSymMatrix<double> matrix(n, A);

  // Defining what we need: the four eigenvectors of A with largest magnitude.

  ARluSymStdEig<double> dprob(4, matrix);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(matrix, dprob);

} // main.
