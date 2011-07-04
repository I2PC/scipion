/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE DNSymGSh.cc.
   Example program that illustrates how to solve a real nonsymmetric
   dense generalized eigenvalue problem in real shift and invert mode 
   using the ARluNonSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in real shift
      and invert mode, where A and B are derived from the finite element
      discretization of the 1-dimensional convection-diffusion operator
                        (d^2u / dx^2) + rho*(du/dx)
      on the interval [0,1] with zero Dirichlet boundary conditions
      using linear elements.

   2) Data structure used to represent matrices A and B:

      Although A and B are very sparse in this example, they are 
      stored here as dense matrices (by columns).

   3) Library called by this example:

      The LAPACK package is called by ARluNonSymGenEig to solve
      some linear systems involving (A-sigma*B).

   4) Included header files:

      File             Contents
      -----------      -------------------------------------------
      dnmatrxb.h       StiffnessMatrix, a function that generates
                       matrix A.
      dnmatrxc.h       MassMatrix, a function tha generates 
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

#include "dnmatrxb.h"
#include "dnmatrxc.h"
#include "ardnsmat.h"
#include "ardgnsym.h"
#include "lnsymsol.h"


int main()
{

  // Defining variables;

  int     n;       // Dimension of the problem.
  double  rho;     // Parameter used to define A.
  double* valA;    // pointer to an array that stores the elements of A.
  double* valB;    // pointer to an array that stores the elements of B.

  // Creating matrices A and B.

  n   = 100;
  rho = 10.0;
  StiffnessMatrix(n, rho, valA);
  ARdsNonSymMatrix<double, double> A(n, valA);

  MassMatrix(n, valB);
  ARdsNonSymMatrix<double, double> B(n, valB);

  // Defining what we need: the four eigenvectors nearest to 0.0.

  ARluNonSymGenEig<double> dprob(4L, A, B, 0.0);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, B, dprob);

} // main.

