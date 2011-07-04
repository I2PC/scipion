/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE CompReg.cc.
   Example program that illustrates how to solve a complex
   standard eigenvalue problem in regular mode using the 
   ARCompStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in regular mode,
      where A is obtained from the standard central difference
      discretization of the convection-diffusion operator 
                     (Laplacian u) + rho*(du / dx)
      on the unit squre [0,1]x[0,1] with zero Dirichlet boundary
      conditions.

   2) Data structure used to represent matrix A:

      When using ARCompStdEig, the user is required to provide a
      class that contains a member function which computes the
      the matrix-vector product w = Av. In this example, this class 
      is called CompMatrixA, and MultMv is the function.

   3) Included header files:

      File             Contents
      -----------      -------------------------------------------
      cmatrixa.h       The CompMatrixA class definition.
      arscomp.h        The ARCompStdEig class definition.
      compsol.h        The Solution function.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "arscomp.h"
#include "cmatrixa.h"
#include "compsol.h"


template<class T>
void Test(T type)
{

  // Creating a complex matrix.

  CompMatrixA<T> A(10); // n = 10*10.

  // Defining what we need: the four eigenvectors of A with largest magnitude.
  // A.MultMv is the function that performs the product w <- A.v.

  ARCompStdEig<T, CompMatrixA<T> >
  dprob(A.ncols(), 4, &A, &CompMatrixA<T>::MultMv);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, dprob);

} // Test.


int main()
{

  // Solving a double precision problem with n = 100.

  Test((double)0.0);

  // Solving a single precision problem with n = 100.

#ifndef __SUNPRO_CC

  Test((float)0.0);

#endif

} // main


