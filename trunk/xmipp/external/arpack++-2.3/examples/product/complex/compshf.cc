/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE CompShf.cc.
   Example program that illustrates how to solve a complex
   standard eigenvalue problem in shift and invert mode 
   using the ARCompStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in shift and
      invert mode, where A is derived from the central difference 
      discretization of the 1-dimensional convection-diffusion operator
                         (d^2u/dx^2) + rho*(du/dx)
      on the interval [0,1] with zero Dirichlet boundary conditions.

   2) Data structure used to represent matrix A:

      When using ARCompStdEig, the user is required to provide a
      class that contains a member function which computes the
      matrix-vector product w = OPv, where OP = inv[A - sigma*I].
      In this example, this class is called CompMatrixB, and 
      MultOPv is the function.

   3) Included header files:

      File             Contents
      -----------      -------------------------------------------
      cmatrixb.h       The CompMatrixB class definition.
      arscomp.h        The ARCompStdEig class definition.
      compsol.h        The Solution function.
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
#include "arscomp.h"
#include "cmatrixb.h"
#include "compsol.h"


template<class T>
void Test(T type)
{

  // Creating a complex matrix (n = 100, shift = 0, rho = 10).

  CompMatrixB<T> A(100, arcomplex<T>(0.0,0.0), arcomplex<T>(10.0,0.0));

  // Defining what we need: the four eigenvectors of F nearest to 0.0.
  // A.MultOPv is the function that performs the product w <- OPv.

  ARCompStdEig<T, CompMatrixB<T> >
    dprob(A.ncols(), 4, &A, &CompMatrixB<T>::MultOPv, arcomplex<T>(0.0, 0.0));

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, dprob);

} // Test.


int main()
{

  // Solving a single precision problem with n = 100.

#ifndef __SUNPRO_CC

  Test((float)0.0);

#endif

  // Solving a double precision problem with n = 100.

  Test((double)0.0);

} // main

