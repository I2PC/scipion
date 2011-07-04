/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE CompGReg.cc.
   Example program that illustrates how to solve a complex
   generalized eigenvalue problem in regular mode using 
   the ARCompGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in regular
      mode, where A and B are derived from the finite element 
      discretization of the 1-dimensional convection-diffusion operator
                          (d^2u/dx^2) + rho*(du/dx)
      on the interval [0,1] with zero boundary conditions using
      piecewise linear elements. 


   2) Data structure used to represent matrix A:

      When using ARCompGenEig, the user is required to provide some
      classes that contain as member functions the matrix-vector
      products w = OPv = inv(B)*A*v and w = B*v. 
      In this example, ComplexGenProblemA is a class that contains a 
      member function, MultOPv(v,w), that takes a vector v and returns 
      the product OPv in w. It also contains an object, B, that stores
      matrix B data. The product Bv is performed by MultMv, a member
      function of B.

   3) Included header files:

      File             Contents
      -----------      -------------------------------------------
      cgenprba.h       The ComplexGenProblemA class definition.
      argcomp.h        The ARCompGenEig class definition.
      compgsl.h        The Solution function.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "cgenprba.h"
#include "compgsol.h"
#include "argcomp.h"


template<class T>
void Test(T type)
{

  // Creating a complex problem with n = 100.

  ComplexGenProblemA<T> P(100); // n = 100.

  // Defining what we need: the four eigenvectors with largest magnitude.
  // P.MultOPv is the function that performs the product w <- OPv.
  // P.B.MultMv is the function that performs the product w <- Bv.

  ARCompGenEig<T, ComplexGenProblemA<T>, ComplexMatrixD<T> >
    dprob(P.A.ncols(), 4, &P, &ComplexGenProblemA<T>::MultOPv,
          &P.B, &ComplexMatrixD<T>::MultMv);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(P.A, P.B, dprob);

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


