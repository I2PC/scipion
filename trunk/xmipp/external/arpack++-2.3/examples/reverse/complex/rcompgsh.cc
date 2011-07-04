/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE RCompGSh.cc.
   Example program that illustrates how to solve a complex
   generalized eigenvalue problem in shift and invert mode 
   using the ARrcCompGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in shift and
      invert mode, where A and B are derived from a finite element 
      discretization of a 1-dimensional convection-diffusion operator
                         (d^2u/dx^2) + rho*(du/dx)
      on the interval [0,1] with zero boundary conditions using 
      piecewise linear elements.

   2) Data structure used to represent matrix A:

      ARrcCompGenEig is a class thar requires the user to provide a
      way to perform the matrix-vector products w = OP*Bv =
      inv(A-sigma*B)*B*v and w = B*v, where sigma is the adopted shift.
      In this example a class called ComplexGenProblemB was created 
      with this purpose. ComplexGenProblemB contains a member function,
      MultOPv(v,w), that takes a vector v and returns the product OPv
      in w. It also contains an object, B, that stores matrix B data.
      The product Bv is performed by MultMv, a member function of B.

   3) The reverse communication interface:

      This example uses the reverse communication interface, which
      means that the desired eigenvalues cannot be obtained directly
      from an ARPACK++ class.
      Here, the overall process of finding eigenvalues by using the
      Arnoldi method is splitted into two parts. In the first, a
      sequence of calls to a function called TakeStep is combined
      with matrix-vector products in order to find an Arnoldi basis.
      In the second part, an ARPACK++ function like FindEigenvectors
      (or EigenValVectors) is used to extract eigenvalues and
      eigenvectors.

   4) Included header files:

      File             Contents
      -----------      -------------------------------------------
      cgenprbb.h       The ComplexGenProblemB class definition.
      arrgcomp.h       The ARrcCompGenEig class definition.
      rcompgsl.h       The Solution function.
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
#include "cgenprbb.h"
#include "rcompgsl.h"
#include "arrgcomp.h"


template<class T>
void Test(T type)
{

  // Defining a temporary vector.

  arcomplex<T> temp[101];

  // Defining a complex pencil with n = 100, rho = 10, sigma = 1.

  ComplexGenProblemB<T> P(100, arcomplex<T>(10.0,0.0), arcomplex<T>(1.0,0.0));

  // Creating a complex eigenvalue problem and defining what we need:
  // the four eigenvectors nearest to 1.0.

  ARrcCompGenEig<T> prob(P.A.ncols(), 4L, arcomplex<T>(1.0,0.0));

  // Finding an Arnoldi basis.

  while (!prob.ArnoldiBasisFound()) {

    // Calling ARPACK FORTRAN code. Almost all work needed to
    // find an Arnoldi basis is performed by TakeStep.

    prob.TakeStep();

    switch (prob.GetIdo()) {
    case -1:

      // Performing w <- OP*B*v for the first time.
      // This product must be performed only if GetIdo is equal to
      // -1. GetVector supplies a pointer to the input vector, v,
      // and PutVector a pointer to the output vector, w.

      P.B.MultMv(prob.GetVector(), temp);
      P.MultOPv(temp, prob.PutVector());
      break;

    case  1:

      // Performing w <- OP*B*v when Bv is available.
      // This product must be performed whenever GetIdo is equal to
      // 1. GetProd supplies a pointer to the previously calculated
      // product Bv and PutVector a pointer to the output vector w.

      P.MultOPv(prob.GetProd(), prob.PutVector());
      break;

    case  2:

      // Performing w <- B*v.

      P.B.MultMv(prob.GetVector(), prob.PutVector());

    }
  }

  // Finding eigenvalues and eigenvectors.

  prob.FindEigenvectors();

  // Printing solution.

  Solution(prob);

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

