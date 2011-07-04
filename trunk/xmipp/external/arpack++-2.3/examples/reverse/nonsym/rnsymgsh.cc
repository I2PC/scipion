/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE RNSymGSh.cc.
   Example program that illustrates how to solve a real
   nonsymmetric generalized eigenvalue problem in real shift
   and invert mode using the ARrcNonSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in shift and
      invert mode, where A and B are derived from the finite element
      discretization of the 1-dimensional convection-diffusion operator
                         (d^2u / dx^2) + rho*(du/dx)
      on the interval [0,1] with zero Dirichlet boundary conditions
      using linear elements.
      The shift sigma is a real number.

   2) Data structure used to represent matrix A:

      ARrcNonSymGenEig is a class that requires the user to provide a
      way to perform the matrix-vector products w = OP*Bv =
      inv(A-sigma*B)*B*v and w = B*v, where sigma is the adopted shift.
      In this example a class called NonSymGenProblemB was created with
      this purpose. NonSymGenProblemB contains a member function,
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
      ngenprbb.h       The NonSymGenProblemB class definition.
      arrgnsym.h       The ARrcNonSymGenEig class definition.
      rnsymgsl.h       The Solution function.

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "ngenprbb.h"
#include "rnsymgsl.h"
#include "arrgnsym.h"


template<class T>
void Test(T type)
{

  // Defining a temporary vector.

  T temp[101];

  // Creating a pencil.

  NonSymGenProblemB<T> P(100, 10, 1.0); // n = 100, rho = 10, sigma = 1.

  // Creating a nonsymmetric eigenvalue problem and defining what we need:
  // the four eigenvectors with largest magnitude.

  ARrcNonSymGenEig<T> prob(P.A.ncols(), 4L, 1.0);

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

  // Solving a double precision problem with n = 100.

  Test((double)0.0);

  // Solving a single precision problem with n = 100.

  Test((float)0.0);

} // main

