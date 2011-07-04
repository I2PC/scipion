/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE RSymGCay.cc.
   Example program that illustrates how to solve a real symmetric
   generalized eigenvalue problem in Cayley mode using the
   ARrcSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in Cayley mode,
      where A and B are obtained from the finite element
      discretrization of the 1-dimensional discrete Laplacian
                                d^2u / dx^2
      on the interval [0,1] with zero Dirichlet boundary conditions
      using piecewise linear elements.

   2) Data structure used to represent matrices A and B:

      ARrcSymGenEig is a class that requires the user to provide a
      way to perform the matrix-vector products w = OPv =
      inv(A-sigma*B)*v, w = A*v and w = B*v, where sigma is the
      adopted shift. In this example a class called SymGenProblemB was
      created with this purpose. SymGenProblemB contains a member
      function, MultOPv(v,w), that takes a vector v and returns the
      product OPv in w. It also contains two objects, A and B, that
      store matrices A and B, respectively. The product Bv is
      performed by MultMv, a member function of B, and Av is obtained
      by calling A.MultMv.

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
      sgenprbb.h       The SymGenProblemB class definition.
      arrgsym.h        The ARrcSymGenEig class definition.
      rsymgsol.h       The Solution function.
      blas1c.h         Some blas1 functions.

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "blas1c.h"
#include "sgenprbb.h"
#include "rsymgsol.h"
#include "arrgsym.h"

template<class T>
void Test(T type)
{

  // Defining two temporary vectors.

  T tempA[101], tempB[101];

  // Creating a pencil.

  SymGenProblemB<T> P(100, 150.0);

  // Creating a symmetric eigenvalue problem. 'S' indicates that
  // we will use the shift and invert mode. P.A.ncols() furnishes
  // the dimension of the problem. 4 is the number of eigenvalues
  // sought and 150.0 is the shift.

  ARrcSymGenEig<T> prob('C', P.A.ncols(), 4L, 150.0);

  // Finding an Arnoldi basis.

  while (!prob.ArnoldiBasisFound()) {

    // Calling ARPACK FORTRAN code. Almost all work needed to
    // find an Arnoldi basis is performed by TakeStep.

    prob.TakeStep();

    switch (prob.GetIdo()) {
    case -1:

      // Performing w <- OP*(A+sigma*B)*v for the first time.
      // This product must be performed only if GetIdo is equal to
      // -1. GetVector supplies a pointer to the input vector, v,
      // and PutVector a pointer to the output vector, w.

      P.A.MultMv(prob.GetVector(), tempA);
      P.B.MultMv(prob.GetVector(), tempB);
      axpy(P.A.ncols(), prob.GetShift(), tempB, 1, tempA, 1);
      P.MultOPv(tempA, prob.PutVector());
      break;

    case  1:

      // Performing w <- OP*(A+sigma*B)*v when Bv is available.
      // This product must be performed whenever GetIdo is equal to
      // 1. GetProd supplies a pointer to the previously calculated
      // product Bv, GetVector a pointer to the input vector, v,
      // and PutVector a pointer to the output vector, w.

      P.A.MultMv(prob.GetVector(), tempA);
      axpy(P.A.ncols(), prob.GetShift(), prob.GetProd(), 1, tempA, 1);
      P.MultOPv(tempA, prob.PutVector());
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

