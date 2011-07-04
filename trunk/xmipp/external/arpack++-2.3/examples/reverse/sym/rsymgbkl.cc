/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE RSymGBkl.cc.
   Example program that illustrates how to solve a real symmetric
   generalized eigenvalue problem in buckling mode using the
   ARrcSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in buckling
      mode, where A and B are obtained from the finite element
      discretrization of the 1-dimensional discrete Laplacian
                                d^2u / dx^2
      on the interval [0,1] with zero Dirichlet boundary conditions
      using piecewise linear elements.

   2) Data structure used to represent matrices A and B:

      ARrcSymGenEig is a class that requires the user to provide a
      way to perform the matrix-vector products w = OP*Av =
      inv(A-sigma*B)*A*v and w = A*v, where sigma is the adopted shift.
      In this example a class called SymGenProblemB was created with
      this purpose. SymGenProblemB contains a member function,
      MultOPv(v,w), that takes a vector v and returns the product OPv
      in w. It also contains an object, A, that stores matrix A data.
      The product Av is performed by MultMv, a member function of A.

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

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "sgenprbb.h"
#include "rsymgsol.h"
#include "arrgsym.h"


template<class T>
void Test(T type)
{

  // Defining a temporary vector.

  T temp[101];

  // Creating a pencil.

  SymGenProblemB<T> P(100, 1.0);

  // Creating a symmetric eigenvalue problem. 'B' indicates that
  // we will use the buckling mode. P.A.ncols() furnishes the
  // dimension of the problem. 4 is the number of eigenvalues
  // sought and 1.0 is the shift.

  ARrcSymGenEig<T> prob('B', P.A.ncols(), 4L, 1.0);

  // Finding an Arnoldi basis.

  while (!prob.ArnoldiBasisFound()) {

    // Calling ARPACK FORTRAN code. Almost all work needed to
    // find an Arnoldi basis is performed by TakeStep.

    prob.TakeStep();

    switch (prob.GetIdo()) {
    case -1:

      // Performing w <- OP*A*v for the first time.
      // This product must be performed only if GetIdo is equal to
      // -1. GetVector supplies a pointer to the input vector, v,
      // and PutVector a pointer to the output vector, w.

      P.A.MultMv(prob.GetVector(), temp);
      P.MultOPv(temp, prob.PutVector());
      break;

    case  1:

      // Performing w <- OP*A*v.
      // This product must be performed whenever GetIdo is equal to
      // 1. GetProd supplies a pointer to the previously calculated
      // product Av and PutVector a pointer to the output vector w.

      P.MultOPv(prob.GetProd(), prob.PutVector());
      break;

    case  2:

      // Performing w <- A*v.

      P.A.MultMv(prob.GetVector(), prob.PutVector());

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

