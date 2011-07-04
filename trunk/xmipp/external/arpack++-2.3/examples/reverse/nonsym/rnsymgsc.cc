/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE RNSymGSC.cc.
   Example program that illustrates how to solve a real
   nonsymmetric generalized eigenvalue problem in complex shift
   and invert mode using the ARrcNonSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in shift and
      invert mode, where A is the tridiagonal matrix with 2 on the
      diagonal, -2 on the subdiagonal and 3 on the superdiagonal, and
      B is the tridiagonal matrix with 4 on the diagonal and 1 on the
      off-diagonals.
      The shift sigma is a complex number.

   2) Data structure used to represent matrix A:

      To obtain the eigenvalues of the above problem, the user is
      required to provide a way to perform the matrix-vector products
      w = OP*Bv = real{inv(A-sigma*B)}*B*v, w = A*v and w = B*v. In
      this example, a class called NonSymGenProblemC was created with
      this purpose. NonSymGenProblemC contains a member function,
      MultOPv(v,w), that takes a vector v and returns the product OPv
      in w. It also contains two objects, A and B, that
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
      ngenprbc.h       The NonSymGenProblemC class definition.
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

#include "ngenprbc.h"
#include "rnsymgsl.h"
#include "arrgnsym.h"


template<class ARFLOAT>
void RecoverEigenvalues(long nconv, NonSymGenProblemC<ARFLOAT>& P,
                        ARFLOAT EigVec[], ARFLOAT EigValR[], ARFLOAT EigValI[])
/*
   This function is used to recover the eigenvalues of the original
   problem A.x = B.x.lambda, after calling ARPACK++.
   When a complex shift is used to solve a nonsymmetric problem
   defined by the ARRCNonSymGenEig class, the Rayleigh quotient
   lambda = x'Ax/x'Bx must be formed by the user to obtain the
   desired eigenvalues.
   The Rayleigh quotient cannot be calculated automatically by
   ARPACK++ because ARRCNonSymGenEig do not handle matrix information.
   Other classes such as ARNonSymGenEig and ARLUNonSymGenEig do
   not require the user to define this eigenvalue transformation.
*/

{

  int      j, n, ColJ, ColJp1;
  ARFLOAT  numr, numi, denr, deni;
  ARFLOAT* Ax;

  n  = P.A.ncols();
  Ax = new ARFLOAT[n];

  for (j=0; j<nconv; j++) {

    ColJ   = j*n;
    ColJp1 = ColJ+n;

    if (EigValI[j] == (ARFLOAT)0.0) {

      // Eigenvalue is real. Computing EigVal = x'(Ax)/x'(Mx).

      P.A.MultMv(&EigVec[ColJ], Ax);
      numr = dot(n, &EigVec[ColJ], 1, Ax, 1);
      P.B.MultMv(&EigVec[ColJ], Ax);
      denr = dot(n, &EigVec[ColJ], 1, Ax, 1);
      EigValR[j] =  numr / denr;

    }
    else {

      // Eigenvalue is complex.

      // Computing x'(Ax).

      P.A.MultMv(&EigVec[ColJ], Ax);
      numr = dot(n, &EigVec[ColJ], 1, Ax, 1);
      numi = dot(n, &EigVec[ColJp1], 1, Ax, 1);
      P.A.MultMv(&EigVec[ColJp1], Ax);
      numr = numr + dot(n, &EigVec[ColJp1], 1, Ax, 1);
      numi = -numi + dot(n, &EigVec[ColJ], 1, Ax, 1);

      // Computing x'(Mx).

      P.B.MultMv(&EigVec[ColJ], Ax);
      denr = dot(n, &EigVec[ColJ], 1, Ax, 1);
      deni = dot(n, &EigVec[ColJp1], 1, Ax, 1);
      P.B.MultMv(&EigVec[ColJp1], Ax);
      denr = denr + dot(n, &EigVec[ColJp1], 1, Ax, 1);
      deni = -deni + dot(n, &EigVec[ColJ], 1, Ax, 1);

      // Computing the first eigenvalue of the conjugate pair.

      EigValR[j] = (numr*denr+numi*deni) / lapy2(denr, deni);
      EigValI[j] = (numi*denr-numr*deni) / lapy2(denr, deni);

      // Getting the second eigenvalue of the conjugate pair by taking
      // the conjugate of the first.

      EigValR[j+1] = EigValR[j];
      EigValI[j+1] = -EigValI[j];
      j++;

    }

  }

  delete[] Ax;

} // RecoverEigenvalues.


template<class T>
void Test(T type)
{

  long nconv;

  // Defining a temporary vector.

  T temp[101];

  // Creating a pencil.

  NonSymGenProblemC<T> P(100, 0.4, 0.6); // n = 100, sigma = (0.4, 0.6).

  // Creating a nonsymmetric eigenvalue problem. 'R' indicates that
  // we will use the real part of OPv. P.A.ncols() furnishes
  // the dimension of the problem. 4 is the number of eigenvalues
  // sought and 0.4 + 0.6I is the shift.

  ARrcNonSymGenEig<T> prob(P.A.ncols(), 4L, 'R', 0.4, 0.6);

  // Finding an Arnoldi basis.

  while (!prob.ArnoldiBasisFound()) {

    // Calling ARPACK FORTRAN code. Almost all work needed to
    // find an Arnoldi basis is performed by TakeStep.

    prob.TakeStep();

    switch (prob.GetIdo()) {
    case -1:

      // Performing w <- Re{OP*B*v} for the first time.
      // This product must be performed only if GetIdo is equal to
      // -1. GetVector supplies a pointer to the input vector, v,
      // and PutVector a pointer to the output vector, w.

      P.B.MultMv(prob.GetVector(), temp);
      P.MultOPvRe(temp, prob.PutVector());
      break;

    case  1:

      // Performing w <- Real{OP*B*v} when Bv is available.
      // This product must be performed whenever GetIdo is equal to
      // 1. GetProd supplies a pointer to the previously calculated
      // product Bv and PutVector a pointer to the output vector w.

      P.MultOPvRe(prob.GetProd(), prob.PutVector());
      break;

    case  2:

      // Performing w <- B*v.

      P.B.MultMv(prob.GetVector(), prob.PutVector());

    }
  }

  // Finding eigenvalues and eigenvectors.

  nconv = prob.FindEigenvectors();

  // Recovering eigenvalues of the original problem
  // using the Rayleigh quotient.

  RecoverEigenvalues(nconv, P, prob.RawEigenvectors(),
                     prob.RawEigenvalues(), prob.RawEigenvaluesImag());

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

