/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE BSVD.cc.
   Example program that illustrates how to determine the largest 
   singular values of a matrix using arpack++.

   1) Problem description:

      In this example, Arpack++ is called to solve the symmetric problem:

                             (A'*A)*v = sigma*v

      where A is an n by n real band matrix.
      This formulation is appropriate when m >= n.
      The roles of A and A' must be reversed in the case that m < n.

   2) Data structure used to represent the matrix:

      {ndiagL, ndiagU, A}: matrix A data in band format. The columns
      of A are stored sequentially in vector A. ndiagL and ndiagU
      supply the lower and upper bandwidth of A, respectively. 

   3) Included header files:

      File             Contents
      -----------      --------------------------------------------
      bnmatrxw.h       MatrixW, a function that generates matrix A
                       in band format.
      arbnsmat.h       The ARbdNonSymMatrix class definition.
      arssym.h         The ARSymStdEig class definition.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/
#include <cmath>
#include <iostream>

#include "arssym.h"
#include "bnmatrxw.h"
#include "arbnsmat.h"


int main()
{

  // Defining variables;

  int     i;
  int     n;          // Number of columns in A.
  int     nl;         // Lower bandwidth of A.
  int     nu;         // Upper bandwidth of A.
  double* valA;       // Pointer to an array that stores the elements of A.
  double* svalue = new double[4];

  // Creating a band matrix with n = 100.

  n  = 100;
  nl = 6;
  nu = 3;
  MatrixW(n, nl, nu, valA);

  // Using ARluNonSymMatrix to store matrix information and to
  // perform the product A'Ax (LU decomposition is not used).

  ARbdNonSymMatrix<double, double> A(n, nl, nu, valA);

  // Defining what we need: eigenvalues with largest magnitude.

  ARSymStdEig<double, ARbdNonSymMatrix<double, double> >
    dprob(n, 4L, &A, &ARbdNonSymMatrix<double, double>::MultMtMv);

  // Finding eigenvalues.

  dprob.Eigenvalues(svalue);

  // Calculating singular values.

  for (i=0; i<4; i++) {
    svalue[i] = sqrt(svalue[i]);
  }

  // Printing some information about the problem.

  std::cout << std::endl << "Testing ARPACK++ class ARSymStdEig" << std::endl;
  std::cout << "Obtaining singular values by solving (A'*A)*v = sigma*v" << std::endl;
  std::cout << std::endl << "greatest singular values: " << std::endl; 
  for (i=0; i<4; i++) {
    std::cout << "  sigma [" << i+1 << "]: " << svalue[i] << std::endl;
  }

} // main.

