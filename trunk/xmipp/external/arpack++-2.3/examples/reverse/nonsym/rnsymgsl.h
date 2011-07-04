/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE RNSymGSl.h
   Printing eigenvalues of a nonsymmetric generalized problem
   (ARrcNonSymGenEig version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef RNSYMGSL_H
#define RNSYMGSL_H

#include "arrgnsym.h"

template<class ARFLOAT>
void Solution(ARrcNonSymGenEig<ARFLOAT> &Prob)
/*
  Prints eigenvalues on standard "cout" stream.
*/

{

  int i, n, nconv, mode;

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  std::cout << std::endl << std::endl << "Testing ARPACK++ class ARrcNonSymGenEig" << std::endl;
  std::cout << "Real nonsymmetric generalized eigenvalue problem: A*x - B*x*lambda" << std::endl;
  switch (mode) {
  case 2:
    std::cout << "Regular mode";
    break;
  case 3:
    std::cout << "Shift and invert mode (using real part of OP)";
    break;
  case 4:
    std::cout << "Shift and invert mode (using imaginary part of OP)";
  }
  std::cout << std::endl << std::endl;

  std::cout << "Dimension of the system            : " << n              << std::endl;
  std::cout << "Number of 'requested' eigenvalues  : " << Prob.GetNev()  << std::endl;
  std::cout << "Number of 'converged' eigenvalues  : " << nconv          << std::endl;
  std::cout << "Number of Arnoldi vectors generated: " << Prob.GetNcv()  << std::endl;
  std::cout << "Number of iterations taken         : " << Prob.GetIter() << std::endl;
  std::cout << std::endl;

  if (Prob.EigenvaluesFound()) {
    std::cout << "Eigenvalues:" << std::endl;
    for (i=0; i<nconv; i++) {
      std::cout << "  lambda[" << (i+1) << "]: " << Prob.Eigenvalue(i) << std::endl;
    }
    std::cout << std::endl;
  }

} // Solution


#endif // RNSYMGSL_H

