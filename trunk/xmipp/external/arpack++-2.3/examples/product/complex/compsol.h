/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE CompSol.h
   Template function that exemplify how to print information 
   about complex standard eigenvalue problems.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef COMPSOL_H
#define COMPSOL_H

#include "arcomp.h"
#include <math.h>
#include "blas1c.h"
#include "lapackc.h"
#include "matprod.h"
#include "arscomp.h"

template<class ARMATRIX, class ARFLOAT>
void Solution(ARMATRIX &A, ARCompStdEig<ARFLOAT, ARMATRIX> &Prob)
/*
  Prints eigenvalues and eigenvectors of complex eigen-problems
  on standard "cout" stream.
*/

{

  int                i, n, nconv, mode;
  arcomplex<ARFLOAT> *Ax;
  ARFLOAT            *ResNorm;

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  std::cout << std::endl << std::endl << "Testing ARPACK++ class ARCompStdEig \n";
  std::cout << "Complex eigenvalue problem: A*x - lambda*x" << std::endl;
  switch (mode) {
  case 1:
    std::cout << "Regular mode" << std::endl << std::endl;
    break;
  case 3:
    std::cout << "Shift and invert mode" << std::endl << std::endl;
  }

  std::cout << "Dimension of the system            : " << n              << std::endl;
  std::cout << "Number of 'requested' eigenvalues  : " << Prob.GetNev()  << std::endl;
  std::cout << "Number of 'converged' eigenvalues  : " << nconv          << std::endl;
  std::cout << "Number of Arnoldi vectors generated: " << Prob.GetNcv()  << std::endl;
  std::cout << "Number of iterations taken         : " << Prob.GetIter() << std::endl;
  std::cout << std::endl;

  if (Prob.EigenvaluesFound()) {

    // Printing eigenvalues.

    std::cout << "Eigenvalues:" << std::endl;
    for (i=0; i<nconv; i++) {
      std::cout << "  lambda[" << (i+1) << "]: " << Prob.Eigenvalue(i) << std::endl;
    }
    std::cout << std::endl;
  }

  if (Prob.EigenvectorsFound()) {

    // Printing the residual norm || A*x - lambda*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new arcomplex<ARFLOAT>[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      A.MultMv(Prob.RawEigenvector(i),Ax);
      axpy(n, -Prob.Eigenvalue(i), Prob.RawEigenvector(i), 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/
                   lapy2(real(Prob.Eigenvalue(i)),imag(Prob.Eigenvalue(i)));
    }

    for (i=0; i<nconv; i++) {
      std::cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      std::cout << ")*x(" << (i+1) << ")||: " << ResNorm[i] << std::endl;
    }
    std::cout << "\n";

    delete[] Ax;
    delete[] ResNorm;

  }

} // Solution


#endif // COMPSOL_H
