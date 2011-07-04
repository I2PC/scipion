/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE CompGSol.h
   Template function that exemplify how to print information 
   about complex generalized eigenvalue problems.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef COMPGSOL_H
#define COMPGSOL_H

#include "arcomp.h"
#include <math.h>
#include "blas1c.h"
#include "lapackc.h"
#include "matprod.h"
#include "argcomp.h"

template<class MATRA, class MATRB, class ARFOP, class ARFB, class ARFLOAT>
void Solution(MATRA &A, MATRB &B, ARCompGenEig<ARFLOAT, ARFOP, ARFB> &Prob)
/*
  Prints eigenvalues and eigenvectors of complex generalized
  eigen-problems on standard "cout" stream.
*/

{

  int                i, n, nconv, mode;
  ARFLOAT            *ResNorm;
  arcomplex<ARFLOAT> *Ax, *Bx;


  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  std::cout << std::endl << std::endl << "Testing ARPACK++ class ARCompGenEig" << std::endl;
  std::cout << "Complex generalized eigenvalue problem: A*x - lambda*B*x" << std::endl;
  switch (mode) {
  case 2:
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

    // Printing the residual norm || A*x - lambda*B*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new arcomplex<ARFLOAT>[n];
    Bx      = new arcomplex<ARFLOAT>[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      A.MultMv(Prob.RawEigenvector(i),Ax);
      B.MultMv(Prob.RawEigenvector(i),Bx);
      axpy(n, -Prob.Eigenvalue(i), Bx, 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/
                   lapy2(real(Prob.Eigenvalue(i)),imag(Prob.Eigenvalue(i)));
    }

    for (i=0; i<nconv; i++) {
      std::cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      std::cout << ")*B*x(" << (i+1) << ")||: " << ResNorm[i] << "\n";
    }
    std::cout << std::endl;

    delete[] Ax;
    delete[] Bx;
    delete[] ResNorm;

  }

} // Solution


#endif // COMPGSOL_H
