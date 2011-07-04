/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LSymSol.h
   Template functions that exemplify how to print information
   about symmetric standard and generalized eigenvalue problems.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LSYMSOL_H
#define LSYMSOL_H

#include <math.h>
#include "blas1c.h"
#include "lapackc.h"
#ifdef ARLSMAT_H
#include "arlssym.h"
#include "arlgsym.h"
#elif defined ARUSMAT_H
#include "arussym.h"
#include "arugsym.h"
#elif defined ARDSMAT_H
#include "ardssym.h"
#include "ardgsym.h"
#else
#include "arbssym.h"
#include "arbgsym.h"
#endif


template<class ARMATRIX, class ARFLOAT>
void Solution(ARMATRIX &A, ARluSymStdEig<ARFLOAT> &Prob)
/*
  Prints eigenvalues and eigenvectors of symmetric eigen-problems
  on standard "std::cout" stream.
*/

{

  int     i, n, nconv, mode;
  ARFLOAT *Ax;
  ARFLOAT *ResNorm;

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  std::cout << std::endl << std::endl << "Testing ARPACK++ class ARluSymStdEig \n";
  std::cout << "Real symmetric eigenvalue problem: A*x - lambda*x" << std::endl;
  switch (mode) {
  case 1:
    std::cout << "Regular mode" << std::endl;
    break;
  case 3:
    std::cout << "Shift and invert mode" << std::endl;
  }
  std::cout << std::endl;

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

    Ax      = new ARFLOAT[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      A.MultMv(Prob.RawEigenvector(i), Ax);
      axpy(n, -Prob.Eigenvalue(i), Prob.RawEigenvector(i), 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/fabs(Prob.Eigenvalue(i));
    }

    for (i=0; i<nconv; i++) {
      std::cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      std::cout << ")*x(" << (i+1) << ")||: " << ResNorm[i] << "\n";
    }
    std::cout << "\n";

    delete[] Ax;
    delete[] ResNorm;

  }

} // Solution


template<class MATRA, class MATRB, class ARFLOAT>
void Solution(MATRA &A, MATRB &B, ARluSymGenEig<ARFLOAT> &Prob)
/*
  Prints eigenvalues and eigenvectors of symmetric generalized
  eigen-problems on standard "std::cout" stream.
*/

{

  int     i, n, nconv, mode;
  ARFLOAT *Ax, *Bx, *ResNorm;

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  std::cout << std::endl << std::endl << "Testing ARPACK++ class ARluSymGenEig \n";
  std::cout << "Real symmetric generalized eigenvalue problem: A*x - lambda*B*x";
  std::cout << std::endl;
  switch (mode) {
  case 2:
    std::cout << "Regular mode" << std::endl;
    break;
  case 3:
    std::cout << "Shift and invert mode" << std::endl;
    break;
  case 4:
    std::cout << "Buckling mode" << std::endl;
    break;
  case 5:
    std::cout << "Cayley mode" << std::endl;
  }
  std::cout << std::endl;

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

    Ax      = new ARFLOAT[n];
    Bx      = new ARFLOAT[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      A.MultMv(Prob.RawEigenvector(i), Ax);
      B.MultMv(Prob.RawEigenvector(i), Bx);
      axpy(n, -Prob.Eigenvalue(i), Bx, 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/fabs(Prob.Eigenvalue(i));
    }

    for (i=0; i<nconv; i++) {
      std::cout << "||A*x(" << i << ") - lambda(" << i;
      std::cout << ")*B*x(" << i << ")||: " << ResNorm[i] << "\n";
    }
    std::cout << "\n";

    delete[] Ax;
    delete[] Bx;
    delete[] ResNorm;

  }

} // Solution


#endif // LSYMSOL_H

