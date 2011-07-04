/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE NSymGSol.h
   Template functions that exemplify how to print information 
   about nonsymmetric generalized eigenvalue problems.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef NSYMGSOL_H
#define NSYMGSOL_H

#include <math.h>
#include "blas1c.h"
#include "lapackc.h"
#include "matprod.h"
#include "argnsym.h"

template<class MATRA, class MATRB, class ARFOP, class ARFB, class ARFLOAT>
void Solution(MATRA &A, MATRB &B, ARNonSymGenEig<ARFLOAT, ARFOP, ARFB> &Prob)
/*
  Prints eigenvalues and eigenvectors of nonsymmetric generalized
  eigen-problems on standard "cout" stream.
*/

{

  int     i, n, nconv, mode;
  ARFLOAT *Ax;
  ARFLOAT *Bx, *Bx1;
  ARFLOAT *ResNorm;

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  std::cout << std::endl << std::endl << "Testing ARPACK++ class ARNonSymGenEig" << std::endl;
  std::cout << "Real nonsymmetric generalized eigenvalue problem: A*x - lambda*B*x";
  std::cout << std::endl;
  switch (mode) {
  case 2:
    std::cout << "Regular mode" << std::endl;
    break;
  case 3:
    std::cout << "Shift and invert mode (using real part of OP)" << std::endl;
    break;
  case 4:
    std::cout << "Shift and invert mode (using imaginary part of OP)" << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Dimension of the system            : " << n              << std::endl;
  std::cout << "Number of 'requested' eigenvalues  : " << Prob.GetNev()  << std::endl;
  std::cout << "Number of 'converged' eigenvalues  : " << nconv          << std::endl;
  std::cout << "Number of Arnoldi vectors generated: " << Prob.GetNcv()  << std::endl;
  std::cout << "Number of iterations taken         : " << Prob.GetIter() << std::endl;
  std::cout << std::endl;

  if (Prob.EigenvaluesFound()) {

    // Eigenvalues.

    std::cout << "Eigenvalues:" << std::endl;
    for (i=0; i<nconv; i++) {
      std::cout << "  lambda[" << (i+1) << "]: " << Prob.EigenvalueReal(i);
      if (Prob.EigenvalueImag(i)>=0.0) {
        std::cout << " + " << Prob.EigenvalueImag(i) << " I" << std::endl;
      }
      else {
        std::cout << " - " << fabs(Prob.EigenvalueImag(i)) << " I" << std::endl;
      }
    }
    std::cout << std::endl;
  }

  if (Prob.EigenvectorsFound()) {

    // Printing the residual norm || A*x - lambda*B*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new ARFLOAT[n];
    Bx      = new ARFLOAT[n];
    Bx1     = new ARFLOAT[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {

      if (Prob.EigenvalueImag(i)==0.0) { // Eigenvalue is real.

        A.MultMv(Prob.RawEigenvector(i), Ax);
        B.MultMv(Prob.RawEigenvector(i), Bx);
        axpy(n, -Prob.EigenvalueReal(i), Bx, 1, Ax, 1);
        ResNorm[i] = nrm2(n, Ax, 1)/fabs(Prob.EigenvalueReal(i));

      }
      else {                 // Eigenvalue is complex.

        A.MultMv(Prob.RawEigenvector(i), Ax);
        B.MultMv(Prob.RawEigenvector(i), Bx);
        B.MultMv(Prob.RawEigenvector(i+1), Bx1);
        axpy(n, -Prob.EigenvalueReal(i), Bx, 1, Ax, 1);
        axpy(n, Prob.EigenvalueImag(i), Bx1, 1, Ax, 1);
        ResNorm[i] = nrm2(n, Ax, 1);
        A.MultMv(Prob.RawEigenvector(i+1), Ax);
        axpy(n, -Prob.EigenvalueImag(i), Bx, 1, Ax, 1);
        axpy(n, -Prob.EigenvalueReal(i), Bx1, 1, Ax, 1);
        ResNorm[i] = lapy2(ResNorm[i], nrm2(n, Ax, 1))/
                     lapy2(Prob.EigenvalueReal(i), Prob.EigenvalueImag(i));
        ResNorm[i+1] = ResNorm[i];
        i++;

      }
    }

    for (i=0; i<nconv; i++) {
      std::cout << "||A*x(" << i << ") - lambda(" << i;
      std::cout << ")*B*x(" << i << ")||: " << ResNorm[i] << "\n";
    }
    std::cout << "\n";

    delete[] Ax;
    delete[] Bx;
    delete[] Bx1;
    delete[] ResNorm;

  }

} // Solution


#endif // NSYMGSOL_H
