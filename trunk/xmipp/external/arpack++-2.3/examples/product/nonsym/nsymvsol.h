/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE NSymVSol.h
   Template functions that exemplify how to print information 
   about the singular value decomposition obtained using the
   ARNonSymStdEig function.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef NSYMVSOL_H
#define NSYMVSOL_H

#include <math.h>
#include "blas1c.h"
#include "lapackc.h"
#include "matprod.h"
#include "arsnsym.h"

template<class ARMATRIX, class ARFLOAT>
void Solution(ARMATRIX &A, ARNonSymStdEig<ARFLOAT, ARMATRIX> &Prob)
/*
  Prints singular values and singular vectors of nonsymmetric 
  real matrices on standard "cout" stream.
*/

{

  int     i, m, n, nconv;
  ARFLOAT ResNorm, t;
  ARFLOAT *value;
  ARFLOAT *Ax;
  ARFLOAT *u;

  n     = Prob.GetN();
  m     = A.nrows();
  nconv = Prob.ConvergedEigenvalues();

  std::cout << std::endl << std::endl << "Testing ARPACK++ class ARNonSymStdEig \n";
  std::cout << "Singular value decomposition problem: (A'*A)*x - lambda*x" << std::endl;

  std::cout << "Dimension of the system              : " << n             << std::endl;
  std::cout << "Number of 'requested' singular values: " << Prob.GetNev() << std::endl;
  std::cout << "Number of 'converged' singular values: " << nconv         << std::endl;
  std::cout << "Number of Arnoldi vectors generated  : " << Prob.GetNcv() << std::endl;
  std::cout << "Number of iterations taken           : " << Prob.GetIter() << std::endl;
  std::cout << std::endl;

  value = new ARFLOAT[nconv];
  Ax    = new ARFLOAT[m];
  u     = new ARFLOAT[m];

  if (Prob.EigenvaluesFound()) {

    // Printing singular values.

    std::cout << "Singular values:" << std::endl;
    for (i=0; i<nconv; i++) {
      value[i] = sqrt(Prob.EigenvalueReal(i));
      std::cout << "  sigma[" << (i+1) << "]: " << value[i] << std::endl;
    }
    std::cout << std::endl;
  }

  if (Prob.EigenvectorsFound()) {

    // Printing the residual norm || A*v - sigma*u ||,
    // where v is a right and u a left singular vector.

    for (i=0; i<nconv; i++) {

      // Computing the left singular vectors.

      A.MultMv(Prob.RawEigenvector(i), Ax);
      copy(m, Ax, 1, u, 1);
      t = 1.0/nrm2(m, u, 1);
      scal(m, t, u, 1);

      // Computing the residual norm.

      axpy(n, -value[i], u, 1, Ax, 1);
      ResNorm = nrm2(n, Ax, 1);

      // Printing the norm.

      std::cout << "||A*v(" << (i+1) << ") - sigma(" << (i+1);
      std::cout << ")*u(" << (i+1) << ")||: " << ResNorm << std::endl;
    }

  }

  delete[] value;
  delete[] Ax;
  delete[] u;

} // Solution


#endif // NSYMVSOL_H
