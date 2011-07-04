/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE SymSimp.cc.
   Example program that illustrates how to solve a very simple
   real symmetric standard eigenvalue problem in regular mode 
   using the ARSymStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in regular mode
      where A is the one dimensional discrete Laplacian on the interval 
      [0,1] with zero Dirichlet boundary conditions. 

  2) Data structure used to represent matrices A and B:

      When using ARSymGenEig, the user is required to provide a 
      class that contain as member function the matrix-vector
      product w = Av. In this example, this class is called SymMatrix.

   3) Included header files:

      File             Contents
      -----------      -------------------------------------------
      arssym.h         The ARSymStdEig class definition.
      blas1.h          Some blas 1 routines.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "blas1c.h"
#include "arssym.h"

template<class T>
class SymMatrix {
/*

  This simple class exemplifies how to create a matrix class 
  that can be used by ARPACK++.
  Basically, SymMatrix is required to have a member function 
  that calculates the matrix-vector product SymMatrix*v, 
  where v is a vector with elements of type T.

*/

 private:

  int  n;

 public:

  int  ncols() { return n; }

  void MultMv(T* v, T* w)

  /*
    Matrix vector subroutine
    where the matrix is the one dimensional discrete Laplacian on
    the interval [0,1] with zero Dirichlet boundary conditions. 
  */

  {
  
    int  j;
    T    h2;

    const T two = 2.0;

    w[0] =  two*v[0] - v[1];
    for (j=1; j<n-1; j++) {
      w[j] = - v[j-1] + two*v[j] - v[j+1];
    }
    w[n-1] = - v[n-2] + two*v[n-1];

    // Scaling the vector w by (1 / h^2) using fortran blas routine scal.

    h2 = T((n+1)*(n+1));
    scal(n, h2, w, 1L);

  } // MultMv

  SymMatrix(int  nval) { n = nval; }

}; // SymMatrix.


template<class FLOAT, class EIGPROB>
void Solution(SymMatrix<FLOAT> &A, EIGPROB &Prob)
/*
  This function prints eigenvalues and eigenvetors on standard "cout" 
  stream and exemplifies how to retrieve information from ARPACK++ classes.
*/

{

  int   i, n, nconv, mode;
  FLOAT *Ax;
  FLOAT *ResNorm;

  /*
     ARPACK++ includes some functions that provide information
     about the problem. For example, GetN furnishes the dimension
     of the problem and ConvergedEigenvalues the number of 
     eigenvalues that attained the required accuracy. GetMode 
     indicates if the problem was solved in regular, 
     shift-and-invert or other mode.
  */

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  std::cout << std::endl << std::endl << "Testing ARPACK++ class ARSymEig \n";
  std::cout << "Real symmetric eigenvalue problem: A*x - lambda*x" << std::endl;
  switch (mode) {
  case 1:
    std::cout << "Regular mode" << std::endl << std::endl;
    break;
  case 3: 
    std::cout << "Shift and invert mode" << std::endl << std::endl;
  }

  std::cout << "Dimension of the system            : " << n             << std::endl;
  std::cout << "Number of 'requested' eigenvalues  : " << Prob.GetNev() << std::endl;
  std::cout << "Number of 'converged' eigenvalues  : " << nconv         << std::endl;
  std::cout << "Number of Arnoldi vectors generated: " << Prob.GetNcv() << std::endl;
  std::cout << std::endl;

  /*
    EigenvaluesFound is a boolean function that indicates
    if the eigenvalues were found or not. Eigenvalue can be
    used to obtain one of the "converged" eigenvalues. There
    are other functions that return eigenvectors elements,
    Schur vectors elements, residual vector elements, etc.
  */

  if (Prob.EigenvaluesFound()) {
    std::cout << "Eigenvalues:" << std::endl;
    for (i=0; i<nconv; i++) {
      std::cout << "  lambda[" << (i+1) << "]: " << Prob.Eigenvalue(i) << std::endl;
    }
    std::cout << std::endl;
  }

  /*
    EigenvectorsFound indicates if the eigenvectors are
    available. RawEigenvector is one of the functions that
    provide raw access to ARPACK++ output data. Other functions
    of this type include RawEigenvalues, RawEigenvectors, 
    RawSchurVector, RawResidualVector, etc.
  */

  if (Prob.EigenvectorsFound()) {

    // Printing the residual norm || A*x - lambda*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new FLOAT[n];
    ResNorm = new FLOAT[nconv];

    for (i=0; i<nconv; i++) {
      A.MultMv(Prob.RawEigenvector(i),Ax);
      axpy(n, -Prob.Eigenvalue(i), Prob.RawEigenvector(i), 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1) / fabs(Prob.Eigenvalue(i));
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


void SimpleTest()
{

  int  nconv;

  // Creating a double precision 100x100 matrix.

  SymMatrix<double> A(100);

  // Defining what we need: the four eigenvectors of A with smallest magnitude.

  ARSymStdEig<double, SymMatrix<double> >
    dprob(A.ncols(), 4, &A, &SymMatrix<double>::MultMv, "SM");

  /*
    It is possible to pass other parameters directly to the constructor
    of class ARSymStdEig in order to define a problem. The list
    of parameters includes, among other values, the maximum number of 
    iterations allowed and the the relative accuracy used to define the
    stopping criterion. Alternatively, it is also possible to use function 
    DefineParameters to set ARPACK++ variables after declaring dprob as an 
    object of ARSymStdEig class using the default constructor.
  */

  // Finding eigenvectors.

  nconv = dprob.FindEigenvectors();

  /*
    FindEigenvalues, FindArnoldiBasis and FindSchurVectors
    are alternatives for FindEigenvectors. However, FindSchurVectors
    can only be used with non-symmetric real and complex problems.
    Some other functions like Eigenvectors, Eigenvalues and
    EigenValVectors could also be used. These functions store the
    desired data in user supplied vectors and matrices.
  */

  // Printing solution.

  Solution(A, dprob);

} // SimpleTest


int main()
{

  // Solving a double precision problem with n = 100.

  SimpleTest();

} // main


