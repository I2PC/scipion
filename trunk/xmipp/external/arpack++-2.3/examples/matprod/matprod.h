/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE MatProd.h
   Matrix template that exemplify the data structure that
   can be used with ARPACK++. MatrixWithProduct is a base
   class for all classes in the complex, nonsym and sym
   subdirectories. 

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef MATPROD_H
#define MATPROD_H

template<class ART>
class MatrixWithProduct {

 private:

  int m, n; // Number of rows and columns.

 public:

  int nrows() { return m; }

  int ncols() { return n; }

  virtual void MultMv(ART* v, ART* w) = 0;
  // Matrix-vector product: w = M*v.

  MatrixWithProduct(int nrows, int ncols = 0)
  // Constructor.
  {
    m = nrows;
    n = (ncols?ncols:nrows);
  } // Constructor.

}; // MatrixWithProduct

#endif // MATPROD_H

