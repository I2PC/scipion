/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or   
 * (at your option) any later version.                                 
 *                                                                     
 * This program is distributed in the hope that it will be useful,     
 * but WITHOUT ANY WARRANTY; without even the implied warranty of      
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
 * GNU General Public License for more details.                        
 *                                                                     
 * You should have received a copy of the GNU General Public License   
 * along with this program; if not, write to the Free Software         
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA            
 * 02111-1307  USA                                                     
 *                                                                     
 *  All comments concerning this program package may be sent to the    
 *  e-mail address 'xmipp@cnb.uam.es'                                  
 ***************************************************************************/
/* ------------------------------------------------------------------------- */
/* MATRICES                                                                  */
/* ------------------------------------------------------------------------- */
#ifndef _XMIPPMATRICES_HH
#   define _XMIPPMATRICES_HH

/* ************************************************************************* */
/* INCLUDES                                                                  */
/* ************************************************************************* */

#include <stdlib.h>
#include <iostream>
#include <string>
#include <math.h>
#include <complex>

#include "xmippMatrices1D.hh"
#include "xmippFuncs.hh"
#include "Bilib/types/tsplinebasis.h"
#include "Bilib/types/tboundaryconvention.h"
#include "Bilib/headers/linearalgebra.h"
#include "Bilib/headers/changebasis.h"
#include "Bilib/headers/kernel.h"
#include "Bilib/headers/pyramidfilters.h"
#include "Bilib/headers/pyramidtools.h"

#define maT  matrix2D<T>
#define maT1 matrix2D<T1>
#undef  maTC
#define maTC matrix2D< complex<double> >

/* ************************************************************************* */
/* FORWARD DEFINITIONS                                                       */
/* ************************************************************************* */
#include "Src/MultidimFriends.inc"

template <class T>
   mT mul_elements(const mT &op1, const mT &op2);
template <class T>
   void mul_matrix(const mT &op1, const mT &op2, mT &result);
template <class T>
   void solve (const mT &A, const vT &b, vT &result);
template <class T>
   void solve (const mT &A, const mT &b, mT &result);
// Use this solve_by_svd function to solve linear systems by singular value decomposition
// This is indicated for ill-conditioned systems and hard-to-solve problems
// (But it's slower)
template <class T>
   void solve_by_svd(const mT &A, const vT &b, matrix1D<double> &result,double tolerance);
template <class T>
   void ludcmp(const mT &A, mT &LU, matrix1D<int> &indx, T &d);
template <class T>
   void lubksb(const mT &LU,matrix1D<int> &indx, vT &b);

template <class T>
void svdcmp(const matrix2D<T> &a,matrix2D<double> &u,
               matrix1D<double> &w, matrix2D<double> &v);
void svbksb(matrix2D<double> &u,matrix1D<double> &w,matrix2D<double> &v,
             matrix1D<double> &b,matrix1D<double> &x);
			   
template <class T>
   void apply_geom(mT &m2, matrix2D<double> A, const mT &m1, bool inv,
      bool wrap, T outside=(T)0);

template <class T>
   void apply_geom_Bspline(mT &m2, matrix2D<double> A, const mT &m1,
      int Splinedegree, bool inv, bool wrap, T outside=(T)0);

matrix2D<double> rot2D_matrix(double ang);
matrix2D<double> translation2D_matrix(const matrix1D<double> v);
int best_prec(float F, int _width);

/* ************************************************************************* */
/* CLASS DEFINITION AND PROTOTYPES                                           */
/* ************************************************************************* */
/**@name Xmipp Matrices*/
//@{
/* Speed up ---------------------------------------------------------------- */
/**@name Speed up macros
   This macros are defined to allow high speed in critical parts of
   your program. They shouldn't be used systematically as usually there
   is no checking on the correctness of the operation you are performing.
   Speed comes from three facts: first, they are macros and no function
   call is performed (although most of the critical functions are
   inline functions), there is no checking on the correctness of the
   operation (it could be wrong and you are not warned of it), and
   destination vectors are not returned saving time in the copy
   constructor and in the creation/destruction of temporary vectors.*/
//@{
/**@name Size and shape
    Although they are not defined here you can also use STARTINGX and
    FINISHINGX (defined for matrix1D)*/
//@{
/** TRUE if both arrays have the same shape.
    Two arrays have the same shape if they have the same size and the
    same starting point. Be aware that this is a macro which simplifies to
    a boolean. */
#define SAME_SHAPE2D(v1,v2) \
    (XSIZE(v1)==XSIZE(v2) && \
     YSIZE(v1)==YSIZE(v2) && \
     STARTINGX(v1)==STARTINGX(v2) && \
     STARTINGY(v1)==STARTINGY(v2))

/** Returns the first valid logical Y index.
    \\Ex: int orgY=STARTINGY(m);*/
#define STARTINGY(m)  ((m).yinit)

/** Returns the last valid logical Y index.
    \\Ex: int finY=FINISHINGY(m);*/
#define FINISHINGY(m) ((m).yinit+(m).ydim-1)

/** Access to Y dimension (size).
    This is a macro equivalent to \Ref{RowNo()}
    \\ Ex:
    \begin{verbatim}
    // Set to 0 1 element out of 4
    for (int i=0; i<YSIZE(m); i +=2)
        for (int j=0; j<XSIZE(m); j +=2)
            DIRECT_MAT_ELEM(v,i,j)=0;
   \end{verbatim}*/
#define YSIZE(m) ((m).ydim)

/** For all elements in the array.
    This macro is used to generate loops for the matrix in an easy way.
    It defines internal indexes 'i' and 'j' which ranges the matrix using
    its mathematical definition (ie, logical access).
    \\Ex:
    \begin{verbatim}
    FOR_ALL_ELEMENTS_IN_MATRIX2D(m) {
       cout << m(i,j) << " ";
    }
    \end{verbatim} */
#define FOR_ALL_ELEMENTS_IN_MATRIX2D(m) \
    for (int i=STARTINGY(m); i<=FINISHINGY(m); i++) \
       for (int j=STARTINGX(m); j<=FINISHINGX(m); j++)

/** For all elements in the array between corners.
    This macro is used to generate loops for a matrix in an easy manner.
    It needs an externally defined matrix1D<double> r(2).    
    Then YY(r) and XX(r) range from
    (int) YY(corner1) to (int)YY(corner2), (int) XX(corner1) to
    (int) XX(corner2) (included limits) respectively. Notice that corner1
    and corner2 need only be matrix1D. 
    \\Ex:
    \begin{verbatim}
    matrix1D<double> corner1(2), corner2(2);
    matrix1D<int> r(2);
    XX(corner1)=-1; XX(corner2)=1;
    YY(corner1)=-2; YY(corner2)=2;
    FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1,corner2) {
       cout << v(r) << " ";
    }
    \end{verbatim}*/
#define FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1,corner2) \
    for (YY(r)=YY((corner1)); YY(r)<=YY((corner2)); YY(r)++) \
       for (XX(r)=XX((corner1)); XX(r)<=XX((corner2)); XX(r)++)

/** For all elements in common.
    This macro is used to generate loops for all the elements logically
    in common between two images in an easy manner.
    Then i and j (locally defined) range from
    MAX(STARTINGY(V1),STARTINGY(V2)) to MIN(FINISHINGY(V1),FINISHINGY(V2)),
    MAX(STARTINGX(V1),STARTINGX(V2)) to MIN(FINISHINGX(V1),FINISHINGX(V2))
    (included limits) respectively. You need to define SPEED_UP_temps.
    \\Ex:
    \begin{verbatim}
    matrix2D<double> m1(10,10), m2(20,20);
    m1.set_Xmipp_origin();
    m2.set_Xmipp_origin();
    FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(m1,m2) {
       ...
    }
    \end{verbatim}*/
#define FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(m1,m2) \
    ispduptmp2=MAX(STARTINGY(m1), STARTINGY(m2)); \
    ispduptmp3=MIN(FINISHINGY(m1),FINISHINGY(m2)); \
    ispduptmp4=MAX(STARTINGX(m1), STARTINGX(m2)); \
    ispduptmp5=MIN(FINISHINGX(m1),FINISHINGX(m2)); \
    for (int i=ispduptmp2; i<=ispduptmp3; i++) \
       for (int j=ispduptmp4; j<=ispduptmp5; j++)

/** For all elements in the array, accessed physically.
    This macro is used to generate loops for the matrix in an easy way
    using physical indexes.
    It defines internal indexes 'i' and 'j' which ranges the matrix using
    its physical definition.
    \\Ex:
    \begin{verbatim}
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(m) {
       cout << DIRECT_MAT_ELEM(m,i,j) << " ";
    }
    \end{verbatim} */
#define FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(m) \
    for (int i=0; i<YSIZE(m); i++) \
       for (int j=0; j<XSIZE(m); j++)

//@}

/**@name Memory access*/
//@{
/** Matrix element: Logical access.
    \\ Ex: MAT_ELEM(m,-2,1)=1;
    \\ Ex: val=MAT_ELEM(m,-2,1);*/
#define MAT_ELEM(mat,i,j) \
   DIRECT_MAT_ELEM(mat,(i)-STARTINGY(mat),(j)-STARTINGX(mat))

/** Matrix element: Physical access.
    Be careful because this is physical access, usually matrices follow
    the C convention of starting index==0 (X and Y). This function should
    not be used as it goes against the vector library philosophy unless you
    explicitly want to access directly to any value in the matrix
    without taking into account its logical position
    \\ Ex: DIRECT_MAT_ELEM(m,0,0)=1;
    \\ Ex: val=DIRECT_MAT_ELEM(m,0,0);*/
#define DIRECT_MAT_ELEM(mat,i,j) (mat).__m[(i)*XSIZE(mat)+(j)]

/** Short alias for the previous function. */
#define dMij(M,i,j) DIRECT_MAT_ELEM(M,i,j)

/** Array access.
    This macro gives you access to the array (T *).
    \\ Ex: cout << "This is an int *" << MAT_ARRAY(m) << endl; */
#define MAT_ARRAY(m) MULTIDIM_ARRAY(m)
//@}

/**@name Arithmethic operations
   The vectors and matrices involved in these macros should be created
   with the correct size before entering in them. These macros allow a
   fast operation on R2 and R3 vectors, and small size matrices.
   These macros need some temporary variables. You must "load" them
   by calling the macro SPEED_UP_temps at the beginning of your
   function*/
//@{
/** Matrix (3x3) by vector (3x1) (a=M*b).
    You must "load" the temporary variables, and create the result
    vector with the appropiate size. You can reuse the vector b to
    store the results (that is, M3x3_BY_V3x1(b,M,b);, is allowed).
    \\ Ex: 
    \begin{verbatim}
           double example {
              SPEED_UP_temps;
              matrix1D<double> a(3), b(3);
              matrix2D<double> M(3,3);
              
              M.init_random(0,1);
              b.init_random(0,1);
              M3x3_BY_V3x1(a,M,b);
              
              return a.sum();
           }
    \end{verbatim}*/
#define M3x3_BY_V3x1(a,M,b) {\
    spduptmp0=dMij(M,0,0)*XX(b)+dMij(M,0,1)*YY(b)+dMij(M,0,2)*ZZ(b);\
    spduptmp1=dMij(M,1,0)*XX(b)+dMij(M,1,1)*YY(b)+dMij(M,1,2)*ZZ(b);\
    spduptmp2=dMij(M,2,0)*XX(b)+dMij(M,2,1)*YY(b)+dMij(M,2,2)*ZZ(b);\
    XX(a)=spduptmp0; YY(a)=spduptmp1; ZZ(a)=spduptmp2;}

/** Matrix (3x3) by Matrix (3x3) (A=B*C).
    You must "load" the temporary variables, and create the result
    vector with the appropiate size. You can reuse any of the
    multiplicands to
    store the results (that is, M3x3_BY_M3x3(A,A,B);, is allowed). */
#define M3x3_BY_M3x3(A,B,C) {\
    spduptmp0=dMij(B,0,0)*dMij(C,0,0)+dMij(B,0,1)*dMij(C,1,0)+dMij(B,0,2)*dMij(C,2,0);\
    spduptmp1=dMij(B,0,0)*dMij(C,0,1)+dMij(B,0,1)*dMij(C,1,1)+dMij(B,0,2)*dMij(C,2,1);\
    spduptmp2=dMij(B,0,0)*dMij(C,0,2)+dMij(B,0,1)*dMij(C,1,2)+dMij(B,0,2)*dMij(C,2,2);\
    spduptmp3=dMij(B,1,0)*dMij(C,0,0)+dMij(B,1,1)*dMij(C,1,0)+dMij(B,1,2)*dMij(C,2,0);\
    spduptmp4=dMij(B,1,0)*dMij(C,0,1)+dMij(B,1,1)*dMij(C,1,1)+dMij(B,1,2)*dMij(C,2,1);\
    spduptmp5=dMij(B,1,0)*dMij(C,0,2)+dMij(B,1,1)*dMij(C,1,2)+dMij(B,1,2)*dMij(C,2,2);\
    spduptmp6=dMij(B,2,0)*dMij(C,0,0)+dMij(B,2,1)*dMij(C,1,0)+dMij(B,2,2)*dMij(C,2,0);\
    spduptmp7=dMij(B,2,0)*dMij(C,0,1)+dMij(B,2,1)*dMij(C,1,1)+dMij(B,2,2)*dMij(C,2,1);\
    spduptmp8=dMij(B,2,0)*dMij(C,0,2)+dMij(B,2,1)*dMij(C,1,2)+dMij(B,2,2)*dMij(C,2,2);\
    dMij(A,0,0)=spduptmp0; dMij(A,0,1)=spduptmp1; dMij(A,0,2)=spduptmp2; \
    dMij(A,1,0)=spduptmp3; dMij(A,1,1)=spduptmp4; dMij(A,1,2)=spduptmp5; \
    dMij(A,2,0)=spduptmp6; dMij(A,2,1)=spduptmp7; dMij(A,2,2)=spduptmp8;}

/** Matrix (2x2) by vector (2x1) (a=M*b).
    You must "load" the temporary variables, and create the result
    vector with the appropiate size. You can reuse the vector b to
    store the results (that is, M2x2_BY_V2x1(b,M,b);, is allowed).
    \\ Ex: 
    \begin{verbatim}
           double example {
              SPEED_UP_temps;
              matrix1D<double> a(2), b(2);
              matrix2D<double> M(2,2);
              
              M.init_random(0,1);
              b.init_random(0,1);
              M2x2_BY_V2x1(a,M,b);
              
              return a.sum();
           }
    \end{verbatim}*/
#define M2x2_BY_V2x1(a,M,b) {\
    spduptmp0=dMij(M,0,0)*XX(b)+dMij(M,0,1)*YY(b);\
    spduptmp1=dMij(M,1,0)*XX(b)+dMij(M,1,1)*YY(b);\
    XX(a)=spduptmp0; YY(a)=spduptmp1;}

/** Matrix (2x2) by constant (M2=M1*k).
    You must create the result
    matrix with the appropiate size. You can reuse the matrix M1 to
    store the results (that is, M2x2_BY_CT(M,M,k);, is allowed).
*/
#define M2x2_BY_CT(M2,M1,k) {\
    dMij(M2,0,0)=dMij(M1,0,0)*k; \
    dMij(M2,0,1)=dMij(M1,0,1)*k; \
    dMij(M2,1,0)=dMij(M1,1,0)*k; \
    dMij(M2,1,1)=dMij(M1,1,1)*k;}

/** Matrix (3x3) by constant (M2=M1*k).
    You must create the result
    matrix with the appropiate size. You can reuse the matrix M1 to
    store the results (that is, M2x2_BY_CT(M,M,k);, is allowed).
*/
#define M3x3_BY_CT(M2,M1,k) {\
    dMij(M2,0,0)=dMij(M1,0,0)*k; \
    dMij(M2,0,1)=dMij(M1,0,1)*k; \
    dMij(M2,0,2)=dMij(M1,0,2)*k; \
    dMij(M2,1,0)=dMij(M1,1,0)*k; \
    dMij(M2,1,1)=dMij(M1,1,1)*k; \
    dMij(M2,1,2)=dMij(M1,1,2)*k; \
    dMij(M2,2,0)=dMij(M1,2,0)*k; \
    dMij(M2,2,1)=dMij(M1,2,1)*k; \
    dMij(M2,2,2)=dMij(M1,2,2)*k;}

/** Inverse of a matrix (2x2).
    Input and output matrix cannot be the same one. The output is
    supposed to be already resized. */
#define M2x2_INV(Ainv,A) {\
   spduptmp0=1.0/(dMij(A,0,0)*dMij(A,1,1)-dMij(A,0,1)*dMij(A,1,0)); \
   dMij(Ainv,0,0)= dMij(A,1,1);       dMij(Ainv,0,1)=-dMij(A,0,1); \
   dMij(Ainv,1,0)=-dMij(A,1,0);       dMij(Ainv,1,1)= dMij(A,0,0); \
   M2x2_BY_CT(Ainv,Ainv,spduptmp0);}
//@}
//@}

/// Template class for Xmipp matrices
#include "MultidimCommon.hh"
template <class T> class matrix2D {
#include "Src/MultidimBasic.hh"
/**@name Common functions to all multidimensional arrays
   A set of methods are always the same for any multidimensional array.
   Have a look on the more detailed structure. */
public:
//@{
//@Include: Src/MultidimBasic2.hh
//@}
   // This file contains several definitions for statistics and arithmetic
   // operations. To use it we have redirected the internal type maT
   // (multidimensional array<T>) to matrix<T>. These definitions are
   // outside because in this way we can reuse the module for other
   // libraries

/* Structure =============================================================== */
// Although the structure is defined as public it should not be used by
// the library user, there are functions enough to handle everything. This is
// done so to speed up the library
public:
   int        ydim,xdim;         // dimensions of array [0...ydim-1]
                                 //                     [0...xdim-1]
   int        yinit,xinit;       // indexes of array  [yinit...yinit+ydim-1]
                                 //                   [xinit...xinit+xdim-1]

/* Procedures ============================================================== */
protected:
   /* Operation related ---------------------------------------------------- */
   // Matrix multiplication in a algebraic way
   friend void mul_matrix<>(const mT &op1, const mT &op2, mT &result);

   // The following functions are directly an interface to the Numerical
   // Recipes ones with the same name, but with a better interface
   // LU decomposition
   friend void ludcmp<>(const mT &A, mT &LU, matrix1D<int> &indx, T &d);

   // Solving an equation system based on LU. Remember to free LU outside
   friend void lubksb<>(const mT &LU, matrix1D<int> &indx, vT &b);
          
public:
   /* Constructors/Destructor ---------------------------------------------- */
   /**@name Constructors*/
   //@{
   /** Empty constructor.
       The empty constructor creates a matrix with no memory associated,
       origin=0, size=0, no statistics, ...
       \\Ex: matrix2D<double> m1;*/
   matrix2D() {core_init(); init_shape(); __spcdim=2;}

   /** Dimension constructor.
       The dimension constructor creates a matrix with memory associated
       (but not assigned to anything, could be full of garbage)
       origin=0, size=the given one
       Be careful that first number is the Y dimension (number of rows),
       and the second the X dimension (number of columns).
       \\Ex: matrix2D<double> v1(6,3);*/
   matrix2D(int Ydim, int Xdim) {
      core_init();
      init_shape();
      __m = new T[Ydim * Xdim];

      if (__m == NULL)
         REPORT_ERROR(1001, "Resize: no memory left");

      xdim = Xdim;
      ydim = Ydim;
      __dim=xdim*ydim;
      __spcdim = 2;
      for (long int i=0; i<__dim; i++) __m[i]=0;
   }

   /** Copy constructor.
       The created matrix is a perfect copy of the input matrix but
       with a different memory assignment.
       \\Ex: matrix2D<double> m2(m1); */
   matrix2D(const mT &v) {core_init(); init_shape(); *this=v;}

   // Destructor
   ~matrix2D() {core_deallocate();}
   //@}

   /* Initialisation ------------------------------------------------------- */
   /**@name Initialisation*/
   //@{
   /** Identity matrix of current size.
       If actually the matrix is not squared then an identity matrix
       is generated of size (Xdim x Xdim).
       \\Ex: m.init_identity();*/
   void init_identity() {init_identity(XSIZE(*this));}

   /** Identity matrix of a given size.
       A (dim x dim) identity matrix is generated.
       \\Ex: m.init_identity(3);*/
   void init_identity(int dim) {init_identity(dim, dim);}
   
   /** Identity matrix of a given size.
       A (dimX x dimY) identity matrix is generated. That is, any element i,j 
       of the matrix such that i = j is equal to 1.
       \\Ex: m.init_identity(2,3);*/
   void init_identity(int Ydim, int Xdim) {
      if (Xdim==0 || Ydim==0) {clear(); return;}
      resize(Ydim,Xdim);
      FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
         DIRECT_MAT_ELEM(*this,i,j)=(T)(i==j);
   }
   
   /** Zero initialisation with a new dimension.
       Be careful to the size order (Ydim, Xdim).
       \\Ex: v1.init_zeros(6,3);*/         
   void init_zeros(int Ydim, int Xdim)
         {resize(Ydim,Xdim); init_constant((T)0);}

   /** Makes a matrix from a vector.
       The origin of the matrix is set such that it has one of the
       index origins (X or Y) to the same value as the vector, and
       the other set to 0 according to the shape.
       \\Ex: matrix2D<double> m=from_vector(v);*/
   void from_vector(const vT &op1) {
      // Null vector => Null matrix
      if (XSIZE(op1)==0) {clear(); return;}

      // Look at shape and copy values
      if (op1.isRow()) {
         resize(1,XSIZE(op1));
         for (int j=0; j<XSIZE(op1); j++)
            DIRECT_MAT_ELEM(*this,0,j)=DIRECT_VEC_ELEM(op1,j);
         STARTINGX(*this)=STARTINGX(op1);
         STARTINGY(*this)=0;
      } else {
         resize(XSIZE(op1),1);
         for (int i=0; i<XSIZE(op1); i++)
            DIRECT_MAT_ELEM(*this,i,0)=DIRECT_VEC_ELEM(op1,i);
         STARTINGX(*this)=0;
         STARTINGY(*this)=STARTINGX(op1);
      }
   }
   
   /** Makes a vector from a matrix.
       An exception is thrown if the matrix is not a single row or
       a single column. The origin of the vector is set according to
       the one of the matrix.
       \\Ex: matrix1D<double> v; m.to_vector(v);*/
   void to_vector(vT &op1) const {
      // Null matrix => Null vector
      if (XSIZE(*this)==0 || YSIZE(*this)==0) {op1.clear(); return;}

      // If matrix is not a vector, produce an error
      if (XSIZE(*this)!=1 && (YSIZE(*this)!=1))
         REPORT_ERROR(1102, "To_vector: Matrix cannot be converted to vector");

      // Look at shape and copy values
      if (YSIZE(*this)==1) {
         // Row vector
         op1.resize(XSIZE(*this));
         for (int j=0; j<XSIZE(*this); j++)
             DIRECT_VEC_ELEM(op1,j)=DIRECT_MAT_ELEM(*this,0,j);
         op1.setRow();
         STARTINGX(op1)=STARTINGX(*this);
      } else {
         // Column vector
         op1.resize(YSIZE(*this));
         for (int i=0; i<YSIZE(*this); i++)
            DIRECT_VEC_ELEM(op1,i)=DIRECT_MAT_ELEM(*this,i,0);
         op1.setCol();
         STARTINGX(op1)=STARTINGY(*this);
      }
   }
   //@}

   /* Memory related ------------------------------------------------------- */
   /**@name Size and shape
      The shape of a matrix is defined by its origin and its size.
      The size is clear, and the origin
      is the logical position of the first real position of the array. For
      instance, if we have a matrix of dimension (5,3)=(Ydim, Xdim) and
      origin (-2,-1), this means
      that the array is representing the logical positions
      \begin{verbatim}
      [(-2,-1) (-2,0) (-2,1)
       (-1,-1) (-1,0) (-1,1)
       ( 0,-1) ( 0,0) ( 0,1)
       ( 1,-1) ( 1,0) ( 1,1)
       ( 2,-1) ( 2,0) ( 2,1)]
      \end{verbatim}
      we could access to any of these positions (Ex: v(-2,1)=3;) and actually
      any try to access to a position related to 5 (Ex: v(4,1)=3;), although
      it physically exists, is not logically correct and hence it will
      throw an exception. The startingX and finishingX positions for this
      sample vector are -1 and 1 respectively, while for Y are -2 and 2.
      The "for" iterations through the matrix should include these two
      values if you want to cover the whole matrix.

      \begin{verbatim}
         for (int i=STARTINGY(m); i<=FINISHINGY(m); i++)
            for (int j=STARTINGX(m); j<=FINISHINGX(m); j++)
               MAT_ELEM(m,i,j) += 1;
      \end{verbatim}*/
      
   //@{
   /** Init shape.
       ydim,xdim=0, startingy,startingx=0.*/
   void init_shape()
      {xinit=yinit=0; xdim=ydim=0;}

   /** Copy shape.
       Copy shape variables from a pattern AND THE ARRAY IS RESIZED */
   template <class T1>
      void copy_shape(const maT1 &v)
         {if (XSIZE(*this)!=XSIZE(v) || YSIZE(*this)!=YSIZE(v))
             resize(YSIZE(v), XSIZE(v));
          STARTINGX(*this)=STARTINGX(v);
          STARTINGY(*this)=STARTINGY(v);}

   /** Resize to a given size.
       This function resize the actual array to the given size. The origin
       is not modified. If the actual array is larger than the pattern
       then the values outside the new size are lost, if it is smaller
       then 0's are added. An exception is thrown if there is no memory.
       \\Ex: m1.resize(3,2);
   */
   void resize(int Ydim, int Xdim) {
      if (Xdim==XSIZE(*this) && Ydim==YSIZE(*this)) return;
      if (Xdim<=0 || Ydim<=0)                       {clear(); return;}

      // Ask for memory   
      T* new_m=new T [Ydim*Xdim];
      if (new_m==NULL) REPORT_ERROR(1001,"Resize: no memory left");

      // Copy needed elements, fill with 0 if necessary
      for (int i=0; i<Ydim; i++)
         for (int j=0; j<Xdim; j++) {
            T val;
            if      (i>=YSIZE(*this)) val=0;
            else if (j>=XSIZE(*this)) val=0;
            else                      val=DIRECT_MAT_ELEM(*this,i,j);
            new_m[Xdim*i+j]=val;
         }

      // deallocate old vector
      core_deallocate();

      // assign *this vector to the newly created
      MULTIDIM_ARRAY(*this)=new_m;
      XSIZE(*this)=Xdim;
      YSIZE(*this)=Ydim;
      __dim=Ydim*Xdim;
   }

   /** Produce an array suitable for working with Numerical Recipes.
       This function must be used only as
       a preparation for routines which need that the first physical
       index is 1 and not 0 as it usually is in C. New memory is needed
       to hold the new double pointer array.
       Click here to see an \URL[example]{../Extra_Docs/examples.html#LU} */
   T ** adapt_for_numerical_recipes() const
      {T **m=NULL; ask_Tmatrix(m, 1, YSIZE(*this), 1, XSIZE(*this));
      FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(*this)
         m[i+1][j+1]=DIRECT_MAT_ELEM(*this,i,j);
      return m;}

   /** Produce a pointer suitable for working with Numerical Recipes (2).
       This function meets the same goal as the one before, however this one
       work with 2D arrays as a single pointer. The first element of
       the array is pointed by result[1*Xdim+1], and in general
       result[i*Xdim+j] */
   T * adapt_for_numerical_recipes2() const {return __m-1-XSIZE(*this);}

   /** Load from numerical recipes result. */
   void load_from_numerical_recipes(T **m, int Ydim, int Xdim) {
      resize(Ydim,Xdim);
      for (int i=1; i<=Ydim; i++)
         for (int j=1; j<=Xdim; j++)
            (*this)(i-1,j-1)=m[i][j];
   }

   /** Kill an array produced for numerical recipes.
       The allocated memory is freed. */
   void kill_adaptation_for_numerical_recipes(T **m) const
      {free_Tmatrix(m, 1, YSIZE(*this), 1, XSIZE(*this));}

   /** Kill an array produced for numerical recipes, 2.
       Nothing needs to be done in fact. */
   void kill_adaptation_for_numerical_recipes2(T **m) const {}

   /** Intersects.
       TRUE if this array intersects with the rectangle defined by the
       arguments (x0 is the starting X).*/
   bool intersects(double x0, double y0, double xdim, double ydim) const;

   /** Outside.
       TRUE if the logical index given is outside the definition region
       of this array. */
   bool outside(int i, int j) const;

   /** isBorder.
       TRUE if the logical index given belong to the border of the matrix */
   bool isBorder(int i,int j);

   /** Set logical origin in Xmipp fashion.
       This function adjust the starting points in the matrix such that
       the center of the matrix is defined in the Xmipp fashion
       (see \Ref{Conventions}, \Ref{FIRST_XMIPP_INDEX}).
       \\Ex: m1.rigin();*/
   void set_Xmipp_origin()
      {yinit=FIRST_XMIPP_INDEX(ydim);
       xinit=FIRST_XMIPP_INDEX(xdim);}

   /** Move origin to.
       This function adjust logical indexes such that the Xmipp origin
       of the array moves to the specified position. For instance, an array
       whose x indexes go from -1 to 1, if we move the origin to 4, then
       the x indexes go from 3 to 5. This is very useful for convolution
       operations where you only need to move the logical starting of the
       array. 
       See \Ref{FIRST_XMIPP_INDEX} */
   void move_origin_to(int i, int j) {
      yinit=i+FIRST_XMIPP_INDEX(ydim);
      xinit=j+FIRST_XMIPP_INDEX(xdim);}

   /** Sets the Y origin.
       The logical position of the first physical Y position is set with
       this function. By default the origin is 0 that is the standard
       convention in C.
       \\Ex: m.startingY(-2); */
   int& startingY() {return yinit;}

   /** Another function for setting the Y origin. */
   void set_startingY(int _yinit) {yinit=_yinit;}

   /** Returns the first valid logical Y index.
       \\Ex: int orgY=m.startingY();*/
   int startingY() const {return yinit;}

   /** Returns the last valid logical Y index.
       \\Ex: int finY=m.finishingY();*/
   int finishingY() const {return yinit+ydim-1;}

   /** Sets the X origin.
       The logical position of the first physical X position is set with
       this function. By default the origin is 0 that is the standard
       convention in C.
       \\Ex: m.startingX(-1); */
   int& startingX() {return xinit;}

   /** Another function for setting the X origin. */
   void set_startingX(int _xinit) {xinit=_xinit;}

   /** Returns the first valid logical X index.
       \\Ex: int orgX=m.startingX();*/
   int startingX() const {return xinit;}
   
   /** Returns the last valid logical X index.
       \\Ex: int finX=m.finishingX();*/
   int finishingX() const {return xinit+xdim-1;}

   /** Returns the matrix dimension.
       Pay attention to the dimension order (Y,X).
       \\ Ex: m.get_dim(Ydim,Xdim);*/
   void get_dim(int &Ydim, int &Xdim) const
         {Xdim=xdim; Ydim=ydim;}
   
   /** Returns Y dimension.
       \\Ex: int Ydim=m.RowNo();*/
   int  RowNo() const {return ydim;}

   /** Returns X dimension.
       \\Ex: int Xdim=m.ColNo();*/
   int  ColNo() const {return xdim;}

   /** Same shape.
       Returns true if this object has got the same shape (origin and
       size) than the argument*/
   bool same_shape(const maT &op) const
      {return SAME_SHAPE2D(*this,op);}
   //@}
      
   /* Information Extraction ----------------------------------------------- */
   /**@name Memory access
      This functions allows you to access to the matrix elements.*/
   //@{
   /** Matrix element access via index.
       Returns the value of a matrix logical position. In our example we could
       access from v(-2,-1) to v(2,1). The elements can be used either by value
       or by reference. An exception is thrown if the index is outside the
       logical range. The first argument is the Y position and the second
       the X position.
       \\ Ex: m(-2,1)=1;
       \\ Ex: val=m(-2,1);*/
   T&   operator () (int i, int j) const {
        if (i<yinit || i>=yinit+ydim)
           REPORT_ERROR(1103,"Matrix subscript (i) out of range");
        if (j<xinit || j>=xinit+xdim)
           REPORT_ERROR(1103,"Matrix subscript (j) out of range");
        return MAT_ELEM(*this,i,j);}

   /** Get the pixel at (i,j). Logical access */
   T get_pixel(int i,int j) const {return (*this)(i,j);}

   /** Set the pixel at (i,j). Logical access */
   void set_pixel(int i,int j, T val) {(*this)(i,j)=val;}

   /** Matrix element access via double vector.
       Returns the value of a matrix logical position, but this time the
       element position is determined by a R2 vector.
       The elements can be used either by value or by reference.
       An exception is thrown if the index is outside the logical range.
       Pay attention in the following example that we are accessing the
       same element as in the previous function but, now we have to give
       first the X position instead of the Y one because we are building
       first a vector of the form (x,y).
       \\ Ex: m(vector_R2(1,-2))=1;
       \\ Ex: val=m(vector_R2(1,-2));*/
   T&   operator () (const matrix1D<double> &v) const
        {return MAT_ELEM(*this,ROUND(YY(v)),ROUND(XX(v)));}

   /** Matrix element access via intger vector.*/
   T&   operator () (const matrix1D<int> &v) const
        {return MAT_ELEM(*this,YY(v),XX(v));}

   /** Interpolates the value of the 2D matrix M at the point (x,y).
       Bilinear interpolation.
       (x,y) are in logical coordinates.*/
   T   interpolated_elem(double x, double y, T outside_value=(T)0) const {
       int x0 = FLOOR(x); double fx = x - x0; int x1 = x0 + 1;
       int y0 = FLOOR(y); double fy = y - y0; int y1 = y0 + 1;

       T d00 = outside(y0,x0) ? outside_value : MAT_ELEM(*this,y0,x0);
       T d10 = outside(y1,x0) ? outside_value : MAT_ELEM(*this,y1,x0);
       T d11 = outside(y1,x1) ? outside_value : MAT_ELEM(*this,y1,x1);
       T d01 = outside(y0,x1) ? outside_value : MAT_ELEM(*this,y0,x1);

       double d0 = (T) LIN_INTERP(fx, (double)d00, (double)d01);
       double d1 = (T) LIN_INTERP(fx, (double)d10, (double)d11);    
       return (T) LIN_INTERP(fy, d0, d1);
   }

   /** Interpolates the value of the 2D matrix M at the point (x,y) knowing
       that this image is a set of B-spline coefficients.
       (x,y) are in logical coordinates.*/
   T   interpolated_elem_as_Bspline(double x, double y, int SplineDegree=3) const {
       int SplineDegree_1=SplineDegree-1;

       // Logical to physical
       y-=STARTINGY(*this);
       x-=STARTINGX(*this);

       int lmax = XSIZE(*this);
       int mmax = YSIZE(*this);
       int l1 = CLIP(CEIL(x - SplineDegree_1),0,lmax-1);
       int l2 = CLIP(l1 + SplineDegree,0,lmax-1);
       int m1 = CLIP(CEIL(y - SplineDegree_1),0,mmax-1);
       int m2 = CLIP(m1 + SplineDegree,0,mmax-1);
       double columns = 0.0;
       for (int m=m1; m<=m2; m++) {
           int row_m = XSIZE(*this)*m;
      	   double rows = 0.0;
           for (int l=l1; l<=l2; l++) {
               double xminusl = x-(double)l;
               double Coeff = (double) __m[row_m+l];
	       switch (SplineDegree) {
	          case 2: rows += Coeff * Bspline02(xminusl); break; 
                  case 3: rows += Coeff * Bspline03(xminusl); break;
                  case 4: rows += Coeff * Bspline04(xminusl); break;
                  case 5: rows += Coeff * Bspline05(xminusl); break;
                  case 6: rows += Coeff * Bspline06(xminusl); break;
                  case 7: rows += Coeff * Bspline07(xminusl); break;
                  case 8: rows += Coeff * Bspline08(xminusl); break;
                  case 9: rows += Coeff * Bspline09(xminusl); break;
	       }
           }
           double yminusm=y-(double)m;
	   switch (SplineDegree) {
              case 2: columns += rows * Bspline02(yminusm); break;
              case 3: columns += rows * Bspline03(yminusm); break;
              case 4: columns += rows * Bspline04(yminusm); break;
              case 5: columns += rows * Bspline05(yminusm); break;
              case 6: columns += rows * Bspline06(yminusm); break;
              case 7: columns += rows * Bspline07(yminusm); break;
              case 8: columns += rows * Bspline08(yminusm); break;
              case 9: columns += rows * Bspline09(yminusm); break;
	   }
       }
       return (T)columns;
   }

   /** Extracts the profile between two points.
       Given two logical indexes, this function returns samples of
       the line that joins them. This is done by bilinear interpolation.
       The number of samples in the line is N. */
   void profile(int x0, int y0, int xF, int yF, int N,
      matrix1D<double> &profile) const {
      profile.init_zeros(N);
      double tx_step=(double)(xF-x0)/(N-1);
      double ty_step=(double)(yF-y0)/(N-1);
      double tx=x0, ty=y0;
      for (int i=0; i<N; i++) {
         profile(i)=interpolated_elem(tx,ty);
         tx+=tx_step;
         ty+=ty_step;
      }
   }

   /** Logical to physical index translation.
       This function returns the physical position of a logical one. See
       \URL[Conventions]{../../../Extra_Docs/Conventions.html}
       for more information about these two different accesses.
       \\ Ex: m.logical2physical(i_log,j_log,i_phys,j_phys); */
   void logical2physical(int i_log, int j_log, int &i_phys, int &j_phys) const
       {i_phys=i_log-yinit; j_phys=j_log-xinit;}
   
   /** Physical to logical index translation.
       This function returns the logical position of a physical one. See
       \URL[Conventions]{../../../Extra_Docs/Conventions.html}
       for more information about these two different accesses.
       \\ Ex: m.physical2logical(i_phys,j_phys,i_log,j_log); */
   void physical2logical(int i_phys, int j_phys, int &i_log, int &j_log) const
       {i_log=i_phys+yinit; j_log=j_phys+xinit;}
   
   /** Get row.
       This function returns a row vector corresponding to the choosen
       row inside matrix, the numbering of the rows is also logical not
       physical.
       \\Ex: vector<double> v; m.getRow(-2,v);*/
   void getRow(int i, vT &v) const {
      if (XSIZE(*this)==0 || YSIZE(*this)==0) {v.clear(); return;}
      if (i<STARTINGY(*this) || i>FINISHINGY(*this))
         REPORT_ERROR(1103,"getRow: Matrix subscript (i) greater than matrix dimension");

      v.resize(XSIZE(*this));
      STARTINGX(v)=STARTINGX(*this);
      for (int j=STARTINGX(*this); j<=FINISHINGX(*this); j++)
          VEC_ELEM(v,j)=MAT_ELEM(*this,i,j);
      v.setRow();
   }
   
   /** Return row. The same as previous. */
   vT Row(int i) const {vT aux; getRow(i,aux); return aux;}       

   /** Get Column.
       This function returns a column vector corresponding to the choosen
       column inside matrix, the numbering of the column is also logical not
       physical.
       \\Ex: vector<double> v; m.getCol(-1, v);*/
   void getCol(int j, vT &v) const {
      if (XSIZE(*this)==0 || YSIZE(*this)==0) {v.clear(); return;}
      if (j<STARTINGX(*this) || j>FINISHINGX(*this))
         REPORT_ERROR(1103,"getCol: Matrix subscript (j) greater than matrix dimension");

      v.resize(YSIZE(*this));
      STARTINGX(v)=STARTINGY(*this);
      for (int i=STARTINGY(*this); i<=FINISHINGY(*this); i++)
          VEC_ELEM(v,i)=MAT_ELEM(*this,i,j);
      v.setCol();
   }
   
   /** Return Column. The same as previous. */
   vT Col(int i) const {vT aux; getCol(i,aux); return aux;}       

   /** Set Row.
       This function sets a row vector corresponding to the choosen
       row inside matrix, the numbering of the rows is also logical not
       physical.
       \\Ex: m.setRow(-2,m.row(1));
       \\--> Copies row 1 in row -2*/
   void setRow(int i, const vT &v) {
      if (XSIZE(*this)==0 || YSIZE(*this)==0)
         REPORT_ERROR(1,"setRow: Target matrix is empty");
      if (i<yinit || i>=yinit+ydim)
         REPORT_ERROR(1103,"setRow: Matrix subscript (i) out of range");
      if (v.get_dim()!=xdim)
         REPORT_ERROR(1102,"setRow: Vector dimension different from matrix one");
      if (!v.isRow())
         REPORT_ERROR(1107,"setRow: Not a row vector in assignment");

      i=i-STARTINGY(*this);
      for (int j=0; j<XSIZE(*this); j++)
          DIRECT_MAT_ELEM(*this,i,j)=DIRECT_VEC_ELEM(v,j);
   }

   /** Set Column.
       This function sets a column vector corresponding to the choosen
       column inside matrix, the numbering of the column is also logical not
       physical.
       \\Ex: m.setCol(-1,(m.row(1)).transpose());
       \\--> Copies row 1 in column -1*/
   void setCol(int j, const vT &v) {
      if (XSIZE(*this)==0 || YSIZE(*this)==0)
         REPORT_ERROR(1,"setCol: Target matrix is empty");
      if (j<xinit || j>=xinit+xdim)
         REPORT_ERROR(1103,"setCol: Matrix subscript (j) out of range");
      if (v.get_dim()!=ydim)
         REPORT_ERROR(1102,"setCol: Vector dimension different from matrix one");
      if (!v.isCol())
         REPORT_ERROR(1107,"setCol: Not a column vector in assignment");

      j=j-STARTINGX(*this);
      for (int i=0; i<YSIZE(*this); i++)
          DIRECT_MAT_ELEM(*this,i,j)=DIRECT_VEC_ELEM(v,i);
   }
   //@}

   /* Other utilities ------------------------------------------------------ */
   /**@name Utilities*/
   //@{
   /** This function must take two vectors of the same size, and operate
      element by element according to the operation required. This is the
      function which really implements the operations. Simple calls to it
      perform much faster than calls to the corresponding operators.
      Although it is supposed to be a hidden function not useable by
      normal programmers.
      It must be implemented in every Matrix module, this is so because
      of the Matrix2D, for which the multiplication is not a component
      by component multiplication but an algebraic one.*/
   friend void array_by_array(const maT &op1, const maT &op2, maT &result,
      char operation) {
         if (operation=='*') mul_matrix(op1, op2, result);
         else {
            if (operation=='x') operation='*';
            if (!op1.same_shape(op2))
               REPORT_ERROR(1007,
                  (string)"Array_by_array: different shapes ("+operation+")");
            result.resize(op1);
            core_array_by_array(op1, op2, result, operation);
         }
      }

   /** Algebraic transpose of matrix.
       You can use the transpose in as complex expressions as you like. The
       origin of the matrix is not changed.
       \\Ex: m2=m1.transpose();*/
   mT transpose() const {
      T aux;
      mT result(XSIZE(*this),YSIZE(*this));
      FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(result)
         DIRECT_MAT_ELEM(result,i,j)=DIRECT_MAT_ELEM(*this,j,i);
      STARTINGX(result)=STARTINGX(*this);
      STARTINGY(result)=STARTINGY(*this);
      return result;
   }
   
   /** Reverse matrix values over X axis.
       Maybe better with an example:

       \begin{verbatim}
         [1 2 3          [7 8 9
          4 5 6  ----->   4 5 6
          7 8 9]          1 2 3]
       \end{verbatim}
       \\Ex: m2=m1.reverseX();*/
   mT reverseX() const
      {mT temp(*this); temp.self_reverseX(); return temp;}
   
   /** Reverse matrix values over X axis, keep in this object. */
   void self_reverseX() {
      T aux;
      int jmax=(int)(XSIZE(*this)-1)/2;
      for (int i=0; i<YSIZE(*this); i++)
         for (int j=0; j<=jmax; j++) {
            SWAP(DIRECT_MAT_ELEM(*this,i,j),
               DIRECT_MAT_ELEM(*this,i,XSIZE(*this)-1-j), aux);
         }
   }
   
   /** Reverse matrix values over Y axis.
       Maybe better with an example:

       \begin{verbatim}
         [1 2 3          [3 2 1
          4 5 6  ----->   6 5 4
          7 8 9]          9 8 7]
       \end{verbatim}
       \\Ex: m2=m1.reverseY();*/
   mT reverseY() const
      {mT temp(*this); temp.self_reverseY(); return temp;}

   /** Reverse matrix values over Y axis, keep in this object. */
   void self_reverseY() {
      T aux;
      int imax=(int)(YSIZE(*this)-1)/2;
      for (int i=0; i<=imax; i++)
         for (int j=0; j<xdim; j++) {
            SWAP(DIRECT_MAT_ELEM(*this,i,j),
               DIRECT_MAT_ELEM(*this,YSIZE(*this)-1-i,j), aux);
         }
   }
   
   /** Determinant of a matrix.
       An exception is thrown if the matrix is not squared or it is empty.
       \\Ex: double det=m.det();*/
   T det() const {
      // (see Numerical Recipes, Chapter 2 Section 5)
      if (xdim==0)
         REPORT_ERROR(1108, "determinant: Matrix is empty");
      if (xdim!=ydim)
         REPORT_ERROR(1109, "determinant: Matrix is not squared");

      for (int i=0; i<YSIZE(*this); i++) {
         bool all_zeros=true;
         for (int j=0; j<XSIZE(*this); j++)
            if (ABS(DIRECT_MAT_ELEM(*this,i,j))>XMIPP_EQUAL_ACCURACY) {
	       all_zeros=false; break;
	    }
         if (all_zeros) return 0;
      }

      // Perform decomposition
      matrix1D<int> indx;
      T d;
      mT LU;
      ludcmp(*this, LU, indx, d);

      // Calculate determinant
      for (int i=0; i<XSIZE(*this); i++) d *= (T)LU(i,i);
      return d;
   }

   /** Inverse of a matrix.
       The matrix is inverted using a SVD decomposition.
       In fact the pseudoinverse is returned.
       \\Ex: matrix2D<double> m1_inv; m1.inv(m1_inv);*/
   void inv(mT &result) const  {
      if (xdim==0)
         REPORT_ERROR(1108, "Inverse: Matrix is empty");

      // Perform SVD decomposition
      matrix2D<double> u,v;
      matrix1D<double> w;
      svdcmp(*this,u,w,v); //*this=U*W*V^t

      double tol=compute_max()*MAX(XSIZE(*this),YSIZE(*this))*1e-14;
      result.init_zeros(XSIZE(*this),YSIZE(*this));

      // Compute W^-1
      bool invertible=false;
      for (int i=0; i<XSIZE(w); i++) {
         if (ABS(DIRECT_VEC_ELEM(w,i))>tol) {
            DIRECT_VEC_ELEM(w,i)=1.0/DIRECT_VEC_ELEM(w,i);
            invertible=true;
         } else DIRECT_VEC_ELEM(w,i)=0.0;
      }
      if (!invertible) return;

      // Compute V*W^-1
      FOR_ALL_ELEMENTS_IN_MATRIX2D(v)
         DIRECT_MAT_ELEM(v,i,j)*=DIRECT_VEC_ELEM(w,j);

      // Compute Inverse
      for (int i=0; i<XSIZE(*this); i++)
         for (int j=0; j<YSIZE(*this); j++)
            for (int k=0; k<XSIZE(*this); k++)
               DIRECT_MAT_ELEM(result,i,j)+= (T)(
                  DIRECT_MAT_ELEM(v,i,k)*DIRECT_MAT_ELEM(u,j,k));
   }

   /** Inverse of a matrix.*/
   mT inv() const {mT result; inv(result); return result;}

   /** Solve equation system.
       The equation system is defined by Ax=b, it is solved for x.
       Exceptions are thrown if the equation system is not solvable.
       \\Ex: matrix1D<double> x=m.inv();*/
   friend void solve<>(const mT &A, const vT &b, vT &result);

   /** Solve equation system.
       The same as before but now, b is a matrix and so, x is also a matrix.
       A must be a square matrix.*/ 
   friend void solve<>(const mT &A, const mT &b, mT &result);

   /** Put a window to matrix.
       The matrix is windowed within the two positions given to this function.
       Indexes always refer to logical indexes. If a position is outside the
       actual matrix range then the matrix is padded with init_value until the
       new position is reached. In the following examples suppose that m1
       is the following and that the origin is (-1,-1).

       \begin{verbatim}
               [1 2 3               [1 2 3 0
          m1 =  4 5 6    --->  m1 =  4 5 6 0
                7 8 9]               7 8 9 0]
       \end{verbatim}
       \\Ex: m1.window(-1,-1,1,2);
       */
   void window(int y0, int x0, int yF, int xF, T init_value=0) {
      mT result(yF-y0+1, xF-x0+1);
      STARTINGY(result)=y0;
      STARTINGX(result)=x0;

      FOR_ALL_ELEMENTS_IN_MATRIX2D(result)
         if (j>=STARTINGX(*this) && j<=FINISHINGX(*this) &&
             i>=STARTINGY(*this) && i<=FINISHINGY(*this))
             MAT_ELEM(result,i,j)=MAT_ELEM(*this,i,j);
          else
              MAT_ELEM(result,i,j)=init_value;	  
      *this=result;
   }
   //@}

   /* Geometrical transformations ------------------------------------------ */
   /**@name Geometrical Transformations
      In all geometrical transformations a periodic extension of the matrix
      is supposed, ie, if a pixel goes out on the left, it is entering on
      the right, ...*/
   //@{
   /** Applies a geometrical transformation.
       Any geometrical transformation defined by the matrix A (double (3x3)!!
       ie, in homogeneous R2 coordinates) is applied to the matrix M1.
       The result is stored in m2 (it cannot be the same as the input
       matrix). An exception is thrown if the
       transformation matrix is not 3x3.
       The result matrix is resized to the same dimensions as M1 if M2 is empty
       (0x0) at the beginning, if it is not, ie, if M2 has got some size
       then only those values in the matrix are filled, this is very
       useful for resizing the matrix, then you manually resize the output
       matrix to the desired size and then call this routine.

       The relationship between the output coordinates and the input ones are
       \begin{verbatim}
           out = A * in
         (x,y) = A * (x',y')
       \end{verbatim}

       This function works independently from the logical indexing of each
       matrix, it sets the logical center and the physical center of the image
       and work with these 2 coordinate spaces. At the end the original logical
       indexing of each matrix is kept.

       The procedure followed goes from coordinates in the output matrix
       to the ones in the input one, so the inverse of the A matrix is
       needed. There is a flag telling if the given transformation matrix
       is already
       the inverse one or the normal one. If it is the normal one internally
       the matrix is inversed. If you are to do many "rotations" then
       some time is spent in inverting the matrix. Normally the
       matrix is the normal one.
       
       There is something else to tell about the geometrical tranformation.
       The value of the pixel in the output matrix is computed via
       bilinear interpolation in the input matrix. If any of the pixels
       participating in the interpolation falls outside the input matrix,
       then automatically the corresponding output pixel is set to 0, unless
       that the wrap flag has been set to 1. In this case if the pixel
       falls out by the right hand then it is "wrapped" and the corresponding
       pixel in the left hand is used. The same is appliable to top-bottom.
       Usually wrap mode is off. Wrap mode is interesting for translations
       but not for rotations, for example.
       
       The inverse mode and wrapping mode should be taken by default by the
       routine, g++ seems to have problems with template functions outside
       a class with default parameters. So, I'm sorry, you will have to
       put them always. The usual combination is
       apply_geom(...,IS_NOT_INV,DONT_WRAP). Although you can also use the
       constants IS_INV, or WRAP.
       
       m2 cannot be the same matrix as m1.
       \\Ex: matrix2D<double> A(3,3); A.init_identity; apply_geom(m2,A,m1);*/
   friend void apply_geom<>(mT &m2, matrix2D<double> A,
       const mT &m1, bool inv, bool wrap, T outside);

   /** Apply geom with B-spline interpolation. */
   friend void apply_geom_Bspline<>(mT &m2, matrix2D<double> A,
       const mT &m1, int Splinedegree, bool inv, bool wrap, T outside);

   /** Self apply geom.
       Same as the previous one, but the result is kept in this object */
   void self_apply_geom(matrix2D<double> A, bool inv, bool wrap, T outside=(T)0)
      {mT aux; apply_geom(aux,A,*this,inv,wrap,outside); *this=aux;}

   /** Self apply geom with Bspline interpolation.*/
   void self_apply_geom_Bspline(matrix2D<double> A, int SplineDegree,
      bool inv, bool wrap, T outside=(T)0)
      {mT aux; apply_geom_Bspline(aux,A,*this,SplineDegree,inv,wrap,outside);
       *this=aux;}

   #define IS_INV     true
   #define IS_NOT_INV false
   #define DONT_WRAP  false
   #define WRAP       true
   
   /** Rotate matrix.
       The angle must be in degrees. The result cannot be this object.
       \\Ex: m1.rotate(60,m2);*/
   void rotate(double ang, mT &result, bool wrap=DONT_WRAP) const
      {matrix2D<double> temp=rot2D_matrix(ang);
       apply_geom(result,temp,*this,IS_NOT_INV,wrap);}
   
   /** Rotate matrix.
       Same as the previous one. */
   mT rotate(double ang, bool wrap=DONT_WRAP) const
      {mT aux; rotate(ang,aux,wrap); return aux;}

   /** Rotate matrix.
       Same as the previous one, but the result is kept in this object */
   void self_rotate(double ang, bool wrap=DONT_WRAP)
      {mT aux; rotate(ang,aux,wrap); *this=aux;}

   /** Rotate matrix (using Bspline interpolation).
       The angle must be in degrees. The result cannot be this object.
       \\Ex: m1.rotate(60,m2);*/
   void rotate_Bspline(int Splinedegree,
      double ang, mT &result, bool wrap=DONT_WRAP) const
      {matrix2D<double> temp=rot2D_matrix(ang);
       apply_geom_Bspline(result,temp,*this,Splinedegree,IS_NOT_INV,wrap);}
   
   /** Rotate matrix  (using Bspline interpolation).
       Same as the previous one. */
   mT rotate_Bspline(int Splinedegree, double ang, bool wrap=DONT_WRAP) const
      {mT aux; rotate_Bspline(Splinedegree,ang,aux,wrap); return aux;}

   /** Rotate matrix  (using Bspline interpolation).
       Same as the previous one, but the result is kept in this object */
   void self_rotate_Bspline(int Splinedegree, double ang, bool wrap=DONT_WRAP)
      {mT aux; rotate_Bspline(Splinedegree,ang,aux,wrap); *this=aux;}

   /** Translate matrix.
       The displacement is given as a R2 vector of the form (shift_X,shift_Y).
       The result cannot be this object.
       \\Ex: m2=m1.translate(vector_R2(0,2));
       \\--> m1 is shifted 2 pixels down and stored in m2*/
   void translate(const matrix1D<double> &v, mT &result, bool wrap=WRAP) const
      {matrix2D<double> temp=translation2D_matrix(v);
       apply_geom(result,temp,*this,IS_NOT_INV,wrap);}
   
   /** Translate matrix.
       Same as the previous one. */
   mT translate(const matrix1D<double> &v, bool wrap=WRAP) const
      {mT aux; translate(v,aux,wrap); return aux;}
   
   /** Translate matrix.
       Same as the previous one, but the result is kept in this object */
   void self_translate(const matrix1D<double> &v, bool wrap=WRAP)
      {mT aux; translate(v,aux,wrap); *this=aux;}

   /** Translate matrix (using Bspline interpolation).
       The displacement is given as a R2 vector of the form (shift_X,shift_Y).
       The result cannot be this object.
       \\Ex: m2=m1.translate(vector_R2(0,2));
       \\--> m1 is shifted 2 pixels down and stored in m2*/
   void translate_Bspline(int Splinedegree,
      const matrix1D<double> &v, mT &result, bool wrap=WRAP) const
      {matrix2D<double> temp=translation2D_matrix(v);
       apply_geom_Bspline(result,temp,*this,Splinedegree,IS_NOT_INV,wrap);}
   
   /** Translate matrix (using Bspline interpolation).
       Same as the previous one. */
   mT translate_Bspline(int Splinedegree,
      const matrix1D<double> &v, bool wrap=WRAP) const
      {mT aux; translate_Bspline(Splinedegree,v,aux,wrap); return aux;}
   
   /** Translate matrix (using Bspline interpolation).
       Same as the previous one, but the result is kept in this object */
   void self_translate_Bspline(int Splinedegree,
      const matrix1D<double> &v, bool wrap=WRAP)
      {mT aux; translate_Bspline(Splinedegree,v,aux,wrap); *this=aux;}

   /** Translate center of mass to center.
       If the input has very high values, it is better to rescale it to
       be between 0 and 1. */
   void self_translate_center_of_mass_to_center(bool wrap=WRAP) {
      set_Xmipp_origin();
      matrix1D<double> center;
      center_of_mass(center);
      center*=-1;
      self_translate(center,wrap);
   }

   /** Translate center of mass to center (using Bspline interpolation).
       If the input has very high values, it is better to rescale it to
       be between 0 and 1. */
   void self_translate_center_of_mass_to_center_Bspline(
      int Splinedegree, bool wrap=WRAP) {
      set_Xmipp_origin();
      matrix1D<double> center;
      center_of_mass(center);
      center*=-1;
      self_translate_Bspline(Splinedegree,center,wrap);
   }

   /** Scales to a new size.
       The matrix is scaled (resampled) to fill a new size. It is not the
       same as "window" in this same class. The size can be larger or 
       smaller than the actual one. But the result matrix cannot be
       this object.
       \\Ex: m1.scale_to_size(128,128);*/
   void scale_to_size(int Ydim,int Xdim, mT &result) const {
      matrix2D<double> temp(3,3);
      result.resize(Ydim,Xdim);
      temp.init_identity();
      DIRECT_MAT_ELEM(temp,0,0)=(double)Xdim/(double)XSIZE(*this);
      DIRECT_MAT_ELEM(temp,1,1)=(double)Ydim/(double)YSIZE(*this);
      apply_geom(result,temp,*this,IS_NOT_INV,WRAP);
   }

   /** Scales to a new size.
       Same as the previous one. */
   mT scale_to_size(int Ydim,int Xdim) const
      {mT aux; scale_to_size(Ydim, Xdim, aux); return aux;}

   /** Scales to a new size.
       Same as the previous one, but the result is kept in this object. */
   void self_scale_to_size(int Ydim,int Xdim)
      {mT aux; scale_to_size(Ydim, Xdim, aux); *this=aux;}

   /** Scales to a new size (using Bspline interpolation).
       The matrix is scaled (resampled) to fill a new size. It is not the
       same as "window" in this same class. The size can be larger or 
       smaller than the actual one. But the result matrix cannot be
       this object.
       \\Ex: m1.scale_to_size(128,128);*/
   void scale_to_size_Bspline(int Splinedegree,
      int Ydim,int Xdim, mT &result) const {
      matrix2D<double> temp(3,3);
      result.resize(Ydim,Xdim);
      temp.init_identity();
      DIRECT_MAT_ELEM(temp,0,0)=(double)Xdim/(double)XSIZE(*this);
      DIRECT_MAT_ELEM(temp,1,1)=(double)Ydim/(double)YSIZE(*this);
      apply_geom_Bspline(result,temp,*this,Splinedegree,IS_NOT_INV,WRAP);
   }

   /** Scales to a new size (using Bspline interpolation).
       Same as the previous one. */
   mT scale_to_size_Bspline(int Splinedegree,int Ydim,int Xdim) const
      {mT aux; scale_to_size_Bspline(Splinedegree,Ydim, Xdim, aux); return aux;}

   /** Scales to a new size (using Bspline interpolation).
       Same as the previous one, but the result is kept in this object. */
   void self_scale_to_size_Bspline(int Splinedegree,int Ydim,int Xdim)
      {mT aux; scale_to_size_Bspline(Splinedegree,Ydim, Xdim, aux); *this=aux;}

   /** Superpixel reduce.
       This function reduces the given image averaging in the superpixel
       of the given size. The last columns and rows are removed if the
       size of the original image is not an exact multiple of the given
       superpixel size */
   void superpixel_reduce(mT &result, int size=2) const {
      result.init_zeros(YSIZE(*this)/size,XSIZE(*this)/size);
      int size2=size*size;
      FOR_ALL_ELEMENTS_IN_MATRIX2D(result) {
         for (int ii=0; ii<size; ii++)
	    for (int jj=0; jj<size; jj++)
	       DIRECT_MAT_ELEM(result,i,j)+=
	          DIRECT_MAT_ELEM(*this,size*i+ii,size*j+jj);
      	 DIRECT_MAT_ELEM(result,i,j)/=size2;
      }
   }

   /** Superpixel expand.
       This function copies each pixel to a new image as many times
       as the size of the superpixel. */
   void superpixel_expand(mT &result, int size=2) const {
      result.init_zeros(YSIZE(*this)*size,XSIZE(*this)*size);
      FOR_ALL_ELEMENTS_IN_MATRIX2D(*this) {
         for (int ii=0; ii<size; ii++)
	    for (int jj=0; jj<size; jj++)
	       DIRECT_MAT_ELEM(result,size*i+ii,size*j+jj)=
	          DIRECT_MAT_ELEM(*this,i,j);
      }
   }

   /** Reduce the image by 2 using a BSpline pyramid. */
   void pyramid_reduce(matrix2D<double> &result, int levels=1) const {
      matrix2D<double> aux, aux2;
      produce_spline_coeffs(aux,3);
      for (int i=0; i<levels; i++) {
         aux.reduce_Bspline(aux2,3);
         aux=aux2;
      }
      aux2.produce_image_from_spline_coeffs(result,3);
   }

   /** Expand the image by 2 using a BSpline pyramid. */
   void pyramid_expand(matrix2D<double> &result, int levels=1) const {
      matrix2D<double> aux, aux2;
      produce_spline_coeffs(aux,3);
      cout << levels << endl;
      for (int i=0; i<levels; i++) {
         aux.expand_Bspline(aux2,3);
         aux=aux2;
      }
      aux2.produce_image_from_spline_coeffs(result,3);
   }

   #ifndef DBL_EPSILON
      #define DBL_EPSILON 1e-50
   #endif
   /** Produce spline coefficients.*/
   void produce_spline_coeffs(matrix2D<double> &coeffs, int SplineDegree=3)
      const {
      coeffs.init_zeros(YSIZE(*this),XSIZE(*this));
      STARTINGX(coeffs)=STARTINGX(*this);
      STARTINGY(coeffs)=STARTINGY(*this);
      int Status;
      matrix2D<double> aux;
      type_cast(*this,aux);
      ChangeBasisVolume(MULTIDIM_ARRAY(aux),MULTIDIM_ARRAY(coeffs),
          XSIZE(*this), YSIZE(*this), 1,
          CardinalSpline, BasicSpline, SplineDegree,
          MirrorOffBounds, DBL_EPSILON, &Status);
      if (Status)
         REPORT_ERROR(1,"matrix2D::produce_spline_coeffs: Error");
   }

   /** Produce image from B-spline coefficients. */
   void produce_image_from_spline_coeffs(
      matrix2D<double> &img, int SplineDegree=3) const {
      img.init_zeros(YSIZE(*this),XSIZE(*this));
      STARTINGX(img)=STARTINGX(*this);
      STARTINGY(img)=STARTINGY(*this);
      int Status;
      matrix2D<double> aux;
      type_cast(*this,aux);
      ChangeBasisVolume(MULTIDIM_ARRAY(aux),MULTIDIM_ARRAY(img),
          XSIZE(*this), YSIZE(*this), 1,
          BasicSpline, CardinalSpline, SplineDegree,
          MirrorOnBounds, DBL_EPSILON, &Status);
      if (Status)
         REPORT_ERROR(1,"matrix2D::produce_spline_img: Error");
   }
#undef DBL_EPSILON

   /** Expand a set of B-spline coefficients.
       Knowing that this matrix is a set of B-spline coefficients,
       produce the expanded set of B-spline coefficients using the
       two-scale relationship. */
   void expand_Bspline(matrix2D<double> &expanded, int SplineDegree=3) const {
      double   g[200];    /* Coefficients of the reduce filter */
      long     ng;         /* Number of coefficients of the reduce filter */
      double   h[200];    /* Coefficients of the expansion filter */
      long     nh;         /* Number of coefficients of the expansion filter */
      short    IsCentered; /* Equal TRUE if the filter is a centered spline, FALSE otherwise */

      // Get the filter
      if (GetPyramidFilter( "Centered Spline", SplineDegree,
            g, &ng, h, &nh, &IsCentered))
	 REPORT_ERROR(1,"Unable to load the filter coefficients");

      matrix2D<double> aux;
      type_cast(*this,aux);
      expanded.resize(2*YSIZE(aux),2*XSIZE(aux));
      Expand_2D(MULTIDIM_ARRAY(aux), XSIZE(aux), YSIZE(aux),
         MULTIDIM_ARRAY(expanded), h, nh, IsCentered);
   }

   /** Reduce a set of B-spline coefficients.
       Knowing that this matrix is a set of B-spline coefficients,
       produce the reduced set of B-spline coefficients using the
       two-scale relationship. */
   void reduce_Bspline(matrix2D<double> &reduced, int SplineDegree=3) const {
      double   g[200];    /* Coefficients of the reduce filter */
      long     ng;         /* Number of coefficients of the reduce filter */
      double   h[200];    /* Coefficients of the expansion filter */
      long     nh;         /* Number of coefficients of the expansion filter */
      short    IsCentered; /* Equal TRUE if the filter is a centered spline, FALSE otherwise */

      // Get the filter
      if (GetPyramidFilter( "Centered Spline", SplineDegree,
            g, &ng, h, &nh, &IsCentered))
	 REPORT_ERROR(1,"Unable to load the filter coefficients");

      matrix2D<double> aux;
      type_cast(*this,aux);
      if (XSIZE(aux)%2!=0 && YSIZE(aux)%2!=0)
          aux.resize(YSIZE(aux)-1,XSIZE(aux)-1);
      else if (YSIZE(aux)%2!=0)
          aux.resize(YSIZE(aux)-1,XSIZE(aux));
      else if (XSIZE(aux)%2!=0)
          aux.resize(YSIZE(aux)  ,XSIZE(aux)-1);
      reduced.resize(YSIZE(aux)/2,XSIZE(aux)/2);
      Reduce_2D(MULTIDIM_ARRAY(aux), XSIZE(aux), YSIZE(aux),
         MULTIDIM_ARRAY(reduced), g, ng, IsCentered);
   }

   //@}

   /* Iterators ------------------------------------------------------------ */
   /**@name Iterators*/
   //@{
   /** Apply the same scalar function to all rows.
       This function must take a row vector and return a single value, a column
       vector with these values is returned.
       \\Ex:T vector_sum(vT &v) {return v.sum();};
            v1=m.for_all_rows(&vector_sum);*/
   vT for_all_rows (T (*f)(vT&)) const  {
      vT temp;
      if (XSIZE(*this)==0 || YSIZE(*this)==0) return temp;
      temp.resize(YSIZE(*this));
      STARTINGX(temp)=STARTINGY(*this);
      temp.setCol();
      for (int i=STARTINGY(*this); i<=FINISHINGY(*this); i++) {
         vT aux;
         getRow(i,aux);
         VEC_ELEM(temp,i)=(*f)(aux);
      }
      return temp;
   }

   /** Apply the same scalar function to all columns.
       This function must take a column vector and return a single value, a row
       vector with these values is returned.
       \\Ex:T vector_sum(vT &v) {return v.sum();};
            v1=m.for_all_cols(&vector_sum);*/
   vT for_all_cols (T (*f)(vT&)) const {
      vT temp;
      if (XSIZE(*this)==0 || YSIZE(*this)==0) return temp;
      temp.resize(XSIZE(*this));
      STARTINGX(temp)=STARTINGX(*this);
      temp.setRow();
      for (int j=STARTINGX(*this); j<=FINISHINGX(*this); j++) {
         vT aux;
         getCol(j,aux);
         VEC_ELEM(temp,j)=(*f)(aux);
      }
      return temp;
   }
   
   /** Apply the same vectorial function to all rows.
       This function must take a row vector and return a row vector (of the
       same size as the input one), a new
       matrix with these transformed rows is returned.
       \\Ex:vT vector_norm(vT &v) {return v.normalize();};
            m2=m.for_all_rows(&vector_norm);*/
   mT for_all_rows (vT (*f)(vT&)) const
      {mT aux(*this); aux.for_all_rows(f); return aux;}

   /** Apply a vectorial function to all rows, keep in this object. */
   void for_all_rows (vT (*f)(vT&)) {
      if (XSIZE(*this)==0 || YSIZE(*this)==0) return;
      for (int i=STARTINGY(*this); i<=FINISHINGY(*this); i++) {
         vT aux;
         getRow(i,aux);
         setRow(i,(*f)(aux));
      }
   }

   /** Apply the same vectorial function to all columns.
       This function must take a columnm vector and return a column vector
       (of the same size as the input one), a
       new matrix with these transformed columns is returned.
       \\Ex:vT vector_norm(vT &v) {return v.normalize();};
            m2=m.for_all_cols(&vector_norm);*/
   mT for_all_cols (vT (*f)(vT&)) const
      {mT aux(*this); aux.for_all_cols(f); return aux;}

   /** Apply a vectorial function to all rows, keep in this object. */
   void for_all_cols (vT (*f)(vT&)) {
      if (XSIZE(*this)==0 || YSIZE(*this)==0) return;
      for (int j=STARTINGX(*this); j<=FINISHINGX(*this); j++) {
         vT aux;
         getCol(j,aux);
         setCol(j,(*f)(aux));
      }
   }

   //@}

   /* Operation related ---------------------------------------------------- */
   /**@name Algebraic operations
      NOTICE!!!: the matrix by matrix multiplier operator (*) has been
      redefined in this class and it represents the true matrix by
      matrix algebraic multiplication, if you want an element-wise
      multiplication you have to use the function defined in this
      section "mul_elements" */
   //@{
   /** Matrix multiplication element by element.
      \\Ex: m3=mul_elements(m1,m2); */
   friend void mul_elements(const mT &op1, const mT &op2, mT &result)
      {array_by_array(op1,op2,result,'x');}

   /** The same as before but the result is returned */
   friend mT mul_elements(const mT &op1, const mT &op2)
      {mT temp; mul_elements(op1,op2,temp); return temp;}

   /** Matrix by vector multiplication.
      \\Ex: v2=A*v1;*/
   vT operator * (const vT &op1) const  {
      vT result;
      if (XSIZE(*this)!=XSIZE(op1))
         REPORT_ERROR(1102,"Not compatible sizes in matrix by vector");
      if (!op1.isCol())
         REPORT_ERROR(1102,"Vector is not a column");

      result.init_zeros(YSIZE(*this));
      for (int i=0; i<YSIZE(*this); i++)
         for (int j=0; j<XSIZE(op1); j++)
               DIRECT_VEC_ELEM(result,i) += DIRECT_MAT_ELEM(*this,i,j)*
                  DIRECT_VEC_ELEM(op1,j);
      result.setCol();
      STARTINGX(result)=STARTINGY(*this);
      return result;
   }
   //@}

   /* Matrix types --------------------------------------------------------- */
   /**@name Matrix types*/
   //@{
   /** True if matrix is a single row or a single column.
       \\Ex: if (m.IsVector()) cout << "The matrix is like a vector\n";*/
   bool IsVector () const {return (XSIZE(*this)==1) || (YSIZE(*this)==1);}
   
   /** True if the matrix is square.
       \\Ex: if (m.IsSquare()) cout << "The matrix is square\n";*/
   bool IsSquare () const {return (XSIZE(*this)==YSIZE(*this));}
   
   /** True if the matrix is singular (det()==0).
       \\Ex: if (m.IsSingular()) cout << "The matrix is singular\n";*/
   bool IsSingular() const {return ABS(det())<XMIPP_EQUAL_ACCURACY;}
   
   /** True if the matrix is diagonal.
       \\Ex: if (m.IsDiagonal()) cout << "The matrix is diagonal\n";*/
   bool IsDiagonal() const
      {if (XSIZE(*this)!=YSIZE(*this)) return false;
       FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
             if (i!=j && ABS(DIRECT_MAT_ELEM(*this,i,j))>XMIPP_EQUAL_ACCURACY)
                return false;
       return true;}

   /** True if the matrix is scalar.
       A scalar matrix is diagonal and all the values at the diagonal are
       the same*/
   bool IsScalar() const
      {if (!IsDiagonal()) return false;
       for (int i=1; i<YSIZE(*this); i++)
          if (ABS(DIRECT_MAT_ELEM(*this,i,i)-DIRECT_MAT_ELEM(*this,0,0))>
              XMIPP_EQUAL_ACCURACY) return false;
       return true;}

   /** True if the matrix is symmetric.
       \\Ex: if (m.IsSymmetric()) cout << "The matrix is symmetric\n";*/
   bool IsSymmetric() const
      {if (XSIZE(*this)!=YSIZE(*this)) return false;
       for (int i=0; i<YSIZE(*this); i++)
          for (int j=i+1; j<XSIZE(*this); j++)
             if (ABS(DIRECT_MAT_ELEM(*this,i,j)-DIRECT_MAT_ELEM(*this,j,i))>
                 XMIPP_EQUAL_ACCURACY) return false;
       return true;}

   /** True if the matrix is skew-symmetric (anti-symmetric).
       \\Ex: if (m.IsSkewSymmetric()) cout << "The matrix is skewsymmetric\n";*/
   bool IsSkewSymmetric() const
      {if (XSIZE(*this)!=YSIZE(*this)) return false;
       for (int i=0; i<YSIZE(*this); i++)
          for (int j=i+1; j<XSIZE(*this); j++)
             if (ABS(DIRECT_MAT_ELEM(*this,i,j)+DIRECT_MAT_ELEM(*this,j,i))>
                 XMIPP_EQUAL_ACCURACY) return false;
       return true;}

   /** True if the matrix is upper-triangular.
       \\Ex: if (m.IsUpperTriangular()) cout << "The matrix is upper triangular\n";*/
   bool IsUpperTriangular() const
      {if (XSIZE(*this)!=YSIZE(*this)) return false;
       for (int i=1; i<YSIZE(*this); i++)
          for (int j=0; j<i-1; j++)
             if (ABS(DIRECT_MAT_ELEM(*this,i,j))>XMIPP_EQUAL_ACCURACY)
                return false;
       return true;}

   /** True if the matrix is lower-triangular.
       \\Ex: if (m.IsLowerTriangular()) cout << "The matrix is lower triangular\n";*/
   bool IsLowerTriangular() const
      {if (XSIZE(*this)!=YSIZE(*this)) return false;
       for (int i=1; i<YSIZE(*this); i++)
          for (int j=i+1; j<XSIZE(*this); j++)
             if (ABS(DIRECT_MAT_ELEM(*this,i,j))>XMIPP_EQUAL_ACCURACY)
                return false;
       return true;}

   /** True if the matrix is identity.
       \\Ex: if (m.IsIdent()) cout << "The matrix is identity\n";*/
   bool IsIdent() const {return IsScalar() &&
      ABS(DIRECT_MAT_ELEM(*this,0,0)-(T)1)<XMIPP_EQUAL_ACCURACY;}
   
   /** True if the matrix is null.
       \\Ex: if (m.IsZero()) cout << "The matrix is null\n";*/
   bool IsZero() const {return IsScalar() && 
      ABS(DIRECT_MAT_ELEM(*this,0,0)-(T)0)<XMIPP_EQUAL_ACCURACY;}

   /** Maximum element.
       This function returns the index of the maximum element of an array.
       array(i,j). Returns -1 if the array is empty*/
   void max_index(int &imax, int &jmax) const {
      if (XSIZE(*this)==0) {imax=jmax=-1; return;}
      imax=jmax=0;
      T   max=MAT_ELEM(*this,imax,jmax);
      FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
         if (MAT_ELEM(*this,i,j)>max)
            {max=MAT_ELEM(*this,i,j); imax=i; jmax=j;}
   }

   /** Minimum element.
       This function returns the index of the minimum element of an array.
       array(i,j). Returns -1 if the array is empty*/
   void min_index(int &imin, int &jmin) const {
      if (XSIZE(*this)==0) {imin=jmin=-1; return;}
      imin=jmin=0;
      T   min=MAT_ELEM(*this,imin,jmin);
      FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
         if (MAT_ELEM(*this,i,j)>min)
            {min=MAT_ELEM(*this,i,j); imin=i; jmin=j;}
   }
   //@}
};

#include "Src/MultidimFriends_implementation.hh"
// Special case for complex numbers
template <>
complex<double> matrix2D< complex<double> >::interpolated_elem(
   double x, double y, complex<double> outside_value) const;

template <>
void core_array_by_scalar< complex<double> >(const maTC &op1,
   const complex<double> &op2, maTC &result, char operation);

template <>
void core_scalar_by_array< complex<double> >(const complex<double> &op1,
   const maTC &op2, maTC &result, char operation);

template <>
void core_array_by_array< complex<double> >(const maTC &op1, const maTC &op2,
   maTC &result, char operation);

template <>
inline void matrix2D<complex<double> >::produce_spline_coeffs(
   matrix2D<double> &coeffs, int SplineDegree) const {
   // *** STILL TO DO
   cerr << "Spline coefficients of a complex matrix is not implemented\n";
}

/**@name Related functions
   These functions are not methods of matrix1D */
//@{

/* Geometry ---------------------------------------------------------------- */
/**@name Geometry with matrices */
//@{
/** Creates a rotational matrix (3x3) for images.
    The rotation angle is in degrees.
    \\Ex: m=rot2D_matrix(60);*/
matrix2D<double> rot2D_matrix(double ang);

/** Creates a translational matrix (3x3) for images.
    The shift is given as a R2 vector (shift_X, shift_Y);
    An exception is thrown if the displacement is not a R2 vector.
    \\Ex: m=translation2D_matrix(vector_R2(1,0));
    \\--> Displacement of 1 pixel to the right */
matrix2D<double> translation2D_matrix(const matrix1D<double> v);

/** Creates a rotational matrix (4x4) for volumes around system axis.
    The rotation angle is in degrees, and the rotational axis is
    either 'X', 'Y' or 'Z'. An exception is thrown if the axis given
    is not one of these.
    The returned matrices are respectively
    alpha degrees around Z
    \begin{verbatim}
    [ cos(A) -sin(A)     0   ]
    [ sin(A)  cos(A)     0   ]
    [   0       0        1   ]
    \end{verbatim}
    alpha degrees around Y
    \begin{verbatim}
    [ cos(A)    0    -sin(A) ]
    [   0       1       0    ]
    [ sin(A)    0     cos(A) ]
    \end{verbatim}
    alpha degrees around X
    \begin{verbatim}
    [   1       0       0    ]
    [   0     cos(A) -sin(A) ]
    [   0     sin(A)  cos(A) ]
    \end{verbatim}
    \\Ex: m=rot3D_matrix(60,'X');*/
matrix2D<double> rot3D_matrix(double ang, char axis);

/** Creates a rotational matrix (4x4) for volumes around any axis.
    The rotation angle is in degrees, and the rotational axis is
    given as a R3 vector. An exception is thrown if the axis is not a
    R3 vector. The axis needs not to be unitary.
    \\Ex: m=rot3D_matrix(60,vector_R3(1,1,1));*/
matrix2D<double> rot3D_matrix(double ang, const matrix1D<double> &axis);

/** Matrix which transforms the given axis into Z.
    A geometrical transformation matrix (4x4) is returned such that the
    given axis is rotated until it is aligned with the Z axis. This is
    very useful in order to produce rotational matrices, for instance,
    around any axis.
    \\Ex:
    \begin{verbatim}
    matrix2D<double> A=align_with_Z(axis);
    return A.transpose() * rot3D_matrix(ang,'Z') * A;
    \end{verbatim}
    
    The returned matrix is such that A*axis=Z, where Z and axis are
    column vectors.*/
matrix2D<double> align_with_Z(const matrix1D<double> &axis);

/** Creates a translational matrix (4x4) for volumes.
    The shift is given as a R3 vector (shift_X, shift_Y, shift_Z);
    An exception is thrown if the displacement is not a R3 vector.
    \\Ex: m=translation3D_matrix(vector_R3(0,0,2));
    \\--> Displacement of 2 pixels down */
matrix2D<double> translation3D_matrix(const matrix1D<double> &v);

/** Creates a scaling matrix (4x4) for volumes.
    The scaling factors for the different axis must be given as a vector.
    So that, XX(sc)=scale for X axis, YY(sc)=...*/
matrix2D<double> scale3D_matrix(const matrix1D<double> &sc);

//@}



/* Other useful functions -------------------------------------------------- */
/**@name Miscellaneous*/
//@{
/** Reduce both matrices to a common size.
    Search the range of logical indexes for which both matrices have got
    valid values, and cut both to that size, the corresponding
    origin is automatically computed.
    \\Ex: matrix2D<double> m1(5,3); m1.startingX()=-2; m1.startingY()=-2;
    \\matrix2D<double> m2(4,4); m2.startingX()=0; m2.startingY()=0;
    \\cut_to_common_size(m1,m2);
    \\--> m1 and m2 range from (0,0)=(y,x) to (2,0) */
template <class T>
   void cut_to_common_size(mT &m1, mT &m2) {
      int y0=MAX(STARTINGY(m1) ,STARTINGY(m2));
      int yF=MIN(FINISHINGY(m1),FINISHINGY(m2));
      int x0=MAX(STARTINGX(m1) ,STARTINGX(m2));
      int xF=MIN(FINISHINGX(m1),FINISHINGX(m2));
      m1.window(y0,x0,yF,xF);
      m2.window(y0,x0,yF,xF);
   }

/** Does a radial average of a matrix, around the pixel where is the origin.
    A vector is returned where:
	 - the first element is the mean of the pixels whose
	   distance to the origin is (0-1),
	 - the second element is the mean of the pixels
       whose distance to the origin is (1-2)
	 - and so on.
    A second vector radial_count is returned containing the number of
    pixels over which each radial average was calculated.
    if rounding=true, element=round(distance);
         - so the first element is the mean of the voxels whose
           distance to the origin is (0.5-1.5),
         - the second element is the mean of the voxels
       whose distance to the origin is (1.5-2.5)
         - and so on. */
template <class T>
   void radial_average(const matrix2D<T> &m, const matrix1D<int> &center_of_rot,
      matrix1D<T> &radial_mean, matrix1D<int> &radial_count,
      const bool &rounding=false) {
      matrix1D<double> idx(2);

      /* First determine the maximum distance that one should expect,
         to set the dimension of the radial average vector */
      matrix1D<int> distances(4);
      double y=STARTINGY(m)-YY(center_of_rot);
      double x=STARTINGX(m)-XX(center_of_rot);
      distances(0)=(int)floor(sqrt(x*x+y*y));
      x=FINISHINGX(m)-XX(center_of_rot);
      y=STARTINGY(m)-YY(center_of_rot);
      distances(1)=(int)floor(sqrt(x*x+y*y));
      x=STARTINGX(m)-XX(center_of_rot);
      y=FINISHINGY(m)-YY(center_of_rot);
      distances(2)=(int)floor(sqrt(x*x+y*y));
      x=FINISHINGX(m)-XX(center_of_rot);
      y=FINISHINGY(m)-YY(center_of_rot);
      distances(3)=(int)floor(sqrt(x*x+y*y));
      int dim=(int)CEIL(distances.compute_max())+1;
      if (rounding) dim++;

      // Define the vectors
      radial_mean.resize(dim);
      radial_mean.init_zeros();
      radial_count.resize(dim);
      radial_count.init_zeros();   

      /* Perform the radial sum and count pixels that contribute to
         every distance */
      FOR_ALL_ELEMENTS_IN_MATRIX2D(m)
      {
         YY(idx)=i-YY(center_of_rot);
	     XX(idx)=j-XX(center_of_rot);
	     // Determine distance to the center
	     int distance;
	     if (rounding) distance=(int)ROUND(idx.module());
	     else distance=(int)floor(idx.module());

         // Sum te value to the pixels with the same distance  
	     radial_mean(distance)+=m(i,j);
	     // Count the pixel
	     radial_count(distance)++;      	  	  
      }

      // Perfor the mean
      FOR_ALL_ELEMENTS_IN_MATRIX1D(radial_mean)
      {
         radial_mean(i)/=(T)radial_count(i);
      }
   }


/** Solve equation system, nonnegative solution.
    The equation system is defined by Ax=b, it is solved for x.
    x is forced to be nonnegative. It is designed to cope with
    large equation systems. This function is borrowed from
    LAPACK nnls.
    
    The norm of the vector Ax-b is returned.*/
double solve_nonneg(const matrix2D<double> &A, const matrix1D<double> &b,
   matrix1D<double> &result);

/** Solve equation system, symmetric positive-definite matrix.
    The equation system is defined by Ax=b, it is solved for x.
    This method can only be applied if A is positive-definite matrix
    and symmetric. It applies a Cholesky factorization and 
    backsubstitution (see Numerical Recipes). */
void solve_via_Cholesky(const matrix2D<double> &A, const matrix1D<double> &b,
   matrix1D<double> &result);

/** Evaluate quadratic form.
    Given x, c and H this function returns the value of the quadratic
    form val=c^t*x+0.5*x^t*H^t*H*x and the gradient of the quadratic form at x
    grad=c+H*x.
    
    Exceptions are thrown if the vectors and matrices do not have consistent
    dimensions.*/
void eval_quadratic(const matrix1D<double> &x, const matrix1D<double> &c,
   const matrix2D<double> &H, double &val, matrix1D<double> &grad);

/** Solves Quadratic programming subproblem. 
    \begin{verbatim}
       min 0.5*x'Cx + d'x   subject to:  A*x <= b 
        x                                Aeq*x=beq
      	                                 bl<=x<=bu
    \end{verbatim}
    */
void quadprog(const matrix2D<double> &C, const matrix1D<double> &d,
   const matrix2D<double> &A,   const matrix1D<double> &b,
   const matrix2D<double> &Aeq, const matrix1D<double> &beq,
         matrix1D<double> &bl,        matrix1D<double> &bu,
         matrix1D<double> &x);
	 
	 
/** Solves the least square problem.
    \begin{verbatim}
      min 0.5*(Norm(C*x-d))   subject to:  A*x <= b 
      x                                    Aeq*x=beq  
      	             	      	           bl<=x<=bu 
    \end{verbatim}
    */
void lsqlin(const matrix2D<double> &C, const matrix1D<double> &d,
   const matrix2D<double> &A,   const matrix1D<double> &b,
   const matrix2D<double> &Aeq, const matrix1D<double> &beq,
         matrix1D<double> &bl,        matrix1D<double> &bu,
         matrix1D<double> &x);
//@}
//@}
//@}

/* Implementations of common routines -------------------------------------- */
/* Print shape ------------------------------------------------------------- */
template <class T>
   void mT::print_shape(ostream &out) const {
   out << "Size(Y,X): " << YSIZE(*this) << "x" << XSIZE(*this)
       << " i=[" << STARTINGY(*this) << ".." << FINISHINGY(*this) << "]"
       << " j=[" << STARTINGX(*this) << ".." << FINISHINGX(*this) << "]";
}

/* Get size--- ------------------------------------------------------------- */
template <class T>
   void mT::get_size(int *size) const
   {size[0]=xdim; size[1]=ydim; size[2]=1;}

/* Outside ----------------------------------------------------------------- */
template <class T>
   bool mT::outside(const matrix1D<double> &v) const {
   if (XSIZE(v)<2)
      REPORT_ERROR(1,"Outside: index vector has got not enough components");
   return (XX(v)<STARTINGX(*this) || XX(v)>FINISHINGX(*this) ||
           YY(v)<STARTINGY(*this) || YY(v)>FINISHINGY(*this));
}

template <class T>
   bool mT::outside(int i, int j) const {
   return (j<STARTINGX(*this) || j>FINISHINGX(*this) ||
           i<STARTINGY(*this) || i>FINISHINGY(*this));
}

/* Intersects -------------------------------------------------------------- */
template <class T>
   bool mT::intersects(const mT &m) const
      {return intersects(STARTINGX(m), STARTINGY(m), XSIZE(m)-1, YSIZE(m)-1);}

template <class T>
   bool mT::intersects(const matrix1D<double> &corner1,
      const matrix1D<double> &corner2) const {
       if (XSIZE(corner1)!=2 || XSIZE(corner2)!=2)
          REPORT_ERROR(1002,"intersects 1D: corner sizes are not 1");
       return intersects(XX(corner1),YY(corner1),
          XX(corner2)-XX(corner1),YY(corner2)-YY(corner1));
}

template <class T>
   bool mT::intersects(double x0, double y0, double xdim, double ydim) const {
   SPEED_UP_temps;
   spduptmp0=MAX(STARTINGY(*this), y0);
   spduptmp1=MIN(FINISHINGY(*this),y0+ydim);
   if (spduptmp0>spduptmp1) return false;

   spduptmp0=MAX(STARTINGX(*this), x0);
   spduptmp1=MIN(FINISHINGX(*this),x0+xdim);
   if (spduptmp0>spduptmp1) return false;
   return true;
}

/* Corner ------------------------------------------------------------------ */
template <class T>
   bool mT::isCorner(const matrix1D<double> &v) {
   if (XSIZE(v)<2)
      REPORT_ERROR(1,"isCorner: index vector has got not enough components");
   return ((XX(v)==STARTINGX(*this)  && YY(v)==STARTINGY(*this))  ||
           (XX(v)==STARTINGX(*this)  && YY(v)==FINISHINGY(*this)) ||
           (XX(v)==FINISHINGX(*this) && YY(v)==STARTINGY(*this))  ||
           (XX(v)==FINISHINGX(*this) && YY(v)==FINISHINGY(*this)));
}

/* Border ------------------------------------------------------------------ */
template <class T>
   bool mT::isBorder(const matrix1D<int> &v) 
{
   if (XSIZE(v)<2)
      REPORT_ERROR(1,"isBorder: index vector has got not enough components");
   return  isBorder(YY(v),XX(v));
}

template <class T>
   bool mT::isBorder(int i,int j) 
{
   return (j==STARTINGX(*this)  || j==FINISHINGX(*this)  ||
           i==STARTINGY(*this)  || i==FINISHINGY(*this));
}

/* Patch ------------------------------------------------------------------- */
template <class T>
   void mT::patch(const mT &patch_array, char operation) {
      SPEED_UP_temps;
      FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(patch_array,*this)
         switch (operation) {
            case '=': MAT_ELEM(*this,i,j) =MAT_ELEM(patch_array,i,j); break;
            case '+': MAT_ELEM(*this,i,j)+=MAT_ELEM(patch_array,i,j); break;
            case '-': MAT_ELEM(*this,i,j)-=MAT_ELEM(patch_array,i,j); break;
            case '*': MAT_ELEM(*this,i,j)*=MAT_ELEM(patch_array,i,j); break;
            case '/': MAT_ELEM(*this,i,j)/=MAT_ELEM(patch_array,i,j); break;
         }
}

/* Statistics in region ---------------------------------------------------- */
template <class T>
void mT::compute_stats(double &avg, double &stddev, T &min_val, T &max_val,
   const matrix1D<double> &corner1, const matrix1D<double> &corner2) const {
   min_val=max_val=(*this)(corner1);
   matrix1D<double> r(3);
   double N=0, sum=0, sum2=0;
   FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1,corner2) {
      sum+=(*this)(r); sum2+=(*this)(r)*(*this)(r); N++;
      if      ((*this)(r)<min_val) min_val=(*this)(r);
      else if ((*this)(r)>max_val) max_val=(*this)(r);
   }
   if (N!=0) {
      avg=sum/N;
      stddev=sqrt(sum2/N-avg*avg);
   } else {avg=stddev=0;}
}

/* Minimum and maximum in region ------------------------------------------- */
template <class T>
void mT::compute_double_minmax(double &min_val, double &max_val,
   const matrix1D<double> &corner1, const matrix1D<double> &corner2) const {
   min_val=max_val=(*this)(corner1);
   matrix1D<double> r(2);
   FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1,corner2) {
      if      ((*this)(r)<min_val) min_val=(*this)(r);
      else if ((*this)(r)>max_val) max_val=(*this)(r);
   }
}

/* Center of mass ---------------------------------------------------------- */
template <class T>
   void mT::center_of_mass(matrix1D<double> &center, void * mask) {
      center.init_zeros(2);
      double mass=0;
      matrix2D<int> *imask=(matrix2D<int> *) mask;
      FOR_ALL_ELEMENTS_IN_MATRIX2D(*this) {
         if (imask==NULL || MAT_ELEM(*imask,i,j)) {
            XX(center)+=j*MAT_ELEM(*this,i,j);
            YY(center)+=i*MAT_ELEM(*this,i,j);
	    mass+=MAT_ELEM(*this,i,j);
      	 }
      }
      if (mass!=0) center/=mass;
   }

/* Output stream ----------------------------------------------------------- */
template <class T>
   ostream& operator << (ostream& ostrm, const mT& v) {
   if (XSIZE(v)==0 || YSIZE(v)==0)
      ostrm << "NULL matrix\n";
   else {
      ostrm << endl;
      double max_val=ABS(MULTIDIM_ELEM(v,0));
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(v)
         max_val=MAX(max_val,ABS(MULTIDIM_ELEM(v,i)));
      int prec=best_prec(max_val,10);

      for (int i=STARTINGY(v); i<=FINISHINGY(v); i++) {
         for (int j=STARTINGX(v); j<=FINISHINGX(v); j++) {
            ostrm << FtoA((double)MAT_ELEM(v,i,j),10,prec) << ' ';
         }
         ostrm << endl;
      }
   }

   return ostrm;
}

/* Solve Ax=b -------------------------------------------------------------- */
// (see Numerical Recipes, Chapter 2 Section 3)
template <class T>
   void solve(const mT &A, const vT &b, vT &result) {
   if (A.xdim==0)
      REPORT_ERROR(1108, "Solve: Matrix is empty");
   if (A.xdim!=A.ydim)
      REPORT_ERROR(1109, "Solve: Matrix is not squared");
   if (A.xdim!=b.get_dim())
      REPORT_ERROR(1102, "Solve: Different sizes of Matrix and Vector");
   if (b.isRow())
      REPORT_ERROR(1107, "Solve: Not correct vector shape");

   // Perform LU decomposition and then solve
   matrix1D<int> indx;
   T d;
   mT LU;
   ludcmp(A, LU, indx, d);
   result=b;
   lubksb(LU, indx, result);
}

template <class T>
   void solve_by_svd(const matrix2D<T> &A, const matrix1D<T> &b,
                     matrix1D<double> &result,double tolerance) {
   if (A.xdim==0)
      REPORT_ERROR(1108, "Solve: Matrix is empty");
   if (A.xdim!=A.ydim)
      REPORT_ERROR(1109, "Solve: Matrix is not squared");
   if (A.xdim!=b.get_dim())
      REPORT_ERROR(1102, "Solve: Different sizes of Matrix and Vector");
   if (b.isRow())
      REPORT_ERROR(1107, "Solve: Not correct vector shape");

   // First perform de single value decomposition
   // Xmipp interface that calls to svdcmp of numerical recipes
   matrix2D<double> u,v;
   matrix1D<double> w;
   svdcmp(A,u,w,v);
   
   // Here is checked if eigenvalues of the svd decomposition are acceptable
   // If a value is lower than tolerance, the it's zeroed, as this increases
   // the precision of the routine. 
   FOR_ALL_ELEMENTS_IN_MATRIX1D(w)
      if(w(i)<tolerance) w(i)=0;
	  
   // Set size of matrices
   result.resize(b.get_dim());

   // Xmipp interface that calls to svdksb of numerical recipes
   matrix1D<double> bd;
   type_cast(b,bd);
   svbksb(u,w,v,bd,result);  
}

// (see Numerical Recipes, Chapter 2 Section 1)
template <class T>
   void solve(const mT &A, const mT &b, mT &result) {
   if (A.xdim==0)
      REPORT_ERROR(1108, "Solve: Matrix is empty");
   if (A.xdim!=A.ydim)
      REPORT_ERROR(1109, "Solve: Matrix is not squared");
   if (A.ydim!=b.ydim)
      REPORT_ERROR(1102, "Solve: Different sizes of A and b");

   // Solve
   result=b;
   mT Aux=A;
   gaussj(Aux.adapt_for_numerical_recipes2(),Aux.ydim,
          result.adapt_for_numerical_recipes2(),b.xdim);
}

/* Interface to numerical recipes: ludcmp ---------------------------------- */
template <class T>
   void ludcmp(const mT &A, mT &LU, matrix1D<int> &indx, T &d) {
   LU=A;
   indx.resize(XSIZE(A));
   ludcmp (LU.adapt_for_numerical_recipes2(),XSIZE(A),
      indx.adapt_for_numerical_recipes(),&d);
}

/* Interface to Numerical Recipes: lubksb  --------------------------------- */
template <class T>
   void lubksb(const mT &LU, matrix1D<int> &indx, vT &b) {
   lubksb(LU.adapt_for_numerical_recipes2(),XSIZE(indx),
      indx.adapt_for_numerical_recipes(),
      b.adapt_for_numerical_recipes());
}

/* Interface to Numerical Recipes: svdcmp  --------------------------------- */
#define VIA_BILIB
template <class T>
void svdcmp(const matrix2D<T> &a,matrix2D<double> &u,
               matrix1D<double> &w, matrix2D<double> &v) {
   // svdcmp only works with double
   type_cast(a,u);
   
   // Set size of matrices
   w.init_zeros(u.ColNo());
   v.init_zeros(u.ColNo(),u.ColNo());
   
   // Call to the numerical recipes routine
#ifdef VIA_NR
   svdcmp(MULTIDIM_ARRAY(u),
          u.RowNo(),u.ColNo(),
          MULTIDIM_ARRAY(w),
          MULTIDIM_ARRAY(v));
#endif
#ifdef VIA_BILIB
   int status;
   SingularValueDecomposition (MULTIDIM_ARRAY(u),
          u.RowNo(),u.ColNo(),
          MULTIDIM_ARRAY(w),
          MULTIDIM_ARRAY(v),
	  5000, &status);
#endif
}
#undef VIA_NR
#undef VIA_BILIB

/* Apply a geometrical transformation -------------------------------------- */
// It is translated from function transforma in Lib/varios2.c
// We should check which one performs better.
//#define DEBUG_APPLYGEO
template <class T>
   void apply_geom(mT &M2,matrix2D<double> A, const mT &M1, bool inv,
   bool wrap, T outside) {
   int m1, n1, m2, n2;
   double x, y, xp, yp;
   double minxp, minyp, maxxp, maxyp;
   int   cen_x, cen_y, cen_xp, cen_yp;
   double wx, wy; // Weights in X,Y directions for bilinear interpolation
   int Xdim, Ydim;
   
   if ((XSIZE(A)!=3) || (YSIZE(A)!=3))
      REPORT_ERROR(1102,"Apply_geom: geometrical transformation is not 3x3");
   if (A.IsIdent())   {M2=M1; return;}
   if (XSIZE(M1)==0)  {M2.clear(); return;}
   
   if (!inv) A=A.inv();

   // For scalings the output matrix is resized outside to the final
   // size instead of being resized inside the routine with the
   // same size as the input matrix
   if (XSIZE(M2)==0) M2.resize(M1);
   if (outside!=0.) {
   // Initialise output matrix with value=outside
     FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D (M2) { DIRECT_MAT_ELEM(M2,i,j)=outside;}
   }
  
   // Find center and limits of image
   cen_y  = (int)(YSIZE(M2)/2);
   cen_x  = (int)(XSIZE(M2)/2);
   cen_yp = (int)(YSIZE(M1)/2);
   cen_xp = (int)(XSIZE(M1)/2);
   minxp  = -cen_xp;
   minyp  = -cen_yp;
   maxxp  = XSIZE(M1)-cen_xp-1;
   maxyp  = YSIZE(M1)-cen_yp-1;
   Xdim   = XSIZE(M1);
   Ydim   = YSIZE(M1);
   // Now we go from the output image to the input image, ie, for any pixel
   // in the output image we calculate which are the corresponding ones in
   // the original image, make an interpolation with them and put this value
   // at the output pixel
   //#define DEBUG_APPLYGEO
   #ifdef DEBUG_APPLYGEO
      cout << "A\n"    << A     << endl
	   << "(cen_x ,cen_y )=(" << cen_x  << "," << cen_y  << ")\n"
	   << "(cen_xp,cen_yp)=(" << cen_xp << "," << cen_yp << ")\n"
	   << "(min_xp,min_yp)=(" << minxp  << "," << minyp  << ")\n"
	   << "(max_xp,max_yp)=(" << maxxp  << "," << maxyp  << ")\n"
      ;
   #endif
   //#undef DEBUG_APPLYGEO
   for (int i=0; i<YSIZE(M2); i++) {
      // Calculate position of the beginning of the row in the output image
      x= -cen_x;
      y=i-cen_y;
      // For OldXmipp origins with even XSIZE & YSIZE:
      //      x= -cen_x+0.5;
      //      y=i-cen_y+0.5;

      // Calculate this position in the input image according to the
      // geometrical transformation
      // they are related by
      // coords_output(=x,y) = A * coords_input (=xp,yp)
      xp=x*dMij(A,0,0) + y*dMij(A,0,1) + dMij(A,0,2);
      yp=x*dMij(A,1,0) + y*dMij(A,1,1) + dMij(A,1,2);
     
      for (int j=0; j<XSIZE(M2); j++) {
         bool interp;
         T tmp;
         
         #ifdef DEBUG_APPLYGEO
         cout << "Computing (" << i << "," << j << ")\n";
         cout << "   (y, x) =(" << y << "," << x << ")\n"
              << "   before wrapping (y',x')=(" << yp << "," << xp << ") " << endl;
         #endif

         // If the point is outside the image, apply a periodic extension
         // of the image, what exits by one side enters by the other
         interp=true;
         if (wrap) {
            if (xp<minxp-XMIPP_EQUAL_ACCURACY ||
	        xp>maxxp+XMIPP_EQUAL_ACCURACY)
		xp=realWRAP(xp, minxp-0.5, maxxp+0.5);
            if (yp<minyp-XMIPP_EQUAL_ACCURACY ||
	        yp>maxyp+XMIPP_EQUAL_ACCURACY) yp=realWRAP(yp, minyp-0.5, maxyp+0.5);
         } else {
            if (xp<minxp-XMIPP_EQUAL_ACCURACY ||
	        xp>maxxp+XMIPP_EQUAL_ACCURACY) interp=false;
            if (yp<minyp-XMIPP_EQUAL_ACCURACY ||
	        yp>maxyp+XMIPP_EQUAL_ACCURACY) interp=false;
         }

         #ifdef DEBUG_APPLYGEO
         cout << "   after wrapping (y',x')=(" << yp << "," << xp << ") " << endl;
         cout << "   Interp = " << interp << endl;
         x++;
         #endif

         if (interp) {            
            // Calculate the integer position in input image, be careful
            // that it is not the nearest but the one at the top left corner
            // of the interpolation square. Ie, (0.7,0.7) would give (0,0)
            // Calculate also weights for point m1+1,n1+1
            wx=xp+cen_xp; m1=(int) wx; wx=wx-m1; m2=m1+1;
            wy=yp+cen_yp; n1=(int) wy; wy=wy-n1; n2=n1+1;
            // m2 and n2 can be out by 1 so wrap must be check here
            if (wrap) {
	       if(m2>=Xdim) m2=0;
	       if(n2>=Ydim) n2=0;
            }
            #ifdef DEBUG_APPLYGEO
            cout << "   From (" << n1 << "," << m1 << ") and ("
                 << n2 << "," << m2 << ")\n";
            cout << "   wx= " << wx << " wy= " << wy << endl;
            #endif

            // Perform interpolation
            // if wx == 0 means that the rightest point is useless for this
            // interpolation, and even it might not be defined if m1=xdim-1
            // The same can be said for wy.
               tmp  = (T) ((1-wy)*(1-wx)*dMij(M1,n1,m1));
            if (wx!=0 && m2<M1.xdim)
               tmp += (T) ((1-wy)*   wx *dMij(M1,n1,m2));
            if (wy!=0 && n2<M1.ydim){
               tmp += (T) (   wy *(1-wx)*dMij(M1,n2,m1));
               if (wx!=0 && m2<M1.xdim)
                  tmp += (T) (wy *   wx *dMij(M1,n2,m2));
            }
            dMij(M2,i,j) = tmp;

            #ifdef DEBUG_APPYGEO
            cout << "   val= " << tmp << endl;
            #endif
         }

         // Compute new point inside input image
         xp += dMij(A,0,0);
         yp += dMij(A,1,0);
      }
   }
}
#undef DEBUG_APPLYGEO

/* Apply a geometrical transformation -------------------------------------- */
// It is translated from function transforma in Lib/varios2.c
// We should check which one performs better.
//#define DEBUG
template <class T>
   void apply_geom_Bspline(mT &M2,matrix2D<double> A, const mT &M1,
   int Splinedegree, bool inv, bool wrap, T outside) {
   int m1, n1, m2, n2;
   double x, y, xp, yp;
   double minxp, minyp, maxxp, maxyp;
   int   cen_x, cen_y, cen_xp, cen_yp;
   
   if ((XSIZE(A)!=3) || (YSIZE(A)!=3))
      REPORT_ERROR(1102,"Apply_geom: geometrical transformation is not 3x3");
   if (A.IsIdent())   {M2=M1; return;}
   if (XSIZE(M1)==0)  {M2.clear(); return;}
   
   if (!inv) A=A.inv();

   // For scalings the output matrix is resized outside to the final
   // size instead of being resized inside the routine with the
   // same size as the input matrix
   if (XSIZE(M2)==0) M2.resize(M1);
   // Initialise output matrix with value=outside
   FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D (M2) { DIRECT_MAT_ELEM(M2,i,j)=outside;}
  
   // Find center and limits of image
   cen_y  = (int)(YSIZE(M2)/2);
   cen_x  = (int)(XSIZE(M2)/2);
   cen_yp = (int)(YSIZE(M1)/2);
   cen_xp = (int)(XSIZE(M1)/2);
   minxp  = -cen_xp;
   minyp  = -cen_yp;
   maxxp  = XSIZE(M1)-cen_xp-1;
   maxyp  = YSIZE(M1)-cen_yp-1;
   
   // Build the B-spline coefficients
   matrix2D<double> Bcoeffs;
   M1.produce_spline_coeffs(Bcoeffs,Splinedegree);
   STARTINGX(Bcoeffs)=(int)minxp;
   STARTINGY(Bcoeffs)=(int)minyp;

   // Now we go from the output image to the input image, ie, for any pixel
   // in the output image we calculate which are the corresponding ones in
   // the original image, make an interpolation with them and put this value
   // at the output pixel
   for (int i=0; i<YSIZE(M2); i++) {
      // Calculate position of the beginning of the row in the output image
      x= -cen_x;
      y=i-cen_y;

      // Calculate this position in the input image according to the
      // geometrical transformation
      // they are related by
      // coords_output(=x,y) = A * coords_input (=xp,yp)
      xp=x*dMij(A,0,0) + y*dMij(A,0,1) + dMij(A,0,2);
      yp=x*dMij(A,1,0) + y*dMij(A,1,1) + dMij(A,1,2);
     
      for (int j=0; j<XSIZE(M2); j++) {
         bool interp;
         T tmp;
         
         #ifdef DEBUG_APPLYGEO
         cout << "Computing (" << i << "," << j << ")\n";
         cout << "   (y, x) =(" << y << "," << x << ")\n"
              << "   before wrapping (y',x')=(" << yp << "," << xp << ") " << endl;
         #endif

         // If the point is outside the image, apply a periodic extension
         // of the image, what exits by one side enters by the other
         interp=true;
         if (wrap) {
            if (xp<minxp-XMIPP_EQUAL_ACCURACY ||
	        xp>maxxp+XMIPP_EQUAL_ACCURACY)
		xp=realWRAP(xp, minxp-0.5, maxxp+0.5);
            if (yp<minyp-XMIPP_EQUAL_ACCURACY ||
	        yp>maxyp+XMIPP_EQUAL_ACCURACY) yp=realWRAP(yp, minyp-0.5, maxyp+0.5);
         } else {
            if (xp<minxp-XMIPP_EQUAL_ACCURACY ||
	        xp>maxxp+XMIPP_EQUAL_ACCURACY) interp=false;
            if (yp<minyp-XMIPP_EQUAL_ACCURACY ||
	        yp>maxyp+XMIPP_EQUAL_ACCURACY) interp=false;
         }

         #ifdef DEBUG_APPLYGEO
         cout << "   after wrapping (y',x')=(" << yp << "," << xp << ") " << endl;
         cout << "   Interp = " << interp << endl;
         x++;
         #endif

         if (interp) {            
            dMij(M2,i,j) = (T)Bcoeffs.interpolated_elem_as_Bspline(
	       xp,yp,Splinedegree);

            #ifdef DEBUG_APPLYGEO
            cout << "   val= " << tmp << endl;
            #endif
         }

         // Compute new point inside input image
         xp += dMij(A,0,0);
         yp += dMij(A,1,0);
      }
   }
}
#undef DEBUG_APPLYGEO

template <>
inline void matrix2D< complex<double> >::scale_to_size_Bspline(int Splinedegree,
   int Ydim,int Xdim, matrix2D< complex<double> > &result) const {
   matrix2D<double> re, im, scre, scim;
   re.resize(YSIZE(*this),XSIZE(*this));
   im.resize(YSIZE(*this),XSIZE(*this));
   Complex2RealImag(MULTIDIM_ARRAY(*this),
      MULTIDIM_ARRAY(re),MULTIDIM_ARRAY(im),MULTIDIM_SIZE(*this));
   re.scale_to_size_Bspline(Splinedegree, Ydim, Xdim, scre);
   im.scale_to_size_Bspline(Splinedegree, Ydim, Xdim, scim);
   result.resize(Ydim,Xdim);
   RealImag2Complex(MULTIDIM_ARRAY(re),MULTIDIM_ARRAY(im),
      MULTIDIM_ARRAY(result),MULTIDIM_SIZE(re));
}

template <> bool matrix2D< complex<double> >::IsDiagonal() const;
template <> bool matrix2D< complex<double> >::IsScalar() const;
template <> bool matrix2D< complex<double> >::IsSymmetric() const;
template <> bool matrix2D< complex<double> >::IsSkewSymmetric() const;
template <> bool matrix2D< complex<double> >::IsUpperTriangular() const;
template <> bool matrix2D< complex<double> >::IsLowerTriangular() const;

/* Matrix by matrix multiplication ----------------------------------------- */
// xinit and yinit are taken from op1
template <class T>
   void mul_matrix(const mT &op1, const mT &op2, mT &result) {
   if (XSIZE(op1)!=YSIZE(op2))
      REPORT_ERROR(1102,"Not compatible sizes in matrix multiplication");
   
   result.init_zeros(YSIZE(op1),XSIZE(op2));
   for (int i=0; i<YSIZE(op1); i++)
      for (int j=0; j<XSIZE(op2); j++)
         for (int k=0; k<XSIZE(op1); k++)
            DIRECT_MAT_ELEM(result,i,j) += DIRECT_MAT_ELEM(op1,i,k)*
               DIRECT_MAT_ELEM(op2,k,j);
   STARTINGY(result)=STARTINGY(op1);
   STARTINGX(result)=STARTINGX(op1);
}

/* Vector by matrix multiplication ----------------------------------------- */
template <class T>
   vT vT::operator * (const mT &M) {
   vT result;
   if (XSIZE(*this)!=YSIZE(M))
      REPORT_ERROR(1102,"Not compatible sizes in matrix by vector");
   if (!isRow())
      REPORT_ERROR(1102,"Vector is not a row");
   
   result.init_zeros(XSIZE(M));
   for (int j=0; j<XSIZE(M); j++)
      for (int i=0; i<YSIZE(M); i++)
            DIRECT_VEC_ELEM(result,j) += DIRECT_VEC_ELEM(*this,i) *
               DIRECT_MAT_ELEM(M,i,j);
   result.setRow();
   STARTINGX(result)=STARTINGX(M);
   return result;
}

#undef maT
#undef maT1
#endif
