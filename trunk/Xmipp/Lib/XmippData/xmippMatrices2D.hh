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
#include <iostream.h>
#include <string>
#include <math.h>
#include <complex.h>

#include "xmippMatrices1D.hh"
#include "xmippFuncs.hh"

#define maT  matrix2D<T>
#define maT1 matrix2D<T1>
#undef  maTC
#define maTC matrix2D<double_complex>

/* ************************************************************************* */
/* FORWARD DEFINITIONS                                                       */
/* ************************************************************************* */
#include "Src/MultidimFriends.inc"
template <>
ostream& operator << (ostream & ostrm, const matrix2D<double_complex> &m);

template <class T>
   mT mul_elements(const mT &op1, const mT &op2);
template <class T>
   void mul_matrix(const mT &op1, const mT &op2, mT &result) _THROW;
template <class T>
   void solve (const mT &A, const vT &b, vT &result) _THROW;
template <class T>
   void solve (const mT &A, const mT &b, mT &result) _THROW;
// Use this solve_by_svd function to solve linear systems by singular value decomposition
// This is indicated for ill-conditioned systems and hard-to-solve problems
// (But it's slower)
template <class T>
   void solve_by_svd(const mT &A, const vT &b, matrix1D<double> &result,double tolerance) _THROW;
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
      bool wrap, T outside=(T)0) _THROW;

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
   friend void mul_matrix<>(const mT &op1, const mT &op2, mT &result) _THROW;

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
   matrix2D(int Ydim, int Xdim)
      {core_init(); init_shape(); resize(Ydim, Xdim);  __spcdim=2;}

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
   void init_identity(int Ydim, int Xdim);
   
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
   void from_vector(const vT &op1);
   
   /** Makes a vector from a matrix.
       An exception is thrown if the matrix is not a single row or
       a single column. The origin of the vector is set according to
       the one of the matrix.
       \\Ex: matrix1D<double> v; m.to_vector(v);*/
   void to_vector(vT &op1) const _THROW;
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
   void resize(int Ydim, int Xdim) _THROW {
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
   void load_from_numerical_recipes(T **m, int Ydim, int Xdim);

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
   T&   operator () (int i, int j) const _THROW {
        if (i<yinit || i>=yinit+ydim)
           REPORT_ERROR(1103,"Matrix subscript (i) out of range");
        if (j<xinit || j>=xinit+xdim)
           REPORT_ERROR(1103,"Matrix subscript (j) out of range");
        return MAT_ELEM(*this,i,j);}
   
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

   /** Interpolates the value of the 2D matrix M at the point (x,y) */
   T   interpolated_elem(double x, double y, T outside_value=(T)0);

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
   void getRow(int i, vT &v) const _THROW {
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
   void getCol(int j, vT &v) const _THROW {
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
   void setRow(int i, const vT &v) _THROW {
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
   void setCol(int j, const vT &v) _THROW {
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
      char operation) _THROW {
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
   mT transpose() const;
   
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
   void self_reverseX();
   
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
   void self_reverseY();
   
   /** Determinant of a matrix.
       An exception is thrown if the matrix is not squared or it is empty.
       \\Ex: double det=m.det();*/
   T det() const _THROW;

   /** Inverse of a matrix.
       The matrix is inverted using a SVD decomposition.
       In fact the pseudoinverse is returned.
       \\Ex: matrix2D<double> m1_inv; m1.inv(m1_inv);*/
   void inv(mT &result) const;

   /** Inverse of a matrix.*/
   mT inv() const {mT result; inv(result); return result;}

   /** Solve equation system.
       The equation system is defined by Ax=b, it is solved for x.
       Exceptions are thrown if the equation system is not solvable.
       \\Ex: matrix1D<double> x=m.inv();*/
   friend void solve<>(const mT &A, const vT &b, vT &result) _THROW;

   /** Solve equation system.
       The same as before but now, b is a matrix and so, x is also a matrix.
       A must be a square matrix.*/ 
   friend void solve<>(const mT &A, const mT &b, mT &result) _THROW;

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
   void window(int y0, int x0, int yF, int xF, T init_value=0);
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
       const mT &m1, bool inv, bool wrap, T outside) _THROW;
   /** Self apply geom.
       Same as the previous one, but the result is kept in this object */
   void self_apply_geom(matrix2D<double> A, bool inv, bool wrap, T outside=(T)0)
      {mT aux; apply_geom(aux,A,*this,inv,wrap,outside); *this=aux;}

   #define IS_INV     TRUE
   #define IS_NOT_INV FALSE
   #define DONT_WRAP  FALSE
   #define WRAP       TRUE
   
   /** Rotate matrix.
       The angle must be in degrees. The result cannot be this object.
       \\Ex: m1.rotate(60,m2);*/
   void rotate(double ang, mT &result, bool wrap=DONT_WRAP) const;
   
   /** Rotate matrix.
       Same as the previous one. */
   mT rotate(double ang, bool wrap=DONT_WRAP) const
      {mT aux; rotate(ang,aux,wrap); return aux;}

   /** Rotate matrix.
       Same as the previous one, but the result is kept in this object */
   void self_rotate(double ang, bool wrap=DONT_WRAP)
      {mT aux; rotate(ang,aux,wrap); *this=aux;}

   /** Translate matrix.
       The displacement is given as a R2 vector of the form (shift_X,shift_Y).
       The result cannot be this object.
       \\Ex: m2=m1.translate(vector_R2(0,2));
       \\--> m1 is shifted 2 pixels down and stored in m2*/
   void translate(const matrix1D<double> &v, mT &result, bool wrap=WRAP) const;
   
   /** Translate matrix.
       Same as the previous one. */
   mT translate(const matrix1D<double> &v, bool wrap=WRAP) const
      {mT aux; translate(v,aux,wrap); return aux;}
   
   /** Translate matrix.
       Same as the previous one, but the result is kept in this object */
   void self_translate(const matrix1D<double> &v, bool wrap=WRAP)
      {mT aux; translate(v,aux,wrap); *this=aux;}

   /** Translate center of mass to center.
       If the input has very high values, it is better to rescale it to
       be between 0 and 1. */
   void self_translate_center_of_mass_to_center(bool wrap=WRAP);

   /** Scales to a new size.
       The matrix is scaled (resampled) to fill a new size. It is not the
       same as "window" in this same class. The size can be larger or 
       smaller than the actual one. But the result matrix cannot be
       this object.
       \\Ex: m1.scale_to_size(128,128);*/
   void scale_to_size(int Ydim,int Xdim, mT &result) const;

   /** Scales to a new size.
       Same as the previous one. */
   mT scale_to_size(int Ydim,int Xdim) const
      {mT aux; scale_to_size(Ydim, Xdim, aux); return aux;}

   /** Scales to a new size.
       Same as the previous one, but the result is kept in this object. */
   void self_scale_to_size(int Ydim,int Xdim)
      {mT aux; scale_to_size(Ydim, Xdim, aux); *this=aux;}

   /** Superpixel reduce.
       This function reduces the given image averaging in the superpixel
       of the given size. The last columns and rows are removed if the
       size of the original image is not an exact multiple of the given
       superpixel size */
   void superpixel_reduce(mT &result, int size=2) const;

   /** Superpixel expand.
       This function copies each pixel to a new image as many times
       as the size of the superpixel. */
   void superpixel_expand(mT &result, int size=2) const;
   //@}

   /* Iterators ------------------------------------------------------------ */
   /**@name Iterators*/
   //@{
   /** Apply the same scalar function to all rows.
       This function must take a row vector and return a single value, a column
       vector with these values is returned.
       \\Ex:T vector_sum(vT &v) {return v.sum();};
            v1=m.for_all_rows(&vector_sum);*/
   vT for_all_rows (T (*f)(vT&)) const;

   /** Apply the same scalar function to all columns.
       This function must take a column vector and return a single value, a row
       vector with these values is returned.
       \\Ex:T vector_sum(vT &v) {return v.sum();};
            v1=m.for_all_cols(&vector_sum);*/
   vT for_all_cols (T (*f)(vT&)) const;
   
   /** Apply the same vectorial function to all rows.
       This function must take a row vector and return a row vector (of the
       same size as the input one), a new
       matrix with these transformed rows is returned.
       \\Ex:vT vector_norm(vT &v) {return v.normalize();};
            m2=m.for_all_rows(&vector_norm);*/
   mT for_all_rows (vT (*f)(vT&)) const
      {mT aux(*this); aux.for_all_rows(f); return aux;}

   /** Apply a vectorial function to all rows, keep in this object. */
   void for_all_rows (vT (*f)(vT&));

   /** Apply the same vectorial function to all columns.
       This function must take a columnm vector and return a column vector
       (of the same size as the input one), a
       new matrix with these transformed columns is returned.
       \\Ex:vT vector_norm(vT &v) {return v.normalize();};
            m2=m.for_all_cols(&vector_norm);*/
   mT for_all_cols (vT (*f)(vT&)) const
      {mT aux(*this); aux.for_all_cols(f); return aux;}

   /** Apply a vectorial function to all rows, keep in this object. */
   void for_all_cols (vT (*f)(vT&));

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
   vT operator * (const vT &op1) const;
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
   bool IsDiagonal() const;

   /** True if the matrix is scalar.
       A scalar matrix is diagonal and all the values at the diagonal are
       the same*/
   bool IsScalar() const;

   /** True if the matrix is symmetric.
       \\Ex: if (m.IsSymmetric()) cout << "The matrix is symmetric\n";*/
   bool IsSymmetric() const;

   /** True if the matrix is skew-symmetric (anti-symmetric).
       \\Ex: if (m.IsSkewSymmetric()) cout << "The matrix is skewsymmetric\n";*/
   bool IsSkewSymmetric() const;

   /** True if the matrix is upper-triangular.
       \\Ex: if (m.IsUpperTriangular()) cout << "The matrix is upper triangular\n";*/
   bool IsUpperTriangular() const;

   /** True if the matrix is lower-triangular.
       \\Ex: if (m.IsLowerTriangular()) cout << "The matrix is lower triangular\n";*/
   bool IsLowerTriangular() const;

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
   void max_index(int &i, int &j) const;

   /** Minimum element.
       This function returns the index of the minimum element of an array.
       array(i,j). Returns -1 if the array is empty*/
   void min_index(int &i, int &j) const;
   //@}
};

/**@name Related functions
   These functions are not methods of matrix1D */
//@{

/* Geometry ---------------------------------------------------------------- */
/**@name Geometry
   These functions allow you to create easily R2 and R3 vectors. You must use
   all these functions with IS_NOT_INV in \Ref{apply_geom}*/
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
matrix2D<double> translation2D_matrix(const matrix1D<double> v) _THROW;

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
matrix2D<double> rot3D_matrix(double ang, char axis) _THROW;

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
matrix2D<double> align_with_Z(const matrix1D<double> &axis) _THROW;

/** Creates a translational matrix (4x4) for volumes.
    The shift is given as a R3 vector (shift_X, shift_Y, shift_Z);
    An exception is thrown if the displacement is not a R3 vector.
    \\Ex: m=translation3D_matrix(vector_R3(0,0,2));
    \\--> Displacement of 2 pixels down */
matrix2D<double> translation3D_matrix(const matrix1D<double> &v) _THROW;

/** Creates a scaling matrix (4x4) for volumes.
    The scaling factors for the different axis must be given as a vector.
    So that, XX(sc)=scale for X axis, YY(sc)=...*/
matrix2D<double> scale3D_matrix(const matrix1D<double> &sc) _THROW;

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
   void cut_to_common_size(mT &m1, mT &m2);

/** Does a radial average of a matrix, around the pixel where is the origin.
    A vector is returned where:
	 - the first element is the mean of the pixels whose
	   distance to the origin is (0-1),
	 - the second element is the mean of the pixels 
       whose distance to the origin is (1-2)
	 - and so on. */
template <class T>
void radial_average(const matrix2D<T> &m, const matrix1D<int> &center_of_rot,
                    matrix1D<T> &radial_mean) _THROW;

/** Solve equation system, nonnegative solution.
    The equation system is defined by Ax=b, it is solved for x.
    x is forced to be nonnegative. It is designed to cope with
    large equation systems. This function is borrowed from
    LAPACK nnls.
    
    The norm of the vector Ax-b is returned.
    \\Ex: matrix1D<double> x=m.inv();*/
double solve_nonneg(const matrix2D<double> &A, const matrix1D<double> &b,
   matrix1D<double> &result) _THROW;

/** Evaluate quadratic form.
    Given x, c and H this function returns the value of the quadratic
    form val=c^t*x+0.5*x^t*H^t*H*x and the gradient of the quadratic form at x
    grad=c+H*x.
    
    Exceptions are thrown if the vectors and matrices do not have consistent
    dimensions.*/
void eval_quadratic(const matrix1D<double> &x, const matrix1D<double> &c,
   const matrix2D<double> &H, double &val, matrix1D<double> &grad) _THROW;
//@}
//@}
//@}
#undef maT
#undef maT1
#endif
