/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef MATRIX2D_H_
#define MATRIX2D_H_

#include <fstream>
#include <errno.h>
#include <fcntl.h>
#include "../../external/bilib/headers/linearalgebra.h"
#include "xmipp_macros.h"
#include "xmipp_filename.h"
#include "xmipp_error.h"
#include "matrix1d.h"

#ifdef XMIPP_MMAP
#include <sys/mman.h>
#endif


// Forward declarations
template<typename T>
class Matrix2D;

template<typename T>
void ludcmp(const Matrix2D<T>& A, Matrix2D<T>& LU, Matrix1D< int >& indx, T& d);

template<typename T>
void lubksb(const Matrix2D<T>& LU, Matrix1D< int >& indx, Matrix1D<T>& b);

template<typename T>
void svdcmp(const Matrix2D< T >& a,
            Matrix2D< double >& u,
            Matrix1D< double >& w,
            Matrix2D< double >& v);

void svbksb(Matrix2D< double >& u,
            Matrix1D< double >& w,
            Matrix2D< double >& v,
            Matrix1D< double >& b,
            Matrix1D< double >& x);


/** Cholesky decomposition.
 * Given M, this function decomposes M as M=L*L^t where L is a lower triangular matrix.
 * M must be positive semi-definite.
 */
void cholesky(const Matrix2D<double> &M, Matrix2D<double> &L);

/** Schur decomposition.
 * Given M, this function decomposes M as M = O*T*O' where O is an orthogonal matrix.
 */
void schur(const Matrix2D<double> &M, Matrix2D<double> &O, Matrix2D<double> &T);

/** @defgroup Matrices Matrix2D Matrices
 * @ingroup DataLibrary
 */
//@{
/** @name Matrices speed up macros */
//@{

/** Array access.
 *
 * This macro gives you access to the array (T)
 */
#define MATRIX2D_ARRAY(m) ((m).mdata)

/** For all elements in the array
 *
 * This macro is used to generate loops for the matrix in an easy way. It
 * defines internal indexes 'i' and 'j' which ranges the matrix using its
 * mathematical definition (ie, logical access).
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_MATRIX2D(m)
 * {
 *     std::cout << m(i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_MATRIX2D(m) \
    for (size_t i=0; i<(m).mdimy; i++) \
        for (size_t j=0; j<(m).mdimx; j++)

/** Access to a matrix element
 * v is the array, i and j define the element v_ij.
 *
  * @code
 * MAT_ELEM(m, 0, 0) = 1;
 * val = MAT_ELEM(m, 0, 0);
 * @endcode
 */
#define MAT_ELEM(m,i,j) ((m).mdata[((size_t)i)*(m).mdimx+((size_t)j)])

/** X dimension of the matrix
 */
#define MAT_XSIZE(m) ((m).mdimx)

/** Y dimension of the matrix
 */
#define MAT_YSIZE(m) ((m).mdimy)

/** Total elements of the matrix
 */
#define MAT_SIZE(m) ((m).mdim)

/** Matrix element: Element access
 *
 * This is just a redefinition
 * of the function above
 * @code
 * dMij(m, -2, 1) = 1;
 * val = dMij(m, -2, 1);
 * @endcode
 */
#define dMij(m, i, j)  MAT_ELEM(m, i, j)

/** Matrix element: Element access
 *
 * This is just a redefinition
 * of the function above
 */
#define dMn(m, n)  ((m).mdata[(n)])

/** Matrix (3x3) by vector (3x1) (a=M*b)
 *
 * You must "load" the temporary variables, and create the result vector with
 * the appropiate size. You can reuse the vector b to store the results (that
 * is, M3x3_BY_V3x1(b, M, b);, is allowed).
 *
 * @code
 * double example
 * {
 *     SPEED_UP_temps;
 *
 *     Matrix1D< double > a(3), b(3);
 *     Matrix2D< double > M(3, 3);
 *
 *     M.init_random(0, 1);
 *     b.init_random(0, 1);
 *     M3x3_BY_V3x1(a, M, b);
 *
 *     return a.sum();
 * }
 * @endcode
 */
#define M3x3_BY_V3x1(a, M, b) { \
        spduptmp0 = dMn(M, 0) * XX(b) + dMn(M, 1) * YY(b) + dMn(M, 2) * ZZ(b); \
        spduptmp1 = dMn(M, 3) * XX(b) + dMn(M, 4) * YY(b) + dMn(M, 5) * ZZ(b); \
        spduptmp2 = dMn(M, 6) * XX(b) + dMn(M, 7) * YY(b) + dMn(M, 8) * ZZ(b); \
        XX(a) = spduptmp0; YY(a) = spduptmp1; ZZ(a) = spduptmp2; }

/** Matrix (3x3) by Matrix (3x3) (A=B*C)
 *
 * You must "load" the temporary variables, and create the result vector with
 * the appropiate size. You can reuse any of the multiplicands to store the
 * results (that is, M3x3_BY_M3x3(A, A, B);, is allowed).
 */
#define M3x3_BY_M3x3(A, B, C) { \
        spduptmp0 = dMn(B,0) * dMn(C,0) + dMn(B,1) * dMn(C,3) + dMn(B,2) * dMn(C,6); \
        spduptmp1 = dMn(B,0) * dMn(C,1) + dMn(B,1) * dMn(C,4) + dMn(B,2) * dMn(C,7); \
        spduptmp2 = dMn(B,0) * dMn(C,2) + dMn(B,1) * dMn(C,5) + dMn(B,2) * dMn(C,8); \
        spduptmp3 = dMn(B,3) * dMn(C,0) + dMn(B,4) * dMn(C,3) + dMn(B,5) * dMn(C,6); \
        spduptmp4 = dMn(B,3) * dMn(C,1) + dMn(B,4) * dMn(C,4) + dMn(B,5) * dMn(C,7); \
        spduptmp5 = dMn(B,3) * dMn(C,2) + dMn(B,4) * dMn(C,5) + dMn(B,5) * dMn(C,8); \
        spduptmp6 = dMn(B,6) * dMn(C,0) + dMn(B,7) * dMn(C,3) + dMn(B,8) * dMn(C,6); \
        spduptmp7 = dMn(B,6) * dMn(C,1) + dMn(B,7) * dMn(C,4) + dMn(B,8) * dMn(C,7); \
        spduptmp8 = dMn(B,6) * dMn(C,2) + dMn(B,7) * dMn(C,5) + dMn(B,8) * dMn(C,8); \
        dMn(A, 0) = spduptmp0; \
        dMn(A, 1) = spduptmp1; \
        dMn(A, 2) = spduptmp2; \
        dMn(A, 3) = spduptmp3; \
        dMn(A, 4) = spduptmp4; \
        dMn(A, 5) = spduptmp5; \
        dMn(A, 6) = spduptmp6; \
        dMn(A, 7) = spduptmp7; \
        dMn(A, 8) = spduptmp8; }

/** Matrix (2x2) by vector (2x1) (a=M*b)
 *
 * You must "load" the temporary variables, and create the result vector with
 * the appropiate size. You can reuse the vector b to store the results (that
 * is, M2x2_BY_V2x1(b, M, b);, is allowed).
 *
 * @code
 * double example
 * {
 *     SPEED_UP_temps;
 *
 *     Matrix1D< double > a(2), b(2);
 *     Matrix2D< double > M(2, 2);
 *
 *     M.init_random(0, 1);
 *     b.init_random(0, 1);
 *
 *     M2x2_BY_V2x1(a, M, b);
 *
 *     return a.sum();
 * }
 * @endcode
 */
#define M2x2_BY_V2x1(a, M, b) { \
        spduptmp0 = dMn(M, 0) * XX(b) + dMn(M, 1) * YY(b); \
        spduptmp1 = dMn(M, 2) * XX(b) + dMn(M, 3) * YY(b); \
        XX(a) = spduptmp0; \
        YY(a) = spduptmp1; }

/** Matrix (2x2) by constant (M2=M1*k)
 *
 * You must create the result matrix with the appropriate size. You can reuse
 * the matrix M1 to store the results (that is, M2x2_BY_CT(M, M, k);, is
 * allowed).
 */
#define M2x2_BY_CT(M2, M1, k) { \
        dMn(M2, 0) = dMn(M1, 0) * k; \
        dMn(M2, 1) = dMn(M1, 1) * k; \
        dMn(M2, 2) = dMn(M1, 2) * k; \
        dMn(M2, 3) = dMn(M1, 3) * k; }

/** Matrix (3x3) by constant (M2=M1*k)
 *
 * You must create the result matrix with the appropiate size. You can reuse the
 * matrix M1 to store the results (that is, M2x2_BY_CT(M, M, k);, is allowed).
 */
#define M3x3_BY_CT(M2, M1, k) { \
        dMn(M2, 0) = dMn(M1, 0) * k; \
        dMn(M2, 1) = dMn(M1, 1) * k; \
        dMn(M2, 2) = dMn(M1, 2) * k; \
        dMn(M2, 3) = dMn(M1, 3) * k; \
        dMn(M2, 4) = dMn(M1, 4) * k; \
        dMn(M2, 5) = dMn(M1, 5) * k; \
        dMn(M2, 6) = dMn(M1, 6) * k; \
        dMn(M2, 7) = dMn(M1, 7) * k; \
        dMn(M2, 8) = dMn(M1, 8) * k; }
/** Matrix (4x4) by constant (M2=M1*k)
 *
 * You must create the result matrix with the appropiate size. You can reuse the
 * matrix M1 to store the results (that is, M2x2_BY_CT(M, M, k);, is allowed).
 */
#define M4x4_BY_CT(M2, M1, k) { \
        dMn(M2, 0) = dMn(M1, 0) * k; \
        dMn(M2, 1) = dMn(M1, 1) * k; \
        dMn(M2, 2) = dMn(M1, 2) * k; \
        dMn(M2, 3) = dMn(M1, 3) * k; \
        dMn(M2, 4) = dMn(M1, 4) * k; \
        dMn(M2, 5) = dMn(M1, 5) * k; \
        dMn(M2, 6) = dMn(M1, 6) * k; \
        dMn(M2, 7) = dMn(M1, 7) * k; \
        dMn(M2, 8) = dMn(M1, 8) * k; \
  dMn(M2, 9) = dMn(M1, 9) * k; \
  dMn(M2,10) = dMn(M1,10) * k; \
  dMn(M2,11) = dMn(M1,11) * k; \
  dMn(M2,12) = dMn(M1,12) * k; \
  dMn(M2,13) = dMn(M1,13) * k; \
  dMn(M2,14) = dMn(M1,14) * k; \
  dMn(M2,15) = dMn(M1,15) * k;}

/** Inverse of a matrix (2x2)
 *
 * Input and output matrix cannot be the same one. The output is supposed to be
 * already resized.
 */
#define M2x2_INV(Ainv, A) { \
  spduptmp0 = dMn(A,0) * dMn(A,3) - dMn(A,1) * dMn(A,2);\
  if (spduptmp0==0.0) \
   REPORT_ERROR(ERR_NUMERICAL,"2x2 matrix is not invertible"); \
        spduptmp0 = 1.0 / spduptmp0; \
        dMn(Ainv, 0) =  dMn(A,3); \
        dMn(Ainv, 1) = -dMn(A,1); \
        dMn(Ainv, 2) = -dMn(A,2); \
        dMn(Ainv, 3) =  dMn(A,0); \
        M2x2_BY_CT(Ainv, Ainv, spduptmp0); }

/** Inverse of a matrix (3x3)
 *
 * Input and output matrix cannot be the same one. The output is supposed to be
 * already resized.
 */
#define M3x3_INV(Ainv, A) { \
        dMn(Ainv, 0) =   dMn(A,8)*dMn(A,4)-dMn(A,7)*dMn(A,5); \
        dMn(Ainv, 1) = -(dMn(A,8)*dMn(A,1)-dMn(A,7)*dMn(A,2)); \
        dMn(Ainv, 2) =   dMn(A,5)*dMn(A,1)-dMn(A,4)*dMn(A,2); \
        dMn(Ainv, 3) = -(dMn(A,8)*dMn(A,3)-dMn(A,6)*dMn(A,5)); \
        dMn(Ainv, 4) =   dMn(A,8)*dMn(A,0)-dMn(A,6)*dMn(A,2); \
        dMn(Ainv, 5) = -(dMn(A,5)*dMn(A,0)-dMn(A,3)*dMn(A,2)); \
        dMn(Ainv, 6) =   dMn(A,7)*dMn(A,3)-dMn(A,6)*dMn(A,4); \
        dMn(Ainv, 7) = -(dMn(A,7)*dMn(A,0)-dMn(A,6)*dMn(A,1)); \
        dMn(Ainv, 8) =   dMn(A,4)*dMn(A,0)-dMn(A,3)*dMn(A,1); \
        spduptmp0 = dMn(A,0)*dMn(Ainv,0)+dMn(A,3)*dMn(Ainv,1)+dMn(A,6)*dMn(Ainv,2); \
  if (spduptmp0==0.0) \
   REPORT_ERROR(ERR_NUMERICAL,"3x3 matrix is not invertible"); \
     spduptmp0 = 1.0 / spduptmp0; \
        M3x3_BY_CT(Ainv, Ainv, spduptmp0); }

/** Inverse of a matrix (4x4)
 *
 * Input and output matrix cannot be the same one. The output is supposed to be
 * already resized.
 */
#define M4x4_INV(Ainv, A) { \
        dMn(Ainv, 0) =   dMn(A,5)*(dMn(A,10)*dMn(A,15)-dMn(A,11)*dMn(A,14))+\
             dMn(A,6)*(dMn(A,11)*dMn(A,13)-dMn(A,9) *dMn(A,15))+\
             dMn(A,7)*(dMn(A,9) *dMn(A,14)-dMn(A,10)*dMn(A,13));\
  dMn(Ainv, 1) =   dMn(A,1)*(dMn(A,11)*dMn(A,14)-dMn(A,10)*dMn(A,15))+\
       dMn(A,2)*(dMn(A,9) *dMn(A,15)-dMn(A,11)*dMn(A,13))+\
       dMn(A,3)*(dMn(A,10)*dMn(A,13)-dMn(A,9) *dMn(A,14));\
  dMn(Ainv, 2) =   dMn(A,1)*(dMn(A,6) *dMn(A,15)-dMn(A,7) *dMn(A,14))+\
       dMn(A,2)*(dMn(A,7) *dMn(A,13)-dMn(A,5) *dMn(A,15))+\
       dMn(A,3)*(dMn(A,5) *dMn(A,14)-dMn(A,6) *dMn(A,13));\
  dMn(Ainv, 3) =   dMn(A,1)*(dMn(A,7) *dMn(A,10)-dMn(A,6) *dMn(A,11))+\
       dMn(A,2)*(dMn(A,5) *dMn(A,11)-dMn(A,7) *dMn(A,9))+\
       dMn(A,3)*(dMn(A,6) *dMn(A,9) -dMn(A,5) *dMn(A,10));\
  dMn(Ainv, 4) =   dMn(A,4)*(dMn(A,11)*dMn(A,14)-dMn(A,10)*dMn(A,15))+\
       dMn(A,6)*(dMn(A,8) *dMn(A,15)-dMn(A,11)*dMn(A,12))+\
       dMn(A,7)*(dMn(A,10)*dMn(A,12)-dMn(A,8) *dMn(A,14));\
  dMn(Ainv, 5) =   dMn(A,0)*(dMn(A,10)*dMn(A,15)-dMn(A,11)*dMn(A,14))+\
       dMn(A,2)*(dMn(A,11)*dMn(A,12)-dMn(A,8) *dMn(A,15))+\
       dMn(A,3)*(dMn(A,8) *dMn(A,14)-dMn(A,10)*dMn(A,12));\
  dMn(Ainv, 6) =   dMn(A,0)*(dMn(A,7) *dMn(A,14)-dMn(A,6) *dMn(A,15))+\
       dMn(A,2)*(dMn(A,4) *dMn(A,15)-dMn(A,7) *dMn(A,12))+\
       dMn(A,3)*(dMn(A,6) *dMn(A,12)-dMn(A,4) *dMn(A,14));\
  dMn(Ainv, 7) =   dMn(A,0)*(dMn(A,6) *dMn(A,11)-dMn(A,7) *dMn(A,10))+\
       dMn(A,2)*(dMn(A,7) *dMn(A,8) -dMn(A,4) *dMn(A,11))+\
       dMn(A,3)*(dMn(A,4) *dMn(A,10)-dMn(A,6) *dMn(A,8));\
  dMn(Ainv, 8) =   dMn(A,4)*(dMn(A,9) *dMn(A,15)-dMn(A,11)*dMn(A,13))+\
       dMn(A,5)*(dMn(A,11)*dMn(A,12)-dMn(A,8) *dMn(A,15))+\
       dMn(A,7)*(dMn(A,8) *dMn(A,13)-dMn(A,9) *dMn(A,12));\
  dMn(Ainv, 9) =   dMn(A,0)*(dMn(A,11)*dMn(A,13)-dMn(A,9) *dMn(A,15))+\
       dMn(A,1)*(dMn(A,8) *dMn(A,15)-dMn(A,11)*dMn(A,12))+\
       dMn(A,3)*(dMn(A,9) *dMn(A,12)-dMn(A,8) *dMn(A,13));\
  dMn(Ainv,10) =   dMn(A,0)*(dMn(A,5) *dMn(A,15)-dMn(A,7) *dMn(A,13))+\
       dMn(A,1)*(dMn(A,7) *dMn(A,12)-dMn(A,4) *dMn(A,15))+\
       dMn(A,3)*(dMn(A,4) *dMn(A,13)-dMn(A,5) *dMn(A,12));\
  dMn(Ainv,11) =   dMn(A,0)*(dMn(A,7) *dMn(A,9) -dMn(A,5) *dMn(A,11))+\
       dMn(A,1)*(dMn(A,4) *dMn(A,11)-dMn(A,7) *dMn(A,8))+\
       dMn(A,3)*(dMn(A,5) *dMn(A,8) -dMn(A,4) *dMn(A,9));\
  dMn(Ainv,12) =   dMn(A,4)*(dMn(A,10)*dMn(A,13)-dMn(A,9) *dMn(A,14))+\
       dMn(A,5)*(dMn(A,8) *dMn(A,14)-dMn(A,10)*dMn(A,12))+\
       dMn(A,6)*(dMn(A,9) *dMn(A,12)-dMn(A,8) *dMn(A,13));\
  dMn(Ainv,13) =   dMn(A,0)*(dMn(A,9) *dMn(A,14)-dMn(A,10)*dMn(A,13))+\
       dMn(A,1)*(dMn(A,10)*dMn(A,12)-dMn(A,8) *dMn(A,14))+\
       dMn(A,2)*(dMn(A,8) *dMn(A,13)-dMn(A,9) *dMn(A,12));\
  dMn(Ainv,14) =   dMn(A,0)*(dMn(A,6) *dMn(A,13)-dMn(A,5) *dMn(A,14))+\
       dMn(A,1)*(dMn(A,4) *dMn(A,14)-dMn(A,6) *dMn(A,12))+\
       dMn(A,2)*(dMn(A,5) *dMn(A,12)-dMn(A,4) *dMn(A,13));\
  dMn(Ainv,15) =   dMn(A,0)*(dMn(A,5) *dMn(A,10)-dMn(A,6) *dMn(A,9))+\
       dMn(A,1)*(dMn(A,6) *dMn(A,8) -dMn(A,4) *dMn(A,10))+\
       dMn(A,2)*(dMn(A,4) *dMn(A,9) -dMn(A,5) *dMn(A,8));\
        spduptmp0 = dMn(A,0)*(dMn(A,5)*(dMn(A,10)*dMn(A,15)-dMn(A,11)*dMn(A,14))\
                              +dMn(A,6)*(dMn(A,11)*dMn(A,13)-dMn(A,9) *dMn(A,15))\
                              +dMn(A,7)*(dMn(A,9) *dMn(A,14)-dMn(A,10)*dMn(A,13)))\
                    +dMn(A,1)*(dMn(A,4)*(dMn(A,11)*dMn(A,14)-dMn(A,10)*dMn(A,15))\
                              +dMn(A,6)*(dMn(A,8) *dMn(A,15)-dMn(A,11)*dMn(A,12))\
                              +dMn(A,7)*(dMn(A,10)*dMn(A,12)-dMn(A,8) *dMn(A,14)))\
        +dMn(A,2)*(dMn(A,4)*(dMn(A,9) *dMn(A,15)-dMn(A,11)*dMn(A,13))\
         +dMn(A,5)*(dMn(A,11)*dMn(A,12)-dMn(A,8) *dMn(A,15))\
         +dMn(A,7)*(dMn(A,8) *dMn(A,13)-dMn(A,9) *dMn(A,12)))\
        +dMn(A,3)*(dMn(A,4)*(dMn(A,10)*dMn(A,13)-dMn(A,9) *dMn(A,14))\
         +dMn(A,5)*(dMn(A,8) *dMn(A,14)-dMn(A,10)*dMn(A,12))\
         +dMn(A,6)*(dMn(A,9) *dMn(A,12)-dMn(A,8) *dMn(A,13))); \
  if (spduptmp0==0.0) \
   REPORT_ERROR(ERR_NUMERICAL,"4x4 matrix is not invertible"); \
  spduptmp0 = 1.0 / spduptmp0; \
        M4x4_BY_CT(Ainv, Ainv, spduptmp0); }

/** Matrix2D class */
template<typename T>
class Matrix2D
{
public:
    // The array itself
    T* mdata;

    // Destroy data
    bool destroyData;

    // Mapped data
    bool mappedData;

    // File descriptor for mapped files
    int fdMap;

    // Mapped data original pointer
    char* mdataOriginal;

    // Number of elements in X
    size_t mdimx;

    // Number of elements in Y
    size_t mdimy;

    // Total number of elements
    size_t mdim;
    //@}

    /// @name Constructors
    //@{
    /** Empty constructor
     */
    Matrix2D()
    {
        coreInit();
    }

    Matrix2D(const FileName &fnMappedMatrix, int Ydim, int Xdim, size_t offset=0)
    {
        coreInit(fnMappedMatrix,Ydim,Xdim,offset);
    }

    /** Dimension constructor
     */
    Matrix2D(int Ydim, int Xdim)
    {
        coreInit();
        initZeros(Ydim, Xdim);
    }

    /** Copy constructor
     */
    Matrix2D(const Matrix2D<T>& v)
    {
        coreInit();
        *this = v;
    }

    /** Destructor.
     */
    ~Matrix2D()
    {
        coreDeallocate();
    }

    /** Assignment.
     *
     * You can build as complex assignment expressions as you like. Multiple
     * assignment is allowed.
     *
     * @code
     * v1 = v2 + v3;
     * v1 = v2 = v3;
     * @endcode
     */
    Matrix2D<T>& operator=(const Matrix2D<T>& op1)
    {
        if (&op1 != this)
        {
            if (MAT_XSIZE(*this)!=MAT_XSIZE(op1) ||
                MAT_YSIZE(*this)!=MAT_YSIZE(op1))
                resizeNoCopy(op1);
            memcpy(mdata,op1.mdata,op1.mdim*sizeof(T));
        }

        return *this;
    }
    //@}

    /// @name Core memory operations for Matrix2D
    //@{
    /** Clear.
     */
    void clear()
    {
        coreDeallocate();
        coreInit();
    }

    /** Core init from mapped file.
     * Offset is in bytes. */
    void coreInit(const FileName &fn, int Ydim, int Xdim, size_t offset=0)
    {
#ifdef XMIPP_MMAP

        mdimx=Xdim;
        mdimy=Ydim;
        mdim=mdimx*mdimy;
        destroyData=false;
        mappedData=true;
        fdMap = open(fn.c_str(),  O_RDWR, S_IREAD | S_IWRITE);
        if (fdMap == -1)
            REPORT_ERROR(ERR_IO_NOTOPEN,fn);
        const size_t pagesize=sysconf(_SC_PAGESIZE);
        size_t offsetPages=(offset/pagesize)*pagesize;
        size_t offsetDiff=offset-offsetPages;
        if ( (mdataOriginal = (char*) mmap(0,Ydim*Xdim*sizeof(T)+offsetDiff, PROT_READ | PROT_WRITE, MAP_SHARED, fdMap, offsetPages)) == MAP_FAILED )
            REPORT_ERROR(ERR_MMAP_NOTADDR,(String)"mmap failed "+integerToString(errno));
        mdata=(T*)(mdataOriginal+offsetDiff);
#else

        REPORT_ERROR(ERR_MMAP,"Mapping not supported in Windows");
#endif

    }

    /** Core init.
     * Initialize everything to 0
     */
    void coreInit()
    {
        mdimx=mdimy=mdim=0;
        mdata=NULL;
        mdataOriginal=NULL;
        destroyData=true;
        mappedData=false;
        fdMap=-1;
    }

    /** Core allocate.
     */
    void coreAllocate( int _mdimy, int _mdimx)
    {
        if (_mdimy <= 0 ||_mdimx<=0)
        {
            clear();
            return;
        }

        mdimx=_mdimx;
        mdimy=_mdimy;
        mdim=_mdimx*_mdimy;
        mdata = new T [mdim];
        mdataOriginal = NULL;
        mappedData=false;
        fdMap=-1;
        if (mdata == NULL)
            REPORT_ERROR(ERR_MEM_NOTENOUGH, "coreAllocate: No space left");
    }

    /** Core deallocate.
     * Free all mdata.
     */
    void coreDeallocate()
    {
        if (mdata != NULL && destroyData)
            delete[] mdata;
        if (mappedData)
        {
#ifdef XMIPP_MMAP
            munmap(mdataOriginal,mdimx*mdimy*sizeof(T));
            close(fdMap);
#else

            REPORT_ERROR(ERR_MMAP,"Mapping not supported in Windows");
#endif

        }
        mdata=NULL;
        mdataOriginal=NULL;
    }
    //@}

    /// @name Size and shape of Matrix2D
    //@{
    /** Resize to a given size
     */
    void resize(size_t Ydim, size_t Xdim, bool noCopy=false)
    {

        if (Xdim == mdimx && Ydim == mdimy)
            return;

        if (Xdim <= 0 || Ydim <= 0)
        {
            clear();
            return;
        }

        T * new_mdata;
        size_t YXdim=Ydim*Xdim;

        try
        {
            new_mdata = new T [YXdim];
        }
        catch (std::bad_alloc &)
        {
            REPORT_ERROR(ERR_MEM_NOTENOUGH, "Allocate: No space left");
        }

        // Copy needed elements, fill with 0 if necessary
        if (!noCopy)
        {
        	T zero=0; // Useful for complexes
            for (size_t i = 0; i < Ydim; i++)
                for (size_t j = 0; j < Xdim; j++)
                {
                    T *val=NULL;
                    if (i >= mdimy)
                        val = &zero;
                    else if (j >= mdimx)
                        val = &zero;
                    else
                        val = &mdata[i*mdimx + j];
                    new_mdata[i*Xdim+j] = *val;
                }
        }
        else
            memset(new_mdata,0,YXdim*sizeof(T));

        // deallocate old vector
        coreDeallocate();

        // assign *this vector to the newly created
        mdata = new_mdata;
        mdimx = Xdim;
        mdimy = Ydim;
        mdim = Xdim * Ydim;
        mappedData = false;
    }

    /** Resize according to a pattern.
     *
     * This function resize the actual array to the same size and origin
     * as the input pattern. If the actual array is larger than the pattern
     * then the trailing values are lost, if it is smaller then 0's are
     * added at the end
     *
     * @code
     * v2.resize(v1);
     * // v2 has got now the same structure as v1
     * @endcode
     */
    template<typename T1>
    void resize(const Matrix2D<T1> &v)
    {
        if (mdimx != v.mdimx || mdimy != v.mdimy)
            resize(v.mdimy, v.mdimx);
    }

    /** Resize to a given size (don't copy old elements)
     */
    inline void resizeNoCopy(int Ydim, int Xdim)
    {
        resize(Ydim, Xdim, true);
    }

    /** Resize according to a pattern.
     *  Do not copy old elements.
     */
    template<typename T1>
    inline void resizeNoCopy(const Matrix2D<T1> &v)
    {
        if (mdimx != v.mdimx || mdimy != v.mdimy)
            resize(v.mdimy, v.mdimx, true);
    }

    /** Map to file.
     * The matrix is mapped to a file. The file is presumed to be already created with
     * enough space for a matrix of size Ydim x Xdim.
     * Offset is in bytes.
     */
    void mapToFile(const FileName &fn, int Ydim, int Xdim, size_t offset=0)
    {
        if (mdata!=NULL)
            clear();

#ifdef XMIPP_MMAP

        coreInit(fn,Ydim,Xdim,offset);
#else

        resizeNoCopy(Ydim, Xdim);
#endif

    }

    /** Extract submatrix and assign to this object.
     */
    void submatrix(int i0, int j0, int iF, int jF)
    {
        if (i0 < 0 || j0 < 0 || iF >= MAT_YSIZE(*this) || jF >= MAT_XSIZE(*this))
            REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS,"Submatrix indexes out of bounds");
        Matrix2D<T> result(iF - i0 + 1, jF - j0 + 1);

        FOR_ALL_ELEMENTS_IN_MATRIX2D(result)
        MAT_ELEM(result, i, j) = MAT_ELEM(*this, i+i0, j+j0);

        *this = result;
    }

    /** Same shape.
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument
     */
    template <typename T1>
    bool sameShape(const Matrix2D<T1>& op) const
    {
        return ((mdimx == op.mdimx) && (mdimy == op.mdimy));
    }

    /** X dimension
     *
     * Returns X dimension
     */
    inline size_t Xdim() const
    {
        return mdimx;
    }

    /** Y dimension
     *
     * Returns Y dimension
     */
    inline size_t Ydim() const
    {
        return mdimy;
    }
    //@}

    /// @name Initialization of Matrix2D values
    //@{
    /** Same value in all components.
     *
     * The constant must be of a type compatible with the array type, ie,
     * you cannot  assign a double to an integer array without a casting.
     * It is not an error if the array is empty, then nothing is done.
     *
     * @code
     * v.initConstant(3.14);
     * @endcode
     */
    void initConstant(T val)
    {
        for (size_t j = 0; j < mdim; j++)
            mdata[j] = val;
    }

    /** Initialize to zeros with a given size.
     */
    void initConstant(size_t Ydim, size_t Xdim, T val)
    {
        if (mdimx!=Xdim || mdimy!=Ydim)
            resizeNoCopy(Ydim, Xdim);
        initConstant(val);
    }

    /** Initialize to zeros with current size.
     *
     * All values are set to 0. The current size and origin are kept. It is not
     * an error if the array is empty, then nothing is done.
     *
     * @code
     * v.initZeros();
     * @endcode
     */
    void initZeros()
    {
        memset(mdata,0,mdimx*mdimy*sizeof(T));
    }

    /** Initialize to zeros with a given size.
     */
    void initZeros(size_t Ydim, size_t Xdim)
    {
        if (mdimx!=Xdim || mdimy!=Ydim)
            resizeNoCopy(Ydim, Xdim);
        memset(mdata,0,mdimx*mdimy*sizeof(T));
    }

    /** Initialize to zeros following a pattern.
      *
      * All values are set to 0, and the origin and size of the pattern are
      * adopted.
      *
      * @code
      * v2.initZeros(v1);
      * @endcode
      */
    template <typename T1>
    void initZeros(const Matrix2D<T1>& op)
    {
        if (mdimx!=op.mdimx || mdimy!=op.mdimy)
            resizeNoCopy(op);
        memset(mdata,0,mdimx*mdimy*sizeof(T));
    }

    /** Initialize to random random numbers, uniform or gaussian
     */
    void initRandom(size_t Ydim, size_t Xdim, double op1, double op2, RandomMode mode = RND_UNIFORM)
    {
        if (mdimx!=Xdim || mdimy!=Ydim)
            resizeNoCopy(Ydim, Xdim);
        for (size_t j = 0; j < mdim; j++)
           mdata[j] = static_cast< T > (mode == RND_UNIFORM ? rnd_unif(op1, op2) : rnd_gaus(op1, op2));
    }

    /** Initialize to gaussian numbers */
    void initGaussian(int Ydim, int Xdim, double op1=0., double op2=1.)
    {
      initRandom(Ydim, Xdim, op1, op2, RND_GAUSSIAN);
    }

    /** 2D Identity matrix of current size
     *
     * If actually the matrix is not squared then an identity matrix is
     * generated of size (Xdim x Xdim).
     *
     * @code
     * m.initIdentity();
     * @endcode
     */
    void initIdentity()
    {
        initIdentity(MAT_XSIZE(*this));
    }

    /** 2D Identity matrix of a given size
     *
     * A (dim x dim) identity matrix is generated.
     *
     * @code
     * m.initIdentity(3);
     * @endcode
     */
    void initIdentity(int dim)
    {
        initZeros(dim, dim);
        for (int i = 0; i < dim; i++)
            MAT_ELEM(*this,i,i) = 1;
    }
    /** 2D gaussian matrix of a given size and with a given variance. The amplitude of the Gaussian is set to 1.
     *
     * A (dim x dim) gaussian matrix is generated.
     *
     * @code
     * m.initGaussian(3,1);
     * @endcode
     */
    void initGaussian(int dim, double var)
    {
        double center = ((double)dim)/2;
        initZeros(dim, dim);
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                MAT_ELEM(*this,i,j) = std::exp(-( (i-center)*(i-center)+(j-center)*(j-center) )/(2*var*var));
    }
    //@}

    /// @name Operators for Matrix2D
    //@{

    /** Matrix element access
     */
    T& operator()(int i, int j) const
    {
        return MAT_ELEM((*this),i,j);
    }
    /** Parenthesis operator for phyton
    */
    void setVal(T val,int y, int x)
    {
        MAT_ELEM((*this),y,x)=val;
    }
    /** Parenthesis operator for phyton
    */
    T getVal( int y, int x) const
    {
        return MAT_ELEM((*this),y,x);
    }

    /** v3 = v1 * k.
     */
    Matrix2D<T> operator*(T op1) const
    {
        Matrix2D<T> tmp(*this);
        for (size_t i=0; i < mdim; i++)
            tmp.mdata[i] = mdata[i] * op1;
        return tmp;
    }

    /** v3 = v1 / k.
     */
    Matrix2D<T> operator/(T op1) const
    {
        Matrix2D<T> tmp(*this);
        for (size_t i=0; i < mdim; i++)
            tmp.mdata[i] = mdata[i] / op1;
        return tmp;
    }

    /** v3 = k * v2.
     */
    friend Matrix2D<T> operator*(T op1, const Matrix2D<T>& op2)
    {
        Matrix2D<T> tmp(op2);
        for (size_t i=0; i < op2.mdim; i++)
            tmp.mdata[i] = op1 * op2.mdata[i];
        return tmp;
    }

    /** v3 *= k.
      */
    void operator*=(T op1)
    {
        for (size_t i=0; i < mdim; i++)
            mdata[i] *= op1;
    }

    /** v3 /= k.
      */
    void operator/=(T op1)
    {
        for (size_t i=0; i < mdim; i++)
            mdata[i] /= op1;
    }

    /** Matrix by vector multiplication
     *
     * @code
     * v2 = A*v1;
     * @endcode
     */
    Matrix1D<T> operator*(const Matrix1D<T>& op1) const
    {
        Matrix1D<T> result;

        if (mdimx != VEC_XSIZE(op1))
            REPORT_ERROR(ERR_MATRIX_SIZE, "Not compatible sizes in matrix by vector");

        if (!op1.isCol())
            REPORT_ERROR(ERR_MATRIX, "Vector is not a column");

        result.initZeros(mdimy);
        for (size_t i = 0; i < mdimy; i++)
            for (size_t j = 0; j < mdimx; j++)
                VEC_ELEM(result,i) += MAT_ELEM(*this,i, j) * VEC_ELEM(op1,j);

        result.setCol();
        return result;
    }

    /** Matrix by Matrix multiplication
     *
     * @code
     * C = A*B;
     * @endcode
     */
    Matrix2D<T> operator*(const Matrix2D<T>& op1) const
    {
        Matrix2D<T> result;
        if (mdimx != op1.mdimy)
            REPORT_ERROR(ERR_MATRIX_SIZE, "Not compatible sizes in matrix multiplication");

        result.initZeros(mdimy, op1.mdimx);
        for (size_t i = 0; i < mdimy; i++)
            for (size_t j = 0; j < op1.mdimx; j++)
                for (size_t k = 0; k < mdimx; k++)
                    MAT_ELEM(result,i, j) += MAT_ELEM(*this,i, k) * MAT_ELEM(op1, k, j);
        return result;
    }

    /** Matrix summation
     *
     * @code
     * C = A + B;
     * @endcode
     */
    Matrix2D<T> operator+(const Matrix2D<T>& op1) const
    {
        Matrix2D<T> result;
        if (mdimx != op1.mdimx || mdimy != op1.mdimy)
            REPORT_ERROR(ERR_MATRIX_SIZE, "operator+: Not same sizes in matrix summation");

        result.initZeros(mdimy, mdimx);
        for (size_t i = 0; i < mdimy; i++)
            for (size_t j = 0; j < mdimx; j++)
                result(i, j) = (*this)(i, j) + op1(i, j);

        return result;
    }

    /** Matrix summation
     *
     * @code
     * A += B;
     * @endcode
     */
    void operator+=(const Matrix2D<T>& op1) const
    {
        if (mdimx != op1.mdimx || mdimy != op1.mdimy)
            REPORT_ERROR(ERR_MATRIX_SIZE, "operator+=: Not same sizes in matrix summation");

        for (size_t i = 0; i < mdimy; i++)
            for (size_t j = 0; j < mdimx; j++)
                MAT_ELEM(*this,i, j) += MAT_ELEM(op1, i, j);
    }

    /** Matrix subtraction
     *
     * @code
     * C = A - B;
     * @endcode
     */
    Matrix2D<T> operator-(const Matrix2D<T>& op1) const
    {
        Matrix2D<T> result;
        if (mdimx != op1.mdimx || mdimy != op1.mdimy)
            REPORT_ERROR(ERR_MATRIX_SIZE, "operator-: Not same sizes in matrix summation");

        result.initZeros(mdimy, mdimx);
        for (size_t i = 0; i < mdimy; i++)
            for (size_t j = 0; j < mdimx; j++)
                result(i, j) = (*this)(i, j) - op1(i, j);

        return result;
    }

    /** Matrix substraction
     *
     * @code
     * A -= B;
     * @endcode
     */
    void operator-=(const Matrix2D<T>& op1) const
    {
        if (mdimx != op1.mdimx || mdimy != op1.mdimy)
            REPORT_ERROR(ERR_MATRIX_SIZE, "operator-=: Not same sizes in matrix summation");

        for (size_t i = 0; i < mdimy; i++)
            for (size_t j = 0; j < mdimx; j++)
                MAT_ELEM(*this,i, j) -= MAT_ELEM(op1, i, j);
    }
    /** Equality.
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument and the same values (within accuracy).
     */
    bool equal(const Matrix2D<T>& op,
               double accuracy = XMIPP_EQUAL_ACCURACY) const
    {
        if (!sameShape(op))
            return false;
        for (size_t i = 0; i < mdimy; i++)
            for (size_t j = 0; j < mdimx; j++)
                if (fabs( MAT_ELEM(*this,i,j) - MAT_ELEM(op,i,j) ) > accuracy)
                    return false;
        return true;
    }
    //@}

    /// @name Utilities for Matrix2D
    //@{
    /** Maximum of the values in the array.
      *
      * The returned value is of the same type as the type of the array.
      */
    T computeMax() const
    {
        if (mdim <= 0)
            return static_cast< T >(0);

        T maxval = mdata[0];
        for (size_t n = 0; n < mdim; n++)
            if (mdata[n] > maxval)
                maxval = mdata[n];
        return maxval;
    }

    /** Minimum of the values in the array.
       *
       * The returned value is of the same type as the type of the array.
       */
    T computeMin() const
    {
        if (mdim <= 0)
            return static_cast< T >(0);

        T minval = mdata[0];
        for (size_t n = 0; n < mdim; n++)
            if (mdata[n] < minval)
                minval = mdata[n];
        return minval;
    }

    /** Produce a 2D array suitable for working with Numerical Recipes
    *
    * This function must be used only as a preparation for routines which need
    * that the first physical index is 1 and not 0 as it usually is in C. New
    * memory is needed to hold the new double pointer array.
    */
    T** adaptForNumericalRecipes() const
    {
        T** m = NULL;
        ask_Tmatrix(m, 1, mdimy, 1, mdimx);

        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < mdimx; j++)
                m[i+1][j+1] = mdata[i*mdimx + j];

        return m;
    }

    /** Produce a 1D pointer suitable for working with Numerical Recipes (2)
     *
     * This function meets the same goal as the one before, however this one
     * work with 2D arrays as a single pointer. The first element of the array
     * is pointed by result[1*Xdim+1], and in general result[i*Xdim+j]
     */
    T* adaptForNumericalRecipes2() const
    {
        return mdata - 1 - mdimx;
    }

    /** Load 2D array from numerical recipes result.
     */
    void loadFromNumericalRecipes(T** m, int Ydim, int Xdim)
    {
        if (mdimx!=Xdim || mdimy!=Ydim)
            resizeNoCopy(Ydim, Xdim);

        for (int i = 1; i <= Ydim; i++)
            for (int j = 1; j <= Xdim; j++)
                MAT_ELEM(*this,i - 1, j - 1) = m[i][j];
    }

    /** Kill a 2D array produced for numerical recipes
     *
     * The allocated memory is freed.
     */
    void killAdaptationForNumericalRecipes(T** m) const
    {
        free_Tmatrix(m, 1, mdimy, 1, mdimx);
    }

    /** Kill a 2D array produced for numerical recipes, 2.
     *
     * Nothing needs to be done.
     */
    void killAdaptationForNumericalRecipes2(T** m) const
        {}

    /** Read this matrix from file.
     *  The matrix is assumed to be already resized.
      */
    void read(const FileName &fn)
    {
        std::ifstream fhIn;
        fhIn.open(fn.c_str());
        if (!fhIn)
            REPORT_ERROR(ERR_IO_NOTEXIST,fn);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
        fhIn >> MAT_ELEM(*this,i,j);
        fhIn.close();
    }

    /** Write this matrix to file
      */
    void write(const FileName &fn) const
    {
        std::ofstream fhOut;
        fhOut.open(fn.c_str());
        if (!fhOut)
            REPORT_ERROR(ERR_IO_NOTOPEN,(std::string)"write: Cannot open "+fn+" for output");
        fhOut << *this;
        fhOut.close();
    }

    /** Show matrix
      */
    friend std::ostream& operator<<(std::ostream& ostrm, const Matrix2D<T>& v)
    {
        if (v.Xdim() == 0 || v.Ydim() == 0)
            ostrm << "NULL matrix\n";
        else
        {
            ostrm << std::endl;
            double max_val = v.computeMax();
            int prec = bestPrecision(max_val, 10);

            for (size_t i = 0; i < v.Ydim(); i++)
            {
                for (size_t j = 0; j < v.Xdim(); j++)
                {
                    ostrm << floatToString((double) v(i, j), 10, prec) << ' ';
                }
                ostrm << std::endl;
            }
        }

        return ostrm;
    }

    /** Makes a matrix from a vector
     *
     * The origin of the matrix is set such that it has one of the index origins
     * (X or Y) to the same value as the vector, and the other set to 0
     * according to the shape.
     *
     * @code
     * Matrix2D< double > m = fromVector(v);
     * @endcode
     */
    void fromVector(const Matrix1D<T>& op1)
    {
        // Null vector => Null matrix
        if (op1.size() == 0)
        {
            clear();
            return;
        }

        // Look at shape and copy values
        if (op1.isRow())
        {
            if (mdimy!=1 || mdimx!=VEC_XSIZE(op1))
                resizeNoCopy(1, VEC_XSIZE(op1));

            for (size_t j = 0; j < VEC_XSIZE(op1); j++)
                MAT_ELEM(*this,0, j) = VEC_ELEM(op1,j);
        }
        else
        {
            if (mdimy!=1 || mdimx!=VEC_XSIZE(op1))
                resizeNoCopy(VEC_XSIZE(op1), 1);

            for (size_t i = 0; i < VEC_XSIZE(op1); i++)
                MAT_ELEM(*this, i, 0) = VEC_ELEM(op1,i);
        }
    }

    /** Makes a vector from a matrix
     *
     * An exception is thrown if the matrix is not a single row or a single
     * column. The origin of the vector is set according to the one of the
     * matrix.
     *
     * @code
     * Matrix1D< double > v;
     * m.toVector(v);
     * @endcode
     */
    void toVector(Matrix1D<T>& op1) const
    {
        // Null matrix => Null vector
        if (mdimx == 0 || mdimy == 0)
        {
            op1.clear();
            return;
        }

        // If matrix is not a vector, produce an error
        if (!(mdimx == 1 || mdimy == 1))
            REPORT_ERROR(ERR_MATRIX_DIM,
                         "toVector: Matrix cannot be converted to vector");

        // Look at shape and copy values
        if (mdimy == 1)
        {
            // Row vector
            if (VEC_XSIZE(op1)!=mdimx)
                op1.resizeNoCopy(mdimx);

            memcpy(&VEC_ELEM(op1,0),&MAT_ELEM(*this,0,0),mdimx*sizeof(double));

            op1.setRow();
        }
        else
        {
            // Column vector
            if (VEC_XSIZE(op1)!=mdimy)
                op1.resizeNoCopy(mdimy);

            for (size_t i = 0; i < mdimy; i++)
                VEC_ELEM(op1,i) = MAT_ELEM(*this,i, 0);

            op1.setCol();
        }
    }

    /**Copy matrix to stl::vector
     */
    void copyToVector(std::vector<T> &v)
    {
        v.assign(mdata, mdata+mdim);
    }
    /**Copy stl::vector to matrix
      */
    void copyFromVector(std::vector<T> &v,int Xdim, int Ydim)
    {
        if (mdimx!=Xdim || mdimy!=Ydim)
            resizeNoCopy(Ydim, Xdim);
        copy( v.begin(), v.begin()+v.size(), mdata);
    }

    /** Get row
     *
     * This function returns a row vector corresponding to the choosen
     * row inside the nth 2D matrix, the numbering of the rows is also
     * logical not physical.
     *
     * @code
     * std::vector< double > v;
     * m.getRow(-2, v);
     * @endcode
     */
    void getRow(size_t i, Matrix1D<T>& v) const
    {
        if (mdimx == 0 || mdimy == 0)
        {
            v.clear();
            return;
        }

        if (i < 0 || i >= mdimy)
            REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, "getRow: Matrix subscript (i) greater than matrix dimension");

        if (VEC_XSIZE(v)!=mdimx)
            v.resizeNoCopy(mdimx);
        memcpy(&VEC_ELEM(v,0),&MAT_ELEM(*this,i,0),mdimx*sizeof(double));

        v.setRow();
    }

    /** Get Column
     *
     * This function returns a column vector corresponding to the
     * choosen column.
     *
     * @code
     * std::vector< double > v;
     * m.getCol(-1, v);
     * @endcode
     */
    void getCol(size_t j, Matrix1D<T>& v) const
    {
        if (mdimx == 0 || mdimy == 0)
        {
            v.clear();
            return;
        }

        if (j < 0 || j >= mdimx)
            REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS,"getCol: Matrix subscript (j) greater than matrix dimension");

        if (VEC_XSIZE(v)!=mdimy)
            v.resizeNoCopy(mdimy);
        for (size_t i = 0; i < mdimy; i++)
            VEC_ELEM(v,i) = MAT_ELEM(*this,i, j);

        v.setCol();
    }

    /** Set Row
     *
     * This function sets a row vector corresponding to the choosen row in the 2D Matrix
     *
     * @code
     * m.setRow(-2, m.row(1)); // Copies row 1 in row -2
     * @endcode
     */
    void setRow(size_t i, const Matrix1D<T>& v)
    {
        if (mdimx == 0 || mdimy == 0)
            REPORT_ERROR(ERR_MATRIX_EMPTY, "setRow: Target matrix is empty");

        if (i < 0 || i >= mdimy)
            REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, "setRow: Matrix subscript (i) out of range");

        if (VEC_XSIZE(v) != mdimx)
            REPORT_ERROR(ERR_MATRIX_SIZE,
                         "setRow: Vector dimension different from matrix one");

        if (!v.isRow())
            REPORT_ERROR(ERR_MATRIX_DIM, "setRow: Not a row vector in assignment");

        memcpy(&MAT_ELEM(*this,i,0),&VEC_ELEM(v,0),mdimx*sizeof(double));
    }

    /** Set Column
     *
     * This function sets a column vector corresponding to the choosen column
     * inside matrix.
     *
     * @code
     * m.setCol(0, (m.row(1)).transpose()); // Copies row 1 in column 0
     * @endcode
     */
    void setCol(size_t j, const Matrix1D<T>& v)
    {
        if (mdimx == 0 || mdimy == 0)
            REPORT_ERROR(ERR_MATRIX_EMPTY, "setCol: Target matrix is empty");

        if (j < 0 || j>= mdimx)
            REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, "setCol: Matrix subscript (j) out of range");

        if (VEC_XSIZE(v) != mdimy)
            REPORT_ERROR(ERR_MATRIX_SIZE,
                         "setCol: Vector dimension different from matrix one");

        if (!v.isCol())
            REPORT_ERROR(ERR_MATRIX_DIM, "setCol: Not a column vector in assignment");

        for (size_t i = 0; i < mdimy; i++)
            MAT_ELEM(*this,i, j) = VEC_ELEM(v,i);
    }

    /** Determinant of a matrix
     *
     * An exception is thrown if the matrix is not squared or it is empty.
     *
     * @code
     * double det = m.det();
     * @endcode
     */
    T det() const
    {
        // (see Numerical Recipes, Chapter 2 Section 5)
        if (mdimx == 0 || mdimy == 0)
            REPORT_ERROR(ERR_MATRIX_EMPTY, "determinant: Matrix is empty");

        if (mdimx != mdimy)
            REPORT_ERROR(ERR_MATRIX_SIZE, "determinant: Matrix is not squared");

        for (size_t i = 0; i < mdimy; i++)
        {
            bool all_zeros = true;
            for (size_t j = 0; j < mdimx; j++)
                if (fabs(MAT_ELEM((*this),i, j)) > XMIPP_EQUAL_ACCURACY)
                {
                    all_zeros = false;
                    break;
                }

            if (all_zeros)
                return 0;
        }

        // Perform decomposition
        Matrix1D< int > indx;
        T d;
        Matrix2D<T> LU;
        ludcmp(*this, LU, indx, d);

        // Calculate determinant
        for (size_t i = 0; i < mdimx; i++)
            d *= (T) MAT_ELEM(LU,i , i);

        return d;
    }
    ///determinat of 3x3 matrix
    T det3x3() const
    {
        return (
                   dMij(*this,0,0)*( dMij(*this,2,2)*dMij(*this,1,1)-
                                     dMij(*this,2,1)*dMij(*this,1,2) )-
                   dMij(*this,1,0)*( dMij(*this,2,2)*dMij(*this,0,1)-
                                     dMij(*this,2,1)*dMij(*this,0,2) )+
                   dMij(*this,2,0)*( dMij(*this,1,2)*dMij(*this,0,1)-
                                     dMij(*this,1,1)*dMij(*this,0,2) )
               );
    }
    /** Algebraic transpose of a Matrix
     *
     * You can use the transpose in as complex expressions as you like. The
     * origin of the vector is not changed.
     *
     * @code
     * v2 = v1.transpose();
     * @endcode
     */
    Matrix2D<T> transpose() const
    {
        Matrix2D<T> result(mdimx, mdimy);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(result)
        MAT_ELEM(result,i,j) = MAT_ELEM((*this),j,i);
        return result;
    }

    /** Inverse of a matrix
     *
     * The matrix is inverted using a SVD decomposition. In fact the
     * pseudoinverse is returned.
     *
     * @code
     * Matrix2D< double > m1_inv;
     * m1.inv(m1_inv);
     * @endcode
     */
    void inv(Matrix2D<double>& result) const
    {
        if (mdimx == 0 || mdimy == 0)
            REPORT_ERROR(ERR_MATRIX_EMPTY, "Inverse: Matrix is empty");
        result.initZeros(mdimx, mdimy);
        SPEED_UP_temps0;
        if (mdimx==2)
        {
            M2x2_INV(result,*this);
        }
        else if (mdimx==3)
        {
            M3x3_INV(result,*this);
        }
        else if (mdimx==4)
        {
            M4x4_INV(result,*this);
        }
        else
        {
            // Perform SVD decomposition
            Matrix2D< double > u, v;
            Matrix1D< double > w;
            svdcmp(*this, u, w, v); // *this = U * W * V^t

            double tol = computeMax() * XMIPP_MAX(mdimx, mdimy) * 1e-14;

            // Compute W^-1
            bool invertible = false;
            FOR_ALL_ELEMENTS_IN_MATRIX1D(w)
            {
                if (fabs(VEC_ELEM(w,i)) > tol)
                {
                    VEC_ELEM(w,i) = 1.0 / VEC_ELEM(w,i);
                    invertible = true;
                }
                else
                    VEC_ELEM(w,i) = 0.0;
            }

            if (!invertible)
                return;

            // Compute V*W^-1
            FOR_ALL_ELEMENTS_IN_MATRIX2D(v)
            MAT_ELEM(v,i,j) *= VEC_ELEM(w,j);

            // Compute Inverse
            for (size_t i = 0; i < mdimx; i++)
                for (size_t j = 0; j < mdimy; j++)
                    for (size_t k = 0; k < mdimx; k++)
                        MAT_ELEM(result,i,j) += MAT_ELEM(v,i,k) * MAT_ELEM(u,j,k);
        }
    }

    /** Perform SVD decomposition
     *  *this = U * W * V^t
     * */

    void svd(Matrix2D<double> &U, Matrix1D<double> &W, Matrix2D<double> &V) const
    {
      svdcmp(*this, U, W, V);
    }

    /** Perform SVD decomposition and add an index vector
     * with the descending order of eigenvalues
     */
    void eigs(Matrix2D<double> &U, Matrix1D<double> &W, Matrix2D<double> &V, Matrix1D<int> &indexes) const
    {
      svdcmp(*this, U, W, V);
      indexes.resizeNoCopy(W);
      indexes.enumerate();

      double dAux;
      int iAux;

      FOR_ALL_ELEMENTS_IN_MATRIX1D(W)
      {
          for (int j = i; j > 0 && dMi(W, j) > dMi(W, j-1); --j)
          {
            VEC_SWAP(W, j, j-1, dAux);
            VEC_SWAP(indexes, j, j-1, iAux);
          }
      }
    }

    /** Inverse of a matrix
     */
    Matrix2D<double> inv() const
    {
        Matrix2D<double> result;
        inv(result);

        return result;
    }

    /** Inverse the current matrix
     */
    void selfInverse()
    {
        Matrix2D<T> auxMatrix(*this);
        auxMatrix.inv(*this);
    }

    /** True if the matrix is identity
     *
     * @code
     * if (m.isIdentity())
     *     std::cout << "The matrix is identity\n";
     * @endcode
     */
    bool isIdentity() const
    {
        for (size_t i = 0; i < mdimy; i++)
            for (size_t j = 0; j < mdimx; j++)
                if (i != j)
                {
                    if (MAT_ELEM(*this,i,j)!=0)
                        return false;
                }
                else
                {
                    if (MAT_ELEM(*this,i,j)!=1)
                        return false;
                }
        return true;
    }
    //@}
};

typedef Matrix2D<double> DMatrix;
typedef Matrix2D<int> IMatrix;

template<typename T>
bool operator==(const Matrix2D<T>& op1, const Matrix2D<T>& op2)
{
    return op1.equal(op2);
}


/**@name Matrix Related functions
 * These functions are not methods of Matrix2D
 */
//@{
/** LU Decomposition
 */
template<typename T>
void ludcmp(const Matrix2D<T>& A, Matrix2D<T>& LU, Matrix1D< int >& indx, T& d)
{
    LU = A;
    if (VEC_XSIZE(indx)!=A.mdimx)
        indx.resizeNoCopy(A.mdimx);
    ludcmp(LU.adaptForNumericalRecipes2(), A.mdimx,
           indx.adaptForNumericalRecipes(), &d);
}


/** LU Backsubstitution
 */
template<typename T>
void lubksb(const Matrix2D<T>& LU, Matrix1D< int >& indx, Matrix1D<T>& b)
{
    lubksb(LU.adaptForNumericalRecipes2(), indx.size(),
           indx.adaptForNumericalRecipes(),
           b.adaptForNumericalRecipes());
}

/** SVD Backsubstitution
 */
void svbksb(Matrix2D< double >& u,
            Matrix1D< double >& w,
            Matrix2D< double >& v,
            Matrix1D< double >& b,
            Matrix1D< double >& x);

//#define VIA_NR
#define VIA_BILIB
/** SVD Decomposition
 */
template<typename T>
void svdcmp(const Matrix2D< T >& a,
            Matrix2D< double >& u,
            Matrix1D< double >& w,
            Matrix2D< double >& v)
{
    // svdcmp only works with double
    typeCast(a, u);

    // Set size of matrices
    w.initZeros(u.mdimx);
    v.initZeros(u.mdimx, u.mdimx);

    // Call to the numerical recipes routine
#ifdef VIA_NR

    svdcmp(u.mdata,
           u.mdimy, u.mdimx,
           w.vdata,
           v.mdata);
#endif

#ifdef VIA_BILIB

    int status;
    SingularValueDecomposition(u.mdata,
                               u.mdimy, u.mdimx,
                               w.vdata,
                               v.mdata,
                               5000, &status);
#endif
}
#undef VIA_NR
#undef VIA_BILIB

/** Generalized eigenvector decomposition.
 * Solves the problem Av=dBv.
 * The decomposition is such that A=B P D P^-1. A and B must be square matrices of the same size.
 */
void generalizedEigs(const Matrix2D<double> &A, const Matrix2D<double> &B, Matrix1D<double> &D, Matrix2D<double> &P);

/** First eigenvectors of a real, symmetric matrix.
 * Solves the problem Av=dv.
 * Only the eigenvectors of the largest M eigenvalues are returned as columns of P
 */
void firstEigs(const Matrix2D<double> &A, size_t M, Matrix1D<double> &D, Matrix2D<double> &P);

/** Last eigenvectors of a real, symmetric matrix.
 * Solves the problem Av=dv.
 * Only the eigenvectors of the smallest M eigenvalues are returned as columns of P
 */
void lastEigs(const Matrix2D<double> &A, size_t M, Matrix1D<double> &D, Matrix2D<double> &P);

/** Find connected components of a graph.
 * Assuming that the matrix G represents an undirected graph (of size NxN), this function returns a vector (of size N) that indicates
 * for each element which is the number of its connected component.
 */
void connectedComponentsOfUndirectedGraph(const Matrix2D<double> &G, Matrix1D<int> &component);

/** Conversion from one type to another.
 *
 * If we have an integer array and we need a double one, we can use this
 * function. The conversion is done through a type casting of each element
 * If n >= 0, only the nth volumes will be converted, otherwise all NSIZE volumes
 */
template<typename T1, typename T2>
void typeCast(const Matrix2D<T1>& v1,  Matrix2D<T2>& v2)
{
    if (v1.mdim == 0)
    {
        v2.clear();
        return;
    }

    if (v1.mdimx!=v2.mdimx || v1.mdimy!=v2.mdimy)
        v2.resizeNoCopy(v1);
    for (unsigned long int n = 0; n < v1.mdim; n++)
        v2.mdata[n] = static_cast< T2 > (v1.mdata[n]);
}

/** Conversion from one type to another.
 * In some cases, the two types are the same. So a faster way is simply by assignment.
 */
template<typename T1>
void typeCast(const Matrix2D<T1>& v1,  Matrix2D<T1>& v2)
{
    v2=v1;
}

/** Helper class for solving linear systems */
class PseudoInverseHelper
{
public:
    Matrix2D<double> A, AtA, AtAinv;
    Matrix1D<double> Atb, b, bpredicted;
};

/** Helper class for solving Weighted Least Squares */
class WeightedLeastSquaresHelper: public PseudoInverseHelper
{
public:
    Matrix1D<double> w; //Weights
};

/** Solve Linear system Ax=b with pseudoinverse.
 * A and b must be set inside the PseudoInverseHelper, the rest of the
 * fields in PseudoInverseHelper are used by this routine to avoid
 * several allocation/deallocations */
void solveLinearSystem(PseudoInverseHelper &h, Matrix1D<double> &result);

/** Solve Weighted least square problem Ax=b with pseudoinverse and weights w.
 * A, w and b must be set inside the WeightedLeastSquaresHelper, the rest of the
 * fields in WeightedLeastSquaresHelper are used by this routine to avoid
 * several allocation/deallocations.
 *
 * The normal equations of this problem are A^t W A x = A^t W b,
 * where W is a diagonal matrix whose entries are in the vector w. */
void weightedLeastSquares(WeightedLeastSquaresHelper &h, Matrix1D<double> &result);

/** Solve Weighted least square problem Ax=b and weights w, with RANSAC.
 * Tol is a tolerance value: if an equation is fulfilled with an error smaller than tol,
 * then it is considered to fulfill the model. Niter is the number of RANSAC iterations to perform.
 * The outlier fraction is the fraction of equations that, at maximum, can be considered as outliers.
 */
void ransacWeightedLeastSquares(WeightedLeastSquaresHelper &h, Matrix1D<double> &result,
		double tol, int Niter=10000, double outlierFraction=0.25, int Nthreads=1);

/** Sparse element.
 *  This class is used to create the SparseMatrices. */
class SparseElement
{
public:
    size_t i;
    size_t j;
    double value;
};

/** Function to sort the sparse elements by their i,j position */
inline bool operator< (const SparseElement& _x, const SparseElement& _y)
{
    return ( _x.i < _y.i) || (( _x.i== _y.i) && ( _x.j< _y.j));
}

//@}
//@}

// Implementation of the vector*matrix
// Documented in matrix1D.h
template<typename T>
Matrix1D<T> Matrix1D<T>::operator*(const Matrix2D<T>& M)
{
    Matrix1D<T> result;

    if (VEC_XSIZE(*this) != MAT_YSIZE(M))
        REPORT_ERROR(ERR_MATRIX_SIZE, "Not compatible sizes in matrix by vector");

    if (!isRow())
        REPORT_ERROR(ERR_MATRIX_DIM, "Vector is not a row");

    result.initZeros(MAT_XSIZE(M));
    for (size_t j = 0; j < MAT_XSIZE(M); j++)
        for (size_t i = 0; i < MAT_YSIZE(M); i++)
            VEC_ELEM(result,j) += VEC_ELEM(*this,i) * MAT_ELEM(M,i, j);

    result.setRow();
    return result;
}

#endif /* MATRIX2D_H_ */
