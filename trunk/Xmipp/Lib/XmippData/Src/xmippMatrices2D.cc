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
#include "../xmippMatrices2D.hh"
#include "../xmippArgs.hh"
#include <stdio.h>

/* ************************************************************************* */
/* IMPLEMENTATIONS                                                           */
/* ************************************************************************* */
#define maT matrix2D<T>
#define ma  matrix2D
#include "MultidimBasic.inc"
#undef ma
#undef maT
   // This file contain several general functions for arithmetic operations
   // They are defined for a general multidimensional_array<T>, but
   // we have "redirected" maT to matrix2D<T>, this way they will do OK in
   // this matrix library

/* Print shape ------------------------------------------------------------- */
template <class T>
   void mT::print_shape(ostream &out) const {
   out << "Size(Y,X): " << YSIZE(*this) << "x" << XSIZE(*this)
       << " i=[" << STARTINGY(*this) << ".." << FINISHINGY(*this) << "]"
       << " j=[" << STARTINGX(*this) << ".." << FINISHINGX(*this) << "]";
}

/* Outside ----------------------------------------------------------------- */
template <class T>
   bool mT::outside(const matrix1D<double> &v) const _THROW {
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
      const matrix1D<double> &corner2) const _THROW {
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
   if (spduptmp0>spduptmp1) return FALSE;

   spduptmp0=MAX(STARTINGX(*this), x0);
   spduptmp1=MIN(FINISHINGX(*this),x0+xdim);
   if (spduptmp0>spduptmp1) return FALSE;
   return TRUE;
}

/* Corner ------------------------------------------------------------------ */
template <class T>
   bool mT::isCorner(const matrix1D<double> &v) _THROW {
   if (XSIZE(v)<2)
      REPORT_ERROR(1,"isCorner: index vector has got not enough components");
   return ((XX(v)==STARTINGX(*this)  && YY(v)==STARTINGY(*this))  ||
           (XX(v)==STARTINGX(*this)  && YY(v)==FINISHINGY(*this)) ||
           (XX(v)==FINISHINGX(*this) && YY(v)==STARTINGY(*this))  ||
           (XX(v)==FINISHINGX(*this) && YY(v)==FINISHINGY(*this)));
}

/* Border ------------------------------------------------------------------ */
template <class T>
   bool mT::isBorder(const matrix1D<int> &v) _THROW 
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

// Special case for complex numbers
template <>
ostream& operator << (ostream& ostrm, const matrix2D<double_complex> &v) {
   if (XSIZE(v)==0 || YSIZE(v)==0)
      ostrm << "NULL matrix\n";
   else
      ostrm << endl;
      for (int i=STARTINGY(v); i<=FINISHINGY(v); i++) {
         for (int j=STARTINGX(v); j<=FINISHINGX(v); j++) {
            ostrm << MAT_ELEM(v,i,j) << ' ';
         }
         ostrm << endl;
      }
   return ostrm;
}

/* Non square matrix Identity---------------------------------------------- */
template <class T>
   void mT::init_identity(int Ydim, int Xdim) {
   if (Xdim==0 || Ydim==0) {clear(); return;}
   resize(Ydim,Xdim);
   FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
      DIRECT_MAT_ELEM(*this,i,j)=(T)(i==j);
}

/* Make a matrix from a vector --------------------------------------------- */
template <class T>
   void mT::from_vector(const vT &op1) {
   // Null vector => Null matrix
   if (op1.get_dim()==0) clear(); return;
   
   // Look at shape and copy values
   if (op1.isRow()) {
      resize(1,XSIZE(op1));
      for (int j=0; j<XSIZE(op1); j++)
         DIRECT_MAT_ELEM(*this,0,j)=DIRECT_VEC_ELEM(op1,j);
      STARTINGX(*this)=STARTINGX(op1);
      STARTINGY(*this)=0;
   } else {
      temp.resize(XSIZE(op1),1);
      for (int i=0; i<XSIZE(op1); i++)
         DIRECT_MAT_ELEM(*this,i,0)=DIRECT_VEC_ELEM(op1,i);
      STARTINGX(*this)=0;
      STARTINGY(*this)=STARTINGX(op1);
   }
   
   return temp;
}

/* Load from a numerical recipes matrix ------------------------------------ */
template <class T>
   void mT::load_from_numerical_recipes(T **m, int Ydim, int Xdim) {
   resize(Ydim,Xdim);
   for (int i=1; i<=Ydim; i++)
      for (int j=1; j<=Xdim; j++)
         (*this)(i-1,j-1)=m[i][j];
}

/* Make a vector from a matrix --------------------------------------------- */
template <class T>
   void mT::to_vector(vT &op1) const _THROW {
   // Null matrix => Null vector
   if (XSIZE(*this)==0 || YSIZE(*this)==0) {op1.clear(); return;}
   
   // If matrix is not a vector, produce an error
   if (XSIZE(*this)!=1 && (YSIZE(*this)!=1))
      REPORT_ERROR(1102, "To_vector: Matrix cannot be converted to vector");
   
   // Look at shape and copy values
   if (YSIZE(*this)==1) {
      // Row vector
      op1.resize(XSIZE(*this));
      for (int j=0; j<XSIZE(this); j++)
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

/*
 * "Bi-linear Interpolation"
 */

template <class T>
   T mT::interpolated_elem(double x, double y, T outside_value) {
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

// Special case for complex numbers
double_complex matrix2D<double_complex>::interpolated_elem(
   double x, double y, double_complex outside_value) {
    int x0 = FLOOR(x); double fx = x - x0; int x1 = x0 + 1;
    int y0 = FLOOR(y); double fy = y - y0; int y1 = y0 + 1;
    
    double_complex d00 = outside(y0,x0) ? outside_value : dMij(*this,y0,x0);
    double_complex d10 = outside(y1,x0) ? outside_value : dMij(*this,y1,x0);
    double_complex d11 = outside(y1,x1) ? outside_value : dMij(*this,y1,x1);
    double_complex d01 = outside(y0,x1) ? outside_value : dMij(*this,y0,x1);

    double_complex d0 = LIN_INTERP(fx, d00, d01);
    double_complex d1 = LIN_INTERP(fx, d10, d11);    
    return LIN_INTERP(fy, d0, d1);
}

/* Transpose --------------------------------------------------------------- */
template <class T>
   mT mT::transpose() const {
   T aux;
   mT result(XSIZE(*this),YSIZE(*this));
   FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(result)
      DIRECT_MAT_ELEM(result,i,j)=DIRECT_MAT_ELEM(*this,j,i);
   STARTINGX(result)=STARTINGX(*this);
   STARTINGY(result)=STARTINGY(*this);
   return result;
}

/* Reverse X --------------------------------------------------------------- */
template <class T>
   void mT::self_reverseX() {
   T aux;
   int jmax=(int)(XSIZE(*this)-1)/2;
   for (int i=0; i<YSIZE(*this); i++)
      for (int j=0; j<=jmax; j++) {
         SWAP(DIRECT_MAT_ELEM(*this,i,j),
            DIRECT_MAT_ELEM(*this,i,XSIZE(*this)-1-j), aux);
      }
}

/* Reverse Y --------------------------------------------------------------- */
template <class T>
   void mT::self_reverseY() {
   T aux;
   int imax=(int)(YSIZE(*this)-1)/2;
   for (int i=0; i<=imax; i++)
      for (int j=0; j<xdim; j++) {
         SWAP(DIRECT_MAT_ELEM(*this,i,j),
            DIRECT_MAT_ELEM(*this,YSIZE(*this)-1-i,j), aux);
      }
}

/* Determinant ------------------------------------------------------------- */
// (see Numerical Recipes, Chapter 2 Section 5)
template <class T> T mT::det() const _THROW {
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

/* Solve Ax=b -------------------------------------------------------------- */
// (see Numerical Recipes, Chapter 2 Section 3)
template <class T>
   void solve(const mT &A, const vT &b, vT &result) _THROW {
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
                     matrix1D<double> &result,double tolerance) _THROW 
{
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
   void solve(const mT &A, const mT &b, mT &result) _THROW {
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

/* Solve Cx=d, nonnegative x */
double solve_nonneg(const matrix2D<double> &C, const matrix1D<double> &d,
   matrix1D<double> &result) _THROW {
   if (C.xdim==0)
      REPORT_ERROR(1108, "Solve_nonneg: Matrix is empty");
   if (C.ydim!=d.get_dim())
      REPORT_ERROR(1102, "Solve_nonneg: Different sizes of Matrix and Vector");
   if (d.isRow())
      REPORT_ERROR(1107, "Solve_nonneg: Not correct vector shape");
   
   matrix2D<double> Ct(XSIZE(C),YSIZE(C)); // Ct=C^t
   FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Ct)
      DIRECT_MAT_ELEM(Ct,i,j)=DIRECT_MAT_ELEM(C,j,i);

   result.init_zeros(YSIZE(Ct));
   double rnorm;

   // Watch out that matrix Ct is transformed.
   int success=nnls(MULTIDIM_ARRAY(Ct),XSIZE(Ct),YSIZE(Ct),
        MULTIDIM_ARRAY(d),
        MULTIDIM_ARRAY(result),
        &rnorm, NULL, NULL, NULL);
   if (success==1) cerr << "Warning, too many iterations in nnls\n";
   else if (success==2)
      REPORT_ERROR(1,"Solve_nonneg: Not enough memory");
   return rnorm;
}

/* Inverse of the matrix --------------------------------------------------- */
// See Numerical Recipes Chapter 2, Section 4
template <class T>
   void mT::inv(mT &result) const {
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

//interface to numerical recipes
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

//interface to numerical recipes
template <class T>
void svdcmp(const matrix2D<T> &a,matrix2D<double> &u,
               matrix1D<double> &w, matrix2D<double> &v)
{
   // svdcmp only works with double
   type_cast(a,u);
   
   // Set size of matrices
   w.init_zeros(u.ColNo());
   v.init_zeros(u.ColNo(),u.ColNo());
   
   // Call to the numerical recipes routine
   svdcmp(MULTIDIM_ARRAY(u),
          u.RowNo(),u.ColNo(),
          MULTIDIM_ARRAY(w),
          MULTIDIM_ARRAY(v));
}

//interface to numerical recipes
void svbksb(matrix2D<double> &u,matrix1D<double> &w,matrix2D<double> &v,
             matrix1D<double> &b,matrix1D<double> &x)
{
   // Call to the numerical recipes routine. Results will be stored in X
   svbksb(u.adapt_for_numerical_recipes2(),
          w.adapt_for_numerical_recipes(),
          v.adapt_for_numerical_recipes2(),
          u.RowNo(),u.ColNo(),
          b.adapt_for_numerical_recipes(),
          x.adapt_for_numerical_recipes());
}

/* Window ------------------------------------------------------------------ */
template <class T>
   void mT::window(int y0, int x0, int yF, int xF, T init_value) {
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

/* For all rows ------------------------------------------------------------ */
template <class T>
   vT mT::for_all_rows(T (*f)(vT&)) const {
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

template <class T>
   vT mT::for_all_cols(T (*f)(vT&)) const {
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

template <class T>
   void mT::for_all_rows(vT (*f)(vT&)) {
   if (XSIZE(*this)==0 || YSIZE(*this)==0) return;
   for (int i=STARTINGY(*this); i<=FINISHINGY(*this); i++) {
      vT aux;
      getRow(i,aux);
      setRow(i,(*f)(aux));
   }
}

template <class T>
   void mT::for_all_cols(vT (*f)(vT&)) {
   if (XSIZE(*this)==0 || YSIZE(*this)==0) return;
   for (int j=STARTINGX(*this); j<=FINISHINGX(*this); j++) {
      vT aux;
      getCol(j,aux);
      setCol(j,(*f)(aux));
   }
}

/* Apply a geometrical transformation -------------------------------------- */
// It is translated from function transforma in Lib/varios2.c
// We should check which one performs better.
//#define DEBUG
template <class T>
   void apply_geom(mT &M2,matrix2D<double> A, const mT &M1, bool inv,
   bool wrap, T outside) _THROW {
   int m1, n1, m2, n2;
   double x, y, xp, yp;
   double minxp, minyp, maxxp, maxyp;
   int   cen_x, cen_y, cen_xp, cen_yp;
   double wx, wy; // Weights in X,Y directions for bilinear interpolation
   
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
   
   // Now we go from the output image to the input image, ie, for any pixel
   // in the output image we calculate which are the corresponding ones in
   // the original image, make an interpolation with them and put this value
   // at the output pixel
   //#define DEBUG
   #ifdef DEBUG
      cout << "A\n"    << A     << endl
	   << "(cen_x ,cen_y )=(" << cen_x  << "," << cen_y  << ")\n"
	   << "(cen_xp,cen_yp)=(" << cen_xp << "," << cen_yp << ")\n"
	   << "(min_xp,min_yp)=(" << minxp  << "," << minyp  << ")\n"
	   << "(max_xp,max_yp)=(" << maxxp  << "," << maxyp  << ")\n"
      ;
   #endif
      //#undef DEBUG
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
         
         #ifdef DEBUG
         cout << "Computing (" << i << "," << j << ")\n";
         cout << "   (y, x) =(" << y << "," << x << ")\n"
              << "   before wrapping (y',x')=(" << yp << "," << xp << ") " << endl;
         #endif

         // If the point is outside the image, apply a periodic extension
         // of the image, what exits by one side enters by the other
         interp=TRUE;
         if (wrap) {
            if (xp<minxp-XMIPP_EQUAL_ACCURACY ||
	        xp>maxxp+XMIPP_EQUAL_ACCURACY)
		xp=realWRAP(xp, minxp, maxxp);
            if (yp<minyp-XMIPP_EQUAL_ACCURACY ||
	        yp>maxyp+XMIPP_EQUAL_ACCURACY) yp=realWRAP(yp, minyp, maxyp);
         } else {
            if (xp<minxp-XMIPP_EQUAL_ACCURACY ||
	        xp>maxxp+XMIPP_EQUAL_ACCURACY) interp=FALSE;
            if (yp<minyp-XMIPP_EQUAL_ACCURACY ||
	        yp>maxyp+XMIPP_EQUAL_ACCURACY) interp=FALSE;
         }

         #ifdef DEBUG
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

            #ifdef DEBUG
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

            #ifdef DEBUG
            cout << "   val= " << tmp << endl;
            #endif
         }

         // Compute new point inside input image
         xp += dMij(A,0,0);
         yp += dMij(A,1,0);
      }
   }
}
#undef DEBUG

/* Geometrical operations -------------------------------------------------- */
template <class T>
   void mT::rotate(double ang, mT &result, bool wrap) const
      {matrix2D<double> temp=rot2D_matrix(ang);
       apply_geom(result,temp,*this,IS_NOT_INV,wrap);}

template <class T>
   void mT::translate(const matrix1D<double> &v, mT &result, bool wrap) const
      {matrix2D<double> temp=translation2D_matrix(v);
       apply_geom(result,temp,*this,IS_NOT_INV,wrap);}

template <class T>
   void mT::scale_to_size(int Ydim,int Xdim, mT &result) const {
      matrix2D<double> temp(3,3);
      result.resize(Ydim,Xdim);
      temp.init_identity();
      DIRECT_MAT_ELEM(temp,0,0)=(double)Xdim/(double)XSIZE(*this);
      DIRECT_MAT_ELEM(temp,1,1)=(double)Ydim/(double)YSIZE(*this);
      apply_geom(result,temp,*this,IS_NOT_INV,WRAP);
   }

template <class T>
   void mT::self_translate_center_of_mass_to_center(bool wrap) {
      set_Xmipp_origin();
      matrix1D<double> center;
      center_of_mass(center);
      center*=-1;
      self_translate(center,wrap);
   }

template <class T>
   void mT::superpixel_reduce(mT &result, int size) const {
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

template <class T>
   void mT::superpixel_expand(mT &result, int size) const {
      result.init_zeros(YSIZE(*this)*size,XSIZE(*this)*size);
      FOR_ALL_ELEMENTS_IN_MATRIX2D(*this) {
         for (int ii=0; ii<size; ii++)
	    for (int jj=0; jj<size; jj++)
	       DIRECT_MAT_ELEM(result,size*i+ii,size*j+jj)=
	          DIRECT_MAT_ELEM(*this,i,j);
      }
   }

/* Matrix by matrix multiplication ----------------------------------------- */
// xinit and yinit are taken from op1
template <class T>
   void mul_matrix(const mT &op1, const mT &op2, mT &result) _THROW {
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

/* Matrix by vector multiplication ----------------------------------------- */
template <class T>
   vT mT::operator * (const vT &op1) const {
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

/* Is diagonal ------------------------------------------------------------- */
template <class T>
   bool mT::IsDiagonal() const
      {if (XSIZE(*this)!=YSIZE(*this)) return FALSE;
       FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
             if (i!=j && ABS(DIRECT_MAT_ELEM(*this,i,j))>XMIPP_EQUAL_ACCURACY)
                return FALSE;
       return TRUE;}

// Special case for complex matrices
bool matrix2D<double_complex>::IsDiagonal() const
   {if (XSIZE(*this)!=YSIZE(*this)) return FALSE;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
          if (i!=j && abs(DIRECT_MAT_ELEM(*this,i,j))>XMIPP_EQUAL_ACCURACY)
             return FALSE;
    return TRUE;}

/* Is Scalar --------------------------------------------------------------- */
template <class T>
   bool mT::IsScalar() const
      {if (!IsDiagonal()) return FALSE;
       for (int i=1; i<YSIZE(*this); i++)
          if (ABS(DIRECT_MAT_ELEM(*this,i,i)-DIRECT_MAT_ELEM(*this,0,0))>
              XMIPP_EQUAL_ACCURACY) return FALSE;
       return TRUE;}

// Special case for complex matrices
bool matrix2D<double_complex>::IsScalar() const
   {if (!IsDiagonal()) return FALSE;
    for (int i=1; i<YSIZE(*this); i++)
       if (abs(DIRECT_MAT_ELEM(*this,i,i)-DIRECT_MAT_ELEM(*this,0,0))>
           XMIPP_EQUAL_ACCURACY) return FALSE;
    return TRUE;}

/* Is Symmetric ------------------------------------------------------------ */
template <class T>
   bool mT::IsSymmetric() const
      {if (XSIZE(*this)!=YSIZE(*this)) return FALSE;
       for (int i=0; i<YSIZE(*this); i++)
          for (int j=i+1; j<XSIZE(*this); j++)
             if (ABS(DIRECT_MAT_ELEM(*this,i,j)-DIRECT_MAT_ELEM(*this,j,i))>
                 XMIPP_EQUAL_ACCURACY) return FALSE;
       return TRUE;}

// Special case for complex matrices
bool matrix2D<double_complex>::IsSymmetric() const
   {if (XSIZE(*this)!=YSIZE(*this)) return FALSE;
    for (int i=0; i<YSIZE(*this); i++)
       for (int j=i+1; j<XSIZE(*this); j++)
          if (abs(DIRECT_MAT_ELEM(*this,i,j)-DIRECT_MAT_ELEM(*this,j,i))>
              XMIPP_EQUAL_ACCURACY) return FALSE;
    return TRUE;}

/* Is Skew symmetric ------------------------------------------------------- */
template <class T>
   bool mT::IsSkewSymmetric() const
      {if (XSIZE(*this)!=YSIZE(*this)) return FALSE;
       for (int i=0; i<YSIZE(*this); i++)
          for (int j=i+1; j<XSIZE(*this); j++)
             if (ABS(DIRECT_MAT_ELEM(*this,i,j)+DIRECT_MAT_ELEM(*this,j,i))>
                 XMIPP_EQUAL_ACCURACY) return FALSE;
       return TRUE;}


// Special case for complex matrices
bool matrix2D<double_complex>::IsSkewSymmetric() const
   {if (XSIZE(*this)!=YSIZE(*this)) return FALSE;
    for (int i=0; i<YSIZE(*this); i++)
       for (int j=i+1; j<XSIZE(*this); j++)
          if (abs(DIRECT_MAT_ELEM(*this,i,j)+DIRECT_MAT_ELEM(*this,j,i))>
              XMIPP_EQUAL_ACCURACY) return FALSE;
    return TRUE;}

/* Is Upper triangular ----------------------------------------------------- */
template <class T>
   bool mT::IsUpperTriangular() const
      {if (XSIZE(*this)!=YSIZE(*this)) return FALSE;
       for (int i=1; i<YSIZE(*this); i++)
          for (int j=0; j<i-1; j++)
             if (ABS(DIRECT_MAT_ELEM(*this,i,j))>XMIPP_EQUAL_ACCURACY)
                return FALSE;
       return TRUE;}

// Special case for complex matrices
bool matrix2D<double_complex>::IsUpperTriangular() const
   {if (XSIZE(*this)!=YSIZE(*this)) return FALSE;
    for (int i=1; i<YSIZE(*this); i++)
       for (int j=0; j<i-1; j++)
          if (abs(DIRECT_MAT_ELEM(*this,i,j))>XMIPP_EQUAL_ACCURACY)
             return FALSE;
    return TRUE;}

/* Is Lower triangular ----------------------------------------------------- */
template <class T>
   bool mT::IsLowerTriangular() const
      {if (XSIZE(*this)!=YSIZE(*this)) return FALSE;
       for (int i=1; i<YSIZE(*this); i++)
          for (int j=i+1; j<XSIZE(*this); j++)
             if (ABS(DIRECT_MAT_ELEM(*this,i,j))>XMIPP_EQUAL_ACCURACY)
                return FALSE;
       return TRUE;}

// Special case for complex matrices
bool matrix2D<double_complex>::IsLowerTriangular() const
   {if (XSIZE(*this)!=YSIZE(*this)) return FALSE;
    for (int i=1; i<YSIZE(*this); i++)
       for (int j=i+1; j<XSIZE(*this); j++)
          if (abs(DIRECT_MAT_ELEM(*this,i,j))>XMIPP_EQUAL_ACCURACY)
             return FALSE;
    return TRUE;}

/* Max index --------------------------------------------------------------- */
template <class T>
   void mT::max_index(int &imax, int &jmax) const {
   if (XSIZE(*this)==0) {imax=jmax=-1; return;}
   imax=jmax=0;
   T   max=MAT_ELEM(*this,imax,jmax);
   FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
      if (MAT_ELEM(*this,i,j)>max)
         {max=MAT_ELEM(*this,imax,jmax); imax=i; jmax=j;}
}

/* Min index --------------------------------------------------------------- */
template <class T>
   void mT::min_index(int &imin, int &jmin) const {
   if (XSIZE(*this)==0) {imin=jmin=-1; return;}
   imin=jmin=0;
   T   min=MAT_ELEM(*this,imin,jmin);
   FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
      if (MAT_ELEM(*this,i,j)>min)
         {min=MAT_ELEM(*this,imin,jmin); imin=i; jmin=j;}
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

/* Rotation 2D ------------------------------------------------------------- */
matrix2D<double> rot2D_matrix(double ang) {
   matrix2D<double> result(3,3);
   double cosine, sine;
   
   ang=DEG2RAD(ang);
   cosine=cos(ang);
   sine=sin(ang);

   DIRECT_MAT_ELEM(result,0,0)=cosine;
   DIRECT_MAT_ELEM(result,0,1)=-sine;
   DIRECT_MAT_ELEM(result,0,2)=0;
 
   DIRECT_MAT_ELEM(result,1,0)=sine;
   DIRECT_MAT_ELEM(result,1,1)=cosine;
   DIRECT_MAT_ELEM(result,1,2)=0;
 
   DIRECT_MAT_ELEM(result,2,0)=0;
   DIRECT_MAT_ELEM(result,2,1)=0;
   DIRECT_MAT_ELEM(result,2,2)=1;
   
   return result;
}

/* Translation 2D ---------------------------------------------------------- */
matrix2D<double> translation2D_matrix(matrix1D<double> v) _THROW {
   if (v.get_dim()!=2)
      REPORT_ERROR(1002,"Translation2D_matrix: vector is not in R2");

   matrix2D<double> result(3,3);

   result.init_identity();
   DIRECT_MAT_ELEM(result,0,2)=XX(v);
   DIRECT_MAT_ELEM(result,1,2)=YY(v);

   return result;   
}

/* Rotation 3D around the system axes -------------------------------------- */
matrix2D<double> rot3D_matrix(double ang, char axis) _THROW {
   matrix2D<double> result(4,4);
   double cosine, sine;
   
   ang=DEG2RAD(ang);
   cosine=cos(ang);
   sine=sin(ang);
   
   result.init_zeros();
   DIRECT_MAT_ELEM(result,3,3)=1;
   switch (axis) {
      case 'Z':
         DIRECT_MAT_ELEM(result,0,0)=cosine;
         DIRECT_MAT_ELEM(result,0,1)=-sine;
         DIRECT_MAT_ELEM(result,1,0)=sine;
         DIRECT_MAT_ELEM(result,1,1)=cosine;
         DIRECT_MAT_ELEM(result,2,2)=1;
         break;
      case 'Y':
         DIRECT_MAT_ELEM(result,0,0)=cosine;
         DIRECT_MAT_ELEM(result,0,2)=-sine;
         DIRECT_MAT_ELEM(result,2,0)=sine;
         DIRECT_MAT_ELEM(result,2,2)=cosine;
         DIRECT_MAT_ELEM(result,1,1)=1;
         break;
      case 'X':
         DIRECT_MAT_ELEM(result,1,1)=cosine;
         DIRECT_MAT_ELEM(result,1,2)=-sine;
         DIRECT_MAT_ELEM(result,2,1)=sine;
         DIRECT_MAT_ELEM(result,2,2)=cosine;
         DIRECT_MAT_ELEM(result,0,0)=1;
         break;
      default:
         REPORT_ERROR(1105,"rot3D_matrix: Unknown axis");
   }
   return result;
}

/* Align a vector with Z axis */
matrix2D<double> align_with_Z(const matrix1D<double> &axis) _THROW {
   matrix1D<double>  Axis;
   matrix2D<double>  A(4,4);

   if (axis.get_dim()!=3)
      REPORT_ERROR(1002,"align_with_Z: Axis is not in R3");
      
   // Copy axis and compute length of the projection on YZ plane
   Axis=axis.normalize();
   double proj_mod=sqrt(YY(Axis)*YY(Axis)+ZZ(Axis)*ZZ(Axis));

   A(3,3)=1;
   if (proj_mod>XMIPP_EQUAL_ACCURACY) { // proj_mod!=0
      // Build Matrix A, which makes the turning axis coincident with Z
      A(0,0)= proj_mod;
      A(0,1)=-XX(Axis)*YY(Axis)/proj_mod;
      A(0,2)=-XX(Axis)*ZZ(Axis)/proj_mod;
      A(1,0)= 0;
      A(1,1)= ZZ(Axis)/proj_mod;
      A(1,2)=-YY(Axis)/proj_mod;
      A(2,0)= XX(Axis);
      A(2,1)= YY(Axis);
      A(2,2)= ZZ(Axis);
   } else {
      // I know that the Axis is the X axis
      A(0,0)= 0; A(0,1)=0; A(0,2)=-1;
      A(1,0)= 0; A(1,1)=1; A(1,2)=0;
      A(2,0)= 1; A(2,1)=0; A(2,2)=0;
   }
   return A;
}

/* Rotation 3D around any axis -------------------------------------------- */
matrix2D<double> rot3D_matrix(double ang, const matrix1D<double> &axis) {
   // Compute a matrix which makes the turning axis coincident with Z
   // And turn around this axis
   matrix2D<double> A=align_with_Z(axis);
   return A.transpose() * rot3D_matrix(ang,'Z') * A;
}

/* Translation 3D ---------------------------------------------------------- */
matrix2D<double> translation3D_matrix(const matrix1D<double> &v) _THROW {
   if (XSIZE(v)!=3)
      REPORT_ERROR(1002,"Translation3D_matrix: vector is not in R3");

   matrix2D<double> result(4,4);

   result.init_identity();
   DIRECT_MAT_ELEM(result,0,3)=XX(v);
   DIRECT_MAT_ELEM(result,1,3)=YY(v);
   DIRECT_MAT_ELEM(result,2,3)=ZZ(v);

   return result;   
}

/* Scale 3D ---------------------------------------------------------------- */
matrix2D<double> scale3D_matrix(const matrix1D<double> &sc) _THROW {
   if (XSIZE(sc)!=3)
      REPORT_ERROR(1002,"Scale3D_matrix: vector is not in R3");

   matrix2D<double> result(4,4);

   result.init_identity();
   DIRECT_MAT_ELEM(result,0,0)=XX(sc);
   DIRECT_MAT_ELEM(result,1,1)=YY(sc);
   DIRECT_MAT_ELEM(result,2,2)=ZZ(sc);

   return result;   
}

/* Cut to common size ------------------------------------------------------ */
template <class T>
   void cut_to_common_size(mT &V1, mT &V2) {
   int y0=MAX(STARTINGY(V1) ,STARTINGY(V2));
   int yF=MIN(FINISHINGY(V1),FINISHINGY(V2));
   int x0=MAX(STARTINGX(V1) ,STARTINGX(V2));
   int xF=MIN(FINISHINGX(V1),FINISHINGX(V2));
   V1.window(y0,x0,yF,xF);
   V2.window(y0,x0,yF,xF);
}

/* Radial average ---------------------------------------------------------- */
template <class T>
void radial_average(const matrix2D<T> &m, const matrix1D<int> &center_of_rot,
                    matrix1D<T> &radial_mean) _THROW
{
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
   
   // Define the vectors
   radial_mean.resize(dim);
   matrix1D<int> radial_count(dim);
   
   
   /* Perform the radial sum and count pixels that contribute to
      every distance */
   FOR_ALL_ELEMENTS_IN_MATRIX2D(m)
   {
      YY(idx)=i-YY(center_of_rot);
	  XX(idx)=j-XX(center_of_rot);
	  // Determine distance to the center
	  int distance=(int)floor(idx.module());
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

/* Quadratic form ---------------------------------------------------------- */
void eval_quadratic(const matrix1D<double> &x, const matrix1D<double> &c,
   const matrix2D<double> &H, double &val, matrix1D<double> &grad) _THROW {
   if (XSIZE(x)!=XSIZE(c))
      REPORT_ERROR(1102,"Eval_quadratic: Not compatible sizes in x and c");
   if (XSIZE(H)!=XSIZE(x))
      REPORT_ERROR(1102,"Eval_quadratic: Not compatible sizes in x and H");

   // H*x, store in grad
   grad.init_zeros(XSIZE(x));
   for (int i=0; i<YSIZE(H); i++)
      for (int j=0; j<XSIZE(x); j++)
            DIRECT_VEC_ELEM(grad,i) += DIRECT_MAT_ELEM(H,i,j)*
               DIRECT_VEC_ELEM(x,j);
   
   // Now, compute c^t*x+1/2*x^t*H*x
   // Add c to the gradient
   double quad=0;
   val=0;
   for (int j=0; j<XSIZE(x); j++) {
      quad+=DIRECT_VEC_ELEM(grad,j)*DIRECT_VEC_ELEM(grad,j); // quad=x^t*H^t*H*x
      val+=DIRECT_VEC_ELEM(c,j)*DIRECT_VEC_ELEM(x,j);        // val=c^t*x

      DIRECT_VEC_ELEM(grad,j)+=DIRECT_VEC_ELEM(c,j);         // grad+=c
   }
   val+=0.5*quad;
}

/* ------------------------------------------------------------------------- */
/* INSTANTIATE                                                               */
/* ------------------------------------------------------------------------- */
template <class T>
   T vector_first(vT &x) {return x(0);}

template <class T>
   vT vector_minus_first(vT &x) {return x-x(0);}

template <class T>
   void instantiate_matrix (matrix2D<T> v) {
      matrix2D<T>      a;
      matrix1D<double> r;
      matrix1D<int>    ir;
      double           d;
      T                Taux;
      
      // General functions for multidimensional arrays
      a==a;
      a=1-a;
      a=a-1;
      a=a*a;
      a=a+a;
      a.print_stats();
      a.compute_max();
      a.compute_min();
      a.compute_avg();
      a.compute_stddev();
      a.compute_double_minmax(d,d);
      a.compute_double_minmax(d,d,r,r);
      a.compute_double_minmax(d,d,ir,ir);
      a.compute_stats(d,d,Taux,Taux,ir,ir);
      a.range_adjust(0,1);
      a.statistics_adjust(0,1);
      a.effective_range(99);
      a.init_random(0,1);
      a.add_noise(0,1);      
      a.threshold("abs_above",1,0);
      a.count_threshold("abs_above",1,0);
      a.center_of_mass(r);

      a.print_shape();
      a.outside(r);
      a.outside(0,0);
	  
 	  matrix1D<int> pixel(2);
	  a.isBorder(pixel);
      a.isBorder(0,0);

      a.intersects(a);
      a.intersects(r,r);
      a.isCorner(r);
      a.patch(a);
      cout << a;
      a.window(0,0,1,1);
      cut_to_common_size(a,a);
   
      // Specific for matrices
      matrix1D<T> b;
      matrix2D<T> A;
      a.init_identity(10,10);
      a.interpolated_elem(3.5,3.5);
      a.reverseX();
      a.reverseY();
      a=a.inv()*a.det();
      b=a*b; b=b*a; /****b=b*a;*/
      solve(A,b,b);
      solve(A,A,A);
      matrix2D<double> B;
      apply_geom(a, B, a, IS_NOT_INV,DONT_WRAP);
      a.scale_to_size(32,32,a);
      a.superpixel_reduce(a,2);
      a.superpixel_expand(a,2);
      a.for_all_rows(&vector_minus_first);
      a.for_all_rows(&vector_first);
      a.for_all_cols(&vector_minus_first);
      a.for_all_cols(&vector_first);
      if (a.IsDiagonal() || a.IsScalar() || a.IsSymmetric() ||
          a.IsSkewSymmetric() || a.IsUpperTriangular() ||
          a.IsLowerTriangular()) cout << "Special matrix\n";
      a.translate(vector_R2(0,0));
      a.rotate(60);
      a.scale_to_size(45,45);
      a.self_translate_center_of_mass_to_center();
      int imax;
      a.max_index(imax,imax);
      a.min_index(imax,imax);
	  
	  matrix1D<int> center;
	  radial_average(v,center,b);
      // Singular value decomposition
      svdcmp(v,B,r,B);
      solve_by_svd(v,b,r,d);
}

void instantiate_complex_matrix () {
      matrix2D<double_complex> a;
      matrix1D<double>         r;
      
      // General functions for multidimensional arrays
      a.print_shape();
      a.outside(r);
      a.outside(0,0);
      
      a=1.0-a;
      a=a-1.0;
      a=a*a;
      a.intersects(a);
      a.intersects(r,r);
      a.isCorner(r);
      a.patch(a);
      cout << a;
      a.window(0,0,1,1);
      cut_to_common_size(a,a);
   
      // Specific for matrices
      matrix1D<double_complex> b;
      a.init_identity(10,10);
      a.interpolated_elem(3.5,3.5);
      a.reverseX();
      a.reverseY();
      b=a*b; /**** b=b*a; */
      matrix2D<double> B;
      apply_geom(a, B, a, IS_NOT_INV,DONT_WRAP);
      a.scale_to_size(32,32,a);
      a.superpixel_reduce(a,2);
      a.superpixel_expand(a,2);
      a.for_all_rows(&vector_minus_first);
      a.for_all_rows(&vector_first);
      a.for_all_cols(&vector_minus_first);
      a.for_all_cols(&vector_first);
      if (a.IsDiagonal() || a.IsScalar() || a.IsSymmetric() ||
          a.IsSkewSymmetric() || a.IsUpperTriangular() ||
          a.IsLowerTriangular()) cout << "Special matrix\n";
      a.translate(vector_R2(0,0));
      a.rotate(60);
      a.scale_to_size(45,45);
	  	  
}

void instantiate2D() {
   matrix2D<char>           V0; instantiate_matrix(V0);
   matrix2D<short>          V1; instantiate_matrix(V1);
   matrix2D<int>            V2; instantiate_matrix(V2);
   matrix2D<float>          V3; instantiate_matrix(V3);
   matrix2D<double>         V4; instantiate_matrix(V4);
}
