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
/* VECTORS                                                                   */
/* ------------------------------------------------------------------------- */
#include "../xmippMatrices1D.hh"
#include "../xmippArgs.hh"
#include <strstream.h>
#include <stdio.h>

/* ************************************************************************* */
/* IMPLEMENTATIONS                                                           */
/* ************************************************************************* */
#define maT matrix1D<T>
#define ma  matrix1D
#include "MultidimBasic.inc"
#undef ma
#undef maT
   // This file contain several general functions for arithmetic operations
   // They are defined for a general multidimensional_array<T>, but
   // we have "redirected" maT to matrix1D<T>, this way they will do OK in
   // this vector library

/* Print shape ------------------------------------------------------------- */
template <class T>
   void vT::print_shape(ostream &out) const {
   out << "Size: " << XSIZE(*this)
       << "i=[" << STARTINGX(*this) << ".." << FINISHINGX(*this) << "]";
}

/* Outside ----------------------------------------------------------------- */
template <class T>
   bool vT::outside(const matrix1D<double> &v) const _THROW {
   if (XSIZE(v)<1)
      REPORT_ERROR(1,"Outside: index vector has got not enough components");
   return (XX(v)<STARTINGX(*this) || XX(v)>FINISHINGX(*this));
}

template <class T>
   bool vT::outside(int i) const {
   return (i<STARTINGX(*this) || i>FINISHINGX(*this));
}

/* Intersects -------------------------------------------------------------- */
template <class T>
   bool vT::intersects(const vT &m) const
      {return intersects(STARTINGX(m), XSIZE(m)-1);}

template <class T>
   bool vT::intersects(const matrix1D<double> &corner1,
      const matrix1D<double> &corner2) const _THROW {
       if (XSIZE(corner1)!=1 || XSIZE(corner2)!=1)
          REPORT_ERROR(1002,"intersects 1D: corner sizes are not 1");
       return intersects(XX(corner1),XX(corner2)-XX(corner1));
}

template <class T>
   bool vT::intersects(double x0, double xdim) const {
   SPEED_UP_temps;
   spduptmp0=MAX(STARTINGX(*this), x0);
   spduptmp1=MIN(FINISHINGX(*this),x0+xdim);
   if (spduptmp0>spduptmp1) return FALSE;
   return TRUE;
}

/* Corner ------------------------------------------------------------------ */
template <class T>
   bool vT::isCorner(const matrix1D<double> &v) _THROW {
   if (XSIZE(v)<1)
      REPORT_ERROR(1,"isCorner: index vector has got not enough components");
   return (XX(v)==STARTINGX(*this) || XX(v)==FINISHINGX(*this));
}

/* Border ------------------------------------------------------------------ */
template <class T>
   bool vT::isBorder(const matrix1D<int> &v) _THROW 
{
   if (XSIZE(v)<1)
      REPORT_ERROR(1,"isBorder: index vector has got not enough components");
   return  isBorder(XX(v));
}

template <class T>
   bool vT::isBorder(int i) 
{
   return (i==STARTINGX(*this)  || i==FINISHINGX(*this));
}

/* Patch ------------------------------------------------------------------- */
template <class T>
   void vT::patch(const vT &patch_array, char operation) {
      SPEED_UP_temps;
      FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX1D(patch_array,*this)
         switch (operation) {
            case '=': VEC_ELEM(*this,i) =VEC_ELEM(patch_array,i); break;
            case '+': VEC_ELEM(*this,i)+=VEC_ELEM(patch_array,i); break;
            case '-': VEC_ELEM(*this,i)-=VEC_ELEM(patch_array,i); break;
            case '*': VEC_ELEM(*this,i)*=VEC_ELEM(patch_array,i); break;
            case '/': VEC_ELEM(*this,i)/=VEC_ELEM(patch_array,i); break;
         }
}

/* Output stream ----------------------------------------------------------- */
template <class T>
   ostream& operator << (ostream& out, const vT& v) {
   if (MULTIDIM_SIZE(v)==0) out << "NULL vector\n";
   else {
      // Look for the exponent
      vT aux(v); aux.ABSnD();
      int prec=best_prec(aux.compute_max(),10);

      FOR_ALL_ELEMENTS_IN_MATRIX1D(v) {      
         if (v.row) out << FtoA((double)VEC_ELEM(v,i),10,prec) << ' ';
         else       out << FtoA((double)VEC_ELEM(v,i),10,prec) << '\n';
      }
   }
      
   return out;
}

// Special case for complex numbers
template <>
ostream& operator << (ostream& out, const matrix1D<double_complex>& v) {
   if (MULTIDIM_SIZE(v)==0) out << "NULL vector\n";
   else {
      FOR_ALL_ELEMENTS_IN_MATRIX1D(v) {      
         if (v.row) out << VEC_ELEM(v,i) << ' ';
         else       out << VEC_ELEM(v,i) << '\n';
      }
   }
   return out;
}

/* Linear initialisation --------------------------------------------------- */
// It is not an error that the vector is empty
template <class T>
   void vT::init_linear(T minF, T maxF, int n, const string &mode) _THROW {
   double slope;
   int steps;

   if (mode=="incr") {
      steps=1+(int) FLOOR((double)ABS((maxF-minF))/((double) n));
      slope=n*SGN(maxF-minF);
   } else if (mode=="steps") {
      steps=n;
      slope=(maxF-minF)/(steps-1);
   } else 
      REPORT_ERROR(1005,"Init_linear: Mode not supported ("+mode+")");
   
   if (steps==0) clear();
   else {   
      resize(steps);
      for (int i=0; i<steps; i++)
         VEC_ELEM(*this,i)=(T) ((double)minF+slope*i);
   }
}      

/* Reverse ----------------------------------------------------------------- */
template <class T>
   void vT::self_reverse() {
   int imax=(int)(XSIZE(*this)-1)/2;
   for (int i=0; i<=imax; i++) {
      T aux;
      SWAP(MULTIDIM_ELEM(*this,i),MULTIDIM_ELEM(*this,XSIZE(*this)-1-i),aux);
   }
}

/* Sort vector ------------------------------------------------------------- */
template <class T>
   vT vT::sort() const {
   vT temp;
   matrix1D<double> aux;
   
   if (XSIZE(*this)==0) return temp;
   
   // Initialise data
   type_cast(*this, aux);
   
   // Sort
   double * aux_array=aux.adapt_for_numerical_recipes();
   qcksrt(xdim,aux_array);

   type_cast(aux,temp);
   return temp;
}

/* Index in order ---------------------------------------------------------- */
template <class T>
   matrix1D<int> vT::index_sort() const {
   matrix1D<int>   indx;
   matrix1D<double> temp;
   
   if (XSIZE(*this)==0) return indx;
   if (XSIZE(*this)==1) {
      indx.resize(1); indx(0)=1;
      return indx;
   }
   
   // Initialise data
   indx.resize(xdim);
   type_cast(*this,temp);

   // Sort indexes
   double * temp_array=temp.adapt_for_numerical_recipes();
   int   * indx_array=indx.adapt_for_numerical_recipes();
   indexx(xdim,temp_array,indx_array);

   return indx;
}

/* Window ------------------------------------------------------------------ */
template <class T>
   void vT::window(int x0, int xF, T init_value) {
   vT result(xF-x0+1);
   STARTINGX(result)=x0;
      
   for (int j=x0; j<=xF; j++)
       if (j>=STARTINGX(*this) && j<=FINISHINGX(*this))
               VEC_ELEM(result,j)=VEC_ELEM(*this,j);
       else
               VEC_ELEM(result,j)=init_value;
	        
   *this=result;
}

/* Max index --------------------------------------------------------------- */
template <class T>
   void vT::max_index(int &imax) const {
   if (XSIZE(*this)==0) {imax=-1; return;}
   imax=0;
   T   max=VEC_ELEM(*this,imax);
   FOR_ALL_ELEMENTS_IN_MATRIX1D(*this)
      if (VEC_ELEM(*this,i)>max) {max=VEC_ELEM(*this,i); imax=i;}
}

/* Min index --------------------------------------------------------------- */
template <class T>
   void vT::min_index(int &imin) const {
   if (XSIZE(*this)==0) {imin=-1; return;}
   imin=0;
   T   min=VEC_ELEM(*this,imin);
   FOR_ALL_ELEMENTS_IN_MATRIX1D(*this)
      if (VEC_ELEM(*this,i)<min) {min=VEC_ELEM(*this,i); imin=i;}
}

/* Statistics in region ---------------------------------------------------- */
template <class T>
void vT::compute_stats(double &avg, double &stddev, T &min_val, T &max_val,
   const matrix1D<double> &corner1, const matrix1D<double> &corner2) const {
   min_val=max_val=(*this)(corner1);
   matrix1D<double> r(3);
   double N=0, sum=0, sum2=0;
   FOR_ALL_ELEMENTS_IN_MATRIX1D_BETWEEN(corner1,corner2) {
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
void vT::compute_double_minmax(double &min_val, double &max_val,
   const matrix1D<double> &corner1, const matrix1D<double> &corner2) const {
   min_val=max_val=(*this)(corner1);
   matrix1D<double> r(1);
   FOR_ALL_ELEMENTS_IN_MATRIX1D_BETWEEN(corner1,corner2) {
      if      ((*this)(r)<min_val) min_val=(*this)(r);
      else if ((*this)(r)>max_val) max_val=(*this)(r);
   }
}

/* Vector R2 and R3 -------------------------------------------------------- */
matrix1D<double> vector_R2(double x, double y) {
   matrix1D<double> result(2);
   VEC_ELEM(result,0)=x;
   VEC_ELEM(result,1)=y;
   return result;
}

matrix1D<double> vector_R3(double x, double y, double z) {
   matrix1D<double> result(3);
   VEC_ELEM(result,0)=x;
   VEC_ELEM(result,1)=y;
   VEC_ELEM(result,2)=z;
   return result;
}

matrix1D<int> vector_R3(int x, int y, int z) {
   matrix1D<int> result(3);
   VEC_ELEM(result,0)=x;
   VEC_ELEM(result,1)=y;
   VEC_ELEM(result,2)=z;
   return result;
}

/* Are orthogonal ---------------------------------------------------------- */
int are_orthogonal (matrix1D<double> &v1, matrix1D<double> &v2,
   matrix1D<double> &v3) _THROW {
   if (XSIZE(v1)!=3 || XSIZE(v2)!=3 || XSIZE(v3)!=3)
      REPORT_ERROR(1002,"Are orthogonal: Some vector do not belong to R3");
   try{
      if (dot_product(v1,v2)!=0) return 0;
      if (dot_product(v2,v3)!=0) return 0;
      if (dot_product(v1,v3)!=0) return 0;
   } catch (Xmipp_error) {
      REPORT_ERROR(1007,"Are orthogonal: Vectors are not all of the same shape");
   }
   return 1;
}

/* Are system? ------------------------------------------------------------- */
int are_system (matrix1D<double> &v1, matrix1D<double> &v2,
   matrix1D<double> &v3) _THROW {
   matrix1D<double> aux(3);
   if (XSIZE(v1)!=3 || XSIZE(v2)!=3 || XSIZE(v3)!=3)
     REPORT_ERROR(1002,"Are orthogonal: Some vector do not belong to R3");
   aux=vector_product(v1,v2); if (aux!=v3) return 0;
   aux=vector_product(v2,v3); if (aux!=v1) return 0;
   aux=vector_product(v3,v1); if (aux!=v2) return 0;
   return 1;
}

/* Cut to common size ------------------------------------------------------ */
template <class T>
   void cut_to_common_size(vT &V1, vT &V2) {
   int x0=MAX(STARTINGX(V1) ,STARTINGX(V2));
   int xF=MIN(FINISHINGX(V1),FINISHINGX(V2));
   V1.window(x0,xF);
   V2.window(x0,xF);
}

/* Sort two vectors -------------------------------------------------------- */
template <class T>
   void sort_two_vectors(vT &v1, vT &v2) _THROW {
   T temp;
   if (XSIZE(v1)!=XSIZE(v2) || STARTINGX(v1)!=STARTINGX(v2))
      REPORT_ERROR(1007, "sort_two_vectors: vectors are not of the same shape");
   FOR_ALL_ELEMENTS_IN_MATRIX1D(v1) {
      temp=MIN(VEC_ELEM(v1,i),VEC_ELEM(v2,i));
      VEC_ELEM(v2,i)=MAX(VEC_ELEM(v1,i),VEC_ELEM(v2,i));
      VEC_ELEM(v1,i)=temp;
   }
}

/* Powell's optimizer ------------------------------------------------------ */
void Powell_optimizer(matrix1D<double> &p, int i0, int n,
   double (*f)(double *x), double ftol, double &fret,
   int &iter, const matrix1D<double> &steps, bool show) {
   double **xi=NULL;

   // Adapt indexes of p
   double *pptr=p.adapt_for_numerical_recipes();
   double *auxpptr=pptr+(i0-1);

   // Form direction matrix
   ask_Tmatrix(xi,1,n,1,n);
   for (int i=1; i<=n; i++)
       for (int j=1; j<=n; j++)
           xi[i][j]=(i==j)?steps(i-1):0;
   
   // Optimize
   powell(auxpptr,xi,n,ftol,iter,fret,f, show);
   
   // Exit
   free_Tmatrix(xi,1,n,1,n);
   p.kill_adaptation_for_numerical_recipes(pptr);
}

/* Center of mass ---------------------------------------------------------- */
template <class T>
   void vT::center_of_mass(matrix1D<double> &center, void * mask) {
      center.init_zeros(1);
      double mass=0;
      matrix1D<int> *imask=(matrix1D<int> *) mask;
      FOR_ALL_ELEMENTS_IN_MATRIX1D(*this) {
         if (imask==NULL || VEC_ELEM(*imask,i)) {
            XX(center)+=i*VEC_ELEM(*this,i);
	    mass+=VEC_ELEM(*this,i);
      	 }
      }
      if (mass!=0) center/=mass;
   }

/*
   The multidimensional arrays must be instantiated if we want them to compile
   into a library. All methods of all possible classes must be called. The
   only basic types considered are:
   
      short
      int
      double
      double
      double_complex

   For double_complex arrays not all methods are valid, so they are not in
   the library.

*/

/* ------------------------------------------------------------------------- */
/* INSTANTIATE                                                                   */
/* ------------------------------------------------------------------------- */
template <class T>
   void instantiate_vector (matrix1D<T> v) {
      matrix1D<T>      a;
      matrix1D<double> r;
      matrix1D<int>    ir;
      double           d;
      T                Taux;
      
      // General functions for multidimensional arrays
      a==a;
      a!=a;
      a=1-a;
      a=a-1;
      a=a*a;
      a.print_stats();
      a.compute_double_minmax(d,d);
      a.compute_double_minmax(d,d,r,r);
      a.compute_double_minmax(d,d,ir,ir);
      a.compute_stats(d,d,Taux,Taux,ir,ir);
      a.compute_max();
      a.compute_min();
      a.compute_avg();
      a.compute_stddev();
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
	  a.isBorder(0);
	  matrix1D<int> pixel(1); a.isBorder(pixel);

      a.outside(0);
      a.intersects(a);
      a.intersects(r,r);
      a.isCorner(r);
      a.patch(a);
      cout << a;
      a.window(0,1);
      int imax;
      a.max_index(imax);
      a.min_index(imax);
      cut_to_common_size(a,a);
   
      // Specific for vectors
      matrix1D<int> indx;
      a.init_linear(0,5);
      a.reverse();
      a=a.sort();
      indx=a.index_sort();
      sort_two_vectors(a,a);
}

void instantiate_complex_vector () {
      matrix1D<double_complex> a;
      matrix1D<double>         r;
      
      // General functions for multidimensional arrays
      a=1.0-a;
      a=a-1.0;
      a=a*a;
      a.print_shape();
      a.outside(r);
      a.outside(0);
      a.intersects(a);
      a.intersects(r,r);
      a.isCorner(r);
      a.patch(a);
      cout << a;
      a.window(0,1);
      cut_to_common_size(a,a);
   
      // Specific for vectors
      a.reverse();
}

void instantiate1D() {
   matrix1D<char>           V0; instantiate_vector(V0);
   matrix1D<short>          V1; instantiate_vector(V1);
   matrix1D<int>            V2; instantiate_vector(V2);
   matrix1D<float>          V3; instantiate_vector(V3);
   matrix1D<double>         V4; instantiate_vector(V4);
}
