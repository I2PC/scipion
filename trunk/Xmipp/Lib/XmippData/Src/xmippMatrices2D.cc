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

/* ************************************************************************* */
/* IMPLEMENTATIONS                                                           */
/* ************************************************************************* */
#define maT matrix2D<T>
#define ma  matrix2D
#include "MultidimBasic.inc"
#undef ma
#undef maT
// Special case for complex numbers
template <>
ostream& operator << (ostream& ostrm, const matrix2D< complex<double> > &v) {
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

// Special case for complex numbers
template <>
complex<double> matrix2D< complex<double> >::interpolated_elem(
   double x, double y, complex<double> outside_value) {
    int x0 = FLOOR(x); double fx = x - x0; int x1 = x0 + 1;
    int y0 = FLOOR(y); double fy = y - y0; int y1 = y0 + 1;
    
    complex<double> d00 = outside(y0,x0) ? outside_value : dMij(*this,y0,x0);
    complex<double> d10 = outside(y1,x0) ? outside_value : dMij(*this,y1,x0);
    complex<double> d11 = outside(y1,x1) ? outside_value : dMij(*this,y1,x1);
    complex<double> d01 = outside(y0,x1) ? outside_value : dMij(*this,y0,x1);

    complex<double> d0 = LIN_INTERP(fx, d00, d01);
    complex<double> d1 = LIN_INTERP(fx, d10, d11);    
    return LIN_INTERP(fy, d0, d1);
}

/* Solve Cx=d, nonnegative x */
double solve_nonneg(const matrix2D<double> &C, const matrix1D<double> &d,
   matrix1D<double> &result) {
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

/* Solve Ax=b, A definite positive and symmetric --------------------------- */
void solve_via_Cholesky(const matrix2D<double> &A, const matrix1D<double> &b,
   matrix1D<double> &result) {
   matrix2D<double> Ap=A;
   matrix1D<double> p(XSIZE(A));
   result.resize(XSIZE(A));
   choldc(Ap.adapt_for_numerical_recipes2(),XSIZE(A),
      p.adapt_for_numerical_recipes());
   cholsl(Ap.adapt_for_numerical_recipes2(),XSIZE(A),
      p.adapt_for_numerical_recipes(),b.adapt_for_numerical_recipes(),
      result.adapt_for_numerical_recipes());
}

// Special case for complex numbers
template <>
   void apply_geom_Bspline(matrix2D< complex<double> > &M2,
      matrix2D<double> A, const matrix2D< complex<double> > &M1,
   int Splinedegree, bool inv, bool wrap, complex<double> outside) {
   REPORT_ERROR(1,"apply_geom_Bspline: Not yet implemented for complex matrices\n");
}

/* Interface to numerical recipes: svbksb ---------------------------------- */
void svbksb(matrix2D<double> &u,matrix1D<double> &w,matrix2D<double> &v,
             matrix1D<double> &b,matrix1D<double> &x) {
   // Call to the numerical recipes routine. Results will be stored in X
   svbksb(u.adapt_for_numerical_recipes2(),
          w.adapt_for_numerical_recipes(),
          v.adapt_for_numerical_recipes2(),
          u.RowNo(),u.ColNo(),
          b.adapt_for_numerical_recipes(),
          x.adapt_for_numerical_recipes());
}

template <>
void matrix2D< complex<double> >::scale_to_size_Bspline(int Splinedegree,
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

/* Is diagonal ------------------------------------------------------------- */
template <>
bool matrix2D< complex<double> >::IsDiagonal() const
   {if (XSIZE(*this)!=YSIZE(*this)) return false;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
          if (i!=j && abs(DIRECT_MAT_ELEM(*this,i,j))>XMIPP_EQUAL_ACCURACY)
             return false;
    return true;}

/* Is Scalar --------------------------------------------------------------- */
template <>
bool matrix2D< complex<double> >::IsScalar() const
   {if (!IsDiagonal()) return false;
    for (int i=1; i<YSIZE(*this); i++)
       if (abs(DIRECT_MAT_ELEM(*this,i,i)-DIRECT_MAT_ELEM(*this,0,0))>
           XMIPP_EQUAL_ACCURACY) return false;
    return true;}

/* Is Symmetric ------------------------------------------------------------ */
template <>
bool matrix2D< complex<double> >::IsSymmetric() const
   {if (XSIZE(*this)!=YSIZE(*this)) return false;
    for (int i=0; i<YSIZE(*this); i++)
       for (int j=i+1; j<XSIZE(*this); j++)
          if (abs(DIRECT_MAT_ELEM(*this,i,j)-DIRECT_MAT_ELEM(*this,j,i))>
              XMIPP_EQUAL_ACCURACY) return false;
    return true;}

/* Is Skew symmetric ------------------------------------------------------- */
template <>
bool matrix2D< complex<double> >::IsSkewSymmetric() const
   {if (XSIZE(*this)!=YSIZE(*this)) return false;
    for (int i=0; i<YSIZE(*this); i++)
       for (int j=i+1; j<XSIZE(*this); j++)
          if (abs(DIRECT_MAT_ELEM(*this,i,j)+DIRECT_MAT_ELEM(*this,j,i))>
              XMIPP_EQUAL_ACCURACY) return false;
    return true;}

/* Is Upper triangular ----------------------------------------------------- */
template <>
bool matrix2D< complex<double> >::IsUpperTriangular() const
   {if (XSIZE(*this)!=YSIZE(*this)) return false;
    for (int i=1; i<YSIZE(*this); i++)
       for (int j=0; j<i-1; j++)
          if (abs(DIRECT_MAT_ELEM(*this,i,j))>XMIPP_EQUAL_ACCURACY)
             return false;
    return true;}

/* Is Lower triangular ----------------------------------------------------- */
template <>
bool matrix2D< complex<double> >::IsLowerTriangular() const
   {if (XSIZE(*this)!=YSIZE(*this)) return false;
    for (int i=1; i<YSIZE(*this); i++)
       for (int j=i+1; j<XSIZE(*this); j++)
          if (abs(DIRECT_MAT_ELEM(*this,i,j))>XMIPP_EQUAL_ACCURACY)
             return false;
    return true;}

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
matrix2D<double> translation2D_matrix(matrix1D<double> v) {
   if (v.get_dim()!=2)
      REPORT_ERROR(1002,"Translation2D_matrix: vector is not in R2");

   matrix2D<double> result(3,3);

   result.init_identity();
   DIRECT_MAT_ELEM(result,0,2)=XX(v);
   DIRECT_MAT_ELEM(result,1,2)=YY(v);

   return result;   
}

/* Rotation 3D around the system axes -------------------------------------- */
matrix2D<double> rot3D_matrix(double ang, char axis) {
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
matrix2D<double> align_with_Z(const matrix1D<double> &axis) {
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
matrix2D<double> translation3D_matrix(const matrix1D<double> &v) {
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
matrix2D<double> scale3D_matrix(const matrix1D<double> &sc) {
   if (XSIZE(sc)!=3)
      REPORT_ERROR(1002,"Scale3D_matrix: vector is not in R3");

   matrix2D<double> result(4,4);

   result.init_identity();
   DIRECT_MAT_ELEM(result,0,0)=XX(sc);
   DIRECT_MAT_ELEM(result,1,1)=YY(sc);
   DIRECT_MAT_ELEM(result,2,2)=ZZ(sc);

   return result;   
}

/* Quadratic form ---------------------------------------------------------- */
void eval_quadratic(const matrix1D<double> &x, const matrix1D<double> &c,
   const matrix2D<double> &H, double &val, matrix1D<double> &grad) {
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

/* Quadprog and Lsqlin ----------------------------------------------------- */
/* Structure to pass the objective function and constraints to cfsqp*/
typedef struct {
   matrix2D<double> C;
   matrix2D<double> D;
   matrix2D<double> A;
   matrix2D<double> B;
} CDAB;

/*/////////////////////////////////////////////////////////////////////////
      Internal functions used by the quadprog function
/////////////////////////////////////////////////////////////////////////*/
/* To calculate the value of the objective function */
void quadprog_obj32(int nparam,int j,double* x,double* fj,void* cd)
{
   CDAB* in = (CDAB *)cd;
   matrix2D<double> X;
   XSIZE(X)=1;
   YSIZE(X)=nparam;
   MULTIDIM_ARRAY(X)=x;
   
   matrix2D<double> result;
   result= 0.5*X.transpose()*in->C*X + in->D.transpose()*X;
   
   *fj=result(0,0);
   MULTIDIM_ARRAY(X)=NULL;
   XSIZE(X)=YSIZE(X)=0;
   return;
}

/* To calculate the value of the jth constraint */
void quadprog_cntr32(int nparam,int j,double* x,double* gj,void* cd)
{
   CDAB* in = (CDAB *)cd;
   *gj=0;
   for (int k=0; k<nparam; k++)
      *gj+=in->A(j-1,k)*x[k];
    *gj-=in->B(j-1,0);
    
   return;
}

/* To calculate the value of the derivative of objective function */
void quadprog_grob32(int nparam,int j,double* x,double* gradfj, void (*mydummy)(int, int, double*, double*, void*), void *cd)
{
   CDAB* in = (CDAB *)cd;
   matrix2D<double> X;
   XSIZE(X)=1;
   YSIZE(X)=nparam;
   MULTIDIM_ARRAY(X)=x;

   matrix2D<double> gradient;
   gradient=in->C*X+in->D;
   for (int k=0; k<nparam; k++)
      gradfj[k]=gradient(k,0);   

   MULTIDIM_ARRAY(X)=NULL;
   XSIZE(X)=YSIZE(X)=0;
   return;
}

/* To calculate the value of the derivative of jth constraint */
void quadprog_grcn32(int nparam, int j, double *x, double *gradgj,void (*mydummy)(int, int, double*, double*, void*), void *cd)
{
   CDAB* in = (CDAB *)cd;
   for (int k=0; k<nparam; k++)
      gradgj[k]=in->A(j-1,k);   
   return;
}   

/************************************************************************** 

   Solves Quadratic programming subproblem. 

  min 0.5*x'Cx + d'x   subject to:  A*x <= b 
   x                                Aeq*x=beq
      	             	      	    bl<=x<=bu  

**************************************************************************/
void quadprog(const matrix2D<double> &C, const matrix1D<double> &d,
   const matrix2D<double> &A,   const matrix1D<double> &b,
   const matrix2D<double> &Aeq, const matrix1D<double> &beq,
         matrix1D<double> &bl,        matrix1D<double> &bu,
         matrix1D<double> &x) {
   CDAB prm;
   prm.C=C;
   prm.D.from_vector(d);
   prm.A.init_zeros(YSIZE(A)+YSIZE(Aeq),XSIZE(A));
   prm.B.init_zeros(YSIZE(prm.A),1);

   
   // Copy Inequalities
   for (int i=0; i<YSIZE(A); i++) {
      for (int j=0; j<XSIZE(A); j++)
         prm.A(i,j)=A(i,j);
      prm.B(i,0)=b(i);
   }
   
   // Copy Equalities
   for (int i=0; i<YSIZE(Aeq); i++) {
      for (int j=0; j<XSIZE(Aeq); j++)
         prm.A(i+YSIZE(A),j)=Aeq(i,j);
      prm.B(i+YSIZE(A),0)=beq(i);
   }
   
   double bigbnd=1e30;
   // Bounds
   if (XSIZE(bl)==0) {
      bl.resize(XSIZE(C)); bl.init_constant(-bigbnd);
   }
   if (XSIZE(bu)==0) {
      bu.resize(XSIZE(C)); bu.init_constant( bigbnd);
   }
   
   // Define intermediate variables
   int    mode=100;    // CFSQP mode
   int    iprint=0;    // Debugging 
   int    miter=1000;  // Maximum number of iterations
   double eps=1e-4;   // Epsilon
   double epsneq=1e-4; // Epsilon for equalities
   double udelta=0.e0; // Finite difference approximation
                       // of the gradients. Not used in this function
   int    nparam=XSIZE(C); // Number of variables
   int    nf=1;            // Number of objective functions
   int    neqn=YSIZE(Aeq);          // Number of nonlinear equations
   int    nineqn=YSIZE(A);        // Number of nonlinear inequations
   int    nineq=YSIZE(A);  // Number of linear inequations
   int    neq=YSIZE(Aeq);  // Number of linear equations
   int    inform;
   int    ncsrl=0, ncsrn=0, nfsr=0, mesh_pts[]={0};

   if (XSIZE(x)==0) x.init_zeros(nparam);
   matrix1D<double> f(nf), g(nineq+neq), lambda(nineq+neq+nf+nparam);

   // Call the minimization routine
   cfsqp(nparam,nf,nfsr,nineqn,nineq,neqn,neq,ncsrl,ncsrn,mesh_pts,
         mode,iprint,miter,&inform,bigbnd,eps,epsneq,udelta,
	 MULTIDIM_ARRAY(bl),MULTIDIM_ARRAY(bu),
	 MULTIDIM_ARRAY(x),
         MULTIDIM_ARRAY(f),MULTIDIM_ARRAY(g),
	 MULTIDIM_ARRAY(lambda),
//	 quadprog_obj32,quadprog_cntr32,quadprog_grob32,quadprog_grcn32,
	 quadprog_obj32,quadprog_cntr32,grobfd,grcnfd,
	 (void*)&prm);

   #ifdef DEBUG
      if (inform==0) cout <<"SUCCESSFUL RETURN. \n";
      if (inform==1 || inform==2) cout << "\nINITIAL GUESS INFEASIBLE.\n";
      if (inform==3) printf("\n MAXIMUM NUMBER OF ITERATIONS REACHED.\n");
      if (inform>3) printf("\ninform=%d\n",inform);
   #endif
   
}


/************************************************************************** 

   Solves the least square problem 

  min 0.5*(Norm(C*x-d))   subject to:  A*x <= b 
   x                                Aeq*x=beq  
      	             	      	    bl<=x<=bu 
**************************************************************************/
void lsqlin(const matrix2D<double> &C, const matrix1D<double> &d,
   const matrix2D<double> &A,   const matrix1D<double> &b,
   const matrix2D<double> &Aeq, const matrix1D<double> &beq,
         matrix1D<double> &bl,        matrix1D<double> &bu,
         matrix1D<double> &x) {
   
   // Convert d to matrix2D for multiplication
   matrix2D<double> P;
   P.from_vector(d);
   P=-2*P.transpose()*C;
   P=P.transpose();
   
   //Convert back to vector for passing it to quadprog
   matrix1D<double> newd; 
   P.to_vector(newd);
      
   quadprog(C.transpose()*C,newd,A,b,Aeq,beq,bl,bu,x);
}
