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
/* MATRICES 3D                                                               */
/* ------------------------------------------------------------------------- */
#include "../xmippMatrices3D.hh"
#include "../xmippArgs.hh"
#include "../xmippGeometry.hh"

#include <stdio.h>

/* ************************************************************************* */
/* IMPLEMENTATIONS                                                           */
/* ************************************************************************* */
#define maT matrix3D<T>
#define ma  matrix3D
#include "MultidimBasic.inc"
#undef ma
#undef maT
   // This file contain several general functions for arithmetic operations
   // They are defined for a general multidimensional_array<T>, but
   // we have "redirected" maT to matrix3D<T>, this way they will do OK in
   // this matrix3D library

/* Print shape ------------------------------------------------------------- */
template <class T>
   void VT::print_shape(ostream &out) const {
   out << "Size(Z,Y,X): " << ZSIZE(*this) << "x" << YSIZE(*this) << "x"
       << XSIZE(*this)
       << " k=[" << STARTINGZ(*this) << ".." << FINISHINGZ(*this) << "]"
       << " i=[" << STARTINGY(*this) << ".." << FINISHINGY(*this) << "]"
       << " j=[" << STARTINGX(*this) << ".." << FINISHINGX(*this) << "]";
}

/* Outside ----------------------------------------------------------------- */
template <class T>
   bool VT::outside(const matrix1D<double> &v) const _THROW {
   if (XSIZE(v)<3)
      REPORT_ERROR(1,"Outside: index vector has got not enough components");
   return (XX(v)<STARTINGX(*this) || XX(v)>FINISHINGX(*this) ||
           YY(v)<STARTINGY(*this) || YY(v)>FINISHINGY(*this) ||
           ZZ(v)<STARTINGZ(*this) || ZZ(v)>FINISHINGZ(*this));
}

template <class T>
   bool VT::outside(int k, int i, int j) const {
   return (j<STARTINGX(*this) || j>FINISHINGX(*this) ||
           i<STARTINGY(*this) || i>FINISHINGY(*this) ||
           k<STARTINGZ(*this) || k>FINISHINGZ(*this));
}

/* Intersects -------------------------------------------------------------- */
template <class T>
   bool VT::intersects(const VT &m) const {
      return intersects(STARTINGZ(m),STARTINGY(m),STARTINGX(m),
         XSIZE(m)-1,YSIZE(m)-1,ZSIZE(m)-1);
}

template <class T>
   bool VT::intersects(const matrix1D<double> &corner1,
      const matrix1D<double> &corner2) const _THROW {
       if (XSIZE(corner1)!=2 || XSIZE(corner2)!=2)
          REPORT_ERROR(1002,"intersects 1D: corner sizes are not 1");
       return intersects(XX(corner1),YY(corner1),ZZ(corner1),
          XX(corner2)-XX(corner1),YY(corner2)-YY(corner1),
          ZZ(corner2)-ZZ(corner1));
}


template <class T>
   bool VT::intersects(double x0, double y0, double z0, double xdim, double ydim,
      double zdim) const {
   SPEED_UP_temps;
   spduptmp0=MAX(STARTINGZ(*this), z0);
   spduptmp1=MIN(FINISHINGZ(*this),z0+zdim);
   if (spduptmp0>spduptmp1) return FALSE;
   
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
   bool VT::isCorner(const matrix1D<double> &v) _THROW {
   if (XSIZE(v)<3)
      REPORT_ERROR(1,"isCorner: index vector has got not enough components");
   return ((ZZ(v)==STARTINGZ(*this)  && XX(v)==STARTINGX(*this)  && YY(v)==STARTINGY(*this))  ||
           (ZZ(v)==STARTINGZ(*this)  && XX(v)==STARTINGX(*this)  && YY(v)==FINISHINGY(*this)) ||
           (ZZ(v)==STARTINGZ(*this)  && XX(v)==FINISHINGX(*this) && YY(v)==STARTINGY(*this))  ||
           (ZZ(v)==STARTINGZ(*this)  && XX(v)==FINISHINGX(*this) && YY(v)==FINISHINGY(*this)) ||
           (ZZ(v)==FINISHINGZ(*this) && XX(v)==STARTINGX(*this)  && YY(v)==STARTINGY(*this))  ||
           (ZZ(v)==FINISHINGZ(*this) && XX(v)==STARTINGX(*this)  && YY(v)==FINISHINGY(*this)) ||
           (ZZ(v)==FINISHINGZ(*this) && XX(v)==FINISHINGX(*this) && YY(v)==STARTINGY(*this))  ||
           (ZZ(v)==FINISHINGZ(*this) && XX(v)==FINISHINGX(*this) && YY(v)==FINISHINGY(*this)));
}

/* Border ------------------------------------------------------------------ */
template <class T>
   bool VT::isBorder(const matrix1D<int> &v) _THROW 
{
   if (XSIZE(v)<3)
      REPORT_ERROR(1,"isBorder: index vector has got not enough components");
   return  isBorder(ZZ(v),YY(v),XX(v));
}

template <class T>
   bool VT::isBorder(int k,int i,int j) 
{
   return (j==STARTINGX(*this)  || j==FINISHINGX(*this)  ||
           k==STARTINGZ(*this)  || k==FINISHINGZ(*this)  ||
           i==STARTINGY(*this)  || i==FINISHINGY(*this));
}

/* Patch ------------------------------------------------------------------- */
template <class T>
   void VT::patch(const VT &patch_array, char operation) {
      SPEED_UP_temps;
      FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(patch_array,*this)
         switch (operation) {
            case '=': VOL_ELEM(*this,k,i,j) =VOL_ELEM(patch_array,k,i,j); break;
            case '+': VOL_ELEM(*this,k,i,j)+=VOL_ELEM(patch_array,k,i,j); break;
            case '-': VOL_ELEM(*this,k,i,j)-=VOL_ELEM(patch_array,k,i,j); break;
            case '*': VOL_ELEM(*this,k,i,j)*=VOL_ELEM(patch_array,k,i,j); break;
            case '/': VOL_ELEM(*this,k,i,j)/=VOL_ELEM(patch_array,k,i,j); break;
         }
}

/* Output stream ----------------------------------------------------------- */
template <class T>
   ostream& operator << (ostream& ostrm, const VT& v) {   
   if (v.xdim==0)
      ostrm << "NULL matrix3D\n";
   else
      ostrm << endl;
      double max_val=ABS(MULTIDIM_ELEM(v,0));
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(v)
         max_val=MAX(max_val,ABS(MULTIDIM_ELEM(v,i)));
      int prec=best_prec(max_val,10);
      
      for (int k=STARTINGZ(v); k<=FINISHINGZ(v); k++) {
         ostrm << "Slice No. " << k << endl;
         for (int i=STARTINGY(v); i<=FINISHINGY(v); i++) {
            for (int j=STARTINGX(v); j<=FINISHINGX(v); j++) {
               ostrm << FtoA((double)VOL_ELEM(v,k,i,j),10,prec) << ' ';
            }
            ostrm << endl;
         }
      }
      
   return ostrm;
}

template <>
ostream& operator << (ostream& ostrm, const matrix3D<double_complex>& v) {
   if (v.xdim==0)
      ostrm << "NULL matrix3D\n";
   else
      ostrm << endl;
      
      for (int k=STARTINGZ(v); k<=FINISHINGZ(v); k++) {
         ostrm << "Slice No. " << k << endl;
         for (int i=STARTINGY(v); i<=FINISHINGY(v); i++) {
            for (int j=STARTINGX(v); j<=FINISHINGX(v); j++) {
               ostrm << VOL_ELEM(v,k,i,j) << ' ';
            }
            ostrm << endl;
         }
      }
      
   return ostrm;
}

/* Interpolated element ---------------------------------------------------- */
template <class T>
   T VT::interpolated_elem(double x, double y, double z, T outside_value) {
    int x0 = FLOOR(x); double fx = x - x0; int x1=x0+1;
    int y0 = FLOOR(y); double fy = y - y0; int y1=y0+1;
    int z0 = FLOOR(z); double fz = z - z0; int z1=z0+1;

    T d000 = (outside(z0,y0,x0)==TRUE) ? outside_value : VOL_ELEM(*this,z0,y0,x0);
    T d001 = (outside(z0,y0,x1)==TRUE) ? outside_value : VOL_ELEM(*this,z0,y0,x1);
    T d010 = (outside(z0,y1,x0)==TRUE) ? outside_value : VOL_ELEM(*this,z0,y1,x0);
    T d011 = (outside(z0,y1,x1)==TRUE) ? outside_value : VOL_ELEM(*this,z0,y1,x1);
    T d100 = (outside(z1,y0,x0)==TRUE) ? outside_value : VOL_ELEM(*this,z1,y0,x0);
    T d101 = (outside(z1,y0,x1)==TRUE) ? outside_value : VOL_ELEM(*this,z1,y0,x1);
    T d110 = (outside(z1,y1,x0)==TRUE) ? outside_value : VOL_ELEM(*this,z1,y1,x0);
    T d111 = (outside(z1,y1,x1)==TRUE) ? outside_value : VOL_ELEM(*this,z1,y1,x1);

    double dx00 = LIN_INTERP(fx, (double) d000, (double) d001);
    double dx01 = LIN_INTERP(fx, (double) d100, (double) d101);
    double dx10 = LIN_INTERP(fx, (double) d010, (double) d011);
    double dx11 = LIN_INTERP(fx, (double) d110, (double) d111);
    double dxy0 = LIN_INTERP(fy, (double) dx00, (double) dx10);
    double dxy1 = LIN_INTERP(fy, (double) dx01, (double) dx11);

    return (T) LIN_INTERP(fz, dxy0, dxy1);
}

// Special case for complex numbers
double_complex matrix3D<double_complex>::interpolated_elem(
   double x, double y, double z, double_complex outside_value) {
    int x0 = FLOOR(x); double fx = x - x0; int x1=x0+1;
    int y0 = FLOOR(y); double fy = y - y0; int y1=y0+1;
    int z0 = FLOOR(z); double fz = z - z0; int z1=z0+1;

    double_complex d000 = outside(z0,y0,x0) ? outside_value : dVkij(*this,z0,y0,x0);
    double_complex d001 = outside(z0,y0,x1) ? outside_value : dVkij(*this,z0,y0,x1);
    double_complex d010 = outside(z0,y1,x0) ? outside_value : dVkij(*this,z0,y1,x0);
    double_complex d011 = outside(z0,y1,x1) ? outside_value : dVkij(*this,z0,y1,x1);
    double_complex d100 = outside(z1,y0,x0) ? outside_value : dVkij(*this,z1,y0,x0);
    double_complex d101 = outside(z1,y0,x1) ? outside_value : dVkij(*this,z1,y0,x1);
    double_complex d110 = outside(z1,y1,x0) ? outside_value : dVkij(*this,z1,y1,x0);
    double_complex d111 = outside(z1,y1,x1) ? outside_value : dVkij(*this,z1,y1,x1);

    double_complex dx00 = LIN_INTERP(fx, d000, d001);
    double_complex dx01 = LIN_INTERP(fx, d100, d101);
    double_complex dx10 = LIN_INTERP(fx, d010, d011);
    double_complex dx11 = LIN_INTERP(fx, d110, d111);
    double_complex dxy0 = LIN_INTERP(fy, dx00, dx10);
    double_complex dxy1 = LIN_INTERP(fy, dx01, dx11);

    return LIN_INTERP(fz, dxy0, dxy1);
}

/* Get slice --------------------------------------------------------------- */
template <class T>
   void VT::getSlice(int k, mT &M, char axis) const _THROW {
   if (xdim==0) {M.clear(); return;}   
   switch (axis) {
      case 'Z':
         if (k<zinit || k>=zinit+zdim)
            REPORT_ERROR(1203,"Slice: matrix3D subscript (k) out of range");

         k=k-zinit;
         M.resize(ydim,xdim);
         FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(M)
               DIRECT_MAT_ELEM(M,i,j)=DIRECT_VOL_ELEM(*this,k,i,j);
         STARTINGX(M)=STARTINGX(*this);
         STARTINGY(M)=STARTINGY(*this);
         break;
      case 'Y':
         if (k<yinit || k>=yinit+ydim)
            REPORT_ERROR(1203,"Slice: matrix3D subscript (i) out of range");

         k=k-yinit;
         M.resize(zdim,xdim);
         FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(M)
               DIRECT_MAT_ELEM(M,i,j)=DIRECT_VOL_ELEM(*this,i,k,j);
         STARTINGX(M)=STARTINGX(*this);
         STARTINGY(M)=STARTINGZ(*this);
         break;
      case 'X':
         if (k<xinit || k>=xinit+xdim)
            REPORT_ERROR(1203,"Slice: matrix3D subscript (j) out of range");

         k=k-xinit;
         M.resize(zdim,ydim);
         FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(M)
               DIRECT_MAT_ELEM(M,i,j)=DIRECT_VOL_ELEM(*this,i,j,k);
         STARTINGX(M)=STARTINGY(*this);
         STARTINGY(M)=STARTINGZ(*this);
         break;
      default:
         REPORT_ERROR(1205,(string)"Slice: not supported axis "+axis);
   }
}

/* Set slice --------------------------------------------------------------- */
template <class T>
   void VT::setSlice(int k, const mT &v) _THROW {
   if (xdim==0) return;
   if (k<zinit || k>=zinit+zdim)
      REPORT_ERROR(1203,"setSlice: matrix3D subscript (k) out of range");
   if (v.RowNo()!=ydim || v.ColNo()!=xdim)
      REPORT_ERROR(1202,"setSlice: matrix3D dimensions different from the matrix ones");

   k=k-zinit;
   FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(v)
         DIRECT_VOL_ELEM(*this,k,i,j)=DIRECT_MAT_ELEM(v,i,j);
}

/* Reverse X --------------------------------------------------------------- */
template <class T>
   void VT::reverseX() {
   for (int k=0; k<zdim; k++)
      for (int i=0; i<ydim; i++)
         for (int j=0; j<=(int)(xdim-1)/2; j++) {
            T aux;
            if (k==0 && i==0)
               cout << "Changing " << j << " " << XSIZE(*this)-1-j << endl;
            SWAP(DIRECT_VOL_ELEM(*this,k,i,j),
                 DIRECT_VOL_ELEM(*this,k,i,XSIZE(*this)-1-j),aux);
         }
   STARTINGX(*this)=-FINISHINGX(*this);
}

/* Reverse Y --------------------------------------------------------------- */
template <class T>
   void VT::reverseY() {
   for (int k=0; k<zdim; k++)
      for (int i=0; i<=(int)(ydim-1)/2; i++)
         for (int j=0; j<xdim; j++) {
            T aux;
            SWAP(DIRECT_VOL_ELEM(*this,k,i,j),
                 DIRECT_VOL_ELEM(*this,k,YSIZE(*this)-1-i,j),aux);
         }
   STARTINGY(*this)=-FINISHINGY(*this);
}

/* Reverse Z --------------------------------------------------------------- */
template <class T>
   void VT::reverseZ() {
   for (int k=0; k<=(int)(zdim-1)/2; k++)
      for (int i=0; i<ydim; i++)
         for (int j=0; j<xdim; j++) {
            T aux;
            SWAP(DIRECT_VOL_ELEM(*this,k,i,j),
                 DIRECT_VOL_ELEM(*this,ZSIZE(*this)-1-k,i,j),aux);
         }
   STARTINGZ(*this)=-FINISHINGZ(*this);
}

/* Window ------------------------------------------------------------------ */
template <class T>
   void VT::window(int z0, int y0, int x0, int zF, int yF, int xF) {
   VT result(zF-z0+1, yF-y0+1, xF-x0+1);
   result.zinit=z0;
   result.yinit=y0;
   result.xinit=x0;
   
   for (int k=z0; k<=zF; k++)
      for (int i=y0; i<=yF; i++)
         for (int j=x0; j<=xF; j++)
             if ((k>=zinit && k<=zinit+zdim-1) &&
                 (i>=yinit && i<=yinit+ydim-1) &&
                 (j>=xinit && j<=xinit+xdim-1))
                     VOL_ELEM(result,k,i,j)=VOL_ELEM(*this,k,i,j);
             else
                     VOL_ELEM(result,k,i,j)=0;
   *this=result;
}

/* Apply a geometrical transformation -------------------------------------- */
// It is translated from function transforma in Lib/varios2.c
// We should check which one performs better.
template <class T>
   void apply_geom(VT &V2, matrix2D<double> A, const VT &V1, bool inv,
      bool wrap) _THROW {
   int m1, n1, o1, m2, n2, o2;
   double x, y, z, xp, yp, zp;
   double minxp, minyp, maxxp, maxyp, minzp, maxzp;
   int   cen_x, cen_y, cen_z, cen_xp, cen_yp, cen_zp;
   double wx, wy, wz; // Weights in X,Y,Z directions for bilinear interpolation
   
   if ((XSIZE(A)!=4) || (YSIZE(A)!=4))
      REPORT_ERROR(1102,"Apply_geom3D: geometrical transformation is not 4x4");
   if (A.IsIdent()) {V2=V1; return;}
   if (XSIZE(V1)==0)  {V2.clear(); return;}
   
   if (!inv) A=A.inv();
   
   if (XSIZE(V2)==0) V2.resize(V1);
   
   // For scalings the output matrix is resized outside to the final
   // size instead of being resized inside the routine with the
   // same size as the input matrix
   if (XSIZE(V2)==0) V2.resize(V1);
   
   // Find center of matrix3D
   cen_z  =(int)(V2.zdim/2);
   cen_y  =(int)(V2.ydim/2);
   cen_x  =(int)(V2.xdim/2);
   cen_zp =(int)(V1.zdim/2);
   cen_yp =(int)(V1.ydim/2);
   cen_xp =(int)(V1.xdim/2);
   minxp  = -cen_xp;
   minyp  = -cen_yp;
   minzp  = -cen_zp;
   maxxp  = V1.xdim-cen_xp-1;
   maxyp  = V1.ydim-cen_yp-1;
   maxzp  = V1.zdim-cen_zp-1;
    
   // Now we go from the output matrix3D to the input matrix3D, ie, for any voxel
   // in the output matrix3D we calculate which are the corresponding ones in
   // the original matrix3D, make an interpolation with them and put this value
   // at the output voxel
   
   // V2 is not initialised to 0 because all its pixels are rewritten
   for (int k=0; k<V2.zdim; k++)
      for (int i=0; i<V2.ydim; i++) {
         // Calculate position of the beginning of the row in the output matrix3D
         x= -cen_x;
         y=i-cen_y;
         z=k-cen_z;

         // Calculate this position in the input image according to the
         // geometrical transformation
         // they are related by
         // coords_output(=x,y) = A * coords_input (=xp,yp)
         xp=x*dMij(A,0,0) + y*dMij(A,0,1) + z*dMij(A,0,2) + dMij(A,0,3);
         yp=x*dMij(A,1,0) + y*dMij(A,1,1) + z*dMij(A,1,2) + dMij(A,1,3);
         zp=x*dMij(A,2,0) + y*dMij(A,2,1) + z*dMij(A,2,2) + dMij(A,2,3);

         for (int j=0; j<V2.xdim; j++) {
            bool interp;
            T tmp;

            // If the point is outside the volume, apply a periodic extension
            // of the volume, what exits by one side enters by the other
            interp=TRUE;
            if (wrap) {
               if (xp<minxp || xp>maxxp) xp=realWRAP(xp, minxp, maxxp);
               if (yp<minyp || yp>maxyp) yp=realWRAP(yp, minyp, maxyp);
               if (zp<minzp || zp>maxzp) zp=realWRAP(zp, minzp, maxzp);
            } else {
               if (xp<minxp || xp>maxxp) interp=FALSE;
               if (yp<minyp || yp>maxyp) interp=FALSE;
               if (zp<minzp || zp>maxzp) interp=FALSE;
            }

            if (interp) {            
               // Calculate the integer position in input volume, be careful
               // that it is not the nearest but the one at the top left corner
               // of the interpolation square. Ie, (0.7,0.7) would give (0,0)
               // Calculate also weights for point m1+1,n1+1
               wx=xp+cen_xp; m1=(int) wx; wx=wx-m1; m2=m1+1;
               wy=yp+cen_yp; n1=(int) wy; wy=wy-n1; n2=n1+1;
               wz=zp+cen_zp; o1=(int) wz; wz=wz-o1; o2=o1+1;

               // Perform interpolation
               // if wx == 0 means that the rightest point is useless for this
               // interpolation, and even it might not be defined if m1=xdim-1
               // The same can be said for wy.
                                              tmp  = (T) ((1-wz)*(1-wy)*(1-wx)*dVkij(V1,o1,n1,m1));
               if (wx!=0 && m2<V1.xdim)       tmp += (T) ((1-wz)*(1-wy)*   wx *dVkij(V1,o1,n1,m2));
               if (wy!=0 && n2<V1.ydim){      tmp += (T) ((1-wz)*   wy *(1-wx)*dVkij(V1,o1,n2,m1));
                  if (wx!=0 && m2<V1.xdim)    tmp += (T) ((1-wz)*   wy *   wx *dVkij(V1,o1,n2,m2));
               }
               if (wz!=0 && o2<V1.zdim){
                                              tmp += (T) (   wz *(1-wy)*(1-wx)*dVkij(V1,o1,n1,m1));
                  if (wx!=0 && m2<V1.xdim)    tmp += (T) (   wz *(1-wy)*   wx *dVkij(V1,o1,n1,m2));
                  if (wy!=0 && n2<V1.ydim){   tmp += (T) (   wz *   wy *(1-wx)*dVkij(V1,o1,n2,m1));
                     if (wx!=0 && m2<V1.xdim) tmp += (T) (   wz *   wy *   wx *dVkij(V1,o1,n2,m2));
                  }
               }
               dVkij(V2,k,i,j) = tmp;
            }

            // Compute new point inside input image
            xp += dMij(A,0,0);
            yp += dMij(A,1,0);
            zp += dMij(A,2,0);
         }
      }
}

/* Geometrical operations -------------------------------------------------- */
template <class T>
   void VT::scale_to_size(int Zdim, int Ydim, int Xdim, VT &result) const
      {matrix2D<double> temp(4,4);
       temp.init_identity();
       DIRECT_MAT_ELEM(temp,0,0)=(double)Xdim/(double)xdim;
       DIRECT_MAT_ELEM(temp,1,1)=(double)Ydim/(double)ydim;
       DIRECT_MAT_ELEM(temp,2,2)=(double)Zdim/(double)zdim;
       result.resize(Zdim,Ydim,Xdim);
       apply_geom(result,temp,*this,IS_NOT_INV,WRAP);}

template <class T>
   void VT::self_translate_center_of_mass_to_center(bool wrap) {
      set_Xmipp_origin();
      matrix1D<double> center;
      center_of_mass(center);
      center*=-1;
      self_translate(center,wrap);
   }

/* Max index --------------------------------------------------------------- */
template <class T>
   void VT::max_index(int &kmax, int &imax, int &jmax) const {
   if (XSIZE(*this)==0) {kmax=imax=jmax=-1; return;}
   kmax=imax=jmax=0;
   T   max=VOL_ELEM(*this,kmax,imax,jmax);
   FOR_ALL_ELEMENTS_IN_MATRIX3D(*this)
      if (VOL_ELEM(*this,k,i,j)>max)
         {max=VOL_ELEM(*this,kmax,imax,jmax); kmax=k; imax=i; jmax=j;}
}

/* Min index --------------------------------------------------------------- */
template <class T>
   void VT::min_index(int &kmin, int &imin, int &jmin) const {
   if (XSIZE(*this)==0) {kmin=imin=jmin=-1; return;}
   kmin=imin=jmin=0;
   T   min=VOL_ELEM(*this,kmin,imin,jmin);
   FOR_ALL_ELEMENTS_IN_MATRIX3D(*this)
      if (VOL_ELEM(*this,k,i,j)>min)
         {min=VOL_ELEM(*this,kmin,imin,jmin); kmin=k; imin=i; jmin=j;}
}

/* For all slices ---------------------------------------------------------- */
template <class T>
   VT VT::for_all_slices(mT (*f)(mT&)) const {
   VT temp;
   
   if (xdim==0) {temp.clear(); return temp;};
   temp.copy_shape(*this);
    
   for (int k=STARTINGZ(*this); k<=FINISHINGZ(*this); k++) {
      mT aux;
      getSlice(k,aux);
      aux=(*f)(aux);
      temp.setSlice(k,aux);
   }
   return temp;
}

template <class T>
   vT VT::for_all_slices(T (*f)(mT&)) const {
   vT temp;
   
   if (ZSIZE(*this)==0) {temp.clear(); return temp;};
   temp.resize(ZSIZE(*this));
   STARTINGX(temp)=STARTINGZ(*this);
   
   for (int k=STARTINGZ(*this); k<=FINISHINGZ(*this); k++) {
      mT aux;
      getSlice(k,aux);
      VEC_ELEM(temp,k)=(*f)(aux);
   }
   
   return temp;
}

/* Cut to common size ------------------------------------------------------ */
template <class T>
   void cut_to_common_size(VT &V1, VT &V2) {
   int z0=MAX(STARTINGZ(V1) ,STARTINGZ(V2));
   int zF=MIN(FINISHINGZ(V1),FINISHINGZ(V2));
   int y0=MAX(STARTINGY(V1) ,STARTINGY(V2));
   int yF=MIN(FINISHINGY(V1),FINISHINGY(V2));
   int x0=MAX(STARTINGX(V1) ,STARTINGX(V2));
   int xF=MIN(FINISHINGX(V1),FINISHINGX(V2));
   V1.window(z0,y0,x0,zF,yF,xF);
   V2.window(z0,y0,x0,zF,yF,xF);
}


/* Radial average ---------------------------------------------------------- */
template <class T>
void radial_average(const matrix3D<T> &m, const matrix1D<int> &center_of_rot,
                    matrix1D<T> &radial_mean) _THROW
{
   matrix1D<double> idx(3);
   
   /* First determine the maximum distance that one should expect,
      to set the dimension of the radial average vector */
   matrix1D<int> distances(8);
   double z=STARTINGZ(m)-ZZ(center_of_rot);
   double y=STARTINGY(m)-YY(center_of_rot);
   double x=STARTINGX(m)-XX(center_of_rot);
   distances(0)=(int)floor(sqrt(x*x+y*y+z*z));
   x=FINISHINGX(m)-XX(center_of_rot);
   distances(1)=(int)floor(sqrt(x*x+y*y+z*z));
   y=FINISHINGY(m)-YY(center_of_rot);
   distances(2)=(int)floor(sqrt(x*x+y*y+z*z));
   x=STARTINGX(m)-XX(center_of_rot);
   distances(3)=(int)floor(sqrt(x*x+y*y+z*z));
   z=FINISHINGZ(m)-ZZ(center_of_rot);
   distances(4)=(int)floor(sqrt(x*x+y*y+z*z));
   x=FINISHINGX(m)-XX(center_of_rot);
   distances(5)=(int)floor(sqrt(x*x+y*y+z*z));
   y=STARTINGY(m)-YY(center_of_rot);
   distances(6)=(int)floor(sqrt(x*x+y*y+z*z));
   x=STARTINGX(m)-XX(center_of_rot);
   distances(7)=(int)floor(sqrt(x*x+y*y+z*z));

   int dim=(int)CEIL(distances.compute_max())+1;
   
   // Define the vectors
   radial_mean.resize(dim);
   matrix1D<int> radial_count(dim);
   
   /* Perform the radial sum and count pixels that contribute to
      every distance */
   FOR_ALL_ELEMENTS_IN_MATRIX3D(m) {
      ZZ(idx)=k-ZZ(center_of_rot);
      YY(idx)=i-YY(center_of_rot);
      XX(idx)=j-XX(center_of_rot);
      // Determine distance to the center
      int distance=(int)floor(idx.module());
      // Sum te value to the pixels with the same distance  
      radial_mean(distance)+=m(k,i,j);
      // Count the pixel
      radial_count(distance)++;      	  	  
   }
   
   // Perform the mean
   FOR_ALL_ELEMENTS_IN_MATRIX1D(radial_mean)
      radial_mean(i)/=(T)radial_count(i);
}

/* Center of mass ---------------------------------------------------------- */
template <class T>
   void VT::center_of_mass(matrix1D<double> &center) {
      center.init_zeros(3);
      double mass=0;
      FOR_ALL_ELEMENTS_IN_MATRIX3D(*this) {
         XX(center)+=j*VOL_ELEM(*this,k,i,j);
         YY(center)+=i*VOL_ELEM(*this,k,i,j);
         ZZ(center)+=k*VOL_ELEM(*this,k,i,j);
	 mass+=VOL_ELEM(*this,k,i,j);
      }
      center/=mass;
   }

/* ------------------------------------------------------------------------- */
/* INSTANTIATE                                                               */
/* ------------------------------------------------------------------------- */
template <class T>
   T slice_first(mT &x) {return x(0,0);}

template <class T>
   mT slice_minus_first(mT &x) {return x-x(0,0);}

template <class T>
   void instantiate_matrix3D (matrix3D<T> v) {
      matrix3D<T>      a;
      matrix1D<T>      vectorT;
      matrix1D<double> r;
      double           d;
      
      // General functions for multidimensional arrays
      a==a;
      a=1-a;
      a=a-1;
      a=a*a;
      a.print_stats();
      a.compute_max();
      a.compute_min();
      a.compute_avg();
      a.compute_stddev();
      a.compute_double_minmax(d,d);
      a.range_adjust(0,1);
      a.statistics_adjust(0,1);
      a.effective_range(99);
      a.init_random(0,1);
      a.add_noise(0,1);      
      a.threshold("abs_above",1,0);
      a.count_threshold("abs_above",1,0);

      a.print_shape();
      a.outside(r);
      a.outside(0,0,0);
	  a.isBorder(0,0,0);
      matrix1D<int> pixel(3);	  a.isBorder(pixel);
      a.intersects(a);
      a.intersects(r,r);
      a.isCorner(r);
      a.patch(a);
      cout << a;
      a.window(0,0,0,1,1,1);
      cut_to_common_size(a,a);
   
      // Specific for volumes
      matrix2D<T> aux;
      matrix2D<double> A;
      a.interpolated_elem(3.5,3.5,3.5);
      a.getSlice(0,aux);
      a.setSlice(0,aux);
      a.reverseX();
      a.reverseY();
      a.reverseZ();
      apply_geom(a, A, a, IS_NOT_INV,DONT_WRAP);
      int imax;
      a.max_index(imax,imax,imax);
      a.min_index(imax,imax,imax);
      a.scale_to_size(32,32,32,a);
      a.self_translate_center_of_mass_to_center();
      a.for_all_slices(&slice_minus_first);
      a.for_all_slices(&slice_first);
      cut_to_common_size(a,a);
      radial_average(a,pixel,vectorT);
}

void instantiate_complex_matrix3D () {
      matrix3D<double_complex> a;
      matrix1D<double>         r;
      matrix1D<double_complex> vectorT;
      
      // General functions for multidimensional arrays
      a=1.0-a;
      a=a-1.0;
      a=a*a;
      a.print_shape();
      a.outside(r);
      a.outside(0,0,0);
      a.intersects(a);
      a.intersects(r,r);
      a.isCorner(r);
	  a.isBorder(0,0,0);
      matrix1D<int> pixel(3);	  a.isBorder(pixel);
      a.patch(a);
      cout << a;
      a.window(0,0,0,1,1,1);
      cut_to_common_size(a,a);
   
      // Specific for volumes
      matrix2D<double_complex> aux;
      matrix2D<double> A;
      a.interpolated_elem(3.5,3.5,3.5);
      a.getSlice(0,aux);
      a.setSlice(0,aux);
      a.reverseX();
      a.reverseY();
      a.reverseZ();
      apply_geom(a, A, a, IS_NOT_INV,DONT_WRAP);
      a.scale_to_size(32,32,32,a);
      a.for_all_slices(&slice_minus_first);
      a.for_all_slices(&slice_first);
      cut_to_common_size(a,a);
      radial_average(a,pixel,vectorT);
}

void instantiate3D() {
   // matrix3D<char>           V0; instantiate_matrix3D(V0);
   // matrix3D<short>          V1; instantiate_matrix3D(V1);
   matrix3D<int>            V2; instantiate_matrix3D(V2);
   matrix3D<float>          V3; instantiate_matrix3D(V3);
   matrix3D<double>         V4; instantiate_matrix3D(V4);
}
