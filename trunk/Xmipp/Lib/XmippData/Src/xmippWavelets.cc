/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Antonio Jose Rodriguez Sanchez (ajr@cnb.uam.es)
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
/* WAVELETS                                                                  */
/* ------------------------------------------------------------------------- */

#include <XmippData/xmippWavelets.hh>

/* Wavelet ----------------------------------------------------------------- */

// Set the DWT type --------------------------------------------------------
void set_DWT_type(int DWT_type) {
   pwtset(DWT_type);
}

// DWT ---------------------------------------------------------------------

template <class T>
   void DWT(const matrix1D<T> &v,
   matrix1D<double> &result, int isign)
   {
   unsigned long int nn[1];
   unsigned long int *ptr_nn=nn-1;

   type_cast(v, result);
   nn[0] = XSIZE(result);
   double *ptr_result=MULTIDIM_ARRAY(result)-1;
   wtn(ptr_result, ptr_nn, 1, isign, pwt);
   }

template <class T>
   void DWT(const matrix2D<T> &v,
   matrix2D<double> &result, int isign)
   {
   unsigned long int nn[2];
   unsigned long int *ptr_nn=nn-1;

   type_cast(v, result);
   nn[1] = YSIZE(result);
   nn[0] = XSIZE(result);
   double *ptr_result=MULTIDIM_ARRAY(result)-1;
   wtn(ptr_result, ptr_nn, 2, isign, pwt);
   }

template <class T>
   void DWT(const matrix3D<T> &v,
   matrix3D<double> &result, int isign)
   {
   unsigned long int nn[2];
   unsigned long int *ptr_nn=nn-1;

   type_cast(v, result);
   nn[2] = ZSIZE(result);
   nn[1] = YSIZE(result);
   nn[0] = XSIZE(result);
   double *ptr_result=MULTIDIM_ARRAY(result)-1;
   wtn(ptr_result, ptr_nn, 3, isign, pwt);
   }

// IDWT --------------------------------------------------------------------
void IDWT(const matrix1D<double> &v, matrix1D<double> &result) {
   DWT(v,result,-1);
}

void IDWT(const matrix2D<double> &v, matrix2D<double> &result) {
   DWT(v,result,-1);
}

void IDWT(const matrix3D<double> &v, matrix3D<double> &result) {
   DWT(v,result,-1);
}

// Instantiation -----------------------------------------------------------
template <class T>
   void instantiate_XmippWavelets1(T &t) {
      matrix1D<T> v1;
      matrix2D<T> v2;
      matrix3D<T> v3;
      DWT(v1,v1);
      DWT(v2,v2);
      DWT(v3,v3);
   }

void instantiate_XmippWavelets() {
   double d; instantiate_XmippWavelets1(d);
}

// Select block ------------------------------------------------------------
#define DWT_Imin(s,smax,l) (int)((l=='0')?0:pow(2.0,smax-s-1))
#define DWT_Imax(s,smax,l) (int)((l=='0')?pow(2.0,smax-s-1)-1:pow(2.0,smax-s)-1)
void SelectDWTBlock(int scale, int size_x,
const string &quadrant,
int &x1, int &x2) {   
   double Nx = log10(size_x) / log10(2.0);
   x1=DWT_Imin(scale,Nx,quadrant[0]);
   x2=DWT_Imax(scale,Nx,quadrant[0]);
}

void SelectDWTBlock(int scale, int size_x, int size_y,
const string &quadrant,
int &x1, int &x2, int &y1, int &y2) {   
   double Nx = log10(size_x) / log10(2.0);
   double Ny = log10(size_y) / log10(2.0);
   x1=DWT_Imin(scale,Nx,quadrant[0]);
   y1=DWT_Imin(scale,Ny,quadrant[1]);
   x2=DWT_Imax(scale,Nx,quadrant[0]);
   y2=DWT_Imax(scale,Ny,quadrant[1]);
}

void SelectDWTBlock(int scale,
int size_x, int size_y, int size_z,
const string &quadrant,
int &x1, int &x2, int &y1, int &y2, int &z1, int &z2) {   
   double Nx = log10(size_x) / log10(2.0);
   double Ny = log10(size_y) / log10(2.0);
   double Nz = log10(size_z) / log10(2.0);
   x1=DWT_Imin(scale,Nx,quadrant[0]);
   y1=DWT_Imin(scale,Ny,quadrant[1]);
   z1=DWT_Imin(scale,Nz,quadrant[2]);
   x2=DWT_Imax(scale,Nx,quadrant[0]);
   y2=DWT_Imax(scale,Ny,quadrant[1]);
   z2=DWT_Imax(scale,Nz,quadrant[2]);
}

// Provide block -----------------------------------------------------------
#define DWT_Scale(i,smax) ((int)((i==0)?smax-1:(ABS((CEIL(log10(i+1)/log10(2.0))-smax)))))
#define DWT_Quadrant1D(i,s,smax) ((s!=smax-1)?'1':((i==0)?'0':'1'))
#define DWT_QuadrantnD(i,s,sp,smax) \
  ((s!=sp)?'0':DWT_Quadrant1D(i,s,smax))
   
void Get_Scale_Quadrant(int size_x, int x,
   int &scale, string &quadrant) {
   double Nx = log10(size_x) / log10(2.0);
   quadrant="x";
   scale=DWT_Scale(x,Nx);
   quadrant[0]=DWT_Quadrant1D(x,scale,Nx);
}

void Get_Scale_Quadrant(int size_x, int size_y, int x, int y,
   int &scale, string &quadrant) {
   double Nx = log10(size_x) / log10(2.0);
   double Ny = log10(size_y) / log10(2.0);
   quadrant="xy";
   double scalex=DWT_Scale(x,Nx);
   double scaley=DWT_Scale(y,Ny);
   scale = (int)(MIN(scalex,scaley));
   quadrant[1]=DWT_QuadrantnD(y,scaley,scale,Ny);
   quadrant[0]=DWT_QuadrantnD(x,scalex,scale,Nx);
}

void Get_Scale_Quadrant(int size_x, int size_y, int size_z,
   int x, int y, int z,
   int &scale, string &quadrant) {
   double Nx = log10(size_x) / log10(2.0);
   double Ny = log10(size_y) / log10(2.0);
   double Nz = log10(size_z) / log10(2.0);
   quadrant="xyz";
   double scalex=DWT_Scale(x,Nx);
   double scaley=DWT_Scale(y,Ny);
   double scalez=DWT_Scale(z,Nz);
   scale = (int)(MIN(scalez,MIN(scalex,scaley)));
   quadrant[2]=DWT_QuadrantnD(z,scalez,scale,Nz);
   quadrant[1]=DWT_QuadrantnD(y,scaley,scale,Ny);
   quadrant[0]=DWT_QuadrantnD(x,scalex,scale,Nx);
}

// Clean quadrant ----------------------------------------------------------
void clean_quadrant(matrix2D<double> &I, int scale, const string &quadrant) {
   int x1, y1, x2, y2;
   SelectDWTBlock(scale, XSIZE(I), YSIZE(I), quadrant, x1, x2, y1, y2);
   for (int y=y1; y<=y2; y++)
      for (int x=x1; x<=x2; x++)
         DIRECT_MAT_ELEM(I,y,x)=0;
}

void clean_quadrant(matrix3D<double> &I, int scale, const string &quadrant) {
   int x1, y1, z1, x2, y2, z2;
   SelectDWTBlock(scale, XSIZE(I), YSIZE(I), ZSIZE(I), quadrant,
      x1, x2, y1, y2, z1, z2);
   for (int z=z1; z<=z2; z++)
      for (int y=y1; y<=y2; y++)
         for (int x=x1; x<=x2; x++)
            DIRECT_VOL_ELEM(I,z,y,x)=0;
}

// Soft thresholding -------------------------------------------------------
void soft_thresholding(matrix2D<double> &I, double th) {
   FOR_ALL_ELEMENTS_IN_MATRIX2D(I) 
      if (ABS(I(i,j))>th) 
         if (I(i,j)>0) I(i,j)-=th; else I(i,j)+=th;
      else I(i,j)=0;
}

void soft_thresholding(matrix3D<double> &I, double th) {
   FOR_ALL_ELEMENTS_IN_MATRIX3D(I) 
      if (ABS(I(k,i,j))>th) 
         if (I(k,i,j)>0) I(k,i,j)-=th; else I(k,i,j)+=th;
      else I(k,i,j)=0;
}
