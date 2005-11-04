/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Antonio Jose Rodriguez Sanchez (ajr@cnb.uam.es)
 *              Arun Kulshreshth        (arun_2000_iitd@yahoo.com)
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
#ifndef _XMIPPWAVELETS_HH
   #define _XMIPPWAVELETS_HH

#include "xmippMatrices3D.hh"
#include "Src/NumericalRecipes.hh"

/**@name Wavelets */
//@{

#define DAUB4   4
#define DAUB12 12
#define DAUB20 20

/**@name Bilib Wavelets */
//@{
/** B-spline Wavelet transform of a vector.
    The B-spline wavelet transform of the input array is computed.
    The size of the array must be so as to allow the downsampling by 2
    as many times as the number of iterations. For instance,
    if iterations is 1, then it must be a multiple of 2. If iterations is 2,
    then it must be a multiple of 4. If itereations if 3, then it must be
    a multiple of 8. And so on.
    
    If the isign=-1 then the inverse wavelet transform is perdormed*/
void Bilib_DWT(const matrix1D<double> &input,
   matrix1D<double> &result, int iterations, int isign=1);

/** B-spline Wavelet transform of a matrix.*/
void Bilib_DWT(const matrix2D<double> &input,
   matrix2D<double> &result, int iterations, int isign=1);

/** B-spline Wavelet transform of a matrix.*/
void Bilib_DWT(const matrix3D<double> &input,
   matrix3D<double> &result, int iterations, int isign=1);
//@}

/**@name Numerical recipes wavelets */
//@{
/** Set DWT type.
    The DWT type should be set before starting making transforms.
    Valid types are: DAUB4, DAUB12, DAUB20*/
    void set_DWT_type(int DWT_type);

/** DWT of a vector.
    The output vector can be the same as the input one. Previously the
    type of DWT must be set with set_DWT_type. If isign=1 the direct DWT
    is performed, if isign=-1 the inverse DWT is done*/
template <class T>
   void DWT(const matrix1D<T> &v, matrix1D<double> &result, int isign=1) {
      unsigned long int nn[1];
      unsigned long int *ptr_nn=nn-1;

      type_cast(v, result);
      nn[0] = XSIZE(result);
      double *ptr_result=MULTIDIM_ARRAY(result)-1;
      wtn(ptr_result, ptr_nn, 1, isign, pwt);
   }

/** DWT of a array.
    The output array can be the same as the input one. Previously the
    type of DWT must be set with set_DWT_type. If isign=1 the direct DWT
    is performed, if isign=-1 the inverse DWT is done*/
template <class T>
   void DWT(const matrix2D<T> &v, matrix2D<double> &result, int isign=1) {
      unsigned long int nn[2];
      unsigned long int *ptr_nn=nn-1;

      type_cast(v, result);
      nn[1] = YSIZE(result);
      nn[0] = XSIZE(result);
      double *ptr_result=MULTIDIM_ARRAY(result)-1;
      wtn(ptr_result, ptr_nn, 2, isign, pwt);
   }

/** DWT of a volume.
    The output vector can be the same as the input one. Previously the
    type of DWT must be set with set_DWT_type. If isign=1 the direct DWT
    is performed, if isign=-1 the inverse DWT is done*/
template <class T>
   void DWT(const matrix3D<T> &v, matrix3D<double> &result, int isign=1) {
      unsigned long int nn[2];
      unsigned long int *ptr_nn=nn-1;

      type_cast(v, result);
      nn[2] = ZSIZE(result);
      nn[1] = YSIZE(result);
      nn[0] = XSIZE(result);
      double *ptr_result=MULTIDIM_ARRAY(result)-1;
      wtn(ptr_result, ptr_nn, 3, isign, pwt);
   }

/** IDWT of a vector.
    The output vector can be the
    same as the input one. Previously the
    type of DWT must be set with set_DWT_type. */
   void IDWT(const matrix1D<double> &v,
   matrix1D<double> &result);

/** IDWT of an array.
    The output vector can be the
    same as the input one. Previously the
    type of DWT must be set with set_DWT_type.  */
   void IDWT(const matrix2D<double> &v,
   matrix2D<double> &result);

/** IDWT of a volume.
    The output volume can be the
    same as the input one. Previously the
    type of DWT must be set with set_DWT_type. */
   void IDWT(const matrix3D<double> &v,
   matrix3D<double> &result);

/** DWT Low pass versions.
    This function returns the low pass versions at different scales.
    The low pass version of the image at scale s is stored in the 01
    quadrant of that scale.*/
    void DWT_lowpass(const matrix2D<double> &v,
       matrix2D<double> &result);
//@}

#define DWT_Imin(s,smax,l) (int)((l=='0')?0:pow(2.0,smax-s-1))
#define DWT_Imax(s,smax,l) (int)((l=='0')?pow(2.0,smax-s-1)-1:pow(2.0,smax-s)-1)

/** Select Block 1D.
    Given the scale (s=0 is the finest) and the quadrant
    "0"(Lower frequencies) or "1"(Higher frequencies)
    this routine returns the indices that should be explored
    for this block (x1 and x2 should be included in the for). */
    template <class T>
    void SelectDWTBlock(int scale, const matrix1D<T> &I,
       const string &quadrant,
       int &x1, int &x2) {
         double Nx = Get_Max_Scale(XSIZE(I));
         I.physical2logical(DWT_Imin(scale,Nx,quadrant[0]),x1);
         I.physical2logical(DWT_Imax(scale,Nx,quadrant[0]),x2);
    }

/** Select Block 2D.
    Given the scale (s=0 is the finest) and the quadrant
    "xy"="00" (Upper left), "01" (Upper right), "10" (Lower left),
    "11" (Lower right).
    this routine returns the indices that should be explored
    for this block (the extremes should be included in the for). */
    template <class T>
    void SelectDWTBlock(int scale, const matrix2D<T> &I,
       const string &quadrant,
       int &x1, int &x2, int &y1, int &y2) {
         double Nx = Get_Max_Scale(XSIZE(I));
         double Ny = Get_Max_Scale(YSIZE(I));
         x1=DWT_Imin(scale,Nx,quadrant[0]);
         y1=DWT_Imin(scale,Ny,quadrant[1]);
         x2=DWT_Imax(scale,Nx,quadrant[0]);
         y2=DWT_Imax(scale,Ny,quadrant[1]);
         I.physical2logical(y1,x1,y1,x1);
         I.physical2logical(y2,x2,y2,x2);
    }

/** Select Block 3D.
    Given the scale (s=0 is the finest) and the quadrant
    "xyz"="000", "001", "010", "011", "100", "101", "110", "111".
    this routine returns the indices that should be explored
    for this block (the extremes should be included in the for). */
    template <class T>
    void SelectDWTBlock(int scale, const matrix3D<T> &I,
       const string &quadrant,
       int &x1, int &x2, int &y1, int &y2, int &z1, int &z2) {
         double Nx = Get_Max_Scale(XSIZE(I));
         double Ny = Get_Max_Scale(YSIZE(I));
         double Nz = Get_Max_Scale(ZSIZE(I));
         x1=DWT_Imin(scale,Nx,quadrant[0]);
         y1=DWT_Imin(scale,Ny,quadrant[1]);
         z1=DWT_Imin(scale,Nz,quadrant[2]);
         x2=DWT_Imax(scale,Nx,quadrant[0]);
         y2=DWT_Imax(scale,Ny,quadrant[1]);
         z2=DWT_Imax(scale,Nz,quadrant[2]);
         I.physical2logical(z1,y1,x1,z1,y1,x1);
         I.physical2logical(z2,y2,x2,z2,y2,x2);
    }

   /** Get maximum scale.
       This function returns the maximum scale achievable by the 
       DWT transform of a given size.*/
   inline int Get_Max_Scale(int size) {return ROUND(log10((double)size) / log10(2.0));}

   /** Given a quadrant number it returns the string associated to it.
       That is nothing more than its corresponding binary representation */
   string Quadrant2D(int q);

   /** Given a quadrant number it returns the string associated to it.
       That is nothing more than its corresponding binary representation */
   string Quadrant3D(int q);

   /** Get scale and quadrant 1D.
       Given a point and the maximum size of the image,
       this routine returns the scale and quadrant it belongs. */
   void Get_Scale_Quadrant(int size_x, int x,
      int &scale, string &quadrant);

   /** Get scale and quadrant 2D.
       Given a point and the maximum size of the image,
       this routine returns the scale and quadrant it belongs. */
   void Get_Scale_Quadrant(int size_x, int size_y, int x, int y,
      int &scale, string &quadrant);

   /** Get scale and quadrant 3D.
       Given a point and the maximum size of the image,
       this routine returns the scale and quadrant it belongs. */
   void Get_Scale_Quadrant(int size_x, int size_y, int size_z,
      int x, int y, int z,
      int &scale, string &quadrant);

   /**@name Denoising */
   //@{
   /** Remove all information within a quadrant and scale */
   void clean_quadrant(matrix2D<double> &I, int scale, const string &quadrant);

   /** Remove all information within a quadrant and scale */
   void clean_quadrant(matrix3D<double> &I, int scale, const string &quadrant);

   /** Soft thresholding 2D.
       Substract a value from all coefficients, if the the value is
       greater than the absolute value of the coefficient, that coefficient
       is set to 0. */
   void soft_thresholding(matrix2D<double> &I, double th);

   /** Soft thresholding 3D.
       Substract a value from all coefficients, if the the value is
       greater than the absolute value of the coefficient, that coefficient
       is set to 0. */
   void soft_thresholding(matrix3D<double> &I, double th);
   
   /** Adaptive soft thresholding 2D.
       Chang, Yu, Betterli. IEEE Int. Conf. Image Processing*/
   void adaptive_soft_thresholding(matrix2D<double> &I, int scale);

   /** Keep central part 2D.
       Keep those coefficients in a certain radius. */
   void DWT_keep_central_part(matrix2D<double> &I, double R);

   /** Bayesian, Wiener filtering.
       Bijaoui, Signal Processing 2002, 82: 709-712. The denoising procedure
       is applied up to the scale given.
       SNR0 is the smallest SNR and SNRF is the largest SNR.
       
       This function returns the estimated coefficients for S and N at
       each scale. If denoise is set to false, then S and N coefficients
       are estimated but they are not applied to the image.*/
   matrix1D<double> bayesian_wiener_filtering(matrix2D<double> &WI, int allowed_scale,
      double SNR0=0.1, double SNRF=0.2, bool white_noise=false, int tell=0,
      bool denoise=true);

   /** Bayesian, Wiener filtering.
       This is the function that really denoise. */
   void bayesian_wiener_filtering(matrix2D<double> &WI,
      int allowed_scale, matrix1D<double> &estimatedS);

   /** Bayesian, Wiener filtering.
       Bijaoui, Signal Processing 2002, 82: 709-712. The denoising procedure
       is applied up to the scale given.
       SNR0 is the smallest SNR and SNRF is the largest SNR.
       
       This function returns the estimated coefficients for S and N at
       each scale. If denoise is set to false, then S and N coefficients
       are estimated but they are not applied to the image.*/
   matrix1D<double> bayesian_wiener_filtering(matrix3D<double> &WI, int allowed_scale,
      double SNR0=0.1, double SNRF=0.2, bool white_noise=false, int tell=0,
      bool denoise=true);

   /** Bayesian, Wiener filtering.
       This is the function that really denoise. */
   void bayesian_wiener_filtering(matrix3D<double> &WI,
      int allowed_scale, matrix1D<double> &estimatedS);
   //@}
//@}
#endif
