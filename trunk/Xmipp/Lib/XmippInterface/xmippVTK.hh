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
/*****************************************************************************/
/* INTERACTION WITH VTK                                                      */
/*****************************************************************************/

#ifndef _XMIPP_VTK_HH
   #define _XMIPP_VTK_HH

#ifdef _HAVE_VTK
   #include <XmippData/xmippMatrices3D.hh>
   #include <XmippData/xmippImages.hh>
   #include <XmippData/xmippVolumes.hh>
   #include <vtkStructuredPoints.h>
   #include <vtkImageData.h>
/**@name VTK */
//@{
/**@name Xmipp <--> VTK */
//@{
   /** Xmipp Vector --> VTK.
       scalarN is the number of scalars in the VTK array. The xmipp Array only
       has 1 value so the difference with scalarN is filled with 0's*/
   template <class T, class VTKT>
      void xmippArray2VTK(const matrix1D<T> &v, VTKT * &retval, int scalarN=1);
   
   /** Xmipp Matrix --> VTK */
   template <class T, class VTKT>
      void xmippArray2VTK(const matrix2D<T> &v, VTKT * &retval, int scalarN=1);
   
   /** Xmipp Volume --> VTK */
   template <class T, class VTKT>
      void xmippArray2VTK(const matrix3D<T> &v, VTKT * &retval, int scalarN=1);

   /** Shift and wrap a 1D image along direction X or Y.
       Used for adjusting the volume from/for VTK
       when the size is even. Howver this function always shift, no
       matter the size*/
   template <class T>
      void shift_for_VTK(matrix1D<T> &v);

   /** Shift and wrap a 2D image along direction X or Y.
       Used for adjusting the volume from/for VTK
       when the size is even. Howver this function always shift, no
       matter the size. Valid directions are 'x' or 'y'*/
   template <class T>
      void shift_for_VTK(matrix2D<T> &v, char dir);

   /** Shift and wrap a 3D image along direction X or Y.
       Used for adjusting the volume from/for VTK
       when the size is even. Howver this function always shift, no
       matter the size. Valid directions are 'x' or 'y'*/
   template <class T>
      void shift_for_VTK(matrix3D<T> &v, char dir);

   /** FFT_Xmipp -> VTK 1D.
       Converts a Xmipp FFT into a vtkImageData with an FFT.
       An exception is thrown if the input array does not have an even
       dimension on X.*/
   void xmippFFT2VTK(matrix1D <complex <double > > &v, vtkImageData * &retval) _THROW;
   
   /** FFT_Xmipp -> VTK 2D.
       Converts a Xmipp FFT into a vtkImageData with an FFT.
       An exception is thrown if the input array does not have an even
       dimension on X.*/
   void xmippFFT2VTK(FourierImageXmipp &v, vtkImageData * &retval) _THROW;

   /** FFT_Xmipp -> VTK 3D.
       Converts a Xmipp FFT into a vtkImageData with an FFT.
       An exception is thrown if the input array does not have an even
       dimension on X.*/
   void xmippFFT2VTK(FourierVolumeXmipp &v, vtkImageData * &retval) _THROW;

   /** VTK -> FFT_Xmipp.
       Converts a vtkImageData with an FFT to an Xmipp FFT */
   void VTK2xmippFFT(vtkImageData *v, matrix1D< complex <double > > &retval);

   /** VTK -> FFT_Xmipp.
       Converts a vtkImageData with an FFT to an Xmipp FFT */
   void VTK2xmippFFT(vtkImageData *v, FourierImageXmipp &retval);

   /** VTK -> FFT_Xmipp.
       Converts a vtkImageData with an FFT to an Xmipp FFT */
   void VTK2xmippFFT(vtkImageData *v, FourierVolumeXmipp &retval);

   /** Resize a XmippVector after a VTK object.
       An exception is thrown if the dimensionality of the VTK object
       does not fit into the target Xmipp type.*/
   template <class T, class VTKT>
      void xmippArray_resize_VTK(matrix1D<T> &retval, VTKT *v) _THROW;

   /** Resize a XmippMatrix after a VTK object.
       An exception is thrown if the dimensionality of the VTK object
       does not fit into the target Xmipp type.*/
   template <class T, class VTKT>
      void xmippArray_resize_VTK(matrix2D<T> &retval, VTKT *v) _THROW;

   /** Resize a XmippVolume after a VTK object.
       An exception is thrown if the dimensionality of the VTK object
       does not fit into the target Xmipp type.*/
   template <class T, class VTKT>
      void xmippArray_resize_VTK(matrix3D<T> &retval, VTKT *v) _THROW;

   /** VTK --> XmippVector */
   template <class T, class VTKT>
      void VTK2xmippArray(VTKT *v, matrix1D<T> &retval);

   /** VTK --> XmippMatrix */
   template <class T, class VTKT>
      void VTK2xmippArray(VTKT *v, matrix2D<T> &retval);

   /** VTK --> XmippVolume */
   template <class T, class VTKT>
      void VTK2xmippArray(VTKT *v, matrix3D<T> &retval);

   /** "Copy constructor" for VTKImageData.
       The output array can be NULL, or full of data before operation.*/
   template <class VTKT>
      void VTK2VTK(VTKT *v_in, VTKT *&v_out, bool change_type=FALSE,
         int new_type=VTK_FLOAT);

   /** Same shape.
       TRUE if the two VTK objects have the same dimensions and sizes */
   bool same_shape(vtkImageData *v1, vtkImageData *v2);

   /** Center an FFT. */
   void CenterFFT(vtkImageData *&v);

   /** Show VTK object. */
   template <class VTKT>
      void VTK_print (ostream &out, VTKT *v);

   /** Index to frequency.
       This function can be used with vectors of any size (1,2,3).
       The Digital spectrum is limited between -1/2 and 1/2.
       If the vector has got more than 3 coordinates, then an exception
       is thrown*/
   inline void FFT_idx2digfreq(vtkImageData *fft, const matrix1D<int> &idx,
      matrix1D<double> &freq) _THROW {
         if (XSIZE(idx)<1 || XSIZE(idx)>3)
            REPORT_ERROR(1,"FFT_idx2digfreq: Index is not of the correct size");
         freq.resize(XSIZE(idx));

         int *wholeExtent = fft->GetWholeExtent();
	 
         double size[3];
         size[0] = (double)(wholeExtent[0] + wholeExtent[1] + 1);
         size[1] = (double)(wholeExtent[2] + wholeExtent[3] + 1);
         size[2] = (double)(wholeExtent[4] + wholeExtent[5] + 1);

         FOR_ALL_ELEMENTS_IN_MATRIX1D(idx) {
            VEC_ELEM(freq,i)=(VEC_ELEM(idx,i)<size[i]/2)?
               VEC_ELEM(idx,i):-size[i]+VEC_ELEM(idx,i);
            VEC_ELEM(freq,i) /= size[i];
         }
   }
   
   /** Frequency to index.
       This function can be used with vectors of any size (1,2,3).
       The Digital spectrum is limited between -1/2 and 1/2.
       If the vector has got more than 3 coordinates, then an exception
       is thrown*/
   inline void digfreq2FFT_idx(vtkImageData *fft, const matrix1D<double> &freq,
      matrix1D<int> &idx) _THROW {
         if (XSIZE(freq)<1 || XSIZE(freq)>3)
            REPORT_ERROR(1,"digfreq2FFT_idx: freq is not of the correct size");
         idx.resize(XSIZE(freq));

         int *wholeExtent = fft->GetWholeExtent();
	 
         double size[3];
         size[0] = (double)(wholeExtent[0] + wholeExtent[1] + 1);
         size[1] = (double)(wholeExtent[2] + wholeExtent[3] + 1);
         size[2] = (double)(wholeExtent[4] + wholeExtent[5] + 1);

         FOR_ALL_ELEMENTS_IN_MATRIX1D(idx) {
	    if (VEC_ELEM(freq,i)>=0)
	       VEC_ELEM(idx,i)=(int)(size[i]*VEC_ELEM(freq,i));
	    else
	       VEC_ELEM(idx,i)=(int)(size[i]+size[i]*VEC_ELEM(freq,i));
         }
   }
   
   /** Digital to Continuous frequency.
       The pixel size must be given in Amstrongs. The digital frequency is
       between [-1/2,1/2].*/
   inline void digfreq2contfreq(const matrix1D<double> &digfreq,
      matrix1D<double> &contfreq, double pixel_size)
      {contfreq=digfreq/pixel_size;}
   
   /** Continuous to Digital frequency.
       The pixel size must be given in Amstrongs. The digital frequency is
       between [-1/2,1/2].*/
   inline void contfreq2digfreq(const matrix1D<double> &contfreq,
      matrix1D<double> &digfreq, double pixel_size)
      {digfreq=contfreq*pixel_size;}
//@}

/**@name FFT */
//@{
   /** FFT of a VTK vector/image/volume.
       N=no. samples in FFT, if it is 0 then the same as the ones in real
       space are used. By default, the Fourier transform is centered.*/
   void FFT_VTK(vtkImageData *v_in, vtkImageData *&fft_out,
      bool do_not_center=FALSE);

   /** FFT of an XmippVector. */
   template <class T>
      void FFT_VTK(const matrix1D<T> &v, vtkImageData *&fft_out,
      bool do_not_center=FALSE)
         {vtkImageData *vtkI=NULL; xmippArray2VTK(v, vtkI);
          FFT_VTK(vtkI,fft_out,do_not_center); vtkI->Delete();}

   /** FFT of an XmippMatrix. */
   template <class T>
      void FFT_VTK(const matrix2D<T> &v, vtkImageData *&fft_out,
      bool do_not_center=FALSE)
         {vtkImageData *vtkI=NULL; xmippArray2VTK(v, vtkI);
          FFT_VTK(vtkI,fft_out,do_not_center); vtkI->Delete();}

   /** FFT of an XmippVolume. */
   template <class T>
      void FFT_VTK(const matrix3D<T> &v, vtkImageData * &fft_out,
      bool do_not_center=FALSE)
         {vtkImageData *vtkI=NULL; xmippArray2VTK(v, vtkI);
          FFT_VTK(vtkI,fft_out,do_not_center); vtkI->Delete();}

   /** Magnitude of FFT. Valid for 2D and 3D */
   template <class maT>
      void FFT_magnitude(vtkImageData *fft_in, maT &mag);

   /** Phase of FFT.
       An exception is thrown if the operation is not valid. */
   template <class maT>
      void FFT_phase(vtkImageData *fft_in, maT &phase) _THROW;

   /** Inverse FFT of a vector/image/volume */
   void IFFT_VTK(vtkImageData *fft_in, vtkImageData *&v_out, bool
      is_not_centered=FALSE);

   /** IFFT of an XmippVector. */
   template <class T>
      void IFFT_VTK(vtkImageData * fft_in, matrix1D<T> &v,
         bool is_not_centered=FALSE)
         {vtkImageData *vtkI=NULL; IFFT_VTK(fft_in,vtkI,is_not_centered);
          VTK2xmippArray(vtkI, v); vtkI->Delete();}

   /** IFFT of an XmippMatrix. */
   template <class T>
      void IFFT_VTK(vtkImageData * fft_in, matrix2D<T> &v,
         bool is_not_centered=FALSE)
         {vtkImageData *vtkI=NULL; IFFT_VTK(fft_in,vtkI,is_not_centered);
          VTK2xmippArray(vtkI, v); vtkI->Delete();}

   /** IFFT of an XmippVolume.
       Although it is prepared for doing so, you should not use centered
       FFTs to produce the real image, in general, the result is this is
       done looks like the original picture but most values are lost. */
   template <class T>
      void IFFT_VTK(vtkImageData * fft_in, matrix3D<T> &v,
         bool is_not_centered=FALSE)
         {vtkImageData *vtkI=NULL; IFFT_VTK(fft_in,vtkI,is_not_centered);
          VTK2xmippArray(vtkI, v); vtkI->Delete();}

   /** Autocorrelation function of an Xmipp matrix.
       Fast calcuation of the autocorrelation matrix of a given one using
	   Fast Fourier Transform. (Using the correlation theorem) */
  template <class T>
   void auto_correlation_matrix(matrix2D<T> const &Img,matrix2D<double> &R);

   /** Autocorrelation function of an Xmipp matrix.
       Fast calcuation of the correlation matrix on two matrices using
	   Fast Fourier Transform. (Using the correlation theorem) */
  template <class T>
   void correlation_matrix(matrix2D<T> const &m1,matrix2D<T> const &m2,
                           matrix2D<double> &R);


   /** Series convolution function. Gives the convolution of two series
   	   given as Xmipp Vectors. Result is stored in result vector.
       Fast calcuation of the convolution result using
	   Fast Fourier Transform. If FullConvolution, set by default to
	   FALSE, is TRUE the full convolution series is returned. Otherwise
	   the convolution vector refers only to the valid values, whose
	   number is the greater dimension of the two series.
	   Note: Complex numbers are allowed */
template <class T>
void series_convolution(matrix1D<T> &series1,matrix1D<T> &series2,
                        matrix1D<T> &result,bool FullConvolution=FALSE);

   /** Convolution of a series (vector) with a given filter. Result is stored
       in result vector. Use this function when the dimension of the result
	   must be the dimension of the given series (first argument)
	   Note: Complex numbers are allowed */
template <class T>
void convolve(matrix1D<T> &series1,matrix1D<T> &filter, matrix1D<T> &result);

   /** numerical_derivative : This function computes the numerical derivative
      of a matrix in Y direction (rows) or X direction (columns) of a given
	   matrix, using a Savitzky-Golay filter on every row or column, and then
	   convolving.Input matrix is M, result is stored in D, direction can have
	   values of 'x' or 'y'. Order is the derivative order.
	   Window size and polynomial order are parameters for the
	   Savitzky-Golay filter that define the number of points forward (+window_size)
	   and backward (-window_size) considered to calculate the filter, and the degree of
	   the polynomial to interpolate these values, respectively. Default values are
	   window_size=2 and polynomial_order=4, which are equivalent to a
	   5-point algorithm to calculate the derivate and give good results.
	   But they can be changed provided that polynomial_order <= 2*window_size.
	   As one can expect, the values of the matrix in a border of size
	   window_size are not accurate ones, as there aren't enough points to perform
	   the estimation. As a rule of thumb, the greater the window, the more the
	   filtering and the less the precission of the derivatives, and the
	   greater the order of the polynomial, the greater the precission.  */

void numerical_derivative(matrix2D<double> &M,matrix2D<double> &D,
							char direction,int order,
							int window_size=2,int polynomial_order=4);



   /** Some variables used in speed up macros related with VTK */
   #define SPEED_UP_vtk \
      int *inExt, maxX, maxY, maxZ; \
      float *vtkPtr;

   /** Macro to apply to all elements in a VTKImageData.
       You must be careful of moving the pointer (vtkPtr) appropiately
       with the number of scalars. The VTK object is supposed to have
       VTK_FLOATs. Indexes k,i,j are defined corresponding to Z,Y,X.
       The SPEED_UP_vtk variables are needed.*/
   #define FOR_ALL_ELEMENTS_IN_VTK(v) \
      inExt=v->GetExtent(); \
      maxX = inExt[1] - inExt[0]; \
      maxY = inExt[3] - inExt[2]; \
      maxZ = inExt[5] - inExt[4]; \
      vtkPtr = (float *) v->GetScalarPointer(); \
      for (int k = 0; k <= maxZ; k++) \
          for (int i = 0; i <= maxY; i++) \
              for (int j = 0; j <= maxX; j++) \


//@}
//@}
#endif

#endif
