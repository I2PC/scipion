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
   #include <XmippData/xmippFFT.hh>
   #include <vtkStructuredPoints.h>
   #include <vtkImageData.h>
/**@name VTK */
//@{
/* Xmipp --> VTK ----------------------------------------------------------- */
#define SET_VTK_TYPE(T,retval) \
   if      (typeid(T)==typeid(unsigned char)) \
      retval->SetScalarType(VTK_UNSIGNED_CHAR); \
   else if (typeid(T)==typeid(short int)) \
      retval->SetScalarType(VTK_SHORT); \
   else if (typeid(T)==typeid(int)) \
      retval->SetScalarType(VTK_INT); \
   else if (typeid(T)==typeid(float) || typeid(T)==typeid(double)) \
      retval->SetScalarType(VTK_FLOAT);

#define COPY_VALUES_TO_VTK(T,v,retval, scalarN) \
   retval->AllocateScalars(); \
   if (typeid(T)!=typeid(double)) { \
      T *ptr= (T *) retval->GetScalarPointer(); \
      for (int i=0; i<MULTIDIM_SIZE(v); i++) \
      { \
         *ptr++=MULTIDIM_ELEM(v,i); \
	 for (int n=1; n<scalarN; n++) *ptr++=0; \
      } \
   } \
   else \
   { \
      float *ptr= (float *) retval->GetScalarPointer(); \
      for (int i=0; i<MULTIDIM_SIZE(v); i++) { \
         *ptr++=(float)MULTIDIM_ELEM(v,i); \
	 for (int n=1; n<scalarN; n++) *ptr++=0.0f; \
      } \
   }

/**@name Xmipp <--> VTK */
//@{
   /** Xmipp Vector --> VTK.
       scalarN is the number of scalars in the VTK array. The xmipp Array only
       has 1 value so the difference with scalarN is filled with 0's*/
   template <class T, class VTKT>
      void xmippArray2VTK(const matrix1D<T> &v, VTKT * &retval, int scalarN=1) {
         if (retval==NULL) retval=VTKT::New();
         else retval->Initialize();
         retval->SetDimensions(XSIZE(v),1,1);
         retval->SetNumberOfScalarComponents(scalarN);
         SET_VTK_TYPE(T,retval);
         COPY_VALUES_TO_VTK(T,v,retval,scalarN);
         retval->Update();
   }
   
   /** Xmipp Matrix --> VTK */
   template <class T, class VTKT>
      void xmippArray2VTK(const matrix2D<T> &v, VTKT * &retval, int scalarN=1)  {
         if (retval==NULL) retval=VTKT::New();
         else retval->Initialize();
         retval->SetDimensions(XSIZE(v),YSIZE(v),1);
         retval->SetNumberOfScalarComponents(scalarN);
         SET_VTK_TYPE(T,retval);
         COPY_VALUES_TO_VTK(T,v,retval,scalarN);
         retval->Update();
      }
   
   /** Xmipp Volume --> VTK */
   template <class T, class VTKT>
      void xmippArray2VTK(const matrix3D<T> &v, VTKT * &retval, int scalarN=1) {
         if (retval==NULL) retval=VTKT::New();
         else retval->Initialize();
         retval->SetDimensions(XSIZE(v),YSIZE(v),ZSIZE(v));
         retval->SetNumberOfScalarComponents(scalarN);
         SET_VTK_TYPE(T,retval);
         COPY_VALUES_TO_VTK(T,v,retval,scalarN);
         retval->Update();
      }

   /** Shift and wrap a 1D image along direction X or Y.
       Used for adjusting the volume from/for VTK
       when the size is even. However this function always shift, no
       matter the size*/
   template <class T>
      void shift_for_VTK(matrix1D<T> &v) {
         int firstX=STARTINGX(v);
         int finalX=FINISHINGX(v);
         T aux=v(firstX);
         for (int j=firstX; j<finalX; j++)
             v(j)=v(j+1);
         v(finalX)=(T)aux;
      }

   /** Shift and wrap a 2D image along direction X or Y.
       Used for adjusting the volume from/for VTK
       when the size is even. However this function always shift, no
       matter the size. Valid directions are 'x' or 'y'*/
   template <class T>
      void shift_for_VTK(matrix2D<T> &v, char dir) {
         int firstY=STARTINGY(v);
         int firstX=STARTINGX(v);
         int finalY=FINISHINGY(v);
         int finalX=FINISHINGX(v);
         if (dir=='y') {
            for (int j=firstX; j<=finalX; j++) {
   	        T aux=v(firstY,j);
   	        for (int i=firstY; i<finalY; i++)
   	           v(i,j)=v(i+1,j);
   	        v(finalY,j)=(T)aux;
            }
         }
         if (dir=='x') {
            for (int i=firstY; i<=finalY; i++) {
   	        T aux=v(i,firstX);
   	        for (int j=firstX; j<finalX; j++)
   	           v(i,j)=v(i,j+1);
   	        v(i,finalX)=(T)aux;
            }
         }
      }

   /** Shift and wrap a 3D image along direction X or Y.
       Used for adjusting the volume from/for VTK
       when the size is even. However this function always shift, no
       matter the size. Valid directions are 'x' or 'y'*/
   template <class T>
      void shift_for_VTK(matrix3D<T> &v, char dir) {
         int firstZ=STARTINGZ(v);
         int firstY=STARTINGY(v);
         int firstX=STARTINGX(v);
         int finalZ=FINISHINGZ(v);
         int finalY=FINISHINGY(v);
         int finalX=FINISHINGX(v);
         if (dir=='z') {
            for (int i=firstY; i<=finalY; i++)
	        for (int j=firstX; j<=finalX; j++) {
   	            T aux=v(firstZ,i,j);
   	            for (int k=firstZ; k<finalZ; k++)
   	    	       v(k,i,j)=v(k+1,i,j);
   	            v(finalZ,i,j)=(T)aux;
	        }
         }
         if (dir=='y') {
            for (int k=firstZ; k<=finalZ; k++)
	        for (int j=firstX; j<=finalX; j++) {
   	            T aux=v(k,firstY,j);
   	            for (int i=firstY; i<finalY; i++)
   	    	       v(k,i,j)=v(k,i+1,j);
   	            v(k,finalY,j)=(T)aux;
	        }
         }
         if (dir=='x') {
            for (int k=firstZ; k<=finalZ; k++)
	        for (int i=firstY; i<=finalY; i++) {
   	            T aux=v(k,i,firstX);
   	            for (int j=firstX; j<finalX; j++)
   	    	       v(k,i,j)=v(k,i,j+1);
   	            v(k,i,finalX)=(T)aux;
	        }
         }
      }

   /** FFT_Xmipp -> VTK 1D.
       Converts a Xmipp FFT into a vtkImageData with an FFT.
       An exception is thrown if the input array does not have an even
       dimension on X.*/
   void xmippFFT2VTK(matrix1D <complex <double > > &v, vtkImageData * &retval);
   
   /** FFT_Xmipp -> VTK 2D.
       Converts a Xmipp FFT into a vtkImageData with an FFT.
       An exception is thrown if the input array does not have an even
       dimension on X.*/
   void xmippFFT2VTK(FourierImageXmipp &v, vtkImageData * &retval);

   /** FFT_Xmipp -> VTK 3D.
       Converts a Xmipp FFT into a vtkImageData with an FFT.
       An exception is thrown if the input array does not have an even
       dimension on X.*/
   void xmippFFT2VTK(FourierVolumeXmipp &v, vtkImageData * &retval);

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
      void xmippArray_resize_VTK(matrix1D<T> &retval, VTKT *v) {
         int dim[3]; v->GetDimensions(dim);
         if (dim[1]!=1 && dim[2]!=1)
            REPORT_ERROR(1,"VTK2xmippVector: VTK image is not a vector");
         retval.resize(dim[0]);
      }

   /** Resize a XmippMatrix after a VTK object.
       An exception is thrown if the dimensionality of the VTK object
       does not fit into the target Xmipp type.*/
   template <class T, class VTKT>
      void xmippArray_resize_VTK(matrix2D<T> &retval, VTKT *v) {
         int dim[3]; v->GetDimensions(dim);
         if (dim[2]!=1)
            REPORT_ERROR(1,"VTK2xmippMatrix: VTK image is not a matrix");
         retval.resize(dim[1],dim[0]);
      }

   /** Resize a XmippVolume after a VTK object.
       An exception is thrown if the dimensionality of the VTK object
       does not fit into the target Xmipp type.*/
   template <class T, class VTKT>
      void xmippArray_resize_VTK(matrix3D<T> &retval, VTKT *v) {
         int dim[3]; v->GetDimensions(dim);
         retval.resize(dim[2],dim[1],dim[0]);
      }

   /** Resize a VTK object after another VTK object.*/
   template <class VTKT>
      void VTK_resize_VTK(VTKT *v_in, VTKT *&v_out, bool change_type=FALSE,
         int new_type=VTK_FLOAT) {
         if (v_out==NULL)  v_out=VTKT::New();
         else              v_out->PrepareForNewData();
         if (v_in!=NULL) {
            if (!change_type) v_out->SetScalarType(v_in->GetScalarType());
            else              v_out->SetScalarType(new_type);
            v_out->SetNumberOfScalarComponents(v_in->GetNumberOfScalarComponents());
            v_out->CopyStructure(v_in);
            v_out->Update();
         }
      }

#define COPY_VALUES_FROM_VTK(T,v,retval) \
   xmippArray_resize_VTK(retval,v); \
   int VTK_scalar_type=v->GetScalarType(); \
   unsigned char *uptr=NULL; \
   short int     *sptr=NULL; \
   int           *iptr=NULL; \
   float         *fptr=NULL; \
   switch (VTK_scalar_type) { \
      case VTK_UNSIGNED_CHAR: uptr=(unsigned char *) v->GetScalarPointer(); break; \
      case VTK_SHORT: sptr=(short int *) v->GetScalarPointer(); break; \
      case VTK_INT: iptr=(int *) v->GetScalarPointer(); break; \
      case VTK_FLOAT: fptr=(float *) v->GetScalarPointer(); break; \
   } \
   for (int i=0; i<MULTIDIM_SIZE(retval); i++) \
   {\
   	      switch (VTK_scalar_type) \
	      { \
               case VTK_UNSIGNED_CHAR: MULTIDIM_ELEM(retval,i)=(T)*uptr++; break; \
               case VTK_SHORT:         MULTIDIM_ELEM(retval,i)=(T)*sptr++; break; \
      	       case VTK_INT:           MULTIDIM_ELEM(retval,i)=(T)*iptr++; break; \
               case VTK_FLOAT:         MULTIDIM_ELEM(retval,i)=(T)*fptr++; break; \
          } \
   } 

   /** VTK --> XmippVector */
   template <class T, class VTKT>
      void VTK2xmippArray(VTKT *v, matrix1D<T> &retval) {
         retval.clear();
         if (v==NULL) return;
         COPY_VALUES_FROM_VTK(T,v,retval);
      }

   /** VTK --> XmippMatrix */
   template <class T, class VTKT>
      void VTK2xmippArray(VTKT *v, matrix2D<T> &retval) {
         retval.clear();
         if (v==NULL) return;
         COPY_VALUES_FROM_VTK(T,v,retval);
      }

   /** VTK --> XmippVolume */
   template <class T, class VTKT>
      void VTK2xmippArray(VTKT *v, matrix3D<T> &retval) {
         retval.clear();
         if (v==NULL) return;
         COPY_VALUES_FROM_VTK(T,v,retval);
      }

   /** "Copy constructor" for VTKImageData.
       The output array can be NULL, or full of data before operation.*/
   template <class VTKT>
      void VTK2VTK(VTKT *v_in, VTKT *&v_out, bool change_type=FALSE,
         int new_type=VTK_FLOAT) {
         if (v_in==NULL) return;
         VTK_resize_VTK(v_in,v_out,change_type,new_type);
         v_out->CopyAndCastFrom(v_in,v_in->GetExtent());
         v_out->Update();
      }

   /** Same shape.
       TRUE if the two VTK objects have the same dimensions and sizes */
   bool VTK_same_shape(vtkImageData *v1, vtkImageData *v2);

   /** Self center an FFT.
       Use the offset to uncenter properly odd size images*/
   void VTK_CenterFFT(vtkImageData *&v, int zoff=0, int yoff=0, int xoff=0);

   /** Center an FFT.
       Use the offset to uncenter properly odd size images*/
   void VTK_CenterFFT(vtkImageData *&v_in, vtkImageData *&v_out,
      int zoff=0, int yoff=0, int xoff=0);

   /** Show VTK object. */
   template <class VTKT>
      void VTK_print (ostream &out, VTKT *v) {
         if (v==NULL) {out << "NULL vtkdata\n"; return;}
         int maxC = v->GetNumberOfScalarComponents();

         // find the region to loop over
         int *inExt=v->GetExtent();
         int maxX = inExt[1] - inExt[0];
         int maxY = inExt[3] - inExt[2]; 
         int maxZ = inExt[5] - inExt[4];

         // Compute maximum value in the array
         float *inPtr  = (float *) v->GetScalarPointer();
         float max_val=ABS(*inPtr);
         for (int idxZ = 0; idxZ <= maxZ; idxZ++)
             for (int idxY = 0; idxY <= maxY; idxY++)
                 for (int idxX = 0; idxX <= maxX; idxX++)
                     for (int k=0; k<maxC; k++) {
                         max_val=MAX(max_val,ABS(*inPtr++));
                     }

         // Show
         int prec=best_prec(max_val,10);
         inPtr  = (float *) v->GetScalarPointer();
         for (int idxZ = 0; idxZ <= maxZ; idxZ++) {
             out << "Slice " << idxZ << " ----------\n";
             for (int idxY = 0; idxY <= maxY; idxY++) {
                 for (int idxX = 0; idxX <= maxX; idxX++) {
                     out << "(";
                     for (int k=0; k<maxC; k++) {
                        float f;
                        f=(ABS(*inPtr)>XMIPP_EQUAL_ACCURACY)? *inPtr:0;
                        if (k==maxC-1) out << FtoA(f,10,prec);
                        else           out << FtoA(f,10,prec) << ",";
                        inPtr++;
                     }
                     out << ") ";
                 }
                 out << endl;
             }
         }
      }

   /** Index to frequency.
       This function can be used with vectors of any size (1,2,3).
       The Digital spectrum is limited between -1/2 and 1/2.
       If the vector has got more than 3 coordinates, then an exception
       is thrown*/
   inline void VTK_FFT_idx2digfreq(vtkImageData *fft, const matrix1D<int> &idx,
      matrix1D<double> &freq) {
         if (XSIZE(idx)<1 || XSIZE(idx)>3)
            REPORT_ERROR(1,"FFT_idx2digfreq: Index is not of the correct size");
         freq.resize(XSIZE(idx));

         int *wholeExtent = fft->GetWholeExtent();
	 
         double size[3];
         size[0] = (double)(wholeExtent[0] + wholeExtent[1] + 1);
         size[1] = (double)(wholeExtent[2] + wholeExtent[3] + 1);
         size[2] = (double)(wholeExtent[4] + wholeExtent[5] + 1);

         FOR_ALL_ELEMENTS_IN_MATRIX1D(idx)
            FFT_IDX2DIGFREQ(VEC_ELEM(idx,i),size[i],VEC_ELEM(freq,i));
   }

   /** Frequency to index.
       This function can be used with vectors of any size (1,2,3).
       The Digital spectrum is limited between -1/2 and 1/2.
       If the vector has got more than 3 coordinates, then an exception
       is thrown*/
   inline void VTK_digfreq2FFT_idx(vtkImageData *fft, const matrix1D<double> &freq,
      matrix1D<int> &idx) {
         if (XSIZE(freq)<1 || XSIZE(freq)>3)
            REPORT_ERROR(1,"digfreq2FFT_idx: freq is not of the correct size");
         idx.resize(XSIZE(freq));

         int *wholeExtent = fft->GetWholeExtent();
	 
         double size[3];
         size[0] = (double)(wholeExtent[0] + wholeExtent[1] + 1);
         size[1] = (double)(wholeExtent[2] + wholeExtent[3] + 1);
         size[2] = (double)(wholeExtent[4] + wholeExtent[5] + 1);

         FOR_ALL_ELEMENTS_IN_MATRIX1D(idx)
            DIGFREQ2FFT_IDX(VEC_ELEM(freq,i),size[i],VEC_ELEM(idx,i));
   }
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
      void VTK_FFT_magnitude(vtkImageData *fft_in, maT &mag) {
            vtkImageMagnitude *mag_filter=vtkImageMagnitude::New();
            mag_filter->SetInput(fft_in); mag_filter->Update();
            VTK2xmippArray(mag_filter->GetOutput(),mag);
            mag_filter->Delete();
      }

   /** Magnitude of FFT. */
   void VTK_FFT_magnitude(FourierImageXmipp &fft_in,
      matrix2D<double> &mag, bool do_not_center=false);

   /** Phase of FFT.
       An exception is thrown if the operation is not valid. */
   template <class maT>
      void VTK_FFT_phase(vtkImageData *fft_in, maT &phase) {
         // Check validity of operation
         int maxC = fft_in->GetNumberOfScalarComponents();
         if (maxC!=2)
            REPORT_ERROR(1,"FFT_phase: Input array is not a valid FFT");
         if (fft_in->GetScalarType()!=VTK_FLOAT)
            REPORT_ERROR(1,"FFT_phase: Input array is not a valid FFT");

         // Resize Output
         vtkImageData *fft_phase=vtkImageData::New();
         fft_phase->SetScalarType(VTK_FLOAT);
         fft_phase->CopyStructure(fft_in);
         fft_phase->AllocateScalars();

         float *outPtr = (float *) fft_phase->GetScalarPointer();
         SPEED_UP_vtk;
         FOR_ALL_ELEMENTS_IN_VTK(fft_in) {
            *outPtr++=(float) atan2(vtkPtr[1],vtkPtr[0]);
            vtkPtr+=2;
         }

         // Translate to Xmipp
         VTK2xmippArray(fft_phase,phase);
         fft_phase->Delete();
      }

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
