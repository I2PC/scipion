/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Copyright (c) 2000 , CSIC .
 *
 * Permission is granted to copy and distribute this file, for noncommercial
 * use, provided (a) this copyright notice is preserved, (b) no attempt
 * is made to restrict redistribution of this file, and (c) this file is
 * restricted by a compilation copyright.
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 *
 *****************************************************************************/

#ifdef _HAVE_VTK

#include <XmippData/xmippArgs.hh>
#include "../xmippVTK.hh"
#include <typeinfo>
#include <vtkImageFFT.h>
#include <vtkImageRFFT.h>
#include <vtkImageMagnitude.h>
#include <vtkImageFourierCenter.h>
#include <vtkImageExtractComponents.h>

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

template <class T, class VTKT>
   void xmippArray2VTK(const matrix1D<T> &v, VTKT * &retval, int scalarN) {
   if (retval==NULL) retval=VTKT::New();
   else retval->Initialize();
   retval->SetDimensions(XSIZE(v),1,1);
   retval->SetNumberOfScalarComponents(scalarN);
   SET_VTK_TYPE(T,retval);
   COPY_VALUES_TO_VTK(T,v,retval,scalarN);
   retval->Update();
}

template <class T, class VTKT>
   void xmippArray2VTK(const matrix2D<T> &v, VTKT * &retval, int scalarN) {
   if (retval==NULL) retval=VTKT::New();
   else retval->Initialize();
   retval->SetDimensions(XSIZE(v),YSIZE(v),1);
   retval->SetNumberOfScalarComponents(scalarN);
   SET_VTK_TYPE(T,retval);
   COPY_VALUES_TO_VTK(T,v,retval,scalarN);
   retval->Update();
}

template <class T, class VTKT>
   void xmippArray2VTK(const matrix3D<T> &v, VTKT * &retval, int scalarN) {
   if (retval==NULL) retval=VTKT::New();
   else retval->Initialize();
   retval->SetDimensions(XSIZE(v),YSIZE(v),ZSIZE(v));
   retval->SetNumberOfScalarComponents(scalarN);
   SET_VTK_TYPE(T,retval);
   COPY_VALUES_TO_VTK(T,v,retval,scalarN);
   retval->Update();
}

void xmippFFT2VTK(FourierImageXmipp &v, vtkImageData * &retval) _THROW {

   if (retval==NULL) retval=vtkImageData::New();
   else retval->Initialize();
   retval->SetDimensions(XSIZE(v()),YSIZE(v()),1);
   retval->SetNumberOfScalarComponents(2);
   retval->SetScalarType(VTK_FLOAT);
   retval->AllocateScalars();
   retval->Update();
   
   float *ptr= (float *) retval->GetScalarPointer();    
   FOR_ALL_ELEMENTS_IN_MATRIX2D(v())
   {
      *ptr++=(v(i,j)).real();
	  *ptr++=(v(i,j)).imag(); 
   }   
}

void xmippFFT2VTK(matrix1D <complex <double > > &v, vtkImageData * &retval) _THROW {

   if (retval==NULL) retval=vtkImageData::New();
   else retval->Initialize();
   retval->SetDimensions(XSIZE(v),1,1);
   retval->SetNumberOfScalarComponents(2);
   retval->SetScalarType(VTK_FLOAT);
   retval->AllocateScalars(); 
   retval->Update();
   
   float *ptr= (float *) retval->GetScalarPointer();    
   FOR_ALL_ELEMENTS_IN_MATRIX1D(v)
   {
      *ptr++=(v(i)).real();
	  *ptr++=(v(i)).imag(); 
   }   
}

void xmippFFT2VTK(FourierVolumeXmipp &v, vtkImageData * &retval) _THROW
{

   if (retval==NULL) retval=vtkImageData::New();
   else retval->Initialize();
   retval->SetDimensions(XSIZE(v()),YSIZE(v()),ZSIZE(v()));
   retval->SetNumberOfScalarComponents(2);
   retval->SetScalarType(VTK_FLOAT);
   retval->AllocateScalars(); 
   retval->Update();
   
   float *ptr= (float *) retval->GetScalarPointer();    
   FOR_ALL_ELEMENTS_IN_MATRIX3D(v())
   {
      *ptr++=(v(k,i,j)).real();
	  *ptr++=(v(k,i,j)).imag(); 
   }   
}

/* Xmipp resize after VTK -------------------------------------------------- */
template <class T, class VTKT>
   void xmippArray_resize_VTK(matrix1D<T> &retval, VTKT *v) _THROW {
   int dim[3]; v->GetDimensions(dim);
   if (dim[1]!=1 && dim[2]!=1)
      REPORT_ERROR(1,"VTK2xmippVector: VTK image is not a vector");
   /* CO:
   int scalar_components=v->GetNumberOfScalarComponents();
   if (scalar_components!=1)
      REPORT_ERROR(1,"VTK2xmippVector: VTK image does not contain scalars");
   */

   retval.resize(dim[0]);
}

template <class T, class VTKT>
   void xmippArray_resize_VTK(matrix2D<T> &retval, VTKT *v) _THROW {
   int dim[3]; v->GetDimensions(dim);
   if (dim[2]!=1)
      REPORT_ERROR(1,"VTK2xmippMatrix: VTK image is not a matrix");
   /* CO:
   int scalar_components=v->GetNumberOfScalarComponents();
   if (scalar_components!=1)
      REPORT_ERROR(1,"VTK2xmippMatrix: VTK image does not contain scalars");
   */

   retval.resize(dim[1],dim[0]);
}

template <class T, class VTKT>
   void xmippArray_resize_VTK(matrix3D<T> &retval, VTKT *v) _THROW {
   int dim[3]; v->GetDimensions(dim);
   /*
   int scalar_components=v->GetNumberOfScalarComponents();
   if (scalar_components!=1)
      REPORT_ERROR(1,"VTK2xmippVolume: VTK image does not contain scalars");
   */

   retval.resize(dim[2],dim[1],dim[0]);
}

/* VTK --> Xmipp ----------------------------------------------------------- */
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
   

template <class T, class VTKT>
   void VTK2xmippArray(VTKT *v, matrix1D<T> &retval) {
   retval.clear();
   if (v==NULL) return;
   COPY_VALUES_FROM_VTK(T,v,retval);
}

template <class T, class VTKT>
   void VTK2xmippArray(VTKT *v, matrix2D<T> &retval) {
   retval.clear();
   if (v==NULL) return;
   COPY_VALUES_FROM_VTK(T,v,retval);
}

template <class T, class VTKT>
   void VTK2xmippArray(VTKT *v, matrix3D<T> &retval) {
   retval.clear();
   if (v==NULL) return;
   COPY_VALUES_FROM_VTK(T,v,retval);
}

/* VTK -> VTK -------------------------------------------------------------- */
template <class VTKT>
   void VTK2VTK(VTKT *v_in, VTKT *&v_out, bool change_type, int new_type) {
   if (v_out==NULL)  v_out=VTKT::New();
   else              v_out->PrepareForNewData();
   if (v_in!=NULL) {
      if (!change_type) v_out->SetScalarType(v_in->GetScalarType());
      else              v_out->SetScalarType(new_type);
      v_out->SetNumberOfScalarComponents(v_in->GetNumberOfScalarComponents());
      v_out->CopyStructure(v_in);
      v_out->CopyAndCastFrom(v_in,v_in->GetExtent());
      v_out->Update();
   }
}

/* Same shape -------------------------------------------------------------- */
//#define DEBUG
bool same_shape(vtkImageData *v1, vtkImageData *v2) {
   int dim1[3]; v1->GetDimensions(dim1);
   int dim2[3]; v2->GetDimensions(dim2);
   #ifdef DEBUG
      cout << "v1 (ZxYxX):" << dim1[2] << " " << dim1[1] << " " << dim1[0] << endl
    	   << "v2 (ZxYxX):" << dim2[2] << " " << dim2[1] << " " << dim2[0] << endl
    	   << "Scalartypes: " << v1->GetScalarType() << " "
    	   << v2->GetScalarType() << endl
    	   << "ScalarN: " << v1->GetNumberOfScalarComponents() << " "
           << v2->GetNumberOfScalarComponents() << endl;
   #endif
   if (dim1[0]!=dim2[0] || dim1[1]!=dim2[1] || dim1[2]!=dim2[2])
      return FALSE;
   if (v1->GetScalarType()!=v2->GetScalarType()) return FALSE;
   if (v1->GetNumberOfScalarComponents()!=v2->GetNumberOfScalarComponents())
      return FALSE;
   return TRUE;
}

/* Center FFT -------------------------------------------------------------- */
void CenterFFT(vtkImageData *&v) {
   if (v==NULL) return;
   vtkImageFourierCenter *fftcenter=vtkImageFourierCenter::New();
   fftcenter->SetInput(v); fftcenter->Update();
   VTK2VTK(fftcenter->GetOutput(),v);
   fftcenter->Delete();
}



/* VTK -> Xmipp FFT -------------------------------------------------------- */
void VTK2xmippFFT(vtkImageData *v, matrix1D< complex <double> > &retval) {
   
   retval.clear();
   if (v==NULL) return;
   int dim[3]; v->GetDimensions(dim);
   int scalar_components=v->GetNumberOfScalarComponents();
   if (dim[2]!=1 || dim[1]!=1)
      REPORT_ERROR(1,"VTK2xmippFFT: VTK image is not a vector");
   if (scalar_components!=2)
      REPORT_ERROR(1,"VTK2xmippFFT: VTK image does not contain complexes");
   
   retval.resize(dim[0]);
   
   v->SetScalarType(VTK_FLOAT);
   v->AllocateScalars(); 
   float *ptr= (float *) v->GetScalarPointer();  
   FOR_ALL_ELEMENTS_IN_MATRIX1D(retval)
   {
      complex<double> c(*ptr,*(ptr+1));
      retval(i)=c;
	  ptr+=2;
   }                
}


/* VTK -> Xmipp FFT -------------------------------------------------------- */
void VTK2xmippFFT(vtkImageData *v, FourierImageXmipp &retval) {
   
   retval.clear();
   if (v==NULL) return;
   int dim[3]; v->GetDimensions(dim);
   int scalar_components=v->GetNumberOfScalarComponents();
   if (dim[2]!=1)
      REPORT_ERROR(1,"VTK2xmippFFT: VTK image is not a matrix");
   if (scalar_components!=2)
      REPORT_ERROR(1,"VTK2xmippFFT: VTK image does not contain complexes");
   
   retval().resize(dim[1],dim[0]);
   
   v->SetScalarType(VTK_FLOAT);
   v->AllocateScalars(); 
   float *ptr= (float *) v->GetScalarPointer();  
   FOR_ALL_ELEMENTS_IN_MATRIX2D(retval())
   {
      complex<double> c(*ptr,*(ptr+1));
      retval(i,j)=c;
	  ptr+=2;
   }                
}

void VTK2xmippFFT(vtkImageData *v, FourierVolumeXmipp &retval)
{
   retval.clear();
   if (v==NULL) return;
   int dim[3]; v->GetDimensions(dim);
   int scalar_components=v->GetNumberOfScalarComponents();
   if (scalar_components!=2)
      REPORT_ERROR(1,"VTK2xmippFFT: VTK volume does not contain complexes");

   retval().resize(dim[2],dim[1],dim[0]);
   
   v->SetScalarType(VTK_FLOAT);
   v->AllocateScalars(); 
   float *ptr= (float *) v->GetScalarPointer();  
   FOR_ALL_ELEMENTS_IN_MATRIX3D(retval())
   {
      complex<double> c(*ptr,*(ptr+1));
      retval(k,i,j)=c;
	  ptr+=2;
   }             
}

/* Show VTK ---------------------------------------------------------------- */
template <class VTKT>
   void VTK_print(ostream &out, VTKT *v) {
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

/* FFT --------------------------------------------------------------------- */
void FFT_VTK(vtkImageData *v_in, vtkImageData * &fft_out,
   bool do_not_center) {
   vtkImageFFT *fft=vtkImageFFT::New();
   fft->SetInput((vtkStructuredPoints *) v_in); fft->Update();
   VTK2VTK(fft->GetOutput(),fft_out);
   fft->Delete();
   if (!do_not_center) CenterFFT(fft_out);
}

/* Magnitude --------------------------------------------------------------- */
template <class maT>
   void FFT_magnitude(vtkImageData *fft_in, maT &mag) {
      vtkImageMagnitude *mag_filter=vtkImageMagnitude::New();
      mag_filter->SetInput(fft_in); mag_filter->Update();
      VTK2xmippArray(mag_filter->GetOutput(),mag);
      mag_filter->Delete();
}

/* Phase ------------------------------------------------------------------- */
/* This is a very short exercise of implementing a VTK filter. It is
   inspired in vtkImageMagnitude */
template <class maT>
   void FFT_phase(vtkImageData *fft_in, maT &phase) _THROW {
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

/* IFFT -------------------------------------------------------------------- */
void IFFT_VTK(vtkImageData *fft_in, vtkImageData *&v_out,
   bool is_not_centered) {
   vtkImageRFFT *rfft=vtkImageRFFT::New();
   if (is_not_centered) {
      rfft->SetInput(fft_in); rfft->Update();
   } else {
      vtkImageFourierCenter *fftcenter=vtkImageFourierCenter::New();
      fftcenter->SetInput(fft_in); fftcenter->Update();
      rfft->SetInput(fftcenter->GetOutput()); rfft->Update();
      fftcenter->Delete();
   }

   vtkImageExtractComponents *real_part=vtkImageExtractComponents::New();
   real_part->SetInput(rfft->GetOutput());
   real_part->SetComponents(0);
   real_part->Update();
   
   VTK2VTK(real_part->GetOutput(),v_out);
   rfft->Delete();
   real_part->Delete();
}


/* Correlation and Autocorrelation ----------------------------------------- */
template <class T>
void auto_correlation_matrix(const matrix2D<T> &Img,matrix2D<double> &R)
{
   R.resize(Img);
   // Compute the Fourier Transform
   vtkImageData *VTK_FourierTransform=NULL;
   FourierImageXmipp FFT1,FFT2; 
   FFT_VTK(Img, VTK_FourierTransform,TRUE);
   VTK2xmippFFT (VTK_FourierTransform, FFT1);   

   FFT2().resize(FFT1());
   // Compute the complex conjugate
   FOR_ALL_ELEMENTS_IN_MATRIX2D(FFT2())
   {
      complex<double> A(FFT1(i,j).real(),-FFT1(i,j).imag());
	  FFT2(i,j)=A;
   }
   // Multiply the matrices element by element
   mul_elements(FFT1(),FFT2(),FFT1());
   // Invert the product, in order to obtain the correlation image
   xmippFFT2VTK (FFT1, VTK_FourierTransform);
   IFFT_VTK(VTK_FourierTransform, R,TRUE);

   double N=R.RowNo()*R.ColNo();
   FOR_ALL_ELEMENTS_IN_MATRIX2D(R)
      R(i,j)=R(i,j)/N;
	  
   // Center the resulting image to obtain a centered autocorrelation
   int ym=(int)R.RowNo()/2;
   int xm=(int)R.ColNo()/2;
   int xp,yp;
   for(int i=0;i<R.RowNo();i++)   
      for(int j=0;j<R.ColNo()/2;j++)
	  {	  
    	  // Determine the destination row
		  if(i>=ym)
			 yp=i-ym;
		  else 
			 yp=i+ym;
		  // Determine the destination column
    	  xp=j+xm;
		  // Swap origin and destination
		  double tmp;
		  SWAP(R(i,j),R(yp,xp),tmp);	  
      }
}

template <class T>
void correlation_matrix(const matrix2D<T> &m1, const matrix2D<T> &m2,
                              matrix2D<double> &R)
{
   // Compute the Fourier Transform of the images
    vtkImageData *VTK_FourierTransform=NULL;
    FourierImageXmipp FFT1,FFT2; 
    FFT_VTK(m1, VTK_FourierTransform,TRUE);
    VTK2xmippFFT (VTK_FourierTransform, FFT1);
    FFT_VTK(m2, VTK_FourierTransform,TRUE);
    VTK2xmippFFT (VTK_FourierTransform, FFT2);

    // Compute the complex conjugate
    
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(FFT2())
    {
       complex<double> A(DIRECT_MAT_ELEM(FFT2(),i,j).real(),-
                         DIRECT_MAT_ELEM(FFT2(),i,j).imag());
	   DIRECT_MAT_ELEM(FFT2(),i,j)=A;
    }
    // Multiply the matrices element by element
    mul_elements(FFT1(),FFT2(),FFT1());
    // Invert the product, in order to obtain the correlation image
    matrix2D<double> R_aux(R);    

    xmippFFT2VTK (FFT1, VTK_FourierTransform);
    IFFT_VTK(VTK_FourierTransform, R_aux,TRUE);

    int xinit = 0;//R.startingX(); Let's use direct coordinate to speed this up
    int yinit = 0;//R.startingY();
    int xend  = R.ColNo()-1;//last index in direct_mat
    int yend  = R.RowNo()-1;
    
    // Divide by the number of pixels to obtain the correlation value used 
    // in Xmipp
    double N=R.RowNo()*R.ColNo();
    int ii,jj;

    int ym=(int)R.RowNo()/2;
    int xm=(int)R.ColNo()/2;

    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(R)
       {
	ii=intWRAP(i+ym,yinit,yend);
	jj=intWRAP(j+xm,xinit,xend);
	DIRECT_MAT_ELEM (R,ii,jj)= DIRECT_MAT_ELEM (R_aux,i,j)/N;
        }
}

/* Convolution of series --------------------------------------------------- */
template <class T>
void series_convolution(matrix1D<T> &series1,matrix1D<T> &series2,
                        matrix1D<T> &result,bool FullConvolution)
{
	// Store dimension of series
	int dim1=series1.get_dim();
	int dim2=series2.get_dim();
	// Resize series to the size of the resulting series
	// (Zeros are stored in the expanded values)
	series1.resize(dim1+dim2-1);
	series2.resize(dim1+dim2-1);
	result.resize(dim1+dim2-1);
	// Fourier Transform the two series
    vtkImageData *VTK_FourierTransform=NULL;
    matrix1D<complex <double> > FFT1,FFT2;
    FFT_VTK(series1,VTK_FourierTransform,TRUE);
    VTK2xmippFFT (VTK_FourierTransform, FFT1);
    FFT_VTK(series2, VTK_FourierTransform,TRUE);
    VTK2xmippFFT (VTK_FourierTransform, FFT2);
    // Multiply the vectors element by element to do the 
	// convolution in the Fourier space.
    FOR_ALL_ELEMENTS_IN_MATRIX1D(FFT1)
		FFT1(i)=FFT1(i)*FFT2(i);
	// Recover the convolution result by inverse FFT
	result.resize(series1.get_dim()*2-1);	// Resize to store the result
    xmippFFT2VTK (FFT1, VTK_FourierTransform);
    IFFT_VTK(VTK_FourierTransform, result,TRUE);
	// Restore the dimensions
	series1.resize(dim1);
	series2.resize(dim2);
	
	/* If the full convolution is required, nothing more remains to be done.
	Otherwise, if the valid values are required, return the central ones. */
	if(FullConvolution==FALSE)
	{
		// First get the maximum dimension of the series, which is the dimension
		// of the result
		int dim=MAX(dim1,dim2);
		// Determine the number of values to discard
		int discard=result.get_dim()-dim;
		// Divide it by two as we have to discard them in both sides of the vector
		discard=discard/2;  // Integer division is intended
		// copy required values (simple displacement of values)
		for(int i=STARTINGX(result);i<STARTINGX(result)+dim;i++)
			result(i)=result(i+discard);
		// and finally resize to discard not copied values
		result.resize(dim);
	}
}

/* Convolution of a series with a given filter --------------------------- */
template <class T>
void convolve(matrix1D<T> &series1,matrix1D<T> &filter, matrix1D<T> &result)
{
	// Store dimension of series
	int dim1=series1.get_dim();
	int dim2=filter.get_dim();
	// Resize series to the size of the resulting series
	// (Zeros are stored in the expanded values)
	series1.resize(dim1+dim2-1);
	filter.resize(dim1+dim2-1);
	result.resize(dim1+dim2-1);
	// Fourier Transform the two series
    vtkImageData *VTK_FourierTransform=NULL;
    matrix1D<complex <double> > FFT1,FFT2;
    FFT_VTK(series1,VTK_FourierTransform,TRUE);
    VTK2xmippFFT (VTK_FourierTransform, FFT1);
    FFT_VTK(filter, VTK_FourierTransform,TRUE);
    VTK2xmippFFT (VTK_FourierTransform, FFT2);
    // Multiply the vectors element by element to do the 
	// convolution in the Fourier space.
    FOR_ALL_ELEMENTS_IN_MATRIX1D(FFT1)
	{
		FFT1(i)=FFT1(i)*FFT2(i);
	}	
	// Recover the convolution result by inverse FFT
	result.resize(series1.get_dim()*2-1);	// Resize to store the result
    xmippFFT2VTK (FFT1, VTK_FourierTransform);
    IFFT_VTK(VTK_FourierTransform, result,TRUE);
	
	
	// Restore the dimensions
	series1.resize(dim1);
	filter.resize(dim2);
	// Determine the number of values to discard
	int discard=result.get_dim()-dim1;
	// Divide it by two as we have to discard them in both sides of the vector
	discard=discard/2;  // Integer division is intended
	// copy required values (simple displacement of values)
	for(int i=STARTINGX(result);i<STARTINGX(result)+dim1;i++)
		result(i)=result(i+discard);
	// and finally resize to discard not copied values
	result.resize(dim1);
}


/* Numerical derivative of a matrix ----------------------------- */
void numerical_derivative(matrix2D<double> &M,matrix2D<double> &D,
							char direction,int order,
							int window_size,int polynomial_order)
{
    // Set D to be a copy in shape of M
	D.copy_shape(M);
	matrix1D<double> v,rotated; 
	matrix1D<double> ans; // To obtain results 
	// Wrap around version of the Savitzky-Golay coefficients
    rotated.resize(2*window_size+1);

	double *pans=ans.adapt_for_numerical_recipes();
	double *pv=v.adapt_for_numerical_recipes();
	double *protated=rotated.adapt_for_numerical_recipes();
	// Calculate the Savitzky-Golay filter coeficients
	 savgol(protated, 2*window_size+1, window_size,
		     window_size, order, polynomial_order);
     // Savitzky-Golay filter is returned in wrap-around style, so
	 // correct it to use with the convolution routine
	 int dim=rotated.get_dim();
	 matrix1D<double> coeficients(dim);
	 for(int i=0;i<dim;i++)
	 {
	 	int j=i+window_size;
		if(j<dim)
			coeficients(j)=rotated(i);
		else
			coeficients(j-dim)=rotated(i);
	 }
     // Apply the Savitzky-Golay filter to every row or column
	if(direction=='x')
	{
		 // For every row (values in a row are values of the X direction)
		 for(int i=STARTINGY(M);i<=FINISHINGY(M);i++)
		 {
			M.getRow(i,v); 
			convolve(v,coeficients,ans);
			ans.setRow();
			D.setRow(i,ans);			
		 }
	}
	else if(direction=='y')
	{
		 // For every column (values in a column are values of the Y direction)
		 for(int i=STARTINGX(M);i<=FINISHINGX(M);i++)
		 {
			M.getCol(i,v);
			convolve(v,coeficients,ans);
			ans.setCol();
			D.setCol(i,ans);						
		 }
	}	

	ans.kill_adaptation_for_numerical_recipes(pans);
	v.kill_adaptation_for_numerical_recipes(pv);
	coeficients.kill_adaptation_for_numerical_recipes(protated);
}


/* Instantiate ------------------------------------------------------------- */
template <class T, class VTKT>
void instantiate_VTK(T a, VTKT *vtkI) {
   matrix1D<T> v;
   matrix2D<T> m;
   matrix3D<T> V;
   xmippArray2VTK(v,vtkI);
   xmippArray2VTK(m,vtkI);
   xmippArray2VTK(V,vtkI);
   VTK2xmippArray(vtkI,v);
   VTK2xmippArray(vtkI,m);
   VTK2xmippArray(vtkI,V);
   VTK2VTK(vtkI,vtkI);
   VTK_print(cout,vtkI);
   
   matrix2D<double> R;
   int ca;
   auto_correlation_matrix(m,R);
   correlation_matrix(m,m,R);
   series_convolution(v,v,v);
   convolve(v,v,v);
}


void instantiate_VTK() {
   vtkImageData        *vtkID;
   vtkStructuredPoints *vtkSP;
   unsigned char u; instantiate_VTK(u,vtkID); instantiate_VTK(u,vtkSP);
   short int     s; instantiate_VTK(s,vtkID); instantiate_VTK(s,vtkSP);
   int           i; instantiate_VTK(i,vtkID); instantiate_VTK(i,vtkSP);
   float         f; instantiate_VTK(f,vtkID); instantiate_VTK(f,vtkSP);
   double        d; instantiate_VTK(d,vtkID); instantiate_VTK(d,vtkSP);

   matrix2D<double> m; FFT_magnitude(vtkID,m); FFT_phase(vtkID,m);
   matrix3D<double> V; FFT_magnitude(vtkID,V); FFT_phase(vtkID,V);
}
#endif
