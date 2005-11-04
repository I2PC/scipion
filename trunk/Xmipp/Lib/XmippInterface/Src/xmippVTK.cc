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

#ifdef _HAVE_VTK

#include <XmippData/xmippArgs.hh>
#include "../xmippVTK.hh"
#include <typeinfo>
#include <vtkImageFFT.h>
#include <vtkImageRFFT.h>
#include <vtkImageMagnitude.h>
#include <vtkImageExtractComponents.h>

void xmippFFT2VTK(matrix1D <complex <double > > &v, vtkImageData * &retval) {
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

void xmippFFT2VTK(FourierImageXmipp &v, vtkImageData * &retval) {
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
      *ptr++=(float)(v(i,j).real());
         *ptr++=(float)(v(i,j).imag());
   }   
}

void xmippFFT2VTK(FourierVolumeXmipp &v, vtkImageData * &retval)
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
void VTK_CenterFFT(vtkImageData *&v, int zoff, int yoff, int xoff) {
   if (v==NULL) return;
   vtkImageData *v_centered=NULL;
   VTK_CenterFFT(v,v_centered, zoff, yoff, xoff);
   VTK2VTK(v_centered,v);
   v_centered->Delete();
}

void VTK_CenterFFT(vtkImageData *&v_in, vtkImageData *&v_out,
   int zoff, int yoff, int xoff) {
   if (v_in==NULL) return;

   SPEED_UP_vtk;
   VTK_resize_VTK(v_in,v_out); // Copy the structure
   float *v_inPtr =(float *)v_in->GetScalarPointer();
   float *v_outPtr=(float *)v_out->GetScalarPointer();
   
   // Compute the shift
   inExt=v_in->GetExtent();
   maxX = inExt[1] - inExt[0];
   maxY = inExt[3] - inExt[2];
   maxZ = inExt[5] - inExt[4];
   int maxXY=(maxX+1)*(maxY+1);
   int shift_X=(int)((double)(maxX+1)/2.0);
   int shift_Y=(int)((double)(maxY+1)/2.0);
   int shift_Z=(int)((double)(maxZ+1)/2.0);
   // Copy the elements applying the appropiate shift
   int dim=v_in->GetNumberOfScalarComponents();
   for (int k = 0; k <= maxZ; k++) {
       int k0=k*maxXY*dim; // pointer to the slice k
       for (int i = 0; i <= maxY; i++) {
           int i0=k0+i*(maxX+1)*dim; // pointer to the column i
           for (int j = 0; j <= maxX; j++) {
              float *ptr_j=v_inPtr+i0+j*dim; // pointer to the element
              
              // New coordinate
              int kp=intWRAP(k+shift_Z+zoff,0,maxZ);
              int ip=intWRAP(i+shift_Y+yoff,0,maxY);
              int jp=intWRAP(j+shift_X+xoff,0,maxX);
              
              // Compute the shift of this element
              float *ptr_jp=v_outPtr+dim*(kp*maxXY+ip*(maxX+1)+jp);
              
              // Now copy all elements from 1 position to the other
              for (int d=0; d<dim; d++) *ptr_jp++=*ptr_j++;
           }
       }
   }
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

/* FFT --------------------------------------------------------------------- */
void FFT_VTK(vtkImageData *v_in, vtkImageData * &fft_out,
   bool do_not_center) {
   vtkImageFFT *fft=vtkImageFFT::New();
   fft->SetInput((vtkStructuredPoints *) v_in); fft->Update();
   VTK2VTK(fft->GetOutput(),fft_out);
   fft->Delete();
   if (!do_not_center) VTK_CenterFFT(fft_out);
}

/* Magnitude --------------------------------------------------------------- */
void VTK_FFT_magnitude(FourierImageXmipp &fft_in,
   matrix2D<double> &mag, bool do_not_center) {
   vtkImageData *fftI=NULL; xmippFFT2VTK(fft_in,fftI);
   if (!do_not_center) VTK_CenterFFT(fftI);
   VTK_FFT_magnitude(fftI,mag);
   fftI->Delete();
}

/* IFFT -------------------------------------------------------------------- */
void IFFT_VTK(vtkImageData *fft_in, vtkImageData *&v_out,
   bool is_not_centered) {
   vtkImageRFFT *rfft=vtkImageRFFT::New();
   if (is_not_centered) {
      rfft->SetInput(fft_in); rfft->Update();
   } else {
      int dim[3]; fft_in->GetDimensions(dim);
      int zoff=(dim[2]%2==1)?1:0;
      int yoff=(dim[1]%2==1)?1:0;
      int xoff=(dim[0]%2==1)?1:0;
      cout << "Offset: " << zoff << "," << yoff << "," << xoff << endl;
      vtkImageData *fft_not_centered=NULL;
      VTK_CenterFFT(fft_in,fft_not_centered,zoff,yoff,xoff);
      rfft->SetInput(fft_not_centered); rfft->Update();
      fft_not_centered->Delete();
   }

   vtkImageExtractComponents *real_part=vtkImageExtractComponents::New();
   real_part->SetInput(rfft->GetOutput());
   real_part->SetComponents(0);
   real_part->Update();
   
   VTK2VTK(real_part->GetOutput(),v_out);
   rfft->Delete();
   real_part->Delete();
}
#endif
