/***************************************************************************
 *
 * Authors: Pedro Antonio de Alarcón (pedro@cnb.uam.es)
 *          Carlos Oscar S. Sorzano
 *          Alberto Pascual Montano
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

#include "../xmippImages.hh"

/* Input (read) from file specifying its dimensions ------------------------ */
template <class T>
void ImageT<T>::read(FileName name, int Ydim, int Xdim, bool reversed,
  Image_Type image_type) _THROW {
  FILE *fh;
  clear(); 

  if ((fh = fopen(name.c_str(), "rb")) == NULL)
    REPORT_ERROR(1501,"Image::read: File " + name + " not found");
  ImageT<T>::read(fh, Ydim, Xdim, reversed, image_type);
  fn_img=name;
  fclose(fh);
}
  
template <class T>
void ImageT<T>::read(FILE * &fh, int Ydim, int Xdim, bool reversed,
  Image_Type image_type) {  
  img.resize(Ydim,Xdim);
  FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(img)
      switch (image_type) {
         case IBYTE:
            unsigned char u;
            FREAD (&u, sizeof(unsigned char), 1, fh, reversed);
            MULTIDIM_ELEM(img,i)=u;
            break;
         case IFLOAT:
            float f;
            FREAD (&f, sizeof(float), 1, fh, reversed);
            MULTIDIM_ELEM(img,i)=f;
            break;
      }
}

// Specific function to read images with complex numbers in them
void ImageT<complex<double> >::read(FILE * &fh, int Ydim, int Xdim, bool reversed,
  Image_Type image_type) 
{
  img.resize(Ydim,Xdim);
  FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(img)
  {
        float a,b;
	    // read real part of a complex number
        FREAD (&a, sizeof(float), 1, fh, reversed);
	    // read imaginary part of a complex number 
	    FREAD (&b, sizeof(float), 1, fh, reversed);
	    // Assign the number
        complex<double> c(a,b);
        MULTIDIM_ELEM(img,i)=c;
	    
  }
}
/* Output (write) ---------------------------------------------------------- */
template <class T>
void ImageT<T>::write(FileName name, bool reversed, Image_Type image_type) _THROW {
  FILE *fp;
  if (name != "") ImageT<T>::rename(name); 

  if ((fp = fopen(fn_img.c_str(), "wb")) == NULL) {
    REPORT_ERROR(1503,"Image::write: File " + fn_img + " cannot be saved");
  };
  ImageT<T>::write(fp, reversed, image_type);
  fclose(fp);  
}

template <class T>
void ImageT<T>::write(FILE * &fh, bool reversed, Image_Type image_type) {
  if (XSIZE(img)==0 || YSIZE(img)==0) return;
  FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(img)
      switch (image_type) {
         case IBYTE:
            unsigned char u;
            u=(unsigned char) MULTIDIM_ELEM(img,i);
            FWRITE (&u, sizeof(unsigned char), 1, fh, reversed);
            break;
         case IFLOAT:
            float f;
            f=(float) MULTIDIM_ELEM(img,i);
            FWRITE (&f, sizeof(float), 1, fh, reversed);
            break;
      }
}

// Specific function to write images with complex numbers in them
void ImageT<complex<double> >::write(FILE * &fh, bool reversed, Image_Type image_type)
{
  FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(img)
  {
            float a,b;
            a=(float) (MULTIDIM_ELEM(img,i)).real();
	        b=(float) (MULTIDIM_ELEM(img,i)).imag();
            FWRITE (&a, sizeof(float), 1, fh, reversed);
	        FWRITE (&b, sizeof(float), 1, fh, reversed);
  }
}

/* Read Xmipp image -------------------------------------------------------- */
template <class T>
void ImageXmippT<T>::read(FileName name, bool skip_type_check, bool reversed) _THROW {
  FILE *fp; 

  ImageXmippT<T>::rename(name); 

  if ((fp = fopen(fn_img.c_str(), "rb")) == NULL)
    REPORT_ERROR(1501,(string)"ImageXmipp::read: File "+fn_img+" not found");
  // Read header
  if (!header.read(fp, skip_type_check, reversed))
     REPORT_ERROR(1502,"ImageXmipp::read: File " + fn_img +
        " is not a valid Xmipp file");   

  // Read whole image and close file
   ImageT<T>::read(fp, header.iYdim(), header.iXdim(), header.reversed(), IFLOAT);
  fclose(fp);

  /* Apply the geometric transformations in the header to the 
     loaded image.
     This is to be compatible with old Xmipp images taht kept 
     geometric transformation information in the header.
     From now on, this information will never be keep in the header
     but in the pixels themself (image).
  */

  // rotate image fAngle1 degrees if necessary
  // img = img.rotate(header.rotAngle());     

  // scale if necessary (check this with Carlos)
  if ((header.Scale() != 0.) && (header.Scale() != 1.)) {
    header.set_dimension(header.Ydim()*header.Scale(), header.Xdim()*header.Scale());
    img.scale_to_size(header.iYdim(), header.iXdim());
  }; 

  header.set_header();  // Set header in a Xmipp consistent state
}

/* Write Xmipp image ------------------------------------------------------- */
template <class T>
void ImageXmippT<T>::write(const FileName &name, bool force_reversed) _THROW {
  FILE *fp;
  if (name != "") ImageXmippT<T>::rename(name); 
  if ((fp = fopen(fn_img.c_str(), "wb")) == NULL)
    REPORT_ERROR(1503,(string)"ImageXmipp::write: File "+fn_img +
       " cannot be written");
  adjust_header();
  bool reversed=(force_reversed) ? !header.reversed():header.reversed();
  header.write(fp, force_reversed);
  ImageT<T>::write(fp, reversed, IFLOAT);
  fclose(fp);  
}

/* Is Xmipp image? --------------------------------------------------------- */
int Is_ImageXmipp(FileName fn, bool skip_type_check,
   bool reversed) _THROW {
   FILE *fp; 
   int result;
   headerXmipp header(headerXmipp::IMG_XMIPP);

   // Open file
   if ((fp = fopen(fn.c_str(), "rb")) == NULL)
     REPORT_ERROR(1501,"Is_ImageXmipp: File " + fn + " not found");

   // Read header
   result = header.read(fp, skip_type_check, reversed);

   fclose(fp);
   
   return result;
}

/* Is Fourier Xmipp image? ------------------------------------------------- */
int Is_FourierImageXmipp(FileName fn, bool skip_type_check,
   bool reversed) _THROW {
   FILE *fp; 
   int result;
   headerXmipp header(headerXmipp::IMG_FOURIER);

   // Open file
   if ((fp = fopen(fn.c_str(), "rb")) == NULL)
     REPORT_ERROR(1501,"Is_FourierImageXmipp: File " + fn + " not found");

   // Read header
   result = header.read(fp, skip_type_check, reversed);

   fclose(fp);
   
   return result;
}

/* Convert a Xmipp Image into  a Fourier Xmipp Image -------------------------------*/
// A simple copy of real numbers is done 
void ImageXmipp_to_FourierImageXmipp(ImageXmipp &I,FourierImageXmipp &F)
{
    // Adjust the size of the Fourier Image
    F().resize(I().RowNo(),I().ColNo());
    // And copy 
 	FOR_ALL_ELEMENTS_IN_MATRIX2D(I())
	{
	    F(i,j)=I(i,j);
	}
}

/* Convert a Fourier Xmipp Image  into a Xmipp Image -------------------------------*/
// A simple copy of real parts is done 
void FourierImageXmipp_to_ImageXmipp(FourierImageXmipp &F,ImageXmipp &I)
{
    // Adjust the size of the Fourier Image
    I().resize(F().RowNo(),F().ColNo());
    // And copy 
 	FOR_ALL_ELEMENTS_IN_MATRIX2D(I())
	{
	    I(i,j)=F(i,j).real();
	}
}

// Initialise an oversampled image (ready for work) ------------------------
void ImageOver::init(int _vmin, int _vmax, int _vistep,
          int _umin, int _umax, int _uistep) {
   overvmin=_vmin;
   overumin=_umin;
   overvmax=_vmax;
   overumax=_umax;
   vistep=_vistep;
   uistep=_uistep;
//   img.init_zeros((_vmax-_vmin+1)*_vistep,(_umax-_umin+1)*_uistep);
   img.init_zeros((_vmax-_vmin)*_vistep+1,(_umax-_umin)*_uistep+1);
   STARTINGY(img)=0;
   STARTINGX(img)=0;
//   STARTINGY(img)=_vmin*_vistep - (_vistep-1)/2;
//   STARTINGX(img)=_umin*_uistep - (_uistep-1)/2;
}

// Window ------------------------------------------------------------------
void ImageOver::window(int _v0, int _u0, int _vF, int _uF) {
   overvmin=_v0;
   overumin=_u0;
   overvmax=_vF;
   overumax=_uF;
   
   int newYdim=(_vF-_v0)*vistep+1;
   int newXdim=(_uF-_u0)*uistep+1;
   img.set_Xmipp_origin();
   img.window(FIRST_XMIPP_INDEX(newYdim),FIRST_XMIPP_INDEX(newXdim),
      LAST_XMIPP_INDEX(newYdim),LAST_XMIPP_INDEX(newXdim));
   STARTINGY(img)=0;
   STARTINGX(img)=0;
}

// Clear -------------------------------------------------------------------
void ImageOver::clear() {
   overvmin=overvmax=0;
   overumin=overumax=0;
   vistep=uistep=0;
   Image::clear();
}   

// Generate the normal image by averaging ----------------------------------
void ImageOver::downsample(Image *I) const {
IMGMATRIX(*I).resize(overvmax-overvmin+1,overumax-overumin+1);
   for (int i=overvmin; i<=overvmax; i++)
      for (int j=overumin; j<=overumax; j++) {
        IMGPIXEL(*I,i,j)=0;
        for (int v=(i-overvmin)*vistep; v<(i+1-overvmin)*vistep; v++)
           for (int u=(j-overumin)*uistep; u<(j+1-overumin)*uistep; u++) {
              IMGPIXEL(*I,i,j) += IMGPIXEL(*this,u,v);
           }
        IMGPIXEL(*I,i,j) /= vistep*uistep;
      }
}

// Generate the oversample image by interpolation --------------------------
void ImageOver::oversample(Image *I) const {}


/* Instantiation ----------------------------------------------------------- */
template <class T>
   void instantiate_images(T t) {
   ImageT<T> img;
   ImageXmippT<T> imgx;     
   FileName name;
   
   Image_Type image_type;
   
   int Ydim,Xdim;
   bool reversed,skip_type_check,force_reversed;
   FILE *fh;
   
   img.read(name,Ydim,Xdim,reversed,image_type);        
   img.read(fh  ,Ydim,Xdim,reversed,image_type);
   img.write(name,reversed,image_type);     
   img.write(fh,  reversed,image_type);
   imgx.read(name,skip_type_check,reversed);
   imgx.write(name,force_reversed);


}

void instantiate_images() {
//   short           s; instantiate_images(s);
//   int             i; instantiate_images(i);
   float           f; instantiate_images(f);
   double          d; instantiate_images(d);
   complex<double> c; instantiate_images(c);
}
