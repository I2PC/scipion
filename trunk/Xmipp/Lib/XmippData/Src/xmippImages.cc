/***************************************************************************
 *
 * Authors: Pedro Antonio de Alarcón (pedro@cnb.uam.es)
 *          Carlos Oscar S. Sorzano
 *          Alberto Pascual Montano
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Part of this module has been developed by Lorenzo Zampighi and Nelson Tang
 * Dept. Physiology of the David Geffen School of Medicine
 * Univ. of California, Los Angeles.
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
#include "../xmippImagic.hh"
#include "../xmippGeometry.hh"
 
/* Input (read) from file specifying its dimensions ------------------------ */
template <class T>
bool ImageT<T>::read(const FileName &name, float fIform, int Ydim,
  int Xdim, bool reversed, Image_Type image_type) _THROW {
  FILE *fh;
  clear(); 

  if ((fh = fopen(name.c_str(), "rb")) == NULL)
    REPORT_ERROR(1501,"Image::read: File " + name + " not found");
  bool ret;
  if ((ret = ImageT<T>::read(fh, fIform, Ydim, Xdim, reversed, image_type)))
     fn_img=name;
  fclose(fh);
  return (ret);
}
  
template <class T>
bool ImageT<T>::read(FILE * &fh, float fIform, int Ydim, int Xdim,
   bool reversed, Image_Type image_type) {  
  img.resize(Ydim,Xdim);

  FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(img)
      switch (image_type) {
         case IBYTE:
            unsigned char u;
            FREAD (&u, sizeof(unsigned char), 1, fh, reversed);
            MULTIDIM_ELEM(img,i)=u;
            break;
         case I16:
            unsigned short us;
            FREAD (&us, sizeof(unsigned short), 1, fh, reversed);
            MULTIDIM_ELEM(img,i)=us;
            break;
         case IFLOAT:
            float f;
            FREAD (&f, sizeof(float), 1, fh, reversed);
            MULTIDIM_ELEM(img,i)=f;
            break;
      }
   return (TRUE);
}

// Specific function to read images with complex numbers in them
bool ImageT<complex<double> >::read(FILE * &fh, float fIform,
   int Ydim, int Xdim, bool reversed, Image_Type image_type)
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
  return (true);
}

// LoadImage loads an image from file, returning a pointer to the generic
// Image interface.
template <class T>
ImageT<T> *ImageT<T>::LoadImage (FileName name, bool apply_geo)
{
  ImageT<T> *ret = NULL;
  if (name.find (IMAGIC_TAG) != string::npos)
  {
    ImageImagicT<T> *i = new ImageImagicT<T>();
    if (i->read (name))
      ret = i;
    else
      delete i;
  }
  else // For now, assume anything else is ImageXmipp type.
  {
    ImageXmippT<T> *i = new ImageXmippT<T>();
    if (i->read (name,FALSE,FALSE,apply_geo,FALSE))
      ret = i;
    else
      delete i;
  }
  return (ret);
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
  double a,b;
  if (image_type!=IFLOAT) {
     double min_val, max_val;
     (*this)().compute_double_minmax(min_val,max_val);
     if (image_type==IBYTE) a=255;
     else                   a=65535;
     a/=(max_val-min_val);
     b=min_val;
  }
  
  FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(img)
      switch (image_type) {
         case IBYTE:
            unsigned char u;
            u=(unsigned char) ROUND(a*(MULTIDIM_ELEM(img,i)-b));
            FWRITE (&u, sizeof(unsigned char), 1, fh, reversed);
            break;
	 case I16:
	    unsigned short us;	    
            us=(unsigned short) ROUND(a*(MULTIDIM_ELEM(img,i)-b));
            FWRITE (&us, sizeof(unsigned short), 1, fh, reversed);
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

/* Check OldXmipp header location for image offsets ----------------------- */
template <class T>
void ImageXmippT<T>::check_oldxmipp_header() {
  if (header.fXoff()==0. && header.fYoff()==0. && header.fZoff()==0.) {
    // Might be oldXmipp image header
    float ang=header.old_rot();   
    matrix2D<double> mat=header.fGeo_matrix();
    if (ABS(mat(0,0)-COSD(ang))<XMIPP_EQUAL_ACCURACY &&
	ABS(mat(1,1)-COSD(ang))<XMIPP_EQUAL_ACCURACY &&
	ABS(mat(0,1)-SIND(ang))<XMIPP_EQUAL_ACCURACY &&
	ABS(mat(1,0)+SIND(ang))<XMIPP_EQUAL_ACCURACY &&
        mat(2,0)!=0. && mat(2,1)!=0. ) {
      // This indeed seems to be an OldXmipp style header with non-zero offsets
      cerr << "WARNING%% Copying shifts from old to new header location: "<< (*this).name() <<endl;
      header.fXoff()=-(float)mat(2,0);
      header.fYoff()=-(float)mat(2,1);
      if (XSIZE(ImageT<T>::img)%2==0) header.fXoff()+=0.5;
      if (YSIZE(ImageT<T>::img)%2==0) header.fYoff()+=0.5;
    } 
  } 
}

/* Get geometric transformation matrix from 2D-image header ---------------- */
template <class T>
matrix2D<double> ImageXmippT<T>::get_transformation_matrix(bool only_apply_shifts) {
  matrix2D<double>    A(3,3);
  double psi=realWRAP(header.Psi(),-180,180);
  double theta=realWRAP(header.Theta(),-180,180);

  A.init_identity();

  /* This is to be compatible with old Xmipp images, that store image
     translations in another position of the header */
  ImageXmippT<T>::check_oldxmipp_header();

  if (only_apply_shifts) {
    Euler_angles2matrix(0.,0.,0.,A);
    A(0,2)=-header.fXoff();   
    A(1,2)=-header.fYoff();
  } else {

    if (theta==0.) {
      // For untilted images: apply Euler matrix
      Euler_angles2matrix (header.Phi(),0.,header.Psi(),A);
    } else { 
      // For tilted images: only apply Psi 
      // Take another_set into account
      if (theta<0.) {
	theta=-theta;
	psi=realWRAP(psi-180.,-180,180);
      }
      Euler_angles2matrix(0.,0.,psi,A);
    }
    A(0,2)=-header.fXoff();   
    A(1,2)=-header.fYoff();
  }

  // Also for only_apply_shifts: mirror if necessary!
  if (header.Flip()==1) {
    A(0,0)=-A(0,0); 
    A(0,1)=-A(0,1); 
  }

  return A;
}

/* Read Xmipp image -------------------------------------------------------- */
template <class T>
bool ImageXmippT<T>::read(const FileName &name, bool skip_type_check,
  bool reversed, bool apply_geo, bool only_apply_shifts) _THROW {
  FILE *fp;
  bool ret;

  ImageXmippT<T>::rename(name); 
  if ((fp = fopen(ImageT<T>::fn_img.c_str(), "rb")) == NULL)
    REPORT_ERROR(1501,(string)"ImageXmipp::read: File "+ImageT<T>::fn_img+" not found");
  // Read header
  if (!header.read(fp, skip_type_check, reversed))
     REPORT_ERROR(1502,"ImageXmipp::read: File " + ImageT<T>::fn_img +
        " is not a valid Xmipp file");   

  // Read whole image and close file
  if ((ret = ImageT<T>::read(fp, header.fIform(), header.iYdim(), header.iXdim(), header.reversed(), IFLOAT)))
  {

    if (apply_geo || only_apply_shifts) {
      // Apply the geometric transformations in the header to the loaded image.
      // Transform image without wrapping, set new values to first element in the matrix
      T  outside=DIRECT_MAT_ELEM(ImageT<T>::img,0,0);
      ImageT<T>::img.self_apply_geom(ImageXmippT<T>::get_transformation_matrix(only_apply_shifts),IS_INV,DONT_WRAP,outside);
    }

    // scale if necessary (check this with Carlos)
    if ((header.Scale() != 0.) && (header.Scale() != 1.)) {
      header.set_dimension(header.Ydim()*header.Scale(), header.Xdim()*header.Scale());
      ImageT<T>::img.scale_to_size(header.iYdim(), header.iXdim());
    }; 

    header.set_header();  // Set header in a Xmipp consistent state

  }
  fclose(fp);
  return (ret);
}

/* Write Xmipp image ------------------------------------------------------- */
template <class T>
void ImageXmippT<T>::write(const FileName &name, bool force_reversed) _THROW {
  FILE *fp;
  if (name != "") ImageXmippT<T>::rename(name); 
  if ((fp = fopen(ImageT<T>::fn_img.c_str(), "wb")) == NULL)
    REPORT_ERROR(1503,(string)"ImageXmipp::write: File "+ImageT<T>::fn_img +
       " cannot be written");
  adjust_header();
  bool reversed=(force_reversed) ? !header.reversed():header.reversed();
  header.write(fp, force_reversed);
  ImageT<T>::write(fp, reversed, IFLOAT);
  fclose(fp);  
}

/* Is Xmipp image? --------------------------------------------------------- */
int Is_ImageXmipp(const FileName &fn, bool skip_type_check,
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
int Is_FourierImageXmipp(const FileName &fn, bool skip_type_check,
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

/* Get Image size ---------------------------------------------------------- */
void GetXmippImageSize(const FileName &fn, int &Ydim, int &Xdim) {
   FILE *fp; 
   int result;
   headerXmipp header(headerXmipp::IMG_XMIPP);

   // Open file
   if ((fp = fopen(fn.c_str(), "rb")) == NULL)
     REPORT_ERROR(1501,"Is_ImageXmipp: File " + fn + " not found");

   // Read header
   result = header.read(fp, FALSE, FALSE);

   fclose(fp);
   if (result) {
      Ydim=header.iYdim();
      Xdim=header.iXdim();
   } else {
      Ydim=Xdim=-1;
   }
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
   float fIform;
   FILE *fh;
   
   img.read(name,fIform,Ydim,Xdim,reversed,image_type);        
   img.read(fh  ,fIform,Ydim,Xdim,reversed,image_type);
   img.LoadImage(name);
   img.write(name,reversed,image_type);     
   img.write(fh,  reversed,image_type);
   imgx.read(name,skip_type_check,reversed);
   imgx.write(name,force_reversed);
   imgx.get_transformation_matrix(reversed);
   imgx.check_oldxmipp_header();


}

void instantiate_images() {
//   short           s; instantiate_images(s);
//   int             i; instantiate_images(i);
   float           f; instantiate_images(f);
   double          d; instantiate_images(d);
   complex<double> c; instantiate_images(c);
}
