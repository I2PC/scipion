/***************************************************************************
 *
 * Authors: Pedro Antonio de Alarcón (pedro@cnb.uam.es)
 *          Carlos Oscar S. Sorzano
 *          Alberto Pascual Montano
 *          Roberto Marabini
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

#include "../xmippVolumes.hh"

/* Input (read) from file specifying its dimensions ------------------------ */
template <class T> void VolumeT<T>::read(FileName name,
  int Zdim, int Ydim, int Xdim, bool reversed, Volume_Type volume_type) _THROW {
  FILE *fh;
  clear(); 
  fn_img=name;
  if ((fh = fopen(fn_img.c_str(), "rb")) == NULL)
    REPORT_ERROR(1501,"Volume::read: File " + fn_img + " not found");
  read(fh, Zdim, Ydim, Xdim, reversed, volume_type);

  fclose(fh);
}
  
template <class T> void VolumeT<T>::read(FILE *fh,
  int Zdim, int Ydim, int Xdim, bool reversed, Volume_Type volume_type) {
  img.resize(Zdim, Ydim,Xdim);
  FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(img)
      switch (volume_type) {
         case VBYTE:
            unsigned char u;
            FREAD (&u, sizeof(unsigned char), 1, fh, reversed);
            MULTIDIM_ELEM(img,i)=(T)u;
            break;
         case VINT:
            int ii;
            FREAD (&i, sizeof(int), 1, fh, reversed);
            MULTIDIM_ELEM(img,i)=(T)ii;
            break;
            //NOTE integers and floats need to be reversed in identical way
         case VFLOAT:
            float f;
            FREAD (&f, sizeof(float), 1, fh, reversed);
            MULTIDIM_ELEM(img,i)=(T)f;
            break;
      }
}


void VolumeT<complex<double> >::read(FILE *fh,
  int Zdim, int Ydim, int Xdim, bool reversed, Volume_Type volume_type)
{
  img.resize(Zdim, Ydim,Xdim);
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
template <class T> void VolumeT<T>::write(FileName name, bool reversed,
   Volume_Type volume_type) _THROW {
   FILE *fp;
   if (name != "") rename(name);  

   if ((fp = fopen(fn_img.c_str(), "wb")) == NULL) {
     REPORT_ERROR(1503,"Volume::write: File " + fn_img + " cannot be saved");
   };
   write(fp, reversed, volume_type);
   fclose(fp);  
}

template <class T> void VolumeT<T>::write(FILE *fh, bool reversed,
   Volume_Type volume_type) {
   if (XSIZE(img)==0 || YSIZE(img)==0 || ZSIZE(img)==0) return;
   FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(img)
       switch (volume_type) {
          case VBYTE:
             unsigned char u;
             u=(unsigned char) MULTIDIM_ELEM(img,i);
             FWRITE (&u, sizeof(unsigned char), 1, fh, reversed);
             break;
          case VFLOAT:
             float f;
             f=(float) MULTIDIM_ELEM(img,i);
             FWRITE (&f, sizeof(float), 1, fh, reversed);
             break;
          case VINT:
             int ii;
             ii=(int) MULTIDIM_ELEM(img,i);
             FWRITE (&ii, sizeof(int), 1, fh, reversed);
             break;
       }
}

// Specific function to write volumes with complex numbers in them
void VolumeT<complex<double> >::write(FILE *fh, bool reversed,
   Volume_Type volume_type) 
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
template <class T> void VolumeXmippT<T>::read(const FileName &name,
  bool skip_type_check, bool force_reversed) _THROW {
  FILE *fp; 

  rename(name); 
  if ((fp = fopen(fn_img.c_str(), "rb")) == NULL)
    REPORT_ERROR(1501,(string)"VolumeXmipp::read: File "+fn_img+" not found");

  // Read header
  if (!header.read(fp, skip_type_check, force_reversed))
     REPORT_ERROR(1502,"VolumeXmipp::read: File " + fn_img +
        " is not a valid Xmipp file");   

  // Read whole image and close file
  VolumeT<T>::read(fp, header.iSlices(), header.iYdim(), header.iXdim(),
     header.reversed(), VFLOAT);
  fclose(fp);

  header.set_header();  // Set header in a Xmipp consistent state
}

/* Write Xmipp Volume  ------------------------------------------------------- */
template <class T> void VolumeXmippT<T>::write(const FileName &name,
  bool force_reversed) _THROW {
  FILE *fp;
  if (name != "") rename(name); 
  if ((fp = fopen(fn_img.c_str(), "wb")) == NULL)
    REPORT_ERROR(1503,(string)"VolumeXmipp::write: File "+fn_img +
       " cannot be written");
  adjust_header();  
  header.write(fp, force_reversed);
  VolumeT<T>::write(fp, header.reversed(), VFLOAT);
  fclose(fp);  
}

/* Is Xmipp image? --------------------------------------------------------- */
int Is_VolumeXmipp(FileName fn, bool skip_type_check, bool force_reversed)
   _THROW {
   FILE *fp; 
   int result=0;
   headerXmipp header(headerXmipp::VOL_XMIPP);
   headerXmipp header1(headerXmipp::VOL_INT);

   // Open file
   if ((fp = fopen(fn.c_str(), "rb")) == NULL)
        REPORT_ERROR(1501,"Is_VolumeXmipp: File " + fn + " not found");

   // Read header
   result = header.read(fp, skip_type_check, force_reversed);
   if(result==1) {fclose(fp); return (headerXmipp::VOL_XMIPP);}


   result = header1.read(fp, skip_type_check, force_reversed);
   fclose(fp);
   return (headerXmipp::VOL_INT * result);
}

/* Is Xmipp image? --------------------------------------------------------- */
int Is_FourierVolumeXmipp(FileName fn, bool skip_type_check, bool force_reversed)
   _THROW {
   FILE *fp; 
   int result=0;
   headerXmipp header(headerXmipp::VOL_FOURIER);

   // Open file
   if ((fp = fopen(fn.c_str(), "rb")) == NULL)
        REPORT_ERROR(1501,"Is_FourierVolumeXmipp: File " + fn + " not found");

   // Read header
   result = header.read(fp, skip_type_check, force_reversed);
   fclose(fp);
   return result;
}

template <class T>
   void instantiate_Volume (VolumeT<T> v) {
   VolumeT<int>    Vint;
   VolumeT<double> Vdouble;
   
   VolumeT<T>      Va;             // Empty constructor
   VolumeT<T>      Vb(1,1,1);      // Constructor with size
   VolumeT<T>      Vc(Vint);       // Copy constructor from int.
   VolumeT<T>      Vcc(Vdouble);    // Copy constructor from double.
   matrix3D<double>         r;
   FileName        aux_fn_img;
   FILE *fh;
   VolumeT<T>      Vd(aux_fn_img); // Constructor using filename.
   
   Va=Vc; // Assignment from int
   Va=r;  // Assignment from matrix3D.
   Vb.rename(aux_fn_img);
   Vb.clear();
   Vb.move_origin_to_center();
   Vb.adapt_to_size(1,1,1);
   Vb().resize(1,1,1);// Matrix access
   Vb(0,0,0)=Va(0,0,0);// Voxel access.
   Vb.name();//Name access
   cout << Va; //
   Vb.read(aux_fn_img, 1, 1, 1, (bool) 0,  (Volume_Type) 0);
   Vb.write(aux_fn_img);// the other read and write are called by this ones
}

template <class T>
   void instantiate_VolumeXmipp (VolumeXmippT<T> Vx) {
   VolumeXmippT<T> VXa;        // Empty constructorx
   VolumeXmippT<T> VXb(1,1,1); //  Constructor with size.
   VolumeXmippT<T> VX("1");    // Constructor with filename, read from disk.
   FileName        aux_fn_img;
   VolumeXmippT<int>    VXint;
   VolumeXmippT<double> VXdouble;
   VolumeXmippT<int>    Vc;    
   VolumeXmippT<double> Vcc;   
   VolumeXmippT<T>      Vd(Vx);  // Copy constructor
   VolumeT<int>         Vint;
   VolumeT<double>      Vdouble;

   VXa.clear(); 
   cout << VXa;
   VXa =VXa;
   VolumeT<double> V; VXb=V;
   matrix3D<float> m; VXb=m;
// I do not know what to do with this
//   VXa.assign_from(Vint);
//   VXa.assign_from(Vdouble);
   VXa.read(aux_fn_img);
   VXa.write(aux_fn_img);
   VXa.adjust_header();
   VXa.clear_header();
   VXa.rename("");
   VXa.reversed();
}

void instantiateVolume() {
   VolumeT<double>           V1;instantiate_Volume(V1);
   VolumeT<int>              V2;instantiate_Volume(V2);
   VolumeT<complex<double> > V3;instantiate_Volume(V3);
}

void instantiateVolumeXmipp() {
   VolumeXmippT<double>      V1; instantiate_VolumeXmipp(V1);
   VolumeXmippT<int>         V2; instantiate_VolumeXmipp(V2);
   VolumeXmippT<complex<double> > V3; instantiate_VolumeXmipp(V3);
}
