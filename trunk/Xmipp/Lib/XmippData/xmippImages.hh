/***************************************************************************
 *
 * Authors: Alberto Pascual Montano (pascual@cnb.uam.es)
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

#ifndef _XMIPPIMAGES_H
   #define _XMIPPIMAGES_H

/* ************************************************************************* */
/* INCLUDES                                                                  */
/* ************************************************************************* */

#include "xmippFuncs.hh"
#include "xmippMatrices2D.hh"
#include "xmippHeader.hh"
#include <typeinfo> /* because typeid is used */
#include <complex>
/* ************************************************************************* */
/* FORWARD DEFINITIONS                                                       */
/* ************************************************************************* */
// This forward definitions are needed for defining operators functions that
// use other clases type




/* ************************************************************************* */
/* IMAGES                                                                    */
/* ************************************************************************* */
/**@name Images*/
//@{
typedef enum {IBYTE=1, IFLOAT=2, I16=3} Image_Type;

/** Basic Image Class.
    The image class is a general class which only contains 
    the image itself and a filename for it. It has got a "float" matrix2D as
    member, and basically all operations between images are based on that
    class.

    This class is the usual one when you want to operate images in memory.
    Images belonging to this class cannot be saved, instead you could
    use the class ImageXmipp which inherits from this class, so you have
    the whole funcitonality of the class Image plus the possibility
    of saving at the end.

    See \Ref{Image logical and physical access} for a detailed information
    about the pixel accessing, and the conventions used in the image
    definition. 
*/

template <class T> class ImageT {
protected:
   FileName          fn_img;              // name of the image
public:
   matrix2D<T>  img;                 // matrix with the image
public:
   // Constructors/Destructor ..............................................
   /**@name Constructors*/
   //@{
   /** Empty constructor.
       An empty image with size 0x0 is created.
       \\Ex: Image I;*/
   ImageT() {fn_img = "";}

   /** Constructor with size.
       A blank image (0.0 filled) is created with the given size. Pay
       attention to the dimension order: Y and then X.
       \\Ex: Image I(64,64);*/
   ImageT(int Ydim, int Xdim){
     img.resize(Ydim,Xdim);
     fn_img="";
   };            
   
   /** Constructor using filename.
       The name for the image is assigned but no load is performed.
       The image content is empty, ie, with size 0x0.
       Ex: Image I("g0ta0001"); */
   ImageT(FileName _name){fn_img=_name;};
   
   /** Copy constructor.
       \\Ex: Image I2(I1); */
   ImageT(const ImageT &I){
     img=I.img;
     fn_img=I.fn_img;
   };
   //@}
   
   /**@name Some operations*/
   //@{
   /** Assignment.*/
   ImageT& operator = (const ImageT &I)
      {if (this!=&I) {fn_img=I.fn_img;img=I.img;}
       return *this;}

   /** Assignment from matrix.*/
   template <class T1>
      ImageT& operator = (const matrix2D<T1> &m)
         {if (&img!=(matrix2D<T> *) &m) {fn_img=""; type_cast(m,img);}
          return *this;}

   /** Rename image.
       Give a new name to the image.
       \\ Ex: I.rename("new_name"); */
   virtual void rename(FileName newName)   {   fn_img = newName;   }

   /** Empty image.
       This function clears the image to a 0x0 image without name.
       \\ Ex: I.clear(); */
   virtual void clear() {fn_img = ""; img.clear();}

   /** Sets the origin of the image at its center.
       The exact formula for the center of an image is defined in
       the Conventions Section, this function modify the indexes of
       the image in such a way that now the logical origin is at its
       center and the indexes may take negative values. Use startingX,
       startingY, finishingX and finishingY to setup loops.
       \\ Ex: I.move_origin_to_center();
       @see startingX, finishingX */
   void move_origin_to_center() {
      img.startingY()=FIRST_XMIPP_INDEX(img.ydim);
      img.startingX()=FIRST_XMIPP_INDEX(img.xdim);
   }

   /** Fill with 0 and move origin to center.
       This function resizes the image to the given size, fill it with
       0.0 and then moves the image logical origin to its center. See
       previous function.
       \\Ex: I.adapt_to_size(64,64)
       @see move_origin_to_center */
   void adapt_to_size(int Ydim, int Xdim)
      {img.init_zeros(Ydim,Xdim); move_origin_to_center();}
   //@}
   
   /**@name Image access*/
   //@{
   /** 2D Matrix access.
       This operator can be used to access the matrix, and the
       matrix operations defined in matrix2D. In this way we could
       resize an image just by resizing its associated matrix or
       we could add two images by adding their matrices.
       \\Ex: I().resize(128,128);
       \\Ex: I2()=I1()+I2(); */
   virtual matrix2D<T>&  operator () () {return img;}
   virtual const matrix2D<T>&  operator () () const {return img;}

   /** Pixel access.
       This operator is used to access a pixel within the image.
       This is a logical access, so you could access to negative
       positions if the image has been defined so (see the general
       explanation for the class).
       \\Ex: cout << "Grey level of pixel (-3,-3) of the image = "
                  << I(-3,-3) << endl;
       \\Ex: I(-3,-3)=I(-3,-2); */
   T& operator () (int i, int j) const {return img(i,j);}

   /** Name access.
       This function is used to know the name of the image. It
       cannot be used to assign a new one.
       \\ Ex: cout << "Image name: " << I.name() << endl; */
   const FileName name() const {return fn_img;}

   /** Cout << Image.
       Shows name and size */
   friend ostream & operator << (ostream &out, const ImageT &I)
      {out << "Image Name   : " << I.fn_img << endl
           << "dimensions   : " << I.img.RowNo() << " x " << I.img.ColNo()
           << "  (rows x columns)" << endl;
      return out;}
   //@}
   
   /**@name I/O functions
      All these functions work with the image written in raw floats.*/
   //@{
   /** Read image from disk, given the image's dimensions.
       If the image doesn't exist at the given path then an exception is
       thrown.
       
       The reversed flag indicates if the elements in element_size
       must be read in a reversed way.
       
       Elements are supposed to be in the following order (y,x)=
       (0,0)(0,1)(0,2), ..., (0,Xdim-1), (1,0), (1,1), ...
       
       The element size can be adjusted so that raw images of bytes (IBYTE),
       unsigned ints woth 16 bits (I16) and floats (IFLOAT) can be read. 
       \\ Ex: I.read(65,65,"g0ta0001.raw");*/
   bool read(const FileName &name, float fIform, int Ydim, int Xdim,
      bool reversed=FALSE, Image_Type image_type=IBYTE) _THROW;

   /** Read image from disk using a file pointer.
       This is the core routine of the previous one. */
   bool read(FILE * &fh, float fIform, int Ydim, int Xdim, bool reversed, Image_Type image_type);

  // LoadImage is a static function that loads an image file depending on
  // what kind of image type (e.g., Xmipp, Imagic) it is; it returns a
  // pointer to the base Image class.  Caller is responsible for deleting
  // the memory for the object.  Returns NULL on error.
  static ImageT<T> *LoadImage (FileName name, bool apply_geo=FALSE);

   /** Write image to disk.
       If there is any problem in the writing, an exception is thrown.
       You can give a name to the written volume different from the one used
       when it was read. From this point the filename of the volume has
       changed. This is somehow like the "Save as ..." and "Save".
       Valid types are IBYTE, I16, IFLOAT.
       \\ Ex: I.write() ---> Save
       \\ Ex: I.write("g0ta0001.raw") ---> Save as */
   void write(FileName name = "", bool reversed=FALSE,
      Image_Type image_type=IBYTE) _THROW;

   /** Write image to disk using a file pointer.
       This is the core routine of the previous one. */
   void write(FILE * &fh, bool reversed, Image_Type image_type);
   //@}
};

/**@name Image Speed up macros*/
//@{
/** Matrix access.
    This macro does the same as the normal matrix access but
    in a faster way as no function call is generated. This macro can be
    used with any of the derived classes from the Image class.
    \\ Ex: IMGMATRIX(I).resize(128,128);
    \\ Ex: IMGMATRIX(I2)=IMGMATRIX(I1)+IMGMATRIX(I2); */
#define IMGMATRIX(I) ((I).img)

/** Array access.
    This macro allows you to access to the bidimensional array behind
    the image (double *). */
#define IMGARRAY(I) MAT_ARRAY((I).img)

/** Pixel access.
    This macro does the same as the normal pixel access (remember,
    logical access) but in a faster way as no function call is
    generated. This macro can be
    used with any of the derived classes from the Image class,
    except for the Oversampled Images (see \Ref{ImageOver}) as their 
    logical positions are related in a special way to the physical ones.
    \\Ex: cout << "Grey level of pixel (-3,-3) of the image = "
               << IMGPIXEL(I,-3,-3) << endl;
    \\Ex: IMGPIXEL(I,-3,-3)=IMGPIXEL(I,-3,-2); */
#define IMGPIXEL(I,i,j) MAT_ELEM(((I).img),(i),(j))

/** Physical pixel access.
    The physical pixel access gives you access to a pixel by
    its physical position and not by its logical one. This access
    shouldn't be used as a custom, use instead the logical access,
    but there might be cases in which this access might be interesting.
    Physical positions start at index 0 in C.
    This macro can be used by any of the derived classes from the
    Image class.
    \\Ex: cout << "This is the first pixel stored in the image " <<
          DIRECT_IMGPIXEL(I,0,0) << endl;
    @see DIRECT_MAT_ELEM */
#define DIRECT_IMGPIXEL(I,i,j) DIRECT_MAT_ELEM(((I).img),(i),(j))

//@}

/* ************************************************************************* */
/* XMIPP IMAGES                                                                    */
/* ************************************************************************* */
/** Xmipp 2D Images.
    The Xmipp image is a normal image (inherited from Image class) plus a
    Spider header. This is the appropiate class to work with images in
    memory which we want to save later.
    The data in the header is not directly accesible from
    the programmer and must be set through object functions, in this way
    the coherence in the header is assured. See File Formats for
    more information about the Spider format.
    
    In principle, the image starts at pixel (0,0) but this can be modified
    for any other logical access. See class \Ref{Image} for more information.
    @see Image
    @see Euler Angles
*/

template <class T> class ImageXmippT: public ImageT<T> {
protected:
   headerXmipp header;         		               // Declares a header

public:
   // Constructors .........................................................
   /**@name Constructors*/
   //@{
   /** Empty constructor.
       Creates an empty (0x0) image with no information in the header.
       \\ Ex: ImageXmipp IX;*/
   ImageXmippT():ImageT<T>() {
     if (typeid(T)==typeid(complex<double>))
        header.headerType() = headerXmipp::IMG_FOURIER; // Sets header of type Image_Fourier (Complex) 
     else if (typeid(T)==typeid(double))
        header.headerType() = headerXmipp::IMG_XMIPP;  // Sets header of type Image_XMipp
     header.Slices() = 1;                 	       // Sets header Slices 
   }; 			       
                                           
   /** Constructor with size.
       Creates a 0.0 filled image of size Ydim x Xdim.
       \\ Ex: ImageXmipp IX(64,64); */
   ImageXmippT (int Ydim, int Xdim):ImageT<T>(Ydim, Xdim) {
     if (typeid(T)==typeid(complex<double>))
        header.headerType() = headerXmipp::IMG_FOURIER; // Sets header of type Image_Fourier (Complex) 
     else if (typeid(T)==typeid(double))
        header.headerType() = headerXmipp::IMG_XMIPP;  // Sets header of type Image_XMipp
     header.set_dimension(Ydim, Xdim);                 // Sets header dimensions     
     header.Slices() = 1;                 	       // Sets header Slices 
     header.set_header(); 			       // Initialize header
     header.set_time();				       // Set time and date	
     header.set_date();
   };

   /** Constructor with filename, read from disk.
       The filename given must exist, then the file is loaded in the
       ImageXmipp class structure. You have loaded the image at the
       declaration time.
       \\ Ex: ImageXmipp IX("g1ta0001.spd"); */
   ImageXmippT (FileName _name, bool apply_geo=FALSE):ImageT<T>(_name){
     if (typeid(T)==typeid(complex<double>))
        header.headerType() = headerXmipp::IMG_FOURIER; // Sets header of type Image_Fourier (Complex) 
     else if (typeid(T)==typeid(double))
        header.headerType() = headerXmipp::IMG_XMIPP;  // Sets header of type Image_XMipp
        read(_name,FALSE,FALSE,apply_geo,FALSE);       // Read image from file
   }
   
   /** Copy constructor.
       \\ Ex: ImageXmipp IX2(IX1); */
   ImageXmippT (const ImageXmippT<T> &I): ImageT<T>(I) {header = I.header;}

   /** Empty image.
       All information is cleared.
       \\Ex: IX.clear();*/
   void clear() {clear_header(); ImageT<T>::clear();}
   //@}
   
   // Overload some operators ..............................................
   /**@name Some operators*/
   //@{
   /** Show the header information of a Xmipp image.
       \\ Ex: cout << IX; */
   friend ostream& operator << (ostream& out, const ImageXmippT<T> &I)
      {out << (ImageT<T> &) I << I.header; return out;}
   
   /** Assignment from another Xmipp image.
       \\ Ex: ImageXmipp IX1, IX2; IX2=IX1; */
   ImageXmippT<T>& operator= (const ImageXmippT<T> &op1)
      {if (&op1!=this) {
          this->ImageT<T>::operator = (op1);
          header = op1.header;
       }
       return *this;}
   
   /** Assignment from a generic image.
       The Euler angles are set to 0.
       \\Ex: Image I; ImageXmipp IX; IX=I;*/
   ImageXmippT<T>& operator= (const ImageT<T> &op1) {
      if (this!=&op1) {
         this->ImageT<T>::operator = (op1);
         clear_header(); adjust_header();
      }
      return *this;
   }

   /** Assignment from a matrix.
       The Euler angles are set to 0 and the image filename is set to "".
       \\Ex: matrix2D<int> m; ImageXmipp IX; IX=m; */
   template <class T1>
      ImageXmippT<T>& operator=(const matrix2D<T1> &op1) {
         if (&img!=(matrix2D<T> *) &op1) {
            this->ImageT<T>::operator = (op1);
            clear_header(); adjust_header();
         }
         return *this;
      }
   
   /** Assignment from any kind of image.
       \\Ex:ImageOver IO; ImageXmipp IX; IX.assign_from(&IO); */
   void assign_from(ImageT<T> *v) {*this=*v;}
   //@}

   // Input/Output .........................................................
   /**@name Input/Output*/
   //@{
   /** Read Xmipp image from disk.
       If the image doesn't exist at the given path then an exception is
       thrown.
       
       If only_apply_shifts is TRUE, then only the shifts of the header are
       applied (even if apply_geo is FALSE).
       If apply_geo is TRUE and only_apply_shifts is FALSE, then shifts
       and in-plane rotation are applied. 
       
       If skip_type_check is FALSE, then the routine automatically checks
       the endianness of thee file, and reversed should be set to FALSE.
       
       If skip_type_check is TRUE, then the programmer must know whether
       the endian is reversed or not, and provide this information in reversed.
       
       \\ Ex: IX.read("g1ta0002.spd");*/
   virtual bool read(const FileName &_name, bool skip_type_check=FALSE,
      bool reversed=FALSE, bool apply_geo=FALSE,
      bool only_apply_shifts=FALSE) _THROW;
   
   /** Write Xmipp image to disk.
       If there is any problem in the writing, an exception is thrown.
       You can give a name to the written image different from the one used
       when it was read. From this point the filename of the image has
       changed. This is somehow like the "Save as ..." and "Save".
       \\ Ex: IX.write() ---> Save
       \\ Ex: IX.write("g1ta0002.spd") ---> Save as */
   virtual void write(const FileName &_name = "", bool force_reversed=FALSE)
      _THROW;
   //@}

   // Header operations interface ..........................................
   // These are the only access to the header allowed 
   /**@name Header access*/
   //@{
   /** Resets header. */
   void clear_header() {header.clear();}

   /** Adjust header.
       Force header to have the dimensions of the image, time, date updated */
   void adjust_header() {
      if (typeid(T)==typeid(complex<double>))
         header.headerType() = headerXmipp::IMG_FOURIER; // Sets header of type Image_Fourier (Complex) 
      else if (typeid(T)==typeid(double))
         header.headerType() = headerXmipp::IMG_XMIPP;  // Sets header of type Image_XMipp
      header.set_dimension(YSIZE(img), XSIZE(img)); // Sets header dimensions     
      header.Slices() = 1;                 	    // Sets header Slices 
      header.set_time();			    // Set time and date 
      header.set_date();
      header.set_title(fn_img);                     // Set title
      header.set_header(); 			    // Initialize header
   }

   /** Change Filename.
       \\Ex: IX.rename("newName.spd"); */
   void rename(const FileName &newName) 
      {ImageT<T>::rename(newName); header.set_title(newName);}


   /** Reversed status.
       This is used for the little/big endian process. */
   bool reversed() const {return header.reversed();}

  /** Get geometric transformation matrix from 2D-image header.  */
  matrix2D<double> get_transformation_matrix(bool only_apply_shifts=FALSE);

   /** Check OldXmipp header location for image translation offsets
       If image appears to be in OldXmipp format, copy offsets to NewXmipp header location */
   void check_oldxmipp_header();

   /**@name Origin Offsets
      Origin Offsets Xoff and Yoff are in pixels and do not have to be integer values.*/
   /**@name As a block*/
   //@{
   /** Set origin Offsets 
       \\Ex: IX.set_originOffsets(1.51,-3.43); */
      void set_originOffsets(float _Xoff, float _Yoff)
        {header.set_originOffsets( _Xoff, _Yoff);}

   /** Get origin Offsets 
       \\Ex: IX.get_originOffsets(Xoff,Yoff); */
      void get_originOffsets(float &_Xoff, float &_Yoff) const
        {header.get_originOffsets(_Xoff, _Yoff);}
   //@}
   /**@name Independently*/
   //@{
   /**@name origin Offset in X-direction */
   //@{
   /** Set Xoff
       \\Ex: IX.Xoff()=3.50;*/
   float& Xoff() {return header.fXoff();}
   
   /** Get Xoff.
       \\Ex: cout << "Origin Offset in X-direction: " << IX.Xoff() << endl;*/
   float  Xoff() const {return header.fXoff();}
   //@}
   /**@name origin Offset in Y-direction */
   //@{
   /** Set Yoff
       \\Ex: IX.Yoff()=3.50;*/
   float& Yoff() {return header.fYoff();}
   
   /** Get Yoff.
       \\Ex: cout << "Origin Offset in Y-direction: " << IX.Yoff() << endl;*/
   float  Yoff() const {return header.fYoff();}
   //@}
   //@}

   /**@name Euler angles
       The angles must be in degrees and can be either positive or negative.
       \\Phi is the first Euler angle. A rotation around Z.
       \\Theta is the second Euler angle. A rotation around Y.
       \\Psi is the third Euler angle. A rotation around Z again.
       \\See the Angle Convention for more information about the Euler angles.*/
   //@{
   /**@name As a block*/
   //@{
   /** Set angles. (default position)
       \\Ex: IX.set_eulerAngles(30,-10,350); */
   template <class T1>
      void set_eulerAngles(T1 _Phi, T1 _Theta, T1 _Psi)
        {header.set_eulerAngles((float) _Phi, (float) _Theta, (float) _Psi);}
   
   /** Set angles. Spider has three places in which Euler angles
   can be stored. First alternative position
       \\Ex: IX.set_eulerAngles1(30,-10,350); */
   template <class T1>
      void set_eulerAngles1(T1 _Phi1, T1 _Theta1, T1 _Psi1)
        {header.set_eulerAngles1((float) _Phi1, (float) _Theta1, (float) _Psi1);}
   
   /** Set angles.Spider has three places in which Euler angles
   can be stored. Second alternative position
       \\Ex: IX.set_eulerAngles2(30,-10,350); */
   template <class T1>
      void set_eulerAngles2(T1 _Phi2, T1 _Theta2, T1 _Psi2)
        {header.set_eulerAngles2((float) _Phi2, (float) _Theta2, (float) _Psi2);}
   /** Clears fFlag flag. The number of triads of Euler angles stored in the header (up to three) is stored here. set_eulerAngles2 makes fFlag=2, set_eulerAngles1 makes fFlag=max(fFlag,1), set_eulerAngles does not change fFlag
        \\Ex: IX.clear_fFlag_flag(); */
      void clear_fFlag_flag()
        {header.clear_fFlag_flag();}
   
   /** Get angles. 
       \\Ex: IX.get_eulerAngles(phi,theta,psi); */
   template <class T1>
      void get_eulerAngles(T1 &_Phi, T1 &_Theta, T1 &_Psi) const
        {header.get_eulerAngles(_Phi, _Theta, _Psi);}
   
   /** Get angles. Spider has three places in which Euler angles
   can be stored. First alternative position
       \\Ex: IX.get_eulerAngles1(phi,theta,psi); */
   template <class T1>
      void get_eulerAngles1(T1 &_Phi1, T1 &_Theta1, T1 &_Psi1) const 
        {header.get_eulerAngles1(_Phi1, _Theta1, _Psi1);}
   
   /** Get angles. Spider has three places in which Euler angles
   can be stored. Second alternative position
       \\Ex: IX.get_eulerAngles2(phi,theta,psi); */
   template <class T1>
      void get_eulerAngles2(float &_Phi2, float &_Theta2, float &_Psi2) const
        {header.get_eulerAngles2(_Phi2, _Theta2, _Psi2);}

   /** How many Euler angles are set?*/
   float Is_flag_set(void) { return(header.Is_flag_set());}
   //@}
   
   /**@name Independently*/
   //@{
   /**@name Old Rotational angle (fAngle1) */
   //@{
   /** Set old rotational angle.
       \\Ex: IX.old_rot()=30;*/
   float& old_rot() {return header.old_rot();}
   
   /** Get old rotational angle.
       \\Ex: cout << "First Euler angle " << IX.old_rot() << endl;*/
   float old_rot() const {return header.old_rot();}
   //@}

   /**@name Rotational angle */
   //@{
   /** Set Phi.
       \\Ex: IX.Phi()=30;*/
   float& Phi() {return header.Phi();}
   
   /** Get Phi.
       \\Ex: cout << "First Euler angle " << IX.Phi() << endl;*/
   float  Phi() const {return header.Phi();}

   /** Set Rotational angle.
       \\Ex: IX.rot()=30;*/
   float& rot() {return header.Phi();}
   
   /** Get rot.
       \\Ex: cout << "First Euler angle " << IX.rot() << endl;*/
   float  rot() const {return header.Phi();}


   /** Set Phi1. First alternative phi angle
       \\Ex: IX.Phi1()=30;*/
   float& Phi1() {return header.Phi1();}
   
   /** Get Phi1. First alternative phi angle
       \\Ex: cout << "First Euler angle " << IX.Phi1() << endl;*/
   float  Phi1() const {return header.Phi1();}

   /** Set 1st Rotational angle. First alternative phi angle
       \\Ex: IX.rot1()=30;*/
   float& rot1() {return header.Phi1();}
   
   /** Get rot1. First alternative phi angle
       \\Ex: cout << "First Euler angle " << IX.rot1() << endl;*/
   float  rot1() const {return header.Phi1();}


   /** Set Phi2. First alternative phi angle
       \\Ex: IX.Phi()=30;*/
   float& Phi2() {return header.Phi2();}
   
   /** Get Phi2. Second alternative phi angle
       \\Ex: cout << "First Euler angle " << IX.Phi2() << endl;*/
   float  Phi2() const {return header.Phi2();}

   /** Set 2sd Rotational angle. First alternative phi angle
       \\Ex: IX.rot()=30;*/
   float& rot2() {return header.Phi2();}
   
   /** Get rot2. Second alternative phi angle
       \\Ex: cout << "First Euler angle " << IX.rot2() << endl;*/
   float  rot2() const {return header.Phi2();}
   //@}
   
   /**@name Tilting */
   //@{
   /** Set Theta.
       \\Ex: IX.Theta()=-10;*/
   float& Theta() {return header.Theta();}

   /** Get Theta.
       \\Ex: cout << "Second Euler angle " << IX.Theta() << endl;*/
   float  Theta() const {return header.Theta();}

   /** Set Tilting angle.
       \\Ex: IX.tilt()=-10;*/
   float& tilt() {return header.Theta();}

   /** Get Tilting angle.
       \\Ex: cout << "Second Euler angle " << IX.tilt() << endl;*/
   float  tilt() const {return header.Theta();}


   /** Set Theta1.
       \\Ex: IX.Theta1()=-10;*/
   float& Theta1() {return header.Theta1();}

   /** Get Theta1.
       \\Ex: cout << "Second Euler angle " << IX.Theta1() << endl;*/
   float  Theta1() const {return header.Theta1();}

   /** Set 1st Tilting angle.
       \\Ex: IX.tilt1()=-10;*/
   float& tilt1() {return header.Theta1();}

   /** Get 1st Tilting angle.
       \\Ex: cout << "Second Euler angle " << IX.tilt1() << endl;*/
   float  tilt1() const {return header.Theta1();}

   /** Set Theta2.
       \\Ex: IX.Theta2()=-10;*/
   float& Theta2() {return header.Theta2();}

   /** Get Theta2.
       \\Ex: cout << "Second Euler angle " << IX.Theta2() << endl;*/
   float  Theta2() const {return header.Theta2();}

   /** Set 2sd Tilting angle.
       \\Ex: IX.tilt2()=-10;*/
   float& tilt2() {return header.Theta2();}

   /** Get 2sd Tilting angle.
       \\Ex: cout << "Second Euler angle " << IX.tilt2() << endl;*/
   float  tilt2() const {return header.Theta2();}
   //@}

   /**@name Psi angle */
   //@{
   /** Set Psi.
       \\Ex: IX.Psi()=350;*/
   float& Psi() {return header.Psi();}
   /** Get Psi.
       \\Ex: cout << "Third Euler angle " << IX.psi() << endl;*/
   float  Psi() const {return header.Psi();}
   /** Set psi.
       \\Ex: IX.psi()=350;*/
   float& psi() {return header.Psi();}
   /** Get psi.
       \\Ex: cout << "Third Euler angle " << IX.psi() << endl;*/
   float  psi() const {return header.Psi();}

   /** Set Psi1.
       \\Ex: IX.Psi1()=350;*/
   float& Psi1() {return header.Psi1();}
   /** Get Psi1.
       \\Ex: cout << "Third Euler angle " << IX.psi1() << endl;*/
   float  Psi1() const {return header.Psi1();}
   /** Set psi1.
       \\Ex: IX.psi1()=350;*/
   float& psi1() {return header.Psi1();}
   /** Get psi1.
       \\Ex: cout << "Third Euler angle " << IX.psi1() << endl;*/
   float  psi1() const {return header.Psi1();}

   /** Set Psi2.
       \\Ex: IX.Psi2()=350;*/
   float& Psi2() {return header.Psi2();}
   /** Get Psi2.
       \\Ex: cout << "Third Euler angle " << IX.psi2() << endl;*/
   float  Psi2() const {return header.Psi2();}
   /** Set psi2.
       \\Ex: IX.psi2()=350;*/
   float& psi2() {return header.Psi2();}
   /** Get psi2.
       \\Ex: cout << "Third Euler angle " << IX.psi2() << endl;*/
   float  psi2() const {return header.Psi2();}
   //@}
   //@}
   //@}
   //@}
         
};

/**@name Related functions */
//@{
/** True if the given volume is an Xmipp image. */
int Is_ImageXmipp(const FileName &fn, bool skip_type_check=FALSE,
      bool reversed=FALSE) _THROW;

/** True if the given volume is a Fourier Xmipp image. */
int Is_FourierImageXmipp(const FileName &fn, bool skip_type_check=FALSE,
      bool reversed=FALSE) _THROW;

/** Get size of an image.
    It returns -1 if the given image is not an Xmipp file*/
void GetXmippImageSize(const FileName &fn, int &Ydim, int &Xdim);
//@}

/* ************************************************************************* */
/* TYPE DEFINITIONS                                                          */
/* ************************************************************************* */
typedef ImageT<double> Image;
typedef ImageXmippT<double> ImageXmipp;
typedef ImageT< complex<double> > FourierImage;
typedef ImageXmippT< complex<double> > FourierImageXmipp;

/**@name Conversion between Image Types */
//@{
/** Converts a Xmipp Image (real numbers) into  a Fourier Xmipp Image (complex numbers)  */
void ImageXmipp_to_FourierImageXmipp(ImageXmipp &I,FourierImageXmipp &F);

/** Converts a Fourier Xmipp Image  into a Xmipp Image. 
    A simple copy of real parts is done */
void FourierImageXmipp_to_ImageXmipp(FourierImageXmipp &F,ImageXmipp &I);

//@}

/* ************************************************************************* */
/* OVERSAMPLED IMAGES                                                        */
/* ************************************************************************* */
/** Oversampled images.
    The oversampled images are images which allow a more accurate treatment
    of information by oversampling all pixels. The idea of this class
    is to have an image with a logical size smaller than its physical one,
    for this reason you could use non integer logical positions to access
    to different samples of the "same" pixel. Let's set an example,
    blob footprints are of size 5x5, for instance, with the center
    at physical position 3x3. It's convenient to have this footprint
    defined between -2 and 2 with the origin at the center of the image.
    But it's also convenient to have more than 25 values of the footprint,
    this could be done by sampling finer each pixel, say 51 samples per
    pixel, this would result in an image of 255x255 samples (with the
    center at [127][127]). But we don't to access to the pixel at position
    (-120,-120) but to the one at logical position (-2.35,-2.35) which
    we know is close to the border. The oversampled image class allows
    you this kind of image access. You have to say from where to where
    your image is defined for each axis (in this case -2...2 for both) and
    which are the sampling rates for both axes (51 for both). This is
    the initialisation process. From now on you can work with your
    image in the way formerly described.
    
    Pay attention to two points:
    \begin{enumerate}
    \item The sampling rate must be an odd number, this is because
    so the logical oversampled pixels (-2,-2), (-2, -1) ... are exactly
    on one cell of the underlying 2D matrix.
    \item As the oversampled pixels are centered with respect to the 
    "superpixels" defined by the 2D matrix, the extent of the oversampled
    image is a little larger than from (-2,-2) to (2,2), ie, from (-2.5,-2.5)
    to (2.5,2.5)
    \end{enumerate}
    
    Oversampled images are normal
    images except for pixel access which can be done in two fashions:
    either using the normal Image procedures (in this case, you are
    restricted to the macro IMGPIXEL), or using the fractional indexes
    defined for this class. Have
    a look at the following example of how to compute the blob footprint
    with this kind of images. Pay special attention to the pixel access
    at the end of the loop. Notice that x and y moves at values between
    -2 and 2, which keep logical meaning of what we are doing while
    all the burden of computing where to put these values at the image
    is done by the library.
    
    \begin{verbatim}
    void footprint_blob( 
       ImageOver &blobprint,         // blob foot_print table
       const struct blobtype &blob,  // blob description
       int   istep,                  // number of foot-print samples per one sample 
                                     // on projection plane in u,v directions
       int   normalise)              // set to 1 if you want normalise. Usually
                                     // you may omit it and no normalisation is performed
    {  
    // Resize output image and redefine the origin of it
       int footmax=CEIL(blob.radius);
       blobprint.init(-footmax, footmax, istep, -footmax, footmax, istep);

    // Run for indexes in the Image class
       for (int i = STARTINGY(blobprint()); i <= FINISHINGY(blobprint()); i++)
          for (int j = STARTINGX(blobprint()); j <= FINISHINGX(blobprint()); j++) {
             // Compute oversampled index, and then its value
             double vi, ui; IMG2OVER(blobprint,i,j,vi,ui);
             double r=sqrt(vi*vi + ui*ui);
             IMGPIXEL(blobprint,i,j)=blob_val(r,blob);
          }

    // Adjust the footprint structure
       if (normalise) blobprint()=blobprint()/blobprint().sum();
    }
    \end{verbatim}
    Note: for this class the X axis has been renamed as U, and Y as V.
*/
class ImageOver: public Image {
public:
   int   uistep, vistep;                  // number of samples per
                                          // one sample on normal image (50)
                                          // in u,v directions
   int   overumax, overvmax;              // table borders in normal units (-2,2)
   int   overumin, overvmin;              // table borders in normal units (-2,2)
                                          // They should be an integer number
public:
   /**@name Initialisation */
   //@{
   /** Prepare the Oversampled Image for work.
       This function defines the oversampled image, it's very
       important to call it before starting to work with the image.
       If the sampling rate in any of the directions is an even number, 
       then it is substituted by the next odd number (+1). This is
       done in order to force the logical oversampled pixels
       to be exactly in one cell of the underlying 2D matrix.
       \\ Ex: IO.init(-2,2,51,-2,2,51); */
   void init(int _vmin, int _vmax, int _vistep,
             int _umin, int _umax, int _uistep);

   /** Window.
       Set a window in the logical space. */
   void window(int _v0, int _u0, int _vF, int _uF);

   /** Empty image.
       A 0x0 image with no information about oversampling is produced.
       Remember to re-initialise the oversampled image before
       starting to use it again.
       \\ Ex: IO.clear();*/
   void clear();
   //@}

   /**@name Index conversions */
   //@{
   /** Returns the exact position of a non-integer location.
       The existing indexes (iu,iv) within the overimage are returned
       according to a certain (u,v) and the Oversampled Image definition.
       Usually this function is not needed in normal oversampled images
       operation as you can access to pixels with fractional indexes.
       In our example of footprints, this funciton would make
       the conversion between the oversampled position (-1,1) to the
       image index (-51,51) (what is returned). Don't be mislead by the
       physical position of the pixel (-51,51) which is (76,178).
       \\ Ex: IO.over2img(y,x,iy,ix); */
   void over2img(double v, double u,int &iv, int &iu) const _THROW
      { if (v<overvmin || v>overvmax)
           REPORT_ERROR(1505,"ImgeOver::over2img: v out of range");
        if (u<overumin || u>overumax)
           REPORT_ERROR(1505,"ImgeOver::over2img: u out of range");
        iu=(int)ROUND((((u)-overumin)*uistep));
        iv=(int)ROUND((((v)-overvmin)*vistep));}

   /** Speed up pixel index macro.
       This macro is the same as the function over2img but faster
       due to no function call is performed.
       \\ Ex: OVER2IMG(IO,y,x,iy,ix); */
   #define OVER2IMG(IO,v,u,iv,iu) \
       iu=(int)ROUND((((u)-(IO).overumin)*(IO).uistep));\
       iv=(int)ROUND((((v)-(IO).overvmin)*(IO).vistep));

   /** Returns the logical position of a "physical" location.
       This function makes the conversion from image logical
       location to oversampled logical position. In our example of footprints
       the conversion between the oversampled position (-51,51) to the
       image index (-1,1) (what is returned). Notice that the
       "physical" position of the pixel (-51,51) is not the true physical
       one, that would be (76,178).
       \\ Ex: IO.img2over(iv,iu,v,u) */
   void img2over(int iv, int iu, double &v, double &u) const _THROW
      { if (iu<0 || iu>XSIZE(img))
           REPORT_ERROR(1505,"ImageOver::img2over: iu out of range");
        if (iv<0 || iv>YSIZE(img))
           REPORT_ERROR(1505,"ImageOver::img2over: iv out of range");
        u=(double)overumin+iu/(double)uistep;
        v=(double)overvmin+iv/(double)vistep;}
       
   /** Speed up pixel index macro.
       This macro is the same as the function img2over but faster
       due to no function call is performed.
       \\ Ex: IMG2OVER(IO,iy,ix,y,x); */
   #define IMG2OVER(IO,iv,iu,v,u) \
       u=(double)(IO).overumin+(iu)/(double)((IO).uistep); \
       v=(double)(IO).overvmin+(iv)/(double)((IO).vistep);
   //@}

   /**@name Image access */
   //@{
   /** Constant pixel access with fractional indexes.
       This function returns you a constant pixel access referred with
       a fractional index. If you want to access to the pixels in the
       classic Image fashion then you should use the macro \Ref{IMGPIXEL}.
       An exception is thrown if you try to access an image outside the
       image.
       \\ Ex: cout << IO(1.34,-0.56) << endl; */
   double operator () (double v, double u) const _THROW
      { if (v<overvmin || v>overvmax)
           REPORT_ERROR(1505,"ImgeOver::over2img: v out of range");
        if (u<overumin || u>overumax)
           REPORT_ERROR(1505,"ImgeOver::over2img: u out of range");
       int iu, iv;
       OVER2IMG(*this,v,u,iv,iu);
       return MAT_ELEM(img,iv,iu);}

   /** Pixel access with fractional indexes.
       This function allows you to set pixels referred with
       a fractional index. If you want to access to the pixels in the
       classic Image fashion then you should use the macro \Ref{IMGPIXEL}.
       \\ Ex: IO(1.34,-0.56)=1; */
   double& operator () (double v, double u)
      {int iu, iv; OVER2IMG(*this,v,u,iv,iu); return MAT_ELEM(img,iv,iu);}
   
   /** Speed up pixel access macro.
       This macro is a speeded up version of the pixel access routines
       for oversampled images.
       \\ Ex:OVERPIXEL(IO,1.34,-0.56)=1; */
   #define OVERPIXEL(IO,y,x) IMGPIXEL((IO), \
       (int)ROUND(((u)*(IO).uistep)), \
       (int)ROUND(((v)*(IO).vistep)))

   // The following two functions have been redefined (they should be
   // inherited from Image) because the compiler doesn't admit this
   // inheritance and the redefinition of the fractional index pixel
   // access.
   /** Matrix access.
       This function allows you to access the matrix2D over which the
       complex data type is based.
       \\ Ex: cout << IO(); */
   matrix2D<double>& operator () () {return img;}
   const matrix2D<double>& operator () () const {return img;}

   //@}

   /**@name Size and shape information */
   //@{
   /** Maximum value in not sampled units on U axis.
       In our example: 2.
       \\ Ex: int Umax=IO.umax(); */
   int umax() const {return overumax;}

   /** Maximum value in not sampled units on V axis.
       In our example: 2.
       \\ Ex: int Vmax=IO.vmax(); */
   int vmax() const {return overvmax;}

   /** Minimum value in not sampled units on U axis.
       In our example: -2.
       \\ Ex: int Umin=IO.umin(); */
   int umin() const {return overumin;}

   /** Minimum value in not sampled units on V axis.
       In our example: -2.
       \\ Ex: int Vmin=IO.vmin(); */
   int vmin() const {return overvmin;}

   /** Sampling rate in U axis.
       In our example: 51.
       \\ Ex: int Usampling=IO.Ustep(); */
   int Ustep() const {return uistep;}

   /** Sampling rate in V axis.
       In our example: 51.
       \\ Ex: int Vsampling=IO.Vstep(); */
   int Vstep() const {return vistep;}
   //@}

   /**@name Tools */
   //@{
   /** Generate the normal image by averaging.
       This function passes from an oversampled image to a normal
       one, in our example from the 250x250 image to a 5x5. A
       pointer to the normal image must be given.
       \\ Ex: IO.down_sample(&I); */
   void downsample(Image *I) const;
   
   /** Generate an oversampled image from a normal one.
       This function passes from a normal image to an oversampled one,
       in our example from 5x5 to 250x250. First you have to initialise
       the parameters of the oversampled image and then the oversampling
       is done by bilinear interpolation.
       \\ Ex: oversample(&I)
       NOT IMPLEMENTED */
   void oversample(Image *I) const;
   //@}
};

//@Include: xmippImagic.hh

//@}
#endif
