/***************************************************************************
 *
 * Authors:      Alberto Pascual Montano (pascual@cnb.uam.es)
 *               Carlos Oscar Sanchez Sorzano
 *               Roberto Marabini
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

#ifndef _XMIPPVOLUMES_H
#define _XMIPPVOLUMES_H

/* ************************************************************************* */
/* INCLUDES                                                                  */
/* ************************************************************************* */

#include <iostream.h>
#include <stdio.h>
#include <typeinfo>
#include "xmippFuncs.hh"
#include "xmippMatrices3D.hh"
#include "xmippHeader.hh"

/* ************************************************************************* */
/* VOLUMES                                                                   */
/* ************************************************************************* */
/**@name Volumes*/
//@{
typedef enum {VBYTE=1, VFLOAT=2, VINT=3, VUCHAR=4, V16=5} Volume_Type;

/** Basic Volume Class.
    The volume class is a general class which only contains 
    the volume itself and a filename for it. It has got a float matrix3D as
    member, and basically all operations between volumes are based on that
    class.

    This class is the usual one when you want to operate volumes in memory.
    volumes belonging to this class cannot be saved, instead you could
    use the class volumeXmipp which inherits from this class, so you have
    the whole funcitonality of the class volume plus the possibility
    of saving at the end.

    See \Ref{Logical and physical access} for a detailed information
    about the voxel accessing, and the conventions used in the volume
    definition. 
*/
template <class T> class VolumeT {
public:
   FileName          fn_img;              // name of the image
   matrix3D<T>       img;                 // 3D matrix with the image
public:
   // Constructors/Destructor ..............................................
   /**@name Constructors*/
   //@{
   /** Empty constructor.
       An empty Volume with size 0x0x0 is created.
       \\Ex: Volume<double> V;*/
   VolumeT() {fn_img = "";}

   /** Constructor with size.
       A blank Volume (0.0 filled) is created with the given size. Pay
       attention to the dimension order: Z, Y and X.
       \\Ex: VolumeT<double> V(64,64,64);*/
   VolumeT(int Zdim, int Ydim, int Xdim){
     img.resize(Zdim,Ydim,Xdim);
     fn_img="";
   };            
   
   /** Constructor using filename.
       The name for the Volume is assigned but no load is performed.
       The Volume content is empty, ie, with size 0x0x0.
       Ex: VolumeT<double> V("art0001.vol"); */
   VolumeT(FileName _name) {fn_img=_name;}
   
   /** Copy constructor.
       \\Ex: VolumeT<double> V2(V1); */
   template <class Type>
      VolumeT(const VolumeT<Type> &I) {*this=I;}
   //@}
   
   /**@name Some operations*/
   //@{
   /** Assignment.
       \\ Ex: V2=V1; */
   template <class Type>
      VolumeT& operator = (const VolumeT<Type> &V)
      {if (this!=(VolumeT *) &V) {
          type_cast(V.img,img);
          fn_img=V.fn_img;
       }
       return *this;}

   /** Assignment from matrix3D.*/
   template <class Type>
      VolumeT& operator = (const matrix3D<Type> &m)
         {if (&img!=(matrix3D<T> *) &m) {fn_img=""; type_cast(m,img);}
          return *this;}

   /** Rename Volume.
       Give a new name to the Volume.
       \\ Ex: V.rename("new_name"); */
   virtual void rename(FileName newName) {fn_img = newName;}

   /** Empty Volume.
       This function clears the Volume to a 0x0x0 Volume without name.
       \\ Ex: V.clear(); */
   virtual void clear(){
     fn_img = "";
     img.clear();
   }

   /** Sets the origin of the Volume at its center.
       The exact formula for the center of a Volume is defined in
       the Conventions Section, this function modify the indexes of
       the Volume in such a way that now the logical origin is at its
       center and the indexes may take negative values. Use startingX,
       finishingX, ... to setup loops.
       \\ Ex: V.move_origin_to_center();
       @see startingX, finishingX */
   void move_origin_to_center() {
      img.startingZ()=FIRST_XMIPP_INDEX(img.zdim);
      img.startingY()=FIRST_XMIPP_INDEX(img.ydim);
      img.startingX()=FIRST_XMIPP_INDEX(img.xdim);
   }

   /** Fill with 0 and move origin to center.
       This function resizes the Volume to the given size, fill it with
       0.0 and then moves the Volume logical origin to its center. See
       previous function.
       \\Ex: V.adapt_to_size(64,64)
       @see move_origin_to_center */
   void adapt_to_size(int Zdim, int Ydim, int Xdim) {
      img.init_zeros(Zdim,Ydim,Xdim);
      move_origin_to_center();
   }
   //@}
   
   /**@name Image access*/
   //@{
   /** 3D Matrix access.
       This operator can be used to access the 3D matrix, and the
       matrix operations defined in matrix3D. In this way we could
       resize a Volume just by resizing its associated 3D matrix or
       we could add two Volumes by adding their 3D matrices.
       \\Ex: V().resize(128,128,128);
       \\Ex: V2()=V1()+V2(); */
   matrix3D<T>&  operator () () {return img;}
   const matrix3D<T>&  operator () () const {return img;}
   
   /** Voxel access.
       This operator is used to access a voxel within the Volume.
       This is a logical access, so you could access to negative
       positions if the Volume has been defined so (see the general
       explanation for the class).
       \\Ex: cout << "Grey level of voxel (2,-3,-3) of the Volume = "
                  << V(2,-3,-3) << endl;
       \\Ex: V(2,-3,-3)=V(2,-3,-2); */
   T& operator () (int z, int y, int x) const {return img(z,y,x);}


   /** Name access.
       This function is used to know the name of the Volume. It
       cannot be used to assign a new one.
       \\ Ex: cout << "Volume name: " << V.name() << endl; */
   FileName name() const {return (FileName)fn_img;}

   /** Cout << Volume.
       Shows name and size */
   friend ostream & operator << (ostream &out, const VolumeT<T> &V)
      {out << "Volume Name   : " << V.fn_img << endl
           << "dimensions   : " << V.img.SliNo() << " x " 
           << V.img.RowNo() << " x " << V.img.ColNo()
           << "  (slices x rows x columns)" << endl;
      return out;}
   //@}

   /**@name I/O functions
      All these functions work with the image written in raw floats.*/
   //@{
   /** Read Volume from disk, given the Volume's dimensions.
       If the Volume doesn't exist at the given path then an exception is
       thrown.
       
       The reversed flag indicates if the elements in element_size
       must be read in a reversed way.
       
       Elements are supposed to be in the following order (y,x)=
       (0,0)(0,1)(0,2), ..., (0,Xdim-1), (1,0), (1,1), ...
       
       The element size can be adjusted so that raw images of bytes (VBYTE),
       unsigned ints of 16 bits (V16) and floats (VFLOAT) can be read. 
       \\ Ex: V.read(65,65,65,"art0001.raw");*/

   void read(FileName name, int Zdim, int Ydim, int Xdim, bool reversed=FALSE,
      Volume_Type volume_type=VBYTE) _THROW;

   /** Read image from disk using a file pointer.
       This is the core routine of the previous one. */
   void read(FILE *fh, int Zdim, int Ydim, int Xdim, bool reversed,
      Volume_Type volume_type);

   /** Write Volume to disk.
       If there is any problem in the writing, an exception is thrown.
       You can give a name to the written Volume different from the one used
       when it was read. From this point the filename of the Volume has
       changed. This is somehow like the "Save as ..." and "Save".
       \\ Ex: V.write() ---> Save
       \\ Ex: V.write("art0002.raw") ---> Save as */
   void write(FileName name = "", bool reversed=FALSE,
      Volume_Type volume_type=VBYTE) _THROW;

   /** Write image to disk using a file pointer.
       This is the core routine of the previous one. */
   void write(FILE *fh, bool reversed, Volume_Type volume_type);
   //@}
};

/**@name Speed up macros*/
//@{
/** 3D Matrix access.
    This macro does the same as the normal 3D matrix access but
    in a faster way as no function call is generated.
    \\ Ex: VOLMATRIX(V).resize(128,128,128);
    \\ Ex: VOLMATRIX(V2)=VOLMATRIX(V1)+VOLMATRIX(V2); */
#define VOLMATRIX(V) ((V).img)

/** Array access.
    This macro allows you to access to the tridimensional array behind
    the image (float ***). */
#define VOLARRAY(V) (VOL_ARRAY(V).img)

/** Voxel access.
    This macro does the same as the normal voxel access (remember,
    logical access) but in a faster way as no function call is
    generated.
    \\Ex: cout << "Grey level of voxel (2,-3,-3) of the Volume = "
               << VOLVOXEL(V,2,-3,-3) << endl;
    \\Ex: VOLVOXEL(I,2,-3,-3)=VOLVOXEL(I,2,-3,-2); */
#define VOLVOXEL(V,k,i,j) VOL_ELEM(((V).img),(k),(i),(j))

/** Physical voxel access.
    The physical voxel access gives you access to a voxel by
    its physical position and not by its logical one. This access
    shouldn't be used as a custom, use instead the logical access,
    but there might be cases in which this access might be interesting.
    Physical positions start at index 0 in C.
    \\Ex: cout << "This is the first voxel stored in the Volume " <<
          DIRECT_VOLVOXEL(V,0,0,0) << endl;
    @see DIRECT_VOL_ELEM */
#define DIRECT_VOLVOXEL(V,k,i,j) DIRECT_VOL_ELEM(((V).img),(k),(i),(j))
//@}

/* ************************************************************************* */
/* XMIPP VOLUMES                                                                    */
/* ************************************************************************* */
/** Xmipp 3D Volumes.
    The Xmipp volume is a normal volume (inherited from volume class) plus a
    Spider header. This is the appropiate class to work with volumes in
    memory which we want to save later.
    The data in the header is not directly accesible from
    the programmer and must be set through object functions, in this way
    the coherence in the header is assured. See File Formats for
    more information about the Spider format.
    
    In principle, the volume starts at voxel (0,0) but this can be modified
    for any other logical access. See class \Ref{Volume} for more information.

    The Euler angles are useless in a Xmipp volume, and although the
    Spider header has got space for them they are not used and cannot be
    accessed.
    @see Volume
*/

template <class T> class VolumeXmippT: public VolumeT<T> {
protected:
   headerXmipp header;				       // Declares a header

public:
   // Constructors .........................................................
   /**@name Constructors*/
   //@{
   /** Empty constructor.
       Creates an empty (0x0x0) image with no information in the header.
       \\ Ex: VolumeXmipp VX;*/
   VolumeXmippT():VolumeT<T>(){
     if(typeid(T) == typeid(double) )
          header.headerType() = headerXmipp::VOL_XMIPP;
                                        // Sets header of type Image_XMipp
     else if(typeid(T) == typeid(int) )
          header.headerType() = headerXmipp::VOL_INT;
                                        // Sets header of type Image_XMipp
     else if(typeid(T) == typeid(complex<double>) )
          header.headerType() = headerXmipp::VOL_FOURIER;
                                        // Sets header of type Image_XMipp (complex)
     else{
     cout << "\nError: VolumeXmipp should be,complex<double>, double or integer\n";
     exit(0);
     }                                   
   }
                                           
   /** Constructor with size.
       Creates a 0.0 filled volume of size Zdim x Ydim x Xdim.
       \\ Ex: VolumeXmipp<double> VX(64,64,64); */
   VolumeXmippT (int Zdim, int Ydim, int Xdim):VolumeT<T>(Zdim, Ydim, Xdim) {
     if(  typeid(T) == typeid(double) )
          header.headerType() = headerXmipp::VOL_XMIPP;
                                        // Sets header of type Image_XMipp
     else if (typeid(T) == typeid(int) )   
          header.headerType() = headerXmipp::VOL_INT;
                                        // Sets header of type Image int
     else if (typeid(T) == typeid(complex<double>) )  
          header.headerType() = headerXmipp::VOL_FOURIER;
                                        // Sets header of type Image_XMipp (complex)
     else{
     cout << "\nError: VolumeXmipp should be double or integer\n";
     exit(0);
     }                                   
     
     header.set_dimension(Ydim, Xdim);                 // Sets header dimensions     
     header.Slices() = Zdim;                 	       // Sets header Slices 
     header.set_header(); 			       // Initialize header
     header.set_time();				       // Set time and date	
     header.set_date();
   };

   /** Constructor with filename, read from disk.
       The filename given must exist, then the file is loaded in the
       VolumeXmipp class structure. You have loaded the volume at the
       declaration time.
       \\ Ex: VolumeXmipp<double> VX("art0001.vol"); */
   VolumeXmippT(FileName _name):VolumeT<T>(_name){
     if(  typeid(T) == typeid(double) )
          header.headerType() = headerXmipp::VOL_XMIPP;
                                        // Sets header of type Image_XMipp
     else if (typeid(T) == typeid(int) )   
          header.headerType() = headerXmipp::VOL_INT;
                                        // Sets header of type Image int
     else if (typeid(T) == typeid(complex<double>) )   
          header.headerType() = headerXmipp::VOL_FOURIER;
                                        // Sets header of type Image_XMipp (complex)
     else{
       cout << "\nError: VolumeXmipp should be double or integer\n";
       exit(0);
     }                                   
     read(_name);  			  	       // Read image from file
   }
   
   /** Copy constructor.
       \\ Ex: VolumeXmipp<double> VX2(VX1); */
   VolumeXmippT (VolumeXmippT &I): VolumeT<T>(I) {header = I.header;}

   /** Empty image.
       All information is cleared.
       \\Ex: VX.clear();*/
   void clear() {clear_header(); VolumeT<T>::clear();}
   //@}
   
   // Overload some operators ..............................................
   /**@name Some operators*/
   //@{
   /** Show the header information of a Xmipp volume.
       \\ Ex: cout << VX; */
   friend ostream& operator << (ostream& out, const VolumeXmippT<T> &V)
      {out << (VolumeT<T> &) V << V.header; return out;}

   /** Assignment from another Xmipp volume.
       \\ Ex: VolumeXmipp VX1, VX2; VX2=VX1; */
   VolumeXmippT & operator= (const VolumeXmippT<T> &op1)
      {if (&op1!=this) {
          this->VolumeT<T>::operator = (op1);
          header = op1.header;
       }
       return *this;}

   /** Assignment from a generic image.
       \\Ex: Volume<double> V; VolumeXmipp<double> VX; VX=V;*/
   VolumeXmippT & operator= (const VolumeT<T> &op1) {
      if (this!=&op1) {
         this->VolumeT<T>::operator = (op1);
         clear_header(); adjust_header();
      }
      return *this;
   }

   /** Assignment from a 3D matrix.
       \\Ex: matrix3D<float> m; VolumeXmipp VX; VX=m; */
   template <class Type>
      VolumeXmippT& operator=(const matrix3D<Type> &op1) {
         if (&img!=(matrix3D<T> *) &op1) {
            this->VolumeT<T>::operator = (op1);
            clear_header(); adjust_header();
         }
         return *this;
      }

   /** Assignment from any kind of volume.*/
   //Ex:Volume<double> VO; VolumeXmipp VX; VX.assign_from(&VO); 
   template <class Type>
      void assign_from(VolumeT<Type> *v) {*this=*v;}
   //@}

   // Input/Output .........................................................
   /**@name Input/Output*/
   //@{
   /** Read Xmipp volume from disk.
       If the volume doesn't exist at the given path then an exception is
       thrown.        
       The type check is a test performed on input image to check
       if it comes from a big or little endian machine. It is done
       over the \Ref{xmippHeader} field fIform. Sometimes, this value
       is corrupted although the whole image is still valid. You can
       skip this check and provide the reversed status via force_reversed.
       \\ Ex: VX.read("art0001.vol");*/
   void read(const FileName &_name, bool skip_type_check=FALSE,
      bool force_reversed=FALSE) _THROW;

   /** Write Xmipp volume to disk.
       If there is any problem in the writing, an exception is thrown.
       You can give a name to the written volume different from the one used
       when it was read. From this point the filename of the volume has
       changed. This is somehow like the "Save as ..." and "Save".
       \\ Ex: VX.write() ---> Save
       \\ Ex: VX.write("art0002.vol") ---> Save as 
       If force_reversed is TRUE then image is saved in reversed mode,
       if not it is saved in the same mode as it was loaded.*/
   void write(const FileName &_name = "", bool force_reversed=FALSE) _THROW;
   //@}

   // Header operations interface ..........................................
   // These are the only access to the header allowed 
   /**@name Header access*/
   //@{
   /** Adjust header.
       Force header to have the dimensions of the image, time, date updated */
   void adjust_header() {
      if(typeid(T) == typeid(double) )
           header.headerType() = headerXmipp::VOL_XMIPP;
                                         // Sets header of type Image_XMipp
      else if (typeid(T) == typeid(int) )   
           header.headerType() = headerXmipp::VOL_INT;
                                         // Sets header of type Image int
      else if (typeid(T) == typeid(complex<double>) )   
           header.headerType() = headerXmipp::VOL_FOURIER;

      header.set_dimension(YSIZE(img), XSIZE(img)); // Sets header dimensions     
      header.Slices() = ZSIZE(img);                 // Sets header Slices 
      header.set_time();			    // Set time and date 
      header.set_date();
      header.set_title(fn_img);                     // Set title
      header.set_header(); 			    // Initialize header
   }

   /** Resets header. */
   void clear_header() {header.clear();}

   /** Change Filename.
       \\Ex: IX.rename("newName.spd"); */
   void rename(FileName newName)
      {VolumeT<T>::rename(newName); header.set_title(newName);}

   /** Reversed status.
       This is used for the little/big endian process. */
   bool reversed() const {return header.reversed();}
   //@}
};


/* ************************************************************************* */
/* TYPE DEFINITIONS                                                          */
/* ************************************************************************* */
typedef VolumeT<double> Volume;
typedef VolumeXmippT<double> VolumeXmipp;
typedef VolumeT<complex<double> > FourierVolume;
typedef VolumeXmippT<complex<double> > FourierVolumeXmipp;


/**@name Related functions */
//@{
/** True if the given volume is an Xmipp volume. See \Ref{volumeXmipp::read}
    for an explanation of skip_type_check and force_reversed.*/
int Is_VolumeXmipp(const FileName &fn, bool skip_type_check=FALSE,
   bool force_reversed=FALSE) _THROW;

/** True if the given volume is a Fourier Xmipp volume. See \Ref{volumeXmipp::read}
    for an explanation of skip_type_check and force_reversed.*/
int Is_FourierVolumeXmipp(const FileName &fn, bool skip_type_check=FALSE,
   bool force_reversed=FALSE) _THROW;

/** Get size of a volume.
    It returns -1 if the file is not an Xmipp volume.*/
void GetXmippVolumeSize(const FileName &fn, int &Zdim, int &Ydim, int &Xdim);
//@}

//@}
#endif
