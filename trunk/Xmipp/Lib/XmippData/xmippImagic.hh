/***************************************************************************
 *
 * Authors: Lorenzo Zampighi and Nelson Tang
 *
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
/*
 * xmippImagic.hh header file defines a subclass of Image for reading
 * Imagic-format files.
 */

#ifndef _XMIPP_IMAGIC_HH
#define _XMIPP_IMAGIC_HH

#include <vector>
#include "xmippImages.hh"
                                                                               
#define GCC_VERSION (__GNUC__ * 10000 \
    + __GNUC_MINOR__ * 100 \
    + __GNUC_PATCHLEVEL__)
/* Test for GCC > 3.3.0 */
#if GCC_VERSION >= 30300
    #include <sstream>
#else
    #include <strstream.h>
#endif


/**@name IMAGIC Images */
//@{

/* Extensions for the Imagic header and image files */
static const char *IMAGIC_HEADER_EXT = "hed";
static const char *IMAGIC_IMAGE_EXT  = "img";

/* String that identifies a selline as an Imagic image */
static const char *IMAGIC_TAG = "imagic:";
/* Length of the tag */
static const short unsigned IMAGIC_TAG_LEN = 7;
/* Separator character */
static const char IMAGIC_TAG_SEP = ':';

/** Types of Imagic files.
    Valid types are IMAGIC_REAL, IMAGIC_INTG, IMAGIC_PACK, IMAGIC_COMP.
*/
enum ImageImagicType { IMAGIC_REAL, IMAGIC_INTG, IMAGIC_PACK, IMAGIC_COMP };

/* structure that holds Imagic image information */
struct ImageImagicInfo
{
  unsigned int num_img;
  unsigned int xsize, ysize;
  vector<ImageImagicType> img_types;
};

/** Imagic Image class */
template <class T> class ImageImagicT: public ImageT<T>
{
public:
  /** Empty constructor. */
  ImageImagicT() : ImageT<T>() { name_parsed = false; };
  /** Constructor with size */
  ImageImagicT(int Ydim, int Xdim) :
    ImageT<T>(Ydim, Xdim) { name_parsed = false; };
  /** Constructor with image name */
  ImageImagicT(FileName _name) : ImageT<T>(_name) { name_parsed = false; };

  /** Copy constructor */
  ImageImagicT(const ImageImagicT &I) : ImageT<T>(I)
  {
    if ((name_parsed = I.name_parsed))
    {
      hedfname = I.hedfname;
      imgfname = I.imgfname;
      imgnum = I.imgnum;
    }
  };

  /** Rename. */
  virtual void rename (FileName newName)
  {
    if (newName != fn_img)
    {
      ImageT<T>::rename (newName);
      name_parsed = false;
    }
  };
  
  /** Clear */
  virtual void clear() { ImageT<T>::clear(); name_parsed = false; };

  /** Assignment operator. */
  ImageImagicT<T>& operator= (const ImageImagicT<T> &I)
  {
    if (this != &I)
    {
      this->ImageT<T>::operator= (I);
      if ((name_parsed = I.name_parsed))
      {
	hedfname = I.hedfname;
	imgfname = I.imgfname;
	imgnum = I.imgnum;
      }
    }
    return (*this);
  };

  /** Assignment operator. */
  ImageImagicT<T>& operator= (const ImageT<T> &I)
  {
    if (this != &I)
    {
      this->ImageT<T>::operator= (I);
      name_parsed = false;
    }
    return (*this);
  };

  /** Assignment operator */
  template <class T1>
  ImageImagicT<T>& operator= (const matrix2D<T1> &m)
  {
    if (&img != (matrix2D<T> *) &m)
    {
      this->ImageT<T>::operator= (m);
      name_parsed = false;
    }
    return (*this);
  };

  /** Show */
  friend ostream& operator<< (ostream& out, const ImageImagicT<T> &I)
  {
    out << (ImageT<T>&) I;
    out << "IMAGIC header fname: " << get_hedfname() << endl;
    out << "IMAGIC image fname:  " << get_imgfname() << endl;
    out << "image number:        " << get_imgnum() << endl;
    return (out);
  };

  /** Read file. */
  bool read (const FileName &name) _THROW;

  /** Write file. Not implemented. */
  void write (const FileName &name="", bool reversed=FALSE,
	      Image_Type image_type=IBYTE) _THROW
  {
    REPORT_ERROR (1503, "ImageImagic::write: can't directly save");
  };
  
  /** Low level write. Not implemented. */
  void write (FILE * &fh, bool reversed, Image_Type image_type)
  {
    REPORT_ERROR (1503, "ImageImagic::write: can't directly save");
  };

  /** Get Header name. */
  const FileName& getHedFname() { parseFname(); return (hedfname); };
  /** Get Image name */
  const FileName& getImgFname() { parseFname(); return (imgfname); };
  /** Get image number */
  const int getImgNum() { parseFname(); return (imgnum); };
  
  virtual void parseFname(); // Not meant to be called directly, but needs to
                             // be kept public so the template can be
                             // instantiated.

protected:
  bool name_parsed;
  FileName hedfname, imgfname;
  int imgnum;  
};

/** Alias for ImageImagic type */
typedef ImageImagicT<double> ImageImagic;

// Utility functions

/** Creates a string of the format 'imagic:hedfname:num'.
    suitable for use as an image name. */
inline string ImagicMakeName (const char *hed_fname, unsigned int imgnum)
{
#if GCC_VERSION < 30300
       char aux[15];
       ostrstream ss(aux,sizeof(aux));
#else
       ostringstream ss;
#endif
  ss << IMAGIC_TAG << hed_fname << IMAGIC_TAG_SEP << imgnum;
  return (ss.str());
};

/** Looks at an IMAGIC header file for information about
    the images it contains.*/
const ImageImagicInfo ImagicGetImgInfo (const FileName &hed_fname);

/** Creates a new Imagic header/image file pair,
    filling it with the data pointed to by the vector parameter. */
template <class T>
bool ImagicWriteImagicFile (const FileName &hed_fname,
			    const vector<ImageT<T> *> &imgs,
			    ImageImagicType img_type = IMAGIC_REAL);

//@}
#endif
