/***************************************************************************
 *
 * Authors: Carlos Oscar (coss@cnb.uam.es)
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

#ifndef _XMIPPMICROGRAPH_H
   #define _XMIPPMICROGRAPH_H

/* ************************************************************************* */
/* INCLUDES                                                                  */
/* ************************************************************************* */
#include <vector>
#include "xmippFuncs.hh"
#include "xmippImages.hh"

/* ************************************************************************* */
/* FORWARD DEFINITIONS                                                       */
/* ************************************************************************* */
// This forward definitions are needed for defining operators functions that
// use other clases type

/* ************************************************************************* */
/* MICROGRAPHY                                                               */
/* ************************************************************************* */
/**@name Micrographs*/
//@{
/** Particle coordinates.
    This structure stores the X,Y position of the particle. */
struct Particle_coords {
   /// Label
   int label;
   /// X position
   int X;
   /// Y position
   int Y;
   /// Valid
   bool valid;
};

/** Micrography class.
    This class manages a large micrograph on disk. The image is not loaded
    into memory, that should avoid memory problems 
*/
class Micrograph {
protected:
   /* This image will contain a single particle from the micrograph,
      this is done to avoid asking/freeing memory all time. */
   Image                   single_particle;
   vector<Particle_coords> coords;
   FileName                fn_coords;
   FileName                fn_micrograph;
   FileName                fn_inf;
   int                     X_window_size;
   int                     Y_window_size;
   int                     Xdim;
   int                     Ydim;
   int                     __depth;
   bool                    __reversed;
   unsigned char           *m8;
   short int               *m16;
   float                   *m32;
   bool                    __scaling_valid;
   float                   __a;
   float                   __b;
   /* bool                    __in_core; */
   int                     fh_micrograph;
   vector<string>          labels;
public:
   /** Constructor */
   Micrograph() {clear();}

   /** Clear */
   void clear();

   /** Get micrograph depth. */
   int depth() const {return __depth;}

   /** Is the file reversed in disk */
   bool reversed() {return __reversed;}

   /** Open micrograph.
       An exception is thrown if the file is not valid. */
   void open_micrograph(const FileName &fn_micrograph, /*bool in_core=FALSE,*/
      bool reversed=FALSE) _THROW;

   /** Close micrograpgh.
       After working with the file, you must close it. */
   void close_micrograph();
   
   /** Compute scaling for 8 bits */
   void compute_8_bit_scaling();

   /** Micrograph filename. */
   string micrograph_name() { return( fn_micrograph ); }

   /** Save coordinates to disk. */
   void write_coordinates(int label, const FileName &fn_coords="") _THROW;

   /** Read coordinates from disk.
       Coordinates are read into the selected family, the rest of
       families are untouched as well as the coordinates already belonging
       to this family */
   void read_coordinates(int label, const FileName &fn_coords) _THROW;

   /** Particle number.
       Number of particles in the coordinate list */
   int ParticleNo() const {return coords.size();}

   /** Particle.
       Return the list of particles. */
   vector<Particle_coords> & Particles() {return coords;}

   /** Set window size.
       This window is set upon each coordinate and is used to cut all
       images. */
   void set_window_size(int _X_window_size, int _Y_window_size)
      {X_window_size=_X_window_size; Y_window_size=_Y_window_size;}

   /** Scissor.
       The single particle is selected by an index within the particle
       coordinate list. If the index is beyond the number of particles
       \Ref{ParticleNo}, or the window size is not set (\Ref{set_window_size})
       an exception is thrown.
       
       Make sure that index n represents a valid particle before cutting it
       
       The scale affects the particle position, such that the position cut
       is pos*scale, but not the window size.
       
       If only check is true then the particle is not scissored, but
       the routine only checks if it can be done.
       
       Returns 0 if an error ocurred and 1 if everything is all right*/
   int scissor(const Particle_coords &P, Image &result,
      double scaleX=1, double scaleY=1, bool only_check=false) _THROW;

   /** Access to array of 8 bits. */
   unsigned char * array8() const {return m8;}

   /** Access to array of 16 bits. */
   short int * array16() const {return m16;}

   /** Access to array of 32 bits. */
   float * array32() const {return m32;}

   /** Pixel access for reading.
       These coordinates follow the physical Xmipp \URL[convention]
       {../../../Extra_Docs/Conventions.html} for coordinates */
   float operator ()(int x, int y) const _THROW {
      if (y<0 || y>=Ydim || x<0 || x>=Xdim)
         REPORT_ERROR(1,"Micrograph::(): index out of range");
      if      (__depth== 8) {
         return m8[y*Xdim+x];
      } else if (__depth==16) {
         unsigned short int retval=m16[y*Xdim+x];
         if (__reversed) {
	    unsigned char *ptr=(unsigned char *)&retval, temp;
	    SWAP(*ptr,*(ptr+1),temp);
         }
         return retval;
      } else if (__depth==32) {
         float retval=m32[y*Xdim+x];
         if (__reversed) {
	    unsigned char *ptr=(unsigned char *)&retval, temp;
	    SWAP(*ptr,*(ptr+3),temp);
	    SWAP(*(ptr+1),*(ptr+2),temp);
         }
         return retval;
      } else REPORT_ERROR(1,"Micrograph::(): depth is not 8, 16 or 32");
   }

   /** Pixel access for writing. */
   void set_val(int x, int y, double new_val) _THROW {
      if (y<0 || y>=Ydim || x<0 || x>=Xdim)
         REPORT_ERROR(1,"Micrograph::set_val: index out of range");
      if      (__depth== 8) m8[y*Xdim+x]=(unsigned char) new_val;
      else if (__depth==16) m16[y*Xdim+x]=(short int) new_val;
      else if (__depth==32) m32[y*Xdim+x]=(float) new_val;
      else REPORT_ERROR(1,"Micrograph::set_val: depth is not 8, 16 or 32");
   }

   /** Pixel value with 8 bits. */
   unsigned char val8(int x, int y) const
       {if (!__scaling_valid) return (unsigned char) (*this)(x,y);
        else return (unsigned char) (__a*(*this)(x,y)+__b);}

   /** Produce all single particle images.
       The file fn_micrograph+".sel" is also generated. The angle is the angle
       from the Y axis to the tilt axis, angles are positive clockwise.
       Images are rotated by -ang.
       If this angle is 0 no rotation is applied.*/
   void produce_all_images(int label, const FileName &fn_root,
      int starting_index=1, const FileName &fn_image="", double ang=0,
      double gamma=0., double psi=0.) _THROW;

   /** Search coordinate near a position.
       By default the precission is set to 3 pixels. The index of the coordinate
       within the list is returned. Returns -1 if none. */
   int search_coord_near(int x, int y, int prec=3) const;

   /** Remove a coordinate from the coordinate list.
       An exception is thrown if the index is out of range within the
       coordinate list */
   void invalidate_coord(int n) _THROW;

   /** Add coordinate. */
   void add_coord(int x, int y, int label);
   
   /** Move last coordinate to this position. */
   void move_last_coord_to(int x, int y);

   /** Access to coordinate structure.
       If the index is out of range then an exception is thrown. */
   Particle_coords & coord(int n) _THROW {
      if (n<0 || n>ParticleNo())
         REPORT_ERROR(1,"Micrograph::coord(): index out of range");
      return coords[n];
   }

   /** Add label.
       The index assigned to the label is returned */
   int add_label(const string &label)
      {labels.push_back(label); return labels.size()-1;}

   /** Number of labels. */
   int LabelNo() {return labels.size();}
   
   /** Get a label. 
       An exception is thrown if the index is greater than the
       number of labels */
   string & get_label(int n) _THROW {
      if (n<0 || n>LabelNo())
         REPORT_ERROR(1,"Micrograph::get_label(): index out of range");
      return labels[n];
   }
   
   /** Return micrograph size */
   void size(int &_Xdim, int &_Ydim) const {_Xdim=Xdim; _Ydim=Ydim;}
};

/** Downsample.
    The mask must be normalized to have energy=1. Two runs are done,
    in the first one the input and output ranges are studied. In the
    second one, the input image is effectively downsampled and rescaled
    such that input and output images have the same range. Output
    grey values always start at 0. If the input and output images have
    different bit size, then the range is scaled by the bit difference, ie,
    if the input ranges 0-255, the output will range between 0 and 65535 */
void downsample(const Micrograph &M, int Xstep, int Ystep,
   const matrix2D<double> &kernel, Micrograph &Mp);

/**@name Normalization
   This functions implement the normalization of a single image. They should
   be called with all images in the corresponding SelFile. In the
   following documentation m(x) is the mean of x, v(x) is the variance,
   bg(x) is its background.
   
   The original image is X and is supposed to be related to each single
   projection by I=a*(X+n)+b where a and b are different for every projection.
   
   Noise is assumed to follow a gaussian distribution N(0,sqrt(v(n)))
   
   In general the background mask is used only to compute statistics, while
   the mask is the one which is really applied to the image. Supply NULL if
   you don't want any mask to be applied
   
   When b is used
   it is measured as the mean in the background, while a*sqrt(v(n)) is the
   standard deviation in the same area.
   */
//@{
   /** OldXmipp normalization.
      Formula:
      \begin{verbatim}
      I'=(I-m(I))/sqrt(v(I))
      \end{verbatim}
      Properties:
      \begin{verbatim}
      m(I')=0                             m(bg(I'))=-m(x)/sqrt((v(X)+v(n)))
      v(I')=1                             v(bg(I'))=v(n)/(v(X)+v(n))
      \end{verbatim}
      Comments: it's not bad but positivity constraints cannot be imposed
   */
   void normalize_OldXmipp(Image *I);

   /** Near_OldXmipp normalization.
      Formula:
      \begin{verbatim}
      I'=(I-m(I))/a*v(n)
      \end{verbatim}
      Properties:
      \begin{verbatim}
      m(I')=0                             m(bg(I'))=-m(x)/sqrt(v(n))
      v(I')=(v(X)+v(n))/v(n)              v(bg(I'))=1
      \end{verbatim}
      Comments: it's not bad but positivity constraints cannot be imposed
   */
   void normalize_Near_OldXmipp(Image *I, const matrix2D<int> &bg_mask);

   /** OldXmipp decomposition.
      Formula:
      \begin{verbatim}
      I'=(I-b)/a*sqrt(v(n))
      I''=I'*mask
      I'''=(I''-m(I''))/sqrt(v(I''))
      \end{verbatim}
      Properties:
      \begin{verbatim}
      m(I')=m(X)/sqrt(v(n))               m(bg(I'))=0
      v(I')=(v(X)+v(n))/v(n)              v(bg(I'))=1
      \end{verbatim}
      Comments: it's not bad but positivity constraints cannot be imposed.
         If no mask is applied, then this formula is a beautiful decomposition
         of the OldXmipp method in two steps.
   */
   void normalize_OldXmipp_decomposition(Image *I,
      const matrix2D<int> &bg_mask, const matrix2D<double> *mask=NULL);

   /** Michael's normalization.
      Formula:
      \begin{verbatim}
      I'=(I-b)/b
      \end{verbatim}
      Properties:
      \begin{verbatim}
      m(I')=0                             m(bg(I'))=-a*m(x)/b
      v(I')=a^2*(v(X)+v(n))/b^2           v(bg(I'))=a^2*v(n)/b^2
      \end{verbatim}
      Comments: it's not bad but positivity constraints cannot be imposed and
         the statistical properties are not so good.
   */
   void normalize_Michael(Image *I, const matrix2D<int> &bg_mask);

   /** NewXmipp's normalization.
      Formula:
      \begin{verbatim}
      I'=(I-b)/a*sqrt(v(n))
      // I''=(I'>0)? I':0
      // I'''=I''-a*sqrt(v(n)/2*PI)
      I''''=I'''*mask
      \end{verbatim}
      Properties:
      \begin{verbatim}
      m(I')=m(X)/sqrt(v(n))               m(bg(I'))=0
      v(I')=(v(X)+v(n))/v(n)              v(bg(I'))=1
      \end{verbatim}
      Comments: In general, we cannot assure that mass projects into positive
         numbers, so the "denoising" capability directly on the images is
	 disabled. However, a positivity constraint can be applied on the 3D
	 volume.
   */
   void normalize_NewXmipp(Image *I, const matrix2D<int> &bg_mask);

   /** NewXmipp 2's normalization.
      Formula:
      \begin{verbatim}
      I'=(I-m(bg))/(m(I)-m(bg))
      \end{verbatim}
      Properties:
      \begin{verbatim}
      m(I')=m(X)/sqrt(v(n))               m(bg(I'))=0
      v(I')=(v(X)+v(n))/v(n)              v(bg(I'))=1
      \end{verbatim}
      Comments: In general, we cannot assure that mass projects into positive
         numbers, so the "denoising" capability directly on the images is
	 disabled. However, a positivity constraint can be applied on the 3D
	 volume.
   */
   void normalize_NewXmipp2(Image *I, const matrix2D<int> &bg_mask);
//@}
//@}
#endif
