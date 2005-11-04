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

#ifndef _XMIPPTOMOGRAM_H
   #define _XMIPPTOMOGRAM_H

/* ************************************************************************* */
/* INCLUDES                                                                  */
/* ************************************************************************* */
#include "xmippFuncs.hh"
#include "xmippMatrices3D.hh"

/* ************************************************************************* */
/* TOMOGRAM                                                               */
/* ************************************************************************* */
/**@name Tomograms*/
//@{

/** Tomogram class.
    This class manages a large tomogram on disk. The volume is not loaded
    into memory, that should avoid memory problems 
*/
class Tomogram {
protected:
   FileName                fn_tomogram;
   FileName                fn_inf;
   long                    Xdim;
   long                    Ydim;
   long                    Zdim;
   long                    XYdim;
   int                     __depth;
   int                     __offset;
   bool                    __reversed;
   bool                    __is_signed;
   unsigned char           *m8;
   short int               *m16;
   unsigned short int      *um16;
   float                   *m32;
   int                     fh_tomogram;
public:
   /** Constructor */
   Tomogram() {clear();}

   /** Clear */
   void clear();

   /** Get tomogram depth. */
   int depth() const {return __depth;}

   /** Is the file reversed in disk */
   bool reversed() {return __reversed;}

   /** Open tomogram.
       An exception is thrown if the file is not valid. */
   void open_tomogram(const FileName &fn_tomogram, 
      bool reversed=false);

   /** Close tomogram.
       After working with the file, you must close it. */
   void close_tomogram();
   
   /** tomogram filename. */
   string tomogram_name() { return( fn_tomogram ); }

   /** Access to array of 8 bits. */
   unsigned char * array8() const {return m8;}

   /** Access to array of 16 bits. */
   short int * array16() const {return m16;}

   /** Access to unsigned array of 16 bits. */
   unsigned short int * arrayU16() const {return um16;}

   /** Access to array of 32 bits. */
   float * array32() const {return m32;}

   /** Pixel access for reading.
       These coordinates follow the physical Xmipp \URL[convention]
       {../../../Extra_Docs/Conventions.html} for coordinates */
   float operator ()(int x, int y, int z) const {
      if (y<0 || y>=Ydim || x<0 || x>=Xdim || z<0 || z>=Zdim)
         REPORT_ERROR(1,"tomogram::(): index out of range");
      if      (__depth== 8) {
         return m8[z*XYdim+y*Xdim+x];
      } else if (__depth==16) {
        if (__is_signed) {
            short int retval=m16[z*XYdim+y*Xdim+x];
            if (__reversed) {
               unsigned char *ptr=(unsigned char *)&retval, temp;
               SWAP(*ptr,*(ptr+1),temp);
            }
            return retval;
         } else {
            unsigned short int retval=um16[z*XYdim+y*Xdim+x];
            if (__reversed) {
	       unsigned char *ptr=(unsigned char *)&retval, temp;
	       SWAP(*ptr,*(ptr+1),temp);
            }
            return retval;
         }
      } else if (__depth==32) {
         float retval=m32[z*XYdim+y*Xdim+x];
         if (__reversed) {
	    unsigned char *ptr=(unsigned char *)&retval, temp;
	    SWAP(*ptr,*(ptr+3),temp);
	    SWAP(*(ptr+1),*(ptr+2),temp);
         }
         return retval;
      } else REPORT_ERROR(1,"tomogram::(): depth is not 8, 16 or 32");
   }

   /** Pixel access for writing. */
   void set_val(int x, int y, int z, double new_val) {
      if (y<0 || y>=Ydim || x<0 || x>=Xdim || z<0 || z>=Zdim)
         REPORT_ERROR(1,"tomogram::set_val: index out of range");
      if      (__depth== 8) m8[z*XYdim+y*Xdim+x]=(unsigned char) new_val;
      else if (__depth==16 && __is_signed) m16[z*XYdim+y*Xdim+x]=
                                                   (short int) new_val;
      else if (__depth==16 && !__is_signed) um16[z*XYdim+y*Xdim+x]=
                                                   (unsigned short int) new_val;
      else if (__depth==32) m32[z*XYdim+y*Xdim+x]=(float) new_val;
      else REPORT_ERROR(1,"tomogram::set_val: depth is not 8, 16 or 32");
   }

   /** Return tomogram size */
   void size(int &_Xdim, int &_Ydim, int &_Zdim) const
      {_Xdim=Xdim; _Ydim=Ydim; _Zdim=Zdim;}
   
   /** Get piece of tomogram.
       The suggested initial point (r0) and length are provided. If the piece
       fits into the tomogram, the suggested initial point and length are
       kept. If the final point (rF=r0+length) falls outside the tomogram,
       then the initial point is shifted until it fits. The new initial point
       is returned. If the suggested length is bigger than the tomogram size,
       then an exception is thrown. */
    void get_piece(matrix1D<int> &r0, matrix1D<int> &length,
       matrix3D<double> &piece);

    /** Set piece in tomogram.
       This routine does the opposite function than the previous one. The
       initial point is assumed to come from it, therefore no check is done
       on its validity. */
    void set_piece(matrix1D<int> &r0, matrix1D<int> &length,
       matrix3D<double> &piece);
};
//@}
#endif
