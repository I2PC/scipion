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
#ifndef _PROG_RANGE_ADJUST_HH
#  define _PROG_RANGE_ADJUST_HH

#include <XmippData/xmippImages.hh>
#include <XmippData/xmippVolumes.hh>

/**@name Adjust grey level range */
//@{
/* Range_adjust Program Parameters ----------------------------------------- */
/** Parameter class for the project program */
class Prog_Range_adjust_Parameters {
public:
   /// min_noise in %
   double min_sigma;
   /// max_noise in %
   double max_sigma;
public:
   /** Read from a command line.
       An exception might be thrown by any of the internal conversions,
       this would mean that there is an error in the command line and you
       might show a usage message. */
   void read(int argc, char **argv) _THROW;

   /** Usage message.
       This function shows the way of introducing this parameters. */
   void usage();

   /** Show parameters. */
   void show();

   /** Range adjust of an image. The input image is modified. */
   void apply(Image *I);

   /** Range adjust of a volume. The input volume is modified. */
   void apply(Volume *V);
};
//@}
#endif
