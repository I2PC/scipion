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
#ifndef _PROG_ADJUST_VOLUME_HH
#  define _PROG_ADJUST_VOLUME_HH

#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippVolumes.hh>
#include <XmippData/xmippSelFiles.hh>

/**@name Adjust volume program */
//@{
/* Adjust volume Program Parameters ------------------------------------------- */
/** Parameter class for the project program */
class Prog_Adjust_Volume_Parameters {
public:
   /// Filename with the input volume
   FileName fn_vol;
   /// Filename with the input projections
   FileName fn_sel;
   /** Filename of the output volume.
       If empty the input one is used. */
   FileName fn_out;
   /// Optimize
   bool optimize;
   /// Probability of being evaluated
   double probb_eval;
public:
   // Input volume
   matrix3D<double> V;
   // SelFile
   SelFile SF;   
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

   /** Produce side information. */
   void produce_side_info();

   /** Mismatching.
       This function returns the overall mismatiching between the
       experimental projections and the theoretical projections of the current
       volume. */
   double mismatching(double a, double b);

   /** Apply.
       This is the function that really does the job */
   void apply(matrix3D<double> &output_volume);

   /** Run.
       Calls apply and save the result. */
   void run();
};
//@}
#endif
