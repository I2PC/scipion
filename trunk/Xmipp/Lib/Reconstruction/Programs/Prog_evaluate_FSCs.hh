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
#ifndef _PROG_EVALUATE_FSCs_HH
#  define _PROG_EVALUATE_FSCs_HH

#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippVolumes.hh>
#include <XmippData/xmippSelFiles.hh>
#include "../phantom.hh"

/**@name Evaluation FSCs program */
//@{
/* Evaluation FSCs Program Parameters --------------------------------------- */
/** Parameter class for the evaluate FSCs program */
class Prog_Evaluate_FSCs_Parameters {
public:
   /// Phantom filename: either Xmipp volume
   FileName fn_phantom;

   /// Reconstruction filename: Xmipp volume or Selfile
   FileName fn_recons;

   /// Reconstruction filename 2: Selfile
   FileName fn_recons2;

   #define ESTIMATE_SINGLE_FSC         0
   #define ESTIMATE_AVERAGE_RESOLUTION 1
   #define ESTIMATE_AVERAGE_FSC        2
   #define COMPARE_TWO_SETS            3
   /** Action to perform.
       Valid actions: ESTIMATE_SINGLE_FSC, ESTIMATE_AVERAGE_RESOLUTION,
       ESTIMATE_AVERAGE_FSC, COMPARE_TWO_SETS. */
   int action;
   
   /// Sampling rate
   double sampling_rate;

// Side information .........................................................
public:
   // Phantom volume
   VolumeXmipp phantom;
   
   // Reconstructed volume
   VolumeXmipp reconstruction;
   
   // Selfile with all reconstructions
   SelFile SF_recons;
   
   // Selfile with the second set of reconstructions
   SelFile SF_recons2;
   
   // Output file (not always necessary)
   FileName fn_out;

public:
   /** Read from a command line.
       An exception might be thrown by any of the internal conversions,
       this would mean that there is an error in the command line and you
       might show a usage message. */
   void read(int argc, char **argv);

   /** Usage message.
       This function shows the way of introducing this parameters. */
   void usage();

   /** cout << prm; */
   friend ostream & operator << (ostream &out,
      const Prog_Evaluate_FSCs_Parameters &prm);

   /** Produce side information.
       An exception might be thrown if some data cannot be accessed.*/
   void produce_side_info() _THROW;

   /** Compute average and stddev of the resolution of a selfile.
       Given a selfile, what is the average resolution (FSC<0.5) in (1/A).
       The vector of all resolutions is provided to check if there is
       any outlier. */
   void compute_average_resolution(double &avg_resol, double &stddev_resol,
      matrix1D<double> &resol);

   /** Compute the average, min and max FSC of a selfile.*/
   void compute_average_FSC(matrix1D<double> &frequency,
      matrix1D<double> &avg_FSC, matrix1D<double> &min_FSC,
      matrix1D<double> &max_FSC);

   /** Compare the FSCs of two selfiles.
       The average of the difference (F2-F1) is returned as well as
       its stddev at each frequency.
   */
   void compare_two_sets(matrix1D<double> &frequency,
      matrix1D<double> &avg_diff_FSC, matrix1D<double> &stddev_diff_FSC);
};

/** Main Evaluation FSCs routine.
    This is the main function of the program Evaluate. It takes a
    set of evaluation parameters and returns the results in a structure. */
void ROUT_Evaluate_FSCs(Prog_Evaluate_FSCs_Parameters &prm) _THROW;

//@}
#endif
