/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2002)
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
#ifndef _PROG_ANGULAR_DISTANCE
   #define _PROG_ANGULAR_DISTANCE

#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippDocFiles.hh>
#include <Reconstruction/symmetries.hh>

/**@name Angular Distance */
//@{
/** Angular Distance parameters. */
class Prog_angular_distance_prm {
public:
   /** Filename angle doc 1 */
   FileName fn_ang1;
   /** Filename angle doc 2 */
   FileName fn_ang2;
   /** Filename symmetry file */
   FileName fn_sym;
   /** Filename of output file with merging */
   FileName fn_ang_out;
   /** Check mirrors for Spider APMQ */
   bool check_mirrors;
public:
   // DocFile 1
   DocFile DF1;
   // DocFile 2
   DocFile DF2;
   // Symmetry List
   SymList SL;
public:
   /// Read argument from command line
   void read(int argc, char **argv) _THROW;

   /// Show
   void show();

   /// Usage
   void usage();

   /** Produce side info.
       Read all document files and symmetry list if any.
       
       An exception is thrown if both files are not of the same length*/
   void produce_side_info() _THROW;
   
   /** Second angle set.
       Given two sets of angles, this function modifies set 2 so that
       the difference with set 1 are minimized searching in the second
       way of expressing the angles.
       
       The sum of the absolute values of the differences among angles.
       If the projdir_mode is false then the angular distance is measured
       among all axes. If it is true, then it is only measured in the
       projection direction.*/
   double Prog_angular_distance_prm::second_angle_set(
      double rot1, double tilt1, double psi1,
      double &rot2, double &tilt2, double &psi2, bool projdir_mode=false);
   
   /** Check symmetries.
       Given two sets of angles, this function modifies set 2 so that
       the difference with set 1 are minimized searching in the symmetry
       list and the second set. Return the angle distance.
       
       See the method second_angle_set of this class to understand
       projdir_mode*/
   double Prog_angular_distance_prm::check_symmetries(
      double rot1, double tilt1, double psi1,
      double &rot2, double &tilt2, double &psi2, bool projdir_mode=false);
};
//@}
#endif
