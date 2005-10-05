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
#ifndef _PROG_SYMMETRIZE_HH
#  define _PROG_SYMMETRIZE_HH

#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippVolumes.hh>
#include <XmippData/xmippMasks.hh>
#include "../symmetries.hh"

/**@name Symmetrize Program */
//@{
/* Test parameters --------------------------------------------------------- */
/// Symmetrize Parameters
class Symmetrize_Parameters {
public:
   /// input file
   FileName        fn_in;
   /// output file
   FileName        fn_out;
   /// symmetry file
   FileName        fn_sym;
   /// Do not generate subgroup
   bool            do_not_generate_subgroup;
   /// wrap or don't wrap input file during symmetrisation
   bool            wrap;
public:
   /** Read parameters from command line. */
   void read(int argc, char **argv);
   
   /** Usage */
   void usage();

   /** Show parameters */
   friend ostream & operator << (ostream &out, const Symmetrize_Parameters
      &prm);
};

/** Really symmetrize.*/
void symmetrize(const SymList &SL, VolumeXmipp &V_in, VolumeXmipp &V_out,
   bool wrap=TRUE, bool show_progress=FALSE);

/** Really symmetrize using Bsplines */
void symmetrize_Bspline(const SymList &SL, VolumeXmipp &V_in, VolumeXmipp &V_out,
			int Splinedegree, bool wrap, bool do_outside_avg);

/** Main program */
void ROUT_symmetrize(const Symmetrize_Parameters &prm);
//@}

#endif
