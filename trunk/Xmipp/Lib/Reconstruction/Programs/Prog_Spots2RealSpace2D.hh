/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@mipg.upenn.edu)
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
#ifndef _PROG_SPOTS2REALSPACE2D_HH
#  define _PROG_SPOTS2REALSPACE2D_HH

#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippMatrices2D.hh>
#include <XmippData/xmippProjection.hh>
#include <XmippInterface/xmippAPH.hh>

/**@name Spots-->Real Space 2D program */
//@{
/* Spot-->Real Space Program Parameters ------------------------------------ */
/** Spot --> Real Space parameters.
    Here are all the parameters needed to produce the real space projections
    from the 2D Fourier Spots.
*/
class Spot2RealSpace2D_Parameters {
public:
   /** Input file (MRC aph)
   */
   FileName        fnaph_in;    
   /** Reference Image (tilt=0) */
   FileName        fnaph_ref;
   /** output file (Spider real space)
   */
   FileName        fn_out;
   /** Maximum spot quality factor taken into account in the Fourier transform.
   */
   int             maxIQ;
   /** No. of unit cells along the X axis.
   */
   int             NoCells_X;       
   /** No. of unit cells along the Y axis.
   */
   int             NoCells_Y;
   /** Move the phase origin (degrees)
   */
   matrix1D<double> Phase_Shift; 
// Not longer needed
//   /** Tilt sign (+1, -1) */
//   int             tilt_sign;
   /** Keep contrast.
       If FALSE then the image is contrast is reversed.
   */
   int             KeepContrast;  
   /** number of samples in Real space  (X)
   */
   int             CellXdim; 
   /** number of samples in Real space  (Y)
   */
   int             CellYdim;
   /** Amplitud factor. The magnitude of the Fourier transform will be divided
   by this guy.  It can be calculated from the output of the MRC program 
   origtilt
   */   
   float           Scale_Factor;
   /** Scale between the magnification in the different micrographies */
   double SamplingScale;
   /** Align A axis with x.
       This is used for phantoms, where usually it is */
   bool            align_a_axis_with_x;
   /** Generate symmetrical reflections.
       By default, no */
   bool            generate_symmetrical_reflections;
   /** Symmetry group */
   string          str_symmetry_group;
public:
   /* Side information */
   /** This image APH */
   APHFile2D       aph_file;
   /** Reference APH */
   APHFile2D       aph_ref;
   /** Lattice vectors in the volume */
   matrix2D<double> vol_latt_vec;
   /** Lattice vectors in the volume */
   matrix2D<double> proj_latt_vec;
   /** Rotational angle */
   double          rot;
   /** Tilt angle */
   double          tilt;
   /** Psi angle */
   double          psi;
   /* Symmetry group code */
   int             symmetry_group;
//   /** Mirror correction phase shift in X */
//   double          mirror_phase_X;
//   /** Mirror correction phase shift in Y */
//   double          mirror_phase_Y;
public:
   /** This routine reads the parameters, supplied by the user, from a file. 
   */
   void read_from_file(const FileName &fnprm);
   /** Show parameters. */
   friend ostream& operator << (ostream &o, const Spot2RealSpace2D_Parameters &prm);
   /** Produce Side Information */
   void produce_SideInfo() _THROW;
};

   void ROUT_Spots2RealSpace(Spot2RealSpace2D_Parameters &prm,
   Projection &prj);
   
//@}
#endif
