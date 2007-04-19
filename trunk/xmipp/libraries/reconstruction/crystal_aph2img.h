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

#include <data/funcs.h>
#include <data/matrix2d.h>
#include <data/projection.h>
#include <data/selfile.h>
#include <data/matrix1d.h>
#include <interface/aph_origmerg.h>

/**@name Spots <-->Real Space 2D program */
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
   /** Sampling Rate in the different micrographies */
   double SamplingRate;
   /** Generate symmetrical reflections.
       By default, no */
   bool            generate_symmetrical_reflections;
   /** Symmetry group */
   string          str_symmetry_group;
#ifdef NEVERDEFINED
   /** vector perpendicular to projection plane */
   matrix1D<double> v_perpendicular_proj_plane;
#endif
public:
   /* Side information */
   /** This image APH */
   APHFileorigmerg       aph_file;
   /** Matrix that stores the number of spots for a given (k,h) */
   matrix2D<int> Counter;
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
   /** taxa angle,  CONVENTION FOR MEASURING TILT AXIS TO ASTAR IS
    THAT THE ANGLE IS FROM TILTAXIS TO ASTAR IN THE DIRECTION GIVEN
    BY ASTAR TO BSTAR BEING POSITIVE.*/
   double          taxa;
   /** angle between a and b (real space and degrees) */
   double          a_b_ang;//2*180,10*90,5*60
   /** a module(real space in A) */
   double          a_mag;
   /** b module(real space in A) */
   double          b_mag;
   /** Tilt angle following MRC conventions*/
   double          mrc_tilt;
   /* Symmetry group code */
   int             symmetry_group;
   /** MRC micrograph label. If label= -1 -> wildcard */
   int             mrc_label;
public:
   /** This routine reads the parameters, supplied by the user, from a file.
   */
   void read_from_file(const FileName &fnprm);
   /** Show parameters. */
   friend ostream& operator << (ostream &o, const Spot2RealSpace2D_Parameters &prm);
   /** Produce Side Information */
   void produce_SideInfo();
#ifdef NEVERDEFINED
   /** Constructor */
   Spot2RealSpace2D_Parameters()
   {
   v_perpendicular_proj_plane.resize(3);
   }
#endif
};

   void ROUT_Spots2RealSpace(Spot2RealSpace2D_Parameters &prm,
   Projection &prj);
//@}

//@{
/* Real Space --> Spot Program Parameters ------------------------------------ */
/** Real Space --> Spot  parameters.
    Here are all the parameters needed to produce from the real space projections
    2D Fourier Spots .
*/
class RealSpace2Spots2D_Parameters {
public:
   /** Input file (Spider real space)
   */
   FileName        fn_in;
   /** Output file (MRC aph)
   */
   FileName        fnaph_out;
public:
   /* Side information */
   /** This image APH */
   APHFileorigmerg       aph_file;
   /** crystal vector_a */
   matrix1D<double> vector_a;
   /** crystal vector_b */
   matrix1D<double> vector_b;
   /** Rotational angle */
   double          rot;
   /** Tilt angle */
   double          tilt;
   /** Psi angle */
   double          psi;
   /** MRC micrograph label. If label= -1 -> wildcard */
   int             mrc_label;
   /** taxa angle,  CONVENTION FOR MEASURING TILT AXIS TO ASTAR IS
    THAT THE ANGLE IS FROM TILTAXIS TO ASTAR IN THE DIRECTION GIVEN
    BY ASTAR TO BSTAR BEING POSITIVE.*/
   double          taxa;
   /** Tilt angle following MRC conventions*/
   double          mrc_tilt;
public:
   /** This routine reads the parameters, supplied by the user, from a file.
   */
   void read_from_file(const FileName &fnprm);
   /** Show parameters. */
   friend ostream& operator << (ostream &o, const RealSpace2Spots2D_Parameters &prm);
   /** Produce Side Information */
   void produce_SideInfo();

   /** Constructor */
   RealSpace2Spots2D_Parameters()
   {
    vector_a.resize(2);
    vector_b.resize(2);
   }

};

   void ROUT_RealSpace2Spots(RealSpace2Spots2D_Parameters &prm,
   Projection &prj);
//@}
/** Discrete inverse, but not fast Fourier transform
*/
void IDFT(const matrix2D< complex<double> > &FT, matrix2D<double> &I,
   int ydim, int xdim);

/** Discrete direct, but not fast Fourier transform
*/
void DFT(const matrix2D<double> &I,  matrix2D< complex<double> > &FT);

#endif
