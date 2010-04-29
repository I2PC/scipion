/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@mipg.upenn.edu)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Copyright (c) 1999 , CSIC .
 *
 * Permission is granted to copy and distribute this file, for noncommercial
 * use, provided (a) this copyright notice is preserved, (b) no attempt
 * is made to restrict redistribution of this file, and (c) this file is
 * restricted by a compilation copyright.
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 *
 *****************************************************************************/
#ifndef _PROG_SPOTS2REALSPACE3D_HH
#  define _PROG_SPOTS2REALSPACE3D_HH

#include <data/funcs.h>
#include <data/matrix3d.h>
#include <interface/aph3d.h>
#include <data/volume.h>
#include <data/types.h>

/**@defgroup spots2realspace3d spots2realspace3d (Convert crystal 3DSpots to Real Space volume)
   @ingroup ReconsLibraryPrograms */
//@{
/* Spot(3D)-->Real Space Program Parameters ------------------------------------ */
/** Spot(3D) --> Real Space parameters.
    Here are all the parameters needed to produce the real space projections
    from the 3D Fourier Spots.
*/
class Spot2RealSpace3D_Parameters
{
public:
    /** Input file (CCP4 prepmklcf.log)
    */
    FileName        fnaph_in;
    /** output file (Spider real space)
    */
    FileName        fn_out;
    /** No. of unit cells along the X,Y,Z axis.
    */
    Matrix1D<int> NoCells;
    /** Move the phase origin (degrees)
    */
    Matrix1D<double> Phase_Shift;
    /** Keep contrast.
        If FALSE then the image is contrast is reversed.
    */
    int             KeepContrast;
    int             SpaceGroup;
    /** number of samples in Real space  (X,Y,Z)
    */
    Matrix1D<int>             Celldim;
    /** This image APH */
    APHFile3D       aph_file;
public:
    /** This routine reads the parameters, supplied by the user, from a file.
    */
    void read_from_file(const FileName &fnprm);
    /** Show parameters. */
    friend std::ostream& operator << (std::ostream &o, const Spot2RealSpace3D_Parameters &prm);
};
/** Inverse Fourier trnsform */
void IDFT_3D(const Matrix3D< std::complex<double> > &FT, Matrix3D<double> &I);
void ROUT_Spots2RealSpace_3D(Spot2RealSpace3D_Parameters &prm,
                             VolumeXmipp &V1);

/** Computing remaining reflections from those of asymmetric unit (P1). */
void symmetrize_P1(Matrix3D< std::complex<double> > &FT,
                   Spot2RealSpace3D_Parameters &prm);

/** Computing remaining reflections from those of asymmetric unit. The
 asymmetric unit involves the reflections with H>=0 and K >= H. The remaining
 ones are the ones which are computed by the routine.

 It is possible to impose the symmetry in the phase in the case of special
 reflections. In P4212 the special reflections are:

               Real               Imaginary
              ------             -----------
              (2n,0,L)            (2n+1,0,L)
              (0,2n,L)            (0,2n+1,L)
              (H,K,0)
              (H,H,L)

 This is enabled/disabled by means of the argument (impose).

 */

void symmetrize_P4212(Matrix3D< std::complex<double> > &FT,
                      Spot2RealSpace3D_Parameters &prm);

/** Brings the input reflection  H K L into the asymmetric unit according to
 the symmetry P4212. The routine returns the new values of the indexes
 H' K' L' within the asymmetric unit, as well as a flag indicating whether
 the reflection is special. And, in such a case, returns its phase value
 (0 or 90 degrees) corresponding to its character real or imaginary.

 The asymmetric unit in P4212 involves H,K,Z >=0 and H <= K. When H < 0,
 or H==0 and K < 0, or H==K==0 and Z < 0 the Conjugate Symmetry Property
 of the DFT is applied in order to bring the reflection into the H,K >= 0
 and make easier the subsequent processes.

                                                                           */

void symmetrize_P4(Matrix3D< std::complex<double> > &FT,
                   Spot2RealSpace3D_Parameters &prm);


void AsymmUnitP4212(int *ih, int *ik, int *il, int *ip1, int *ip2,
                    int *spec, int *iptest);
/**
  Does matrix multiplication to bring reflections into the asymmetric unit.

      (H' K' Z' AMP' PHS') = (H K Z AMP PHS) <A>

  where <A> has form:

         A[0]  A[2]   0     0    A[5]
         A[1]  A[3]   0     0    A[6]
          0     0    A[4]   0     0
          0     0     0     1     0
          0     0     0     0    A[7]

   for all cases.*/
void AsymmUnitP4(int *ih, int *ik, int *il, int *ip1, int *ip2,
                 int *spec, int *iptest);

void MatrixMult(int A[], int *ih, int *ik, int *il, int *ip1, int *ip2);

/**
   Checks if the current reflection (H K L) is a special reflection, i.e.,
   whose phase must be either real (0 or PI) or imaginary (PI/2 or 3*PI/2).

   It returns spec=1 if the reflection is special. Otherwise, 0.
   It returns in iptest 0 if real and 90 is imaginary.

   Conditions checked are:

      - H=0 special
      - K=0 special
      - Z=0 special
      - H=K special
      - If for H=0 or K=0 K+H is odd, it indicates an imaginary value for the
        reflection.

      Except the last condition, all other special reflections are real.

      Summary of special reflections.

                 Real               Imaginary
                ------             -----------
                (2n,0,L)            (2n+1,0,L)
                (0,2n,L)            (0,2n+1,L)
                (H,K,0)
                (H,H,L)
                                                                              */
void CheckSpecP4212(int ih, int ik, int il, int *spec, int *iptest);
void CheckSpecP4(int ih, int ik, int il, int *spec, int *iptest);
//@}
#endif
