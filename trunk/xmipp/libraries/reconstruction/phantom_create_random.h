/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
//Sun Nov 14 22:07:48 EST 1999: added binary option (R. Marabini)

#ifndef _PROG_RANDOM_PHANTOM_HH
#  define _PROG_RANDOM_PHANTOM_HH

#include <data/funcs.h>

#include <data/phantom.h>

/**@defgroup RandomPhantom phantom_create_random (Generate phantoms randomly)
   @ingroup ReconsLibraryPrograms */
//@{

/* Random Phantom parameters ----------------------------------------------- */
/** Random Phantom program parameters.
    This class contains all information needed about the Random phantom
    genration  process.
    See the user help on Random Phantom Creation 
    for more information. */
class Prog_Random_Phantom_Parameters
{
public:
    /// Input Phantom description file
    FileName fn_random;
    /// Output Phantom description file
    FileName fn_output;
    /// Minimum volume for the randomly generated features
    double    min_vol;
    /// Use only discrete (integer) values for the features's density
    int  discrete;
    /// center of the features should be at a distance greater than distance
    double RPP_distance;
    /// center of the features should be inside an sphere of radius radius
    double RPP_radius;
    /// Compute volume statistics, and only this
    int N_stats;
    /// Filename with the CTF
    FileName fn_CTF;
    /// Xdim of projections
    int Xdim;
    /// Ydim of projections
    int Ydim;
    /// Target SNR
    double target_SNR;
public:
    /// Empty constructor
    Prog_Random_Phantom_Parameters();

    /** Read parameters from a command line.
        This function reads the parameters from a command line
        defined by argc and argv. An exception might be thrown by any
        of the internal conversions, this would mean that there is
        an error in the command line and you might show a usage message. */
    void read(int argc, char **argv);

    /** Usage message.
        This function shows the way of introdustd::cing this parameters. */
    void usage();
};

/* Random Phantom Side information ----------------------------------------- */
/** Random Phantom Side Information.
    The Side information is useful information that needn't be computed
    more than once at the beginning from the Random Phantom parameters and
    that is used all over the program. */
class Random_Phantom_Side_Info
{
public:
    /// True if in voxel mode
    bool voxel_mode;
    /// Random specification phantom
    Phantom         Random;
    /// Voxel phantom
    VolumeXmipp     VoxelPhantom;
public:
    /** Produce Initial and Side information for Random Phantom.
        This function loads the phantom specification defined in the
        Random Phantom parameters. It also check that the phantom file meets
        the conditions to be a Random Phantom specification (see
        the user help on random phantom creation
        for more information). */
    void produce_Side_Info(const Prog_Random_Phantom_Parameters &prm);
};

/* Generate a realization of the random phantom  --------------------------- */
/** Generate a realization of the random phantom.
    This function returns a realization of the Random phantom in the
    Side Information. It takes care that created features are at least
    of a volume superior to the minimum volume of the parameters. */
void generate_realization_of_random_phantom(
    const Prog_Random_Phantom_Parameters &prm,
    Random_Phantom_Side_Info &side, Phantom &Realization);

/* Main Random Phantom Routine --------------------------------------------- */
/** Main Random Phantom Routine.
    This function takes an input Random Phantom specification and writes
    to disk the generated phantom realization (if any output filename has
    been given). The realization is also
    returned in the Phantom variable. */
void ROUT_random_phantom(const Prog_Random_Phantom_Parameters &prm,
                         Phantom &Realization);
//@}
#endif
