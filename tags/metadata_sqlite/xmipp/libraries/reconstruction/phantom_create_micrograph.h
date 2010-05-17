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
#ifndef _PROG_PHANTOM_CREATE_MICROGRAPH_HH
#define _PROG_PHANTOM_CREATE_MICROGRAPH_HH

#include <data/funcs.h>
#include <data/volume.h>
#include "phantom_simulate_microscope.h"

/**@defgroup MicrographProgram phantom_create_micrograph (Microscope simulation)
   @ingroup ReconsLibraryPrograms */
//@{
/* Create micrograph Program Parameters ------------------------------------ */
/** Parameter class for the create micrograph program */
class Prog_Phantom_Create_Micrograph_Parameters
{
public:
    /// Volume
    FileName fn_vol;
    
    /// Rootname for the micrographs
    FileName fn_root;

    /// Size of the micrograph
    int Xdim;
    
    /// Number of micrographs
    int Nmicrographs;

    /// Density of particles (a percentage<100)
    double density;
    
    /// CTF
    Prog_Microscope_Parameters microscope;
public:
    // Volume
    VolumeXmipp V;
    
    // Number of projections per micrograph
    int Nproj;
public:
    /** Empty constructor */
    Prog_Phantom_Create_Micrograph_Parameters();

    /** Read from a command line.
        An exception might be thrown by any of the internal conversions,
        this would mean that there is an error in the command line and you
        might show a usage message. */
    void read(int argc, char **argv);

    /** Usage message.
        This function shows the way of introducing these parameters. */
    void usage();

    /** Show parameters. */
    void show();

    /** Produce side information. */
    void produce_side_info();

    /** Run */
    void run();
};
//@}
#endif
