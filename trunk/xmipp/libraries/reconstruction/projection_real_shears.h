/***************************************************************************
 *
 * Authors:     Slavica JONIC (slavica.jonic@impmc.jussieu.fr, slavica.jonic@a3.epfl.ch)
 *              Jean-Noël PIOCHE (jnp95@hotmail.com)
 *
 * Biomedical Imaging Group, EPFL (Lausanne, Suisse).
 * Structures des Assemblages Macromoléculaires, IMPMC UMR 7590 (Paris, France).
 * IUT de Reims-Châlons-Charleville (Reims, France).
 *
 * Last modifications by JNP the 27/05/2009 15:52:45
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

#ifndef __PROJECTION_REAL_SHEARS_H__
#define __PROJECTION_REAL_SHEARS_H__

/*****************************************************************************
 *    System includes
 ****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <iostream>

/*****************************************************************************
 *    Toolbox defines
 ****************************************************************************/
#include <external/bilib/configs.h>
#include <external/bilib/headers/messagedisplay.h>
#include <data/metadata.h>

/*****************************************************************************
 *    Toolbox includes
 ****************************************************************************/
#include <external/bilib/headers/linearalgebra.h>
#include <external/bilib/headers/getput.h>
#include <external/bilib/headers/getputd.h>
#include <external/bilib/headers/changebasis.h>

/*****************************************************************************
 *    New toolbox includes
 ****************************************************************************/
#include <data/metadata.h>
#include <data/projection.h>
#include <data/image.h>


#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503130808472633361816000000000e-16
#endif

/**@name ProjRealShears Projection library program using Real-Shears */
//@{
/// Structure for holding a volume
class VolumeStruct
{
public:
    const MultidimArray<double> *volume;
    int Xdim;
    Matrix1D<double> angles; // Angles
    Matrix1D<double> shifts; //Shifts

    // Constructor
    VolumeStruct(const MultidimArray<double> &V);
};

/// Make projection
void project_Volume(VolumeStruct &Data, Projection &P, int Ydim, int Xdim,
    double rot, double tilt, double psi);

///Main class of this program
class Projection_real_shears
{
    //--------------- Fields ---------------
    public :
        /// Filename of the projection parameters file (input).
        FileName fn_proj_param;

        /// Selection file with all projections (output).
        FileName fn_sel_file;

        ///Tell if the displaying is active
        bool display;

        ///X Translation parameter
        double shiftX;
        ///Y Translation parameter
        double shiftY;

        ///Projection Xdim
        int proj_Xdim;
        ///Projection Ydim
        int proj_Ydim;

        ///Volume file
        FileName fnPhantom;

        ///Seed name of the projections files
        std::string fnProjectionSeed;

        ///Starting number for the projections files
        int starting;
        ///Extension name of the projections files
        std::string fn_projection_extension;
        ///Angular file
        FileName fn_angle;

        // Volume to project
        Image<double> V;
        ///Content of the angle file
        MetaData DF;
        ///Basics Data for projections
        VolumeStruct *Data;

    //------------------ Functions -------------------
    public :
        ///Main function
        void ROUT_project_real_shears();

        ///Does starting instructions
        void produceSideInfo();

        ///Read input and output file parameters only
        void read(int argc, char **argv);

        /// Usage message. This function shows the way of introducing this parameters.
        void usage();

        ///"Overloaded" function in order to use translation parameters
        void read(const FileName &fn_proj_param); //read_withShift [...]
};
//@}
#endif
