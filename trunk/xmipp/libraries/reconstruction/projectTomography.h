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
#ifndef _PROG_PROJECTION_TOMOGRAPHY_HH
#define _PROG_PROJECTION_TOMOGRAPHY_HH

#include <data/funcs.h>
#include <data/docfile.h>
#include <data/selfile.h>
#include <data/projection.h>

/**@defgroup ProjectionTomographyProgram projectTomograpy (project for tilt series)
   @ingroup ReconsLibraryPrograms */
//@{
/* Projection Tomogrpahy Program Parameters -------------------------------- */
/** Parameter class for the project program */
class Prog_Project_Tomography_Parameters
{
public:
    /// Filename with the Projection_Parameters.
    FileName fn_proj_param;
    /// Selection file with all projections
    FileName fn_sel_file;
    /// Only create angles, do not project
    bool only_create_angles;

#define TELL_SHOW_ANGLES 0x1
    /** Debugging variable.
        This is a bitwise flag with the following valid labels:
        \\TELL_SHOW_ANGLES: the program shows the angles for each image.*/
    int tell;
public:
    /** Read from a command line.
        An exception might be thrown by any of the internal conversions,
        this would mean that there is an error in the command line and you
        might show a usage message. */
    void read(int argc, char **argv);

    /** Usage message.
        This function shows the way of introdustd::cing this parameters. */
    void usage();
};

/* Projection parameters --------------------------------------------------- */
/** Projecting parameters.
    This class reads a set of projection parameters in a file (see
    xmipp_project_tomography
    for more information about the file structure) and extract the
    useful information from it.*/
class Projection_Tomography_Parameters
{
public:
    /** Phantom filename.
        It must be a Xmipp volume. */
    FileName fnPhantom;
    /// Starting name for all projections
    std::string   fnProjectionSeed;
    /// First projection number. By default, 1.
    int      starting;
    /// Extension for projection filenames. This is optional
    std::string   fn_projection_extension;

    /// Projection Xdim
    int      proj_Xdim;
    /// Projection Ydim
    int      proj_Ydim;

    /// Debugging level. See \ref Prog_Project_Parameters::tell
    int tell;

    /// Rotational angle of the tilt axis
    double axisRot;
    /// Tilt angle of the tilt axis
    double axisTilt;
    /// Offset of the tilt axis
    Matrix1D<double> raxis;
    /// Minimum tilt around this axis
    double tilt0;
    /// Maximum tilt around this axis
    double tiltF;
    /// Step in tilt
    double tiltStep;

    /// Bias to be applied to each pixel grey value */
    double    Npixel_avg;
    /// Standard deviation of the noise to be added to each pixel grey value
    double    Npixel_dev;

    /// Bias to apply to the image center
    double    Ncenter_avg;
    /// Standard deviation of the image center
    double    Ncenter_dev;

    /// Bias to apply to the angles
    double    Nangle_avg;
    /// Standard deviation of the angles
    double    Nangle_dev;
public:
    /** Read projection parameters from a file.
        An exception is thrown if the file is not found or any of the
        parameters is not found in the right place.*/
    void read(const FileName &fn_proj_param);
};

/** Project program Side information.
    This class contains side necessary information for the Project program.
    This information can be obtained from the parameters and is basically
    the Xmipp volume or phantom description plus a flag saying which
    of the two is valid. */
class PROJECT_Tomography_Side_Info
{
public:
    /// Document File for the projecting angles. Order: rot, tilt, psi
    DocFile        DF;
    /// Phantom Xmipp volume
    VolumeXmipp    phantomVol;
public:
    /** Produce Project Side information.
        This function produce the side information from the project
        program parameters. Basically it loads the phantom.*/
    void produce_Side_Info(const Projection_Tomography_Parameters &prm,
                           const Prog_Project_Tomography_Parameters &prog_prm);
};

/* Effectively project ----------------------------------------------------- */
/** Effectively project.
    This is the routine which effectively projects, it needs the projection
    parameters and the side information. The Projection field will keep
    at the end the last projection, this is useful in case you want
    to project only one image, although it is also written to disk.
    The returned number is the total number of projections generated.
    A selection file with all images is also returned.*/
int PROJECT_Tomography_Effectively_project(
    const Projection_Tomography_Parameters &prm,
    PROJECT_Tomography_Side_Info &side, Projection &proj, SelFile &SF);

/* Main routine ------------------------------------------------------------ */
/** Main Project routine.
    Generate a set of projections given the projection parameters.
    This is the main projecting routine. This function generates a set
    of projections according to the projecting parameters defined.
    The projections are written to disk.

    The Projection field will keep
    at the end the last projection, this is useful in case you want
    to project only one image, although it is also written to disk.
    The returned number is the total number of projections generated.
    A selection file with all images is also returned (and saved if any
    name has been given in the parameters).*/
int ROUT_Tomography_project(Prog_Project_Tomography_Parameters &prm,
    Projection &proj, SelFile &SF);
//@}
#endif
