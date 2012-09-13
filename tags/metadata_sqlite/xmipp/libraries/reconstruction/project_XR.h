/***************************************************************************
 * Authors:     AUTHOR_NAME (joton@cnb.csic.es)
 *
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

#ifndef _PROJECTXR_H_
#define _PROJECTXR_H_


#include <data/funcs.h>
#include <data/metadata.h>
#include <data/projection.h>
#include <data/psf_xr.h>

#include <data/transformations.h>


/**@defgroup ProjectionXRProgram project_xr (project for tilt series)
   @ingroup ReconsLibraryPrograms */
//@{
/* Projection XR Program Parameters -------------------------------- */
/** Parameter class for the project program */
class Prog_Project_XR_Parameters
{
public:
    /// Filename with the Projection_Parameters.
    FileName fn_proj_param;
    /// Selection file with all projections
    FileName fn_sel_file;
    /// Filename with the Microscope Parameters.
    FileName fn_psf_xr;
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
        This function shows the way of introducing this parameters. */
    void usage();
};

/* Projection parameters --------------------------------------------------- */
/** Projecting parameters.
    This class reads a set of projection parameters in a file (see
    xmipp_project_xr
    for more information about the file structure) and extract the
    useful information from it.*/
class Projection_XR_Parameters
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
class PROJECT_XR_Side_Info
{
public:
    /// Document File for the projecting angles. Order: rot, tilt, psi
    MetaData        DF;
    /// Phantom Xmipp volume
    Image<double>    phantomVol;
public:
    /** Produce Project Side information.
        This function produce the side information from the project
        program parameters. Basically it loads the phantom.*/
    void produce_Side_Info(const Projection_XR_Parameters &prm,
                           const Prog_Project_XR_Parameters &prog_prm);
};

/* Effectively project ----------------------------------------------------- */
/** Effectively project.
    This is the routine which effectively projects, it needs the projection
    parameters and the side information. The Projection field will keep
    at the end the last projection, this is useful in case you want
    to project only one image, although it is also written to disk.
    The returned number is the total number of projections generated.
    A selection file with all images is also returned.*/
int PROJECT_XR_Effectively_project(
    const Projection_XR_Parameters &prm,
    PROJECT_XR_Side_Info &side, Projection &proj,XmippXRPSF &psf, MetaData &SF);

/** From voxel volumes, off-centered tilt axis.
    This routine projects a volume that is rotating (angle) degrees
    around the axis defined by the two angles (axisRot,axisTilt) and
    that passes through the point raxis. The projection can be futher
    inplane rotated and shifted through the parameters
    (inplaneRot) and (rinplane).

    All vectors involved must be 3D.

    The projection model is rproj=H Rinplane Raxis r +
                                  Rinplane (I-Raxis) raxis + rinplane

    Where Raxis is the 3D rotation matrix given by the axis and
    the angle.

    Off-centered not implemented. Rotations are around volume center
*/
void project_xr_Volume_offCentered(MultidimArray<double> &V,XmippXRPSF &psf, Projection &P,
   int Ydim, int Xdim, double axisRot, double axisTilt,
   const Matrix1D<double> &raxis, double angle, double inplaneRot,
   const Matrix1D<double> &rinplane);

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
int ROUT_XR_project(Prog_Project_XR_Parameters &prm,
    Projection &proj, MetaData &SF);
//@}


#endif /* _PROJECTXR_H_ */
