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
#ifndef _PROG_PROJECTION_TOMOGRAPHY_HH
#define _PROG_PROJECTION_TOMOGRAPHY_HH

#include <data/xmipp_funcs.h>
#include <data/metadata.h>
#include <data/projection.h>
#include <data/xmipp_program.h>

/**@defgroup ProjectionTomographyProgram projectTomograpy (project for tilt series)
   @ingroup ReconsLibrary */
//@{
/* Projection Tomography Program Class-------------------------------- */
class ProgProjectTomography: public XmippProgram
{
public:
    /// Filename with the Projection_Parameters.
    FileName fn_proj_param;
    /// Selection file with all projections
    FileName fn_sel_file;
    /// Projection parameters
    ParametersProjectionTomography projParam;

protected:
    virtual void defineParams();
    virtual void readParams();

public:
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
    virtual void run();
};

/** Project program Side information.
    This class contains side necessary information for the Project program.
    This information can be obtained from the parameters and is basically
    the Xmipp volume or phantom description plus a flag saying which
    of the two is valid. */
class TomoProjectSideInfo
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
    void produceSideInfo(const ParametersProjectionTomography &prm);
};
//@}
#endif
