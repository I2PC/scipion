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

#include <data/psf_xr.h>
#include <data/xmipp_program.h>
#include <data/xmipp_threads.h>

/**@defgroup ProjectionXRProgram project_xr (project for tilt series)
   @ingroup ReconsLibrary */
//@{

/// Data struct to be passed to threads
struct XrayThread
{
    XRayPSF                *psf;
    MultidimArray<double> *vol;
    Image<double>         *imOut;
    ParallelTaskDistributor * td;
};

/** Project program Side information.
    This class contains side necessary information for the Project program.
    This information can be obtained from the parameters and is basically
    the Xmipp volume or phantom description plus a flag saying which
    of the two is valid. */
class XrayProjPhantom
{
public:
    /// Phantom Xmipp volume
    ImageGeneric          iniVol;
    MultidimArray<double> rotVol;

public:
    /** Produce Project Side information.
        This function produce the side information from the project
        program parameters. Basically it loads the phantom.*/
    void read(const ParametersProjectionTomography &prm);
};

/* Projection XR Program -------------------------------- */
/** Program class for the project program */
class ProgXrayProject: public virtual XmippProgram
{

public:
    /// Filename with the Projection_Parameters.
    FileName fn_proj_param;
    /// Selection file with all projections
    FileName fn_sel_file;
    /// Filename with the Microscope Parameters.
    FileName fn_psf_xr;
    // Projection parameters
    ParametersProjectionTomography projParam;
    // Microscope optics parameters
    XRayPSF psf;
    // Input volume sampling
    double dxo;
    /// Number of threads;
    int nThr;

    XrayProjPhantom phantom;
    Projection   proj;
    MetaData     projMD;
    XrayThread * mainDataThread;
    ParallelTaskDistributor * td;

protected:
    virtual void defineParams();
    virtual void readParams();

public:

    virtual void run();

protected:

    void preRun();
    void postRun();
};

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
void XrayProjectVolumeOffCentered(XrayProjPhantom &side, XRayPSF &psf, Projection &P,
                                  int Ydim, int Xdim, int  idxSlice = 1);

/// Thread Job to generate an X-ray microscope projection
void threadXrayProject(ThreadArgument &thArg);
//@}

#endif /* _PROJECTXR_H_ */
