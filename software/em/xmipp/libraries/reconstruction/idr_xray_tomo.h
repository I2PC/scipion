/***************************************************************************
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
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

#ifndef IDR_XRAY_TOMO_H_
#define IDR_XRAY_TOMO_H_

#include "project_xray.h"
#include "recons.h"

/**@defgroup IDRXrayTomo IDR Xray Tomography
   @ingroup ReconsLibrary */

//@{


class ProgIDRXrayTomo: public virtual XmippProgram
{
public:

    /// Metadafile with angles and projection file names
    FileName fnInputProjMD;
    /// Rootname for intermediate/exchange projections metadatas and files
    FileName fnRootInter;
    FileName fnInterProjs;
    FileName fnInterProjsMD;
    FileName fnInterAngles;
    /// Reconstructed output volume file name
    FileName fnOutVol;
    /// Initial volume
    FileName fnStart;

    // Microscope optics parameters
    /// Xray Microscope PSF Parameters
    FileName fnPSF;
    /// Xray PSF Object
    XRayPSF psf;
    /// threshold for psfSlabs
    double psfThr;
    // Input projections sampling
    double sampling;
    /// Number of threads;
    int nThr;

    // IDR Params
    /// Number of iterations
    int itNum;
    /// Relaxation parameter
    Matrix1D<double> lambda_list;


    // Reconstruction method
    enum
    {
        RECONS_ART,
        RECONS_FOURIER,
        RECONS_TOMO3D
    } reconsMethod;

    ProgReconsBase *reconsProgram;

    Image<double> muVol;
    MetaData     projMD;
    MetaData     interProjMD;
    XrayProjPhantom phantom;
    Projection   proj;
    ParallelTaskDistributor * td;


protected:

    virtual void defineParams();
    virtual void readParams();
    void preRun();
    void postRun();

public:

    virtual void run();
    void reconstruct(const FileName &fnProjs, const FileName &fnVol);
    ProgReconsBase * createReconsProgram(const FileName &input, const FileName &output);
};


/** Reconstruct tomogram projections using external tomo3D
 *
 * @param MD Includes tilt angle values and projections file names
 * @param fnOut Reconstructed output volume name
 * @param params Other parameters to be passed to tomo3D
 * @return True if external system call ran right. Otherwise False
 */
int reconsTomo3D(const MetaData& MD, const FileName& fnOut, const String& params = "");

/** Reconstruct tomogram projections using external tomo3D
 *
 * @param fnAngles Text file with angles sequence
 * @param fnProjs  MRC stack file with projections
 * @param fnOut    Reconstructed output volume name
 * @param params   Other parameters to be passed to tomo3D
 * @return True if external system call ran right. Otherwise False
 */
int reconsTomo3D(const FileName& fnAngles, const FileName& fnProjs,
                 const FileName& fnOut, const String& params = "");


//@}
#endif /* IDR_XRAY_TOMO_H_ */
