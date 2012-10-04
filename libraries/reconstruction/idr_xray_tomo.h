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


/**@defgroup IDRXrayTomo IDR Xray Tomography
   @ingroup ReconsLibrary */

//@{


class ProgIDRXrayTomo: public virtual XmippProgram
{
public:

	/// Metadafile with angles and projection file names
	FileName fnInputProj;
	/// Rootname for intermediate/exchange projections metadatas and files
	FileName fnIntermProjs;
	/// Reconstructed output volume file name
	FileName fnOutVol;

    // Microscope optics parameters
    XRayPSF psf;
    /// threshold for psfSlabs
    double psfThr;
    // Input volume sampling
    double dxo;
    /// Number of threads;
    int nThr;

    /// IDR Params
    /// Number of iterations
    int itNum;
    /// Relaxation parameter
    Matrix1D<double> lambda_list;


    XrayProjPhantom phantom;
    Projection   proj;
    MetaData     projMD;
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


//@}
#endif /* IDR_XRAY_TOMO_H_ */
