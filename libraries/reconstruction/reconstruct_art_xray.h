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

#ifndef RECONSTRUCT_ART_XRAY_H_
#define RECONSTRUCT_ART_XRAY_H_

#include "project_xray.h"
#include "recons.h"


/**@defgroup Reconstruction Program ART
   @ingroup ReconsLibrary */
//@{

/* The user interface program should make a call to the run routine.
  */
class ProgReconsXrayART: public ProgReconsBase
{
protected:
    bool            isMpi; // True if prog is mpi version
public:

    /// Selfile with the input images
    FileName fnDoc;
    /// Output filename
    FileName fnOut;
    /// psf Filename
    FileName fnPSF;
    /// Start Volume Filename
    FileName fnStart;
    /// Lambda
    double lambdaART;
    /// Number of iterations
    int Nit;
    /// threshold for psfSlabs
    double psfThr;
    /// Sampling rate
    double sampling;
    /// Metadata with projections info
    MetaData MDin;
    /// Basis function. By default, blobs
    Basis basis;
    /// Microscope parameters
    XRayPSF psf;
    /// Vol with the Igeometrical distribution along specimen volume
    MultidimArray<double> IgeoVol;
    /// Projection X dimension
    int             projXdim;
    /// Projection Y dimension
    int             projYdim;
    // Number of threads to use
    int nThreads;
    ThreadManager * thMgr;

    ProgReconsXrayART()
    {}
    ~ProgReconsXrayART()
    {}

    ///Functions of common reconstruction interface
    void setIO(const FileName &fn_in, const FileName &fn_out)
    {}
    void defineParams();
    void readParams();
    void show();
    void run();
    void preProcess(Image<double> &volVoxels);
    void postProcess();
    //    void preIteration();
    //    void postIteration();
    /** ART single step */
    double singleStep(MultidimArray<double> &muVol, const Projection &projExp,
                      double rot, double tilt, double psi);


}
;
//@}

#endif /* RECONSTRUCT_ART_XRAY_H_ */
