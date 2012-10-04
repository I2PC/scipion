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

#include "idr_xray_tomo.h"


void ProgIDRXrayTomo::defineParams()
{
    addUsageLine("Iterative Data Refinement applied to the reconstruction of X-ray tomography projections.");
    addSeeAlsoLine("xray_psf_create, xray_import, xray_project");
    //Params
    //    projParam.defineParams(this); // Projection parameters
    addParamsLine("-i <md_file>  : Metadata file with input projections");
    addParamsLine("alias --input");
    addParamsLine("--oroot <rootname=\"idr_xray\">  : Rootname used for refined projections");
    addParamsLine("-o <Vol_file>      : Filename for refined output volume");
    addParamsLine("alias --output");
    addParamsLine(" [-s <Ts>]          : Sampling rate (nm)");
    addParamsLine("alias --sampling_rate;");
    addParamsLine("[--thr <threads=1>]           : Number of concurrent threads.");

    addParamsLine("== Iterations options == ");
    addParamsLine("  [-l <...>]                  : Relaxation factor, by default 0.01 (recommended range 0.0 - 0.1). ");
    addParamsLine("                              : A list of lambda values is also accepted as \"-l lambda0 lambda1 ...\"");
    addParamsLine("  [-n <noit=1>]               : Number of iterations");

    addParamsLine("== Xray options == ");
    addParamsLine("[--psf <psf_param_file=\"\">] : XRay-Microscope parameters file. If not set, then default parameters are chosen.");
    addParamsLine("[--threshold <thr=0.0>]       :+ Normalized threshold relative to maximum of PSF to reduce the volume into slabs");
    //Examples
    addExampleLine("Generating a set of projections using a projection parameter file and two threads:", false);
    addExampleLine("xmipp_xray_project -i volume.vol --oroot images --params projParams.xmd -s 10 --psf psf_560.xmd --thr 2");
    addExampleLine("Generating a single projection at 45 degrees around X axis:", false);
    addExampleLine("xmipp_xray_project -i volume.vol -o image.spi --angles 45 0 90 -s 10 --psf psf_560.xmd");
    addExampleLine("Generating a single projection at 45 degrees around Y axis with default PSF:", false);
    addExampleLine("xmipp_xray_project -i volume.vol -o image.spi --angles 45 90 90 -s 10");
    //Example projection file
    addExampleLine("In the following link you can find an example of projection parameter file:",false);
    addExampleLine(" ",false);
    addExampleLine("http://newxmipp.svn.sourceforge.net/viewvc/newxmipp/trunk/testXmipp/input/tomoProjection.param",false);
    addExampleLine(" ",false);
    addExampleLine("In the following link you can find an example of X-ray microscope parameters file:",false);
    addExampleLine(" ",false);
    addExampleLine("http://newxmipp.svn.sourceforge.net/viewvc/newxmipp/trunk/testXmipp/input/xray_psf.xmd",false);
}

/* Read from command line ================================================== */
void ProgIDRXrayTomo::readParams()
{
    //    fn_proj_param = getParam("-i");
    //    fn_sel_file   = getParam("-o");
//    projParam.readParams(this);

    fn_psf_xr = getParam("--psf");
    dxo  = (checkParam("-s"))? getDoubleParam("-s")*1e-9 : -1 ;
    psfThr = getDoubleParam("--threshold");
    nThr = getIntParam("--thr");

    psf.read(fn_psf_xr);
    psf.verbose = verbose;
    psf.nThr = nThr;
}//readParams


void ProgIDRXrayTomo::run()
{}


