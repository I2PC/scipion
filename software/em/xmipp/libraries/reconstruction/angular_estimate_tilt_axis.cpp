/***************************************************************************
 *
 * Authors:    Carlos Oscar Sorzano          (coss@cnb.csic.es)
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
#include "angular_estimate_tilt_axis.h"
#include <data/micrograph.h>
#include <data/xmipp_filename.h>
#include <string>

//#define DEBUG

//Define Program parameters
void ProgAngularEstimateTiltAxis::defineParams()
{
    //Usage
    addUsageLine("Estimate the tilt axis from a set of corresponding particles");

    // Params
    addParamsLine("  --untilted <posfile>       : Input coordinates in the untilted micrograph");
    addParamsLine("  --tilted <posfile>         : Input coordinates in the tilted micrograph");
    addParamsLine("  -o <metadata>              : Output metadata file with tilt axis angles");
}

//Read params
void ProgAngularEstimateTiltAxis::readParams()
{
    fnUntilted = getParam("--untilted");
    fnTilted = getParam("--tilted");
    fnOut = getParam("-o");
}

// Main program  ===============================================================
void ProgAngularEstimateTiltAxis::run()
{
	MetaData mdU, mdT;
	mdU.read(fnUntilted);
	mdT.read(fnTilted);
	TiltPairAligner aligner;
	MetaData mdOut;
	size_t id;
	int count;
	double alphaU, alphaT, gamma;
	if (mdU.size()==0)
	{
		id=mdOut.addObject();
		mdOut.setValue(MDL_ANGLE_Y,id);
		mdOut.setValue(MDL_ANGLE_Y2,id);
		mdOut.setValue(MDL_ANGLE_TILT,id);
		mdOut.write(fnOut);
		return;
	}
	FOR_ALL_OBJECTS_IN_METADATA2(mdU,mdT)
	{
		int xu, yu, xt, yt;
		mdU.getValue(MDL_XCOOR,xu,__iter.objId);
		mdU.getValue(MDL_YCOOR,yu,__iter.objId);
		mdT.getValue(MDL_XCOOR,xt,__iter2.objId);
		mdT.getValue(MDL_YCOOR,yt,__iter2.objId);

        aligner.addCoordinatePair(xu,yu,xt,yt);
	}

	aligner.calculatePassingMatrix();
	aligner.computeGamma();

	gamma=aligner.gamma;
	aligner.computeAngles(alphaU, alphaT, gamma);

	if (fnOut.exists())
		mdOut.read(fnOut);

	id=mdOut.addObject();
	mdOut.setValue(MDL_ANGLE_Y,alphaU,id);
	mdOut.setValue(MDL_ANGLE_Y2,alphaT,id);
	mdOut.setValue(MDL_ANGLE_TILT,gamma,id);
//
//	mdU.getValue(MDL_IMAGE, imname,id);
//	mdOut.setValue(MDL_IMAGE, imname,id);

	mdOut.write(fnOut);


}

