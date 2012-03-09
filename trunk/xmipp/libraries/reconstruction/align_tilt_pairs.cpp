/***************************************************************************
 *
 * Authors:    Sjors Scheres                 (scheres@cnb.csic.es)
 *             Carlos Oscar Sorzano          (coss@cnb.csic.es)
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
#include "align_tilt_pairs.h"

//Define Program parameters
void ProgAlignTiltPairs::defineParams()
{
    //Usage
    addUsageLine("Center the tilted images of all tilted-untilted image pairs.");
    addUsageLine("+This program receives as input a metadata with two sets of images, untilted and tilted.");
    addUsageLine("+Untilted images must have been previously classified and aligned. Then, the tilted images");
    addUsageLine("+are aligned to their untilted versions. The alignment parameters from the classification plus");
    addUsageLine("+the alignment parameters from the tilted-untilted comparison are used to deduce the");
    addUsageLine("+3D alignment parameters for the tilted images.");
    addUsageLine("+* If the particles under study are disc-shaped and lying flat on the micrograph (i.e. the have a relatively");
    addUsageLine("+  low height (Z) and are more-or-less spherically shaped in XY), this program is expected to give good results");
    addUsageLine("+  using the default options.");
    addUsageLine("+* If the particles have a different shape, the centering may be poor, since especially in the X direction the");
    addUsageLine("+  shift becomes poorly defined. In these cases, three alternative options may be tried:");
    addUsageLine("+    --skip_stretching  will skip the cosine-stretching prior to centering.");
    addUsageLine("+    --force_x_zero will force the shift in the X direction to be zero, and will only center the images in the Y direction");
    addUsageLine("+    --skip_centering will skip the entire centering, so that only the Psi angle of the tilted images will be modified.");

    //Examples
    addExampleLine("To center tilted images allowing a maximum shift of 10 pixels:",false);
    addExampleLine("xmipp_align_tilt_pairs -i class_00001@classes.xmd -o alignedImages.xmd --max_shift 10");

    // Params
    addParamsLine("  -i <metadata>              : Input metadata with untilted and tilted images");
    addParamsLine("  -o <metadata>              : Output metadata file with rotations & translations for 3D reconstruction.");
    addParamsLine(" [--max_shift <value=0.0>]   : Discard images which shift more (in pixels).");
    addParamsLine(" [--force_x_zero]            : Force x-shift to be zero.");
    addParamsLine(" [--skip_stretching]         : Default action is to stretch the tilted images in X direction by");
    addParamsLine("                             : 1/cos(tilt_angle), before centering them. Use this option to skip it.");
    addParamsLine(" [--skip_centering]          : Default action is to center tilted images based on cross-correlation with");
    addParamsLine("                             : the untilted pairs. Use this option to skip it.");
}

//Read params
void ProgAlignTiltPairs::readParams()
{
    fnIn = getParam("-i");
    fnOut = getParam("-o");
    max_shift = getDoubleParam("--max_shift");
    force_x_zero = checkParam("--force_x_zero");
    do_center = !checkParam("--skip_centering");
    do_stretch = !checkParam("--skip_stretching");
}

// Show ====================================================================
void ProgAlignTiltPairs::show()
{
    std::cerr
    << "Input metadata: " <<  fnIn << std::endl
    << "Output metadata: " <<  fnOut << std::endl;
    if (max_shift != 0)
        std::cerr << "Discard images that shift more than: " << max_shift << std::endl;
    if (force_x_zero)
        std::cerr << "Force x-shift to be zero " << std::endl;
    if (!do_stretch)
        std::cerr << "Skip cosine stretching " << std::endl;
    if (!do_center)
        std::cerr << "Skip centering " << std::endl;
}

// Center one tilted image  =====================================================
bool ProgAlignTiltPairs::centerTiltedImage(const MultidimArray<double> &imgU,
		double tilt, MultidimArray<double> &imgT,
		double alphaT,
		double &shiftX, double &shiftY, double &crossCorrelation)
{
    Image<double> save;
    save()=imgU;
    save.write("PPPuntilted.xmp");
    save()=imgT;
    save.write("PPPtilted.xmp");

    // Cosine stretching, store stretched image in imgTaux
    double cos_tilt = COSD(tilt);
    MultidimArray<double> imgTaux;
    if (do_stretch)
    {
        Matrix2D<double> A(3, 3);
        A.initIdentity();
        A(0, 0) = cos_tilt;
    	std::cout << "Tilt=" << tilt << " A=\n" << A << std::endl;
        applyGeometry(LINEAR, imgTaux, imgT, A, IS_INV, WRAP);
    }
    else
    	imgTaux = imgT;

    save()=imgTaux;
    save.write("PPPtiltedStretched.xmp");

    // Calculate cross-correlation
    MultidimArray<double> Mcorr;
    CorrelationAux auxCorr;
    correlation_matrix(imgU, imgTaux, Mcorr, auxCorr);

    // Look for maximum
    if (force_x_zero)
        FOR_ALL_ELEMENTS_IN_ARRAY2D(Mcorr)
            if (j != 0)
                A2D_ELEM(Mcorr, i, j) = 0.;

    int imax,jmax;
    Mcorr.maxIndex(imax, jmax);
    crossCorrelation = A2D_ELEM(Mcorr, imax, jmax);
    shiftX = jmax;
    shiftY = imax;

    // Readjust shift
    if (do_stretch)
    	shiftX *= cos_tilt;
    double shift = sqrt(shiftX * shiftX + shiftY * shiftY);
    std::cout << "shiftX=" << shiftX << " shiftY=" << shiftY << std::endl;
    std::cout << "Press any key\n";
    char c; std::cin >> c;

    if ((max_shift==0) || (shift < max_shift))
    {
    	// Accept shift
        double cos_alphaT = COSD(alphaT), sin_alphaT = SIND(alphaT);
        double correctedShiftX=(-shiftX * cos_alphaT - shiftY * sin_alphaT);
        double correctedShiftY=( shiftX * sin_alphaT - shiftY * cos_alphaT);
        shiftX=correctedShiftX;
        shiftY=correctedShiftY;
        return true;
    }
    else
    	// Reject shift
    	return false;
}

// Main program  ===============================================================
void ProgAlignTiltPairs::run()
{
    Image<double>     imgU, imgT;
    MultidimArray<double> Maux;
    Matrix2D<double>  A(3, 3);

    MetaData MDin, MDout;
    MDin.read(fnIn);
    size_t stepBar;
    if (verbose)
    {
    	size_t N=MDin.size();
        init_progress_bar(N);
        stepBar = XMIPP_MAX(1, (int)(1 + (N/60)));
    }

    size_t imgno = 0;
    FileName fnTilted, fnUntilted;

    double alphaU, alphaT, gamma;
    MultidimArray<double> verticalU, verticalT;
    size_t nDiscarded=0;
	Matrix2D<double> E, Tu, Tup;
	Matrix1D<double> vShift(3);
	vShift.initZeros();
    FOR_ALL_OBJECTS_IN_METADATA(MDin)
    {
    	bool flip;
    	MDin.getValue(MDL_FLIP,flip,__iter.objId);
    	if (flip)
    		continue;

    	// Read untilted and tilted images
    	double inPlaneU;
    	MDin.getValue(MDL_ANGLEPSI,inPlaneU,__iter.objId);
    	MDin.getValue(MDL_IMAGE,fnUntilted,__iter.objId);
    	imgU.read(fnUntilted);
    	imgU().setXmippOrigin();
    	MDin.getValue(MDL_IMAGE_TILTED,fnTilted,__iter.objId);
    	imgT.read(fnTilted);
    	imgT().setXmippOrigin();
    	std::cout << "Reading " << fnTilted << std::endl;

    	// Get alignment parameters
    	double shiftXu, shiftYu;
    	MDin.getValue(MDL_ANGLE_Y,alphaU,__iter.objId);
    	MDin.getValue(MDL_ANGLE_Y2,alphaT,__iter.objId);
    	MDin.getValue(MDL_ANGLETILT,gamma,__iter.objId);
    	MDin.getValue(MDL_SHIFTX,shiftXu,__iter.objId);
    	MDin.getValue(MDL_SHIFTX,shiftYu,__iter.objId);

    	// Correct untilted alignment
    	Euler_angles2matrix(-inPlaneU,gamma,alphaT,E,true);
    	XX(vShift)=shiftXu;
    	YY(vShift)=shiftYu;
    	translation3DMatrix(vShift,Tu);
    	Tup=E*Tu*E.inv();

    	// Align tilt axis with Y
    	rotate(BSPLINE3,verticalU,imgU(),-alphaU,'Z',WRAP);
    	rotate(BSPLINE3,verticalT,imgT(),-alphaT,'Z',WRAP);

    	// Align tilt and untilted projections
    	double shiftX=-MAT_ELEM(Tup,0,3), shiftY=-MAT_ELEM(Tup,1,3), crossCorrelation=0;
    	bool enable=true;
        if (do_center)
        	enable=centerTiltedImage(verticalU, gamma, verticalT, alphaT, shiftX, shiftY, crossCorrelation);
        if (!enable)
        	++nDiscarded;

        // Write results
        size_t idOut=MDout.addObject();
        MDout.setValue(MDL_IMAGE,fnTilted,idOut);
        MDout.setValue(MDL_ANGLEROT,-inPlaneU,idOut);
        MDout.setValue(MDL_ANGLETILT,gamma,idOut);
        MDout.setValue(MDL_ANGLEPSI,alphaT,idOut);
        MDout.setValue(MDL_SHIFTX,shiftX,idOut);
        MDout.setValue(MDL_SHIFTY,shiftY,idOut);

        if (++imgno % stepBar == 0)
            progress_bar(imgno);
    }

    if (verbose && max_shift > 0)
    {
        progress_bar(MDin.size());
        std::cout << "  Discarded " << nDiscarded << " tilted images that shifted too much" << std::endl;
    }

    // Write out selfile
    MDout.write(fnOut);
}

