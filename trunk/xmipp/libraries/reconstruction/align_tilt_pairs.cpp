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
#include "align_tilt_pairs.h"

//Define Program parameters
void ProgAlignTiltPairs::defineParams() {
	//Usage
	addUsageLine("Center the tilted images of all tilted-untilted image pairs.");
	addUsageLine("+This program receives as input a metadata with two sets of images, untilted and tilted.");
	addUsageLine("+Untilted images must have been previously classified and aligned. Then, the tilted images");
	addUsageLine("+are aligned to their untilted versions. The alignment parameters from the classification plus");
	addUsageLine("+the alignment parameters from the tilted-untilted comparison are used to deduce the");
	addUsageLine("+3D alignment parameters for the tilted images.");

	//Examples
	addExampleLine("To center tilted images allowing a maximum shift of 10 pixels:",false);
	addExampleLine("xmipp_align_tilt_pairs -i class_00001@classes.xmd -o alignedImages.xmd --max_shift 10");

	// Params
	addParamsLine("  -i <metadata>              : Input metadata with untilted and tilted images");
	addParamsLine("  -o <metadata>              : Output metadata file with rotations & translations for 3D reconstruction.");
	addParamsLine("  --ref <file>               : 2D average of the untilted images");
	addParamsLine(" [--max_shift <value=0.0>]   : Discard images shifting more than a given threshold (in pixels).");
	addParamsLine("                             :+Set it to 0 for no shift estimate between tilted and untilted images");
	addParamsLine(" [--do_stretch]              : Stretch tilted image to fit into the untilted one.");
	addParamsLine("                             :+Do it only with thin particles");
	addParamsLine(" [--do_not_align_tilted]     : Do not align tilted images to untilted ones.");
	addParamsLine("                             :+Do not align if the quality of the tilted images is low.");
}

//Read params
void ProgAlignTiltPairs::readParams() {
	fnIn = getParam("-i");
	fnOut = getParam("-o");
	fnRef = getParam("--ref");
	max_shift = getDoubleParam("--max_shift");
	do_stretch = checkParam("--do_stretch");
	do_not_align_tilted = checkParam("--do_not_align_tilted");
}

// Show ====================================================================
void ProgAlignTiltPairs::show() {
	std::cout << "Input metadata:  " << fnIn << std::endl << "Output metadata: "
			<< fnOut << std::endl << "Reference image: " << fnRef << std::endl;
	if (max_shift != 0)
		std::cout << "Discard images that shift more than: " << max_shift
				<< std::endl;
	if (do_stretch)
		std::cout << "Stretching tilted images\n";
}

// Center one tilted image  =====================================================
bool ProgAlignTiltPairs::centerTiltedImage(const MultidimArray<double> &imgU,
		double inPlaneU, double alphaT, double alphaU,
		double tilt, MultidimArray<double> &imgT, double &shiftX,
		double &shiftY, CorrelationAux &auxCorr) {
	// Cosine stretching, store stretched image in imgTaux
	double cos_tilt = COSD(tilt);
	MultidimArray<double> imgTaux;
	Matrix2D<double> E, E2D;
	if (!do_stretch)
		tilt=0.;
	Euler_angles2matrix(-alphaU-inPlaneU,tilt,alphaT,E,true);
	E2D.initIdentity(3);
	MAT_ELEM(E2D,0,0)=MAT_ELEM(E,0,0);
	MAT_ELEM(E2D,0,1)=MAT_ELEM(E,0,1);
	MAT_ELEM(E2D,1,0)=MAT_ELEM(E,1,0);
	MAT_ELEM(E2D,1,1)=MAT_ELEM(E,1,1);
	applyGeometry(LINEAR, imgTaux, imgT, E2D, IS_INV, WRAP);

	// Calculate best shift
	CorrelationAux aux;
	bestShift(imgU, imgTaux, shiftX, shiftY, auxCorr);
	std::cout << shiftX << " " << shiftY << std::endl;
	Matrix1D<double> vShift(3);
	XX(vShift)=-shiftX;
	YY(vShift)=-shiftY;
	ZZ(vShift)=0;
	Matrix2D<double> Tt, Ttp;
	translation3DMatrix(vShift, Tt);
	Ttp=E.inv()*Tt*E;

	// Readjust shift
	shiftX=MAT_ELEM(Ttp,0,3);
	shiftY=MAT_ELEM(Ttp,1,3);
	double shift = sqrt(shiftX * shiftX + shiftY * shiftY);

	Image<double> save;
	save() = imgT;
	save.write("PPPtilted.xmp");
	save() = imgU;
	save.write("PPPuntiltedRef.xmp");
	save() = imgTaux;
	save.write("PPPtiltedAdjusted.xmp");
	std::cout << shiftX << " " << shiftY << std::endl;
	std::cout << "Press any key\n";
	char c;
	//std::cin >> c;

	return (shift < max_shift);
}

// Main program  ===============================================================
void ProgAlignTiltPairs::run() {
	Image<double> imgU, imgT;
	MultidimArray<double> Maux;
	Matrix2D<double> A(3, 3);

	MetaData MDin, MDout;
	MDin.read(fnIn);
	MDin.removeDisabled();
	size_t stepBar;
	if (verbose) {
		size_t N = MDin.size();
		init_progress_bar(N);
		stepBar = XMIPP_MAX(1, (int)(1 + (N/60)));
	}

	Image<double> imgRef;
	imgRef.read(fnRef);
	imgRef().setXmippOrigin();

	size_t imgno = 0;
	FileName fnTilted, fnUntilted;

	double alphaU, alphaT, gamma;
	size_t nDiscarded = 0;
	Matrix2D<double> E, Tu, Tup;
	Matrix1D<double> vShift(3);
	vShift.initZeros();
	CorrelationAux auxCorr;
	FOR_ALL_OBJECTS_IN_METADATA(MDin) {
		bool flip;
		MDin.getValue(MDL_FLIP, flip, __iter.objId);

		// Read untilted and tilted images
		double inPlaneU;
		MDin.getValue(MDL_ANGLEPSI, inPlaneU, __iter.objId);
		MDin.getValue(MDL_IMAGE, fnUntilted, __iter.objId);
		imgU.read(fnUntilted);
		imgU().setXmippOrigin();
		MDin.getValue(MDL_IMAGE_TILTED, fnTilted, __iter.objId);
		imgT.read(fnTilted);
		imgT().setXmippOrigin();

		Image<double> save;
		save() = imgU();
		save.write("PPPuntilted.xmp");
		Matrix2D<double> E;
		rotation2DMatrix(-inPlaneU,E);
		applyGeometry(LINEAR, save(), imgU(), E, IS_INV, WRAP);
		save.write("PPPuntiltedAligned.xmp");

		// Get alignment parameters
		double shiftXu, shiftYu;
		MDin.getValue(MDL_ANGLE_Y, alphaU, __iter.objId);
		MDin.getValue(MDL_ANGLE_Y2, alphaT, __iter.objId);
		MDin.getValue(MDL_ANGLETILT, gamma, __iter.objId);
		MDin.getValue(MDL_SHIFTX, shiftXu, __iter.objId);
		MDin.getValue(MDL_SHIFTX, shiftYu, __iter.objId);

		// Correct untilted alignment
		if (flip)
			Euler_angles2matrix(-inPlaneU, gamma, -alphaT, E, true);
		else
			Euler_angles2matrix(-inPlaneU, gamma, alphaT, E, true);
		XX(vShift) = shiftXu;
		YY(vShift) = shiftYu;
		translation3DMatrix(vShift, Tu);
		Tup = E * Tu * E.inv();

		// Align tilt and untilted projections
		double shiftX = 0, shiftY = 0;
		bool enable = true;
		if (max_shift > 0)
			enable = centerTiltedImage(imgRef(), inPlaneU, alphaT, alphaU, gamma, imgT(), shiftX, shiftY, auxCorr);
		if (!enable) {
			++nDiscarded;
			shiftX = shiftY = 0;
		}

		// Write results
		size_t idOut = MDout.addObject();
		MDout.setValue(MDL_IMAGE, fnTilted, idOut);
		if (flip) {
			MDout.setValue(MDL_ANGLEROT, -inPlaneU, idOut);
			MDout.setValue(MDL_ANGLETILT, gamma + 180, idOut);
			MDout.setValue(MDL_ANGLEPSI, alphaT, idOut);
			MDout.setValue(MDL_SHIFTX, -MAT_ELEM(Tup,0,3) + shiftX, idOut);
			MDout.setValue(MDL_SHIFTY, -MAT_ELEM(Tup,1,3) + shiftY, idOut);
		} else {
			MDout.setValue(MDL_ANGLEROT, -inPlaneU, idOut);
			MDout.setValue(MDL_ANGLETILT, gamma, idOut);
			MDout.setValue(MDL_ANGLEPSI, alphaT, idOut);
			MDout.setValue(MDL_SHIFTX, -MAT_ELEM(Tup,0,3) + shiftX, idOut);
			MDout.setValue(MDL_SHIFTY, -MAT_ELEM(Tup,1,3) + shiftY, idOut);
		}
		MDout.setValue(MDL_ENABLED, (int) enable, idOut);

		if (++imgno % stepBar == 0)
			progress_bar(imgno);
	}

	if (verbose && max_shift > 0) {
		progress_bar(MDin.size());
		if (nDiscarded>0)
			std::cout << "  Discarded " << nDiscarded
			<< " tilted images that shifted too much" << std::endl;
	}

	// Write out selfile
	MDout.write(fnOut);
}

