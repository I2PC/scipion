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

//#define DEBUG

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

    //Examples
    addExampleLine("To center tilted images allowing a maximum shift of 10 pixels:",false);
    addExampleLine("xmipp_align_tilt_pairs -i class_00001@classes.xmd -o alignedImages.xmd --max_shift 10");

    // Params
    addParamsLine("  -i <metadata>              : Input metadata with untilted and tilted images");
    addParamsLine("  -o <metadata>              : Output metadata file with rotations & translations for 3D reconstruction.");
    addParamsLine("  --ref <file>               : 2D average of the untilted images");
    addParamsLine(" [--max_shift <value=10>]    : Discard images shifting more than a given threshold (in percentage of the image size).");
    addParamsLine("                             :+Set it to 0 for no shift estimate between tilted and untilted images");
    addParamsLine(" [--do_stretch]              : Stretch tilted image to fit into the untilted one.");
    addParamsLine("                             :+Do it only with thin particles");
    addParamsLine(" [--do_not_align_tilted]     : Do not align tilted images to untilted ones.");
    addParamsLine("                             :+Do not align if the quality of the tilted images is low.");
}

//Read params
void ProgAlignTiltPairs::readParams()
{
    fnIn = getParam("-i");
    fnOut = getParam("-o");
    fnRef = getParam("--ref");
    max_shift = getDoubleParam("--max_shift");
    do_stretch = checkParam("--do_stretch");
    do_not_align_tilted = checkParam("--do_not_align_tilted");
}

// Show ====================================================================
void ProgAlignTiltPairs::show()
{
    std::cout
    << "Input metadata:  " << fnIn << std::endl
    << "Output metadata: " << fnOut << std::endl
    << "Reference image: " << fnRef << std::endl;
    if (max_shift != 0)
        std::cout
        << "Discard images that shift more than: " << max_shift << std::endl;
    if (do_stretch)
        std::cout << "Stretching tilted images\n";
}

// Center one tilted image  =====================================================
//#define DEBUG
bool ProgAlignTiltPairs::centerTiltedImage(const MultidimArray<double> &imgRefU,
        bool flip,
        double inPlaneU, double shiftXu, double shiftYu,
        double alphaT, double alphaU,
        double tilt, MultidimArray<double> &imgT, double &shiftX,
        double &shiftY, CorrelationAux &auxCorr)
{
    // Cosine stretching, store stretched image in imgT2DClass
    MultidimArray<double> imgT2DClass;
    Matrix2D<double> E, E2D, Mu2D, A2D;

    if (!do_stretch)
        tilt = flip ? 180. : 0.;

    Euler_angles2matrix(flip ? alphaU : -alphaU, tilt, alphaT, E, true);

    rotation2DMatrix(inPlaneU, Mu2D, true);
    MAT_ELEM(Mu2D,0,2) = flip ? shiftXu : -shiftXu;
    MAT_ELEM(Mu2D,1,2) = -shiftYu;

    if (flip)
    {
        // We need to multiply the first column by -1, because
        // in Xmipp the Mu2D order is: first flip (if necessary), then rot, then shift
        //MAT_ELEM(Mu2D,0,0)*=-1;

        MAT_ELEM(Mu2D,1,0)*=-1;
        MAT_ELEM(Mu2D,2,0)*=-1;

      // We need to multiply the first row by -1, as needed by the RCT theory
        MAT_ELEM(Mu2D,0,1)*=-1;
        MAT_ELEM(Mu2D,0,2)*=-1;
    }

    E2D.initIdentity(3);
    MAT_ELEM(E2D,0,0)=MAT_ELEM(E,0,0);
    MAT_ELEM(E2D,0,1)=MAT_ELEM(E,0,1);
    MAT_ELEM(E2D,1,0)=MAT_ELEM(E,1,0);
    MAT_ELEM(E2D,1,1)=MAT_ELEM(E,1,1);
    A2D = Mu2D * E2D.inv();

    applyGeometry(LINEAR, imgT2DClass, imgT, A2D, IS_NOT_INV, WRAP);

    // Calculate best shift
    int max_shift_pixels=(int)(max_shift/100.0*XSIZE(imgT));
    CorrelationAux aux;
    double corr=bestShift(imgRefU, imgT2DClass, shiftX, shiftY, auxCorr, NULL, max_shift_pixels);
    // double corr=bestShiftRealSpace(imgRefU, imgT2DClass, shiftX, shiftY, NULL, max_shift_pixels, 0.5);

#ifdef DEBUG

    std::cout << "alphaU=" << alphaU << " inplaneU=" << inPlaneU << " tilt=" << tilt << " alphaT=" << alphaT << std::endl;
    std::cout << "Best shift= " << shiftX << " " << shiftY << " corr=" << corr << std::endl;
#endif

    Matrix1D<double> vShift(2);
    XX(vShift)=shiftX;
    YY(vShift)=shiftY;
    Matrix2D<double> Tt, Tt2D;
    translation2DMatrix(vShift, Tt2D,true);
    Tt=A2D.inv()*Tt2D.inv()*A2D;

    // Readjust shift
    shiftX=MAT_ELEM(Tt,0,2);
    shiftY=MAT_ELEM(Tt,1,2);
    double shift = sqrt(shiftX * shiftX + shiftY * shiftY);

//#define DEBUG_2
#ifdef DEBUG_2

    Image<double> save;
    save() = imgT;
    save.write("PPPtilted.xmp");
    save() = imgRefU;
    save.write("PPPuntiltedRef.xmp");
    save() = imgT2DClass;
    save.write("PPPtiltedAdjusted.xmp");
    std::cout << "Corrected shift= " << shiftX << " " << shiftY << std::endl;
    std::cout << "Press any key\n";
    char c = getchar();
#endif

    return (shift < (double)max_shift_pixels) || corr<0;
}
//#undef DEBUG

// Main program  ===============================================================
void ProgAlignTiltPairs::run()
{
    Image<double> imgT;
    MultidimArray<double> Maux;
    Matrix2D<double> A(3, 3);

    MetaData mdIn, mdOut;
    mdIn.read(fnIn);
    mdIn.removeDisabled();

    initProgress(mdIn.size());

    Image<double> imgRef;
    imgRef.read(fnRef);
    imgRef().setXmippOrigin();

    size_t imgno = 0;
    FileName fnTilted, fnUntilted;

    double alphaU, alphaT, tilt, inPlaneU, minusInPlaneU;
    size_t nDiscarded = 0;
    Matrix2D<double> E, Tu, Tup;
    Matrix1D<double> vShift(3);
    vShift.initZeros();
    CorrelationAux auxCorr;
    MDRow row, rowOut;
    bool flip;

    FOR_ALL_OBJECTS_IN_METADATA(mdIn)
    {
        mdIn.getRow(row, __iter.objId);

        row.getValue(MDL_FLIP, flip);
        // Read untilted and tilted images
        row.getValue(MDL_ANGLE_PSI, inPlaneU);
        row.getValue(MDL_IMAGE, fnUntilted);
        row.getValue(MDL_IMAGE_TILTED, fnTilted);
        imgT.read(fnTilted);
        imgT().setXmippOrigin();


#ifdef DEBUG
        {
            Image<double> imgU;
            imgU.read(fnUntilted);
            imgU().setXmippOrigin();
            Image<double> save;
            save() = imgU();
            save.write("PPPuntilted.xmp");
            Matrix2D<double> E;
            save.readApplyGeo(fnUntilted, row);
            save.write("PPPuntiltedAligned.xmp");
            geo2TransformationMatrix(row, E);
            //rotation2DMatrix(-inPlaneU,E);
            //save() = imgU();
            //applyGeometry(LINEAR, save(), imgU(), E, IS_NOT_INV, WRAP);
            //save.write("PPPuntiltedAligned2.xmp");
        }
#endif

        // Get alignment parameters
        double shiftXu, shiftYu;
        row.getValue(MDL_ANGLE_Y, alphaU);
        row.getValue(MDL_ANGLE_Y2, alphaT);
        row.getValue(MDL_ANGLE_TILT, tilt);
        row.getValue(MDL_SHIFT_X, shiftXu);
        row.getValue(MDL_SHIFT_Y, shiftYu);

        if (flip)
        {
            tilt += 180;
            minusInPlaneU = inPlaneU+alphaU;
        }
        else
        	minusInPlaneU=-(inPlaneU+alphaU);

        // Correct untilted alignment
        Euler_angles2matrix( minusInPlaneU, tilt, alphaT, E, true);
        XX(vShift) = shiftXu;
        YY(vShift) = shiftYu;
        translation3DMatrix(vShift, Tu);
        Tup = E * Tu * E.inv();

        // Align tilt and untilted projections
        double shiftX = 0, shiftY = 0;
        bool enable = true;
        if (max_shift > 0)
            enable = centerTiltedImage(imgRef(), flip, inPlaneU, shiftXu, shiftYu,
                                       alphaT, alphaU, tilt, imgT(), shiftX, shiftY, auxCorr);
        if (!enable)
        {
            ++nDiscarded;
            shiftX = shiftY = 0;
        }

        // Write results
        rowOut.setValue(MDL_IMAGE, fnTilted);
        rowOut.setValue(MDL_ANGLE_ROT, minusInPlaneU);
        rowOut.setValue(MDL_ANGLE_TILT, tilt);
        rowOut.setValue(MDL_ANGLE_PSI, alphaT);
        rowOut.setValue(MDL_SHIFT_X, -MAT_ELEM(Tup,0,3) + shiftX);
        rowOut.setValue(MDL_SHIFT_Y, -MAT_ELEM(Tup,1,3) + shiftY);
        rowOut.setValue(MDL_ENABLED, (int) enable);
        mdOut.addRow(rowOut);

        setProgress(++imgno);
    }

    if (verbose && nDiscarded)
        std::cout << "  Discarded " << nDiscarded
        << " tilted images that shifted too much" << std::endl;

    // Write out selfile
    mdOut.write(fnOut);
}

