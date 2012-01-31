/***************************************************************************
 *
 * Authors:    Sjors Scheres                 (scheres@cnb.csic.es)
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
    addUsageLine("+This program receives as input a metadata with attwo sets of images, untilted and tilted. In the first one, a metadata");
    addUsageLine("+is required with values for: *angleRot*,the angle between y-axis and tilt axis taken from micrographs.");
    addUsageLine("+Also the in-plane translations (*shiftX* and *shiftY*) and rotation(*anglePsi*) are expected from a");
    addUsageLine("+previous alignment step. In the tilted metadata should be provided the angle between y-axis and tilt");
    addUsageLine("+axis, *anglePsi* and the tilt angle taken from the experiment.");
    addUsageLine("+As output, the first euler angle *angleRot* of tilted images is set with *anglePsi* of untilted ones.");
    addUsageLine("+The in-plane translations *shiftX* and *shiftY* of tilted images are calculated using cross-correlation");
    addUsageLine("+yielding also the cross-correlation coeficient *maxCC*");
    addUsageLine("+* If the particles under study are disc-shaped and lying flat on the micrograph (i.e. the have a relatively");
    addUsageLine("+  low height (Z) and are more-or-less spherically shaped in XY), this program is expected to give good results");
    addUsageLine("+  using the default options.");
    addUsageLine("+* If the particles have a different shape, the centering may be poor, since especially in the X direction the");
    addUsageLine("+  shift becomes poorly defined. In these cases, three alternative options may be tried:");
    addUsageLine("+    --skip_stretching  will skip the cosine-stretching prior to centering.");
    addUsageLine("+    --force_x_zero will force the shift in the X direction to be zero, and will only center the images in the Y direction");
    addUsageLine("+    --skip_centering will skip the entire centering, so that only the Psi angle of the tilted images will be modified.");

    //Examples
    addExampleLine("To center tilted images g0t.sel allowing a maximum shift of 10 pixels:",false);
    addExampleLine("xmipp_align_tilt_pairs -u g0u.sel -t g0t.sel -o g0t_aligned.doc --max_shift 10");
    // Params
    addParamsLine(" -u <metadata>                : Input metadata with untilted images");
    addParamsLine(" -t <metadata>                : Input metadata with tilted images. Both metadatas should");
    addParamsLine("                              : have same number of entries, and the same order is assumed for image pairs).");
    addParamsLine(" [-o <metadata>]              : Output metadata file with rotations & translations.");
    addParamsLine(" alias --odoc;");
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
    FileName fn_sel;

    // metaData File with untilted images
    fn_sel = getParam("-u");
    mdU.read(fn_sel);
    // metaData file with tilted images
    fn_sel = getParam("-t");
    mdT.read(fn_sel);

    if (mdU.size() != mdT.size())
        REPORT_ERROR(ERR_MD_OBJECTNUMBER, "Unequal number of active images in untilted and tilted metadata files");
    // Write out document file?
    mdOut = getParam("-o");
    // Maximum shift (discard images that shift more in last iteration)
    max_shift = getDoubleParam("--max_shift");
    // Force x-shift to be zero?
    force_x_zero = checkParam("--force_x_zero");
    // Perform centering?
    do_center = !checkParam("--skip_centering");
    // Perform cosine stretching?
    do_stretch = !checkParam("--skip_stretching");
}

// Show ====================================================================
void ProgAlignTiltPairs::show()
{
    std::cerr << " MetaData file with untilted images   : " <<  mdU.getFilename() << std::endl;
    std::cerr << " MetaData file with tilted images     : " <<  mdT.getFilename() << std::endl;
    if (!mdOut.empty())
        std::cerr << " Output document file (tilted images) : " << mdOut << std::endl;
    if (max_shift != 0)
        std::cerr << " Discard images that shift more than  : " << max_shift << std::endl;
    if (force_x_zero)
        std::cerr << " Force x-shift to be zero " << std::endl;
    if (!do_stretch)
        std::cerr << " Skip cosine stretching " << std::endl;
    if (!do_center)
        std::cerr << " Skip centering " << std::endl;

}

// Center one tilted image  =====================================================
bool ProgAlignTiltPairs::centerTiltedImage(const Image<double> &imgU, Image<double> &imgT)
{

    Matrix2D<double> A(3, 3);
    MultidimArray<double> Maux(imgT()), Mcorr(imgT());
    double maxcorr, aux;
    float xshift, yshift, shift;
    int imax, jmax;

    double cos_tilt, imgT_psi;
    rowT.getValue(MDL_ANGLETILT, cos_tilt);
    cos_tilt = COSD(cos_tilt);
    rowT.getValue(MDL_ANGLEPSI, imgT_psi);

    if (do_stretch)
    {
        // Cosine stretching, store stretched image in Maux
        A.initIdentity();

        A(0, 0) = cos_tilt;
        Maux.initZeros();
        applyGeometry(LINEAR, Maux, imgT(), A, IS_INV, DONT_WRAP);
    }
    else
        Maux = imgT();

    // Calculate cross-correlation
    CorrelationAux auxCorr;
    correlation_matrix(Maux, imgU(), Mcorr, auxCorr);

    if (force_x_zero)
    {
        FOR_ALL_ELEMENTS_IN_ARRAY2D(Mcorr)
        {
            if (j != 0)
                A2D_ELEM(Mcorr, i, j) = 0.;
        }
    }

    Mcorr.maxIndex(imax, jmax);
    maxcorr = A2D_ELEM(Mcorr, imax, jmax);
    xshift = (float) jmax;
    yshift = (float) imax;

    // Calculate correlation coefficient
    A.initIdentity();
    A(0, 2) = xshift;
    A(1, 2) = yshift;
    selfApplyGeometry(LINEAR, Maux, A, IS_INV, DONT_WRAP);
    Maux.setXmippOrigin();
    // ccf = correlationIndex(imgU(), Maux);

    if (do_stretch)
        xshift *= cos_tilt;
    shift = sqrt(xshift * xshift + yshift * yshift);

    if ((max_shift < XMIPP_EQUAL_ACCURACY) || (shift < max_shift))
    {
        // Store shift in the header of the image
        // Take rotation into account for shifts
        double cos_psi = COSD(imgT_psi), sin_psi = SIND(imgT_psi);
        rowT.setValue(MDL_SHIFTX, (-xshift * cos_psi - yshift * sin_psi));
        rowT.setValue(MDL_SHIFTY, (xshift * sin_psi - yshift * cos_psi));
        aux = correlationIndex(imgU(), Maux);
        rowT.setValue(MDL_MAXCC, aux);
        return true;
    }

    return false;
}

// Main program  ===============================================================
void ProgAlignTiltPairs::run()
{

    FileName          fn;
    Image<double>     imgU, imgT;
    MultidimArray<double> Maux;
    Matrix2D<double>  A(3, 3);
    double            ccf, outside, aux;
    bool              OK;
    int               imgno, barf, n_images, n_discarded;


    n_images = mdT.size();
    n_discarded = 0;
    if (verbose)
    {
        std::cout << "  Centering of " << n_images << " tilted images" << std::endl;
        init_progress_bar(n_images);
        barf = XMIPP_MAX(1, (int)(1 + (n_images / 60)));
    }

    imgno = 0;
    FileName file_name;

    FOR_ALL_OBJECTS_IN_METADATA2(mdU, mdT)
    {
        idU = __iter.objId;
        idT = __iter2.objId;
        mdU.getRow(rowU, idU);
        mdT.getRow(rowT, idT);

        // Read in untilted image and apply shifts (center) and Phi (align tilt-axis with y-axis)
        rowU.getValue(MDL_IMAGE, fn);
        imgU.read(fn);
        imgU().setXmippOrigin();
        rowU.getValue(MDL_ANGLEROT, aux);
        Euler_angles2matrix(aux, 0., 0., A);
        rowU.getValue(MDL_SHIFTX, aux);
        A(0, 2) = -aux;
        rowU.getValue(MDL_SHIFTY, aux);
        A(1, 2) = -aux;
        outside = dAij(imgU(), 0, 0);
        selfApplyGeometry(LINEAR, imgU(), A, IS_INV, DONT_WRAP, outside);

        // Read in tilted image and apply Psi (align tilt-axis with y-axis) and shifts if present
        rowT.getValue(MDL_IMAGE, fn);
        imgT.read(fn);
        rowT.getValue(MDL_ANGLEPSI, aux);
        Euler_angles2matrix(0., 0., aux, A);
        outside = dAij(imgT(), 0, 0);
        selfApplyGeometry(LINEAR, imgT(), A, IS_INV, DONT_WRAP, outside);
        imgT().setXmippOrigin();

        //Always set Psi of tilted image with Rot of untilted
        rowU.getValue(MDL_ANGLEROT, aux);
        rowT.setValue(MDL_ANGLEPSI, aux);
        rowT.setValue(MDL_MAXCC, 1.);
        rowT.setValue(MDL_ENABLED, 1);

        if (do_center && !centerTiltedImage(imgU, imgT))//something went wrong with this image
        {
            rowT.setValue( MDL_ENABLED, -1);
            rowT.setValue( MDL_SHIFTX, 0. );
            rowT.setValue( MDL_SHIFTY, 0. );
            rowT.setValue( MDL_MAXCC, 0. );
            n_discarded++;
        }

        mdT.setRow(rowT, idT);

        if (++imgno % barf == 0)
            progress_bar(imgno);
    }


    if (verbose && max_shift > 0)
    {
        progress_bar(n_images);
        std::cout << "  Discarded " << n_discarded << " tilted images that shifted too much" << std::endl;
    }

    // Write out selfile
    mdT.write(mdOut);
}

