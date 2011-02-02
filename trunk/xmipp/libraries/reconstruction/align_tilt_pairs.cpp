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
    addUsageLine("Center the tilted images of all tilted-untilted image pairs. Besides, the in-plane rotation");
    addUsageLine("of the untilted images Psi, is copied to the first-rotation Phi of the tilted ones. Therefore,");
    addUsageLine("the untilted images should be aligned (using the align2d program), prior to applying centilt.");
    //Examples
    addExampleLine("Typical use of this program:",false);
    addExampleLine("lign_tilt_pairs -u g0u.sel -t g0t.sel -doc centilt.doc -max_shift 10 -oext cen");
    addExampleLine("The input images and selfile are not overwritten, instead new images called imagename.cen are written out, with a corresponding selfile called g0t_cen.sel. This selfile could be used in a subsequent volume reconstruction using ART. Any image that would shift more than 10 pixels will be discarded, i.e. it will not be included in the output selfile, and no .cen image is written out.",false);
    addExampleLine(" It is always a good idea to output document files like centilt.doc, since it allows to restore the angles and shifts from a previous point in the image processing pathway, using the program headerinfo.");
    // Params
    addParamsLine(" -u <untilted_metadata_file> : Input metadata with untilted images.");
    addParamsLine(" -t <tilted_metadata_file>   : Input metadata with tilted images (Untilted and tilted selfile should have");
    addParamsLine("                             : same number of entries, and the same order is assumed for image pairs).");
    addParamsLine(" [-o <output_metadata_file>] : Output metadata file with rotations & translations.");
    addParamsLine(" alias --odoc;");
    addParamsLine(" [--oext <extension>]        : For output tilted images; if not to overwrite input.");
    addParamsLine(" [--max_shift <value=0.0>]   : Discard images which shift more (in pixels).");
    addParamsLine(" [--skip_stretching]         : Default action is to stretch the tilted images in X direction by");
    addParamsLine("                             : 1/cos(tilt_angle), before centering them. Use this option to skip it.");
    addParamsLine(" [--force_x_zero]            : Force x-shift to be zero.");
    addParamsLine(" [-skip_centering]           : Default action is to center tilted images based on cross-correlation with");
    addParamsLine("                             : the untilted pairs. Use this option to skip it.");
}

//Read params
void ProgAlignTiltPairs::readParams()
{
  FileName fn_sel;

  // metaData File with untilted images
  fn_sel = getParam("-u");
  SFu.read(fn_sel);
  // metaData file with tilted images
  fn_sel = getParam("-t");
  SFt.read(fn_sel);
  if (SFu.size() != SFt.size())
      REPORT_ERROR(ERR_MD_OBJECTNUMBER, "Unequal number of active images in untilted and tilted metadata files");
  // Extension if not to overwrite input images
  oext = getParam("--oext");
  // Write out document file?
  fn_doc = getParam("-o");
  // Maximum shift (discard images that shift more in last iteration)
  max_shift = getCountParam("--max_shift");
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
    std::cerr << " MetaData file with untilted images   : " <<  SFu.getFilename() << std::endl;
    std::cerr << " MetaData file with tilted images     : " <<  SFt.getFilename() << std::endl;
    if (oext != "")
        std::cerr << " Output extension for tilted images   : " << oext << std::endl;
    if (fn_doc != "")
        std::cerr << " Output document file (tilted images) : " << fn_doc << std::endl;
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
bool ProgAlignTiltPairs::center_tilted_image(const Image<double> &Iu, Image<double> &It, double &ccf)
{

    Matrix2D<double> A(3, 3);
    MultidimArray<double> Maux(It()), Mcorr(It());
    int              n_max = -1;
    bool             neighbourhood = true;
    int              imax, jmax, i_actual, j_actual, x_zero = 1;
    double           maxcorr, xmax, ymax, sumcorr;
    float            xshift, yshift, shift;

    if (do_stretch)
    {
        // Cosine stretching, store stretched image in Maux
        A.initIdentity();
        A(0, 0) = COSD(It.rot());
        Maux.initZeros();
        applyGeometry(LINEAR, Maux, It(), A, IS_INV, DONT_WRAP);
    }
    else
        Maux = It();

    // Calculate cross-correlation
    correlation_matrix(Maux, Iu(), Mcorr);
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

    if (force_x_zero)
    {
        x_zero = 0;
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
    ccf = correlation_index(Iu(), Maux);

    if (do_stretch)
        xshift *= COSD(It.rot());
    shift = sqrt(xshift * xshift + yshift * yshift);

    if ((max_shift < XMIPP_EQUAL_ACCURACY) || (shift < max_shift))
    {
        // Store shift in the header of the image
        // Take rotation into account for shifts
        It.setShifts(-xshift*COSD(It.psi()) - yshift*SIND(It.psi()),
                     xshift*SIND(It.psi()) - yshift*COSD(It.psi()));
        It.setRot(Iu.psi());
        return true;
    }
    else
        return false;
}

// Main program  ===============================================================
void ProgAlignTiltPairs::run()
{

    FileName          fn_img;
    Image<double>     Iu, It;
    MultidimArray<double> Maux;
    Matrix2D<double>  A(3, 3);
    double            ccf, outside;
    bool              OK;
    int               imgno, barf, n_images, n_discarded;

    n_images = SFt.size();
    n_discarded = 0;
    std::cerr << "  Centering of " << n_images << " tilted images" << std::endl;
    init_progress_bar(n_images);
    barf = XMIPP_MAX(1, (int)(1 + (n_images / 60)));
    imgno = 0;

    //while (imgno < n_images)

    FileName file_name;

    if (SFu.size() != SFt.size())
        REPORT_ERROR(ERR_MD_OBJECTNUMBER, "Error: different metadata lengths in 'centilt'");

    FOR_ALL_OBJECTS_IN_METADATA2(SFu, SFt)
    {
        // Read in untilted image and apply shifts (center) and Phi (align tilt-axis with y-axis)
        SFu.getValue(MDL_IMAGE, file_name, __iter.objId);
        Iu.read(file_name);

        Iu().setXmippOrigin();
        Euler_angles2matrix(Iu.rot(), 0., 0., A);
        A(0, 2) = -Iu.Xoff();
        A(1, 2) = -Iu.Yoff();
        outside = dAij(Iu(), 0, 0);
        selfApplyGeometry(LINEAR, Iu(), A, IS_INV, DONT_WRAP, outside);

        // Read in tilted image and apply Psi (align tilt-axis with y-axis) and shifts if present
        SFt.getValue( MDL_IMAGE, file_name, __iter2.objId);
        It.read(file_name );
        // Store original matrix for later output
        Maux = It();
        Euler_angles2matrix(0., 0., It.psi(), A);
        outside = dAij(It(), 0, 0);
        selfApplyGeometry(LINEAR, It(), A, IS_INV, DONT_WRAP, outside);
        It().setXmippOrigin();

        if (do_center)
            OK = center_tilted_image(Iu, It, ccf);
        else
        {
            OK = true;
            ccf = 1.;
            It.setRot(Iu.psi());
        }
        if (OK)
        {
            fn_img = It.name();
            if (oext != "")
            {
                fn_img = fn_img.withoutExtension() + "." + oext;
            }
            SFt.setValue( MDL_ANGLEROT,It.rot(), __iter2.objId );
            SFt.setValue( MDL_ANGLETILT,It.tilt(), __iter2.objId );
            SFt.setValue( MDL_ANGLEPSI,It.psi(), __iter2.objId );
            SFt.setValue( MDL_SHIFTX,It.Xoff(), __iter2.objId );
            SFt.setValue( MDL_SHIFTY,It.Yoff(), __iter2.objId );
            SFt.setValue( MDL_MAXCC,ccf , __iter2.objId);
            SFt.setValue( MDL_IMAGE,fn_img, __iter2.objId );
            SFt.setValue( MDL_ENABLED, 1, __iter2.objId);

            // Re-store original matrix & write out tilted image
            It() = Maux;
            It.write(fn_img);
        }
        else
        {
            SFt.setValue( MDL_ENABLED, -1, __iter2.objId);
            SFt.setValue( MDL_SHIFTX, 0., __iter2.objId);
            SFt.setValue( MDL_SHIFTY, 0., __iter2.objId);
            SFt.setValue( MDL_ENABLED, -1, __iter2.objId);
            SFt.setValue( MDL_SHIFTX,0. , __iter2.objId);
            SFt.setValue( MDL_SHIFTY,0., __iter2.objId );
            SFt.setValue( MDL_MAXCC,0. , __iter2.objId);
            n_discarded++;
        }

        imgno++;
        if (imgno % barf == 0)
            progress_bar(imgno);
    }

    progress_bar(n_images);
    if (max_shift > 0)
        std::cerr << "  Discarded " << n_discarded << " tilted images that shifted too much" << std::endl;


    // Write out selfile
    fn_img = SFt.getFilename();
    if (oext != "")
        fn_img = fn_img.insertBeforeExtension("_" + oext);
    SFt.write(fn_img);
    if (fn_doc != "")
    {
        //delete all discarted images
        SFt.removeObjects(MDValueEQ(MDL_ENABLED, -1));
        SFt.write(fn_doc);
    }
}

