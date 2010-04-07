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

#ifndef METADATA
// Read arguments ==========================================================
void Prog_centilt_prm::read(int argc, char **argv)
{

    FileName fn_sel;

    // Selfile with untilted images
    fn_sel = getParameter(argc, argv, "-u");
    SFu.read(fn_sel);
    // Selfile with tilted images
    fn_sel = getParameter(argc, argv, "-t");
    SFt.read(fn_sel);
    if (SFu.ImgNo() != SFt.ImgNo())
        REPORT_ERROR(1, "Unequal number of active images in untilted and tilted selfiles");
    // Extension if not to overwrite input images
    oext = getParameter(argc, argv, "-oext", "");
    // Write out document file?
    fn_doc = getParameter(argc, argv, "-doc", "");
    // Maximum shift (discard images that shift more in last iteration)
    max_shift = textToFloat(getParameter(argc, argv, "-max_shift", "0"));
    // Force x-shift to be zero?
    force_x_zero = checkParameter(argc, argv, "-force_x_zero");
    // Perform centering?
    do_center = !checkParameter(argc, argv, "-skip_centering");
    // Perform cosine stretching?
    do_stretch = !checkParameter(argc, argv, "-skip_stretching");
}

// Show ====================================================================
void Prog_centilt_prm::show()
{
    std::cerr << " Selfile untilted images              : " <<  SFu.name() << std::endl;
    std::cerr << " Selfile tilted images                : " <<  SFt.name() << std::endl;
    if (oext != "")
        std::cerr << " Output extension for tilted images   : " << oext << std::endl;
    if (max_shift != 0)
        std::cerr << " Discard images that shift more than  : " << max_shift << std::endl;
    if (fn_doc != "")
        std::cerr << " Output document file (tilted images) : " << fn_doc << std::endl;
    if (force_x_zero)
        std::cerr << " Force x-shift to be zero " << std::endl;
    if (!do_stretch)
        std::cerr << " Skip cosine stretching " << std::endl;
    if (!do_center)
        std::cerr << " Skip centering " << std::endl;

}

// usage ===================================================================
void Prog_centilt_prm::usage()
{
    std::cerr << "Usage:  " << std::endl;
    std::cerr << "  centilt [options]" << std::endl;
    std::cerr << "   -u <selfile>             : Selfile containing untilted images \n"
    << "   -t <selfile>             : Selfile containing tilted images \n"
    << " [ -oext <extension> ]      : For output tilted images; if not to overwrite input\n"
    << " [ -max_shift <float> ]     : Discard images that shift more [pix]\n"
    << " [ -doc <docfile> ]         : write output document file with rotations & translations \n"
    << " [ -force_x_zero ]          : Force x-shift to be zero \n"
    << " [ -skip_stretching ]       : do not perform cosine stretching \n"
    << " [ -skip_centering ]        : do not perform centering, i.e. only modify angles \n"
    << std::endl;
}

// Center one tilted image  =====================================================
bool Prog_centilt_prm::center_tilted_image(const ImageXmipp &Iu, ImageXmipp &It, double &ccf)
{

    Matrix2D<double> A(3, 3), Maux(It()), Mcorr(It());
    int              n_max = -1;
    bool             neighbourhood = true;
    int              imax, jmax, i_actual, j_actual, x_zero = 1;
    double           maxcorr, xmax, ymax, sumcorr;
    float            xshift, yshift, shift;

    if (do_stretch)
    {
        // Cosine stretching, store stretched image in Maux
        A.initIdentity();
        A(0, 0) = COSD(It.Theta());
        Maux.initZeros();
        applyGeometry(Maux, A, It(), IS_INV, DONT_WRAP);
    }
    else Maux = It();

    // Calculate cross-correlation
    correlation_matrix(Maux, Iu(), Mcorr);
    if (force_x_zero)
    {
        FOR_ALL_ELEMENTS_IN_MATRIX2D(Mcorr)
        {
            if (j != 0) MAT_ELEM(Mcorr, i, j) = 0.;
        }
    }
    Mcorr.maxIndex(imax, jmax);
    maxcorr = MAT_ELEM(Mcorr, imax, jmax);

    if (force_x_zero)
    {
        x_zero = 0;
        FOR_ALL_ELEMENTS_IN_MATRIX2D(Mcorr)
        {
            if (j != 0) MAT_ELEM(Mcorr, i, j) = 0.;
        }
    }
    Mcorr.maxIndex(imax, jmax);
    maxcorr = MAT_ELEM(Mcorr, imax, jmax);
    xshift = (float) jmax;
    yshift = (float) imax;

    // Calculate correlation coefficient
    A.initIdentity();
    A(0, 2) = xshift;
    A(1, 2) = yshift;
    Maux.selfApplyGeometry(A, IS_INV, DONT_WRAP);
    Maux.setXmippOrigin();
    ccf = correlation_index(Iu(), Maux);

    if (do_stretch) xshift *= COSD(It.Theta());
    shift = sqrt(xshift * xshift + yshift * yshift);

    if ((max_shift < XMIPP_EQUAL_ACCURACY) || (shift < max_shift))
    {
        // Store shift in the header of the image
	// Take rotation into account for shifts
	It.set_Xoff(-xshift*COSD(It.psi()) - yshift*SIND(It.psi()));
	It.set_Yoff( xshift*SIND(It.psi()) - yshift*COSD(It.psi()));
        It.set_Phi(Iu.psi());
        return true;
    }
    else return false;
}

// Main program  ===============================================================
void Prog_centilt_prm::centilt()
{

    DocFile           DFo;
    FileName          fn_img;
    ImageXmipp        Iu, It;
    Matrix2D<double>  Maux, A(3, 3);
    Matrix1D<double>  dataline(6);
    double            ccf, outside;
    bool              OK;
    int               imgno, barf, n_images, n_discarded;

    if (fn_doc != "")
    {
        DFo.reserve(SFt.ImgNo()*2 + 1);
        DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Corr (6)");
    }

    n_images = SFt.ImgNo();
    n_discarded = 0;
    std::cerr << "  Centering of " << n_images << " tilted images" << std::endl;
    init_progress_bar(n_images);
    barf = XMIPP_MAX(1, (int)(1 + (n_images / 60)));
    imgno = 0;

    while (imgno < n_images)
    {

        SFu.go_beginning();
        SFt.go_beginning();
        SFu.jump(imgno);
        SFt.jump(imgno - n_discarded);
        // Read in untilted image and apply shifts (center) and Phi (align tilt-axis with y-axis)
        Iu.read(SFu.get_current_file());
        Iu().setXmippOrigin();
        Euler_angles2matrix(Iu.Phi(), 0., 0., A);
        A(0, 2) = -Iu.Xoff();
        A(1, 2) = -Iu.Yoff();
        outside = dMij(Iu(), 0, 0);
        Iu().selfApplyGeometry(A, IS_INV, DONT_WRAP, outside);

        // Read in tilted image and apply Psi (align tilt-axis with y-axis) and shifts if present
        It.read(SFt.get_current_file());
        // Store original matrix for later output
        Maux.resize(It());
        Maux = It();
        Euler_angles2matrix(0., 0., It.psi(), A);
        outside = dMij(It(), 0, 0);
        It().selfApplyGeometry(A, IS_INV, DONT_WRAP, outside);
        It().setXmippOrigin();

        if (do_center) OK = center_tilted_image(Iu, It, ccf);
        else
        {
            OK = true;
            ccf = 1.;
            It.set_Phi(Iu.psi());
        }
        if (OK)
        {
            fn_img = It.name();
            if (oext != "")
            {
                fn_img = fn_img.without_extension() + "." + oext;
            }
            // Add line to document file
            if (fn_doc != "")
            {
                dataline(0) = It.Phi();
                dataline(1) = It.Theta();
                dataline(2) = It.psi();
                dataline(3) = It.Xoff();
                dataline(4) = It.Yoff();
                dataline(5) = ccf;
                DFo.append_comment(fn_img);
                DFo.append_data_line(dataline);
            }
            SFt.set_current_filename(fn_img);
            // Re-store original matrix & write out tilted image
            It() = Maux;
            It.write(fn_img);
        }
        else
        {
            SFt.remove(It.name());
            n_discarded++;
        }

        imgno++;
        if (imgno % barf == 0) progress_bar(imgno);
    }
    progress_bar(n_images);
    if (max_shift > 0) std::cerr << "  Discarded " << n_discarded << " tilted images that shifted too much" << std::endl;

    // Write out document file
    if (fn_doc != "") DFo.write(fn_doc);

    // Write out selfile
    fn_img = SFt.name();
    if (oext != "") fn_img = fn_img.insert_before_extension("_" + oext);
    SFt.write(fn_img);


}


#else


// Read arguments ==========================================================
void Prog_centilt_prm::read(int argc, char **argv)
{

    FileName fn_sel;

    // metaData File with untilted images
    fn_sel = getParameter(argc, argv, "-u");
    SFu.read(fn_sel,NULL);
    // metaData file with tilted images
    fn_sel = getParameter(argc, argv, "-t");
    SFt.read(fn_sel,NULL);
    if (SFu.size() != SFt.size())
        REPORT_ERROR(1, "Unequal number of active images in untilted and tilted metadata files");
    // Extension if not to overwrite input images
    oext = getParameter(argc, argv, "-oext", "");
    // Write out document file?
    fn_doc = getParameter(argc, argv, "-doc", "");
    // Maximum shift (discard images that shift more in last iteration)
    max_shift = textToFloat(getParameter(argc, argv, "-max_shift", "0"));
    // Force x-shift to be zero?
    force_x_zero = checkParameter(argc, argv, "-force_x_zero");
    // Perform centering?
    do_center = !checkParameter(argc, argv, "-skip_centering");
    // Perform cosine stretching?
    do_stretch = !checkParameter(argc, argv, "-skip_stretching");
}

// Show ====================================================================
void Prog_centilt_prm::show()
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

// usage ===================================================================
void Prog_centilt_prm::usage()
{
    std::cerr << "Usage:  " << std::endl;
    std::cerr << "  centilt [options]" << std::endl;
    std::cerr << "   -u <metadatafile>             : Metadata File  containing untilted images \n"
  		  	  << "   -t <metadatafile>        : Metadata File  containing tilted images \n"
			  << " [ -oext <extension> ]      : For output tilted images; if not to overwrite input\n"
			  << " [ -doc <metadata> ]        : write output document file with rotations & translations \n"
			  << " [ -max_shift <float> ]     : Discard images that shift more [pix]\n"
			  << " [ -force_x_zero ]          : Force x-shift to be zero \n"
			  << " [ -skip_stretching ]       : do not perform cosine stretching \n"
			  << " [ -skip_centering ]        : do not perform centering, i.e. only modify angles \n"
    << std::endl;
}

// Center one tilted image  =====================================================
bool Prog_centilt_prm::center_tilted_image(const ImageXmipp &Iu, ImageXmipp &It, double &ccf)
{

    Matrix2D<double> A(3, 3), Maux(It()), Mcorr(It());
    int              n_max = -1;
    bool             neighbourhood = true;
    int              imax, jmax, i_actual, j_actual, x_zero = 1;
    double           maxcorr, xmax, ymax, sumcorr;
    float            xshift, yshift, shift;

    if (do_stretch)
    {
        // Cosine stretching, store stretched image in Maux
        A.initIdentity();
        A(0, 0) = COSD(It.Theta());
        Maux.initZeros();
        applyGeometry(Maux, A, It(), IS_INV, DONT_WRAP);
    }
    else Maux = It();

    // Calculate cross-correlation
    correlation_matrix(Maux, Iu(), Mcorr);
    if (force_x_zero)
    {
        FOR_ALL_ELEMENTS_IN_MATRIX2D(Mcorr)
        {
            if (j != 0) MAT_ELEM(Mcorr, i, j) = 0.;
        }
    }
    Mcorr.maxIndex(imax, jmax);
    maxcorr = MAT_ELEM(Mcorr, imax, jmax);

    if (force_x_zero)
    {
        x_zero = 0;
        FOR_ALL_ELEMENTS_IN_MATRIX2D(Mcorr)
        {
            if (j != 0) MAT_ELEM(Mcorr, i, j) = 0.;
        }
    }
    Mcorr.maxIndex(imax, jmax);
    maxcorr = MAT_ELEM(Mcorr, imax, jmax);
    xshift = (float) jmax;
    yshift = (float) imax;

    // Calculate correlation coefficient
    A.initIdentity();
    A(0, 2) = xshift;
    A(1, 2) = yshift;
    Maux.selfApplyGeometry(A, IS_INV, DONT_WRAP);
    Maux.setXmippOrigin();
    ccf = correlation_index(Iu(), Maux);

    if (do_stretch) xshift *= COSD(It.Theta());
    shift = sqrt(xshift * xshift + yshift * yshift);
    
    if ((max_shift < XMIPP_EQUAL_ACCURACY) || (shift < max_shift))
    {
        // Store shift in the header of the image
	// Take rotation into account for shifts
	It.set_Xoff(-xshift*COSD(It.psi()) - yshift*SIND(It.psi()));
	It.set_Yoff( xshift*SIND(It.psi()) - yshift*COSD(It.psi()));
        It.set_Phi(Iu.psi());
        return true;
    }
    else return false;
}

// Main program  ===============================================================
void Prog_centilt_prm::centilt()
{

    FileName          fn_img;
    ImageXmipp        Iu, It;
    Matrix2D<double>  Maux, A(3, 3);
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
    SFu.firstObject();
    SFt.firstObject();
    
    FileName file_name;
        
    do
    {
        // Read in untilted image and apply shifts (center) and Phi (align tilt-axis with y-axis)
        SF.getValue( MDL_IMAGE, file_name);
        Iu.read( file_name);

        Iu().setXmippOrigin();
        Euler_angles2matrix(Iu.Phi(), 0., 0., A);
        A(0, 2) = -Iu.Xoff();
        A(1, 2) = -Iu.Yoff();
        outside = dMij(Iu(), 0, 0);
        Iu().selfApplyGeometry(A, IS_INV, DONT_WRAP, outside);

        // Read in tilted image and apply Psi (align tilt-axis with y-axis) and shifts if present
        SFt.getValue( MDL_IMAGE, file_name);
        It.read( file_name );
        // Store original matrix for later output
        Maux.resize(It());
        Maux = It();
        Euler_angles2matrix(0., 0., It.psi(), A);
        outside = dMij(It(), 0, 0);
        It().selfApplyGeometry(A, IS_INV, DONT_WRAP, outside);
        It().setXmippOrigin();

        if (do_center) OK = center_tilted_image(Iu, It, ccf);
        else
        {
            OK = true;
            ccf = 1.;
            It.set_Phi(Iu.psi());
        }
        if (OK)
        {
            fn_img = It.name();
            if (oext != "")
            {
                fn_img = fn_img.without_extension() + "." + oext;
            }
			SFt.setValue( MDL_ANGLEROT,It.Phi() );
			SFt.setValue( MDL_ANGLETILT,It.Theta() );
			SFt.setValue( MDL_ANGLEPSI,It.psi() );
			SFt.setValue( MDL_SHIFTX,It.Xoff() );
	        SFt.setValue( MDL_SHIFTY,It.Yoff() );
            SFt.setValue( MDL_MAXCC,ccf );
            SFt.setValue( MDL_IMAGE,fn_img );

            // Re-store original matrix & write out tilted image
            It() = Maux;
            It.write(fn_img);
        }
        else
        {
            SFt.setValue( MDL_ENABLED, -1);
            SFt.setValue( MDL_SHIFTX, 0.);
            SFt.setValue( MDL_SHIFTY, 0.);
            SFt.setValue( MDL_ENABLED, -1);
            
        	SFt.setEnabled(-1);
			SFt.setShiftX( 0. );
			SFt.setShiftY( 0. );
			SFt.setMaxCC( 0.);
			n_discarded++;
        }

        imgno++;
        if (imgno % barf == 0) progress_bar(imgno);
        SFu.nextObject();
    }
    while (SFt.nextObject()!= MetaData::NO_MORE_OBJECTS);

    progress_bar(n_images);
    if (max_shift > 0) std::cerr << "  Discarded " << n_discarded << " tilted images that shifted too much" << std::endl;


    // Write out selfile
    fn_img = SFt.getFilename();
    if (oext != "") fn_img = fn_img.insert_before_extension("_" + oext);
    SFt.write(fn_img);
    if (fn_doc != "")
        {
            //delete all discarted images
			SFt.removeObjects(MDL_ENABLED, -1);
			SFt.write(fn_doc);
        }
}

#endif
